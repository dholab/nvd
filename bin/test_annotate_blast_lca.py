"""
Characterization tests for annotate_blast_lca.py.

These tests capture behavior before/after refactoring to use py_nvd.taxonomy.
Focus on:
- _find_lca_list() LCA calculation
- filter_blast_hits() quality thresholds
- End-to-end output format
"""

import sys
from pathlib import Path

import polars as pl
import polars.exceptions
import pytest

from annotate_blast_lca import (
    LcaParams,
    _find_lca_list,
    filter_blast_hits,
    main,
    parse_args,
)
from py_nvd import taxonomy


@pytest.fixture
def test_taxonomy_sqlite(tmp_path: Path) -> Path:
    """
    Create minimal taxonomy SQLite for testing.

    Taxonomy tree:
        1 (root)
        └── 131567 (cellular organisms)
            └── 2759 (Eukaryota)
                └── 33208 (Metazoa)
                    └── 9604 (Hominidae)
                        └── 207598 (Homininae)
                            ├── 9605 (Homo)
                            │   └── 9606 (Homo sapiens)
                            └── 9596 (Pan)
                                └── 9598 (Pan troglodytes)
    """
    taxdump_dir = tmp_path / "taxdump"
    taxdump_dir.mkdir()

    nodes_content = """\
1\t|\t1\t|\tno rank\t|\t
131567\t|\t1\t|\tno rank\t|\t
2759\t|\t131567\t|\tsuperkingdom\t|\t
33208\t|\t2759\t|\tkingdom\t|\t
9604\t|\t33208\t|\tfamily\t|\t
207598\t|\t9604\t|\tsubfamily\t|\t
9605\t|\t207598\t|\tgenus\t|\t
9606\t|\t9605\t|\tspecies\t|\t
9596\t|\t207598\t|\tgenus\t|\t
9598\t|\t9596\t|\tspecies\t|\t
"""
    (taxdump_dir / "nodes.dmp").write_text(nodes_content)

    names_content = """\
1\t|\troot\t|\t\t|\tscientific name\t|
131567\t|\tcellular organisms\t|\t\t|\tscientific name\t|
2759\t|\tEukaryota\t|\t\t|\tscientific name\t|
33208\t|\tMetazoa\t|\t\t|\tscientific name\t|
9604\t|\tHominidae\t|\t\t|\tscientific name\t|
207598\t|\tHomininae\t|\t\t|\tscientific name\t|
9605\t|\tHomo\t|\t\t|\tscientific name\t|
9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|
9596\t|\tPan\t|\t\t|\tscientific name\t|
9598\t|\tPan troglodytes\t|\t\t|\tscientific name\t|
"""
    (taxdump_dir / "names.dmp").write_text(names_content)

    merged_content = """\
12345\t|\t9606\t|
"""
    (taxdump_dir / "merged.dmp").write_text(merged_content)

    sqlite_path = taxonomy._build_sqlite_from_dmp(taxdump_dir)
    return sqlite_path


@pytest.fixture
def mock_taxonomy_open(test_taxonomy_sqlite: Path, monkeypatch: pytest.MonkeyPatch):
    """Patch taxonomy.open() to use our test SQLite."""
    taxdump_dir = test_taxonomy_sqlite.parent
    monkeypatch.setattr(taxonomy, "_ensure_taxdump", lambda _=None: taxdump_dir)


class TestFindLcaList:
    """Tests for _find_lca_list() LCA calculation."""

    def test_single_taxid(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Single taxid should return itself."""
        with taxonomy.open() as tax:
            tx = tax.taxopy_db
            result = _find_lca_list([9606], tx)
            assert result == 9606

    def test_same_taxid_twice(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Same taxid twice should return itself."""
        with taxonomy.open() as tax:
            tx = tax.taxopy_db
            result = _find_lca_list([9606, 9606], tx)
            assert result == 9606

    def test_human_chimp_lca(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Human and Chimp should return Homininae."""
        with taxonomy.open() as tax:
            tx = tax.taxopy_db
            result = _find_lca_list([9606, 9598], tx)
            assert result == 207598  # Homininae

    def test_empty_list(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Empty list should return None."""
        with taxonomy.open() as tax:
            tx = tax.taxopy_db
            result = _find_lca_list([], tx)
            assert result is None

    def test_none_input(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """None input should return None."""
        with taxonomy.open() as tax:
            tx = tax.taxopy_db
            result = _find_lca_list(None, tx)
            assert result is None


class TestFilterBlastHits:
    """Tests for filter_blast_hits() quality thresholds."""

    def test_adds_smax_column(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Filter should add Smax column."""
        # Create minimal LazyFrame
        df = pl.LazyFrame(
            {
                "task": ["megablast"],
                "sample": ["sample1"],
                "qseqid": ["contig1"],
                "length": [100],
                "pident": [95.0],
                "evalue": [1e-20],
                "bitscore": [200.0],
                "staxids": [9606],
            },
        )

        params = LcaParams()
        result = filter_blast_hits(df, params).collect()

        assert "Smax" in result.columns
        assert "near_tie" in result.columns


class TestEndToEnd:
    """End-to-end tests for annotate_blast_lca.py."""

    def test_produces_output_with_adjusted_taxid(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Test that script produces output with adjusted_taxid column."""
        # Create minimal BLAST input TSV
        input_file = tmp_path / "blast_results.txt"
        input_content = """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig1\t500\tref|NC_001\tHuman gene\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
"""
        input_file.write_text(input_content)

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "-i",
                str(input_file),
                "-o",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        assert output_file.exists()
        content = output_file.read_text()
        # Should have header with adjusted_taxid column
        assert "adjusted_taxid" in content

    def test_help_works(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Test that --help doesn't crash."""
        original_argv = sys.argv
        try:
            sys.argv = ["annotate_blast_lca.py", "--help"]
            with pytest.raises(SystemExit) as exc_info:
                parse_args()
            assert exc_info.value.code == 0
        finally:
            sys.argv = original_argv


class TestLCAGoldenResults:
    """
    Golden tests for LCA calculation with known expected results.

    These tests use realistic BLAST-like data and verify that the LCA
    calculation produces the expected taxonomy assignments. The expected
    values are hand-calculated based on the test taxonomy tree:

        1 (root)
        └── 131567 (cellular organisms)
            └── 2759 (Eukaryota)
                └── 33208 (Metazoa)
                    └── 9604 (Hominidae)
                        └── 207598 (Homininae)
                            ├── 9605 (Homo)
                            │   └── 9606 (Homo sapiens)
                            └── 9596 (Pan)
                                └── 9598 (Pan troglodytes)
    """

    @pytest.fixture
    def golden_blast_data(self, tmp_path: Path) -> Path:
        """
        Create BLAST-like TSV with multiple test scenarios.

        Test cases:
        1. contig1: All hits to same species (9606) -> 9606 (dominant)
        2. contig2: Human (9606) + Chimp (9598) equal scores -> 207598 (LCA)
        3. contig3: Human dominant (90% bitscore) -> 9606 (dominant)
        4. contig4: Near-tie Human/Chimp -> 207598 (LCA fallback)
        5. contig5: Single hit -> 9606 (trivial)
        6. contig6: Invalid taxid mixed with valid -> 9606 (skip invalid)
        """
        blast_file = tmp_path / "golden_blast.txt"
        # Format: task sample qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids
        content = """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig1\t500\tref|NC_001\tHuman gene 1\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
megablast\tsample1\tcontig1\t500\tref|NC_002\tHuman gene 2\t448\t99.3\t1e-99\t795\tHomo sapiens\t9606
megablast\tsample1\tcontig1\t500\tref|NC_003\tHuman gene 3\t445\t99.0\t1e-98\t790\tHomo sapiens\t9606
megablast\tsample1\tcontig2\t500\tref|NC_004\tHuman gene\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
megablast\tsample1\tcontig2\t500\tref|NC_005\tChimp gene\t450\t99.5\t1e-100\t800\tPan troglodytes\t9598
megablast\tsample1\tcontig3\t500\tref|NC_006\tHuman gene\t450\t99.5\t1e-100\t900\tHomo sapiens\t9606
megablast\tsample1\tcontig3\t500\tref|NC_007\tChimp gene\t200\t85.0\t1e-20\t100\tPan troglodytes\t9598
megablast\tsample1\tcontig4\t500\tref|NC_008\tHuman gene\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
megablast\tsample1\tcontig4\t500\tref|NC_009\tChimp gene\t448\t99.3\t1e-99\t795\tPan troglodytes\t9598
megablast\tsample1\tcontig5\t500\tref|NC_010\tHuman gene\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
megablast\tsample1\tcontig6\t500\tref|NC_011\tHuman gene\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
megablast\tsample1\tcontig6\t500\tref|NC_012\tUnknown\t200\t85.0\t1e-20\t100\tUnknown\t999999999
"""
        blast_file.write_text(content)
        return blast_file

    # Expected results for each contig
    # Format: {qseqid: expected_adjusted_taxid}
    EXPECTED_TAXIDS = {
        "contig1": 9606,  # All Human -> Human (dominant)
        "contig2": 207598,  # Human + Chimp equal -> Homininae (LCA)
        "contig3": 9606,  # Human dominant (900 vs 100 bitscore) -> Human
        "contig4": 207598,  # Human + Chimp near-tie -> Homininae (LCA)
        "contig5": 9606,  # Single hit -> Human
        "contig6": 9606,  # Human + invalid -> Human (invalid skipped)
    }

    def test_golden_lca_results(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        golden_blast_data: Path,
        tmp_path: Path,
    ):
        """Verify LCA produces expected results for all golden test cases."""
        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "-i",
                str(golden_blast_data),
                "-o",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        assert output_file.exists(), "Output file was not created"

        # Parse output and verify each contig's adjusted_taxid
        df = pl.read_csv(output_file, separator="\t")

        assert "adjusted_taxid" in df.columns, "Output missing adjusted_taxid column"
        assert "qseqid" in df.columns, "Output missing qseqid column"

        # Get unique qseqid -> adjusted_taxid mappings
        results = dict(df.select(["qseqid", "adjusted_taxid"]).unique().iter_rows())

        # Verify each expected result
        for qseqid, expected_taxid in self.EXPECTED_TAXIDS.items():
            assert qseqid in results, f"Missing result for {qseqid}"
            actual_taxid = results[qseqid]
            assert actual_taxid == expected_taxid, (
                f"Wrong taxid for {qseqid}: expected {expected_taxid}, got {actual_taxid}"
            )

    def test_dominant_hit_threshold(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """
        Verify dominant hit detection respects the >=80% bitscore threshold.

        When one taxid has >=80% of the total bitscore, it should be used
        directly instead of calculating LCA.

        Note: The threshold is >= (not >), so exactly 80% uses dominant path.
        """
        # Create data where Human has exactly 80% of bitscore (boundary case)
        blast_file = tmp_path / "threshold_test.txt"
        content = """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_80pct\t500\tref|NC_001\tHuman\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
megablast\tsample1\tcontig_80pct\t500\tref|NC_002\tChimp\t200\t85.0\t1e-20\t200\tPan troglodytes\t9598
megablast\tsample1\tcontig_79pct\t500\tref|NC_003\tHuman\t450\t99.5\t1e-100\t790\tHomo sapiens\t9606
megablast\tsample1\tcontig_79pct\t500\tref|NC_004\tChimp\t200\t85.0\t1e-20\t210\tPan troglodytes\t9598
"""
        blast_file.write_text(content)

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "-i",
                str(blast_file),
                "-o",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        df = pl.read_csv(output_file, separator="\t")
        results = dict(df.select(["qseqid", "adjusted_taxid"]).unique().iter_rows())

        # 80% exactly: 800/(800+200) = 0.80 -> uses dominant (threshold is >=)
        assert results["contig_80pct"] == 9606, "80% threshold should use dominant hit"

        # 79%: 790/(790+210) = 0.79 -> ALSO uses dominant due to delta_s_window
        # The near-tie window (default 5.0) means both hits are considered,
        # but the dominant calculation still applies per-taxid bitscore sums.
        # This documents actual behavior - may need review if unexpected.
        assert results["contig_79pct"] == 9606, (
            "79% also uses dominant - delta_s_window affects which hits are considered"
        )

    def test_merged_taxid_in_blast_results(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """
        Document behavior with merged (deprecated) taxids.

        IMPORTANT FINDING: The LCA script does NOT resolve merged taxids.
        It passes through the original taxid from the BLAST results.

        The test taxonomy has 12345 -> 9606 (Homo sapiens) as a merged taxid,
        but the script returns 12345 as-is.

        This may be intentional (preserve original BLAST output) or a gap
        that should be addressed. Documenting current behavior.
        """
        blast_file = tmp_path / "merged_test.txt"
        content = """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_merged\t500\tref|NC_001\tOld Human ID\t450\t99.5\t1e-100\t800\tHomo sapiens\t12345
"""
        blast_file.write_text(content)

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "-i",
                str(blast_file),
                "-o",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        df = pl.read_csv(output_file, separator="\t")
        results = dict(df.select(["qseqid", "adjusted_taxid"]).unique().iter_rows())

        # CURRENT BEHAVIOR: Merged taxid is NOT resolved
        # The script returns 12345, not 9606
        # TODO: Consider whether this should be changed
        assert results["contig_merged"] == 12345, (
            "Current behavior: merged taxids are NOT resolved in LCA script"
        )

    def test_semicolon_delimited_taxids_not_supported(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """
        Document that semicolon-delimited taxids are NOT handled by this script.

        IMPORTANT: The LCA script expects staxids to already be integers.
        Semicolon-delimited taxids (e.g., "9606;9598") should be exploded
        upstream by select_top_blast_hits.py before reaching this script.

        Pipeline flow:
        1. Raw BLAST output (may have semicolon-delimited taxids)
        2. select_top_blast_hits.py (explodes semicolons into separate rows)
        3. annotate_blast_results.py (adds lineage)
        4. annotate_blast_lca.py (THIS SCRIPT - expects integer taxids)

        This test documents that passing semicolon-delimited taxids directly
        to this script will cause a schema error.
        """
        blast_file = tmp_path / "semicolon_test.txt"
        content = """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_multi\t500\tref|NC_001\tPrimate gene\t450\t99.5\t1e-100\t800\tHomo sapiens;Pan troglodytes\t9606;9598
"""
        blast_file.write_text(content)

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "-i",
                str(blast_file),
                "-o",
                str(output_file),
            ]
            # This should fail with a schema error because staxids contains
            # semicolon-delimited strings, not integers
            with pytest.raises(polars.exceptions.InvalidOperationError):
                main()
        finally:
            sys.argv = original_argv

    def test_empty_input_creates_empty_output(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Empty input file creates empty output file without error."""
        # Create empty input file
        blast_file = tmp_path / "empty.txt"
        blast_file.touch()

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "-i",
                str(blast_file),
                "-o",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        # Output file should exist and be empty
        assert output_file.exists()
        assert output_file.stat().st_size == 0
