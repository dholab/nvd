"""
Characterization tests for annotate_blast_lca.py.

These tests capture behavior before/after refactoring to use py_nvd.taxonomy.
Focus on:
- relative bitscore selection for assignment references
- dominant versus LCA assignment
- exact taxid/name/rank correspondence
- valid empty-output sentinels
"""

import csv
import sys
from pathlib import Path

import polars as pl
import pytest
from annotate_blast_lca import (
    OUTPUT_COLUMNS,
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
                        └── 5000 (unranked test group)
                            ├── 5001 (Test species one)
                            └── 5002 (Test species two)
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
5000\t|\t9604\t|\tno rank\t|\t
5001\t|\t5000\t|\tspecies\t|\t
5002\t|\t5000\t|\tspecies\t|\t
6001\t|\t1\t|\tspecies\t|\t
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
5000\t|\tUnranked test group\t|\t\t|\tscientific name\t|
5001\t|\tTest species one\t|\t\t|\tscientific name\t|
5002\t|\tTest species two\t|\t\t|\tscientific name\t|
6001\t|\tSparse-lineage test species\t|\t\t|\tscientific name\t|
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

    def test_taxdump_dir(*_args: object, **_kwargs: object) -> Path:
        return taxdump_dir

    monkeypatch.setenv("NVD_TAXONOMY_DB", str(taxdump_dir))
    monkeypatch.setattr(
        taxonomy,
        "_ensure_taxdump",
        test_taxdump_dir,
    )


class TestEndToEnd:
    """End-to-end tests for annotate_blast_lca.py."""

    def test_taxonomy_policy_args_parse(self) -> None:
        """Taxonomy mode and max-age arguments are accepted by the CLI parser."""
        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "-i",
                "input.tsv",
                "-o",
                "output.tsv",
                "--taxonomy-mode",
                "read_only",
                "--taxonomy-max-age-days",
                "7",
            ]
            args = parse_args()
        finally:
            sys.argv = original_argv

        assert args.taxonomy_mode == "read_only"
        assert args.taxonomy_max_age_days == 7

    def test_close_scoring_taxids_resolve_to_their_lca(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ) -> None:
        """References scoring at least 95% of the best hit determine the LCA."""
        blast_file = tmp_path / "close_scoring_taxids.tsv"
        blast_file.write_text(
            """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_edge\t216\tref|human-best\tHuman best\t214\t100\t1e-100\t396\tHomo sapiens\t9606
megablast\tsample1\tcontig_edge\t216\tref|chimp-close\tChimp close\t214\t99.5\t1e-99\t390\tPan troglodytes\t9598
megablast\tsample1\tcontig_edge\t216\tref|sparse-close\tSparse close\t10\t20.0\t1\t390\tSparse-lineage test species\t6001
megablast\tsample1\tcontig_edge\t216\tref|test-close\tTest close\t214\t99.0\t1e-98\t385\tTest species one\t5001
megablast\tsample1\tcontig_edge\t216\tref|human-distant\tHuman distant\t214\t98.0\t1e-90\t374\tHomo sapiens\t9606
""",
        )
        output_file = tmp_path / "output.tsv"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "--input-file",
                str(blast_file),
                "--output-file",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        output = pl.read_csv(output_file, separator="\t")
        assert output.height == 5
        assert set(output["sseqid"]) == {
            "ref|human-best",
            "ref|chimp-close",
            "ref|sparse-close",
            "ref|test-close",
            "ref|human-distant",
        }
        [assignment] = (
            output.select("adjusted_taxid", "adjustment_method")
            .unique()
            .iter_rows(named=True)
        )
        assert assignment == {
            "adjusted_taxid": 1,
            "adjustment_method": "lca",
        }


class TestTaxonomicConsensus:
    """Behavioral tests for taxonomic consensus output."""

    def test_assignment_score_boundary_controls_dominant_versus_lca(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """A reference at 95% participates while one below 95% does not."""
        blast_file = tmp_path / "threshold_test.txt"
        content = """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_at_boundary\t500\tref|NC_001\tHuman\t450\t99.5\t1e-100\t100\tHomo sapiens\t9606
megablast\tsample1\tcontig_at_boundary\t500\tref|NC_002\tChimp\t450\t99.0\t1e-90\t95\tPan troglodytes\t9598
megablast\tsample1\tcontig_below_boundary\t500\tref|NC_003\tHuman\t450\t99.5\t1e-100\t100\tHomo sapiens\t9606
megablast\tsample1\tcontig_below_boundary\t500\tref|NC_004\tChimp\t450\t99.0\t1e-90\t94.9\tPan troglodytes\t9598
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

        results = {
            row["qseqid"]: (row["adjusted_taxid"], row["adjustment_method"])
            for row in (
                pl.read_csv(output_file, separator="\t")
                .select("qseqid", "adjusted_taxid", "adjustment_method")
                .unique()
                .iter_rows(named=True)
            )
        }

        assert results == {
            "contig_at_boundary": (207598, "lca"),
            "contig_below_boundary": (9606, "dominant"),
        }

    def test_merged_taxid_in_blast_results(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Merged taxids are canonicalized before method selection."""
        blast_file = tmp_path / "merged_test.txt"
        content = """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_merged\t500\tref|NC_001\tOld Human ID\t450\t99.5\t1e-100\t800\tHomo sapiens\t12345
megablast\tsample1\tcontig_merged\t500\tref|NC_001\tOld Human ID\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
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

        output = pl.read_csv(output_file, separator="\t")
        assert output.height == 1
        assert output["staxids"].to_list() == [9606]
        [assignment] = (
            output
            .select("adjusted_taxid", "adjustment_method")
            .unique()
            .iter_rows(named=True)
        )
        assert assignment == {
            "adjusted_taxid": 9606,
            "adjustment_method": "dominant",
        }

    def test_invalid_taxids_do_not_create_false_lca_assignments(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ) -> None:
        """Method selection uses only taxids available in the taxonomy."""
        blast_file = tmp_path / "invalid_taxid.tsv"
        blast_file.write_text(
            """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_invalid\t500\tref|valid\tHuman\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
megablast\tsample1\tcontig_invalid\t500\tref|invalid\tUnknown\t450\t99.5\t1e-100\t800\tUnknown\t999999999
""",
        )
        output_file = tmp_path / "output.tsv"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "--input-file",
                str(blast_file),
                "--output-file",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        output = pl.read_csv(output_file, separator="\t")
        assert output.filter(pl.col("sseqid") == "ref|valid")["staxids"].item() == 9606
        assert output.filter(pl.col("sseqid") == "ref|invalid")[
            "staxids"
        ].item() is None
        [assignment] = (
            output
            .select("adjusted_taxid", "adjustment_method")
            .unique()
            .iter_rows(named=True)
        )
        assert assignment == {
            "adjusted_taxid": 9606,
            "adjustment_method": "dominant",
        }

    def test_unranked_lca_keeps_matching_taxid_name_and_rank(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ) -> None:
        """Strict-LCA metadata describes the exact LCA taxid."""
        blast_file = tmp_path / "unranked_lca.tsv"
        blast_file.write_text(
            """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_unranked\t500\tref|one\tOne\t450\t99.5\t1e-100\t800\tTest species one\t5001
megablast\tsample1\tcontig_unranked\t500\tref|two\tTwo\t450\t99.5\t1e-100\t800\tTest species two\t5002
""",
        )
        output_file = tmp_path / "output.tsv"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "--input-file",
                str(blast_file),
                "--output-file",
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        [assignment] = (
            pl.read_csv(output_file, separator="\t")
            .select(
                "adjusted_taxid",
                "adjusted_taxid_name",
                "adjusted_taxid_rank",
                "adjustment_method",
            )
            .unique()
            .iter_rows(named=True)
        )
        assert assignment == {
            "adjusted_taxid": 5000,
            "adjusted_taxid_name": "Unranked test group",
            "adjusted_taxid_rank": "no rank",
            "adjustment_method": "lca",
        }

    def test_query_with_no_canonical_assignment_taxid_fails(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ) -> None:
        blast_file = tmp_path / "unavailable_taxids.tsv"
        blast_file.write_text(
            """\
task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
megablast\tsample1\tcontig_unavailable\t500\tref|invalid\tUnknown\t450\t99.5\t1e-100\t800\tUnknown\t999999999
""",
        )
        output_file = tmp_path / "output.tsv"

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_lca.py",
                "--input-file",
                str(blast_file),
                "--output-file",
                str(output_file),
            ]
            with pytest.raises(RuntimeError, match="contig_unavailable"):
                main()
        finally:
            sys.argv = original_argv

        assert not output_file.exists()

    def test_empty_input_creates_empty_output(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Empty input file creates a header-only output without error."""
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

        assert output_file.exists()
        with output_file.open(newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))

        assert rows == [OUTPUT_COLUMNS]

    def test_header_only_input_creates_header_only_output(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Header-only merged BLAST tables are valid no-hit sentinels."""
        blast_file = tmp_path / "header_only.txt"
        blast_file.write_text(
            "task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids\trank\n",
        )

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

        assert output_file.exists()
        with output_file.open(newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))

        assert rows == [OUTPUT_COLUMNS]
