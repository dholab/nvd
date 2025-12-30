"""
Characterization tests for annotate_blast_results.py.

These tests capture current behavior before refactoring to use py_nvd.taxonomy.
The goal is to verify that refactoring preserves output format and lineage strings.

Test strategy:
- Build a minimal SQLite taxonomy database matching gettax.sqlite schema
- Pass the SQLite path to the script (current implementation needs it)
- After refactoring, the --sqlite_cache arg will be ignored but tests still pass
"""

from pathlib import Path

import pytest
from py_nvd import taxonomy


@pytest.fixture
def test_taxonomy_sqlite(tmp_path: Path) -> Path:
    """
    Create minimal taxonomy SQLite for testing.

    Builds a SQLite database with the same schema as gettax.sqlite,
    containing a small taxonomy tree:
        1 (root)
        └── 131567 (cellular organisms)
            └── 2759 (Eukaryota)
                └── 33208 (Metazoa)
                    └── 9604 (Hominidae)
                        └── 207598 (Homininae)
                            └── 9605 (Homo)
                                └── 9606 (Homo sapiens)

    Returns path to the SQLite file.
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
"""
    (taxdump_dir / "names.dmp").write_text(names_content)

    merged_content = """\
12345\t|\t9606\t|
"""
    (taxdump_dir / "merged.dmp").write_text(merged_content)

    # Build SQLite from .dmp files
    sqlite_path = taxonomy._build_sqlite_from_dmp(taxdump_dir)
    return sqlite_path


@pytest.fixture
def mock_taxonomy_open(test_taxonomy_sqlite: Path, monkeypatch: pytest.MonkeyPatch):
    """
    Patch taxonomy.open() to use our test SQLite.

    This fixture is used AFTER refactoring. The test_taxonomy_sqlite fixture
    provides the SQLite path for the CURRENT implementation.
    """
    taxdump_dir = test_taxonomy_sqlite.parent
    monkeypatch.setattr(taxonomy, "_ensure_taxdump", lambda _=None: taxdump_dir)


@pytest.fixture
def blast_input_file(tmp_path: Path) -> Path:
    """Create minimal BLAST input TSV."""
    input_file = tmp_path / "blast_results.txt"
    # Header + one data row with taxid 9606 (Homo sapiens)
    content = """\
qseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
contig_1\t500\tref|NC_001\tHuman gene\t450\t99.5\t1e-100\t800\tHomo sapiens\t9606
"""
    input_file.write_text(content)
    return input_file


@pytest.fixture
def empty_blast_input(tmp_path: Path) -> Path:
    """Create empty BLAST input file."""
    input_file = tmp_path / "empty_blast.txt"
    input_file.write_text("")
    return input_file


class TestAnnotateBlastResults:
    """Tests for annotate_blast_results.py behavior."""

    def test_output_contains_lineage_string(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        blast_input_file: Path,
        tmp_path: Path,
    ):
        """Test that output contains expected lineage strings for known taxids."""
        output_file = tmp_path / "output.tsv"

        # Import and run the main function
        # We need to simulate command-line args
        import sys

        from annotate_blast_results import main

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_results.py",
                "--input_file",
                str(blast_input_file),
                "--output_file",
                str(output_file),
                "--sample_name",
                "test_sample",
                "--task",
                "megablast",
            ]
            main()
        finally:
            sys.argv = original_argv

        # Verify output
        assert output_file.exists()
        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Should have header + 1 data row
        assert len(lines) == 2

        # Check header
        header = lines[0].split("\t")
        assert header[0] == "task"
        assert header[1] == "sample"
        assert header[-1] == "rank"  # lineage column

        # Check data row
        data = lines[1].split("\t")
        assert data[0] == "megablast"
        assert data[1] == "test_sample"

        # Check lineage contains expected ranks
        lineage = data[-1]
        assert "species:Homo sapiens" in lineage
        assert "genus:Homo" in lineage
        assert "family:Hominidae" in lineage

    def test_empty_input_produces_empty_output(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        empty_blast_input: Path,
        tmp_path: Path,
    ):
        """Test that empty input produces empty output file."""
        output_file = tmp_path / "output.tsv"

        import sys

        from annotate_blast_results import main

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_results.py",
                "--input_file",
                str(empty_blast_input),
                "--output_file",
                str(output_file),
                "--sample_name",
                "test_sample",
                "--task",
                "megablast",
            ]
            main()
        finally:
            sys.argv = original_argv

        # Empty input should produce empty output
        assert output_file.exists()
        assert output_file.read_text() == ""

    def test_multiple_taxids_semicolon_separated(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Test handling of multiple taxids in staxids column."""
        input_file = tmp_path / "multi_taxid.txt"
        # Two taxids separated by semicolon
        content = """\
qseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids
contig_1\t500\tref|NC_001\tMultiple hits\t450\t99.5\t1e-100\t800\tHomo sapiens;Homo\t9606;9605
"""
        input_file.write_text(content)

        output_file = tmp_path / "output.tsv"

        import sys

        from annotate_blast_results import main

        original_argv = sys.argv
        try:
            sys.argv = [
                "annotate_blast_results.py",
                "--input_file",
                str(input_file),
                "--output_file",
                str(output_file),
                "--sample_name",
                "test_sample",
                "--task",
                "megablast",
            ]
            main()
        finally:
            sys.argv = original_argv

        content = output_file.read_text()
        lines = content.strip().split("\n")
        data = lines[1].split("\t")
        lineage = data[-1]

        # Should have lineages for both taxids joined by "; "
        # Each taxid gets its own lineage string, then they're joined
        assert "Homo sapiens" in lineage
        assert "Homo" in lineage
