"""
Characterization tests for hits_to_report.py.

These tests capture behavior after refactoring to use py_nvd.taxonomy.
Focus on:
- deduce_tax_id() LCA logic (now uses tax.find_lca())
- End-to-end output format
"""

from pathlib import Path

import pytest
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


class TestDeduceTaxId:
    """Tests for deduce_tax_id() LCA logic."""

    def test_single_taxid_returns_itself(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
    ):
        """Single taxid should return itself."""
        from hits_to_report import deduce_tax_id

        with taxonomy.open() as tax:
            result = deduce_tax_id({9606}, tax)
            assert result == 9606

    def test_same_taxid_twice_returns_itself(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
    ):
        """Same taxid appearing twice should return itself."""
        from hits_to_report import deduce_tax_id

        with taxonomy.open() as tax:
            # Sets deduplicate, so {9606, 9606} == {9606}
            result = deduce_tax_id({9606}, tax)
            assert result == 9606

    def test_human_chimp_returns_homininae(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
    ):
        """Human (9606) and Chimp (9598) should return Homininae (207598)."""
        from hits_to_report import deduce_tax_id

        with taxonomy.open() as tax:
            result = deduce_tax_id({9606, 9598}, tax)
            # LCA of Human and Chimp is Homininae (207598)
            assert result == 207598

    def test_homo_and_human_returns_homo(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
    ):
        """Homo (9605) and Homo sapiens (9606) should return Homo (9605)."""
        from hits_to_report import deduce_tax_id

        with taxonomy.open() as tax:
            result = deduce_tax_id({9605, 9606}, tax)
            # 9606 is child of 9605, so LCA is 9605
            assert result == 9605


class TestEndToEnd:
    """End-to-end tests for hits_to_report.py."""

    def test_produces_formatted_output(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Test that script produces expected tree output format."""
        import sys

        from hits_to_report import main

        # Create input file with hits
        input_file = tmp_path / "hits.txt"
        # Format: spot_id<tab>taxid
        input_file.write_text("spot1\t9606\nspot2\t9606\nspot3\t9598\n")

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "hits_to_report.py",
                str(input_file),
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        assert output_file.exists()
        content = output_file.read_text()
        # Should contain taxonomy names
        assert "Homo sapiens" in content or "Homininae" in content

    def test_empty_input_produces_empty_tree(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Test that empty input produces empty output."""
        import sys

        from hits_to_report import main

        input_file = tmp_path / "empty.txt"
        input_file.write_text("")

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "hits_to_report.py",
                str(input_file),
                str(output_file),
            ]
            main()
        finally:
            sys.argv = original_argv

        assert output_file.exists()
        # Empty input produces empty output
        content = output_file.read_text()
        assert content == ""
