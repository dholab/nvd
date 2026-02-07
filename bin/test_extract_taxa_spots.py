"""
Characterization tests for extract_taxa_spots.py.

These tests capture behavior before/after refactoring to use py_nvd.taxonomy.
Focus on:
- is_in_target_lineage() membership checks
- resolve_taxa() name-to-taxid resolution
- End-to-end filtering
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


class TestIsInTargetLineage:
    """Tests for is_in_target_lineage() membership checks."""

    def test_exact_match(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Taxid exactly matching target should return True."""
        from extract_taxa_spots import is_in_target_lineage

        with taxonomy.open() as tax:
            # 9606 is in target set {9606}
            result = is_in_target_lineage(
                tax_id=9606,
                target_taxa={9606},
                tax=tax,
                include_children=False,
            )
            assert result is True

    def test_no_match(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Taxid not in target should return False."""
        from extract_taxa_spots import is_in_target_lineage

        with taxonomy.open() as tax:
            # 9598 (chimp) is not in target set {9606}
            result = is_in_target_lineage(
                tax_id=9598,
                target_taxa={9606},
                tax=tax,
                include_children=False,
            )
            assert result is False

    def test_include_children_ancestor_match(
        self, test_taxonomy_sqlite: Path, mock_taxonomy_open,
    ):
        """With include_children, taxid under target ancestor should match."""
        from extract_taxa_spots import is_in_target_lineage

        with taxonomy.open() as tax:
            # 9606 (human) is under 207598 (Homininae)
            result = is_in_target_lineage(
                tax_id=9606,
                target_taxa={207598},  # Homininae
                tax=tax,
                include_children=True,
            )
            assert result is True

    def test_include_children_no_match(
        self, test_taxonomy_sqlite: Path, mock_taxonomy_open,
    ):
        """With include_children, taxid not under target should not match."""
        from extract_taxa_spots import is_in_target_lineage

        with taxonomy.open() as tax:
            # 9606 (human) is not under 9596 (Pan genus)
            result = is_in_target_lineage(
                tax_id=9606,
                target_taxa={9596},  # Pan
                tax=tax,
                include_children=True,
            )
            assert result is False


class TestResolveTaxa:
    """Tests for resolve_taxa() name-to-taxid resolution."""

    def test_resolves_valid_name(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Valid scientific name should resolve to taxid."""
        from extract_taxa_spots import resolve_taxa

        with taxonomy.open() as tax:
            result = resolve_taxa(["Homo sapiens"], tax)
            assert 9606 in result

    def test_unknown_name_skipped(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Unknown name should be skipped with warning."""
        from extract_taxa_spots import resolve_taxa

        with taxonomy.open() as tax:
            result = resolve_taxa(["NotARealTaxon"], tax)
            assert len(result) == 0

    def test_multiple_names(self, test_taxonomy_sqlite: Path, mock_taxonomy_open):
        """Multiple names should all be resolved."""
        from extract_taxa_spots import resolve_taxa

        with taxonomy.open() as tax:
            result = resolve_taxa(["Homo sapiens", "Pan troglodytes"], tax)
            assert 9606 in result
            assert 9598 in result


class TestEndToEnd:
    """End-to-end tests for extract_taxa_spots.py."""

    def test_filters_matching_spots(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Spots with matching taxa should be included in output."""
        import sys

        from extract_taxa_spots import main

        # Create input file with hits
        input_file = tmp_path / "hits.txt"
        # Format: spot_id<tab>taxid (or taxidxcount)
        # spot1 has only human hits, spot2 has only chimp hits
        input_file.write_text("spot1\t9606\nspot2\t9598\n")

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "extract_taxa_spots.py",
                "--gettax_sqlite_path",
                "ignored",
                "--hits_file",
                str(input_file),
                "--output_file",
                str(output_file),
                "--taxa",
                "Homo sapiens",
                "--stringency",
                "0.5",
                "--include_children",
            ]
            main()
        finally:
            sys.argv = original_argv

        assert output_file.exists()
        content = output_file.read_text()
        # spot1 should be included (human matches "Homo sapiens")
        assert "spot1" in content
        # spot2 should NOT be included (chimp doesn't match "Homo sapiens")
        assert "spot2" not in content

    def test_empty_input_handled(
        self,
        test_taxonomy_sqlite: Path,
        mock_taxonomy_open,
        tmp_path: Path,
    ):
        """Empty input should produce empty output."""
        import sys

        from extract_taxa_spots import main

        input_file = tmp_path / "empty.txt"
        input_file.write_text("")

        output_file = tmp_path / "output.txt"

        original_argv = sys.argv
        try:
            sys.argv = [
                "extract_taxa_spots.py",
                "--gettax_sqlite_path",
                "ignored",
                "--hits_file",
                str(input_file),
                "--output_file",
                str(output_file),
                "--taxa",
                "Homo sapiens",
            ]
            main()
        finally:
            sys.argv = original_argv

        # Should complete without error
        # Output file may or may not exist depending on implementation
