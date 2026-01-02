"""
Tests for py_nvd.taxonomy module.

Uses minimal test fixtures to avoid downloading real NCBI data during tests.
"""

import sqlite3
from pathlib import Path

import pytest

from py_nvd import taxonomy
from py_nvd.db import get_taxdump_dir


@pytest.fixture
def minimal_taxdump(tmp_path: Path) -> Path:
    """
    Create minimal taxdump for testing.

    Creates a small taxonomy tree:
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

    Also includes a merged taxid: 12345 -> 9606
    """
    taxdump_dir = tmp_path / "taxdump"
    taxdump_dir.mkdir()

    # nodes.dmp: tax_id | parent_tax_id | rank | ...
    # Format: tax_id<tab>|<tab>parent_tax_id<tab>|<tab>rank<tab>|<tab>...
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

    # names.dmp: tax_id | name | unique_name | name_class
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

    # merged.dmp: old_tax_id | new_tax_id
    merged_content = """\
12345\t|\t9606\t|
"""
    (taxdump_dir / "merged.dmp").write_text(merged_content)

    return taxdump_dir


@pytest.fixture
def taxdump_with_merge_chain(tmp_path: Path) -> Path:
    """
    Create taxdump with multi-level merge chain for testing.

    Same taxonomy tree as minimal_taxdump, but with merge chain:
        11111 -> 22222 -> 9606 (Homo sapiens)

    Also includes a merge to non-existent taxid:
        99999 -> 88888 (88888 doesn't exist)
    """
    taxdump_dir = tmp_path / "taxdump_chain"
    taxdump_dir.mkdir()

    # Same nodes as minimal_taxdump
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

    # Same names as minimal_taxdump
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

    # Multi-level merge chain + merge to non-existent
    merged_content = """\
11111\t|\t22222\t|
22222\t|\t9606\t|
99999\t|\t88888\t|
"""
    (taxdump_dir / "merged.dmp").write_text(merged_content)

    return taxdump_dir


@pytest.fixture
def taxdump_with_cyclic_merge(tmp_path: Path) -> Path:
    """
    Create taxdump with cyclic merge for testing cycle detection.

    Merge cycle: 11111 -> 22222 -> 11111
    """
    taxdump_dir = tmp_path / "taxdump_cycle"
    taxdump_dir.mkdir()

    # Minimal nodes
    nodes_content = """\
1\t|\t1\t|\tno rank\t|\t
9606\t|\t1\t|\tspecies\t|\t
"""
    (taxdump_dir / "nodes.dmp").write_text(nodes_content)

    names_content = """\
1\t|\troot\t|\t\t|\tscientific name\t|
9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|
"""
    (taxdump_dir / "names.dmp").write_text(names_content)

    # Cyclic merge
    merged_content = """\
11111\t|\t22222\t|
22222\t|\t11111\t|
"""
    (taxdump_dir / "merged.dmp").write_text(merged_content)

    return taxdump_dir


@pytest.fixture
def test_taxonomy(minimal_taxdump: Path, monkeypatch: pytest.MonkeyPatch):
    """Provide TaxonomyDB using minimal test data."""
    # Build SQLite from test .dmp files
    taxonomy._build_sqlite_from_dmp(minimal_taxdump)

    # Patch _ensure_taxdump to return our test directory
    monkeypatch.setattr(taxonomy, "_ensure_taxdump", lambda _=None: minimal_taxdump)

    with taxonomy.open() as tax:
        yield tax


class TestTaxonomyDB:
    """Tests for TaxonomyDB class."""

    def test_get_taxon_basic(self, test_taxonomy):
        """Test basic taxon lookup."""
        taxon = test_taxonomy.get_taxon(9606)
        assert taxon is not None
        assert taxon.tax_id == 9606
        assert taxon.scientific_name == "Homo sapiens"
        assert taxon.rank == "species"
        assert taxon.parent_tax_id == 9605

    def test_get_taxon_not_found(self, test_taxonomy):
        """Test lookup of non-existent taxid returns None."""
        taxon = test_taxonomy.get_taxon(999999999)
        assert taxon is None

    def test_get_taxon_resolves_merged(self, test_taxonomy):
        """Test that merged taxids resolve to current taxid."""
        taxon = test_taxonomy.get_taxon(12345)  # Merged -> 9606
        assert taxon is not None
        assert taxon.tax_id == 9606
        assert taxon.scientific_name == "Homo sapiens"

    def test_get_name(self, test_taxonomy):
        """Test get_name convenience method."""
        name = test_taxonomy.get_name(9606)
        assert name == "Homo sapiens"

    def test_get_name_not_found(self, test_taxonomy):
        """Test get_name returns None for non-existent taxid."""
        name = test_taxonomy.get_name(999999999)
        assert name is None

    def test_get_rank(self, test_taxonomy):
        """Test get_rank convenience method."""
        rank = test_taxonomy.get_rank(9606)
        assert rank == "species"

    def test_get_lineage(self, test_taxonomy):
        """Test lineage retrieval."""
        lineage = test_taxonomy.get_lineage(9606)
        assert len(lineage) > 0
        # Should go from root toward leaf (but not include root)
        assert lineage[-1].tax_id == 9606  # Leaf
        assert lineage[-1].scientific_name == "Homo sapiens"
        # Check some ancestors are present
        tax_ids = [t.tax_id for t in lineage]
        assert 9605 in tax_ids  # Homo
        assert 207598 in tax_ids  # Homininae
        assert 9604 in tax_ids  # Hominidae

    def test_get_lineage_caching(self, test_taxonomy):
        """Test that lineage results are cached."""
        # First call
        lineage1 = test_taxonomy.get_lineage(9606)
        # Second call should return cached result
        lineage2 = test_taxonomy.get_lineage(9606)
        # Should be the same object (cached)
        assert lineage1 is lineage2

    def test_get_lineage_ids(self, test_taxonomy):
        """Test get_lineage_ids returns list of tax_ids."""
        lineage_ids = test_taxonomy.get_lineage_ids(9606)
        assert isinstance(lineage_ids, list)
        assert all(isinstance(tid, int) for tid in lineage_ids)
        assert 9606 in lineage_ids  # Leaf
        assert 9605 in lineage_ids  # Homo
        assert 207598 in lineage_ids  # Homininae

    def test_get_lineage_string(self, test_taxonomy):
        """Test get_lineage_string formats correctly."""
        lineage_str = test_taxonomy.get_lineage_string(9606)
        assert isinstance(lineage_str, str)
        assert "species:Homo sapiens" in lineage_str
        assert "genus:Homo" in lineage_str
        assert "subfamily:Homininae" in lineage_str

    def test_get_lineage_string_excludes_unranked_by_default(self, test_taxonomy):
        """Test that unranked taxa are excluded by default."""
        lineage_str = test_taxonomy.get_lineage_string(9606)
        # "cellular organisms" has rank "no rank" in our test data
        # but we set it as "no rank" which should be excluded
        assert "no rank" not in lineage_str

    def test_get_many(self, test_taxonomy):
        """Test batch lookup."""
        result = test_taxonomy.get_many([9606, 9598])
        assert len(result) == 2
        assert 9606 in result
        assert 9598 in result
        assert result[9606].scientific_name == "Homo sapiens"
        assert result[9598].scientific_name == "Pan troglodytes"

    def test_get_many_empty_list(self, test_taxonomy):
        """Test batch lookup with empty list."""
        result = test_taxonomy.get_many([])
        assert result == {}


class TestGetManyMergedBehavior:
    """Tests documenting get_many() behavior with merged taxids."""

    def test_get_many_does_not_resolve_merged_taxids(self, test_taxonomy):
        """Document that get_many() does NOT resolve merged taxids (unlike get_taxon)."""
        # 12345 is merged to 9606 in our test data
        result = test_taxonomy.get_many([12345])
        # get_many queries taxons table directly, so merged taxids are not found
        assert result == {}

    def test_get_many_mixed_valid_and_merged_taxids(self, test_taxonomy):
        """Valid taxids are returned, merged ones are silently omitted."""
        result = test_taxonomy.get_many([9606, 12345, 9598])
        # Only valid taxids are returned
        assert set(result.keys()) == {9606, 9598}
        assert result[9606].scientific_name == "Homo sapiens"
        assert result[9598].scientific_name == "Pan troglodytes"

    def test_get_many_vs_get_taxon_merged_behavior_differs(self, test_taxonomy):
        """Explicitly document the behavioral difference between get_many and get_taxon."""
        # get_taxon resolves merged taxids
        taxon = test_taxonomy.get_taxon(12345)
        assert taxon is not None
        assert taxon.tax_id == 9606

        # get_many does NOT resolve merged taxids
        result = test_taxonomy.get_many([12345])
        assert result == {}


class TestLineageEdgeCases:
    """Tests for lineage methods with edge cases."""

    def test_get_lineage_invalid_taxid_returns_empty_list(self, test_taxonomy):
        """get_lineage returns empty list for non-existent taxid."""
        lineage = test_taxonomy.get_lineage(999999999)
        assert lineage == []
        assert isinstance(lineage, list)

    def test_get_lineage_ids_invalid_taxid_returns_empty_list(self, test_taxonomy):
        """get_lineage_ids returns empty list for non-existent taxid."""
        lineage_ids = test_taxonomy.get_lineage_ids(999999999)
        assert lineage_ids == []

    def test_get_lineage_string_invalid_taxid_returns_empty_string(self, test_taxonomy):
        """get_lineage_string returns empty string for non-existent taxid."""
        lineage_str = test_taxonomy.get_lineage_string(999999999)
        assert lineage_str == ""

    def test_get_lineage_caches_invalid_taxid_result(self, test_taxonomy):
        """Invalid taxid lookup is also cached to avoid repeated DB queries."""
        lineage1 = test_taxonomy.get_lineage(999999999)
        lineage2 = test_taxonomy.get_lineage(999999999)
        # Should be the same object (cached)
        assert lineage1 is lineage2


class TestFindLCA:
    """Tests for LCA (Lowest Common Ancestor) functionality."""

    def test_find_lca_same_species(self, test_taxonomy):
        """Test LCA of same species returns that species."""
        lca = test_taxonomy.find_lca([9606, 9606])
        assert lca == 9606

    def test_find_lca_human_chimp(self, test_taxonomy):
        """Test LCA of human and chimp is Homininae."""
        lca = test_taxonomy.find_lca([9606, 9598])
        assert lca == 207598  # Homininae

    def test_find_lca_single_taxid(self, test_taxonomy):
        """Test LCA of single taxid returns that taxid."""
        lca = test_taxonomy.find_lca([9606])
        assert lca == 9606

    def test_find_lca_empty_list(self, test_taxonomy):
        """Test LCA of empty list returns None."""
        lca = test_taxonomy.find_lca([])
        assert lca is None

    def test_find_lca_invalid_taxids_skipped(self, test_taxonomy):
        """Test that invalid taxids are skipped in LCA calculation."""
        # 999999999 doesn't exist, should be skipped
        lca = test_taxonomy.find_lca([9606, 999999999])
        # Should return 9606 since it's the only valid one
        assert lca == 9606


class TestMergedTaxidResolution:
    """Tests for merged taxid resolution including chains and edge cases."""

    def test_get_taxon_resolves_merge_chain(
        self, taxdump_with_merge_chain: Path, monkeypatch
    ):
        """Multi-level merge resolution: 11111 -> 22222 -> 9606."""
        taxonomy._build_sqlite_from_dmp(taxdump_with_merge_chain)
        monkeypatch.setattr(
            taxonomy, "_ensure_taxdump", lambda _=None: taxdump_with_merge_chain
        )

        with taxonomy.open() as tax:
            taxon = tax.get_taxon(11111)
            assert taxon is not None
            assert taxon.tax_id == 9606
            assert taxon.scientific_name == "Homo sapiens"

    def test_get_taxon_merge_to_nonexistent_returns_none(
        self, taxdump_with_merge_chain: Path, monkeypatch
    ):
        """Merge pointing to non-existent taxid returns None."""
        taxonomy._build_sqlite_from_dmp(taxdump_with_merge_chain)
        monkeypatch.setattr(
            taxonomy, "_ensure_taxdump", lambda _=None: taxdump_with_merge_chain
        )

        with taxonomy.open() as tax:
            # 99999 -> 88888, but 88888 doesn't exist
            taxon = tax.get_taxon(99999)
            assert taxon is None

    def test_get_taxon_cyclic_merge_returns_none(
        self, taxdump_with_cyclic_merge: Path, monkeypatch
    ):
        """Cyclic merge (A -> B -> A) returns None instead of infinite loop."""
        taxonomy._build_sqlite_from_dmp(taxdump_with_cyclic_merge)
        monkeypatch.setattr(
            taxonomy, "_ensure_taxdump", lambda _=None: taxdump_with_cyclic_merge
        )

        with taxonomy.open() as tax:
            # 11111 -> 22222 -> 11111 (cycle)
            # Should return None due to cycle detection, not infinite loop
            taxon = tax.get_taxon(11111)
            assert taxon is None


class TestTaxdumpManagement:
    """Tests for taxdump download and SQLite building."""

    def test_build_sqlite_from_dmp(self, minimal_taxdump: Path):
        """Test that SQLite is built correctly from .dmp files."""
        sqlite_path = taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        assert sqlite_path.exists()
        assert sqlite_path.name == "taxonomy.sqlite"

    def test_is_taxdump_stale_missing(self, tmp_path: Path):
        """Test that missing taxdump is considered stale."""
        taxdump_dir = tmp_path / "nonexistent"
        assert taxonomy._is_taxdump_stale(taxdump_dir) is True

    def test_is_taxdump_stale_fresh(self, minimal_taxdump: Path):
        """Test that fresh taxdump is not stale."""
        # The fixture just created the files, so they're fresh
        assert taxonomy._is_taxdump_stale(minimal_taxdump) is False

    def test_is_taxdump_stale_old(self, minimal_taxdump: Path):
        """Test that 31-day-old taxdump is considered stale."""
        import os
        import time

        nodes_dmp = minimal_taxdump / "nodes.dmp"
        # Set mtime to 31 days ago
        old_time = time.time() - (31 * 24 * 60 * 60)
        os.utime(nodes_dmp, (old_time, old_time))

        assert taxonomy._is_taxdump_stale(minimal_taxdump) is True

    def test_is_taxdump_stale_29_days_is_fresh(self, minimal_taxdump: Path):
        """Test that 29-day-old taxdump is NOT stale (boundary test)."""
        import os
        import time

        nodes_dmp = minimal_taxdump / "nodes.dmp"
        # Set mtime to 29 days ago
        old_time = time.time() - (29 * 24 * 60 * 60)
        os.utime(nodes_dmp, (old_time, old_time))

        assert taxonomy._is_taxdump_stale(minimal_taxdump) is False


class TestEnsureTaxdump:
    """Tests for _ensure_taxdump() control flow paths."""

    def test_ensure_taxdump_fresh_dmp_but_missing_sqlite_rebuilds(
        self, minimal_taxdump: Path, monkeypatch
    ):
        """Path 2: Fresh .dmp files but missing SQLite triggers rebuild only."""
        # Ensure SQLite doesn't exist
        sqlite_path = minimal_taxdump / "taxonomy.sqlite"
        if sqlite_path.exists():
            sqlite_path.unlink()

        # Track if download was called (it shouldn't be)
        download_called = []
        original_download = taxonomy._download_and_extract_taxdump

        def mock_download(taxdump_dir):
            download_called.append(True)
            original_download(taxdump_dir)

        monkeypatch.setattr(taxonomy, "_download_and_extract_taxdump", mock_download)
        monkeypatch.setattr(taxonomy, "get_taxdump_dir", lambda _=None: minimal_taxdump)

        # Call _ensure_taxdump
        result = taxonomy._ensure_taxdump()

        # Download should NOT have been called (files are fresh)
        assert download_called == []
        # SQLite should now exist
        assert sqlite_path.exists()
        assert result == minimal_taxdump

    def test_ensure_taxdump_fresh_with_sqlite_does_nothing(
        self, minimal_taxdump: Path, monkeypatch
    ):
        """Path 3: Fresh taxdump with SQLite returns immediately."""
        # Build SQLite first
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        sqlite_path = minimal_taxdump / "taxonomy.sqlite"
        original_mtime = sqlite_path.stat().st_mtime

        # Track if either function was called
        download_called = []
        build_called = []

        def mock_download(taxdump_dir):
            download_called.append(True)

        def mock_build(taxdump_dir):
            build_called.append(True)
            return taxdump_dir / "taxonomy.sqlite"

        monkeypatch.setattr(taxonomy, "_download_and_extract_taxdump", mock_download)
        monkeypatch.setattr(taxonomy, "_build_sqlite_from_dmp", mock_build)
        monkeypatch.setattr(taxonomy, "get_taxdump_dir", lambda _=None: minimal_taxdump)

        # Call _ensure_taxdump
        result = taxonomy._ensure_taxdump()

        # Neither should have been called
        assert download_called == []
        assert build_called == []
        # SQLite mtime should be unchanged
        assert sqlite_path.stat().st_mtime == original_mtime
        assert result == minimal_taxdump


class TestContextManagers:
    """Tests for context manager behavior."""

    def test_open_returns_taxonomy_db(self, minimal_taxdump: Path, monkeypatch):
        """Test that open() returns a TaxonomyDB instance."""
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        monkeypatch.setattr(taxonomy, "_ensure_taxdump", lambda _=None: minimal_taxdump)

        with taxonomy.open() as tax:
            assert isinstance(tax, taxonomy.TaxonomyDB)

    def test_open_connection_closed_after_context(
        self,
        minimal_taxdump: Path,
        monkeypatch,
    ):
        """Test that database connection is closed after context exits."""
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        monkeypatch.setattr(taxonomy, "_ensure_taxdump", lambda _=None: minimal_taxdump)

        with taxonomy.open() as tax:
            conn = tax._conn

        # After context, connection should be closed
        # Attempting to use it should raise an error
        with pytest.raises(sqlite3.ProgrammingError):
            conn.execute("SELECT 1")

    def test_wal_mode_is_enabled(self, minimal_taxdump: Path, monkeypatch):
        """Verify the database is opened in WAL mode for concurrent access."""
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        monkeypatch.setattr(taxonomy, "_ensure_taxdump", lambda _=None: minimal_taxdump)

        with taxonomy.open() as tax:
            result = tax._conn.execute("PRAGMA journal_mode;").fetchone()
            assert result[0] == "wal"


class TestConcurrentAccess:
    """Tests for concurrent database access."""

    def test_concurrent_readers_do_not_block(self, minimal_taxdump: Path, monkeypatch):
        """Multiple simultaneous readers can query without blocking."""
        import threading
        import time

        taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        monkeypatch.setattr(taxonomy, "_ensure_taxdump", lambda _=None: minimal_taxdump)

        results = []
        errors = []
        barrier = threading.Barrier(5)  # Wait for all threads to start

        def reader_thread(thread_id):
            try:
                barrier.wait(timeout=2)  # Synchronize start
                with taxonomy.open() as tax:
                    for _ in range(10):
                        taxon = tax.get_taxon(9606)
                        if taxon is None or taxon.scientific_name != "Homo sapiens":
                            errors.append(f"Thread {thread_id}: unexpected result")
                            return
                results.append(thread_id)
            except Exception as e:
                errors.append(f"Thread {thread_id}: {e}")

        threads = [threading.Thread(target=reader_thread, args=(i,)) for i in range(5)]
        start = time.time()

        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=5)

        elapsed = time.time() - start

        # All threads should complete successfully
        assert len(errors) == 0, f"Errors: {errors}"
        assert len(results) == 5
        # Should complete quickly (not blocked)
        assert elapsed < 3.0, f"Took too long: {elapsed}s"


@pytest.mark.slow
@pytest.mark.network
class TestNCBIIntegration:
    """
    Integration tests that download real data from NCBI.

    These tests are slow and require network access. They are skipped by default.
    Run with: pytest -m "slow and network"
    """

    def test_download_and_build_taxonomy(self, tmp_path):
        """
        Verify we can download and build taxonomy from NCBI.

        This test:
        1. Downloads taxdump.tar.gz from NCBI FTP (~50MB compressed)
        2. Extracts nodes.dmp, names.dmp, merged.dmp
        3. Builds taxonomy.sqlite
        4. Verifies known taxids resolve correctly
        """
        import os

        # Use a fresh temp directory to force download
        state_dir = tmp_path / "ncbi_test"
        state_dir.mkdir()
        os.environ["NVD_STATE_DIR"] = str(state_dir)

        try:
            with taxonomy.open() as tax:
                # Verify Homo sapiens exists and has correct data
                human = tax.get_taxon(9606)
                assert human is not None, "Homo sapiens (9606) not found"
                assert human.scientific_name == "Homo sapiens"
                assert human.rank == "species"

                # Verify Pan troglodytes (Chimpanzee)
                chimp = tax.get_taxon(9598)
                assert chimp is not None, "Pan troglodytes (9598) not found"
                assert chimp.scientific_name == "Pan troglodytes"

                # Verify lineage works
                lineage = tax.get_lineage(9606)
                assert len(lineage) > 5, "Lineage should have many ancestors"

                # Verify some expected ancestors are present
                lineage_names = [t.scientific_name for t in lineage]
                assert "Homo sapiens" in lineage_names
                assert "Hominidae" in lineage_names or "Primates" in lineage_names

                # Verify LCA works
                lca = tax.find_lca([9606, 9598])  # Human and Chimp
                assert lca is not None, "LCA calculation failed"
                # LCA should be Homininae (207598) or a parent thereof
                lca_taxon = tax.get_taxon(lca)
                assert lca_taxon is not None

        finally:
            if "NVD_STATE_DIR" in os.environ:
                del os.environ["NVD_STATE_DIR"]

    @pytest.mark.slow
    def test_taxonomy_sqlite_is_cached(self, tmp_path):
        """
        Verify that taxonomy.sqlite is cached and reused.

        Second open() should be fast because it uses cached database.
        """
        import os
        import time

        state_dir = tmp_path / "cache_test"
        state_dir.mkdir()
        os.environ["NVD_STATE_DIR"] = str(state_dir)

        try:
            # First open - downloads and builds (slow)
            start1 = time.time()
            with taxonomy.open() as tax:
                _ = tax.get_taxon(9606)
            elapsed1 = time.time() - start1

            # Second open - should use cache (fast)
            start2 = time.time()
            with taxonomy.open() as tax:
                _ = tax.get_taxon(9606)
            elapsed2 = time.time() - start2

            # Second should be much faster (at least 10x)
            assert elapsed2 < elapsed1 / 10, (
                f"Cache not working: first={elapsed1:.2f}s, second={elapsed2:.2f}s"
            )

        finally:
            if "NVD_STATE_DIR" in os.environ:
                del os.environ["NVD_STATE_DIR"]


class TestOfflineMode:
    """Tests for offline mode behavior."""

    def test_offline_with_existing_sqlite(self, minimal_taxdump):
        """Offline mode works when SQLite database exists."""
        # Build SQLite first
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)

        # Should work in offline mode
        with taxonomy.open(state_dir=minimal_taxdump.parent, offline=True) as tax:
            taxon = tax.get_taxon(9606)
            assert taxon is not None
            assert taxon.scientific_name == "Homo sapiens"

    def test_offline_builds_sqlite_from_dmp(self, minimal_taxdump):
        """Offline mode builds SQLite if .dmp files exist but SQLite doesn't."""
        # .dmp files exist but no SQLite yet
        sqlite_path = minimal_taxdump / "taxonomy.sqlite"
        assert not sqlite_path.exists()

        # Should build SQLite from .dmp files
        with taxonomy.open(state_dir=minimal_taxdump.parent, offline=True) as tax:
            taxon = tax.get_taxon(9606)
            assert taxon is not None

        # SQLite should now exist
        assert sqlite_path.exists()

    def test_offline_raises_when_no_data(self, tmp_path):
        """Offline mode raises TaxonomyOfflineError when no data exists."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        with pytest.raises(taxonomy.TaxonomyOfflineError) as exc_info:
            with taxonomy.open(state_dir=empty_dir, offline=True):
                pass

        assert "offline mode is enabled" in str(exc_info.value)
        assert str(empty_dir) in str(exc_info.value)

    def test_offline_env_var(self, minimal_taxdump, monkeypatch):
        """NVD_TAXONOMY_OFFLINE env var enables offline mode."""
        # Build SQLite first
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)

        # Set env var
        monkeypatch.setenv("NVD_TAXONOMY_OFFLINE", "1")

        # Should use offline mode (won't try to download)
        with taxonomy.open(state_dir=minimal_taxdump.parent) as tax:
            taxon = tax.get_taxon(9606)
            assert taxon is not None

    def test_offline_env_var_true(self, minimal_taxdump, monkeypatch):
        """NVD_TAXONOMY_OFFLINE=true also works."""
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        monkeypatch.setenv("NVD_TAXONOMY_OFFLINE", "true")

        with taxonomy.open(state_dir=minimal_taxdump.parent) as tax:
            taxon = tax.get_taxon(9606)
            assert taxon is not None

    def test_explicit_offline_overrides_env_var(self, minimal_taxdump, monkeypatch):
        """Explicit offline=False overrides env var."""
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)

        # Env var says offline, but explicit param says online
        monkeypatch.setenv("NVD_TAXONOMY_OFFLINE", "1")

        # offline=False should override env var
        # (This won't try to download because data isn't stale)
        with taxonomy.open(state_dir=minimal_taxdump.parent, offline=False) as tax:
            taxon = tax.get_taxon(9606)
            assert taxon is not None


class TestTaxonomyDbEnvVar:
    """Tests for NVD_TAXONOMY_DB environment variable override."""

    def test_env_var_overrides_default(self, minimal_taxdump, monkeypatch):
        """NVD_TAXONOMY_DB points to custom taxonomy location."""
        # Build SQLite in the minimal_taxdump directory
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)

        # Set env var to point directly to the taxdump directory
        monkeypatch.setenv("NVD_TAXONOMY_DB", str(minimal_taxdump))

        # Should use the env var path, not the default
        resolved = get_taxdump_dir()
        assert resolved == minimal_taxdump

    def test_env_var_works_with_taxonomy_open(self, minimal_taxdump, monkeypatch):
        """taxonomy.open() respects NVD_TAXONOMY_DB."""
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)
        monkeypatch.setenv("NVD_TAXONOMY_DB", str(minimal_taxdump))

        # Should find taxonomy at the env var location
        with taxonomy.open(offline=True) as tax:
            taxon = tax.get_taxon(9606)
            assert taxon is not None
            assert taxon.scientific_name == "Homo sapiens"

    def test_env_var_with_shared_readonly_taxonomy(
        self, minimal_taxdump, tmp_path, monkeypatch
    ):
        """Simulates cluster setup: shared taxonomy, separate state dir."""
        # Build taxonomy in "shared" location
        taxonomy._build_sqlite_from_dmp(minimal_taxdump)

        # Point to shared taxonomy
        monkeypatch.setenv("NVD_TAXONOMY_DB", str(minimal_taxdump))

        # Use separate state dir (simulating user's home on cluster)
        user_state = tmp_path / "user_state"
        user_state.mkdir()
        monkeypatch.setenv("NVD_STATE_DIR", str(user_state))

        # Should use shared taxonomy, not look in user_state
        with taxonomy.open(offline=True) as tax:
            taxon = tax.get_taxon(9606)
            assert taxon is not None

        # Verify taxonomy was NOT copied to user state
        assert not (user_state / "taxdump" / "taxonomy.sqlite").exists()
