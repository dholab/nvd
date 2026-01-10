"""
Unified taxonomy access layer.

All taxonomy lookups go through this module. It manages a local cache of NCBI
taxonomy data downloaded from taxdump.tar.gz and builds a SQLite database for
fast lookups.

Features:
    - Automatic download and extraction of NCBI taxdump
    - SQLite database built from .dmp files for fast queries
    - Merged taxid resolution (deprecated taxids map to current ones)
    - Lineage caching for performance during batch operations
    - taxopy integration for LCA (Least Common Ancestor) calculations

State directory resolution follows the same hierarchy as db.py:
    1. Explicit path argument
    2. NVD_STATE_DIR environment variable
    3. Default: ~/.nvd/
"""

from __future__ import annotations

import builtins
import os
import sqlite3
import tarfile
import urllib.request
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from typing import TYPE_CHECKING

import taxopy
import taxopy.utilities
from taxopy.exceptions import TaxidError

from py_nvd.db import get_taxdump_dir
from py_nvd.models import Taxon

# Save builtin open before we shadow it with our context manager
_open = builtins.open

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

# NCBI taxdump URL and configuration
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
REQUIRED_DMP_FILES = {"nodes.dmp", "names.dmp", "merged.dmp"}
MAX_AGE = timedelta(days=30)


def _is_taxdump_stale(taxdump_dir: Path) -> bool:
    """Check if taxdump needs refresh (missing or >30 days old)."""
    nodes_dmp = taxdump_dir / "nodes.dmp"
    if not nodes_dmp.exists():
        return True
    mtime = datetime.fromtimestamp(nodes_dmp.stat().st_mtime, tz=timezone.utc)
    age = datetime.now(tz=timezone.utc) - mtime
    return age > MAX_AGE


def _download_and_extract_taxdump(taxdump_dir: Path) -> None:
    """
    Download taxdump.tar.gz and extract required files.

    Only extracts nodes.dmp, names.dmp, and merged.dmp (~500MB total).
    Uses extractfile() for safety against path traversal attacks.
    """
    taxdump_dir.mkdir(parents=True, exist_ok=True)
    tar_path = taxdump_dir / "taxdump.tar.gz"

    # Download
    urllib.request.urlretrieve(TAXDUMP_URL, tar_path)  # noqa: S310

    # Extract only needed files using extractfile() for safety
    # This avoids the deprecated tar.extract() behavior in Python 3.12+
    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            if member.name in REQUIRED_DMP_FILES:
                file_obj = tar.extractfile(member)
                if file_obj is not None:
                    dest_path = taxdump_dir / member.name
                    with _open(dest_path, "wb") as dest:
                        dest.write(file_obj.read())

    # Clean up tarball
    tar_path.unlink()


def _build_sqlite_from_dmp(taxdump_dir: Path) -> Path:
    """
    Build taxonomy.sqlite from .dmp files.

    Creates a SQLite database with the same schema as gettax.sqlite,
    plus a 'merged' table for deprecated taxid resolution.

    Returns path to the built database.
    """
    sqlite_path = taxdump_dir / "taxonomy.sqlite"

    # Remove existing to rebuild
    if sqlite_path.exists():
        sqlite_path.unlink()

    conn = sqlite3.connect(sqlite_path)
    # Use DELETE journal mode for network filesystem compatibility.
    # WAL mode requires shared memory (mmap) which doesn't work on NFS/GPFS/Lustre.
    conn.execute("PRAGMA journal_mode=DELETE")

    # Create schema (matches gettax.sqlite + merged table)
    conn.executescript("""
        CREATE TABLE taxons (
            tax_id INTEGER PRIMARY KEY,
            parent_tax_id INTEGER,
            rank TEXT,
            scientific_name TEXT
        );

        CREATE TABLE merged (
            old_tax_id INTEGER PRIMARY KEY,
            new_tax_id INTEGER
        );

        CREATE INDEX idx_taxons_parent ON taxons(parent_tax_id);
    """)

    # Parse nodes.dmp: tax_id | parent_tax_id | rank | ...
    # Format: tax_id<tab>|<tab>parent_tax_id<tab>|<tab>rank<tab>|<tab>...
    nodes: dict[int, tuple[int, str]] = {}
    with _open(taxdump_dir / "nodes.dmp") as f:
        for line in f:
            parts = line.split("\t|\t")
            tax_id = int(parts[0])
            parent_tax_id = int(parts[1])
            rank = parts[2].strip()
            nodes[tax_id] = (parent_tax_id, rank)

    # Parse names.dmp: tax_id | name | unique_name | name_class
    # We only want "scientific name" entries
    names: dict[int, str] = {}
    with _open(taxdump_dir / "names.dmp") as f:
        for line in f:
            parts = line.split("\t|\t")
            tax_id = int(parts[0])
            name = parts[1]
            name_class = parts[3].rstrip("\t|\n")
            if name_class == "scientific name":
                names[tax_id] = name

    # Insert taxons
    taxon_rows = [
        (tax_id, parent_tax_id, rank, names.get(tax_id, ""))
        for tax_id, (parent_tax_id, rank) in nodes.items()
    ]
    conn.executemany(
        "INSERT INTO taxons (tax_id, parent_tax_id, rank, scientific_name) VALUES (?, ?, ?, ?)",
        taxon_rows,
    )

    # Parse merged.dmp: old_tax_id | new_tax_id
    merged_rows = []
    with _open(taxdump_dir / "merged.dmp") as f:
        for line in f:
            parts = line.split("\t|\t")
            old_tax_id = int(parts[0])
            new_tax_id = int(parts[1].rstrip("\t|\n"))
            merged_rows.append((old_tax_id, new_tax_id))

    conn.executemany(
        "INSERT INTO merged (old_tax_id, new_tax_id) VALUES (?, ?)",
        merged_rows,
    )

    conn.commit()
    conn.close()

    return sqlite_path


def _ensure_taxdump(state_dir: Path | str | None = None) -> Path:
    """
    Ensure taxdump is downloaded, extracted, and SQLite is built.

    Returns path to taxdump directory.
    """
    taxdump_dir = get_taxdump_dir(state_dir)
    sqlite_path = taxdump_dir / "taxonomy.sqlite"

    if _is_taxdump_stale(taxdump_dir):
        _download_and_extract_taxdump(taxdump_dir)
        _build_sqlite_from_dmp(taxdump_dir)
    elif not sqlite_path.exists():
        # .dmp files exist but SQLite needs rebuilding
        _build_sqlite_from_dmp(taxdump_dir)

    return taxdump_dir


def ensure_taxonomy_available(state_dir: Path | str | None = None) -> Path:
    """
    Ensure taxonomy database is downloaded and ready for use.

    Downloads NCBI taxdump if not present or stale (>30 days old),
    extracts required .dmp files, and builds the SQLite database.

    This is the public API for pre-downloading taxonomy data before
    distributed workers need it. Call this once on the submit node
    (e.g., in CHECK_RUN_STATE) to avoid multiple workers racing to
    download simultaneously.

    Args:
        state_dir: Optional state directory override. If None, uses
                   NVD_STATE_DIR env var or ~/.nvd/

    Returns:
        Path to the taxdump directory containing taxonomy.sqlite
    """
    return _ensure_taxdump(state_dir)


class TaxonomyDB:
    """
    Helper class for taxonomy lookups.

    Provides methods for taxon lookup, lineage traversal, and LCA calculation.
    Automatically resolves merged (deprecated) taxids to their current equivalents.
    """

    def __init__(self, conn: sqlite3.Connection, taxdump_dir: Path) -> None:
        self._conn = conn
        self._taxdump_dir = taxdump_dir
        self._lineage_cache: dict[int, list[Taxon]] = {}
        self._taxopy_db: taxopy.TaxDb | None = None

    def get_taxon(self, tax_id: int, _seen: set[int] | None = None) -> Taxon | None:
        """
        Get full taxon info for a tax ID.

        Automatically resolves merged (deprecated) taxids to their
        current equivalents using the merged.dmp data. Handles multi-level
        merges (A → B → C) and detects cycles to prevent infinite recursion.

        Args:
            tax_id: The taxonomy ID to look up
            _seen: Internal parameter for cycle detection (do not use)

        Returns:
            Taxon if found (possibly after resolving merges), None otherwise
        """
        # Cycle detection for corrupt merged.dmp data
        if _seen is None:
            _seen = set()
        if tax_id in _seen:
            return None  # Cycle detected, bail out
        _seen.add(tax_id)

        row = self._conn.execute(
            "SELECT tax_id, scientific_name, rank, parent_tax_id FROM taxons WHERE tax_id = ?",
            (tax_id,),
        ).fetchone()

        if row is None:
            # Check if this is a merged (deprecated) taxid
            merged = self._conn.execute(
                "SELECT new_tax_id FROM merged WHERE old_tax_id = ?",
                (tax_id,),
            ).fetchone()
            if merged:
                return self.get_taxon(merged["new_tax_id"], _seen)  # Recursive lookup
            return None

        return Taxon(
            tax_id=row["tax_id"],
            scientific_name=row["scientific_name"],
            rank=row["rank"],
            parent_tax_id=row["parent_tax_id"],
        )

    def get_name(self, tax_id: int) -> str | None:
        """Get the scientific name for a tax ID."""
        taxon = self.get_taxon(tax_id)
        return taxon.scientific_name if taxon else None

    def get_rank(self, tax_id: int) -> str | None:
        """Get the taxonomic rank for a tax ID."""
        taxon = self.get_taxon(tax_id)
        return taxon.rank if taxon else None

    def get_lineage(self, tax_id: int) -> list[Taxon]:
        """
        Get full lineage for a tax ID (from root to taxon).

        Results are cached for performance during batch operations.
        """
        if tax_id in self._lineage_cache:
            return self._lineage_cache[tax_id]

        lineage = []
        current_id: int | None = tax_id

        while current_id and current_id != 1:  # 1 is root
            taxon = self.get_taxon(current_id)
            if not taxon:
                break
            lineage.append(taxon)
            current_id = taxon.parent_tax_id

        lineage.reverse()
        self._lineage_cache[tax_id] = lineage
        return lineage

    def get_lineage_ids(self, tax_id: int) -> list[int]:
        """
        Get lineage as list of tax IDs from root to leaf.

        Used by hits_to_report.py for LCA calculation.

        Example:
            >>> tax.get_lineage_ids(9606)  # Homo sapiens
            [131567, 2759, ..., 9605, 9606]
        """
        return [t.tax_id for t in self.get_lineage(tax_id)]

    def get_lineage_string(self, tax_id: int, include_unranked: bool = False) -> str:
        """
        Get lineage as formatted string "rank:name; rank:name; ...".

        Used by annotate_blast_results.py for TSV output.

        Args:
            tax_id: The taxonomy ID to look up
            include_unranked: If False (default), skip taxa with empty/no rank

        Example:
            >>> tax.get_lineage_string(9606)
            "superkingdom:Eukaryota; clade:Opisthokonta; kingdom:Metazoa; ..."
        """
        lineage = self.get_lineage(tax_id)
        parts = []
        for taxon in lineage:
            # Skip unranked taxa (rank is None, empty, or "no rank") unless requested
            has_rank = taxon.rank and taxon.rank != "no rank"
            if has_rank or include_unranked:
                rank = taxon.rank or "no rank"
                parts.append(f"{rank}:{taxon.scientific_name}")
        return "; ".join(parts)

    def get_many(self, tax_ids: list[int]) -> dict[int, Taxon]:
        """
        Batch lookup for multiple tax IDs.

        Returns dict mapping tax_id to Taxon. Tax IDs not found in the database
        are silently omitted from the result.

        Note:
            Unlike get_taxon(), this method does NOT resolve merged (deprecated)
            taxids. If you pass a merged taxid, it will be missing from the
            result. Use get_taxon() for individual lookups that need merged
            taxid resolution.
        """
        if not tax_ids:
            return {}

        placeholders = ",".join("?" * len(tax_ids))
        rows = self._conn.execute(
            f"SELECT tax_id, scientific_name, rank, parent_tax_id "  # noqa: S608
            f"FROM taxons WHERE tax_id IN ({placeholders})",
            tax_ids,
        ).fetchall()

        return {
            row["tax_id"]: Taxon(
                tax_id=row["tax_id"],
                scientific_name=row["scientific_name"],
                rank=row["rank"],
                parent_tax_id=row["parent_tax_id"],
            )
            for row in rows
        }

    @property
    def taxopy_db(self) -> taxopy.TaxDb:
        """
        Lazy-load taxopy TaxDb for LCA and advanced operations.

        First access takes 2-3 seconds to load ~100MB into memory.
        Subsequent accesses return the cached instance.

        Example:
            >>> tx = tax.taxopy_db
            >>> taxon = taxopy.Taxon(9606, tx)
            >>> taxon.name
            'Homo sapiens'
        """
        if self._taxopy_db is None:
            self._taxopy_db = taxopy.TaxDb(
                nodes_dmp=str(self._taxdump_dir / "nodes.dmp"),
                names_dmp=str(self._taxdump_dir / "names.dmp"),
                merged_dmp=str(self._taxdump_dir / "merged.dmp"),
            )
        return self._taxopy_db

    def find_lca(self, tax_ids: list[int]) -> int | None:
        """
        Find the Lowest Common Ancestor of multiple tax IDs.

        Uses taxopy's battle-tested LCA algorithm. Invalid or merged
        taxids are silently skipped.

        Args:
            tax_ids: List of taxonomy IDs to find LCA for

        Returns:
            The tax_id of the LCA, or None if no valid taxids provided

        Example:
            >>> tax.find_lca([9606, 9598])  # Human and Chimp
            207598  # Homininae
        """

        if not tax_ids:
            return None
        if len(tax_ids) == 1:
            return tax_ids[0]

        valid_taxons = []
        for tid in tax_ids:
            try:
                valid_taxons.append(taxopy.Taxon(tid, self.taxopy_db))
            except (TaxidError, KeyError, ValueError):
                # Skip invalid taxids. Note: merged taxids should be auto-resolved
                # by taxopy since we initialize TaxDb with merged_dmp.
                continue

        if not valid_taxons:
            return None
        if len(valid_taxons) == 1:
            return valid_taxons[0].taxid

        return taxopy.utilities.find_lca(valid_taxons, self.taxopy_db).taxid


class TaxonomyOfflineError(Exception):
    """Raised when taxonomy database is unavailable in offline mode."""


@contextmanager
def open(  # noqa: A001
    state_dir: Path | str | None = None,
    offline: bool | None = None,
) -> Generator[TaxonomyDB, None, None]:
    """
    Context manager for taxonomy database access.

    Downloads taxdump from NCBI if missing or stale (>30 days).
    Builds SQLite from .dmp files for fast lookups.
    Provides access to taxopy for LCA calculations.

    Args:
        state_dir: Optional explicit state directory. If None, resolves
                   via NVD_STATE_DIR env var, then ~/.nvd/
        offline: If True, never attempt to download taxonomy data.
                 If None (default), checks NVD_TAXONOMY_OFFLINE env var.
                 Set to True for air-gapped environments with pre-cached data.

    Raises:
        TaxonomyOfflineError: If offline=True and no cached taxonomy exists.

    Usage:
        with taxonomy.open() as tax:
            taxon = tax.get_taxon(9606)
            print(taxon.scientific_name)  # "Homo sapiens"

            # Convenience methods for script compatibility
            lineage_str = tax.get_lineage_string(9606)
            lineage_ids = tax.get_lineage_ids(9606)

            # LCA calculation (uses taxopy internally)
            lca = tax.find_lca([9606, 9598])  # Human and Chimp -> Homininae

        # Offline mode (for air-gapped environments)
        with taxonomy.open(offline=True) as tax:
            # Uses cached data only, never downloads
            taxon = tax.get_taxon(9606)
    """
    # Resolve offline mode from env var if not explicitly set
    if offline is None:
        offline = os.environ.get("NVD_TAXONOMY_OFFLINE", "").lower() in ("1", "true")

    taxdump_dir = get_taxdump_dir(state_dir)
    sqlite_path = taxdump_dir / "taxonomy.sqlite"

    if offline and not sqlite_path.exists():
        # Offline mode with no SQLite - try to build from .dmp files
        nodes_dmp = taxdump_dir / "nodes.dmp"
        if nodes_dmp.exists():
            _build_sqlite_from_dmp(taxdump_dir)
        else:
            raise TaxonomyOfflineError(
                f"Taxonomy database not found at {sqlite_path} and offline mode is enabled. "
                f"Either disable offline mode to download from NCBI, or copy taxonomy data to {taxdump_dir}"
            )
    elif not offline:
        # Normal mode: download if stale
        taxdump_dir = _ensure_taxdump(state_dir)
        sqlite_path = taxdump_dir / "taxonomy.sqlite"

    conn = sqlite3.connect(sqlite_path, timeout=30.0)
    conn.row_factory = sqlite3.Row

    try:
        yield TaxonomyDB(conn, taxdump_dir)
    finally:
        conn.close()
