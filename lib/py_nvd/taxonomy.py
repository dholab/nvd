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

Taxonomy cache resolution is intentionally small:
    1. Explicit taxonomy_dir argument
    2. NVD_TAXONOMY_DB environment variable
    3. NVD_CONFIG_DIR/taxdump, defaulting to ~/.config/nvd/taxdump
"""

from __future__ import annotations

import builtins
import os
import sqlite3
import sys
import tarfile
import urllib.error
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

import taxopy
import taxopy.utilities
from taxopy.exceptions import TaxidError

from py_nvd.models import Taxon
from py_nvd.paths import get_taxdump_dir

# Save builtin open before we shadow it with our context manager
_open = builtins.open

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

# NCBI taxdump URL and configuration
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
REQUIRED_DMP_FILES = {"nodes.dmp", "names.dmp", "merged.dmp"}
DEFAULT_TAXONOMY_MAX_AGE_DAYS = 90
MAX_AGE = timedelta(days=DEFAULT_TAXONOMY_MAX_AGE_DAYS)


class TaxonomyMode(StrEnum):
    """Pipeline-facing taxonomy availability modes."""

    READ_ONLY = "read_only"
    MISSING = "missing"


class TaxonomyRefresh(StrEnum):
    """Admin-facing taxonomy preparation refresh modes."""

    MISSING = "missing"
    STALE = "stale"
    FORCE = "force"


class TaxonomyAction(StrEnum):
    """Concrete action selected by taxonomy planning."""

    USE_EXISTING = "use_existing"
    BUILD_SQLITE = "build_sqlite"
    DOWNLOAD_AND_BUILD = "download_and_build"
    FAIL = "fail"


@dataclass(frozen=True)
class TaxonomyStatus:
    """Read-only inspection of a taxonomy directory."""

    taxdump_dir: Path
    missing_dmp_files: tuple[str, ...]
    sqlite_exists: bool
    age_days: int | None

    @property
    def has_sources(self) -> bool:
        """Whether all required NCBI taxdump source files are present."""
        return not self.missing_dmp_files

    @property
    def is_prepared(self) -> bool:
        """Whether taxonomy can be opened without mutating the directory."""
        return self.has_sources and self.sqlite_exists

    def is_stale(self, max_age_days: int = DEFAULT_TAXONOMY_MAX_AGE_DAYS) -> bool:
        """Whether taxonomy is older than the configured freshness window."""
        return self.age_days is None or self.age_days > max_age_days


@dataclass(frozen=True)
class TaxonomyPlan:
    """Decision about how taxonomy should be made available."""

    status: TaxonomyStatus
    action: TaxonomyAction
    reason: str


@dataclass(frozen=True)
class TaxonomyPolicy:
    """Resolved pipeline taxonomy availability policy."""

    mode: TaxonomyMode = TaxonomyMode.MISSING
    max_age_days: int = DEFAULT_TAXONOMY_MAX_AGE_DAYS

    @classmethod
    def from_env_and_params(
        cls,
        *,
        mode: str | TaxonomyMode | None,
        max_age_days: int | None,
        offline: bool,
    ) -> TaxonomyPolicy:
        """Resolve explicit params plus legacy offline mode into one policy."""
        if mode is None and offline:
            resolved_mode = TaxonomyMode.READ_ONLY
        elif mode is None:
            resolved_mode = TaxonomyMode.MISSING
        else:
            resolved_mode = TaxonomyMode(mode)
        return cls(
            mode=resolved_mode,
            max_age_days=max_age_days or DEFAULT_TAXONOMY_MAX_AGE_DAYS,
        )


def _missing_required_dmp_files(taxdump_dir: Path) -> list[Path]:
    """Return required NCBI taxdump files absent from a taxonomy directory."""
    return sorted(
        (
            taxdump_dir / filename
            for filename in REQUIRED_DMP_FILES
            if not (taxdump_dir / filename).exists()
        ),
        key=lambda path: path.name,
    )


def _required_dmp_files_missing(taxdump_dir: Path) -> bool:
    """Whether any required NCBI taxdump source files are missing."""
    return bool(_missing_required_dmp_files(taxdump_dir))


def _is_taxdump_stale(taxdump_dir: Path, max_age: timedelta = MAX_AGE) -> bool:
    """Check if taxdump is missing or older than the configured freshness window."""
    nodes_dmp = taxdump_dir / "nodes.dmp"
    if not nodes_dmp.exists():
        return True
    mtime = datetime.fromtimestamp(nodes_dmp.stat().st_mtime, tz=UTC)
    age = datetime.now(tz=UTC) - mtime
    return age > max_age


def inspect_taxonomy(taxonomy_dir: Path | str | None = None) -> TaxonomyStatus:
    """Inspect taxonomy availability without creating, downloading, or rebuilding."""
    taxdump_dir = get_taxdump_dir(taxonomy_dir=taxonomy_dir)
    nodes_dmp = taxdump_dir / "nodes.dmp"
    age_days = None
    if nodes_dmp.exists():
        mtime = datetime.fromtimestamp(nodes_dmp.stat().st_mtime, tz=UTC)
        age_days = (datetime.now(tz=UTC) - mtime).days
    return TaxonomyStatus(
        taxdump_dir=taxdump_dir,
        missing_dmp_files=tuple(
            path.name for path in _missing_required_dmp_files(taxdump_dir)
        ),
        sqlite_exists=(taxdump_dir / "taxonomy.sqlite").exists(),
        age_days=age_days,
    )


def plan_pipeline_taxonomy(
    status: TaxonomyStatus,
    policy: TaxonomyPolicy,
) -> TaxonomyPlan:
    """Plan taxonomy availability for pipeline runs."""
    if policy.mode is TaxonomyMode.READ_ONLY:
        if status.is_prepared:
            return TaxonomyPlan(
                status,
                TaxonomyAction.USE_EXISTING,
                "prepared taxonomy exists",
            )
        missing = [*status.missing_dmp_files]
        if not status.sqlite_exists:
            missing.append("taxonomy.sqlite")
        return TaxonomyPlan(
            status,
            TaxonomyAction.FAIL,
            "read-only taxonomy mode requires prepared taxonomy at "
            f"{status.taxdump_dir}; missing: {', '.join(sorted(missing))}",
        )

    if status.is_prepared:
        return TaxonomyPlan(
            status,
            TaxonomyAction.USE_EXISTING,
            "prepared taxonomy exists",
        )
    if status.has_sources:
        return TaxonomyPlan(
            status,
            TaxonomyAction.BUILD_SQLITE,
            "taxonomy.sqlite is missing",
        )
    return TaxonomyPlan(
        status,
        TaxonomyAction.DOWNLOAD_AND_BUILD,
        "required taxdump files are missing",
    )


def plan_admin_taxonomy(
    status: TaxonomyStatus,
    *,
    refresh: TaxonomyRefresh = TaxonomyRefresh.MISSING,
    max_age_days: int = DEFAULT_TAXONOMY_MAX_AGE_DAYS,
) -> TaxonomyPlan:
    """Plan taxonomy preparation for explicit admin commands."""
    if refresh is TaxonomyRefresh.FORCE:
        return TaxonomyPlan(
            status,
            TaxonomyAction.DOWNLOAD_AND_BUILD,
            "force refresh requested",
        )
    if not status.has_sources:
        return TaxonomyPlan(
            status,
            TaxonomyAction.DOWNLOAD_AND_BUILD,
            "required taxdump files are missing",
        )
    if refresh is TaxonomyRefresh.STALE and status.is_stale(max_age_days):
        return TaxonomyPlan(
            status,
            TaxonomyAction.DOWNLOAD_AND_BUILD,
            "taxonomy is stale",
        )
    if not status.sqlite_exists:
        return TaxonomyPlan(
            status,
            TaxonomyAction.BUILD_SQLITE,
            "taxonomy.sqlite is missing",
        )
    return TaxonomyPlan(status, TaxonomyAction.USE_EXISTING, "prepared taxonomy exists")


def execute_taxonomy_plan(plan: TaxonomyPlan) -> Path:
    """Execute a taxonomy plan, mutating only for build/download actions."""
    taxdump_dir = plan.status.taxdump_dir
    if plan.action is TaxonomyAction.USE_EXISTING:
        return taxdump_dir
    if plan.action is TaxonomyAction.BUILD_SQLITE:
        _build_sqlite_from_dmp(taxdump_dir)
        return taxdump_dir
    if plan.action is TaxonomyAction.DOWNLOAD_AND_BUILD:
        _download_and_extract_taxdump(taxdump_dir)
        _build_sqlite_from_dmp(taxdump_dir)
        return taxdump_dir
    raise TaxonomyOfflineError(plan.reason)


def _download_and_extract_taxdump(taxdump_dir: Path) -> None:
    """
    Download taxdump.tar.gz and extract required files.

    Only extracts nodes.dmp, names.dmp, and merged.dmp (~500MB total).
    Uses extractfile() for safety against path traversal attacks.

    Extraction is atomic: files are written to .tmp temporaries and only
    renamed to their final names after all three are fully written. An
    interrupted extraction (job kill, disk full, network hiccup) cannot
    corrupt or truncate existing .dmp files.
    """
    taxdump_dir.mkdir(parents=True, exist_ok=True)
    tar_path = taxdump_dir / "taxdump.tar.gz"

    # Download
    urllib.request.urlretrieve(TAXDUMP_URL, tar_path)  # noqa: S310

    # Extract to temporary files, then rename atomically. This prevents an
    # interrupted extraction from leaving truncated .dmp files behind (the
    # "wb" open mode truncates the target to zero bytes before writing).
    tmp_paths: list[tuple[Path, Path]] = []
    rename_ok = False
    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            for member in tar.getmembers():
                if member.name not in REQUIRED_DMP_FILES:
                    continue
                file_obj = tar.extractfile(member)
                if file_obj is None:
                    continue
                tmp = taxdump_dir / f"{member.name}.tmp"
                with _open(tmp, "wb") as dest:
                    dest.write(file_obj.read())
                tmp_paths.append((tmp, taxdump_dir / member.name))

        for tmp, final in tmp_paths:
            os.replace(tmp, final)
        rename_ok = True
    finally:
        if not rename_ok:
            for tmp, _final in tmp_paths:
                tmp.unlink(missing_ok=True)

    # Clean up tarball
    tar_path.unlink()


def _build_sqlite_from_dmp(taxdump_dir: Path) -> Path:  # noqa: C901, PLR0912, PLR0915
    """
    Build taxonomy.sqlite from .dmp files.

    Creates a SQLite database with the same schema as gettax.sqlite,
    plus a 'merged' table for deprecated taxid resolution.

    The build is atomic: writes to a temporary file and renames on success,
    so a crash mid-build never leaves a corrupt taxonomy.sqlite on disk.
    All .dmp files are parsed into memory before the database is created,
    so parse errors cannot produce a partial database.

    Returns path to the built database.

    Raises:
        TaxonomyBuildError: If .dmp files cannot be parsed or the resulting
            database fails validation.
    """
    sqlite_path = taxdump_dir / "taxonomy.sqlite"
    tmp_path = taxdump_dir / "taxonomy.sqlite.tmp"

    # Clean up any stale temp file from a previous failed build
    if tmp_path.exists():
        tmp_path.unlink()

    # Parse all .dmp files into memory BEFORE creating the database.
    # This ensures parse errors never leave a partial taxonomy.sqlite on disk.

    # Parse nodes.dmp: tax_id | parent_tax_id | rank | ...
    nodes: dict[int, tuple[int, str]] = {}
    nodes_path = taxdump_dir / "nodes.dmp"
    with _open(nodes_path) as f:
        for line_num, line in enumerate(f, start=1):
            try:
                parts = line.split("\t|\t")
                tax_id = int(parts[0])
                parent_tax_id = int(parts[1])
                rank = parts[2].strip()
            except (ValueError, IndexError) as e:
                msg = f"Failed to parse nodes.dmp line {line_num}: {e}\n  Line content: {line.rstrip()!r}"
                raise TaxonomyBuildError(msg) from e
            nodes[tax_id] = (parent_tax_id, rank)

    if not nodes:
        msg = f"nodes.dmp is empty: {nodes_path}"
        raise TaxonomyBuildError(msg)

    # Parse names.dmp: tax_id | name | unique_name | name_class
    names: dict[int, str] = {}
    names_path = taxdump_dir / "names.dmp"
    with _open(names_path) as f:
        for line_num, line in enumerate(f, start=1):
            try:
                parts = line.split("\t|\t")
                tax_id = int(parts[0])
                name = parts[1]
                name_class = parts[3].rstrip("\t|\n")
            except (ValueError, IndexError) as e:
                msg = f"Failed to parse names.dmp line {line_num}: {e}\n  Line content: {line.rstrip()!r}"
                raise TaxonomyBuildError(msg) from e
            if name_class == "scientific name":
                names[tax_id] = name

    if not names:
        msg = f"names.dmp contains no scientific name entries: {names_path}"
        raise TaxonomyBuildError(msg)

    # Parse merged.dmp: old_tax_id | new_tax_id
    merged_rows: list[tuple[int, int]] = []
    merged_path = taxdump_dir / "merged.dmp"
    with _open(merged_path) as f:
        for line_num, line in enumerate(f, start=1):
            try:
                parts = line.split("\t|\t")
                old_tax_id = int(parts[0])
                new_tax_id = int(parts[1].rstrip("\t|\n"))
            except (ValueError, IndexError) as e:
                msg = f"Failed to parse merged.dmp line {line_num}: {e}\n  Line content: {line.rstrip()!r}"
                raise TaxonomyBuildError(msg) from e
            merged_rows.append((old_tax_id, new_tax_id))

    # Build taxon rows, joining nodes with names
    taxon_rows = [
        (tax_id, parent_tax_id, rank, names.get(tax_id, ""))
        for tax_id, (parent_tax_id, rank) in nodes.items()
    ]

    conn = sqlite3.connect(tmp_path)
    try:
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

        conn.executemany(
            "INSERT INTO taxons (tax_id, parent_tax_id, rank, scientific_name) VALUES (?, ?, ?, ?)",
            taxon_rows,
        )
        conn.executemany(
            "INSERT INTO merged (old_tax_id, new_tax_id) VALUES (?, ?)",
            merged_rows,
        )
        conn.commit()

        # Post-build validation
        row = conn.execute("SELECT COUNT(*) FROM taxons").fetchone()
        taxon_count = row[0] if row else 0

        row = conn.execute(
            "SELECT COUNT(*) FROM taxons WHERE scientific_name != ''",
        ).fetchone()
        named_count = row[0] if row else 0
    except Exception:
        conn.close()
        if tmp_path.exists():
            tmp_path.unlink()
        raise
    else:
        conn.close()

    if taxon_count == 0:
        tmp_path.unlink()
        msg = f"Taxonomy build produced an empty database from {taxdump_dir}"
        raise TaxonomyBuildError(msg)

    if named_count == 0:
        tmp_path.unlink()
        msg = (
            f"Taxonomy build produced {taxon_count} taxons but none have scientific names. "
            f"names.dmp may be missing or corrupt: {taxdump_dir / 'names.dmp'}"
        )
        raise TaxonomyBuildError(msg)

    # Atomic replace: os.replace overwrites the destination atomically on POSIX
    # (same filesystem), so a crash can never leave the database missing.
    os.replace(tmp_path, sqlite_path)

    return sqlite_path


def _ensure_taxdump(
    taxonomy_dir: Path | str | None = None,
    *,
    refresh: TaxonomyRefresh | str = TaxonomyRefresh.MISSING,
    max_age_days: int = DEFAULT_TAXONOMY_MAX_AGE_DAYS,
) -> Path:
    """
    Ensure taxdump is downloaded, extracted, and SQLite is built when missing.

    Args:
        taxonomy_dir: Optional explicit taxonomy directory.

    Returns:
        Path to taxdump directory.
    """
    status = inspect_taxonomy(taxonomy_dir=taxonomy_dir)
    plan = plan_admin_taxonomy(
        status,
        refresh=TaxonomyRefresh(refresh),
        max_age_days=max_age_days,
    )
    return execute_taxonomy_plan(plan)


def ensure_taxonomy_available(
    taxonomy_dir: Path | str | None = None,
    *,
    refresh: TaxonomyRefresh | str = TaxonomyRefresh.MISSING,
    max_age_days: int = DEFAULT_TAXONOMY_MAX_AGE_DAYS,
) -> Path:
    """
    Ensure taxonomy database is downloaded and ready for use.

    Downloads NCBI taxdump only when required taxonomy files are missing.
    Existing taxonomy data is reused regardless of age for pipeline
    reproducibility.

    This is the public API for pre-downloading taxonomy data before
    distributed workers need it. Call this once on the submit node
    (e.g., in ENSURE_TAXONOMY) to avoid multiple workers racing to
    download simultaneously.

    Args:
        taxonomy_dir: Optional explicit taxonomy directory.
        refresh: Admin refresh policy for preparation.
        max_age_days: Freshness window used when refresh is ``stale``.

    Returns:
        Path to the taxdump directory containing taxonomy.sqlite
    """
    return _ensure_taxdump(
        taxonomy_dir=taxonomy_dir,
        refresh=refresh,
        max_age_days=max_age_days,
    )


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

    def get_taxon_by_name(self, scientific_name: str) -> Taxon | None:
        """
        Look up a taxon by its exact scientific name.

        Case-sensitive exact match against the NCBI scientific name field.
        Does not resolve merged taxa (merges are keyed by tax_id, not name).

        Args:
            scientific_name: Exact scientific name (e.g., "Adenoviridae")

        Returns:
            Taxon if found, None otherwise
        """
        row = self._conn.execute(
            "SELECT tax_id, scientific_name, rank, parent_tax_id "
            "FROM taxons WHERE scientific_name = ?",
            (scientific_name,),
        ).fetchone()
        if row is None:
            return None
        return Taxon(
            tax_id=row["tax_id"],
            scientific_name=row["scientific_name"],
            rank=row["rank"],
            parent_tax_id=row["parent_tax_id"],
        )

    @property
    def taxon_count(self) -> int:
        """Return the number of taxa in the database."""
        row = self._conn.execute("SELECT COUNT(*) FROM taxons").fetchone()
        return row[0] if row else 0

    def get_lineage(self, tax_id: int) -> list[Taxon]:
        """
        Get lineage for a tax ID in ancestor-to-taxon order.

        The returned lineage excludes NCBI's synthetic root taxid 1. It begins
        with the first taxon below root, such as ``131567`` for cellular
        organisms in the standard NCBI taxonomy.

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
        Get lineage as list of tax IDs from below NCBI root to leaf.

        Used by hits_to_report.py for LCA calculation.

        Example:
            >>> tax.get_lineage_ids(9606)  # Homo sapiens
            [131567, 2759, ..., 9605, 9606]
        """
        return [t.tax_id for t in self.get_lineage(tax_id)]

    def get_lineage_string(self, tax_id: int, *, include_unranked: bool = False) -> str:
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


class TaxonomyBuildError(Exception):
    """Raised when the taxonomy SQLite database cannot be built from .dmp files."""


class TaxonomyOfflineError(Exception):
    """Raised when taxonomy database is unavailable in offline mode."""


class ResourceUnavailableError(Exception):
    """Base class for errors raised when a required resource is unavailable."""

    resource_type = "Resource"
    recovery_hint = "remove the --sync flag"

    def __init__(
        self,
        resource_path: Path | str,
        operation: str,
        reason: str,
        original_error: Exception | None = None,
    ) -> None:
        self.resource_path = Path(resource_path)
        self.operation = operation
        self.reason = reason
        self.original_error = original_error
        super().__init__(
            f"{self.resource_type} unavailable during '{operation}': {reason}\n"
            f"Path: {resource_path}\n"
            f"The --sync flag requires this resource to be available.\n"
            f"To continue without it, {self.recovery_hint}.",
        )


class TaxonomyUnavailableError(ResourceUnavailableError):
    """
    Raised when taxonomy database is unavailable and --sync was requested.

    This error indicates that the taxonomy database could not be accessed
    (due to download failures, permission errors, corruption, etc.) and the
    user explicitly requested synchronization via --sync.

    Without --sync, this condition results in a warning and the script
    attempting a lazy download instead of failing.

    Attributes:
        taxdump_dir: Path to the taxdump directory
        operation: What operation was being attempted
        reason: Human-readable explanation of why taxonomy is unavailable
        original_error: The underlying exception that caused the failure
    """

    resource_type = "Taxonomy database"
    recovery_hint = "remove --sync (may fetch latest taxonomy)"

    def __init__(
        self,
        taxdump_dir: Path | str,
        operation: str,
        reason: str,
        original_error: Exception | None = None,
    ) -> None:
        super().__init__(
            resource_path=taxdump_dir,
            operation=operation,
            reason=reason,
            original_error=original_error,
        )

    @property
    def taxdump_dir(self) -> Path:
        """Alias for resource_path for backward compatibility."""
        return self.resource_path


def format_taxonomy_warning(  # noqa: PLR0913
    *,
    operation: str,
    context: str,
    error: Exception | None,
    taxdump_dir: Path | str,
    will_download: bool = True,
    will_continue: bool = False,
) -> str:
    """
    Format a detailed warning message for taxonomy database issues.

    This produces a "wide event" style warning with full context about what
    was being attempted, what went wrong, and what will happen next.

    Args:
        operation: Short description of the operation (e.g., "Loading taxonomy")
        context: Detailed context about what was being done
        error: The exception that occurred (if any)
        taxdump_dir: Path to the taxdump directory
        will_download: Whether a fresh download will be attempted
        will_continue: Whether NVD will continue with existing taxonomy data.

    Returns:
        Formatted multi-line warning string suitable for logging
    """
    # Validate inputs
    assert operation, "operation must not be empty"
    assert context, "context must not be empty"

    # Determine error details
    if error is not None:
        error_type = type(error).__name__
        error_msg = str(error)

        # Provide more helpful explanations for common errors
        if isinstance(error, PermissionError):
            explanation = "The taxonomy directory or files are not readable/writable by the current user."
        elif isinstance(error, sqlite3.DatabaseError):
            explanation = (
                "The taxonomy.sqlite file appears to be corrupted or is not "
                "a valid SQLite database."
            )
        elif isinstance(error, (OSError, IOError)) and "readonly" in str(error).lower():
            explanation = "The taxonomy directory is on a read-only filesystem."
        elif isinstance(error, urllib.error.URLError):
            explanation = "Failed to download taxonomy data from NCBI (network error)."
        elif isinstance(error, TaxonomyOfflineError):
            explanation = (
                "Taxonomy database not found and offline mode is enabled. "
                "The database must be pre-downloaded."
            )
        else:
            explanation = "An unexpected error occurred while accessing taxonomy data."

        error_section = f"""
Error:
  Type: {error_type}
  Message: {error_msg}
  Explanation: {explanation}
"""
    else:
        error_section = """
Issue:
  The pre-cached taxonomy database was not found or is stale.
"""

    if will_continue:
        consequence = """
Consequence:
  Continuing with the existing taxonomy database.
  WARNING: taxonomy may be older than the configured freshness window, so
  taxids, merged identifiers, and names may differ from current NCBI taxonomy.
  For reproducible pipeline runs, this is usually preferable to mutating a
  shared taxonomy directory during analysis.
"""
    elif will_download:
        consequence = """
Consequence:
  Attempting to download fresh taxonomy data from NCBI.
  WARNING: The latest NCBI taxonomy may include changes that affect results:
    - Taxids may have been merged or deprecated
    - Taxonomic classifications may have changed
    - New taxa may have been added
  For reproducible results, use 'nvd taxonomy ensure' before running the pipeline.
"""
    else:
        consequence = """
Consequence:
  Cannot proceed without taxonomy data.
  The --sync flag requires pre-cached taxonomy to be available.
"""

    return f"""
================================================================================
WARNING: TAXONOMY DATABASE ISSUE
================================================================================

Context:
  Operation: {operation}
  Details: {context}
  Taxdump directory: {taxdump_dir}
{error_section}{consequence}
Resolution:
  To pre-download taxonomy: nvd taxonomy ensure
  To require pre-cached taxonomy (fail if unavailable): add --sync to the command

================================================================================
""".strip()


@contextmanager
def open(  # noqa: A001, C901, PLR0912, PLR0913, PLR0915
    offline: bool | None = None,  # noqa: FBT001
    *,
    taxonomy_dir: Path | str | None = None,
    sync: bool = False,
    taxonomy_mode: TaxonomyMode | str | None = None,
    max_age_days: int | None = None,
    warn_callback: Callable[[str], None] | None = None,
) -> Generator[TaxonomyDB, None, None]:
    """
    Context manager for taxonomy database access.

    By default, downloads taxdump from NCBI only when required taxonomy files
    are missing. Existing taxonomy data is reused even when older than the
    freshness warning window so pipeline runs do not mutate shared references.
    Builds SQLite from .dmp files for fast lookups when needed.
    Provides access to taxopy for LCA calculations.

    Args:
        offline: If True, never attempt to download taxonomy data.
                 If None (default), checks NVD_TAXONOMY_OFFLINE env var.
                 Set to True for air-gapped environments with pre-cached data.
        taxonomy_dir: Optional explicit taxonomy directory.
        sync: If True, require pre-cached taxonomy and fail if unavailable.
              If False (default), attempt lazy download with a warning if
              required taxonomy files are missing. Stale-but-complete taxonomy
              is reused for reproducibility.
        taxonomy_mode: Pipeline taxonomy mode. ``read_only`` never mutates;
                       ``missing`` prepares taxonomy only when files are absent.
                       If unset, NVD_TAXONOMY_OFFLINE=1 selects ``read_only``;
                       otherwise the default is ``missing``.
        max_age_days: Freshness warning window.
        warn_callback: Optional callback function to receive warning messages.
                       If None, warnings are printed to stderr. The callback
                       should accept a single string argument.

    Raises:
        TaxonomyOfflineError: If offline=True and no cached taxonomy exists.
        TaxonomyUnavailableError: If sync=True and taxonomy cannot be loaded.

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

        # Sync mode (require pre-cached taxonomy)
        with taxonomy.open(sync=True) as tax:
            # Fails if taxonomy not pre-cached
            taxon = tax.get_taxon(9606)

        # Explicit taxonomy directory
        with taxonomy.open(taxonomy_dir="/shared/taxonomy") as tax:
            taxon = tax.get_taxon(9606)
    """

    def emit_warning(msg: str) -> None:
        if warn_callback is not None:
            warn_callback(msg)
        else:
            print(msg, file=sys.stderr)  # noqa: T201  # fallback warning output to stderr

    # Resolve offline mode from env var if not explicitly set.
    if offline is None:
        offline = os.environ.get("NVD_TAXONOMY_OFFLINE", "").lower() in ("1", "true")

    policy = TaxonomyPolicy.from_env_and_params(
        mode=taxonomy_mode,
        max_age_days=max_age_days,
        offline=offline or sync,
    )
    status = inspect_taxonomy(taxonomy_dir=taxonomy_dir)
    plan = plan_pipeline_taxonomy(status, policy)
    taxdump_dir = status.taxdump_dir
    sqlite_path = taxdump_dir / "taxonomy.sqlite"

    if plan.action is TaxonomyAction.FAIL:
        if sync:
            raise TaxonomyUnavailableError(
                taxdump_dir=taxdump_dir,
                operation="Loading taxonomy database",
                reason=plan.reason,
            )
        raise TaxonomyOfflineError(plan.reason)

    if plan.action is not TaxonomyAction.USE_EXISTING:
        if sync:
            raise TaxonomyUnavailableError(
                taxdump_dir=taxdump_dir,
                operation="Loading taxonomy database",
                reason=plan.reason,
            )
        warning = format_taxonomy_warning(
            operation="Loading taxonomy database",
            context=plan.reason,
            error=None,
            taxdump_dir=taxdump_dir,
            will_download=plan.action is TaxonomyAction.DOWNLOAD_AND_BUILD,
        )
        emit_warning(warning)

    elif status.is_stale(policy.max_age_days):
        warning = format_taxonomy_warning(
            operation="Loading taxonomy database",
            context=(
                f"Taxonomy is older than the configured freshness window "
                f"({policy.max_age_days} days); continuing with the existing database"
            ),
            error=None,
            taxdump_dir=taxdump_dir,
            will_download=False,
            will_continue=True,
        )
        emit_warning(warning)

    try:
        taxdump_dir = execute_taxonomy_plan(plan)
        sqlite_path = taxdump_dir / "taxonomy.sqlite"
    except Exception as e:
        if sync:
            raise TaxonomyUnavailableError(
                taxdump_dir=taxdump_dir,
                operation="Downloading/building taxonomy database",
                reason=str(e),
                original_error=e,
            ) from e
        warning = format_taxonomy_warning(
            operation="Downloading/building taxonomy database",
            context="Failed to download or build taxonomy database",
            error=e,
            taxdump_dir=taxdump_dir,
            will_download=False,
        )
        emit_warning(warning)
        raise

    # Open in read-only mode to avoid "unable to write to readonly database" errors
    # on distributed workers where the state directory may be read-only.
    # URI mode with ?mode=ro tells SQLite not to acquire write locks or create journals.
    try:
        conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True, timeout=30.0)
        conn.row_factory = sqlite3.Row
    except Exception as e:
        if sync:
            raise TaxonomyUnavailableError(
                taxdump_dir=taxdump_dir,
                operation="Opening taxonomy database",
                reason=str(e),
                original_error=e,
            ) from e
        else:
            warning = format_taxonomy_warning(
                operation="Opening taxonomy database",
                context=f"Failed to open {sqlite_path}",
                error=e,
                taxdump_dir=taxdump_dir,
                will_download=False,
            )
            emit_warning(warning)
            raise

    # Pre-use validation: verify the database is not empty or corrupt.
    # This catches cases where the file exists but is truncated (e.g., stale
    # NFS cache, interrupted copy to a distributed worker, or a build that
    # committed an empty transaction).
    try:
        row = conn.execute("SELECT COUNT(*) FROM taxons").fetchone()
        count = row[0] if row else 0
    except sqlite3.Error as e:
        conn.close()
        msg = f"Taxonomy database at {sqlite_path} appears corrupt: {e}"
        raise TaxonomyUnavailableError(
            taxdump_dir=taxdump_dir,
            operation="Validating taxonomy database",
            reason=msg,
            original_error=e,
        ) from e

    if count == 0:
        conn.close()
        msg = (
            f"Taxonomy database at {sqlite_path} contains 0 taxons. "
            f"The database file may be corrupt or was not fully built. "
            f"Try deleting {taxdump_dir} and re-running to trigger a fresh "
            f"download and rebuild."
        )
        raise TaxonomyUnavailableError(
            taxdump_dir=taxdump_dir,
            operation="Validating taxonomy database",
            reason=msg,
        )

    try:
        yield TaxonomyDB(conn, taxdump_dir)
    finally:
        conn.close()
