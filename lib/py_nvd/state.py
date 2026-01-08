"""
Run state management.

Tracks experiments, samples, database versions, and LabKey sync status.
Provides atomic operations that prevent duplicate processing and uploads.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import date, datetime
from pathlib import Path
from typing import Literal

from py_nvd.db import connect
from py_nvd.models import (
    Database,
    DatabaseResolution,
    DbType,
    Preset,
    ProcessedSample,
    ProcessedSampleStatus,
    Run,
    SampleLock,
    Status,
    TaxonomyVersion,
    Upload,
    UploadType,
)


class RunNotFoundError(Exception):
    """Raised when a run_id is not found in the database."""

    def __init__(self, run_id: str):
        self.run_id = run_id
        super().__init__(f"Run not found: {run_id}")


class SampleNotFoundError(Exception):
    """Raised when a sample is not found in the database."""

    def __init__(self, sample_id: str, sample_set_id: str):
        self.sample_id = sample_id
        self.sample_set_id = sample_set_id
        super().__init__(f"Sample not found: {sample_id} in sample set {sample_set_id}")


def compute_sample_set_id(sample_ids: list[str]) -> str:
    """
    Compute a deterministic ID from a set of sample IDs.

    The ID is a truncated SHA256 hash of the sorted, deduplicated sample IDs.
    This provides idempotency: same samples → same ID → duplicate detection.

    Args:
        sample_ids: List of sample identifiers (order doesn't matter)

    Returns:
        16-character hex string (first 64 bits of SHA256)
    """
    canonical = sorted(set(sample_ids))
    content = "\n".join(canonical)
    return hashlib.sha256(content.encode()).hexdigest()[:16]


def register_run(
    run_id: str,
    sample_set_id: str,
    experiment_id: int | None = None,
    state_dir: Path | str | None = None,
) -> Run | None:
    """
    Atomically register a new run.

    Returns the Run if registered, None if sample_set_id already exists.
    This prevents duplicate processing of the same sample set.

    Args:
        run_id: Nextflow's workflow.runName (unique per execution)
        sample_set_id: Hash of sample IDs (from compute_sample_set_id)
        experiment_id: Optional LabKey experiment ID for linkage
        state_dir: Optional state directory override

    Returns:
        The registered Run, or None if sample_set_id was already processed
    """
    with connect(state_dir) as conn:
        try:
            conn.execute(
                """INSERT INTO runs (run_id, sample_set_id, experiment_id, started_at, status)
                   VALUES (?, ?, ?, datetime('now'), 'running')""",
                (run_id, sample_set_id, experiment_id),
            )
            conn.commit()
            return get_run(run_id, state_dir)
        except sqlite3.IntegrityError:
            return None


def get_run(run_id: str, state_dir: Path | str | None = None) -> Run | None:
    """Get a run by ID."""
    with connect(state_dir) as conn:
        row = conn.execute("SELECT * FROM runs WHERE run_id = ?", (run_id,)).fetchone()
        if row:
            return Run.from_row(row)
        return None


def get_run_by_sample_set(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> Run | None:
    """Get a run by sample_set_id (useful for checking if samples were already processed)."""
    with connect(state_dir) as conn:
        row = conn.execute(
            "SELECT * FROM runs WHERE sample_set_id = ?",
            (sample_set_id,),
        ).fetchone()
        if row:
            return Run.from_row(row)
        return None


def update_run_id(
    sample_set_id: str,
    new_run_id: str,
    state_dir: Path | str | None = None,
) -> Run | None:
    """
    Update the run_id for an existing run identified by sample_set_id.

    Used when resuming a run with a new Nextflow run name. The sample lock
    system uses fingerprint (hostname + username) to allow the same user
    on the same machine to re-acquire locks with a new run_id. This function
    propagates that run_id update to both the runs table and sample_locks
    table for consistency.

    Args:
        sample_set_id: The sample set identifier (hash of sample IDs)
        new_run_id: The new run identifier (e.g., new workflow.runName)
        state_dir: Optional state directory override

    Returns:
        The updated Run, or None if no run exists for this sample_set_id
    """
    with connect(state_dir) as conn:
        # Get the old run_id first so we can update referencing tables
        old_run = conn.execute(
            "SELECT run_id FROM runs WHERE sample_set_id = ?",
            (sample_set_id,),
        ).fetchone()

        if not old_run:
            return None

        old_run_id = old_run[0]

        # Update runs table (primary key change)
        conn.execute(
            "UPDATE runs SET run_id = ? WHERE sample_set_id = ?",
            (new_run_id, sample_set_id),
        )

        # Update all tables that reference runs.run_id to keep foreign keys valid
        conn.execute(
            "UPDATE sample_locks SET run_id = ? WHERE run_id = ?",
            (new_run_id, old_run_id),
        )
        conn.execute(
            "UPDATE processed_samples SET run_id = ? WHERE run_id = ?",
            (new_run_id, old_run_id),
        )
        conn.execute(
            "UPDATE taxonomy_versions SET run_id = ? WHERE run_id = ?",
            (new_run_id, old_run_id),
        )

        conn.commit()
        return get_run(new_run_id, state_dir)


def list_runs(
    status: Status | None = None,
    since: date | datetime | None = None,
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[Run]:
    """
    List runs, most recent first.

    Args:
        status: Filter by run status ('running', 'completed', 'failed')
        since: Filter to runs started on or after this date/datetime
        limit: Maximum number of runs to return
        state_dir: Optional state directory override

    Returns:
        List of Run objects matching the filters, ordered by started_at DESC
    """
    with connect(state_dir) as conn:
        query = "SELECT * FROM runs"
        conditions: list[str] = []
        params: list[str | int] = []

        if status is not None:
            conditions.append("status = ?")
            params.append(status)

        if since is not None:
            # Convert to ISO format string for SQLite comparison
            if isinstance(since, datetime):
                since_str = since.isoformat(sep=" ")
            else:
                # date: compare from start of day
                since_str = since.isoformat()
            conditions.append("started_at >= ?")
            params.append(since_str)

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        query += " ORDER BY started_at DESC"

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = conn.execute(query, params).fetchall()
        return [Run.from_row(row) for row in rows]


def complete_run(
    run_id: str,
    status: Status = "completed",
    state_dir: Path | str | None = None,
    *,
    strict: bool = False,
) -> bool:
    """
    Mark a run as completed or failed.

    Args:
        run_id: The run identifier
        status: New status ("completed" or "failed")
        state_dir: Optional state directory override
        strict: If True, raise RunNotFoundError when run_id doesn't exist

    Returns:
        True if the run was found and updated, False if run_id not found.

    Raises:
        RunNotFoundError: If strict=True and run_id doesn't exist.
    """
    with connect(state_dir) as conn:
        cursor = conn.execute(
            "UPDATE runs SET completed_at = datetime('now'), status = ? WHERE run_id = ?",
            (status, run_id),
        )
        conn.commit()

        if cursor.rowcount == 0:
            if strict:
                raise RunNotFoundError(run_id)
            return False
        return True


def register_database(
    db_type: DbType,
    version: str,
    path: str | Path,
    state_dir: Path | str | None = None,
) -> Database:
    """
    Register a reference database.

    The path is canonicalized (resolved to absolute, symlinks followed)
    before storage to ensure consistent lookups.

    Args:
        db_type: Type of database (blast, stat, gottcha2)
        version: Version string for this database
        path: Path to the database (will be canonicalized)
        state_dir: Optional state directory override

    Returns:
        The registered Database
    """
    # Canonicalize path for consistent storage/lookup
    canonical_path = str(Path(path).resolve())

    with connect(state_dir) as conn:
        conn.execute(
            """INSERT OR REPLACE INTO databases 
               (db_type, version, path, checksum, registered_at)
               VALUES (?, ?, ?, NULL, datetime('now'))""",
            (db_type, version, canonical_path),
        )
        conn.commit()
        result = get_database_by_version(db_type, version, state_dir)
        assert result is not None  # We just inserted it
        return result


def get_database_by_version(
    db_type: DbType,
    version: str | None = None,
    state_dir: Path | str | None = None,
) -> Database | None:
    """
    Get database info by version.

    Args:
        db_type: Type of database (blast, stat, gottcha2)
        version: Version string. If None, returns most recently registered.
        state_dir: Optional state directory override

    Returns:
        Database if found, None otherwise
    """
    with connect(state_dir) as conn:
        if version is not None:
            row = conn.execute(
                "SELECT * FROM databases WHERE db_type = ? AND version = ?",
                (db_type, version),
            ).fetchone()
        else:
            row = conn.execute(
                "SELECT * FROM databases WHERE db_type = ? ORDER BY registered_at DESC LIMIT 1",
                (db_type,),
            ).fetchone()
        if row:
            return Database.from_row(row)
        return None


def get_database_by_path(
    db_type: DbType,
    path: str | Path,
    state_dir: Path | str | None = None,
) -> tuple[Database | None, str | None]:
    """
    Look up database by canonical path.

    When multiple versions are registered at the same path, returns the
    most recently registered version and includes a warning.

    Args:
        db_type: Type of database (blast, stat, gottcha2)
        path: Path to the database (will be canonicalized)
        state_dir: Optional state directory override

    Returns:
        (database, warning) tuple:
        - database: The Database if found, None otherwise
        - warning: Warning message if multiple versions found, None otherwise
    """
    # Canonicalize path (resolve symlinks, make absolute)
    canonical_path = str(Path(path).resolve())

    with connect(state_dir) as conn:
        # Use rowid as tiebreaker when timestamps are equal (same second)
        rows = conn.execute(
            """SELECT * FROM databases 
               WHERE db_type = ? AND path = ? 
               ORDER BY registered_at DESC, rowid DESC""",
            (db_type, canonical_path),
        ).fetchall()

        if not rows:
            return None, None

        # Build Database from first (newest) row
        row = rows[0]
        database = Database.from_row(row)

        # Warn if multiple versions exist for this path
        warning = None
        if len(rows) > 1:
            versions = [r["version"] for r in rows]
            warning = (
                f"Multiple versions registered for {db_type} at {canonical_path}: "
                f"{', '.join(versions)}. Using newest: {database.version}"
            )

        return database, warning


def get_databases_by_path(
    db_type: DbType,
    path: str | Path,
    state_dir: Path | str | None = None,
) -> list[Database]:
    """
    Get all database versions registered at a given path.

    Useful for debugging when multiple versions are registered at the same path.

    Args:
        db_type: Type of database (blast, stat, gottcha2)
        path: Path to the database (will be canonicalized)
        state_dir: Optional state directory override

    Returns:
        List of Database objects, ordered by registered_at DESC (newest first)
    """
    # Canonicalize path
    canonical_path = str(Path(path).resolve())

    with connect(state_dir) as conn:
        # Use rowid as tiebreaker when timestamps are equal (same second)
        rows = conn.execute(
            """SELECT * FROM databases 
               WHERE db_type = ? AND path = ? 
               ORDER BY registered_at DESC, rowid DESC""",
            (db_type, canonical_path),
        ).fetchall()

        return [Database.from_row(row) for row in rows]


def list_databases(state_dir: Path | str | None = None) -> list[Database]:
    """List all registered databases."""
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM databases ORDER BY db_type, registered_at DESC",
        ).fetchall()
        return [Database.from_row(row) for row in rows]


def validate_databases(
    required: list[tuple[DbType, str | None]],
    state_dir: Path | str | None = None,
) -> list[str]:
    """
    Validate that required databases are registered and paths exist.
    Returns list of error messages (empty if all valid).
    """
    errors = []
    for db_type, version in required:
        db = get_database_by_version(db_type, version, state_dir)
        if not db:
            errors.append(f"Database not registered: {db_type} ({version or 'any'})")
        elif not Path(db.path).exists():
            errors.append(f"Database path missing: {db.path}")
    return errors


def _resolve_single_database(
    db_type: DbType,
    path: Path | None,
    version: str | None,
    state_dir: Path | str | None = None,
) -> tuple[str | None, str | None, tuple[DbType, str, str] | None]:
    """
    Resolve version for a single database type.

    Implements the resolution logic:
    - path=None: skip (return None for all)
    - path set, version set: register/update, return provided version
    - path set, version None: lookup, return registered version or warn

    Args:
        db_type: Type of database (blast, stat, gottcha2)
        path: Database path (None means skip)
        version: User-provided version (None means auto-resolve)
        state_dir: Optional state directory override

    Returns:
        (resolved_version, warning, registration) tuple:
        - resolved_version: The version to use (None if unresolvable)
        - warning: Warning message if any (None otherwise)
        - registration: (db_type, version, path) if registered (None otherwise)
    """
    # Case 1: No path provided - nothing to resolve
    if path is None:
        return version, None, None

    canonical_path = str(Path(path).resolve())

    # Look up existing registration for this path
    existing_db, lookup_warning = get_database_by_path(db_type, path, state_dir)

    # Case 2: Path and version both provided - register/update
    if version is not None:
        if existing_db is None:
            # Not registered - auto-register
            register_database(db_type, version, canonical_path, state_dir=state_dir)
            return version, None, (db_type, version, canonical_path)

        if existing_db.version == version:
            # Already registered with same version - no-op
            # But still warn if there were multiple versions at this path
            return version, lookup_warning, None

        # Registered with different version - warn and add new registration
        # (old registration is preserved for history)
        warning = (
            f"{db_type} database at {canonical_path}: "
            f"registering new version '{version}' (previously registered as '{existing_db.version}')"
        )
        register_database(db_type, version, canonical_path, state_dir=state_dir)
        return version, warning, (db_type, version, canonical_path)

    # Case 3: Path provided but no version - auto-resolve from registry
    if existing_db is not None:
        # Found in registry - use registered version
        # Include lookup_warning if there were multiple versions
        return existing_db.version, lookup_warning, None

    # Not in registry and no version provided - warn
    warning = (
        f"{db_type} database path not registered: {canonical_path}. "
        f"Provide --{db_type.replace('_', '-')}-db-version or register with: "
        f"nvd state database register {db_type} {canonical_path} --version <version>"
    )
    return None, warning, None


def resolve_database_versions(
    blast_db: Path | None = None,
    blast_db_version: str | None = None,
    gottcha2_db: Path | None = None,
    gottcha2_db_version: str | None = None,
    stat_index: Path | None = None,
    stat_db_version: str | None = None,
    state_dir: Path | str | None = None,
) -> DatabaseResolution:
    """
    Resolve database versions from the registry based on provided paths.

    For each database type, applies the following logic:
    1. If path is None: skip (no resolution needed)
    2. If path is set and version is set:
       - Look up path in registry
       - If not found: auto-register (path, version)
       - If found with same version: no-op
       - If found with different version: warn and update registry
    3. If path is set and version is None:
       - Look up path in registry
       - If found: use registered version
       - If not found: warn "unregistered path"

    Args:
        blast_db: Path to BLAST database (None to skip)
        blast_db_version: User-provided BLAST version (None to auto-resolve)
        gottcha2_db: Path to GOTTCHA2 database (None to skip)
        gottcha2_db_version: User-provided GOTTCHA2 version (None to auto-resolve)
        stat_index: Path to STAT index file (None to skip)
        stat_db_version: User-provided STAT version (None to auto-resolve)
        state_dir: Optional state directory override

    Returns:
        DatabaseResolution with resolved versions, warnings, and registrations
    """
    resolution = DatabaseResolution()

    # Resolve BLAST
    blast_resolved, blast_warning, blast_reg = _resolve_single_database(
        "blast", blast_db, blast_db_version, state_dir
    )
    resolution.blast_db_version = blast_resolved
    if blast_warning:
        resolution.warnings.append(blast_warning)
    if blast_reg:
        resolution.auto_registered.append(blast_reg)

    # Resolve GOTTCHA2
    gottcha2_resolved, gottcha2_warning, gottcha2_reg = _resolve_single_database(
        "gottcha2", gottcha2_db, gottcha2_db_version, state_dir
    )
    resolution.gottcha2_db_version = gottcha2_resolved
    if gottcha2_warning:
        resolution.warnings.append(gottcha2_warning)
    if gottcha2_reg:
        resolution.auto_registered.append(gottcha2_reg)

    # Resolve STAT (uses stat_index as canonical path)
    stat_resolved, stat_warning, stat_reg = _resolve_single_database(
        "stat", stat_index, stat_db_version, state_dir
    )
    resolution.stat_db_version = stat_resolved
    if stat_warning:
        resolution.warnings.append(stat_warning)
    if stat_reg:
        resolution.auto_registered.append(stat_reg)

    return resolution


def mark_sample_completed(
    sample_id: str,
    sample_set_id: str,
    run_id: str,
    blast_db_version: str | None = None,
    stat_db_version: str | None = None,
    taxonomy_hash: str | None = None,
    state_dir: Path | str | None = None,
) -> ProcessedSample:
    """
    Mark a sample as completed (REGISTER_HITS finished).

    Call this after REGISTER_HITS succeeds. The sample is marked 'completed',
    indicating results exist and are ready for upload.

    This is an INSERT operation - the sample should not already exist in the
    processed_samples table.

    Args:
        sample_id: The sample identifier
        sample_set_id: The sample set this sample belongs to
        run_id: The run identifier (workflow.runName)
        blast_db_version: BLAST database version for provenance
        stat_db_version: STAT database version for provenance
        taxonomy_hash: Taxonomy database hash for provenance
        state_dir: Optional state directory override

    Returns:
        The created ProcessedSample record
    """
    with connect(state_dir) as conn:
        conn.execute(
            """INSERT INTO processed_samples
               (sample_id, sample_set_id, run_id, processed_at, blast_db_version,
                stat_db_version, taxonomy_hash, status)
               VALUES (?, ?, ?, datetime('now'), ?, ?, ?, 'completed')""",
            (
                sample_id,
                sample_set_id,
                run_id,
                blast_db_version,
                stat_db_version,
                taxonomy_hash,
            ),
        )
        conn.commit()
        result = get_processed_sample(sample_id, sample_set_id, state_dir)
        assert result is not None  # We just inserted it
        return result


# Backward compatibility alias
register_processed_sample = mark_sample_completed


def get_processed_sample(
    sample_id: str,
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> ProcessedSample | None:
    """Get a processed sample record by sample_id and sample_set_id."""
    with connect(state_dir) as conn:
        row = conn.execute(
            "SELECT * FROM processed_samples WHERE sample_id = ? AND sample_set_id = ?",
            (sample_id, sample_set_id),
        ).fetchone()
        if row:
            return ProcessedSample.from_row(row)
        return None


def get_samples_for_run(
    run_id: str,
    state_dir: Path | str | None = None,
) -> list[ProcessedSample]:
    """Get all processed samples for a specific run."""
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM processed_samples WHERE run_id = ? ORDER BY sample_id",
            (run_id,),
        ).fetchall()
        return [ProcessedSample.from_row(row) for row in rows]


def list_samples(
    run_id: str | None = None,
    status: ProcessedSampleStatus | None = None,
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[ProcessedSample]:
    """
    List processed samples with optional filters.

    Args:
        run_id: Filter to samples from a specific run
        status: Filter by processing status ('completed', 'uploaded', 'failed')
        limit: Maximum number of samples to return
        state_dir: Optional state directory override

    Returns:
        List of ProcessedSample objects, ordered by processed_at DESC
    """
    with connect(state_dir) as conn:
        query = "SELECT * FROM processed_samples"
        conditions: list[str] = []
        params: list[str | int] = []

        if run_id is not None:
            conditions.append("run_id = ?")
            params.append(run_id)

        if status is not None:
            conditions.append("status = ?")
            params.append(status)

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        query += " ORDER BY processed_at DESC"

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = conn.execute(query, params).fetchall()
        return [ProcessedSample.from_row(row) for row in rows]


def get_sample_history(
    sample_id: str,
    state_dir: Path | str | None = None,
) -> list[ProcessedSample]:
    """
    Get all processing records for a sample across all runs.

    Unlike get_processed_sample (which requires sample_set_id), this returns
    every time this sample_id was processed, enabling cross-run history view.

    Args:
        sample_id: The sample ID to look up
        state_dir: Optional state directory override

    Returns:
        List of ProcessedSample objects, ordered by processed_at DESC
    """
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM processed_samples WHERE sample_id = ? ORDER BY processed_at DESC",
            (sample_id,),
        ).fetchall()
        return [ProcessedSample.from_row(row) for row in rows]


def was_sample_ever_processed(
    sample_id: str,
    state_dir: Path | str | None = None,
) -> bool:
    """Check if a sample was ever processed in any run."""
    with connect(state_dir) as conn:
        row = conn.execute(
            "SELECT 1 FROM processed_samples WHERE sample_id = ? LIMIT 1",
            (sample_id,),
        ).fetchone()
        return row is not None


def get_samples_needing_upload(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> list[ProcessedSample]:
    """
    Get samples that are completed but not yet uploaded.

    These are samples ready for LabKey upload.

    Args:
        sample_set_id: The sample set to query
        state_dir: Optional state directory override

    Returns:
        List of ProcessedSample objects with status='completed'
    """
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM processed_samples "
            "WHERE sample_set_id = ? AND status = 'completed' "
            "ORDER BY sample_id",
            (sample_set_id,),
        ).fetchall()
        return [ProcessedSample.from_row(row) for row in rows]


def get_uploaded_sample_ids(
    sample_ids: list[str],
    state_dir: Path | str | None = None,
) -> set[str]:
    """
    Check which of the given sample IDs have been uploaded (in any sample set).

    Used at run start to warn about/filter already-uploaded samples.

    Args:
        sample_ids: List of sample IDs to check
        state_dir: Optional state directory override

    Returns:
        Set of sample IDs that have status='uploaded'
    """
    if not sample_ids:
        return set()

    placeholders = ",".join("?" * len(sample_ids))
    with connect(state_dir) as conn:
        rows = conn.execute(
            f"SELECT DISTINCT sample_id FROM processed_samples "  # noqa: S608
            f"WHERE sample_id IN ({placeholders}) AND status = 'uploaded'",
            sample_ids,
        ).fetchall()
        return {row["sample_id"] for row in rows}


def complete_sample(
    sample_id: str,
    sample_set_id: str,
    status: ProcessedSampleStatus = "completed",
    state_dir: Path | str | None = None,
    *,
    strict: bool = False,
) -> bool:
    """
    Mark a processed sample as completed or failed.

    Args:
        sample_id: The sample identifier
        sample_set_id: The sample set this sample belongs to
        status: New status ("completed" or "failed")
        state_dir: Optional state directory override
        strict: If True, raise SampleNotFoundError when sample doesn't exist

    Returns:
        True if the sample was found and updated, False if not found.

    Raises:
        SampleNotFoundError: If strict=True and sample doesn't exist.
    """
    with connect(state_dir) as conn:
        cursor = conn.execute(
            "UPDATE processed_samples SET status = ? WHERE sample_id = ? AND sample_set_id = ?",
            (status, sample_id, sample_set_id),
        )
        conn.commit()

        if cursor.rowcount == 0:
            if strict:
                raise SampleNotFoundError(sample_id, sample_set_id)
            return False
        return True


def mark_sample_uploaded(
    sample_id: str,
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> bool:
    """
    Mark a sample as uploaded to LabKey.

    This is the terminal state for samples in LabKey-enabled runs.
    Should be called AFTER record_upload() succeeds.

    Args:
        sample_id: The sample identifier
        sample_set_id: The sample set this sample belongs to
        state_dir: Optional state directory override

    Returns:
        True if updated, False if sample not found
    """
    with connect(state_dir) as conn:
        cursor = conn.execute(
            "UPDATE processed_samples SET status = 'uploaded' "
            "WHERE sample_id = ? AND sample_set_id = ?",
            (sample_id, sample_set_id),
        )
        conn.commit()
        return cursor.rowcount > 0


def record_upload(
    sample_id: str,
    sample_set_id: str,
    upload_type: UploadType,
    upload_target: str,
    content_hash: str,
    target_metadata: dict | None = None,
    state_dir: Path | str | None = None,
) -> Upload:
    """
    Record a successful upload.

    Call this AFTER the upload succeeds to ensure we only record
    uploads that actually completed.

    Args:
        sample_id: The sample identifier
        sample_set_id: The sample set this sample belongs to
        upload_type: Type of upload (blast, blast_fasta, gottcha2, gottcha2_fasta)
        upload_target: Where it was uploaded (labkey, local, globus, etc.)
        content_hash: SHA256 of the uploaded content
        target_metadata: Optional target-specific data (e.g., experiment_id, row_id)
    """
    metadata_json = json.dumps(target_metadata) if target_metadata else None
    with connect(state_dir) as conn:
        conn.execute(
            """INSERT INTO uploads 
               (sample_id, sample_set_id, upload_type, upload_target, content_hash, 
                uploaded_at, target_metadata)
               VALUES (?, ?, ?, ?, ?, datetime('now'), ?)""",
            (
                sample_id,
                sample_set_id,
                upload_type,
                upload_target,
                content_hash,
                metadata_json,
            ),
        )
        conn.commit()
        result = get_upload(
            sample_id,
            sample_set_id,
            upload_type,
            upload_target,
            state_dir,
        )
        assert result is not None  # We just inserted it
        return result


def get_upload(
    sample_id: str,
    sample_set_id: str,
    upload_type: UploadType,
    upload_target: str,
    state_dir: Path | str | None = None,
) -> Upload | None:
    """Get a specific upload record."""
    with connect(state_dir) as conn:
        row = conn.execute(
            """SELECT * FROM uploads 
               WHERE sample_id = ? AND sample_set_id = ? AND upload_type = ? AND upload_target = ?""",
            (sample_id, sample_set_id, upload_type, upload_target),
        ).fetchone()
        if row:
            return Upload.from_row(row)
        return None


def check_upload(
    sample_id: str,
    sample_set_id: str,
    upload_type: UploadType,
    upload_target: str,
    content_hash: str,
    state_dir: Path | str | None = None,
) -> Literal["not_uploaded", "already_uploaded", "content_changed"]:
    """
    Check upload status for a sample before uploading.

    Returns:
        "not_uploaded": Safe to upload
        "already_uploaded": Same content already uploaded, skip
        "content_changed": Different content was uploaded previously
    """
    existing = get_upload(
        sample_id,
        sample_set_id,
        upload_type,
        upload_target,
        state_dir,
    )

    if existing is None:
        return "not_uploaded"
    if existing.content_hash == content_hash:
        return "already_uploaded"
    return "content_changed"


def get_uploads_for_sample(
    sample_id: str,
    state_dir: Path | str | None = None,
) -> list[Upload]:
    """
    Get all uploads for a sample across all runs.

    This is the key function for cross-run deduplication — it answers
    "was this sample ever uploaded?" regardless of which sample set it was in.
    """
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM uploads WHERE sample_id = ? ORDER BY uploaded_at",
            (sample_id,),
        ).fetchall()
        return [Upload.from_row(row) for row in rows]


def list_uploads(
    sample_id: str | None = None,
    upload_type: UploadType | None = None,
    upload_target: str | None = None,
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[Upload]:
    """
    List uploads with optional filters.

    Args:
        sample_id: Filter to uploads for a specific sample
        upload_type: Filter by upload type ('blast', 'blast_fasta', 'gottcha2', 'gottcha2_fasta')
        upload_target: Filter by target ('labkey', 'local', etc.)
        limit: Maximum number of uploads to return
        state_dir: Optional state directory override

    Returns:
        List of Upload objects, ordered by uploaded_at DESC
    """
    with connect(state_dir) as conn:
        query = "SELECT * FROM uploads"
        conditions: list[str] = []
        params: list[str | int] = []

        if sample_id is not None:
            conditions.append("sample_id = ?")
            params.append(sample_id)

        if upload_type is not None:
            conditions.append("upload_type = ?")
            params.append(upload_type)

        if upload_target is not None:
            conditions.append("upload_target = ?")
            params.append(upload_target)

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        query += " ORDER BY uploaded_at DESC"

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = conn.execute(query, params).fetchall()
        return [Upload.from_row(row) for row in rows]


def was_sample_ever_uploaded(
    sample_id: str,
    upload_type: UploadType | None = None,
    upload_target: str | None = None,
    state_dir: Path | str | None = None,
) -> bool:
    """
    Check if a sample was ever uploaded in any run.

    This is the key function for cross-run sample deduplication.

    Args:
        sample_id: The sample to check
        upload_type: Optional filter by upload type
        upload_target: Optional filter by target (e.g., 'labkey')

    Returns:
        True if the sample was uploaded (with optional filters), False otherwise
    """
    with connect(state_dir) as conn:
        query = "SELECT 1 FROM uploads WHERE sample_id = ?"
        params: list[str] = [sample_id]

        if upload_type is not None:
            query += " AND upload_type = ?"
            params.append(upload_type)

        if upload_target is not None:
            query += " AND upload_target = ?"
            params.append(upload_target)

        query += " LIMIT 1"
        row = conn.execute(query, params).fetchone()
        return row is not None


def hash_upload_content(data: dict | list) -> str:
    """
    Generate a stable hash of upload content for idempotency checking.

    Uses JSON serialization with sorted keys for deterministic output.
    """
    content = json.dumps(data, sort_keys=True, default=str)
    return hashlib.sha256(content.encode()).hexdigest()


def record_taxonomy_version(
    run_id: str,
    file_hash: str,
    ncbi_rebuild_timestamp: str | None,
    file_size: int,
    downloaded_at: str,
    state_dir: Path | str | None = None,
) -> None:
    """
    Record which taxonomy database version was used for a run.

    Call this at the start of a run to enable reproducibility tracking.
    """
    with connect(state_dir) as conn:
        conn.execute(
            """INSERT OR REPLACE INTO taxonomy_versions 
               (run_id, file_hash, ncbi_rebuild_timestamp, file_size, downloaded_at, recorded_at)
               VALUES (?, ?, ?, ?, ?, datetime('now'))""",
            (run_id, file_hash, ncbi_rebuild_timestamp, file_size, downloaded_at),
        )
        conn.commit()


def get_taxonomy_version(
    run_id: str,
    state_dir: Path | str | None = None,
) -> TaxonomyVersion | None:
    """Get the taxonomy version used for a specific run."""
    with connect(state_dir) as conn:
        row = conn.execute(
            "SELECT * FROM taxonomy_versions WHERE run_id = ?",
            (run_id,),
        ).fetchone()
        if row:
            return TaxonomyVersion.from_row(row)
        return None


def check_taxonomy_drift(
    run_id_a: str,
    run_id_b: str,
    state_dir: Path | str | None = None,
) -> bool:
    """
    Check if two runs used different taxonomy versions.

    Returns True if versions differ (drift detected), False if same.
    Useful for flagging when comparing results between runs.
    """
    version_a = get_taxonomy_version(run_id_a, state_dir)
    version_b = get_taxonomy_version(run_id_b, state_dir)

    if version_a is None or version_b is None:
        return True  # Can't compare if one is missing

    return version_a.file_hash != version_b.file_hash


def register_preset(
    name: str,
    params: dict[str, str | int | float | bool | None],
    description: str | None = None,
    state_dir: Path | str | None = None,
) -> Preset:
    """
    Register or update a parameter preset.

    If a preset with this name exists, it is updated (params and description
    are replaced, created_at is preserved).

    Args:
        name: Unique name for the preset
        params: Dictionary of parameter names to values
        description: Optional description of the preset's purpose
        state_dir: Optional state directory override

    Returns:
        The created or updated Preset
    """
    params_json = json.dumps(params, sort_keys=True)

    with connect(state_dir) as conn:
        # Check if exists
        existing = conn.execute(
            "SELECT created_at FROM presets WHERE name = ?",
            (name,),
        ).fetchone()

        if existing:
            # Update existing preset
            conn.execute(
                """UPDATE presets 
                   SET params = ?, description = ?, updated_at = datetime('now')
                   WHERE name = ?""",
                (params_json, description, name),
            )
        else:
            # Insert new preset
            conn.execute(
                """INSERT INTO presets (name, params, description, created_at, updated_at)
                   VALUES (?, ?, ?, datetime('now'), datetime('now'))""",
                (name, params_json, description),
            )

        conn.commit()

        # Return the preset we just created/updated
        result = get_preset(name, state_dir)
        assert result is not None  # We just inserted/updated it
        return result


def get_preset(name: str, state_dir: Path | str | None = None) -> Preset | None:
    """
    Get a preset by name.

    Args:
        name: Name of the preset to retrieve
        state_dir: Optional state directory override

    Returns:
        The Preset if found, None otherwise
    """
    with connect(state_dir) as conn:
        row = conn.execute("SELECT * FROM presets WHERE name = ?", (name,)).fetchone()

        if row:
            return Preset.from_row(row)
        return None


def list_presets(state_dir: Path | str | None = None) -> list[Preset]:
    """
    List all presets, most recently created first.

    Args:
        state_dir: Optional state directory override

    Returns:
        List of all Presets, ordered by creation date (newest first)
    """
    with connect(state_dir) as conn:
        rows = conn.execute("SELECT * FROM presets ORDER BY created_at DESC").fetchall()

        return [Preset.from_row(row) for row in rows]


def delete_preset(name: str, state_dir: Path | str | None = None) -> bool:
    """
    Delete a preset by name.

    Args:
        name: Name of the preset to delete
        state_dir: Optional state directory override

    Returns:
        True if the preset was deleted, False if it wasn't found
    """
    with connect(state_dir) as conn:
        cursor = conn.execute("DELETE FROM presets WHERE name = ?", (name,))
        conn.commit()
        return cursor.rowcount > 0


# =============================================================================
# Sample Locks
# =============================================================================

DEFAULT_LOCK_TTL_HOURS = 72


def get_machine_fingerprint() -> tuple[str, str]:
    """
    Return (hostname, username) for lock ownership.

    Used to detect when the same run_id is resumed from different machines,
    which would indicate a conflict (two users resuming the same run).
    """
    import os
    import socket

    hostname = socket.gethostname()
    username = os.environ.get("USER") or os.getlogin()
    return (hostname, username)


def get_lock(
    sample_id: str,
    state_dir: Path | str | None = None,
    *,
    include_expired: bool = False,
) -> SampleLock | None:
    """
    Get lock info for a sample.

    Args:
        sample_id: The sample to check
        state_dir: Optional state directory override
        include_expired: If True, return expired locks too

    Returns:
        SampleLock if locked (and not expired, unless include_expired=True), else None
    """
    with connect(state_dir) as conn:
        row = conn.execute(
            "SELECT * FROM sample_locks WHERE sample_id = ?",
            (sample_id,),
        ).fetchone()

        if row is None:
            return None

        lock = SampleLock.from_row(row)

        if not include_expired and lock.is_expired:
            return None

        return lock


def get_locks_for_run(
    run_id: str,
    state_dir: Path | str | None = None,
    *,
    include_expired: bool = False,
) -> list[SampleLock]:
    """
    Get all locks held by a run.

    Args:
        run_id: The run to query
        state_dir: Optional state directory override
        include_expired: If True, include expired locks

    Returns:
        List of SampleLock objects for this run
    """
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM sample_locks WHERE run_id = ? ORDER BY sample_id",
            (run_id,),
        ).fetchall()

        locks = [SampleLock.from_row(row) for row in rows]

        if not include_expired:
            locks = [lock for lock in locks if not lock.is_expired]

        return locks


def list_locks(
    state_dir: Path | str | None = None,
    *,
    include_expired: bool = False,
) -> list[SampleLock]:
    """
    List all sample locks.

    Args:
        state_dir: Optional state directory override
        include_expired: If True, include expired locks

    Returns:
        List of all SampleLock objects
    """
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM sample_locks ORDER BY locked_at DESC"
        ).fetchall()

        locks = [SampleLock.from_row(row) for row in rows]

        if not include_expired:
            locks = [lock for lock in locks if not lock.is_expired]

        return locks


def acquire_sample_lock(
    sample_id: str,
    run_id: str,
    ttl_hours: int = DEFAULT_LOCK_TTL_HOURS,
    state_dir: Path | str | None = None,
) -> SampleLock | None:
    """
    Attempt to acquire a lock for a single sample.

    Lock acquisition logic (fingerprint = hostname + username):
    - No existing lock → acquire
    - Existing lock, same fingerprint → acquire (same user can always re-acquire,
      handles Nextflow -resume which generates new run_id each time)
    - Existing lock, different fingerprint, expired → acquire (steal expired)
    - Existing lock, different fingerprint, active → blocked

    Args:
        sample_id: The sample to lock
        run_id: The run acquiring the lock
        ttl_hours: Lock TTL in hours (default 72)
        state_dir: Optional state directory override

    Returns:
        None if acquired successfully, or the blocking SampleLock if blocked
    """
    hostname, username = get_machine_fingerprint()
    now = datetime.now().astimezone()
    expires = now + __import__("datetime").timedelta(hours=ttl_hours)

    now_str = now.isoformat()
    expires_str = expires.isoformat()

    with connect(state_dir) as conn:
        existing = conn.execute(
            "SELECT * FROM sample_locks WHERE sample_id = ?",
            (sample_id,),
        ).fetchone()

        if existing is None:
            # No lock - acquire
            conn.execute(
                """INSERT INTO sample_locks 
                   (sample_id, run_id, hostname, username, locked_at, expires_at)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                (sample_id, run_id, hostname, username, now_str, expires_str),
            )
            conn.commit()
            return None

        existing_lock = SampleLock.from_row(existing)

        # Check fingerprint first - same user on same host can always re-acquire
        # This handles Nextflow -resume which generates a new run_id each time
        if existing_lock.fingerprint == (hostname, username):
            # Same user/host - allow acquisition, update run_id and refresh TTL
            conn.execute(
                """UPDATE sample_locks 
                   SET run_id = ?, locked_at = ?, expires_at = ?
                   WHERE sample_id = ?""",
                (run_id, now_str, expires_str, sample_id),
            )
            conn.commit()
            return None

        # Different fingerprint - check if we can steal
        if existing_lock.is_expired:
            # Expired lock from another user - steal it
            conn.execute(
                """UPDATE sample_locks 
                   SET run_id = ?, hostname = ?, username = ?, locked_at = ?, expires_at = ?
                   WHERE sample_id = ?""",
                (run_id, hostname, username, now_str, expires_str, sample_id),
            )
            conn.commit()
            return None
        else:
            # Active lock held by different user/host - blocked
            return existing_lock


def acquire_sample_locks(
    sample_ids: list[str],
    run_id: str,
    ttl_hours: int = DEFAULT_LOCK_TTL_HOURS,
    state_dir: Path | str | None = None,
) -> tuple[list[str], list[SampleLock]]:
    """
    Attempt to acquire locks for multiple samples.

    Args:
        sample_ids: List of samples to lock
        run_id: The run acquiring the locks
        ttl_hours: Lock TTL in hours (default 72)
        state_dir: Optional state directory override

    Returns:
        Tuple of (acquired, conflicts) where:
        - acquired: list of sample_ids successfully locked
        - conflicts: list of SampleLock objects for blocked samples
    """
    acquired: list[str] = []
    conflicts: list[SampleLock] = []

    for sample_id in sample_ids:
        result = acquire_sample_lock(sample_id, run_id, ttl_hours, state_dir)
        if result is None:
            acquired.append(sample_id)
        else:
            conflicts.append(result)

    return acquired, conflicts


def release_sample_lock(
    sample_id: str,
    run_id: str,
    state_dir: Path | str | None = None,
) -> bool:
    """
    Release a lock for a single sample if held by this run.

    Args:
        sample_id: The sample to unlock
        run_id: The run releasing the lock
        state_dir: Optional state directory override

    Returns:
        True if released, False if not held by this run
    """
    with connect(state_dir) as conn:
        cursor = conn.execute(
            "DELETE FROM sample_locks WHERE sample_id = ? AND run_id = ?",
            (sample_id, run_id),
        )
        conn.commit()
        return cursor.rowcount > 0


def release_sample_locks(
    sample_ids: list[str],
    run_id: str,
    state_dir: Path | str | None = None,
) -> int:
    """
    Release locks for multiple samples held by this run.

    Args:
        sample_ids: List of samples to unlock
        run_id: The run releasing the locks
        state_dir: Optional state directory override

    Returns:
        Count of locks released
    """
    if not sample_ids:
        return 0

    placeholders = ",".join("?" * len(sample_ids))
    with connect(state_dir) as conn:
        cursor = conn.execute(
            f"DELETE FROM sample_locks WHERE sample_id IN ({placeholders}) AND run_id = ?",  # noqa: S608
            [*sample_ids, run_id],
        )
        conn.commit()
        return cursor.rowcount


def release_all_locks_for_run(
    run_id: str,
    state_dir: Path | str | None = None,
) -> int:
    """
    Release all locks held by a run.

    Args:
        run_id: The run to release locks for
        state_dir: Optional state directory override

    Returns:
        Count of locks released
    """
    with connect(state_dir) as conn:
        cursor = conn.execute(
            "DELETE FROM sample_locks WHERE run_id = ?",
            (run_id,),
        )
        conn.commit()
        return cursor.rowcount


def cleanup_expired_locks(state_dir: Path | str | None = None) -> int:
    """
    Remove all expired locks.

    Args:
        state_dir: Optional state directory override

    Returns:
        Count of locks removed
    """
    now = datetime.now().astimezone().isoformat()
    with connect(state_dir) as conn:
        cursor = conn.execute(
            "DELETE FROM sample_locks WHERE expires_at < ?",
            (now,),
        )
        conn.commit()
        return cursor.rowcount


# =============================================================================
# Run Status (with TTL-based failure detection)
# =============================================================================


def get_effective_run_status(
    run: Run,
    ttl_hours: int = DEFAULT_LOCK_TTL_HOURS,
) -> Status:
    """
    Get the effective status of a run, treating stale 'running' as 'failed'.

    A run that has been 'running' for longer than the TTL is considered
    effectively failed (crashed/interrupted without cleanup).

    Args:
        run: The Run to check
        ttl_hours: Hours after which a 'running' run is considered failed

    Returns:
        Effective status: 'running', 'completed', or 'failed'
    """
    if run.status != "running":
        return run.status

    started = datetime.fromisoformat(run.started_at.replace("Z", "+00:00"))
    if started.tzinfo is None:
        started = started.replace(tzinfo=__import__("datetime").timezone.utc)

    now = datetime.now(__import__("datetime").timezone.utc)
    age_hours = (now - started).total_seconds() / 3600

    if age_hours > ttl_hours:
        return "failed"

    return "running"


def mark_stale_runs_failed(
    ttl_hours: int = DEFAULT_LOCK_TTL_HOURS,
    state_dir: Path | str | None = None,
) -> int:
    """
    Update 'running' runs older than TTL to 'failed'.

    This persists the derived status to the database for clarity.
    Called by `nvd state cleanup` to clean up stale runs.

    Args:
        ttl_hours: Hours after which a 'running' run is considered failed
        state_dir: Optional state directory override

    Returns:
        Count of runs updated
    """
    cutoff = datetime.now(__import__("datetime").timezone.utc) - __import__(
        "datetime"
    ).timedelta(hours=ttl_hours)
    cutoff_str = cutoff.isoformat()

    with connect(state_dir) as conn:
        cursor = conn.execute(
            """UPDATE runs SET status = 'failed', completed_at = datetime('now')
               WHERE status = 'running' AND started_at < ?""",
            (cutoff_str,),
        )
        conn.commit()
        return cursor.rowcount
