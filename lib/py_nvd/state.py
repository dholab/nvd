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
    DbType,
    Preset,
    ProcessedSample,
    ProcessedSampleStatus,
    Run,
    Status,
    TaxonomyVersion,
    Upload,
    UploadType,
)

# --- Run Management ---


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
            return Run(
                run_id=row["run_id"],
                sample_set_id=row["sample_set_id"],
                started_at=row["started_at"],
                status=row["status"],
                experiment_id=row["experiment_id"],
                completed_at=row["completed_at"],
            )
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
            return Run(
                run_id=row["run_id"],
                sample_set_id=row["sample_set_id"],
                started_at=row["started_at"],
                status=row["status"],
                experiment_id=row["experiment_id"],
                completed_at=row["completed_at"],
            )
        return None


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
        return [
            Run(
                run_id=row["run_id"],
                sample_set_id=row["sample_set_id"],
                started_at=row["started_at"],
                status=row["status"],
                experiment_id=row["experiment_id"],
                completed_at=row["completed_at"],
            )
            for row in rows
        ]


def complete_run(
    run_id: str,
    status: Status = "completed",
    state_dir: Path | str | None = None,
) -> None:
    """Mark a run as completed or failed."""
    with connect(state_dir) as conn:
        conn.execute(
            "UPDATE runs SET completed_at = datetime('now'), status = ? WHERE run_id = ?",
            (status, run_id),
        )
        conn.commit()


# --- Database Registry ---


def register_database(
    db_type: DbType,
    version: str,
    path: str,
    checksum: str | None = None,
    state_dir: Path | str | None = None,
) -> Database:
    """Register a reference database."""
    with connect(state_dir) as conn:
        conn.execute(
            """INSERT OR REPLACE INTO databases 
               (db_type, version, path, checksum, registered_at)
               VALUES (?, ?, ?, ?, datetime('now'))""",
            (db_type, version, path, checksum),
        )
        conn.commit()
        result = get_database(db_type, version, state_dir)
        assert result is not None  # We just inserted it
        return result


def get_database(
    db_type: DbType,
    version: str | None = None,
    state_dir: Path | str | None = None,
) -> Database | None:
    """Get database info. If version is None, returns most recently registered."""
    with connect(state_dir) as conn:
        if version:
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
            return Database(
                db_type=row["db_type"],
                version=row["version"],
                path=row["path"],
                registered_at=row["registered_at"],
                checksum=row["checksum"],
            )
        return None


def list_databases(state_dir: Path | str | None = None) -> list[Database]:
    """List all registered databases."""
    with connect(state_dir) as conn:
        rows = conn.execute(
            "SELECT * FROM databases ORDER BY db_type, registered_at DESC",
        ).fetchall()
        return [
            Database(
                db_type=row["db_type"],
                version=row["version"],
                path=row["path"],
                registered_at=row["registered_at"],
                checksum=row["checksum"],
            )
            for row in rows
        ]


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
        db = get_database(db_type, version, state_dir)
        if not db:
            errors.append(f"Database not registered: {db_type} ({version or 'any'})")
        elif not Path(db.path).exists():
            errors.append(f"Database path missing: {db.path}")
    return errors


# --- Processed Sample Tracking ---


def register_processed_sample(
    sample_id: str,
    sample_set_id: str,
    run_id: str,
    blast_db_version: str | None = None,
    stat_db_version: str | None = None,
    taxonomy_hash: str | None = None,
    state_dir: Path | str | None = None,
) -> ProcessedSample:
    """
    Register a sample as being processed.

    Call this when processing starts for a sample. The sample is initially
    marked as 'processing' and should be updated to 'completed' or 'failed'
    when processing finishes.
    """
    with connect(state_dir) as conn:
        conn.execute(
            """INSERT INTO processed_samples 
               (sample_id, sample_set_id, run_id, processed_at, blast_db_version, 
                stat_db_version, taxonomy_hash, status)
               VALUES (?, ?, ?, datetime('now'), ?, ?, ?, 'processing')""",
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
            return ProcessedSample(
                sample_id=row["sample_id"],
                sample_set_id=row["sample_set_id"],
                run_id=row["run_id"],
                processed_at=row["processed_at"],
                status=row["status"],
                blast_db_version=row["blast_db_version"],
                stat_db_version=row["stat_db_version"],
                taxonomy_hash=row["taxonomy_hash"],
            )
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
        return [
            ProcessedSample(
                sample_id=row["sample_id"],
                sample_set_id=row["sample_set_id"],
                run_id=row["run_id"],
                processed_at=row["processed_at"],
                status=row["status"],
                blast_db_version=row["blast_db_version"],
                stat_db_version=row["stat_db_version"],
                taxonomy_hash=row["taxonomy_hash"],
            )
            for row in rows
        ]


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
        status: Filter by processing status ('processing', 'completed', 'failed')
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
        return [
            ProcessedSample(
                sample_id=row["sample_id"],
                sample_set_id=row["sample_set_id"],
                run_id=row["run_id"],
                processed_at=row["processed_at"],
                status=row["status"],
                blast_db_version=row["blast_db_version"],
                stat_db_version=row["stat_db_version"],
                taxonomy_hash=row["taxonomy_hash"],
            )
            for row in rows
        ]


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
        return [
            ProcessedSample(
                sample_id=row["sample_id"],
                sample_set_id=row["sample_set_id"],
                run_id=row["run_id"],
                processed_at=row["processed_at"],
                status=row["status"],
                blast_db_version=row["blast_db_version"],
                stat_db_version=row["stat_db_version"],
                taxonomy_hash=row["taxonomy_hash"],
            )
            for row in rows
        ]


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


def complete_sample(
    sample_id: str,
    sample_set_id: str,
    status: ProcessedSampleStatus = "completed",
    state_dir: Path | str | None = None,
) -> None:
    """Mark a processed sample as completed or failed."""
    with connect(state_dir) as conn:
        conn.execute(
            "UPDATE processed_samples SET status = ? WHERE sample_id = ? AND sample_set_id = ?",
            (status, sample_id, sample_set_id),
        )
        conn.commit()


# --- Upload Tracking (Target-Agnostic) ---


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
            metadata = (
                json.loads(row["target_metadata"]) if row["target_metadata"] else None
            )
            return Upload(
                sample_id=row["sample_id"],
                sample_set_id=row["sample_set_id"],
                upload_type=row["upload_type"],
                upload_target=row["upload_target"],
                content_hash=row["content_hash"],
                uploaded_at=row["uploaded_at"],
                target_metadata=metadata,
            )
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
        return [
            Upload(
                sample_id=row["sample_id"],
                sample_set_id=row["sample_set_id"],
                upload_type=row["upload_type"],
                upload_target=row["upload_target"],
                content_hash=row["content_hash"],
                uploaded_at=row["uploaded_at"],
                target_metadata=json.loads(row["target_metadata"])
                if row["target_metadata"]
                else None,
            )
            for row in rows
        ]


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
        return [
            Upload(
                sample_id=row["sample_id"],
                sample_set_id=row["sample_set_id"],
                upload_type=row["upload_type"],
                upload_target=row["upload_target"],
                content_hash=row["content_hash"],
                uploaded_at=row["uploaded_at"],
                target_metadata=json.loads(row["target_metadata"])
                if row["target_metadata"]
                else None,
            )
            for row in rows
        ]


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
        params: list = [sample_id]

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


# --- Taxonomy Version Tracking ---


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
            return TaxonomyVersion(
                run_id=row["run_id"],
                file_hash=row["file_hash"],
                ncbi_rebuild_timestamp=row["ncbi_rebuild_timestamp"],
                file_size=row["file_size"],
                downloaded_at=row["downloaded_at"],
                recorded_at=row["recorded_at"],
            )
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


# --- Preset Management ---


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
            return Preset(
                name=row["name"],
                params=json.loads(row["params"]),
                description=row["description"],
                created_at=row["created_at"],
                updated_at=row["updated_at"],
            )
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

        return [
            Preset(
                name=row["name"],
                params=json.loads(row["params"]),
                description=row["description"],
                created_at=row["created_at"],
                updated_at=row["updated_at"],
            )
            for row in rows
        ]


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
