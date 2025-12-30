"""
Data models for py_nvd.

Pydantic dataclasses provide validation, serialization, and type safety.
All data crossing API boundaries uses these models.
"""

from typing import Literal

from pydantic.dataclasses import dataclass

# Type aliases for constrained values
Status = Literal["running", "completed", "failed"]
ProcessedSampleStatus = Literal["processing", "completed", "failed"]
DbType = Literal["blast", "stat", "gottcha2", "hostile"]
Platform = Literal["illumina", "ont", "sra"]
UploadType = Literal["blast", "blast_fasta", "gottcha2", "gottcha2_fasta"]


@dataclass(frozen=True)
class Run:
    """A pipeline run."""

    run_id: str
    sample_set_id: str  # SHA256[:16] of sorted sample IDs (idempotency key)
    started_at: str
    status: Status
    experiment_id: int | None = None  # Optional LabKey linkage
    completed_at: str | None = None


@dataclass(frozen=True)
class ProcessedSample:
    """A sample that has been processed with provenance tracking."""

    sample_id: str
    sample_set_id: str
    run_id: str
    processed_at: str
    status: ProcessedSampleStatus
    blast_db_version: str | None = None
    stat_db_version: str | None = None
    taxonomy_hash: str | None = None


@dataclass(frozen=True)
class Database:
    """A registered reference database."""

    db_type: DbType
    version: str
    path: str
    registered_at: str
    checksum: str | None = None


@dataclass(frozen=True)
class Upload:
    """Record of a sample upload to any target (LabKey, local, Globus, etc.)."""

    sample_id: str
    sample_set_id: str
    upload_type: UploadType
    upload_target: str  # 'labkey', 'local', 'globus', etc.
    content_hash: str  # SHA256 of uploaded data
    uploaded_at: str
    target_metadata: dict | None = (
        None  # Target-specific data (e.g., experiment_id, row_id)
    )


@dataclass(frozen=True)
class TaxonomyVersion:
    """Record of taxonomy database version used for a run."""

    run_id: str
    file_hash: str  # SHA256 of gettax.sqlite file
    downloaded_at: str
    recorded_at: str
    ncbi_rebuild_timestamp: str | None = None
    file_size: int | None = None


@dataclass(frozen=True)
class Taxon:
    """A taxon from the NCBI taxonomy."""

    tax_id: int
    scientific_name: str
    rank: str
    parent_tax_id: int | None = None


@dataclass(frozen=True)
class SraDownload:
    """A cached SRA download."""

    srr_accession: str
    download_path: str
    downloaded_at: str
    file_count: int | None = None
    total_bytes: int | None = None


@dataclass(frozen=True)
class Preset:
    """A named collection of run parameters."""

    name: str
    params: dict[str, str | int | float | bool]  # JSON-serializable param values
    created_at: str
    updated_at: str
    description: str | None = None
