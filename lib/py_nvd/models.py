"""
Data models for py_nvd.

Pydantic dataclasses provide validation, serialization, and type safety.
All data crossing API boundaries uses these models.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator
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


@dataclass
class DatabaseResolution:
    """
    Result of resolving database versions from the registry.

    Used by resolve_database_versions() to return resolved versions,
    warnings about unregistered paths or version mismatches, and
    a list of registrations that were performed.

    Note: This is a mutable dataclass (not frozen) because the resolution
    process builds it incrementally.
    """

    # Resolved versions (None if path not provided or unresolvable)
    blast_db_version: str | None = None
    gottcha2_db_version: str | None = None
    stat_db_version: str | None = None

    # Warnings to display to the user
    warnings: list[str] = Field(default_factory=list)

    # Registrations performed: list of (db_type, version, path)
    auto_registered: list[tuple[DbType, str, str]] = Field(default_factory=list)


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


@dataclass(frozen=True)
class Hit:
    """
    A unique biological sequence identified by BLAST.

    The hit_key is a BLAKE3 hash of the canonical sequence (lexicographically
    smaller of the sequence and its reverse complement), making it strand-agnostic
    and deterministic across runs.
    """

    hit_key: str  # BLAKE3 hash of canonical sequence (32 hex chars)
    sequence_length: int  # Original length in bp
    sequence_compressed: bytes  # 2-bit + zlib compressed (with N positions)
    gc_content: float  # GC fraction (0.0-1.0)
    first_seen_date: str  # ISO8601


@dataclass(frozen=True)
class HitObservation:
    """
    Record of when/where a hit was observed.

    Links a Hit to a specific sample and run, enabling queries like
    "which samples have seen this hit?" and "what hits appeared in this run?"
    """

    hit_key: str  # References Hit.hit_key
    sample_set_id: str  # Links to Run.sample_set_id
    sample_id: str  # Source sample
    run_date: str  # ISO8601
    contig_id: str | None = None  # Original SPAdes contig ID (for traceability)


@dataclass(frozen=True)
class HitStats:
    """
    Summary statistics for the hits collection.

    Provides aggregate metrics for dashboard displays and health checks.
    All fields are None when the database is empty.
    """

    total_hits: int
    total_observations: int
    unique_samples: int
    unique_runs: int
    length_min: int | None
    length_max: int | None
    length_median: float | None
    gc_min: float | None
    gc_max: float | None
    gc_median: float | None
    date_first: str | None  # ISO8601, earliest first_seen_date
    date_last: str | None  # ISO8601, latest first_seen_date


@dataclass(frozen=True)
class RecurringHit:
    """
    A hit that appears in multiple samples or runs.

    Used for identifying sequences that recur across the dataset,
    which may indicate persistent infections, contamination, or
    lab reagent artifacts.
    """

    hit_key: str
    sequence_length: int
    gc_content: float
    first_seen_date: str  # ISO8601
    last_seen_date: str  # ISO8601, most recent observation
    sample_count: int  # Number of distinct samples
    run_count: int  # Number of distinct runs (sample_set_ids)
    observation_count: int  # Total observations


@dataclass(frozen=True)
class TimelineBucket:
    """
    Count of new hits discovered in a time period.

    Used for visualizing discovery rate over time (e.g., ASCII histograms).
    """

    period: str  # e.g., "2024-01" for month, "2024-W05" for week
    new_hits: int


# =============================================================================
# NvdParams: Pipeline Parameter Model
# =============================================================================

# Valid values for constrained fields
VALID_TOOLS = frozenset(
    {
        "stat_blast",
        "nvd",
        "stat",
        "blast",
        "stast",
        "gottcha",
        "all",
        "clumpify",
    }
)

DEFAULT_HUMAN_VIRUS_FAMILIES = [
    "Adenoviridae",
    "Anelloviridae",
    "Arenaviridae",
    "Arteriviridae",
    "Astroviridae",
    "Bornaviridae",
    "Peribunyaviridae",
    "Caliciviridae",
    "Coronaviridae",
    "Filoviridae",
    "Flaviviridae",
    "Hepadnaviridae",
    "Hepeviridae",
    "Orthoherpesviridae",
    "Orthomyxoviridae",
    "Papillomaviridae",
    "Paramyxoviridae",
    "Parvoviridae",
    "Picobirnaviridae",
    "Picornaviridae",
    "Pneumoviridae",
    "Polyomaviridae",
    "Poxviridae",
    "Sedoreoviridae",
    "Retroviridae",
    "Rhabdoviridae",
    "Togaviridae",
    "Kolmioviridae",
]


class NvdParams(BaseModel):
    """
    All parameters for an NVD pipeline run.

    This is the single source of truth for parameter handling. All params
    flow through this class, regardless of source (CLI, preset, params-file).
    Field definitions mirror the params block in nextflow.config.

    Note: Nextflow-native options (profile, config, resume) are NOT included
    here as they are not pipeline parameters—they control Nextflow itself.
    """

    # =========================================================================
    # Core Options
    # =========================================================================
    samplesheet: Path | None = Field(
        None,
        description="Path to samplesheet CSV",
        json_schema_extra={"category": "Core"},
    )
    results: Path | None = Field(
        None,
        description="Results directory",
        json_schema_extra={"category": "Core"},
    )
    tools: str | None = Field(
        None,
        description="Workflow(s) to run: all, stat_blast, gottcha, blast, clumpify",
        json_schema_extra={"category": "Core"},
    )
    experiment_id: str | None = Field(
        None,
        description="Experiment identifier (required for LabKey uploads)",
        json_schema_extra={"category": "Core"},
    )
    max_concurrent_downloads: int = Field(
        3,
        description="Maximum concurrent SRA downloads",
        json_schema_extra={"category": "Core"},
    )
    cleanup: bool | None = Field(
        None,
        description="Clean work directory after success",
        json_schema_extra={"category": "Core"},
    )
    work_dir: Path | None = Field(
        None,
        description="Nextflow work directory",
        json_schema_extra={"category": "Core"},
    )

    # =========================================================================
    # Database Versions
    # =========================================================================
    gottcha2_db_version: str | None = Field(
        None,
        description="GOTTCHA2 database version",
        json_schema_extra={"category": "Databases"},
    )
    blast_db_version: str | None = Field(
        None,
        description="BLAST database version",
        json_schema_extra={"category": "Databases"},
    )
    stat_db_version: str | None = Field(
        None,
        description="STAT database version",
        json_schema_extra={"category": "Databases"},
    )

    # =========================================================================
    # Database Paths
    # =========================================================================
    gottcha2_db: Path | None = Field(
        None,
        description="Path to GOTTCHA2 database directory",
        json_schema_extra={"category": "Databases"},
    )
    nvd_files: Path | None = Field(
        None,
        description="Path to NVD resource files directory",
        json_schema_extra={"category": "Databases"},
    )
    blast_db: Path | None = Field(
        None,
        description="Path to BLAST database directory",
        json_schema_extra={"category": "Databases"},
    )
    blast_db_prefix: str | None = Field(
        None,
        description="BLAST database name prefix",
        json_schema_extra={"category": "Databases"},
    )
    stat_index: Path | None = Field(
        None,
        description="Path to STAT index file",
        json_schema_extra={"category": "Databases"},
    )
    stat_dbss: Path | None = Field(
        None,
        description="Path to STAT DBSS file",
        json_schema_extra={"category": "Databases"},
    )
    stat_annotation: Path | None = Field(
        None,
        description="Path to STAT annotation file",
        json_schema_extra={"category": "Databases"},
    )
    human_virus_taxlist: Path | None = Field(
        None,
        description="Path to human virus taxonomy list file",
        json_schema_extra={"category": "Databases"},
    )

    # =========================================================================
    # Preprocessing
    # =========================================================================
    preprocess: bool = Field(
        False,
        description="Enable all preprocessing steps",
        json_schema_extra={"category": "Preprocessing"},
    )
    merge_pairs: bool | None = Field(
        None,
        description="Merge paired read mates based on overlaps",
        json_schema_extra={"category": "Preprocessing"},
    )
    dedup: bool | None = Field(
        None,
        description="Deduplicate reads",
        json_schema_extra={"category": "Preprocessing"},
    )
    trim_adapters: bool | None = Field(
        None,
        description="Trim Illumina adapters",
        json_schema_extra={"category": "Preprocessing"},
    )
    scrub_host_reads: bool | None = Field(
        None,
        description="Remove host reads with hostile",
        json_schema_extra={"category": "Preprocessing"},
    )
    hostile_index: Path | None = Field(
        None,
        description="Path to local hostile index (for offline use)",
        json_schema_extra={"category": "Preprocessing"},
    )
    hostile_index_name: str = Field(
        "human-t2t-hla.rs-viral-202401_ml-phage-202401",
        description="Standard hostile index name",
        json_schema_extra={"category": "Preprocessing"},
    )
    sra_human_db: Path | None = Field(
        None,
        description="Path to human reads database for SRA submission scrubbing",
        json_schema_extra={"category": "Preprocessing"},
    )
    human_read_scrub: Path | None = Field(
        None,
        description="DEPRECATED: Use sra_human_db instead",
        json_schema_extra={"category": "Preprocessing"},
    )
    filter_reads: bool | None = Field(
        None,
        description="Filter reads by quality/length",
        json_schema_extra={"category": "Preprocessing"},
    )
    min_read_quality_illumina: int = Field(
        20,
        description="Minimum average quality for Illumina reads",
        json_schema_extra={"category": "Preprocessing"},
    )
    min_read_quality_nanopore: int = Field(
        12,
        description="Minimum average quality for Nanopore reads",
        json_schema_extra={"category": "Preprocessing"},
    )
    min_read_length: int = Field(
        50,
        description="Minimum read length",
        json_schema_extra={"category": "Preprocessing"},
    )
    max_read_length: int | None = Field(
        None,
        description="Maximum read length (no limit if not specified)",
        json_schema_extra={"category": "Preprocessing"},
    )

    # =========================================================================
    # Analysis Parameters
    # =========================================================================
    cutoff_percent: float = Field(
        0.001,
        description="Minimum abundance threshold (0-1)",
        json_schema_extra={"category": "Analysis"},
    )
    entropy: float = Field(
        0.9,
        description="Entropy threshold for sequence complexity (0-1)",
        json_schema_extra={"category": "Analysis"},
    )
    min_consecutive_bases: int = Field(
        200,
        description="Minimum consecutive bases required",
        json_schema_extra={"category": "Analysis"},
    )
    qtrim: str = Field(
        "t",
        description="Quality trimming mode",
        json_schema_extra={"category": "Analysis"},
    )
    tax_stringency: float = Field(
        0.7,
        description="Taxonomy stringency (0-1)",
        json_schema_extra={"category": "Analysis"},
    )
    include_children: bool = Field(
        True,
        description="Include child taxa in analysis",
        json_schema_extra={"category": "Analysis"},
    )
    human_virus_families: list[str] = Field(
        default_factory=lambda: DEFAULT_HUMAN_VIRUS_FAMILIES.copy(),
        description="Virus family names to include in human virus analysis",
        json_schema_extra={"category": "Analysis"},
    )
    min_gottcha_reads: int = Field(
        250,
        description="Minimum reads for GOTTCHA2",
        json_schema_extra={"category": "Analysis"},
    )
    max_blast_targets: int = Field(
        100,
        description="Maximum BLAST targets to consider",
        json_schema_extra={"category": "Analysis"},
    )
    blast_retention_count: int = Field(
        5,
        description="Number of top BLAST hits to retain",
        json_schema_extra={"category": "Analysis"},
    )

    # =========================================================================
    # LabKey Integration
    # =========================================================================
    labkey: bool = Field(
        False,
        description="Enable LabKey integration",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_server: str | None = Field(
        None,
        description="LabKey server URL",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_project_name: str | None = Field(
        None,
        description="LabKey project name/path",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_webdav: str | None = Field(
        None,
        description="LabKey WebDAV URL for file uploads",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_schema: str | None = Field(
        None,
        description="LabKey database schema name",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_gottcha_fasta_list: str | None = Field(
        None,
        description="LabKey list name for GOTTCHA2 FASTA results",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_gottcha_full_list: str | None = Field(
        None,
        description="LabKey list name for full GOTTCHA2 results",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_gottcha_blast_verified_full_list: str | None = Field(
        None,
        description="LabKey list name for BLAST-verified GOTTCHA2 results",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_blast_meta_hits_list: str | None = Field(
        None,
        description="LabKey list name for BLAST metagenomic hits",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_blast_fasta_list: str | None = Field(
        None,
        description="LabKey list name for BLAST FASTA results",
        json_schema_extra={"category": "LabKey"},
    )
    labkey_exp_id_guard_list: str | None = Field(
        None,
        description="LabKey list name for experiment ID guard",
        json_schema_extra={"category": "LabKey"},
    )

    # =========================================================================
    # Internal Parameters (not user-configurable)
    # =========================================================================
    date: str | None = Field(
        None,
        description="Run date (computed by Nextflow)",
        json_schema_extra={"category": "Internal"},
    )
    resources: Path | None = Field(
        None,
        description="Directory for workflow resource files",
        json_schema_extra={"category": "Internal"},
    )
    refman_registry: Path | None = Field(
        None,
        description="Refman registry file path",
        json_schema_extra={"category": "Internal"},
    )
    monoimage: str = Field(
        "nrminor/nvd2:v0.1.0",
        description="Container image",
        json_schema_extra={"category": "Internal"},
    )
    state_dir: Path | None = Field(
        None,
        description="State directory for run tracking",
        json_schema_extra={"category": "Internal"},
    )

    # =========================================================================
    # Validators
    # =========================================================================

    @field_validator("tools")
    @classmethod
    def validate_tools(cls, v: str | None) -> str | None:
        """Validate tools option (supports comma-separated values)."""
        if v is None:
            return v
        tools_list = [t.strip() for t in v.split(",")]
        invalid = [t for t in tools_list if t not in VALID_TOOLS]
        if invalid:
            raise ValueError(f"Invalid tools: {invalid}. Valid: {sorted(VALID_TOOLS)}")
        return v

    @field_validator("cutoff_percent", "entropy", "tax_stringency")
    @classmethod
    def validate_zero_to_one(cls, v: float) -> float:
        """Validate that value is in 0-1 range."""
        if not 0 <= v <= 1:
            raise ValueError(f"Must be between 0 and 1, got {v}")
        return v

    @field_validator(
        "min_gottcha_reads",
        "max_blast_targets",
        "blast_retention_count",
        "min_consecutive_bases",
        "min_read_length",
        "max_concurrent_downloads",
    )
    @classmethod
    def validate_positive_int(cls, v: int) -> int:
        """Validate that value is a positive integer (>= 1)."""
        if v < 1:
            raise ValueError(f"Must be >= 1, got {v}")
        return v

    @field_validator("min_read_quality_illumina", "min_read_quality_nanopore")
    @classmethod
    def validate_non_negative_int(cls, v: int) -> int:
        """Validate that value is non-negative (>= 0)."""
        if v < 0:
            raise ValueError(f"Must be >= 0, got {v}")
        return v

    @field_validator("max_read_length")
    @classmethod
    def validate_max_read_length(cls, v: int | None) -> int | None:
        """Validate max_read_length is positive if set."""
        if v is not None and v < 1:
            raise ValueError(f"Must be >= 1, got {v}")
        return v

    # =========================================================================
    # Class Methods
    # =========================================================================

    @classmethod
    def merge(cls, *sources: dict[str, Any] | NvdParams | None) -> NvdParams:
        """
        Merge params from multiple sources with later sources taking precedence.

        Args:
            *sources: Dicts or NvdParams instances, in order of increasing
                      precedence. None sources are skipped. None values within
                      dicts are also skipped (they don't override).

        Returns:
            Validated NvdParams with merged values.

        Example:
            params = NvdParams.merge(
                preset_params,      # lowest precedence
                params_file_dict,   # middle
                cli_args,           # highest precedence
            )
        """
        merged: dict[str, Any] = {}
        for source in sources:
            if source is None:
                continue
            if isinstance(source, NvdParams):
                source = source.model_dump(exclude_none=True)
            for key, value in source.items():
                if value is not None:
                    merged[key] = value
        return cls.model_validate(merged)

    # =========================================================================
    # Instance Methods
    # =========================================================================

    def to_nextflow_args(self, pipeline_root: Path) -> list[str]:
        """
        Convert to Nextflow command-line arguments.

        Generates a command list suitable for subprocess.run(). All params
        are converted to --param-name value format. Booleans become
        "true"/"false" strings. Paths are converted to strings.

        Note: Nextflow-native options (profile, config, resume) must be
        added separately by the caller since they are not params.

        Args:
            pipeline_root: Path to the pipeline directory (for nextflow run)

        Returns:
            Command list starting with ["nextflow", "run", <pipeline_root>]
        """
        cmd = ["nextflow", "run", str(pipeline_root)]

        for field_name, value in self.model_dump(exclude_none=True).items():
            # Nextflow params use underscores (Groovy identifiers), so pass field names as-is.
            # The CLI uses kebab-case (--blast-db) per convention, but Typer handles that
            # conversion automatically. Here we're talking directly to Nextflow.
            param_name = field_name

            # Convert value to string
            if isinstance(value, bool):
                str_value = "true" if value else "false"
            elif isinstance(value, Path):
                str_value = str(value)
            elif isinstance(value, list):
                # Lists need special handling for Nextflow
                # Convert to JSON-like format that Nextflow can parse
                str_value = ",".join(str(item) for item in value)
            else:
                str_value = str(value)

            cmd.extend([f"--{param_name}", str_value])

        return cmd


# =============================================================================
# Param Source Tracking
# =============================================================================


@dataclass(frozen=True)
class ParamSource:
    """
    Tracks the source of a parameter value after merging.

    Used by trace_merge() to show where each value came from and
    what it may have overridden.
    """

    field: str
    value: Any
    source: str  # e.g., "cli", "file", "preset", "default"
    overridden_value: Any | None = None
    overridden_source: str | None = None


@dataclass(frozen=True)
class TracedParams:
    """
    Result of trace_merge(): merged params with source tracking.

    Attributes:
        params: The merged and validated NvdParams
        sources: Per-field source information for non-default values
        defaults_used: Count of fields using default values
    """

    params: NvdParams
    sources: list[ParamSource]
    defaults_used: int


def trace_merge(
    *sources: tuple[str, dict[str, Any] | NvdParams | None],
) -> TracedParams:
    """
    Merge params and track where each value came from.

    Like NvdParams.merge(), but returns source tracking information
    for display purposes (e.g., showing which preset value was overridden).

    Args:
        *sources: (label, params) pairs in precedence order (lowest first).
                  e.g., ("preset", preset_dict), ("file", file_dict)
                  None params are skipped.

    Returns:
        TracedParams with merged params and source annotations.

    Example:
        traced = trace_merge(
            ("preset", {"tools": "all", "cutoff_percent": 0.01}),
            ("file", {"tools": "blast", "samplesheet": Path("/data/s.csv")}),
        )
        # traced.params.tools == "blast"
        # traced.sources includes ParamSource(
        #     field="tools", value="blast", source="file",
        #     overridden_value="all", overridden_source="preset"
        # )
    """
    # Track the winning value and source for each field
    field_values: dict[str, tuple[Any, str]] = {}  # field -> (value, source_label)
    # Track what each field's value was before being overridden
    field_history: dict[
        str, list[tuple[Any, str]]
    ] = {}  # field -> [(value, source), ...]

    for label, params in sources:
        if params is None:
            continue

        if isinstance(params, NvdParams):
            params = params.model_dump(exclude_none=True)

        for field, value in params.items():
            if value is None:
                continue

            # Record history if this field was already set
            if field in field_values:
                if field not in field_history:
                    field_history[field] = []
                field_history[field].append(field_values[field])

            field_values[field] = (value, label)

    # Merge to get validated params (this also applies defaults)
    source_dicts = [p for _, p in sources if p is not None]
    merged = NvdParams.merge(*source_dicts)

    # Build source annotations
    source_list: list[ParamSource] = []
    defaults_used = 0

    for field in NvdParams.model_fields:
        value = getattr(merged, field)

        if field in field_values:
            # Value came from one of the sources
            _, source_label = field_values[field]

            # Check if it overrode something
            overridden_value = None
            overridden_source = None
            if field in field_history and field_history[field]:
                overridden_value, overridden_source = field_history[field][-1]

            source_list.append(
                ParamSource(
                    field=field,
                    value=value,
                    source=source_label,
                    overridden_value=overridden_value,
                    overridden_source=overridden_source,
                )
            )
        else:
            # Value is the default
            if value is not None:
                defaults_used += 1
                source_list.append(
                    ParamSource(
                        field=field,
                        value=value,
                        source="default",
                    )
                )

    return TracedParams(
        params=merged,
        sources=source_list,
        defaults_used=defaults_used,
    )


# =============================================================================
# Field Metadata Helpers
# =============================================================================

# Canonical category order for display
PARAM_CATEGORIES = [
    "Core",
    "Databases",
    "Analysis",
    "Preprocessing",
    "LabKey",
    "Internal",
]


def get_field_category(field_name: str) -> str:
    """
    Get the category for an NvdParams field from its metadata.

    Categories are defined in each field's json_schema_extra. This function
    provides a single source of truth for field categorization, used by
    CLI commands for grouping params in output.

    Args:
        field_name: Name of the NvdParams field

    Returns:
        Category name (e.g., "Core", "Analysis"), or "Internal" if not found.
    """
    field_info = NvdParams.model_fields.get(field_name)
    if field_info and field_info.json_schema_extra:
        extra = field_info.json_schema_extra
        # json_schema_extra can be a dict or callable
        if isinstance(extra, dict):
            return extra.get("category", "Internal")
    return "Internal"
