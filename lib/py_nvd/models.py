"""
Data models for py_nvd.

Pydantic dataclasses provide validation, serialization, and type safety.
All data crossing API boundaries uses these models.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Self, TypeVar

from pydantic import BaseModel, Field, field_validator
from pydantic.dataclasses import dataclass

if TYPE_CHECKING:
    from sqlite3 import Row

logger = logging.getLogger(__name__)

T = TypeVar("T")


def _parse_iso8601(timestamp: str) -> datetime:
    """
    Parse an ISO8601 timestamp string to a timezone-aware datetime.

    Handles both 'Z' suffix and explicit timezone offsets. Naive timestamps
    (no timezone info) are assumed to be UTC.

    Args:
        timestamp: ISO8601 formatted timestamp string.

    Returns:
        Timezone-aware datetime in UTC.

    Raises:
        ValueError: If the timestamp cannot be parsed.
    """
    # Handle Z suffix (fromisoformat doesn't accept 'Z' directly)
    normalized = timestamp[:-1] + "+00:00" if timestamp.endswith("Z") else timestamp

    dt = datetime.fromisoformat(normalized)

    # Normalize naive datetimes to UTC
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=UTC)

    return dt


def _from_row(cls: type[T], row: Row) -> T:
    """
    Construct a dataclass instance from a sqlite3.Row.

    Extra columns in the row are ignored; missing columns raise
    an assertion error with diagnostic information.

    Args:
        cls: A dataclass type to instantiate
        row: A sqlite3.Row with column names matching field names

    Returns:
        An instance of cls populated from row data

    Raises:
        AssertionError: If cls is not a dataclass or row is missing required columns
        ValidationError: If Pydantic validation fails
    """
    assert dataclasses.is_dataclass(cls), f"_from_row requires a dataclass, got {cls.__name__}"

    field_names = [f.name for f in dataclasses.fields(cls)]
    row_keys = set(row.keys())

    missing = [f for f in field_names if f not in row_keys]
    assert not missing, (
        f"Row missing columns for {cls.__name__}: {missing}. Available: {sorted(row_keys)}"
    )

    return cls(**{f: row[f] for f in field_names})


# Type aliases for constrained values
Status = Literal["running", "completed", "failed"]
ProcessedSampleStatus = Literal["completed", "uploaded", "failed"]
DbType = Literal["blast", "stat", "gottcha2"]
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

    @property
    def duration_seconds(self) -> float | None:
        """
        Duration of the run in seconds.

        Returns None if:
        - completed_at is None (run still in progress or failed without completion)
        - Either timestamp cannot be parsed
        - Duration would be negative (indicates data corruption)

        Returns:
            Duration in seconds as a float, or None if unavailable.
        """
        if self.completed_at is None:
            return None

        try:
            start = _parse_iso8601(self.started_at)
            end = _parse_iso8601(self.completed_at)

            duration = (end - start).total_seconds()

            if duration < 0:
                logger.warning(
                    "Run %s has negative duration (completed_at before started_at): %ss",
                    self.run_id,
                    duration,
                )
                return None

            return duration  # return in try is clearer here
        except ValueError as e:
            logger.warning("Failed to parse timestamps for run %s: %s", self.run_id, e)
            return None

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a Run from a sqlite3.Row."""
        return _from_row(cls, row)


@dataclass(frozen=True)
class ProcessedSample:
    """
    A sample that has been processed with provenance tracking.

    Status values:
        completed: BLAST/REGISTER_HITS finished; results exist but not uploaded
        uploaded: Results uploaded to LabKey (terminal state)
        failed: Processing failed (retriable in future runs)

    Samples not in the database have not yet completed processing.
    Nextflow handles retry/resume for in-flight samples.
    """

    sample_id: str
    sample_set_id: str
    run_id: str
    processed_at: str
    status: ProcessedSampleStatus
    blast_db_version: str | None = None
    stat_db_version: str | None = None
    taxonomy_hash: str | None = None

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a ProcessedSample from a sqlite3.Row."""
        return _from_row(cls, row)


@dataclass(frozen=True)
class Database:
    """A registered reference database."""

    db_type: DbType
    version: str
    path: str
    registered_at: str
    checksum: str | None = None

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a Database from a sqlite3.Row."""
        return _from_row(cls, row)


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
    target_metadata: dict | None = None  # Target-specific data (e.g., experiment_id, row_id)

    @field_validator("target_metadata", mode="before")
    @classmethod
    def _parse_json_metadata(cls, v: str | dict | None) -> dict | None:
        """Parse JSON string from database into dict."""
        if isinstance(v, str):
            return json.loads(v)
        return v

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct an Upload from a sqlite3.Row."""
        return _from_row(cls, row)


@dataclass(frozen=True)
class TaxonomyVersion:
    """Record of taxonomy database version used for a run."""

    run_id: str
    file_hash: str  # SHA256 of gettax.sqlite file
    downloaded_at: str
    recorded_at: str
    ncbi_rebuild_timestamp: str | None = None
    file_size: int | None = None

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a TaxonomyVersion from a sqlite3.Row."""
        return _from_row(cls, row)


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
class SampleLock:
    """
    A processing lock on a sample to prevent duplicate work.

    Locks are acquired at run start and released on completion/upload.
    TTL-based expiration handles crashed runs. Machine fingerprint
    (hostname + username) enables conflict detection when the same
    run_id is resumed from different machines.
    """

    sample_id: str
    run_id: str
    hostname: str
    username: str
    locked_at: str  # ISO8601
    expires_at: str  # ISO8601

    @property
    def is_expired(self) -> bool:
        """Check if this lock has expired."""
        expires = _parse_iso8601(self.expires_at)
        return datetime.now(UTC) > expires

    @property
    def fingerprint(self) -> tuple[str, str]:
        """Return (hostname, username) tuple for ownership comparison."""
        return (self.hostname, self.username)

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a SampleLock from a sqlite3.Row."""
        return _from_row(cls, row)


@dataclass(frozen=True)
class Preset:
    """A named collection of run parameters."""

    name: str
    params: dict[str, str | int | float | bool]  # JSON-serializable param values
    created_at: str
    updated_at: str
    description: str | None = None

    @field_validator("params", mode="before")
    @classmethod
    def _parse_json_params(cls, v: str | dict) -> dict:
        """Parse JSON string from database into dict."""
        if isinstance(v, str):
            return json.loads(v)
        return v

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a Preset from a sqlite3.Row."""
        return _from_row(cls, row)


@dataclass(frozen=True)
class Hit:
    """
    A unique biological sequence identified by BLAST.

    The hit_key is a BLAKE3 hash of the canonical sequence (lexicographically
    smaller of the sequence and its reverse complement), making it strand-agnostic
    and deterministic across runs.

    Classification fields (top_taxid, top_taxid_name, top_taxid_rank) contain the
    taxonomic assignment from the most recent observation. When a hit has multiple
    observations, we use the most recent non-null classification. These may be None
    for hits from older runs that predate classification tracking.
    """

    hit_key: str  # BLAKE3 hash of canonical sequence (32 hex chars)
    sequence_length: int  # Original length in bp
    sequence_compressed: bytes  # 2-bit + zlib compressed (with N positions)
    gc_content: float  # GC fraction (0.0-1.0)
    first_seen_date: str  # ISO8601

    # Classification from most recent observation (may be None for old data)
    top_taxid: int | None = None
    top_taxid_name: str | None = None
    top_taxid_rank: str | None = None

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a Hit from a sqlite3.Row."""
        return _from_row(cls, row)


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

    @classmethod
    def from_row(cls, row: Row) -> Self:
        """Construct a HitObservation from a sqlite3.Row."""
        return _from_row(cls, row)


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

    Classification fields (top_taxid, top_taxid_name, top_taxid_rank) contain the
    most common taxonomic assignment across all observations. When a hit has
    multiple different classifications, we use the one that appears most frequently.
    These may be None for hits from older runs that predate classification tracking.
    """

    hit_key: str
    sequence_length: int
    gc_content: float
    first_seen_date: str  # ISO8601
    last_seen_date: str  # ISO8601, most recent observation
    sample_count: int  # Number of distinct samples
    run_count: int  # Number of distinct runs (sample_set_ids)
    observation_count: int  # Total observations

    # Most common classification across all observations (may be None for old data)
    top_taxid: int | None = None
    top_taxid_name: str | None = None
    top_taxid_rank: str | None = None


@dataclass(frozen=True)
class TimelineBucket:
    """
    Count of new hits discovered in a time period.

    Used for visualizing discovery rate over time (e.g., ASCII histograms).
    """

    period: str  # e.g., "2024-01" for month, "2024-W05" for week
    new_hits: int


@dataclass(frozen=True)
class MonthCompactionInfo:
    """
    Information about a single month's compaction.

    Used in CompactionResult to report what was compacted.
    """

    month: str  # YYYY-MM format
    observation_count: int  # Total rows compacted
    unique_hits: int  # Distinct hit_keys
    sample_sets: int  # Distinct sample_set_ids
    source_files: int  # Number of uncompacted files merged


@dataclass(frozen=True)
class CompactionResult:
    """
    Result of a hits compaction operation.

    Reports what was compacted, whether it was a dry run, and any
    months that were processed.
    """

    dry_run: bool
    months: list[MonthCompactionInfo]
    total_observations: int  # Sum of all observation_count
    total_source_files: int  # Sum of all source_files
    errors: list[str] = Field(default_factory=list)  # Any errors encountered


@dataclass(frozen=True)
class DeleteResult:
    """
    Result of deleting hits for a sample set.

    Tracks both uncompacted file deletion (cheap) and compacted month
    rewriting (more expensive). Used by delete_sample_set_hits() and
    reported by `nvd state prune`.
    """

    uncompacted_files_deleted: int
    compacted_months_rewritten: int


@dataclass(frozen=True)
class TaxonSummary:
    """
    Summary of a taxon's presence in a sample set.

    Used by get_top_taxa() to report the most prevalent taxa in a run.
    Ordered by sample_count (how many samples have this taxon) rather
    than hit_count (total contigs) since sample prevalence is more
    meaningful for interpreting findings.

    Note: These are raw BLAST/LCA results and require expert interpretation.
    Common false positives include endogenous retroviruses (ERVs) matching
    HIV/HTLV, and host sequences matching viral databases.
    """

    taxid: int | None  # NCBI taxonomy ID (None if unclassified)
    name: str  # Scientific name (e.g., "Human gammaherpesvirus 4")
    rank: str | None  # Taxonomic rank (e.g., "species", "genus", "no rank")
    sample_count: int  # Number of distinct samples with this taxon
    hit_count: int  # Number of distinct contigs (hit_keys) with this taxon
    avg_contig_length: int  # Mean contig length in bp


@dataclass(frozen=True)
class SampleSummary:
    """
    Summary of hits for a single sample within a run.

    Used by get_sample_summaries() to provide per-sample breakdown of findings.
    Useful for identifying which samples have the most interesting or unusual
    results, and for detailed investigation of specific samples.

    The taxa_detected list contains taxon names (excluding noise taxa like
    'root' and 'Homo sapiens') ordered alphabetically.
    """

    sample_id: str
    total_hits: int  # Total distinct contigs (hit_keys) in this sample
    viral_taxa_count: int  # Number of distinct taxa (excluding noise)
    taxa_detected: tuple[str, ...]  # Taxon names found (immutable for frozen dataclass)


@dataclass(frozen=True)
class ContigQuality:
    """
    Assembly quality metrics for contigs in a run.

    Used by get_contig_quality() to provide QC information about the
    assembled contigs. Useful for troubleshooting runs with poor assembly
    or understanding the characteristics of detected sequences.
    """

    total_contigs: int  # Number of distinct hit_keys
    total_bases: int  # Sum of all contig lengths
    mean_length: float  # Average contig length
    median_length: int  # Median contig length
    max_length: int  # Longest contig
    contigs_500bp_plus: int  # Contigs >= 500 bp
    contigs_1kb_plus: int  # Contigs >= 1000 bp
    mean_gc: float  # Average GC content (0.0-1.0)


@dataclass(frozen=True)
class CategorySummary:
    """
    Summary of findings grouped by virus family.

    Used by get_taxa_by_category() to provide a high-level triage view
    of findings. Categories correspond to virus families from the
    human_virus_families enrichment list (see FAMILY_SEARCH_PATTERNS
    in hits.py for the mapping).

    For specific taxon-level details, use get_top_taxa() instead.

    Note: Taxa that don't match any known family pattern are grouped
    under "Other". This may include novel viruses, misclassified taxa,
    or taxa with unexpected naming conventions.
    """

    category: str  # Virus family name (e.g., "Orthoherpesviridae") or "Other"
    sample_count: int  # Number of distinct samples with taxa in this category
    hit_count: int  # Number of distinct contigs (hit_keys) in this category
    taxa: tuple[str, ...]  # Specific taxon names in this category (alphabetized)


@dataclass(frozen=True)
class RunReport:
    """
    Complete summary of a pipeline run for reporting.

    Combines counts, quality metrics, and findings into a single object
    suitable for Slack notifications, CLI display, or JSON export.

    Used by get_run_report() to provide a comprehensive view of a run's
    results. This is the primary data structure for run-level reporting.

    Note: All findings require expert interpretation. See TaxonSummary
    and CategorySummary docstrings for caveats about false positives.
    """

    sample_set_id: str
    run_date: str | None  # ISO8601, earliest run_date in the sample set

    # Counts
    samples_analyzed: int  # Total distinct samples in the run
    samples_with_viral_hits: int  # Samples with at least one non-noise taxon
    samples_negative: int  # Samples with no viral hits (only noise taxa)
    unique_hits: int  # Distinct hit_keys (contigs)
    viral_taxa_found: int  # Distinct taxa (excluding noise)

    # Quality metrics
    median_contig_length: int  # Median sequence_length across unique hits
    contigs_over_500bp: int  # Count of contigs >= 500 bp

    # Highlights
    top_findings: tuple[TaxonSummary, ...]  # Top 5 taxa by sample prevalence
    sample_summaries: tuple[SampleSummary, ...]  # Per-sample breakdown


@dataclass(frozen=True)
class NovelTaxon:
    """
    A taxon seen for the first time in a given run.

    Used by get_novel_taxa() to identify taxa that have never appeared
    in any previous run. Useful for detecting emerging pathogens or
    new contamination sources.

    Note: "Novel" means not seen in the NVD state database, not necessarily
    novel to science. Historical data availability affects what's considered novel.
    """

    taxid: int | None  # NCBI taxonomy ID (None if unclassified)
    name: str  # Scientific name
    rank: str | None  # Taxonomic rank
    sample_count: int  # Samples in this run with this taxon
    hit_count: int  # Contigs in this run with this taxon


@dataclass(frozen=True)
class TaxonChange:
    """
    A taxon's prevalence change relative to historical baseline.

    Used by get_top_movers() to identify taxa with the largest deviation
    from their historical average prevalence. Positive change_pct indicates
    higher prevalence than usual; negative indicates lower.

    Note: Statistical significance depends on baseline_run_count. Small
    baselines produce noisy comparisons.
    """

    taxid: int | None  # NCBI taxonomy ID (None if unclassified)
    name: str  # Scientific name
    rank: str | None  # Taxonomic rank
    current_prevalence_pct: float  # Prevalence in this run (0-100)
    baseline_prevalence_pct: float  # Historical average prevalence (0-100)
    change_pct: float  # current - baseline
    current_sample_count: int  # Samples in this run with this taxon
    baseline_run_count: int  # Number of runs in baseline calculation


@dataclass(frozen=True)
class RareTaxon:
    """
    A taxon that appears infrequently in historical data.

    Used by get_rare_taxa() to flag taxa in the current run that have
    historically appeared in few runs. May indicate unusual findings
    worth investigation.

    Note: Rarity is relative to the available history. A taxon appearing
    in 5% of runs may be common in absolute terms but rare in this dataset.
    """

    taxid: int | None  # NCBI taxonomy ID (None if unclassified)
    name: str  # Scientific name
    rank: str | None  # Taxonomic rank
    historical_run_count: int  # Number of runs where this taxon appeared
    historical_run_pct: float  # Percentage of runs with this taxon (0-100)
    current_sample_count: int  # Samples in this run with this taxon


@dataclass(frozen=True)
class TaxonRunPresence:
    """
    A taxon's presence in a single run.

    Used as part of TaxonHistory to show the taxon's track record
    across runs. Each entry represents one run where the taxon was detected.
    """

    sample_set_id: str  # Run identifier
    run_date: str  # ISO8601 timestamp
    sample_count: int  # Samples in this run with the taxon
    prevalence_pct: float  # Percentage of samples in this run (0-100)


@dataclass(frozen=True)
class TaxonHistory:
    """
    Complete history of a taxon across all runs.

    Used by get_taxon_history() to provide a detailed track record
    for a specific taxon. Useful for understanding whether a finding
    is typical or unusual for this dataset.
    """

    taxid: int | None  # NCBI taxonomy ID (None if unclassified)
    name: str  # Scientific name
    rank: str | None  # Taxonomic rank
    total_runs_seen: int  # Number of runs where taxon appeared
    total_runs_in_system: int  # Total runs in the database
    overall_prevalence_pct: float  # Percentage of runs with this taxon (0-100)
    first_seen_date: str  # ISO8601, earliest run_date
    last_seen_date: str  # ISO8601, most recent run_date
    run_history: tuple[TaxonRunPresence, ...]  # Per-run details, ordered by date


@dataclass(frozen=True)
class RunComparison:
    """
    How a run compares to historical norms.

    Used by get_run_comparison() to contextualize a run's metrics
    against the historical distribution. Helps identify runs that
    are unusually large, small, or different in character.

    Note: Historical averages are simple means. With few runs in
    the system, these may not be representative.
    """

    sample_set_id: str  # Run identifier
    run_date: str | None  # ISO8601 timestamp

    # Current run metrics
    samples_analyzed: int  # Total samples in this run
    unique_hits: int  # Distinct contigs
    taxa_found: int  # Distinct taxa (excluding noise)
    median_contig_length: int  # Median sequence length

    # Historical baselines
    historical_run_count: int  # Number of prior runs
    avg_samples: float  # Historical average samples per run
    avg_hits: float  # Historical average hits per run
    avg_taxa: float  # Historical average taxa per run
    avg_median_length: float  # Historical average median contig length


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
    },
)

DEFAULT_HUMAN_VIRUS_FAMILIES = (
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
)


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
        default=False,
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
        description="Remove host reads with STAT (requires sra_human_db)",
        json_schema_extra={"category": "Preprocessing"},
    )
    sra_human_db: Path | None = Field(
        None,
        description="Path to human reads STAT database for host scrubbing and SRA submission prep",
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
        default=True,
        description="Include child taxa in analysis",
        json_schema_extra={"category": "Analysis"},
    )
    human_virus_families: list[str] = Field(
        default_factory=lambda: list(DEFAULT_HUMAN_VIRUS_FAMILIES),
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
        default=False,
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
    # Slack Notifications
    # =========================================================================
    slack_enabled: bool = Field(
        default=False,
        description="Enable Slack notifications for run completion",
        json_schema_extra={"category": "Notifications"},
    )
    slack_channel: str | None = Field(
        None,
        description="Slack channel ID for notifications",
        json_schema_extra={"category": "Notifications"},
    )

    # =========================================================================
    # Stateless Mode
    # =========================================================================
    stateless: bool = Field(
        default=False,
        description="Run without state management (disables run tracking, LabKey, Slack)",
        json_schema_extra={"category": "Core"},
    )
    taxonomy_dir: Path | None = Field(
        None,
        description="Taxonomy database directory (required when --stateless)",
        json_schema_extra={"category": "Core"},
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
        "nrminor/nvd:v2.4.0",
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
            msg = f"Invalid tools: {invalid}. Valid: {sorted(VALID_TOOLS)}"
            raise ValueError(msg)
        return v

    @field_validator("cutoff_percent", "entropy", "tax_stringency")
    @classmethod
    def validate_zero_to_one(cls, v: float) -> float:
        """Validate that value is in 0-1 range."""
        if not 0 <= v <= 1:
            msg = f"Must be between 0 and 1, got {v}"
            raise ValueError(msg)
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
            msg = f"Must be >= 1, got {v}"
            raise ValueError(msg)
        return v

    @field_validator("min_read_quality_illumina", "min_read_quality_nanopore")
    @classmethod
    def validate_non_negative_int(cls, v: int) -> int:
        """Validate that value is non-negative (>= 0)."""
        if v < 0:
            msg = f"Must be >= 0, got {v}"
            raise ValueError(msg)
        return v

    @field_validator("max_read_length")
    @classmethod
    def validate_max_read_length(cls, v: int | None) -> int | None:
        """Validate max_read_length is positive if set."""
        if v is not None and v < 1:
            msg = f"Must be >= 1, got {v}"
            raise ValueError(msg)
        return v

    @field_validator("slack_channel")
    @classmethod
    def validate_slack_channel(cls, v: str | None) -> str | None:
        """Validate Slack channel ID format."""
        if v is not None and not re.match(r"^C[A-Z0-9]+$", v):
            msg = (
                f"Invalid Slack channel ID: {v}. "
                "Must match pattern C[A-Z0-9]+ (e.g., 'C0123456789')"
            )
            raise ValueError(msg)
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
                as_dict = source.model_dump(exclude_none=True)
            else:
                as_dict = source
            merged.update({k: v for k, v in as_dict.items() if v is not None})
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


def _build_source_annotations(
    merged: NvdParams,
    field_values: dict[str, tuple[Any, str]],
    field_history: dict[str, list[tuple[Any, str]]],
) -> tuple[list[ParamSource], int]:
    """Build per-field source annotations from merge tracking data."""
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
            if field_history.get(field):
                overridden_value, overridden_source = field_history[field][-1]

            source_list.append(
                ParamSource(
                    field=field,
                    value=value,
                    source=source_label,
                    overridden_value=overridden_value,
                    overridden_source=overridden_source,
                ),
            )
        # Value is the default
        elif value is not None:
            defaults_used += 1
            source_list.append(
                ParamSource(
                    field=field,
                    value=value,
                    source="default",
                ),
            )

    return source_list, defaults_used


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
        str,
        list[tuple[Any, str]],
    ] = {}  # field -> [(value, source), ...]

    for label, params in sources:
        if params is None:
            continue

        as_dict = params.model_dump(exclude_none=True) if isinstance(params, NvdParams) else params

        for field, value in as_dict.items():
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
    source_list, defaults_used = _build_source_annotations(
        merged,
        field_values,
        field_history,
    )

    return TracedParams(
        params=merged,
        sources=source_list,
        defaults_used=defaults_used,
    )


# Canonical category order for display
PARAM_CATEGORIES = [
    "Core",
    "Databases",
    "Analysis",
    "Preprocessing",
    "LabKey",
    "Notifications",
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
