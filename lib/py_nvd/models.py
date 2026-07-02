"""Retained data models for the NVD CLI and active pipeline helpers."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from pydantic import BaseModel, ConfigDict, Field, field_validator
from pydantic.dataclasses import dataclass


@dataclass(frozen=True)
class Taxon:
    """A taxon from the NCBI taxonomy."""

    tax_id: int
    scientific_name: str
    rank: str
    parent_tax_id: int | None = None


@dataclass(frozen=True)
class Preset:
    """A named collection of run parameters."""

    name: str
    params: dict[str, str | int | float | bool]
    created_at: str
    updated_at: str
    description: str | None = None

    @field_validator("params", mode="before")
    @classmethod
    def _parse_json_params(cls, v: str | dict) -> dict:
        """Parse JSON string from storage into a params dict."""
        if isinstance(v, str):
            return json.loads(v)
        return v


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
    """All parameters for an NVD pipeline run."""

    model_config = ConfigDict(extra="forbid")

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
    experimental: bool = Field(
        default=False,
        description="Enable experimental release-candidate features",
        json_schema_extra={"category": "Core"},
    )
    skip_assembly: bool = Field(
        default=False,
        description="Skip SPAdes assembly and all downstream contig classification",
        json_schema_extra={"category": "Core"},
    )
    skip_blast: bool = Field(
        default=False,
        description="Skip MEGABLAST and BLASTN contig search.",
        json_schema_extra={"category": "Core"},
    )

    blast_db_version: str | None = Field(
        None,
        description="BLAST database version",
        json_schema_extra={"category": "Databases"},
    )
    virus_index_version: str | None = Field(
        None,
        description="Virus enrichment index version",
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
    virus_index: Path | None = Field(
        None,
        description="Path to prebuilt vertebrate-infecting virus deacon index (.idx file)",
        json_schema_extra={"category": "Databases"},
    )
    virus_index_url: str | None = Field(
        None,
        description="URL to download a prebuilt vertebrate-infecting virus deacon index",
        json_schema_extra={"category": "Databases"},
    )
    virus_reference_fasta: Path | None = Field(
        None,
        description="Custom vertebrate-infecting virus FASTA for building an enrichment index",
        json_schema_extra={"category": "Databases"},
    )
    virus_kmer_size: int = Field(
        31,
        description="K-mer size for building a custom virus enrichment index",
        json_schema_extra={"category": "Databases"},
    )
    virus_window_size: int = Field(
        1,
        description="Minimizer window size for building a custom virus enrichment index",
        json_schema_extra={"category": "Databases"},
    )
    virus_abs_threshold: int = Field(
        1,
        description="Minimum absolute minimizer hits for virus read enrichment",
        json_schema_extra={"category": "Databases"},
    )
    virus_rel_threshold: float = Field(
        0.0,
        description="Minimum relative proportion of minimizers for virus read enrichment (0.0-1.0)",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_ref_path: Path | None = Field(
        None,
        description="Path to a prebuilt sourmash reference sketch database",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_ref_url: str | None = Field(
        None,
        description="URL to download a prebuilt sourmash reference sketch database",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_ref_fasta: Path | None = Field(
        None,
        description="Local FASTA to sketch as an experimental sourmash reference database",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_lineages_path: Path | None = Field(
        None,
        description="Path to a sourmash taxonomy lineages CSV matching the reference sketch database",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_lineages_url: str | None = Field(
        None,
        description="URL to download a sourmash taxonomy lineages CSV matching the reference sketch database",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_ksize: int = Field(
        31,
        description="K-mer size for experimental sourmash sketching",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_scaled: int = Field(
        50,
        description="Scaled value for experimental sourmash sketching",
        json_schema_extra={"category": "Databases"},
    )
    sourmash_threshold_bp: int = Field(
        50,
        description="Minimum estimated base-pair overlap for experimental sourmash gather",
        json_schema_extra={"category": "Databases"},
    )
    preprocess: bool = Field(
        default=False,
        description="Enable all preprocessing steps",
        json_schema_extra={"category": "Preprocessing"},
    )
    dedup: bool = Field(
        default=False,
        description="Deduplicate reads (umbrella: enables both dedup_seq and dedup_pos)",
        json_schema_extra={"category": "Preprocessing"},
    )
    dedup_seq: bool = Field(
        default=False,
        description="Sequence-based deduplication with clumpify (preprocessing)",
        json_schema_extra={"category": "Preprocessing"},
    )
    dedup_pos: bool = Field(
        default=False,
        description="Positional deduplication with samtools markdup (after alignment)",
        json_schema_extra={"category": "Preprocessing"},
    )
    trim_adapters: bool | None = Field(
        None,
        description="Trim Illumina adapters",
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
    host_index: Path | None = Field(
        None,
        description="Path to prebuilt host/contaminant index (.idx file)",
        json_schema_extra={"category": "Preprocessing"},
    )
    host_index_url: str | None = Field(
        None,
        description="URL to download a prebuilt host/contaminant index",
        json_schema_extra={"category": "Preprocessing"},
    )
    host_contaminants_fasta: Path | None = Field(
        None,
        description="Custom contaminant FASTA to build and union with other host indexes",
        json_schema_extra={"category": "Preprocessing"},
    )
    host_kmer_size: int = Field(
        31,
        description="K-mer size for building a custom host/contaminant index",
        json_schema_extra={"category": "Preprocessing"},
    )
    host_window_size: int = Field(
        15,
        description="Minimizer window size for building a custom host/contaminant index",
        json_schema_extra={"category": "Preprocessing"},
    )
    host_abs_threshold: int = Field(
        2,
        description="Minimum absolute minimizer hits to classify as contaminant",
        json_schema_extra={"category": "Preprocessing"},
    )
    host_rel_threshold: float = Field(
        0.01,
        description="Minimum relative proportion of minimizers (0.0-1.0)",
        json_schema_extra={"category": "Preprocessing"},
    )

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

    taxonomy_dir: Path | None = Field(
        None,
        description="Explicit taxonomy database directory",
        json_schema_extra={"category": "Core"},
    )
    taxonomy_mode: str | None = Field(
        None,
        description="Pipeline taxonomy mode: read_only or missing",
        json_schema_extra={"category": "Core"},
    )
    taxonomy_refresh: str = Field(
        "missing",
        description="Admin taxonomy refresh policy: missing, stale, or force",
        json_schema_extra={"category": "Core"},
    )
    taxonomy_max_age_days: int = Field(
        90,
        description="Freshness threshold for taxonomy refreshes and warnings",
        json_schema_extra={"category": "Core"},
    )

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
        "nrminor/nvd:v3.0.0",
        description="Container image",
        json_schema_extra={"category": "Internal"},
    )

    @field_validator(
        "cutoff_percent",
        "entropy",
        "tax_stringency",
        "host_rel_threshold",
        "virus_rel_threshold",
    )
    @classmethod
    def validate_zero_to_one(cls, v: float) -> float:
        """Validate that value is in the 0-1 range."""
        if not 0 <= v <= 1:
            msg = f"Must be between 0 and 1, got {v}"
            raise ValueError(msg)
        return v

    @field_validator(
        "max_blast_targets",
        "blast_retention_count",
        "min_consecutive_bases",
        "min_read_length",
        "max_concurrent_downloads",
        "host_kmer_size",
        "host_window_size",
        "host_abs_threshold",
        "virus_kmer_size",
        "virus_window_size",
        "virus_abs_threshold",
        "sourmash_ksize",
        "sourmash_scaled",
        "taxonomy_max_age_days",
    )
    @classmethod
    def validate_positive_int(cls, v: int) -> int:
        """Validate that value is a positive integer."""
        if v < 1:
            msg = f"Must be >= 1, got {v}"
            raise ValueError(msg)
        return v

    @field_validator("min_read_quality_illumina", "min_read_quality_nanopore")
    @classmethod
    def validate_non_negative_int(cls, v: int) -> int:
        """Validate that value is non-negative."""
        if v < 0:
            msg = f"Must be >= 0, got {v}"
            raise ValueError(msg)
        return v

    @field_validator("sourmash_threshold_bp")
    @classmethod
    def validate_non_negative_sourmash_threshold(cls, v: int) -> int:
        """Validate sourmash_threshold_bp is non-negative."""
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

    @field_validator("taxonomy_mode")
    @classmethod
    def validate_taxonomy_mode(cls, v: str | None) -> str | None:
        """Validate pipeline taxonomy mode."""
        if v is not None and v not in {"read_only", "missing"}:
            msg = "taxonomy_mode must be one of ['missing', 'read_only']"
            raise ValueError(msg)
        return v

    @field_validator("taxonomy_refresh")
    @classmethod
    def validate_taxonomy_refresh(cls, v: str) -> str:
        """Validate admin taxonomy refresh policy."""
        if v not in {"missing", "stale", "force"}:
            msg = "taxonomy_refresh must be one of ['force', 'missing', 'stale']"
            raise ValueError(msg)
        return v

    @field_validator(
        "sourmash_ref_path",
        "sourmash_ref_fasta",
        "sourmash_lineages_path",
        mode="before",
    )
    @classmethod
    def validate_not_url_path(cls, v: object) -> object:
        """Validate that local sourmash path params are not URLs."""
        if v is None:
            return v
        parsed = urlparse(str(v))
        if parsed.scheme in {"http", "https"}:
            msg = "Use the corresponding *_url param for remote sourmash resources"
            raise ValueError(msg)
        return v

    @field_validator("sourmash_ref_url", "sourmash_lineages_url")
    @classmethod
    def validate_url(cls, v: str | None) -> str | None:
        """Validate that sourmash URL params are HTTP(S) URLs."""
        if v is None:
            return v
        parsed = urlparse(v)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            msg = f"Expected an HTTP(S) URL, got {v}"
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

    @classmethod
    def merge(cls, *sources: dict[str, Any] | NvdParams | None) -> NvdParams:
        """Merge params from multiple sources with later sources taking precedence."""
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

    def to_nextflow_args(self, pipeline_root: Path) -> list[str]:
        """Convert params to a Nextflow command argument list."""
        cmd = ["nextflow", "run", str(pipeline_root)]

        for field_name, value in self.model_dump(exclude_none=True).items():
            if isinstance(value, bool):
                str_value = "true" if value else "false"
            elif isinstance(value, Path):
                str_value = str(value)
            elif isinstance(value, list):
                str_value = ",".join(str(item) for item in value)
            else:
                str_value = str(value)

            cmd.extend([f"--{field_name}", str_value])

        return cmd


@dataclass(frozen=True)
class ParamSource:
    """Tracks the source of a parameter value after merging."""

    field: str
    value: Any
    source: str
    overridden_value: Any | None = None
    overridden_source: str | None = None


@dataclass(frozen=True)
class TracedParams:
    """Result of trace_merge(): merged params with source tracking."""

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
            _, source_label = field_values[field]
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
    """Merge params and track where each value came from."""
    field_values: dict[str, tuple[Any, str]] = {}
    field_history: dict[str, list[tuple[Any, str]]] = {}

    for label, params in sources:
        if params is None:
            continue

        as_dict = (
            params.model_dump(exclude_none=True)
            if isinstance(params, NvdParams)
            else params
        )

        for field, value in as_dict.items():
            if value is None:
                continue

            if field in field_values:
                if field not in field_history:
                    field_history[field] = []
                field_history[field].append(field_values[field])

            field_values[field] = (value, label)

    source_dicts = [p for _, p in sources if p is not None]
    merged = NvdParams.merge(*source_dicts)
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
    """Get the category for an NvdParams field from its metadata."""
    field_info = NvdParams.model_fields.get(field_name)
    if field_info and field_info.json_schema_extra:
        extra = field_info.json_schema_extra
        if isinstance(extra, dict):
            return extra.get("category", "Internal")
    return "Internal"
