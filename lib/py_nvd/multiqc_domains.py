"""Typed orchestration of NVD report sections."""

from __future__ import annotations

from typing import Literal, TypeAlias

from pydantic import BaseModel, ConfigDict, NonNegativeInt, ValidationError

from py_nvd.multiqc_assembly import AssemblyRow, assembly_rows
from py_nvd.multiqc_fastx import (
    FASTX_MAX_POINTS,
    FastxRow,
    collect_fastx_profiles,
    fastx_row,
    histogram_series,
)
from py_nvd.multiqc_packages import DeaconReceipt, ReportPackage


class SectionModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)


class DeaconStats(SectionModel):
    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)

    seqs_in: NonNegativeInt
    seqs_out: NonNegativeInt
    seqs_removed: NonNegativeInt
    seqs_out_proportion: float
    seqs_removed_proportion: float
    bp_in: NonNegativeInt
    bp_out: NonNegativeInt
    bp_removed: NonNegativeInt
    bp_out_proportion: float
    bp_removed_proportion: float


class DeaconRow(SectionModel):
    sample_id: str
    query_class: str | None = None
    status: Literal["observed", "complete_empty", "invalid"]
    reads_in: NonNegativeInt | None = None
    reads_retained: NonNegativeInt | None = None
    reads_removed: NonNegativeInt | None = None
    read_retained_fraction: float | None = None
    read_removed_fraction: float | None = None
    bases_in: NonNegativeInt | None = None
    bases_retained: NonNegativeInt | None = None
    bases_removed: NonNegativeInt | None = None
    base_retained_fraction: float | None = None
    base_removed_fraction: float | None = None
    validation_message: str | None = None


class ConfiguredSkipRow(SectionModel):
    sample_id: str
    status: Literal["skipped"] = "skipped"
    reason: str


ReportRow: TypeAlias = DeaconRow | FastxRow | AssemblyRow | ConfiguredSkipRow


class TableSection(SectionModel):
    section_name: str
    description: str
    plot_type: Literal["table"] = "table"
    rows: tuple[ReportRow, ...]


class LineGraphSection(SectionModel):
    section_name: str
    description: str
    plot_type: Literal["linegraph"] = "linegraph"
    data: dict[str, dict[int, int]]
    xlab: str
    ylab: str


ReportSection: TypeAlias = TableSection | LineGraphSection


def report_row_label(row: ReportRow) -> str:
    parts = [row.sample_id]
    if isinstance(row, DeaconRow) and row.query_class is not None:
        parts.append(row.query_class)
    elif isinstance(row, FastxRow):
        for part in (row.stage, row.query_class, row.producer):
            if part is not None and part not in parts:
                parts.append(part)
    elif isinstance(row, AssemblyRow):
        parts.append(row.producer)
    return " | ".join(parts)


def deacon_row(package: ReportPackage) -> DeaconRow:
    receipt = package.receipt
    if not isinstance(receipt, DeaconReceipt):
        message = "DEACON parser received the wrong receipt"
        raise TypeError(message)
    try:
        stats = DeaconStats.model_validate_json(
            package.payload(receipt.stats).read_text(encoding="utf-8"),
            strict=True,
        )
    except (ValidationError, UnicodeDecodeError) as error:
        return DeaconRow(
            sample_id=receipt.sample_id,
            query_class=receipt.query_class,
            status="invalid",
            validation_message=str(error),
        )
    return DeaconRow(
        sample_id=receipt.sample_id,
        query_class=receipt.query_class,
        status="complete_empty"
        if stats.seqs_in == 0 and stats.bp_in == 0
        else "observed",
        reads_in=stats.seqs_in,
        reads_retained=stats.seqs_out,
        reads_removed=stats.seqs_removed,
        read_retained_fraction=stats.seqs_out_proportion,
        read_removed_fraction=stats.seqs_removed_proportion,
        bases_in=stats.bp_in,
        bases_retained=stats.bp_out,
        bases_removed=stats.bp_removed,
        base_retained_fraction=stats.bp_out_proportion,
        base_removed_fraction=stats.bp_removed_proportion,
    )


def configured_skip_rows(
    sample_ids: tuple[str, ...],
    label: str,
) -> tuple[ConfiguredSkipRow, ...]:
    return tuple(
        ConfiguredSkipRow(
            sample_id=sample_id,
            reason=f"{label} disabled by configuration",
        )
        for sample_id in sample_ids
    )


def build_domain_sections(
    *,
    packages: tuple[ReportPackage, ...],
    sample_ids: tuple[str, ...],
    target_enrichment_enabled: bool,
    depletion_enabled: bool,
    assembly_enabled: bool,
) -> dict[str, ReportSection]:
    sections: dict[str, ReportSection] = {}
    for domain, enabled, section_id, name, description in (
        (
            "target_enrichment",
            target_enrichment_enabled,
            "nvd_target_enrichment",
            "Target Enrichment",
            "Target enrichment reads and bases retained or removed.",
        ),
        (
            "depletion",
            depletion_enabled,
            "nvd_depletion",
            "Depletion",
            "Depletion reads and bases retained or removed.",
        ),
    ):
        rows: tuple[ReportRow, ...] = (
            tuple(
                deacon_row(package)
                for package in packages
                if isinstance(package.receipt, DeaconReceipt)
                and package.receipt.domain == domain
            )
            if enabled
            else configured_skip_rows(sample_ids, name.lower())
        )
        if rows:
            sections[section_id] = TableSection(
                section_name=name,
                description=description,
                rows=rows,
            )

    fastx = collect_fastx_profiles(packages, domain="fastx")
    if fastx:
        sections["nvd_fastx_profiles"] = TableSection(
            section_name="Query Sequence Profiles",
            description="Processed-read and filtered-contig FASTX QC profiles.",
            rows=tuple(fastx_row(item) for item in fastx),
        )
    for quality, section_id, name, xlab in (
        (
            False,
            "nvd_fastx_length_distribution",
            "Query Sequence Length Distribution",
            "Length bin start",
        ),
        (
            True,
            "nvd_fastx_quality_distribution",
            "Read Quality Distribution",
            "Mean quality bin start",
        ),
    ):
        series = histogram_series(fastx, quality=quality)
        if series:
            sections[section_id] = LineGraphSection(
                section_name=name,
                description=f"Producer histograms rebinned to at most {FASTX_MAX_POINTS} presentation points.",
                data=series,
                xlab=xlab,
                ylab="Sequence count",
            )

    assembly: tuple[ReportRow, ...] = (
        assembly_rows(
            packages,
            sample_order={
                sample_id: index for index, sample_id in enumerate(sample_ids)
            },
        )
        if assembly_enabled
        else configured_skip_rows(sample_ids, "assembly")
    )
    if assembly:
        sections["nvd_assembly"] = TableSection(
            section_name="Assembly Decision Metrics",
            description="Assembly eligibility, producer profiles, and union contribution.",
            rows=assembly,
        )
    return sections
