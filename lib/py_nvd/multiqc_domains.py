"""Typed orchestration of NVD report sections."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    PositiveInt,
    ValidationError,
)

from py_nvd.multiqc_assembly import AssemblyRow, assembly_rows
from py_nvd.multiqc_blast import (
    BlastDisabledRow,
    BlastTaskRow,
    blast_task_rows,
    collect_megablast_query_partitions,
)
from py_nvd.multiqc_fastx import (
    FASTX_MAX_POINTS,
    FastxRow,
    collect_fastx_profiles,
    fastx_row,
    histogram_series,
)
from py_nvd.multiqc_packages import DeaconReceipt, ReportPackage
from py_nvd.multiqc_query_preparation import (
    PreparedQueryRow,
    collect_prepared_query_summaries,
    prepared_query_rows,
)
from py_nvd.multiqc_taxonomy import (
    FAINTEST_SIGNALS_MAX_ROWS_PER_SAMPLE,
    STRONGEST_SIGNALS_MAX_ROWS_PER_SAMPLE,
    InvalidTaxonFindingRow,
    TaxonFindingRow,
    collect_taxon_big_tables,
    faintest_taxon_rows,
    higher_risk_taxon_rows,
    strongest_taxon_rows,
)


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


ReportRow: TypeAlias = (
    DeaconRow
    | FastxRow
    | AssemblyRow
    | ConfiguredSkipRow
    | PreparedQueryRow
    | BlastTaskRow
    | BlastDisabledRow
    | TaxonFindingRow
    | InvalidTaxonFindingRow
)


class TableColumn(SectionModel):
    title: str
    description: str | None = None
    hidden: bool = False
    format: str | None = None


class TableSection(SectionModel):
    section_name: str
    description: str
    plot_type: Literal["table"] = "table"
    rows: tuple[ReportRow, ...]
    headers: dict[str, TableColumn] = Field(default_factory=dict)
    first_column_title: str | None = None


class LineGraphSection(SectionModel):
    section_name: str
    description: str
    plot_type: Literal["linegraph"] = "linegraph"
    data: dict[str, dict[int, int]]
    xlab: str
    ylab: str
    smooth_points: PositiveInt


ReportSection: TypeAlias = TableSection | LineGraphSection


@dataclass(frozen=True, slots=True)
class DomainConfiguration:
    experimental_enabled: bool
    target_enrichment_enabled: bool
    depletion_enabled: bool
    assembly_enabled: bool
    blast_enabled: bool


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
    elif (
        isinstance(row, (PreparedQueryRow, BlastTaskRow))
        and row.query_class_code is not None
    ):
        parts.append(row.query_class_code)
    elif isinstance(row, (TaxonFindingRow, InvalidTaxonFindingRow)):
        return f"{row.sample_id} · {row.taxon_key}"
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


def blast_task_breakdown_section(
    packages: tuple[ReportPackage, ...],
    *,
    sample_ids: tuple[str, ...],
    blast_enabled: bool,
) -> TableSection | None:
    description = (
        "Queries are searched first with MEGABLAST. Queries without a MEGABLAST "
        "match are forwarded to BLASTN."
    )
    if not blast_enabled:
        return TableSection(
            section_name="BLAST Task Breakdown",
            description=description,
            rows=(BlastDisabledRow(),),
            headers={"notice": TableColumn(title="Run configuration")},
        )

    partitions = collect_megablast_query_partitions(packages)
    if not partitions:
        return None
    return TableSection(
        section_name="BLAST Task Breakdown",
        description=description,
        rows=blast_task_rows(
            partitions,
            sample_order={
                sample_id: index for index, sample_id in enumerate(sample_ids)
            },
        ),
        headers={
            "query_class": TableColumn(title="Query class"),
            "queries_searched": TableColumn(title="Queries searched"),
            "matched_by_megablast": TableColumn(title="Matched by MEGABLAST"),
            "forwarded_to_blastn": TableColumn(title="Forwarded to BLASTN"),
            "forwarded_to_blastn_percentage": TableColumn(
                title="Forwarded to BLASTN (%)",
                format="{:,.1f}",
            ),
            "problem": TableColumn(title="Problem"),
            "validation_details": TableColumn(
                title="Validation details",
                hidden=True,
            ),
        },
    )


def taxon_table_headers() -> dict[str, TableColumn]:
    return {
        "taxon_name": TableColumn(title="Taxon"),
        "taxon_rank": TableColumn(title="Rank"),
        "support_tier": TableColumn(title="Support"),
        "taxon_crumbs": TableColumn(title="Taxon CRUMBS"),
        "relative_crumbs_percent": TableColumn(
            title="CRUMBS share (%)",
            format="{:,.2f}",
        ),
        "supporting_query_count": TableColumn(title="Supporting queries"),
        "taxid": TableColumn(title="TaxID"),
        "who_risk_group": TableColumn(title="WHO risk group"),
        "total_query_span": TableColumn(title="Total query span"),
        "total_crumbs_score": TableColumn(title="Total CRUMBS score"),
        "strong_query_count": TableColumn(title="Strong queries"),
        "moderate_query_count": TableColumn(title="Moderate queries"),
        "weak_query_count": TableColumn(title="Weak queries"),
        "review_query_count": TableColumn(title="Review queries"),
        "redacted_query_count": TableColumn(title="Redacted queries"),
        "supporting_genome_like_contig_count": TableColumn(
            title="Genome-like contigs",
        ),
        "supporting_long_contig_count": TableColumn(title="Long contigs"),
        "supporting_short_contig_count": TableColumn(title="Short contigs"),
        "supporting_merged_pair_count": TableColumn(title="Merged pairs"),
        "supporting_single_read_count": TableColumn(title="Single reads"),
        "support_tier_rule": TableColumn(title="Support rule"),
        "support_note": TableColumn(title="Support note"),
        "problem": TableColumn(title="Problem"),
        "validation_details": TableColumn(
            title="Validation details",
            hidden=True,
        ),
    }


def taxon_projection_section(
    *,
    section_name: str,
    description: str,
    rows: tuple[TaxonFindingRow | InvalidTaxonFindingRow, ...],
) -> TableSection | None:
    if not rows:
        return None
    return TableSection(
        section_name=section_name,
        description=description,
        rows=rows,
        headers=taxon_table_headers(),
        first_column_title="Finding",
    )


def higher_risk_taxonomic_findings_section(
    packages: tuple[ReportPackage, ...],
) -> TableSection | None:
    section = taxonomy_sections(packages).get("nvd_higher_risk_taxonomic_findings")
    return section if isinstance(section, TableSection) else None


def taxonomy_sections(
    packages: tuple[ReportPackage, ...],
) -> dict[str, ReportSection]:
    tables = collect_taxon_big_tables(packages)
    if not tables:
        return {}

    higher_risk_rows, higher_risk_eligible = higher_risk_taxon_rows(tables)
    strongest_rows, strongest_eligible = strongest_taxon_rows(tables)
    faintest_rows, faintest_eligible = faintest_taxon_rows(tables)
    big_table_guidance = (
        "All taxon signals are detailed in the Taxon Big Table. Queries supporting "
        "each taxon are detailed in the Query Big Table."
    )
    projections = (
        (
            "nvd_higher_risk_taxonomic_findings",
            "Higher-Risk Taxonomic Findings",
            (
                f"Showing {sum(isinstance(row, TaxonFindingRow) for row in higher_risk_rows)} "
                f"of {higher_risk_eligible} WHO risk group 2-4 findings. Rows may "
                f"also appear in the signal-strength views. {big_table_guidance}"
            ),
            higher_risk_rows,
        ),
        (
            "nvd_strongest_taxonomic_signals",
            "Strongest Taxonomic Signals",
            (
                f"Showing {sum(isinstance(row, TaxonFindingRow) for row in strongest_rows)} "
                f"of {strongest_eligible} taxon signals with Taxon CRUMBS greater "
                f"than zero, limited to the {STRONGEST_SIGNALS_MAX_ROWS_PER_SAMPLE} "
                "strongest per sample. "
                "Taxon CRUMBS is normalized support per assigned query base, not "
                f"genome-copy abundance. Rows may appear in multiple views. {big_table_guidance}"
            ),
            strongest_rows,
        ),
        (
            "nvd_faintest_taxonomic_signals",
            "Faintest Taxonomic Signals",
            (
                f"Showing {sum(isinstance(row, TaxonFindingRow) for row in faintest_rows)} "
                f"of {faintest_eligible} taxon signals with Taxon CRUMBS greater "
                f"than zero, limited to the {FAINTEST_SIGNALS_MAX_ROWS_PER_SAMPLE} "
                "faintest per sample. Taxon "
                "CRUMBS is normalized support per assigned query base, not genome-copy "
                f"abundance. Rows may appear in multiple views. {big_table_guidance}"
            ),
            faintest_rows,
        ),
    )
    sections: dict[str, ReportSection] = {}
    for section_id, section_name, description, rows in projections:
        section = taxon_projection_section(
            section_name=section_name,
            description=description,
            rows=rows,
        )
        if section is not None:
            sections[section_id] = section
    return sections


def build_domain_sections(
    *,
    packages: tuple[ReportPackage, ...],
    sample_ids: tuple[str, ...],
    sample_platforms: dict[str, str],
    configuration: DomainConfiguration,
) -> dict[str, ReportSection]:
    sections = taxonomy_sections(packages)
    for domain, enabled, section_id, name, description in (
        (
            "target_enrichment",
            configuration.target_enrichment_enabled,
            "nvd_target_enrichment",
            "Target Enrichment",
            "Target enrichment reads and bases retained or removed.",
        ),
        (
            "depletion",
            configuration.depletion_enabled,
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
    for stage, section_id, name, subject in (
        (
            "single_read",
            "nvd_fastx_single_read_length_distribution",
            "Single Read Length Distribution",
            "single reads",
        ),
        (
            "overlap_merged_pair",
            "nvd_fastx_overlap_merged_pair_length_distribution",
            "Merged Pair Length Distribution",
            "overlap-merged pairs",
        ),
        (
            "filtered_contigs",
            "nvd_fastx_filtered_contigs_length_distribution",
            "Contig Length Distribution",
            "filtered contigs",
        ),
    ):
        series = histogram_series(fastx, quality=False, stage=stage)
        if series:
            sections[section_id] = LineGraphSection(
                section_name=name,
                description=f"Producer length histograms for {subject}, rebinned to at most {FASTX_MAX_POINTS} presentation points.",
                data=series,
                xlab="Length bin start",
                ylab="Sequence count",
                smooth_points=FASTX_MAX_POINTS,
            )

    quality_series = histogram_series(fastx, quality=True)
    if quality_series:
        sections["nvd_fastx_quality_distribution"] = LineGraphSection(
            section_name="Read Quality Distribution",
            description=f"Producer quality histograms rebinned to at most {FASTX_MAX_POINTS} presentation points.",
            data=quality_series,
            xlab="Mean quality bin start",
            ylab="Sequence count",
            smooth_points=FASTX_MAX_POINTS,
        )

    assembly: tuple[ReportRow, ...] = (
        assembly_rows(
            packages,
            sample_order={
                sample_id: index for index, sample_id in enumerate(sample_ids)
            },
        )
        if configuration.assembly_enabled
        else configured_skip_rows(sample_ids, "assembly")
    )
    if assembly:
        sections["nvd_assembly"] = TableSection(
            section_name="Assembly Decision Metrics",
            description="Assembly eligibility, producer profiles, and union contribution.",
            rows=assembly,
        )

    prepared_query_summaries = collect_prepared_query_summaries(packages)
    if prepared_query_summaries:
        rows = prepared_query_rows(
            prepared_query_summaries,
            sample_order={
                sample_id: index for index, sample_id in enumerate(sample_ids)
            },
            sample_platforms=sample_platforms,
            assembly_enabled=configuration.assembly_enabled,
            unassembled_read_querying_enabled=configuration.experimental_enabled,
        )
        sections["nvd_prepared_blast_query_batches"] = TableSection(
            section_name="Prepared BLAST Query Batches",
            description=(
                "Prepared query batches observed for each sample and query class. "
                "BLAST querying of unassembled reads alongside assembly contigs is "
                "currently only available in experimental mode. Resume the run if "
                "a summary is invalid; if the problem recurs, inspect the "
                "SUMMARIZE_BLAST_QUERY_BATCHES process and include its logs when "
                "reporting the issue."
            ),
            rows=rows,
            headers={
                "query_class": TableColumn(title="Query class"),
                "sequence_type": TableColumn(title="Sequence type"),
                "availability": TableColumn(title="Availability"),
                "query_sequences": TableColumn(
                    title="Prepared queries",
                    description="Number of sequences in the prepared BLAST query batch.",
                ),
                "platform": TableColumn(title="Platform", hidden=True),
                "problem": TableColumn(title="Problem"),
                "validation_details": TableColumn(
                    title="Validation details",
                    hidden=True,
                ),
            },
        )

    blast_section = blast_task_breakdown_section(
        packages,
        sample_ids=sample_ids,
        blast_enabled=configuration.blast_enabled,
    )
    if blast_section is not None:
        sections["nvd_blast_task_breakdown"] = blast_section
    return sections
