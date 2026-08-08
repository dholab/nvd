"""Typed Taxon Big Table parsing and bounded findings presentation."""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from typing import Annotated, Literal, TypeAlias

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeFloat,
    NonNegativeInt,
    StringConstraints,
    ValidationError,
)

from py_nvd.multiqc_packages import ReportPackage, TaxonBigTableReceipt

NonEmptyString = Annotated[str, StringConstraints(min_length=1)]
SupportTier: TypeAlias = Literal["strong", "moderate", "review", "weak", "redacted"]
RiskGroup: TypeAlias = Literal["RG1", "RG2", "RG3", "RG4"]
HIGHER_RISK_GROUPS = frozenset(("RG2", "RG3", "RG4"))
HIGHER_RISK_FINDINGS_MAX_ROWS = 100
STRONGEST_SIGNALS_MAX_ROWS_PER_SAMPLE = 10
FAINTEST_SIGNALS_MAX_ROWS_PER_SAMPLE = 5
SUPPORT_TIER_ORDER = {
    "strong": 0,
    "moderate": 1,
    "review": 2,
    "weak": 3,
    "redacted": 4,
}
RISK_GROUP_ORDER = {"RG4": 0, "RG3": 1, "RG2": 2}


class TaxonProducerRow(BaseModel):
    """One typed row from a producer-owned Taxon Big Table."""

    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)

    sample_id: NonEmptyString
    taxon_name: NonEmptyString
    taxon_rank: NonEmptyString
    support_tier: SupportTier
    taxon_crumbs: NonNegativeFloat | None = None
    relative_crumbs_percent: NonNegativeFloat | None = None
    supporting_query_count: NonNegativeInt
    taxid: NonNegativeInt
    who_risk_group: RiskGroup | None = None
    total_query_span: NonNegativeInt
    total_crumbs_score: NonNegativeFloat | None = None
    strong_query_count: NonNegativeInt
    moderate_query_count: NonNegativeInt
    weak_query_count: NonNegativeInt
    review_query_count: NonNegativeInt
    redacted_query_count: NonNegativeInt
    supporting_genome_like_contig_count: NonNegativeInt
    supporting_long_contig_count: NonNegativeInt
    supporting_short_contig_count: NonNegativeInt
    supporting_merged_pair_count: NonNegativeInt
    supporting_single_read_count: NonNegativeInt
    support_tier_rule: NonEmptyString
    support_note: NonEmptyString


class TaxonFindingRow(TaxonProducerRow):
    """One taxonomic finding presented in a bounded MultiQC projection."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    taxon_key: str = Field(exclude=True)
    problem: str | None = None
    validation_details: str | None = None


class InvalidTaxonFindingRow(BaseModel):
    """One Taxon Big Table package that could not be read."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    sample_id: str
    taxon_key: Literal["invalid"] = Field(default="invalid", exclude=True)
    problem: Literal["Taxonomic findings could not be read"] = (
        "Taxonomic findings could not be read"
    )
    validation_details: str


class TaxonBigTableError(ValueError):
    """A Taxon Big Table violates its report-facing contract."""


TaxonBigTableParseError: TypeAlias = (
    ValidationError | UnicodeDecodeError | csv.Error | TaxonBigTableError
)


@dataclass(frozen=True, slots=True)
class ParsedTaxonBigTable:
    receipt: TaxonBigTableReceipt
    rows: tuple[TaxonProducerRow, ...]


@dataclass(frozen=True, slots=True)
class InvalidTaxonBigTable:
    receipt: TaxonBigTableReceipt
    error: TaxonBigTableParseError


TaxonBigTableInput: TypeAlias = ParsedTaxonBigTable | InvalidTaxonBigTable
TaxonReportRow: TypeAlias = TaxonFindingRow | InvalidTaxonFindingRow


def parse_taxon_big_table_package(package: ReportPackage) -> ParsedTaxonBigTable:
    receipt = package.receipt
    if not isinstance(receipt, TaxonBigTableReceipt):
        message = "Taxon Big Table parser received the wrong receipt"
        raise TypeError(message)

    with package.payload(receipt.table).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = set(TaxonProducerRow.model_fields)
        if reader.fieldnames is None or not required_columns.issubset(
            reader.fieldnames,
        ):
            message = "Taxon Big Table is missing required columns"
            raise TaxonBigTableError(message)
        payload_rows = tuple(reader)

    parsed_rows: list[TaxonProducerRow] = []
    seen_taxa: set[tuple[int, str, str]] = set()
    for payload_row in payload_rows:
        normalized = dict(payload_row)
        for field_name in (
            "taxon_crumbs",
            "relative_crumbs_percent",
            "total_crumbs_score",
            "who_risk_group",
        ):
            if normalized[field_name] == "":
                del normalized[field_name]
        row = TaxonProducerRow.model_validate_strings(normalized, strict=True)
        if row.sample_id != receipt.sample_id:
            message = "Taxon Big Table sample does not match its receipt"
            raise TaxonBigTableError(message)
        taxon_key = (row.taxid, row.taxon_name, row.taxon_rank)
        if taxon_key in seen_taxa:
            message = "Taxon Big Table contains a duplicate assigned taxon"
            raise TaxonBigTableError(message)
        seen_taxa.add(taxon_key)
        parsed_rows.append(row)
    return ParsedTaxonBigTable(receipt=receipt, rows=tuple(parsed_rows))


def collect_taxon_big_tables(
    packages: tuple[ReportPackage, ...],
) -> tuple[TaxonBigTableInput, ...]:
    tables: list[TaxonBigTableInput] = []
    for package in packages:
        if not isinstance(package.receipt, TaxonBigTableReceipt):
            continue
        try:
            tables.append(parse_taxon_big_table_package(package))
        except (
            ValidationError,
            UnicodeDecodeError,
            csv.Error,
            TaxonBigTableError,
        ) as error:
            tables.append(InvalidTaxonBigTable(receipt=package.receipt, error=error))
    return tuple(tables)


def higher_risk_sort_key(
    row: TaxonProducerRow,
) -> tuple[int, int, float, float, str, str, int]:
    risk_group = row.who_risk_group
    if risk_group not in RISK_GROUP_ORDER:
        message = "Taxon finding has no higher-risk group"
        raise ValueError(message)
    return (
        RISK_GROUP_ORDER[risk_group],
        SUPPORT_TIER_ORDER[row.support_tier],
        -row.taxon_crumbs if row.taxon_crumbs is not None else float("inf"),
        (
            -row.total_crumbs_score
            if row.total_crumbs_score is not None
            else float("inf")
        ),
        row.sample_id,
        row.taxon_name,
        row.taxid,
    )


def higher_risk_taxon_rows(
    tables: tuple[TaxonBigTableInput, ...],
    *,
    max_rows: int = HIGHER_RISK_FINDINGS_MAX_ROWS,
) -> tuple[tuple[TaxonReportRow, ...], int]:
    if max_rows < 1:
        message = "Higher-risk findings row budget must be positive"
        raise ValueError(message)

    eligible: list[TaxonProducerRow] = []
    invalid: list[InvalidTaxonFindingRow] = []
    for table in tables:
        if isinstance(table, InvalidTaxonBigTable):
            invalid.append(
                InvalidTaxonFindingRow(
                    sample_id=table.receipt.sample_id,
                    validation_details=str(table.error),
                ),
            )
            continue
        eligible.extend(
            row for row in table.rows if row.who_risk_group in HIGHER_RISK_GROUPS
        )

    eligible.sort(key=higher_risk_sort_key)
    selected = finding_rows(eligible[:max_rows])
    return ((*selected, *sorted(invalid, key=lambda row: row.sample_id)), len(eligible))


def finding_rows(rows: list[TaxonProducerRow]) -> tuple[TaxonFindingRow, ...]:
    label_counts = Counter((row.sample_id, row.taxon_name) for row in rows)
    return tuple(
        TaxonFindingRow(
            **row.model_dump(),
            taxon_key=(
                row.taxon_name
                if label_counts[(row.sample_id, row.taxon_name)] == 1
                else f"{row.taxon_name} (TaxID {row.taxid})"
            ),
        )
        for row in rows
    )


def signal_sort_key(
    row: TaxonProducerRow,
    *,
    strongest_first: bool,
) -> tuple[float, float, int, float, int, str, int]:
    direction = -1 if strongest_first else 1
    return (
        direction * float(row.taxon_crumbs or 0.0),
        (
            direction * float(row.total_crumbs_score)
            if row.total_crumbs_score is not None
            else float("inf")
        ),
        direction * int(row.supporting_query_count),
        (
            direction * float(row.relative_crumbs_percent)
            if row.relative_crumbs_percent is not None
            else float("inf")
        ),
        SUPPORT_TIER_ORDER[row.support_tier],
        row.taxon_name,
        row.taxid,
    )


def bounded_signal_rows(
    tables: tuple[TaxonBigTableInput, ...],
    *,
    max_rows_per_sample: int,
    strongest_first: bool,
) -> tuple[tuple[TaxonReportRow, ...], int]:
    if max_rows_per_sample < 1:
        message = "Taxon signal row budget must be positive"
        raise ValueError(message)

    eligible_by_sample: dict[str, list[TaxonProducerRow]] = {}
    invalid: list[InvalidTaxonFindingRow] = []
    for table in tables:
        if isinstance(table, InvalidTaxonBigTable):
            invalid.append(
                InvalidTaxonFindingRow(
                    sample_id=table.receipt.sample_id,
                    validation_details=str(table.error),
                ),
            )
            continue
        for row in table.rows:
            if row.taxon_crumbs is not None and row.taxon_crumbs > 0:
                eligible_by_sample.setdefault(row.sample_id, []).append(row)

    selected: list[TaxonProducerRow] = []
    for sample_id in sorted(eligible_by_sample):
        sample_rows = sorted(
            eligible_by_sample[sample_id],
            key=lambda row: signal_sort_key(
                row,
                strongest_first=strongest_first,
            ),
        )
        selected.extend(sample_rows[:max_rows_per_sample])

    eligible_count = sum(len(rows) for rows in eligible_by_sample.values())
    return (
        (
            *finding_rows(selected),
            *sorted(invalid, key=lambda row: row.sample_id),
        ),
        eligible_count,
    )


def strongest_taxon_rows(
    tables: tuple[TaxonBigTableInput, ...],
    *,
    max_rows_per_sample: int = STRONGEST_SIGNALS_MAX_ROWS_PER_SAMPLE,
) -> tuple[tuple[TaxonReportRow, ...], int]:
    return bounded_signal_rows(
        tables,
        max_rows_per_sample=max_rows_per_sample,
        strongest_first=True,
    )


def faintest_taxon_rows(
    tables: tuple[TaxonBigTableInput, ...],
    *,
    max_rows_per_sample: int = FAINTEST_SIGNALS_MAX_ROWS_PER_SAMPLE,
) -> tuple[tuple[TaxonReportRow, ...], int]:
    return bounded_signal_rows(
        tables,
        max_rows_per_sample=max_rows_per_sample,
        strongest_first=False,
    )
