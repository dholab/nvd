"""Typed prepared-query summary parsing and report presentation."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Annotated, Literal, TypeAlias

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    NonNegativeInt,
    StringConstraints,
    ValidationError,
)

from py_nvd.multiqc_packages import (
    PreparedQueryBatchesReceipt,
    ReportPackage,
)

NonEmptyString = Annotated[str, StringConstraints(min_length=1)]
QueryClass: TypeAlias = Literal[
    "short_assembly_contig",
    "long_assembly_contig",
    "overlap_merged_pair",
    "single_read",
]
QuerySource: TypeAlias = Literal["contig", "read_query"]
Availability: TypeAlias = Literal["enabled", "disabled", "invalid"]

QUERY_CLASSES: tuple[QueryClass, ...] = (
    "short_assembly_contig",
    "long_assembly_contig",
    "overlap_merged_pair",
    "single_read",
)
QUERY_SOURCES: dict[QueryClass, QuerySource] = {
    "short_assembly_contig": "contig",
    "long_assembly_contig": "contig",
    "overlap_merged_pair": "read_query",
    "single_read": "read_query",
}
QUERY_CLASS_LABELS: dict[QueryClass, str] = {
    "short_assembly_contig": "Short assembly contig",
    "long_assembly_contig": "Long assembly contig",
    "overlap_merged_pair": "Overlap-merged pair",
    "single_read": "Single read",
}
SEQUENCE_TYPE_LABELS: dict[QuerySource, str] = {
    "contig": "Assembly contig",
    "read_query": "Unassembled read",
}


class ProducerModel(BaseModel):
    """Immutable projection of a trusted producer artifact."""

    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)


class PreparedQueryBatch(ProducerModel):
    sample_id: NonEmptyString
    platform: NonEmptyString
    query_class: QueryClass
    query_source: QuerySource
    n_query_sequences: NonNegativeInt
    query_fasta_present: bool
    query_lookup_present: bool


class PreparedQueryRow(BaseModel):
    """One user-facing prepared-query row."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    sample_id: str
    query_class_code: QueryClass | None = Field(default=None, exclude=True)
    query_class: str | None = None
    sequence_type: str | None = None
    availability: Availability
    query_sequences: NonNegativeInt | None = None
    platform: str | None = None
    problem: str | None = None
    validation_details: str | None = None


class PreparedQuerySummaryError(ValueError):
    """A prepared-query summary violates its semantic contract."""


PreparedQueryParseError: TypeAlias = (
    ValidationError | UnicodeDecodeError | csv.Error | PreparedQuerySummaryError
)


@dataclass(frozen=True, slots=True)
class ParsedPreparedQuerySummary:
    receipt: PreparedQueryBatchesReceipt
    batches: tuple[PreparedQueryBatch, ...]


@dataclass(frozen=True, slots=True)
class InvalidPreparedQuerySummary:
    receipt: PreparedQueryBatchesReceipt
    error: PreparedQueryParseError


PreparedQuerySummary: TypeAlias = (
    ParsedPreparedQuerySummary | InvalidPreparedQuerySummary
)


def parse_prepared_query_package(
    package: ReportPackage,
) -> ParsedPreparedQuerySummary:
    receipt = package.receipt
    if not isinstance(receipt, PreparedQueryBatchesReceipt):
        message = "Prepared-query parser received the wrong receipt"
        raise TypeError(message)

    summary_path = package.payload(receipt.summary)
    with summary_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(
            PreparedQueryBatch.model_fields,
        ).issubset(
            reader.fieldnames,
        ):
            message = "prepared-query summary is missing required columns"
            raise PreparedQuerySummaryError(message)
        batches = tuple(
            PreparedQueryBatch.model_validate_strings(row, strict=True)
            for row in reader
        )

    by_class: dict[QueryClass, PreparedQueryBatch] = {}
    for batch in batches:
        if batch.sample_id != receipt.sample_id:
            message = "prepared-query summary sample does not match its receipt"
            raise PreparedQuerySummaryError(message)
        if batch.query_class in by_class:
            message = f"duplicate prepared-query class: {batch.query_class}"
            raise PreparedQuerySummaryError(message)
        if batch.query_source != QUERY_SOURCES[batch.query_class]:
            message = (
                f"prepared-query source does not match query class: {batch.query_class}"
            )
            raise PreparedQuerySummaryError(message)
        if batch.query_fasta_present != batch.query_lookup_present:
            message = (
                "prepared-query FASTA/lookup presence differs for query class: "
                f"{batch.query_class}"
            )
            raise PreparedQuerySummaryError(message)
        if batch.n_query_sequences > 0 and not batch.query_fasta_present:
            message = (
                "prepared-query count has no batch artifacts for query class: "
                f"{batch.query_class}"
            )
            raise PreparedQuerySummaryError(message)
        by_class[batch.query_class] = batch

    observed_classes = set(by_class)
    expected_classes = set(QUERY_CLASSES)
    if observed_classes != expected_classes:
        missing = sorted(expected_classes - observed_classes)
        unexpected = sorted(observed_classes - expected_classes)
        message = (
            "prepared-query class membership differs; "
            f"missing={missing}, unexpected={unexpected}"
        )
        raise PreparedQuerySummaryError(message)

    return ParsedPreparedQuerySummary(
        receipt=receipt,
        batches=tuple(by_class[query_class] for query_class in QUERY_CLASSES),
    )


def collect_prepared_query_summaries(
    packages: tuple[ReportPackage, ...],
) -> tuple[PreparedQuerySummary, ...]:
    summaries: list[PreparedQuerySummary] = []
    for package in packages:
        if not isinstance(package.receipt, PreparedQueryBatchesReceipt):
            continue
        try:
            summaries.append(parse_prepared_query_package(package))
        except (
            ValidationError,
            UnicodeDecodeError,
            csv.Error,
            PreparedQuerySummaryError,
        ) as error:
            summaries.append(
                InvalidPreparedQuerySummary(receipt=package.receipt, error=error),
            )
    return tuple(summaries)


def prepared_query_rows(
    summaries: tuple[PreparedQuerySummary, ...],
    *,
    sample_order: dict[str, int],
    sample_platforms: dict[str, str],
    assembly_enabled: bool,
    unassembled_read_querying_enabled: bool,
) -> tuple[PreparedQueryRow, ...]:
    rows: list[PreparedQueryRow] = []
    for summary in summaries:
        if isinstance(summary, InvalidPreparedQuerySummary):
            rows.append(
                PreparedQueryRow(
                    sample_id=summary.receipt.sample_id,
                    availability="invalid",
                    problem="Prepared query summary could not be read",
                    validation_details=str(summary.error),
                ),
            )
            continue

        expected_platform = sample_platforms[summary.receipt.sample_id]
        if any(batch.platform != expected_platform for batch in summary.batches):
            rows.append(
                PreparedQueryRow(
                    sample_id=summary.receipt.sample_id,
                    availability="invalid",
                    problem="Prepared query summary conflicts with the sample roster",
                    validation_details=(
                        "prepared-query platform does not match the resolved roster"
                    ),
                ),
            )
            continue

        disabled_observations: list[QueryClass] = []
        for batch in summary.batches:
            enabled = (
                assembly_enabled
                if batch.query_source == "contig"
                else unassembled_read_querying_enabled
            )
            if not enabled and (
                batch.n_query_sequences > 0 or batch.query_fasta_present
            ):
                disabled_observations.append(batch.query_class)
        if disabled_observations:
            rows.append(
                PreparedQueryRow(
                    sample_id=summary.receipt.sample_id,
                    availability="invalid",
                    problem="Prepared query summary conflicts with the run plan",
                    validation_details=(
                        "disabled query classes have observed batches: "
                        + ", ".join(disabled_observations)
                    ),
                ),
            )
            continue

        for batch in summary.batches:
            enabled = query_class_enabled(
                batch.query_class,
                assembly_enabled=assembly_enabled,
                unassembled_read_querying_enabled=unassembled_read_querying_enabled,
            )
            if enabled and not batch.query_fasta_present:
                continue
            rows.append(
                PreparedQueryRow(
                    sample_id=summary.receipt.sample_id,
                    query_class_code=batch.query_class,
                    query_class=QUERY_CLASS_LABELS[batch.query_class],
                    sequence_type=SEQUENCE_TYPE_LABELS[batch.query_source],
                    availability="enabled" if enabled else "disabled",
                    query_sequences=batch.n_query_sequences if enabled else None,
                    platform=expected_platform,
                ),
            )

    class_order = {
        query_class: index for index, query_class in enumerate(QUERY_CLASSES)
    }
    return tuple(
        sorted(
            rows,
            key=lambda row: (
                sample_order[row.sample_id],
                class_order[row.query_class_code]
                if row.query_class_code is not None
                else len(class_order),
            ),
        ),
    )


def query_class_enabled(
    query_class: QueryClass,
    *,
    assembly_enabled: bool,
    unassembled_read_querying_enabled: bool,
) -> bool:
    return (
        assembly_enabled
        if QUERY_SOURCES[query_class] == "contig"
        else unassembled_read_querying_enabled
    )
