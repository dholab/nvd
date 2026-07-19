"""Typed MEGABLAST partition parsing and report presentation."""

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
    MegablastQueryPartitionReceipt,
    ReportPackage,
)
from py_nvd.multiqc_query_preparation import (
    QUERY_CLASS_LABELS,
    QUERY_CLASSES,
    QueryClass,
)

NonEmptyString = Annotated[str, StringConstraints(min_length=1)]


class ProducerModel(BaseModel):
    """Immutable projection of a trusted producer artifact."""

    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)


class MegablastQueryPartition(ProducerModel):
    sample_id: NonEmptyString
    query_class: QueryClass
    queries_in: NonNegativeInt
    megablast_accounted: NonNegativeInt
    blastn_candidates: NonNegativeInt


class BlastTaskRow(BaseModel):
    """One user-facing BLAST task-breakdown row."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    sample_id: str
    query_class_code: str = Field(exclude=True)
    query_class: str
    queries_searched: NonNegativeInt | None = None
    matched_by_megablast: NonNegativeInt | None = None
    forwarded_to_blastn: NonNegativeInt | None = None
    forwarded_to_blastn_percentage: float | None = None
    problem: str | None = None
    validation_details: str | None = None


class BlastDisabledRow(BaseModel):
    """One run-level notice when BLAST searching is disabled."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    sample_id: Literal["Run configuration"] = "Run configuration"
    notice: Literal["BLAST searching was disabled for this run."] = (
        "BLAST searching was disabled for this run."
    )


class MegablastQueryPartitionError(ValueError):
    """A MEGABLAST partition summary violates its semantic contract."""


MegablastQueryPartitionParseError: TypeAlias = (
    ValidationError | UnicodeDecodeError | csv.Error | MegablastQueryPartitionError
)


@dataclass(frozen=True, slots=True)
class ParsedMegablastQueryPartition:
    receipt: MegablastQueryPartitionReceipt
    partition: MegablastQueryPartition


@dataclass(frozen=True, slots=True)
class InvalidMegablastQueryPartition:
    receipt: MegablastQueryPartitionReceipt
    error: MegablastQueryPartitionParseError


MegablastQueryPartitionInput: TypeAlias = (
    ParsedMegablastQueryPartition | InvalidMegablastQueryPartition
)


def parse_megablast_query_partition_package(
    package: ReportPackage,
) -> ParsedMegablastQueryPartition:
    receipt = package.receipt
    if not isinstance(receipt, MegablastQueryPartitionReceipt):
        message = "MEGABLAST query partition parser received the wrong receipt"
        raise TypeError(message)

    summary_path = package.payload(receipt.summary)
    with summary_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(
            MegablastQueryPartition.model_fields,
        ).issubset(reader.fieldnames):
            message = "MEGABLAST query partition summary is missing required columns"
            raise MegablastQueryPartitionError(message)
        payload_rows = tuple(reader)

    if len(payload_rows) != 1:
        message = "MEGABLAST query partition summary must contain exactly one row"
        raise MegablastQueryPartitionError(message)

    partition = MegablastQueryPartition.model_validate_strings(
        payload_rows[0],
        strict=True,
    )
    if partition.sample_id != receipt.sample_id:
        message = "MEGABLAST query partition sample does not match its receipt"
        raise MegablastQueryPartitionError(message)
    if partition.query_class != receipt.query_class:
        message = "MEGABLAST query partition class does not match its receipt"
        raise MegablastQueryPartitionError(message)
    return ParsedMegablastQueryPartition(receipt=receipt, partition=partition)


def collect_megablast_query_partitions(
    packages: tuple[ReportPackage, ...],
) -> tuple[MegablastQueryPartitionInput, ...]:
    partitions: list[MegablastQueryPartitionInput] = []
    for package in packages:
        if not isinstance(package.receipt, MegablastQueryPartitionReceipt):
            continue
        try:
            partitions.append(parse_megablast_query_partition_package(package))
        except (
            ValidationError,
            UnicodeDecodeError,
            csv.Error,
            MegablastQueryPartitionError,
        ) as error:
            partitions.append(
                InvalidMegablastQueryPartition(
                    receipt=package.receipt,
                    error=error,
                ),
            )
    return tuple(partitions)


def blast_task_rows(
    partitions: tuple[MegablastQueryPartitionInput, ...],
    *,
    sample_order: dict[str, int],
) -> tuple[BlastTaskRow, ...]:
    rows: list[BlastTaskRow] = []
    for item in partitions:
        if isinstance(item, InvalidMegablastQueryPartition):
            rows.append(
                BlastTaskRow(
                    sample_id=item.receipt.sample_id,
                    query_class_code=item.receipt.query_class,
                    query_class=QUERY_CLASS_LABELS.get(
                        item.receipt.query_class,
                        item.receipt.query_class,
                    ),
                    problem="BLAST task breakdown could not be read",
                    validation_details=str(item.error),
                ),
            )
            continue

        partition = item.partition
        percentage = (
            round(100 * partition.blastn_candidates / partition.queries_in, 1)
            if partition.queries_in > 0
            else None
        )
        rows.append(
            BlastTaskRow(
                sample_id=partition.sample_id,
                query_class_code=partition.query_class,
                query_class=QUERY_CLASS_LABELS[partition.query_class],
                queries_searched=partition.queries_in,
                matched_by_megablast=partition.megablast_accounted,
                forwarded_to_blastn=partition.blastn_candidates,
                forwarded_to_blastn_percentage=percentage,
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
                class_order.get(row.query_class_code, len(class_order)),
            ),
        ),
    )
