#!/usr/bin/env python3
"""Run long-read assemblers only when cached profiles show usable reads."""

from __future__ import annotations

import argparse
import csv
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, Literal, TypeAlias

from pydantic import BaseModel, ConfigDict, Field, NonNegativeInt, StringConstraints

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence


ASSESSMENT_FIELDS = (
    "sample_id",
    "sequence_count",
    "total_bases",
    "max_read_length",
)
DECISION_FIELDS = (
    "length_floor",
    "qualifying_reads",
    "qualifying_bases",
    "decision",
    "reason",
)
Number: TypeAlias = int | float
TsvCell: TypeAlias = str | int | float | None
NonEmptyString = Annotated[str, StringConstraints(min_length=1)]


class Assembler(StrEnum):
    METAMDBG = "metamdbg"
    MYLOASM = "myloasm"
    METAFLYE = "metaflye"


ASSEMBLER_THRESHOLDS = {
    Assembler.METAMDBG: "metamdbg_min_read_overlap",
    Assembler.MYLOASM: "myloasm_min_read_length",
    Assembler.METAFLYE: "metaflye_min_read_length",
}
REPORT_FIELDS = [
    *ASSESSMENT_FIELDS,
    *(
        f"{assembler.value}_{field}"
        for assembler in Assembler
        for field in DECISION_FIELDS
    ),
]


class ProducerProfile(BaseModel):
    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)


class LengthProfile(ProducerProfile):
    max: NonNegativeInt | None


class ThresholdProfile(ProducerProfile):
    name: NonEmptyString
    value: Number
    sequence_count_at_or_above: NonNegativeInt


class LengthThresholdProfile(ThresholdProfile):
    axis: Literal["length"]
    bases_at_or_above: NonNegativeInt


class QualityThresholdProfile(ThresholdProfile):
    axis: Literal["quality"]
    bases_at_or_above: None


ParsedThreshold: TypeAlias = Annotated[
    LengthThresholdProfile | QualityThresholdProfile,
    Field(discriminator="axis"),
]


class FastxProfile(ProducerProfile):
    sample_id: NonEmptyString
    sequence_count: NonNegativeInt
    total_bases: NonNegativeInt
    length: LengthProfile
    thresholds: tuple[ParsedThreshold, ...]


class AssessmentModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)


class AssemblerDecision(AssessmentModel):
    producer: Assembler
    length_floor: Number
    qualifying_reads: NonNegativeInt
    qualifying_bases: NonNegativeInt
    decision: Literal["run", "skip"]
    reason: Literal["qualifying_reads_present", "no_reads_meet_length_floor"]


class LongReadAssemblyAssessment(AssessmentModel):
    schema_version: Literal["nvd.long-read-assembly-eligibility/v1"] = (
        "nvd.long-read-assembly-eligibility/v1"
    )
    sample_id: NonEmptyString
    sequence_count: NonNegativeInt
    total_bases: NonNegativeInt
    max_read_length: NonNegativeInt | None
    assemblers: tuple[AssemblerDecision, ...]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assess long-read assembler eligibility from a FASTX profile.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    assess = commands.add_parser("assess")
    assess.add_argument("--profile", type=Path, required=True)
    assess.add_argument("--output", type=Path, required=True)
    assess.add_argument("--summary-json", type=Path)
    assess.add_argument("--metamdbg-marker", type=Path, required=True)
    assess.add_argument("--myloasm-marker", type=Path, required=True)
    assess.add_argument("--metaflye-marker", type=Path, required=True)

    report = commands.add_parser("report")
    report.add_argument("--inputs", type=Path, nargs="*", required=True)
    report.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def display_number(value: Number) -> Number:
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def parse_profile(path: Path) -> FastxProfile:
    return FastxProfile.model_validate_json(
        path.read_text(encoding="utf-8"),
        strict=True,
    )


def thresholds_by_name(profile: FastxProfile) -> dict[str, ParsedThreshold]:
    thresholds: dict[str, ParsedThreshold] = {}
    for threshold in profile.thresholds:
        if threshold.name in thresholds:
            message = f"Duplicate FASTX profile threshold: {threshold.name}"
            raise ValueError(message)
        thresholds[threshold.name] = threshold
    return thresholds


def require_length_threshold(
    thresholds: Mapping[str, ParsedThreshold],
    name: str,
) -> LengthThresholdProfile:
    if name not in thresholds:
        message = f"FASTX profile is missing threshold {name!r}."
        raise ValueError(message)
    threshold = thresholds[name]
    if not isinstance(threshold, LengthThresholdProfile):
        message = f"FASTX profile threshold {name!r} must have axis 'length'."
        raise TypeError(message)
    return threshold


def assess_profile(profile: FastxProfile) -> LongReadAssemblyAssessment:
    thresholds = thresholds_by_name(profile)
    decisions: list[AssemblerDecision] = []
    for assembler, threshold_name in ASSEMBLER_THRESHOLDS.items():
        threshold = require_length_threshold(thresholds, threshold_name)
        qualifying_reads = threshold.sequence_count_at_or_above
        decision = "run" if qualifying_reads > 0 else "skip"
        decisions.append(
            AssemblerDecision(
                producer=assembler,
                length_floor=display_number(threshold.value),
                qualifying_reads=qualifying_reads,
                qualifying_bases=threshold.bases_at_or_above,
                decision=decision,
                reason=(
                    "qualifying_reads_present"
                    if decision == "run"
                    else "no_reads_meet_length_floor"
                ),
            ),
        )
    return LongReadAssemblyAssessment(
        sample_id=profile.sample_id,
        sequence_count=profile.sequence_count,
        total_bases=profile.total_bases,
        max_read_length=profile.length.max,
        assemblers=tuple(decisions),
    )


def assessment_tsv_row(
    assessment: LongReadAssemblyAssessment,
) -> dict[str, TsvCell]:
    row: dict[str, TsvCell] = {
        "sample_id": assessment.sample_id,
        "sequence_count": assessment.sequence_count,
        "total_bases": assessment.total_bases,
        "max_read_length": assessment.max_read_length,
    }
    for assembler in assessment.assemblers:
        prefix = assembler.producer.value
        row.update(
            {
                f"{prefix}_length_floor": assembler.length_floor,
                f"{prefix}_qualifying_reads": assembler.qualifying_reads,
                f"{prefix}_qualifying_bases": assembler.qualifying_bases,
                f"{prefix}_decision": assembler.decision,
                f"{prefix}_reason": assembler.reason,
            },
        )
    return row


def write_assessment_tsv(path: Path, assessment: LongReadAssemblyAssessment) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(assessment_tsv_row(assessment))


def write_markers(
    assessment: LongReadAssemblyAssessment,
    markers: Mapping[Assembler, Path],
) -> None:
    decisions = {
        decision.producer: decision.decision for decision in assessment.assemblers
    }
    for assembler, marker in markers.items():
        marker.unlink(missing_ok=True)
        if decisions[assembler] == "run":
            marker.touch()


def write_summary_json(path: Path, assessment: LongReadAssemblyAssessment) -> None:
    path.write_text(assessment.model_dump_json(indent=2) + "\n", encoding="utf-8")


def assess(args: argparse.Namespace) -> None:
    assessment = assess_profile(parse_profile(args.profile))
    write_assessment_tsv(args.output, assessment)
    if args.summary_json is not None:
        write_summary_json(args.summary_json, assessment)
    write_markers(
        assessment,
        {
            Assembler.METAMDBG: args.metamdbg_marker,
            Assembler.MYLOASM: args.myloasm_marker,
            Assembler.METAFLYE: args.metaflye_marker,
        },
    )


def read_report_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REPORT_FIELDS:
            message = f"Unexpected eligibility report columns in {path}."
            raise ValueError(message)
        rows = list(reader)
    if len(rows) != 1:
        message = f"Expected one eligibility row in {path}, found {len(rows)}."
        raise ValueError(message)
    return rows[0]


def combine_reports(args: argparse.Namespace) -> None:
    rows = sorted(
        (read_report_row(path) for path in args.inputs),
        key=lambda row: row["sample_id"],
    )
    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        message = "Eligibility reports contain duplicate sample IDs."
        raise ValueError(message)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.command == "assess":
        assess(args)
    else:
        combine_reports(args)


if __name__ == "__main__":
    main()
