"""Typed assembly package parsing and report-row derivation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Literal, TypeAlias, TypeVar

from pydantic import BaseModel, ConfigDict, NonNegativeInt, ValidationError

from py_nvd.multiqc_fastx import (
    InvalidFastxProfile,
    ParsedFastxProfile,
    collect_fastx_profiles,
)
from py_nvd.multiqc_packages import (
    EligibilityReceipt,
    LongReadEligibilityReceipt,
    PackageEnvelopeError,
    ReportPackage,
    SafeName,
    UnionReceipt,
)


class ProducerSummary(BaseModel):
    model_config = ConfigDict(extra="ignore", frozen=True, strict=True)


class AssemblerDecisionSummary(ProducerSummary):
    producer: SafeName
    length_floor: NonNegativeInt
    qualifying_reads: NonNegativeInt
    qualifying_bases: NonNegativeInt
    decision: Literal["run", "skip"]
    reason: str


class LongReadAssemblyAssessment(ProducerSummary):
    schema_version: Literal["nvd.long-read-assembly-eligibility/v1"]
    sample_id: str
    sequence_count: NonNegativeInt
    total_bases: NonNegativeInt
    max_read_length: NonNegativeInt | None
    assemblers: tuple[AssemblerDecisionSummary, ...]


class AssemblerUnionSummary(ProducerSummary):
    producer: SafeName
    source_count: NonNegativeInt
    source_bases: NonNegativeInt


class LongReadUnionSummary(ProducerSummary):
    schema_version: Literal["nvd.long-read-union-summary/v1"]
    sample_id: str
    source_contig_count: NonNegativeInt
    source_bases: NonNegativeInt
    representative_count: NonNegativeInt
    emitted_contig_count: NonNegativeInt
    emitted_bases: NonNegativeInt
    exact_duplicate_count: NonNegativeInt
    exact_containment_count: NonNegativeInt
    by_assembler: tuple[AssemblerUnionSummary, ...]


@dataclass(frozen=True, slots=True, order=True)
class AssemblyKey:
    sample_id: str
    producer: str


@dataclass(frozen=True, slots=True)
class AssemblyEligibilityDecision:
    key: AssemblyKey
    decision: Literal["run", "skip"]
    sequence_count: int
    reason: str


@dataclass(frozen=True, slots=True)
class MinimumSequenceCountDecision(AssemblyEligibilityDecision):
    minimum_sequence_count: int


@dataclass(frozen=True, slots=True)
class LongReadAssemblerDecision(AssemblyEligibilityDecision):
    total_bases: int
    max_read_length: int | None
    length_floor: int
    qualifying_reads: int
    qualifying_bases: int


@dataclass(frozen=True, slots=True)
class UnionInput:
    key: AssemblyKey
    source_count: int
    source_bases: int


@dataclass(frozen=True, slots=True)
class UnionAggregate(UnionInput):
    emitted_count: int
    emitted_bases: int
    representative_count: int
    duplicate_count: int
    containment_count: int


@dataclass(frozen=True, slots=True)
class AssemblerUnionContribution(UnionInput):
    pass


@dataclass(frozen=True, slots=True)
class InvalidAssemblyInput:
    key: AssemblyKey
    error: ValidationError | UnicodeDecodeError | PackageEnvelopeError


AssemblyReportInput: TypeAlias = (
    AssemblyEligibilityDecision
    | UnionInput
    | ParsedFastxProfile
    | InvalidFastxProfile
    | InvalidAssemblyInput
)


@dataclass(frozen=True, slots=True)
class AssemblyInputs:
    key: AssemblyKey
    eligibility: AssemblyEligibilityDecision | None
    profile: ParsedFastxProfile | None
    union: UnionInput | None
    input_errors: tuple[str, ...]


class AssemblyStatus(StrEnum):
    OBSERVED = "observed"
    COMPLETE_EMPTY = "complete_empty"
    ROUTED_AWAY = "routed_away"
    EXPECTED_OUTPUT_ABSENT = "expected_output_absent"
    INVALID = "invalid"


class AssemblyRow(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, strict=True)

    sample_id: str
    producer: str
    status: AssemblyStatus
    validation_message: str | None = None
    eligibility_status: Literal["observed"] | None = None
    eligibility_decision: Literal["run", "skip"] | None = None
    eligibility_reason: str | None = None
    eligibility_sequence_count: NonNegativeInt | None = None
    eligibility_threshold: NonNegativeInt | None = None
    eligibility_total_bases: NonNegativeInt | None = None
    eligibility_max_read_length: NonNegativeInt | None = None
    eligibility_length_floor: NonNegativeInt | None = None
    eligibility_qualifying_reads: NonNegativeInt | None = None
    eligibility_qualifying_bases: NonNegativeInt | None = None
    sequence_count: NonNegativeInt | None = None
    total_bases: NonNegativeInt | None = None
    min_length: NonNegativeInt | None = None
    max_length: NonNegativeInt | None = None
    mean_length: int | float | None = None
    n50: NonNegativeInt | None = None
    l50: NonNegativeInt | None = None
    source_count: NonNegativeInt | None = None
    source_bases: NonNegativeInt | None = None
    emitted_contig_count: NonNegativeInt | None = None
    emitted_bases: NonNegativeInt | None = None
    representative_count: NonNegativeInt | None = None
    exact_duplicate_count: NonNegativeInt | None = None
    exact_containment_count: NonNegativeInt | None = None


def parse_long_read_eligibility(
    package: ReportPackage,
) -> tuple[AssemblyReportInput, ...]:
    receipt = package.receipt
    if not isinstance(receipt, LongReadEligibilityReceipt):
        message = "long-read eligibility parser received the wrong receipt"
        raise TypeError(message)
    try:
        assessment = LongReadAssemblyAssessment.model_validate_json(
            package.payload(receipt.summary).read_text(encoding="utf-8"),
            strict=True,
        )
    except (ValidationError, UnicodeDecodeError) as error:
        return (
            InvalidAssemblyInput(
                key=AssemblyKey(receipt.sample_id, receipt.producer),
                error=error,
            ),
        )
    if assessment.sample_id != receipt.sample_id:
        return (
            InvalidAssemblyInput(
                key=AssemblyKey(receipt.sample_id, receipt.producer),
                error=PackageEnvelopeError(
                    "long-read eligibility sample does not match its receipt",
                ),
            ),
        )
    return tuple(
        LongReadAssemblerDecision(
            key=AssemblyKey(assessment.sample_id, assembler.producer),
            decision=assembler.decision,
            sequence_count=assessment.sequence_count,
            reason=assembler.reason,
            total_bases=assessment.total_bases,
            max_read_length=assessment.max_read_length,
            length_floor=assembler.length_floor,
            qualifying_reads=assembler.qualifying_reads,
            qualifying_bases=assembler.qualifying_bases,
        )
        for assembler in assessment.assemblers
    )


def parse_union_summary(package: ReportPackage) -> tuple[AssemblyReportInput, ...]:
    receipt = package.receipt
    if not isinstance(receipt, UnionReceipt):
        message = "union parser received the wrong receipt"
        raise TypeError(message)
    try:
        summary = LongReadUnionSummary.model_validate_json(
            package.payload(receipt.summary).read_text(encoding="utf-8"),
            strict=True,
        )
    except (ValidationError, UnicodeDecodeError) as error:
        return (
            InvalidAssemblyInput(
                key=AssemblyKey(receipt.sample_id, receipt.producer),
                error=error,
            ),
        )
    if summary.sample_id != receipt.sample_id:
        return (
            InvalidAssemblyInput(
                key=AssemblyKey(receipt.sample_id, receipt.producer),
                error=PackageEnvelopeError(
                    "union summary sample does not match its receipt",
                ),
            ),
        )
    aggregate = UnionAggregate(
        key=AssemblyKey(summary.sample_id, receipt.producer),
        source_count=summary.source_contig_count,
        source_bases=summary.source_bases,
        emitted_count=summary.emitted_contig_count,
        emitted_bases=summary.emitted_bases,
        representative_count=summary.representative_count,
        duplicate_count=summary.exact_duplicate_count,
        containment_count=summary.exact_containment_count,
    )
    contributions = tuple(
        AssemblerUnionContribution(
            key=AssemblyKey(summary.sample_id, assembler.producer),
            source_count=assembler.source_count,
            source_bases=assembler.source_bases,
        )
        for assembler in summary.by_assembler
    )
    return (aggregate, *contributions)


def collect_assembly_report_inputs(
    packages: tuple[ReportPackage, ...],
) -> tuple[AssemblyReportInput, ...]:
    report_inputs: list[AssemblyReportInput] = list(
        collect_fastx_profiles(packages, domain="assembly"),
    )
    for package in packages:
        receipt = package.receipt
        if isinstance(receipt, EligibilityReceipt):
            report_inputs.append(
                MinimumSequenceCountDecision(
                    key=AssemblyKey(receipt.sample_id, receipt.producer),
                    decision=receipt.decision,
                    sequence_count=receipt.sequence_count,
                    reason=receipt.reason,
                    minimum_sequence_count=receipt.threshold,
                ),
            )
        elif isinstance(receipt, LongReadEligibilityReceipt):
            report_inputs.extend(parse_long_read_eligibility(package))
        elif isinstance(receipt, UnionReceipt):
            report_inputs.extend(parse_union_summary(package))
    return tuple(report_inputs)


def assembly_input_key(report_input: AssemblyReportInput) -> AssemblyKey:
    if isinstance(report_input, ParsedFastxProfile | InvalidFastxProfile):
        producer = report_input.receipt.producer or "missing_producer"
        return AssemblyKey(report_input.receipt.sample_id, producer)
    return report_input.key


T = TypeVar("T")


def single_input(inputs: list[T], role: str) -> tuple[T | None, str | None]:
    if len(inputs) == 1:
        return inputs[0], None
    if not inputs:
        return None, None
    return None, f"multiple {role} inputs"


def parse_assembly_inputs(
    key: AssemblyKey,
    report_inputs: tuple[AssemblyReportInput, ...],
) -> AssemblyInputs:
    errors = [
        str(report_input.error)
        for report_input in report_inputs
        if isinstance(report_input, InvalidAssemblyInput | InvalidFastxProfile)
    ]
    eligibility_inputs = [
        report_input
        for report_input in report_inputs
        if isinstance(report_input, AssemblyEligibilityDecision)
    ]
    profile_inputs = [
        report_input
        for report_input in report_inputs
        if isinstance(report_input, ParsedFastxProfile)
        and report_input.receipt.stage == "assembly"
    ]
    if any(
        isinstance(report_input, ParsedFastxProfile)
        and report_input.receipt.stage != "assembly"
        for report_input in report_inputs
    ):
        errors.append("assembly profile receipt has a non-assembly stage")
    union_inputs = [
        report_input
        for report_input in report_inputs
        if isinstance(report_input, UnionInput)
    ]

    eligibility, eligibility_error = single_input(
        eligibility_inputs,
        "assembly eligibility decision",
    )
    profile, profile_error = single_input(profile_inputs, "assembly profile")
    union, union_error = single_input(union_inputs, "union summary")
    errors.extend(
        error
        for error in (eligibility_error, profile_error, union_error)
        if error is not None
    )
    if eligibility is None and profile is None and union is None and not errors:
        errors.append("assembly package group has no recognized input")
    return AssemblyInputs(
        key=key,
        eligibility=eligibility,
        profile=profile,
        union=union,
        input_errors=tuple(errors),
    )


def union_output_count(union: UnionInput) -> int:
    if isinstance(union, UnionAggregate):
        return union.emitted_count
    return union.source_count


def union_output_bases(union: UnionInput) -> int:
    if isinstance(union, UnionAggregate):
        return union.emitted_bases
    return union.source_bases


def assembly_conflicts(inputs: AssemblyInputs) -> tuple[str, ...]:
    conflicts = list(inputs.input_errors)
    if (
        inputs.eligibility is not None
        and inputs.eligibility.decision == "skip"
        and (inputs.profile is not None or inputs.union is not None)
    ):
        conflicts.append("eligibility skip conflicts with present output")
    if (
        inputs.profile is not None
        and inputs.union is not None
        and (
            inputs.profile.profile.sequence_count != union_output_count(inputs.union)
            or inputs.profile.profile.total_bases != union_output_bases(inputs.union)
        )
    ):
        conflicts.append("profile conflicts with union summary")
    return tuple(sorted(set(conflicts)))


def assembly_status(
    inputs: AssemblyInputs,
    conflicts: tuple[str, ...],
) -> AssemblyStatus:
    if conflicts:
        return AssemblyStatus.INVALID
    if inputs.profile is not None:
        return (
            AssemblyStatus.COMPLETE_EMPTY
            if inputs.profile.profile.sequence_count == 0
            else AssemblyStatus.OBSERVED
        )
    if inputs.union is not None:
        return (
            AssemblyStatus.COMPLETE_EMPTY
            if union_output_count(inputs.union) == 0
            else AssemblyStatus.OBSERVED
        )
    if inputs.eligibility is not None and inputs.eligibility.decision == "skip":
        return AssemblyStatus.ROUTED_AWAY
    if inputs.eligibility is not None and inputs.eligibility.decision == "run":
        return AssemblyStatus.EXPECTED_OUTPUT_ABSENT
    return AssemblyStatus.INVALID


def derive_assembly_row(inputs: AssemblyInputs) -> AssemblyRow:
    conflicts = assembly_conflicts(inputs)
    eligibility = inputs.eligibility
    minimum_count_decision = (
        eligibility if isinstance(eligibility, MinimumSequenceCountDecision) else None
    )
    long_read_decision = (
        eligibility if isinstance(eligibility, LongReadAssemblerDecision) else None
    )
    profile = inputs.profile.profile if inputs.profile is not None else None
    aggregate = inputs.union if isinstance(inputs.union, UnionAggregate) else None
    return AssemblyRow(
        sample_id=inputs.key.sample_id,
        producer=inputs.key.producer,
        status=assembly_status(inputs, conflicts),
        validation_message="; ".join(conflicts) or None,
        eligibility_status="observed" if eligibility is not None else None,
        eligibility_decision=eligibility.decision if eligibility is not None else None,
        eligibility_reason=eligibility.reason if eligibility is not None else None,
        eligibility_sequence_count=(
            eligibility.sequence_count if eligibility is not None else None
        ),
        eligibility_threshold=(
            minimum_count_decision.minimum_sequence_count
            if minimum_count_decision is not None
            else None
        ),
        eligibility_total_bases=(
            long_read_decision.total_bases if long_read_decision is not None else None
        ),
        eligibility_max_read_length=(
            long_read_decision.max_read_length
            if long_read_decision is not None
            else None
        ),
        eligibility_length_floor=(
            long_read_decision.length_floor if long_read_decision is not None else None
        ),
        eligibility_qualifying_reads=(
            long_read_decision.qualifying_reads
            if long_read_decision is not None
            else None
        ),
        eligibility_qualifying_bases=(
            long_read_decision.qualifying_bases
            if long_read_decision is not None
            else None
        ),
        sequence_count=profile.sequence_count if profile is not None else None,
        total_bases=profile.total_bases if profile is not None else None,
        min_length=profile.length.min if profile is not None else None,
        max_length=profile.length.max if profile is not None else None,
        mean_length=profile.length.mean if profile is not None else None,
        n50=profile.length.n50 if profile is not None else None,
        l50=profile.length.l50 if profile is not None else None,
        source_count=inputs.union.source_count if inputs.union is not None else None,
        source_bases=inputs.union.source_bases if inputs.union is not None else None,
        emitted_contig_count=(
            aggregate.emitted_count if aggregate is not None else None
        ),
        emitted_bases=aggregate.emitted_bases if aggregate is not None else None,
        representative_count=(
            aggregate.representative_count if aggregate is not None else None
        ),
        exact_duplicate_count=(
            aggregate.duplicate_count if aggregate is not None else None
        ),
        exact_containment_count=(
            aggregate.containment_count if aggregate is not None else None
        ),
    )


def assembly_rows(
    packages: tuple[ReportPackage, ...],
    *,
    sample_order: dict[str, int],
) -> tuple[AssemblyRow, ...]:
    grouped: dict[AssemblyKey, list[AssemblyReportInput]] = {}
    for report_input in collect_assembly_report_inputs(packages):
        grouped.setdefault(assembly_input_key(report_input), []).append(report_input)
    keys = sorted(
        grouped,
        key=lambda key: (
            sample_order.get(key.sample_id, len(sample_order)),
            key.producer,
        ),
    )
    return tuple(
        derive_assembly_row(parse_assembly_inputs(key, tuple(grouped[key])))
        for key in keys
    )
