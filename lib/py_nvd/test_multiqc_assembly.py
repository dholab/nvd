"""Tests for typed assembly report inputs and row derivation."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

from py_nvd.multiqc_assembly import (
    AssemblerUnionContribution,
    AssemblyKey,
    AssemblyStatus,
    LongReadAssemblerDecision,
    MinimumSequenceCountDecision,
    derive_assembly_row,
    parse_assembly_inputs,
    parse_long_read_eligibility,
)
from py_nvd.multiqc_packages import (
    LongReadEligibilityReceipt,
    ReportPackage,
    write_package,
)

if TYPE_CHECKING:
    from pathlib import Path


def test_skipped_eligibility_without_output_is_routed_away() -> None:
    key = AssemblyKey("sample_A", "spades")
    decision = MinimumSequenceCountDecision(
        key,
        decision="skip",
        sequence_count=20,
        minimum_sequence_count=100,
        reason="short_read_minimum_sequence_count",
    )
    row = derive_assembly_row(parse_assembly_inputs(key, (decision,)))
    assert row.status is AssemblyStatus.ROUTED_AWAY


def test_skipped_eligibility_with_union_output_is_invalid() -> None:
    key = AssemblyKey("sample_A", "myloasm")
    decision = MinimumSequenceCountDecision(
        key,
        decision="skip",
        sequence_count=0,
        minimum_sequence_count=100,
        reason="no_reads_meet_length_floor",
    )
    union = AssemblerUnionContribution(
        key=key,
        source_count=1,
        source_bases=50,
    )
    row = derive_assembly_row(parse_assembly_inputs(key, (decision, union)))
    assert row.status is AssemblyStatus.INVALID
    assert row.validation_message == "eligibility skip conflicts with present output"


def test_eligibility_parser_projects_only_report_fields(tmp_path: Path) -> None:
    summary = tmp_path / "eligibility.json"
    summary.write_text(
        json.dumps(
            {
                "schema_version": "nvd.long-read-assembly-eligibility/v1",
                "sample_id": "sample_A",
                "sequence_count": 120,
                "total_bases": 120_000,
                "max_read_length": 10_000,
                "producer_metadata": {"revision": "future-addition"},
                "assemblers": [
                    {
                        "producer": "myloasm",
                        "length_floor": 1000,
                        "qualifying_reads": 100,
                        "qualifying_bases": 110_000,
                        "decision": "run",
                        "reason": "eligible",
                        "producer_note": "future-addition",
                    },
                ],
            },
        ),
        encoding="utf-8",
    )
    receipt = LongReadEligibilityReceipt(sample_id="sample_A")
    package_dir = write_package(receipt, (summary,), output_root=tmp_path)

    [decision] = parse_long_read_eligibility(ReportPackage(package_dir, receipt))

    assert isinstance(decision, LongReadAssemblerDecision)
    assert decision.key == AssemblyKey("sample_A", "myloasm")


def test_run_decision_without_output_is_reported_as_absent() -> None:
    key = AssemblyKey("sample_A", "spades")
    decision = MinimumSequenceCountDecision(
        key=key,
        decision="run",
        sequence_count=120,
        reason="short_read_minimum_sequence_count",
        minimum_sequence_count=100,
    )

    row = derive_assembly_row(parse_assembly_inputs(key, (decision,)))

    assert row.status is AssemblyStatus.EXPECTED_OUTPUT_ABSENT


def test_duplicate_eligibility_decisions_are_not_first_wins() -> None:
    key = AssemblyKey("sample_A", "spades")
    decisions = tuple(
        MinimumSequenceCountDecision(
            key=key,
            decision=decision,
            sequence_count=120,
            reason="short_read_minimum_sequence_count",
            minimum_sequence_count=100,
        )
        for decision in ("run", "skip")
    )

    row = derive_assembly_row(parse_assembly_inputs(key, decisions))

    assert row.status is AssemblyStatus.INVALID
    assert row.eligibility_decision is None
    assert row.validation_message == "multiple assembly eligibility decision inputs"


def test_duplicate_union_contributions_are_not_first_wins() -> None:
    key = AssemblyKey("sample_A", "myloasm")
    contributions = (
        AssemblerUnionContribution(key=key, source_count=1, source_bases=100),
        AssemblerUnionContribution(key=key, source_count=2, source_bases=200),
    )

    row = derive_assembly_row(parse_assembly_inputs(key, contributions))

    assert row.status is AssemblyStatus.INVALID
    assert row.source_count is None
    assert row.validation_message == "multiple union summary inputs"
