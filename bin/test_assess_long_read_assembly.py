"""Tests for cached-profile long-read assembler decisions."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING

import pytest
from assess_long_read_assembly import REPORT_FIELDS, main

if TYPE_CHECKING:
    from pathlib import Path


def write_profile(
    path: Path,
    *,
    sample_id: str,
    metamdbg_reads: int,
    myloasm_reads: int,
    metaflye_reads: int,
) -> None:
    path.write_text(
        json.dumps(
            {
                "sample_id": sample_id,
                "sequence_count": 4,
                "total_bases": 2399,
                "length": {"min": 150, "max": 1000, "mean": 599.75},
                "thresholds": [
                    {
                        "name": "metamdbg_min_read_overlap",
                        "axis": "length",
                        "value": 200,
                        "sequence_count_at_or_above": metamdbg_reads,
                        "bases_at_or_above": 2249 if metamdbg_reads else 0,
                    },
                    {
                        "name": "myloasm_min_read_length",
                        "axis": "length",
                        "value": 1000,
                        "sequence_count_at_or_above": myloasm_reads,
                        "bases_at_or_above": 1000 if myloasm_reads else 0,
                    },
                    {
                        "name": "metaflye_min_read_length",
                        "axis": "length",
                        "value": 1001,
                        "sequence_count_at_or_above": metaflye_reads,
                        "bases_at_or_above": 1200 if metaflye_reads else 0,
                    },
                ],
            },
        )
        + "\n",
        encoding="utf-8",
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def assess(tmp_path: Path, *, sample_id: str, counts: tuple[int, int, int]) -> Path:
    profile = tmp_path / f"{sample_id}.fastx_profile.json"
    report = tmp_path / f"{sample_id}.long_read_assembly_eligibility.tsv"
    write_profile(
        profile,
        sample_id=sample_id,
        metamdbg_reads=counts[0],
        myloasm_reads=counts[1],
        metaflye_reads=counts[2],
    )
    main(
        [
            "assess",
            "--profile",
            str(profile),
            "--output",
            str(report),
            "--metamdbg-marker",
            str(tmp_path / f"{sample_id}.metamdbg.run"),
            "--myloasm-marker",
            str(tmp_path / f"{sample_id}.myloasm.run"),
            "--metaflye-marker",
            str(tmp_path / f"{sample_id}.metaflye.run"),
        ],
    )
    return report


def test_assess_uses_threshold_counts_for_run_and_skip_decisions(
    tmp_path: Path,
) -> None:
    report = assess(tmp_path, sample_id="sample-1", counts=(3, 1, 0))

    assert (tmp_path / "sample-1.metamdbg.run").exists()
    assert (tmp_path / "sample-1.myloasm.run").exists()
    assert not (tmp_path / "sample-1.metaflye.run").exists()
    assert read_rows(report) == [
        {
            "sample_id": "sample-1",
            "sequence_count": "4",
            "total_bases": "2399",
            "max_read_length": "1000",
            "metamdbg_length_floor": "200",
            "metamdbg_qualifying_reads": "3",
            "metamdbg_qualifying_bases": "2249",
            "metamdbg_decision": "run",
            "metamdbg_reason": "qualifying_reads_present",
            "myloasm_length_floor": "1000",
            "myloasm_qualifying_reads": "1",
            "myloasm_qualifying_bases": "1000",
            "myloasm_decision": "run",
            "myloasm_reason": "qualifying_reads_present",
            "metaflye_length_floor": "1001",
            "metaflye_qualifying_reads": "0",
            "metaflye_qualifying_bases": "0",
            "metaflye_decision": "skip",
            "metaflye_reason": "no_reads_meet_length_floor",
        },
    ]


def test_report_combines_rows_in_sample_order(tmp_path: Path) -> None:
    second = assess(tmp_path, sample_id="sample-2", counts=(1, 0, 0))
    first = assess(tmp_path, sample_id="sample-1", counts=(1, 1, 1))
    combined = tmp_path / "long_read_assembly_eligibility.tsv"

    main(
        [
            "report",
            "--inputs",
            str(second),
            str(first),
            "--output",
            str(combined),
        ],
    )

    assert [row["sample_id"] for row in read_rows(combined)] == [
        "sample-1",
        "sample-2",
    ]


def test_report_writes_header_when_there_are_no_long_read_samples(
    tmp_path: Path,
) -> None:
    combined = tmp_path / "long_read_assembly_eligibility.tsv"

    main(
        [
            "report",
            "--inputs",
            "--output",
            str(combined),
        ],
    )

    assert read_rows(combined) == []
    assert combined.read_text(encoding="utf-8").splitlines() == [
        "\t".join(REPORT_FIELDS),
    ]


def test_assess_reports_samples_with_no_eligible_assembler(tmp_path: Path) -> None:
    report = assess(tmp_path, sample_id="sample-1", counts=(0, 0, 0))

    assert not (tmp_path / "sample-1.metamdbg.run").exists()
    assert not (tmp_path / "sample-1.myloasm.run").exists()
    assert not (tmp_path / "sample-1.metaflye.run").exists()
    [row] = read_rows(report)
    assert row["metamdbg_decision"] == "skip"
    assert row["myloasm_decision"] == "skip"
    assert row["metaflye_decision"] == "skip"


def test_assess_rejects_missing_required_threshold(tmp_path: Path) -> None:
    profile = tmp_path / "sample-1.fastx_profile.json"
    write_profile(
        profile,
        sample_id="sample-1",
        metamdbg_reads=1,
        myloasm_reads=1,
        metaflye_reads=1,
    )
    payload = json.loads(profile.read_text(encoding="utf-8"))
    payload["thresholds"] = payload["thresholds"][:-1]
    profile.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="metaflye_min_read_length"):
        main(
            [
                "assess",
                "--profile",
                str(profile),
                "--output",
                str(tmp_path / "report.tsv"),
                "--metamdbg-marker",
                str(tmp_path / "metamdbg.run"),
                "--myloasm-marker",
                str(tmp_path / "myloasm.run"),
                "--metaflye-marker",
                str(tmp_path / "metaflye.run"),
            ],
        )
