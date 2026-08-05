#!/usr/bin/env python3
"""Run long-read assemblers only when cached profiles show usable reads."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any


ASSEMBLER_THRESHOLDS = {
    "metamdbg": "metamdbg_min_read_overlap",
    "myloasm": "myloasm_min_read_length",
    "metaflye": "metaflye_min_read_length",
}

REPORT_FIELDS = [
    "sample_id",
    "sequence_count",
    "total_bases",
    "max_read_length",
    "metamdbg_length_floor",
    "metamdbg_qualifying_reads",
    "metamdbg_qualifying_bases",
    "metamdbg_decision",
    "metamdbg_reason",
    "myloasm_length_floor",
    "myloasm_qualifying_reads",
    "myloasm_qualifying_bases",
    "myloasm_decision",
    "myloasm_reason",
    "metaflye_length_floor",
    "metaflye_qualifying_reads",
    "metaflye_qualifying_bases",
    "metaflye_decision",
    "metaflye_reason",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assess long-read assembler eligibility from a FASTX profile.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    assess = commands.add_parser("assess")
    assess.add_argument("--profile", type=Path, required=True)
    assess.add_argument("--output", type=Path, required=True)
    assess.add_argument("--metamdbg-marker", type=Path, required=True)
    assess.add_argument("--myloasm-marker", type=Path, required=True)
    assess.add_argument("--metaflye-marker", type=Path, required=True)

    report = commands.add_parser("report")
    report.add_argument("--inputs", type=Path, nargs="*", required=True)
    report.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def display_number(value: object) -> object:
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return value


def load_profile(path: Path) -> dict[str, Any]:
    profile = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(profile, dict):
        message = f"FASTX profile must contain a JSON object: {path}"
        raise TypeError(message)
    return profile


def thresholds_by_name(profile: dict[str, Any]) -> dict[str, dict[str, Any]]:
    thresholds: dict[str, dict[str, Any]] = {}
    for threshold in profile.get("thresholds", []):
        name = threshold.get("name")
        if not isinstance(name, str):
            message = "FASTX profile threshold is missing a string name."
            raise TypeError(message)
        if name in thresholds:
            message = f"Duplicate FASTX profile threshold: {name}"
            raise ValueError(message)
        thresholds[name] = threshold
    return thresholds


def assessment_row(profile: dict[str, Any]) -> dict[str, object]:
    thresholds = thresholds_by_name(profile)
    row: dict[str, object] = {
        "sample_id": profile["sample_id"],
        "sequence_count": profile["sequence_count"],
        "total_bases": profile["total_bases"],
        "max_read_length": profile["length"]["max"],
    }
    for assembler, threshold_name in ASSEMBLER_THRESHOLDS.items():
        if threshold_name not in thresholds:
            message = f"FASTX profile is missing threshold {threshold_name!r}."
            raise ValueError(message)
        threshold = thresholds[threshold_name]
        qualifying_reads = int(threshold["sequence_count_at_or_above"])
        decision = "run" if qualifying_reads > 0 else "skip"
        row.update(
            {
                f"{assembler}_length_floor": display_number(threshold["value"]),
                f"{assembler}_qualifying_reads": qualifying_reads,
                f"{assembler}_qualifying_bases": threshold["bases_at_or_above"],
                f"{assembler}_decision": decision,
                f"{assembler}_reason": (
                    "qualifying_reads_present"
                    if decision == "run"
                    else "no_reads_meet_length_floor"
                ),
            },
        )
    return row


def write_row(path: Path, row: dict[str, object]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def write_markers(row: dict[str, object], markers: dict[str, Path]) -> None:
    for assembler, marker in markers.items():
        marker.unlink(missing_ok=True)
        if row[f"{assembler}_decision"] == "run":
            marker.touch()


def assess(args: argparse.Namespace) -> None:
    row = assessment_row(load_profile(args.profile))
    write_row(args.output, row)
    write_markers(
        row,
        {
            "metamdbg": args.metamdbg_marker,
            "myloasm": args.myloasm_marker,
            "metaflye": args.metaflye_marker,
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
