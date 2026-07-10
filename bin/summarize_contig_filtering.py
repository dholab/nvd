#!/usr/bin/env python3
"""Render Deacon contig-filter summaries as a routing decision table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--target-summary", type=Path, required=True)
    parser.add_argument("--depletion-summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def summary_row(sample_id: str, phase: str, path: Path) -> dict[str, object]:
    summary = json.loads(path.read_text(encoding="utf-8"))
    return {
        "sample_id": sample_id,
        "phase": phase,
        "sequences_in": summary["seqs_in"],
        "sequences_out": summary["seqs_out"],
        "sequences_removed": summary["seqs_removed"],
        "bases_in": summary["bp_in"],
        "bases_out": summary["bp_out"],
        "bases_removed": summary["bp_removed"],
    }


def main() -> None:
    args = parse_args()
    rows = [summary_row(args.sample_id, "target_filter", args.target_summary)]
    if args.depletion_summary is not None:
        rows.append(
            summary_row(args.sample_id, "host_depletion", args.depletion_summary),
        )
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
