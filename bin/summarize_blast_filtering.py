#!/usr/bin/env python3
"""Summarize query retention across a BLAST taxonomic filter."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--query-class", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def query_ids(path: Path) -> set[str]:
    with path.open(encoding="utf-8") as handle:
        return {
            fields[0]
            for line in handle
            if line.strip() and (fields := line.rstrip("\n").split("\t"))
        }


def main() -> None:
    args = parse_args()
    input_ids = query_ids(args.input)
    output_ids = query_ids(args.output)
    with args.summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "sample_id",
                "query_class",
                "stage",
                "queries_in",
                "queries_retained",
                "queries_removed",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": args.sample_id,
                "query_class": args.query_class,
                "stage": args.stage,
                "queries_in": len(input_ids),
                "queries_retained": len(output_ids),
                "queries_removed": len(input_ids - output_ids),
            },
        )


if __name__ == "__main__":
    main()
