#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///

"""Summarize samtools depth -aa rows into CRUMBS coverage support metrics."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import TextIO


OUTPUT_COLUMNS = [
    "sample_id",
    "qseqid",
    "contig_length",
    "covered_bases_1x",
    "breadth_1x",
    "raw_aligned_bases",
    "mean_depth_full",
    "median_depth_full",
    "median_depth_positive",
    "depth_p95",
    "depth_p99",
    "max_depth",
    "crumbs_p95",
    "crumbs_p99",
]
EXPECTED_DEPTH_COLUMNS = 3


class CoverageSummaryError(ValueError):
    """Raised when depth rows cannot be summarized safely."""


@dataclass
class DepthRow:
    qseqid: str
    position: int
    depth: int


def nearest_rank_percentile(values: Sequence[int], percentile: float) -> int:
    """Return a discrete empirical nearest-rank percentile value."""
    if not values:
        message = "cannot compute a percentile for an empty contig"
        raise CoverageSummaryError(message)
    if not 0 < percentile <= 1:
        message = "percentile must be in the interval (0, 1]"
        raise CoverageSummaryError(message)

    sorted_values = sorted(values)
    rank = math.ceil(percentile * len(sorted_values))
    return sorted_values[rank - 1]


def median(values: Sequence[int]) -> float:
    """Return the conventional median."""
    if not values:
        return 0.0

    sorted_values = sorted(values)
    midpoint = len(sorted_values) // 2
    if len(sorted_values) % 2 == 1:
        return float(sorted_values[midpoint])
    return (sorted_values[midpoint - 1] + sorted_values[midpoint]) / 2


def parse_depth_rows(lines: Iterable[str]) -> Iterable[DepthRow]:
    """Yield samtools depth -aa rows with contig, one-based position, and depth."""
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split("\t")
        if len(parts) != EXPECTED_DEPTH_COLUMNS:
            message = f"line {line_number}: expected 3 tab-delimited columns"
            raise CoverageSummaryError(message)

        qseqid, position_text, depth_text = parts
        try:
            position = int(position_text)
            depth = int(depth_text)
        except ValueError as exc:
            message = f"line {line_number}: position and depth must be integers"
            raise CoverageSummaryError(message) from exc

        if not qseqid:
            message = f"line {line_number}: contig name cannot be empty"
            raise CoverageSummaryError(message)
        if position < 1:
            message = f"line {line_number}: position must be one-based"
            raise CoverageSummaryError(message)
        if depth < 0:
            message = f"line {line_number}: depth cannot be negative"
            raise CoverageSummaryError(message)

        yield DepthRow(qseqid=qseqid, position=position, depth=depth)


def summarize_depths(
    sample_id: str,
    qseqid: str,
    depths: Sequence[int],
) -> dict[str, str]:
    """Summarize full-contig depths for one contig."""
    if not depths:
        message = f"contig {qseqid!r} has no depth rows"
        raise CoverageSummaryError(message)

    contig_length = len(depths)
    covered_bases = sum(1 for depth in depths if depth >= 1)
    raw_aligned_bases = sum(depths)
    positive_depths = [depth for depth in depths if depth > 0]
    depth_p95 = nearest_rank_percentile(depths, 0.95)
    depth_p99 = nearest_rank_percentile(depths, 0.99)

    return {
        "sample_id": sample_id,
        "qseqid": qseqid,
        "contig_length": str(contig_length),
        "covered_bases_1x": str(covered_bases),
        "breadth_1x": str(covered_bases / contig_length),
        "raw_aligned_bases": str(raw_aligned_bases),
        "mean_depth_full": str(raw_aligned_bases / contig_length),
        "median_depth_full": str(median(depths)),
        "median_depth_positive": str(median(positive_depths)),
        "depth_p95": str(depth_p95),
        "depth_p99": str(depth_p99),
        "max_depth": str(max(depths)),
        "crumbs_p95": str(sum(min(depth, depth_p95) for depth in depths)),
        "crumbs_p99": str(sum(min(depth, depth_p99) for depth in depths)),
    }


def summarize_rows(
    sample_id: str,
    rows: Iterable[DepthRow],
) -> Iterable[dict[str, str]]:
    """Summarize sorted samtools depth -aa rows into one output row per contig."""
    current_qseqid: str | None = None
    current_depths: list[int] = []
    previous_position = 0
    completed_qseqids: set[str] = set()

    for row in rows:
        if current_qseqid is None:
            current_qseqid = row.qseqid
        elif row.qseqid != current_qseqid:
            completed_qseqids.add(current_qseqid)
            yield summarize_depths(sample_id, current_qseqid, current_depths)
            current_qseqid = row.qseqid
            current_depths = []
            previous_position = 0

        if row.qseqid in completed_qseqids:
            message = f"contig {row.qseqid!r} appears in multiple blocks"
            raise CoverageSummaryError(message)
        if row.position <= previous_position:
            message = f"contig {row.qseqid!r} positions must increase monotonically"
            raise CoverageSummaryError(message)

        missing_positions = row.position - previous_position - 1
        if missing_positions > 0:
            current_depths.extend([0] * missing_positions)
        current_depths.append(row.depth)
        previous_position = row.position

    if current_qseqid is not None:
        yield summarize_depths(sample_id, current_qseqid, current_depths)


def write_summaries(rows: Iterable[dict[str, str]], output: TextIO) -> None:
    """Write coverage summary rows as a TSV."""
    writer = csv.DictWriter(output, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize samtools depth -aa output into CRUMBS coverage metrics",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument(
        "--depth-tsv",
        dest="depth_tsv",
        type=Path,
        default=Path("-"),
        help="samtools depth -aa TSV path, or '-' for stdin",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("-"),
        help="coverage summary TSV path, or '-' for stdout",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)

    try:
        if str(args.depth_tsv) == "-":
            rows = parse_depth_rows(sys.stdin)
            summaries = summarize_rows(args.sample_id, rows)
            if str(args.output) == "-":
                write_summaries(summaries, sys.stdout)
            else:
                with args.output.open("w", newline="", encoding="utf-8") as output:
                    write_summaries(summaries, output)
        else:
            with args.depth_tsv.open(encoding="utf-8") as input_handle:
                rows = parse_depth_rows(input_handle)
                summaries = summarize_rows(args.sample_id, rows)
                if str(args.output) == "-":
                    write_summaries(summaries, sys.stdout)
                else:
                    with args.output.open("w", newline="", encoding="utf-8") as output:
                        write_summaries(summaries, output)
    except (CoverageSummaryError, OSError) as exc:
        sys.stderr.write(f"{exc}\n")
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
