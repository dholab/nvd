#!/usr/bin/env python3
"""Convert a sourmash compare CSV matrix to a pairwise distance table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import polars as pl

MIN_MATRIX_ROWS = 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert sourmash compare --csv output to pairwise distances",
    )
    parser.add_argument("--metric", required=True, help="Metric label, e.g. abund")
    parser.add_argument(
        "--matrix",
        required=True,
        help="CSV matrix from sourmash compare",
    )
    parser.add_argument("--output", required=True, help="Output pairwise TSV")
    return parser.parse_args()


def load_matrix(path: Path) -> pl.DataFrame:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))

    if len(rows) < MIN_MATRIX_ROWS:
        message = f"Matrix CSV {path} did not contain matrix rows"
        raise ValueError(message)

    sample_ids = rows[0]
    matrix_rows = rows[1:]

    if len(sample_ids) != len(matrix_rows):
        message = (
            f"Matrix CSV {path} is not square: header has {len(sample_ids)} names "
            f"but there are {len(matrix_rows)} rows"
        )
        raise ValueError(
            message,
        )

    pairwise_rows: list[dict[str, str | float]] = []
    for sample_a, similarities in zip(sample_ids, matrix_rows, strict=True):
        if len(similarities) != len(sample_ids):
            message = (
                f"Matrix CSV {path} row for {sample_a!r} has {len(similarities)} "
                f"columns, expected {len(sample_ids)}"
            )
            raise ValueError(
                message,
            )

        for sample_b, similarity in zip(sample_ids, similarities, strict=True):
            pairwise_rows.append(
                {
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "similarity": float(similarity),
                },
            )

    return pl.DataFrame(pairwise_rows)


def add_distance_columns(pairwise: pl.DataFrame, metric: str) -> pl.DataFrame:
    return pairwise.with_columns(
        pl.lit(metric).alias("metric"),
        (pl.lit(1.0) - pl.col("similarity")).alias("distance"),
        (pl.col("sample_a") == pl.col("sample_b")).alias("is_self_pair"),
    ).select(
        [
            "metric",
            "sample_a",
            "sample_b",
            "similarity",
            "distance",
            "is_self_pair",
        ],
    )


def main() -> int:
    args = parse_args()
    pairwise = load_matrix(Path(args.matrix)).pipe(add_distance_columns, args.metric)
    pairwise.sort(["sample_a", "sample_b"]).write_csv(args.output, separator="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
