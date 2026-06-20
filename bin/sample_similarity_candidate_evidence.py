#!/usr/bin/env python3
"""Combine per-metric sample similarity candidates into an evidence table."""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

EXPECTED_METRICS = ("abund", "noabund")
EVIDENCE_COLUMNS = [
    "sample_min",
    "sample_max",
    "candidate_evidence",
    "abund_candidate",
    "abund_similarity",
    "abund_distance",
    "abund_threshold",
    "abund_similarity_method",
    "noabund_candidate",
    "noabund_similarity",
    "noabund_distance",
    "noabund_threshold",
    "noabund_similarity_method",
]
TRUE_VALUE = True
FALSE_VALUE = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine abundance-aware and no-abundance candidate TSVs",
    )
    parser.add_argument(
        "--candidates",
        nargs="+",
        required=True,
        help="Per-metric possible sample mix-up candidate TSVs",
    )
    parser.add_argument("--output", required=True, help="Output evidence TSV")
    return parser.parse_args()


def load_candidate(path: Path) -> pl.DataFrame:
    frame = pl.read_csv(path, separator="\t")
    if frame.is_empty():
        return pl.DataFrame(
            schema={
                "sample_min": pl.String,
                "sample_max": pl.String,
                "metric": pl.String,
                "similarity": pl.Float64,
                "distance": pl.Float64,
                "candidate_similarity_threshold": pl.Float64,
                "similarity_method": pl.String,
                "distance_definition": pl.String,
                "reason": pl.String,
            },
        )
    return frame


def metric_candidate_table(candidates: pl.DataFrame, metric: str) -> pl.DataFrame:
    metric_frame = candidates.filter(pl.col("metric") == metric)
    if metric_frame.is_empty():
        return pl.DataFrame(
            schema={
                "sample_min": pl.String,
                "sample_max": pl.String,
                f"{metric}_candidate": pl.Boolean,
                f"{metric}_similarity": pl.Float64,
                f"{metric}_distance": pl.Float64,
                f"{metric}_threshold": pl.Float64,
                f"{metric}_similarity_method": pl.String,
            },
        )
    return metric_frame.select(
        "sample_min",
        "sample_max",
        pl.lit(TRUE_VALUE).alias(f"{metric}_candidate"),
        pl.col("similarity").alias(f"{metric}_similarity"),
        pl.col("distance").alias(f"{metric}_distance"),
        pl.col("candidate_similarity_threshold").alias(f"{metric}_threshold"),
        pl.col("similarity_method").alias(f"{metric}_similarity_method"),
    )


def candidate_evidence(candidates: pl.DataFrame) -> pl.DataFrame:
    tables = [metric_candidate_table(candidates, metric) for metric in EXPECTED_METRICS]
    evidence = tables[0].join(
        tables[1],
        on=["sample_min", "sample_max"],
        how="full",
        coalesce=True,
    )
    if evidence.is_empty():
        return evidence.with_columns(
            pl.lit(FALSE_VALUE).alias("abund_candidate"),
            pl.lit(FALSE_VALUE).alias("noabund_candidate"),
            pl.lit("none").alias("candidate_evidence"),
        ).select(EVIDENCE_COLUMNS)

    evidence = evidence.with_columns(
        pl.col("abund_candidate").fill_null(FALSE_VALUE),
        pl.col("noabund_candidate").fill_null(FALSE_VALUE),
    )
    return (
        evidence.with_columns(
            pl.when(pl.col("abund_candidate") & pl.col("noabund_candidate"))
            .then(pl.lit("both"))
            .when(pl.col("abund_candidate"))
            .then(pl.lit("abund_only"))
            .when(pl.col("noabund_candidate"))
            .then(pl.lit("noabund_only"))
            .otherwise(pl.lit("none"))
            .alias("candidate_evidence"),
        )
        .select(EVIDENCE_COLUMNS)
        .sort(["candidate_evidence", "sample_min", "sample_max"])
    )


def main() -> int:
    args = parse_args()
    candidates = pl.concat(
        [load_candidate(Path(path)) for path in args.candidates],
        how="diagonal",
    )
    candidate_evidence(candidates).write_csv(args.output, separator="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
