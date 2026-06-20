#!/usr/bin/env python3
"""Report nearest sample-sketch neighbors and high-similarity mix-up candidates."""

from __future__ import annotations

import argparse

import polars as pl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build sample similarity QC reports from pairwise sketch distances",
    )
    parser.add_argument("--metric", required=True, help="Metric label, e.g. abund")
    parser.add_argument("--pairwise", required=True, help="Pairwise distance TSV")
    parser.add_argument(
        "--nearest-neighbors-output",
        required=True,
        help="Output nearest-neighbor TSV",
    )
    parser.add_argument(
        "--possible-mixups-output",
        required=True,
        help="Output high-similarity candidate TSV",
    )
    parser.add_argument(
        "--candidate-similarity-threshold",
        type=float,
        default=0.95,
        help="Minimum non-self similarity to report as a possible mix-up candidate",
    )
    return parser.parse_args()


def non_self_pairs(pairwise: pl.DataFrame) -> pl.DataFrame:
    return pairwise.filter(~pl.col("is_self_pair"))


def nearest_neighbors(pairwise: pl.DataFrame, metric: str) -> pl.DataFrame:
    return (
        non_self_pairs(pairwise)
        .sort(["sample_a", "similarity", "sample_b"], descending=[False, True, False])
        .group_by("sample_a")
        .agg(
            pl.first("sample_b").alias("nearest_sample_id"),
            pl.first("similarity").alias("nearest_similarity"),
            pl.first("distance").alias("nearest_distance"),
        )
        .rename({"sample_a": "sample_id"})
        .with_columns(pl.lit(metric).alias("metric"))
        .select(
            [
                "metric",
                "sample_id",
                "nearest_sample_id",
                "nearest_similarity",
                "nearest_distance",
            ],
        )
        .sort("sample_id")
    )


def possible_mixups(
    pairwise: pl.DataFrame,
    metric: str,
    threshold: float,
) -> pl.DataFrame:
    candidates = (
        non_self_pairs(pairwise)
        .with_columns(
            pl.min_horizontal(["sample_a", "sample_b"]).alias("sample_min"),
            pl.max_horizontal(["sample_a", "sample_b"]).alias("sample_max"),
        )
        .unique(subset=["sample_min", "sample_max"], keep="first")
        .filter(pl.col("similarity") >= threshold)
        .with_columns(
            pl.lit(metric).alias("metric"),
            pl.lit(threshold).alias("candidate_similarity_threshold"),
            pl.lit("high_sketch_similarity_between_distinct_samples").alias("reason"),
        )
    )
    return candidates.select(
        [
            "metric",
            "sample_min",
            "sample_max",
            "similarity",
            "distance",
            "candidate_similarity_threshold",
            "reason",
        ],
    ).sort(["similarity", "sample_min", "sample_max"], descending=[True, False, False])


def main() -> int:
    args = parse_args()
    pairwise = pl.read_csv(args.pairwise, separator="\t")
    nearest_neighbors(pairwise, args.metric).write_csv(
        args.nearest_neighbors_output,
        separator="\t",
    )
    possible_mixups(
        pairwise,
        args.metric,
        args.candidate_similarity_threshold,
    ).write_csv(args.possible_mixups_output, separator="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
