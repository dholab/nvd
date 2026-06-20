#!/usr/bin/env python3
"""Report nearest sample-sketch neighbors and high-similarity mix-up candidates."""

from __future__ import annotations

import argparse

import polars as pl

TRUE_VALUE = True
FALSE_VALUE = False


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
        "--candidate-output",
        required=True,
        help="Output high-similarity candidate TSV",
    )
    parser.add_argument(
        "--top-neighbors",
        type=int,
        default=5,
        help="Number of non-self nearest neighbors to report for each sample",
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


def rank_neighbors(pairwise: pl.DataFrame) -> pl.DataFrame:
    return (
        non_self_pairs(pairwise)
        .sort(["sample_a", "similarity", "sample_b"], descending=[False, True, False])
        .with_columns(
            (pl.col("similarity").rank("ordinal", descending=True).over("sample_a"))
            .cast(pl.Int64)
            .alias("neighbor_rank"),
        )
    )


def add_reciprocal_nearest_neighbor(ranked: pl.DataFrame) -> pl.DataFrame:
    top_neighbors = ranked.filter(pl.col("neighbor_rank") == 1).select(
        pl.col("sample_a").alias("sample_b"),
        pl.col("sample_b").alias("sample_a"),
        pl.lit(TRUE_VALUE).alias("has_reciprocal_top_neighbor"),
    )
    return ranked.join(
        top_neighbors,
        left_on=["sample_a", "sample_b"],
        right_on=["sample_a", "sample_b"],
        how="left",
    ).with_columns(
        (
            (pl.col("neighbor_rank") == 1)
            & pl.col("has_reciprocal_top_neighbor").fill_null(FALSE_VALUE)
        ).alias(
            "is_reciprocal_nearest_neighbor",
        ),
    )


def nearest_neighbors(
    pairwise: pl.DataFrame,
    metric: str,
    top_neighbors: int,
) -> pl.DataFrame:
    return (
        rank_neighbors(pairwise)
        .pipe(add_reciprocal_nearest_neighbor)
        .filter(pl.col("neighbor_rank") <= top_neighbors)
        .rename({"sample_a": "sample_id"})
        .rename({"sample_b": "neighbor_sample_id"})
        .with_columns(pl.lit(metric).alias("metric"))
        .select(
            [
                "metric",
                "similarity_method",
                "distance_definition",
                "sample_id",
                "neighbor_rank",
                "neighbor_sample_id",
                "similarity",
                "distance",
                "is_reciprocal_nearest_neighbor",
            ],
        )
        .sort(["sample_id", "neighbor_rank", "neighbor_sample_id"])
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
            pl.lit("high_similarity_between_distinct_query_metagenome_sketches").alias(
                "reason",
            ),
        )
    )
    return candidates.select(
        [
            "metric",
            "similarity_method",
            "distance_definition",
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
    nearest_neighbors(pairwise, args.metric, args.top_neighbors).write_csv(
        args.nearest_neighbors_output,
        separator="\t",
    )
    possible_mixups(
        pairwise,
        args.metric,
        args.candidate_similarity_threshold,
    ).write_csv(args.candidate_output, separator="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
