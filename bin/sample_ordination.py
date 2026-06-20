#!/usr/bin/env python3
"""Run classical PCoA from pairwise sample sketch distances."""

from __future__ import annotations

import argparse

import numpy as np
import polars as pl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a distance matrix and PCoA coordinates from pairwise TSV",
    )
    parser.add_argument("--metric", required=True, help="Metric label, e.g. abund")
    parser.add_argument("--pairwise", required=True, help="Pairwise distance TSV")
    parser.add_argument("--distance-output", required=True, help="Output matrix TSV")
    parser.add_argument(
        "--coordinates-output",
        required=True,
        help="Output coordinates TSV",
    )
    parser.add_argument("--variance-output", required=True, help="Output variance TSV")
    return parser.parse_args()


def sample_ids_from_pairwise(pairwise: pl.DataFrame) -> list[str]:
    return (
        pairwise.select(pl.concat_list(["sample_a", "sample_b"]).alias("sample_id"))
        .explode("sample_id")
        .unique(maintain_order=True)
        .get_column("sample_id")
        .to_list()
    )


def build_distance_matrix(
    pairwise: pl.DataFrame,
    sample_ids: list[str],
) -> pl.DataFrame:
    matrix = pairwise.select(["sample_a", "sample_b", "distance"]).pivot(
        on="sample_b",
        index="sample_a",
        values="distance",
    )

    sample_order = pl.DataFrame(
        {"sample_a": sample_ids, "row_order": list(range(len(sample_ids)))},
    )

    matrix = (
        matrix.join(sample_order, on="sample_a", how="inner")
        .sort("row_order")
        .drop("row_order")
        .select(["sample_a", *sample_ids])
    )

    diagonal_fixed = [
        pl.when(pl.col("sample_a") == sample_id)
        .then(pl.lit(0.0))
        .otherwise(pl.col(sample_id))
        .alias(sample_id)
        for sample_id in sample_ids
    ]
    return matrix.with_columns(diagonal_fixed)


def validate_matrix(distance_matrix: pl.DataFrame, sample_ids: list[str]) -> np.ndarray:
    values = distance_matrix.select(sample_ids).to_numpy()
    expected_shape = (len(sample_ids), len(sample_ids))
    if values.shape != expected_shape:
        message = f"Distance matrix shape {values.shape} != {expected_shape}"
        raise ValueError(message)
    if np.isnan(values).any():
        message = "Distance matrix contains NaN values"
        raise ValueError(message)
    if not np.allclose(np.diag(values), 0.0):
        message = "Distance matrix diagonal contains non-zero values"
        raise ValueError(message)
    if not np.allclose(values, values.T, atol=1e-12):
        message = "Distance matrix is not symmetric"
        raise ValueError(message)
    return values


def run_pcoa(distance_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_samples = distance_values.shape[0]
    if n_samples == 0:
        return np.empty((0, 0)), np.empty((0,))
    if n_samples == 1:
        return np.zeros((1, 1)), np.array([1.0])

    squared = distance_values**2
    centering = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    gram = -0.5 * centering @ squared @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive = eigenvalues > 0
    positive_values = eigenvalues[positive]
    if positive_values.size == 0:
        return np.zeros((n_samples, 1)), np.array([0.0])

    coordinates = eigenvectors[:, positive] * np.sqrt(positive_values)
    proportions = positive_values / positive_values.sum()
    return coordinates, proportions


def coordinates_table(
    metric: str,
    sample_ids: list[str],
    coordinates: np.ndarray,
) -> pl.DataFrame:
    data: dict[str, list[str] | list[float]] = {"sample_id": sample_ids}
    for axis_index in range(coordinates.shape[1]):
        data[f"axis_{axis_index + 1}"] = coordinates[:, axis_index].tolist()
    return (
        pl.DataFrame(data)
        .with_columns(pl.lit(metric).alias("metric"))
        .select(
            [
                "metric",
                "sample_id",
                *[f"axis_{i + 1}" for i in range(coordinates.shape[1])],
            ],
        )
    )


def variance_table(metric: str, proportions: np.ndarray) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "metric": [metric] * len(proportions),
            "axis": [f"axis_{index}" for index in range(1, len(proportions) + 1)],
            "proportion_explained": proportions.tolist(),
        },
    )


def main() -> int:
    args = parse_args()
    pairwise = pl.read_csv(args.pairwise, separator="\t")
    sample_ids = sample_ids_from_pairwise(pairwise)
    distance_matrix = build_distance_matrix(pairwise, sample_ids)
    distance_values = validate_matrix(distance_matrix, sample_ids)
    coordinates, proportions = run_pcoa(distance_values)

    distance_matrix.rename({"sample_a": "sample_id"}).write_csv(
        args.distance_output,
        separator="\t",
    )
    coordinates_table(args.metric, sample_ids, coordinates).write_csv(
        args.coordinates_output,
        separator="\t",
    )
    variance_table(args.metric, proportions).write_csv(
        args.variance_output,
        separator="\t",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
