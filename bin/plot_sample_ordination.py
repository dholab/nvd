#!/usr/bin/env python3
"""Plot sample sketch PCoA ordination."""

from __future__ import annotations

import argparse

import altair as alt
import polars as pl

METRIC_TITLES = {
    "abund": "Ordination of abundance-weighted FracMinHash angular distances",
    "noabund": "Ordination of FracMinHash Jaccard distances",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot sample sketch PCoA ordination")
    parser.add_argument("--metric", required=True, help="Metric label, e.g. abund")
    parser.add_argument(
        "--coordinates",
        required=True,
        help="Ordination coordinates TSV",
    )
    parser.add_argument("--variance", required=True, help="Ordination variance TSV")
    parser.add_argument("--html-output", required=True, help="Output HTML plot")
    parser.add_argument("--png-output", required=True, help="Output PNG plot")
    return parser.parse_args()


def axis_label(variance: pl.DataFrame, axis: str) -> str:
    row = variance.filter(pl.col("axis") == axis)
    if row.is_empty():
        return axis.replace("_", " ").title()
    proportion = row.get_column("proportion_explained").item()
    return f"{axis.replace('_', ' ').title()} ({proportion * 100:.1f}%)"


def build_chart(
    metric: str,
    coordinates: pl.DataFrame,
    variance: pl.DataFrame,
) -> alt.Chart:
    if "axis_2" not in coordinates.columns:
        coordinates = coordinates.with_columns(pl.lit(0.0).alias("axis_2"))

    return (
        alt.Chart(alt.InlineData(values=coordinates.to_dicts()))
        .mark_circle(size=120, opacity=0.85)
        .encode(
            x=alt.X("axis_1:Q", title=axis_label(variance, "axis_1")),
            y=alt.Y("axis_2:Q", title=axis_label(variance, "axis_2")),
            tooltip=["sample_id:N", "axis_1:Q", "axis_2:Q"],
        )
        .properties(
            title=METRIC_TITLES.get(
                metric,
                f"Ordination of sample sketch distances ({metric})",
            ),
            width=720,
            height=520,
        )
    )


def main() -> int:
    args = parse_args()
    coordinates = pl.read_csv(args.coordinates, separator="\t")
    variance = pl.read_csv(args.variance, separator="\t")
    chart = build_chart(args.metric, coordinates, variance)
    chart.save(args.html_output)
    chart.save(args.png_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
