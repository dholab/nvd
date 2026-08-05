#!/usr/bin/env python3

import argparse
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import altair as alt
import polars as pl

SUMMARY_SUFFIX = ".deacon_filter.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create target enrichment summary tables and visualizations.",
    )
    parser.add_argument(
        "--summaries",
        nargs="+",
        required=True,
        help="Deacon summary JSON files from target read enrichment",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()


def sample_id_from_path(path: Path) -> str:
    if path.name.endswith(SUMMARY_SUFFIX):
        return path.name[: -len(SUMMARY_SUFFIX)]
    return path.stem


def read_summary(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        summary = json.load(handle)

    return {
        "sample_id": sample_id_from_path(path),
        "version": summary["version"],
        "index": summary.get("index", ""),
        "k": summary["k"],
        "w": summary["w"],
        "abs_threshold": summary["abs_threshold"],
        "rel_threshold": summary["rel_threshold"],
        "seqs_in": summary["seqs_in"],
        "seqs_enriched": summary["seqs_out"],
        "seqs_not_enriched": summary["seqs_removed"],
        "seqs_enriched_proportion": summary["seqs_out_proportion"],
        "seqs_not_enriched_proportion": summary["seqs_removed_proportion"],
        "bp_in": summary["bp_in"],
        "bp_enriched": summary["bp_out"],
        "bp_not_enriched": summary["bp_removed"],
        "bp_enriched_proportion": summary["bp_out_proportion"],
        "bp_not_enriched_proportion": summary["bp_removed_proportion"],
        "time_seconds": summary["time"],
        "seqs_per_second": summary["seqs_per_second"],
        "bp_per_second": summary["bp_per_second"],
    }


def load_target_enrichment_summaries(paths: Sequence[Path]) -> pl.DataFrame:
    if not paths:
        message = "At least one target enrichment summary JSON is required"
        raise ValueError(message)

    return pl.DataFrame([read_summary(path) for path in sorted(paths)]).sort(
        "bp_enriched_proportion",
        descending=True,
    )


def add_percent_columns(frame: pl.DataFrame) -> pl.DataFrame:
    return frame.with_columns(
        (pl.col("seqs_enriched_proportion") * 100).alias("seqs_enriched_percent"),
        (pl.col("bp_enriched_proportion") * 100).alias("bp_enriched_percent"),
        (pl.col("seqs_not_enriched_proportion") * 100).alias(
            "seqs_not_enriched_percent",
        ),
        (pl.col("bp_not_enriched_proportion") * 100).alias("bp_not_enriched_percent"),
    )


def write_tsv(frame: pl.DataFrame, path: str | Path) -> pl.DataFrame:
    frame.write_csv(path, separator="\t")
    return frame


def ranked_enriched_bases_chart(frame: pl.DataFrame) -> alt.Chart:
    return (
        alt.Chart(frame)
        .mark_bar()
        .encode(
            x=alt.X(
                "bp_enriched_percent:Q",
                title="Target-enriched bases retained (%)",
            ),
            y=alt.Y("sample_id:N", sort="-x", title="Sample"),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip(
                    "bp_enriched_percent:Q",
                    title="Bases retained by enrichment (%)",
                    format=".3f",
                ),
                alt.Tooltip(
                    "seqs_enriched_percent:Q",
                    title="Reads retained by enrichment (%)",
                    format=".3f",
                ),
                alt.Tooltip(
                    "seqs_enriched:Q",
                    title="Target-enriched reads",
                    format=",",
                ),
                alt.Tooltip(
                    "seqs_not_enriched:Q",
                    title="Reads not retained",
                    format=",",
                ),
                alt.Tooltip("seqs_in:Q", title="Input reads", format=","),
                alt.Tooltip("bp_enriched:Q", title="Target-enriched bases", format=","),
                alt.Tooltip(
                    "bp_not_enriched:Q",
                    title="Bases not retained",
                    format=",",
                ),
                alt.Tooltip("bp_in:Q", title="Input bases", format=","),
            ],
        )
        .properties(
            title="Target enrichment by sample",
            width=700,
            height=alt.Step(22),
        )
    )


def retained_filtered_long(frame: pl.DataFrame) -> pl.DataFrame:
    return (
        frame.select(
            "sample_id",
            "seqs_in",
            "seqs_enriched",
            "seqs_not_enriched",
            "bp_in",
            "bp_enriched",
            "bp_not_enriched",
            pl.col("bp_enriched_percent").alias("retained_by_enrichment"),
            pl.col("bp_not_enriched_percent").alias("not_retained_by_enrichment"),
        )
        .unpivot(
            on=["retained_by_enrichment", "not_retained_by_enrichment"],
            index=[
                "sample_id",
                "seqs_in",
                "seqs_enriched",
                "seqs_not_enriched",
                "bp_in",
                "bp_enriched",
                "bp_not_enriched",
            ],
            variable_name="fraction",
            value_name="percent",
        )
        .with_columns(
            pl.when(pl.col("fraction") == "retained_by_enrichment")
            .then(pl.lit("Retained by target enrichment"))
            .otherwise(pl.lit("Not retained by target enrichment"))
            .alias("fraction_label"),
        )
    )


def retained_vs_filtered_chart(frame: pl.DataFrame) -> alt.Chart:
    long_frame = retained_filtered_long(frame)

    return (
        alt.Chart(long_frame)
        .mark_bar()
        .encode(
            x=alt.X("percent:Q", stack="normalize", title="Fraction of input bases"),
            y=alt.Y("sample_id:N", sort=frame["sample_id"].to_list(), title="Sample"),
            color=alt.Color("fraction_label:N", title="Classification"),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip("fraction_label:N", title="Classification"),
                alt.Tooltip("percent:Q", title="Bases (%)", format=".3f"),
                alt.Tooltip(
                    "seqs_enriched:Q",
                    title="Target-enriched reads",
                    format=",",
                ),
                alt.Tooltip(
                    "seqs_not_enriched:Q",
                    title="Reads not retained",
                    format=",",
                ),
                alt.Tooltip("seqs_in:Q", title="Input reads", format=","),
                alt.Tooltip("bp_enriched:Q", title="Target-enriched bases", format=","),
                alt.Tooltip(
                    "bp_not_enriched:Q",
                    title="Bases not retained",
                    format=",",
                ),
                alt.Tooltip("bp_in:Q", title="Input bases", format=","),
            ],
        )
        .properties(
            title="Target-enriched vs not-retained bases",
            width=700,
            height=alt.Step(22),
        )
    )


def reads_vs_bases_scatter_chart(frame: pl.DataFrame) -> alt.Chart:
    points = (
        alt.Chart(frame)
        .mark_circle(size=80)
        .encode(
            x=alt.X(
                "seqs_enriched_percent:Q",
                title="Reads retained by target enrichment (%)",
            ),
            y=alt.Y(
                "bp_enriched_percent:Q",
                title="Bases retained by target enrichment (%)",
            ),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip(
                    "seqs_enriched_percent:Q",
                    title="Reads retained by enrichment (%)",
                    format=".3f",
                ),
                alt.Tooltip(
                    "bp_enriched_percent:Q",
                    title="Bases retained by enrichment (%)",
                    format=".3f",
                ),
                alt.Tooltip(
                    "seqs_enriched:Q",
                    title="Target-enriched reads",
                    format=",",
                ),
                alt.Tooltip(
                    "seqs_not_enriched:Q",
                    title="Reads not retained",
                    format=",",
                ),
                alt.Tooltip("seqs_in:Q", title="Input reads", format=","),
                alt.Tooltip("bp_enriched:Q", title="Target-enriched bases", format=","),
                alt.Tooltip(
                    "bp_not_enriched:Q",
                    title="Bases not retained",
                    format=",",
                ),
                alt.Tooltip("bp_in:Q", title="Input bases", format=","),
            ],
        )
    )
    diagonal = (
        alt.Chart(pl.DataFrame({"x": [0, 100], "y": [0, 100]}))
        .mark_line(strokeDash=[4, 4], color="gray")
        .encode(x="x:Q", y="y:Q")
    )

    return (diagonal + points).properties(
        title="Read-count vs base-count target enrichment",
        width=550,
        height=450,
    )


def save_chart(chart: alt.Chart, stem: Path) -> None:
    chart.save(stem.with_suffix(".html"))
    chart.save(stem.with_suffix(".png"))


def write_visualizations(frame: pl.DataFrame, outdir: Path) -> None:
    save_chart(
        ranked_enriched_bases_chart(frame),
        outdir / "target_enriched_bases_ranked",
    )
    save_chart(
        retained_vs_filtered_chart(frame),
        outdir / "target_retained_vs_filtered_stacked",
    )
    save_chart(
        reads_vs_bases_scatter_chart(frame),
        outdir / "target_reads_vs_bases_scatter",
    )


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = (
        load_target_enrichment_summaries([Path(path) for path in args.summaries])
        .pipe(add_percent_columns)
        .pipe(write_tsv, path=outdir / "target_enrichment_summary.tsv")
    )
    write_visualizations(summary, outdir)


if __name__ == "__main__":
    main()
