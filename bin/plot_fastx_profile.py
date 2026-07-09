#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#   "altair>=6.1.0",
#   "polars>=1.27.1",
#   "vl-convert-python>=1.9.0.post1",
# ]
# ///
"""Plot FASTX profile length and quality histograms."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import altair as alt
import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any


@dataclass(frozen=True, slots=True)
class ProfileSummary:
    sample_id: str
    stage: str
    sequence_count: int
    total_bases: int
    format: str
    quality_encoding: str | None
    thresholds: tuple[dict[str, Any], ...]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot FASTX profile histograms.")
    parser.add_argument("--profile-json", type=Path, required=True)
    parser.add_argument("--length-histogram", type=Path)
    parser.add_argument("--quality-histogram", type=Path)
    parser.add_argument("--output-prefix")
    args = parser.parse_args(argv)
    if args.length_histogram is None and args.quality_histogram is None:
        parser.error(
            "At least one of --length-histogram or --quality-histogram is required.",
        )
    return args


def load_profile(path: Path) -> ProfileSummary:
    profile = json.loads(path.read_text(encoding="utf-8"))
    return ProfileSummary(
        sample_id=profile["sample_id"],
        stage=profile["stage"],
        sequence_count=profile["sequence_count"],
        total_bases=profile["total_bases"],
        format=profile["format"],
        quality_encoding=profile["quality"].get("encoding"),
        thresholds=tuple(profile.get("thresholds", [])),
    )


def read_histogram(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t").with_columns(
        pl.col("bin_start").cast(pl.Int64),
        pl.col("bin_end").cast(pl.Int64),
        pl.col("sequence_count").cast(pl.Int64),
        (pl.col("bin_end").cast(pl.Int64) + 1).alias("plot_bin_end"),
    )


def threshold_frame(profile: ProfileSummary, axis: str) -> pl.DataFrame:
    rows = [
        threshold for threshold in profile.thresholds if threshold.get("axis") == axis
    ]
    if not rows:
        return pl.DataFrame(
            schema={"name": pl.String, "value": pl.Float64, "label": pl.String},
        )
    return pl.DataFrame(
        {
            "name": [row.get("name", "threshold") for row in rows],
            "value": [row["value"] for row in rows],
            "label": [
                f"{row.get('name', 'threshold')}: {float(row['value']):g}"
                for row in rows
            ],
        },
    )


def histogram_title(profile: ProfileSummary, label: str) -> str:
    return f"{profile.sample_id} {profile.stage} {label} histogram"


def profile_subtitle(profile: ProfileSummary) -> str:
    return f"{profile.sequence_count:,} sequences · {profile.total_bases:,} bases"


def quality_axis_title(profile: ProfileSummary) -> str:
    if profile.quality_encoding == "phred33_integer_average":
        return "Integer average read quality (Phred+33)"
    return "Per-read mean quality (Phred+33)"


def base_histogram_chart(
    frame: pl.DataFrame,
    *,
    x_title: str,
) -> alt.Chart:
    return (
        alt.Chart(frame)
        .mark_bar(color="#9CA3AF", stroke="#6B7280", strokeWidth=0.5)
        .encode(
            x=alt.X("bin_start:Q", title=x_title),
            x2="plot_bin_end:Q",
            y=alt.Y("sequence_count:Q", title="Sequences"),
            y2=alt.Y2(datum=0),
            tooltip=[
                alt.Tooltip("sample_id:N", title="Sample"),
                alt.Tooltip("stage:N", title="Stage"),
                alt.Tooltip("bin_start:Q", title="Bin start"),
                alt.Tooltip("bin_end:Q", title="Bin end"),
                alt.Tooltip("sequence_count:Q", title="Sequences", format=","),
            ],
        )
        .properties(
            width=720,
            height=360,
        )
    )


def threshold_rules(thresholds: pl.DataFrame) -> alt.Chart:
    rules = (
        alt.Chart(thresholds)
        .mark_rule(color="#d62728", strokeDash=[6, 4], size=2)
        .encode(
            x="value:Q",
            tooltip=[
                alt.Tooltip("name:N", title="Threshold"),
                alt.Tooltip("value:Q", title="Value"),
            ],
        )
    )
    labels = (
        alt.Chart(thresholds)
        .mark_text(align="left", baseline="top", dx=5, dy=5, color="#d62728")
        .encode(
            x="value:Q",
            y=alt.value(8),
            text="label:N",
            tooltip=[
                alt.Tooltip("name:N", title="Threshold"),
                alt.Tooltip("value:Q", title="Value"),
            ],
        )
    )
    return rules + labels


def empty_quality_chart(profile: ProfileSummary) -> alt.Chart:
    return (
        alt.Chart(
            pl.DataFrame(
                {
                    "message": [
                        "No FASTQ quality scores were available for this FASTX profile.",
                    ],
                },
            ),
        )
        .mark_text(size=16)
        .encode(text="message:N")
        .properties(
            title=histogram_title(profile, "quality"),
            width=720,
            height=180,
        )
    )


def histogram_chart(
    frame: pl.DataFrame,
    *,
    profile: ProfileSummary,
    axis: str,
    x_title: str,
    title: str,
) -> alt.Chart:
    if frame.is_empty() and axis == "quality":
        return empty_quality_chart(profile)
    chart = base_histogram_chart(
        frame,
        x_title=x_title,
    )
    thresholds = threshold_frame(profile, axis)
    layered = chart if thresholds.is_empty() else chart + threshold_rules(thresholds)
    return layered.properties(
        title=alt.TitleParams(text=title, subtitle=profile_subtitle(profile)),
        width=720,
        height=360,
    )


def output_prefix(profile_json: Path, explicit_prefix: str | None) -> Path:
    if explicit_prefix:
        return Path(explicit_prefix)
    suffix = ".fastx_profile.json"
    return profile_json.with_name(
        profile_json.name[: -len(suffix)]
        if profile_json.name.endswith(suffix)
        else profile_json.stem,
    )


def save_chart(chart: alt.Chart, stem: Path) -> None:
    chart.save(Path(f"{stem}.png"))


def write_plots(
    *,
    profile: ProfileSummary,
    length_histogram: pl.DataFrame | None,
    quality_histogram: pl.DataFrame | None,
    prefix: Path,
) -> None:
    if length_histogram is not None:
        save_chart(
            histogram_chart(
                length_histogram,
                profile=profile,
                axis="length",
                x_title="Sequence length (bp)",
                title=histogram_title(profile, "length"),
            ),
            prefix.with_name(f"{prefix.name}.length_histogram"),
        )
    if quality_histogram is not None:
        save_chart(
            histogram_chart(
                quality_histogram,
                profile=profile,
                axis="quality",
                x_title=quality_axis_title(profile),
                title=histogram_title(profile, "quality"),
            ),
            prefix.with_name(f"{prefix.name}.quality_histogram"),
        )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    write_plots(
        profile=load_profile(args.profile_json),
        length_histogram=read_histogram(args.length_histogram)
        if args.length_histogram
        else None,
        quality_histogram=read_histogram(args.quality_histogram)
        if args.quality_histogram
        else None,
        prefix=output_prefix(args.profile_json, args.output_prefix),
    )


if __name__ == "__main__":
    main()
