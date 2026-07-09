"""Tests for FASTX profile histogram plotting."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING

from plot_fastx_profile import (
    histogram_chart,
    histogram_title,
    load_profile,
    main,
    quality_axis_title,
    read_histogram,
)

if TYPE_CHECKING:
    from pathlib import Path


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "stage", "bin_start", "bin_end", "sequence_count"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_profile(
    path: Path,
    *,
    thresholds: list[dict[str, object]] | None = None,
) -> None:
    path.write_text(
        json.dumps(
            {
                "sample_id": "sample-1",
                "stage": "preprocessed_single_read",
                "format": "fastq",
                "sequence_count": 10,
                "total_bases": 1000,
                "length": {"min": 50, "max": 150, "mean": 100.0, "bin_size": 50},
                "quality": {
                    "mean": 30.0,
                    "bin_size": 1,
                    "encoding": "phred33_integer_average",
                },
                "thresholds": thresholds or [],
            },
        )
        + "\n",
        encoding="utf-8",
    )


def test_plot_fastx_profile_writes_plots_next_to_profile(tmp_path: Path) -> None:
    profile = tmp_path / "sample-1.single_read.filtered.fastx_profile.json"
    length_histogram = tmp_path / "sample-1.single_read.filtered.length_histogram.tsv"
    quality_histogram = tmp_path / "sample-1.single_read.filtered.quality_histogram.tsv"
    write_profile(
        profile,
        thresholds=[
            {"name": "min_read_length", "axis": "length", "value": 75},
            {"name": "min_read_quality", "axis": "quality", "value": 30},
        ],
    )
    write_tsv(
        length_histogram,
        [
            {
                "sample_id": "sample-1",
                "stage": "preprocessed_single_read",
                "bin_start": 50,
                "bin_end": 99,
                "sequence_count": 4,
            },
            {
                "sample_id": "sample-1",
                "stage": "preprocessed_single_read",
                "bin_start": 100,
                "bin_end": 149,
                "sequence_count": 6,
            },
        ],
    )
    write_tsv(
        quality_histogram,
        [
            {
                "sample_id": "sample-1",
                "stage": "preprocessed_single_read",
                "bin_start": 30,
                "bin_end": 30,
                "sequence_count": 10,
            },
        ],
    )

    main(
        [
            "--profile-json",
            str(profile),
            "--length-histogram",
            str(length_histogram),
            "--quality-histogram",
            str(quality_histogram),
        ],
    )

    assert (tmp_path / "sample-1.single_read.filtered.length_histogram.png").exists()
    assert (tmp_path / "sample-1.single_read.filtered.quality_histogram.png").exists()
    assert not (
        tmp_path / "sample-1.single_read.filtered.length_histogram.html"
    ).exists()
    assert not (
        tmp_path / "sample-1.single_read.filtered.quality_histogram.html"
    ).exists()

    profile_summary = load_profile(profile)
    quality_frame = read_histogram(quality_histogram)
    quality_chart = histogram_chart(
        quality_frame,
        profile=profile_summary,
        axis="quality",
        x_title=quality_axis_title(profile_summary),
        title=histogram_title(profile_summary, "quality"),
    ).to_dict()
    assert quality_frame["plot_bin_end"].to_list() == [31]
    bar_layer = quality_chart["layer"][0]
    assert bar_layer["mark"]["color"] == "#9CA3AF"
    assert bar_layer["encoding"]["x2"]["field"] == "plot_bin_end"
    assert bar_layer["encoding"]["y2"]["datum"] == 0
    assert quality_chart["title"]["subtitle"] == "10 sequences · 1,000 bases"
    quality_spec = json.dumps(quality_chart)
    assert '"type": "rule"' in quality_spec
    assert '"type": "text"' in quality_spec
    assert "min_read_quality" in quality_spec
    assert (
        quality_axis_title(profile_summary) == "Integer average read quality (Phred+33)"
    )


def test_plot_fastx_profile_skips_quality_plot_when_quality_histogram_is_absent(
    tmp_path: Path,
) -> None:
    profile = tmp_path / "sample-1.contigs.fastx_profile.json"
    length_histogram = tmp_path / "sample-1.contigs.length_histogram.tsv"
    write_profile(profile)
    write_tsv(
        length_histogram,
        [
            {
                "sample_id": "sample-1",
                "stage": "filtered_contigs",
                "bin_start": 100,
                "bin_end": 149,
                "sequence_count": 1,
            },
        ],
    )

    main(
        [
            "--profile-json",
            str(profile),
            "--length-histogram",
            str(length_histogram),
        ],
    )

    assert (tmp_path / "sample-1.contigs.length_histogram.png").exists()
    assert not (tmp_path / "sample-1.contigs.quality_histogram.png").exists()
    assert not (tmp_path / "sample-1.contigs.quality_histogram.html").exists()
