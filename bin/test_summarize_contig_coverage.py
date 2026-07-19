"""Tests for CRUMBS contig coverage summarization."""

from __future__ import annotations

import csv
import sys
from io import StringIO
from typing import TYPE_CHECKING

import pytest
from summarize_contig_coverage import main

if TYPE_CHECKING:
    from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_samtools_depth_rows_produce_full_contig_crumbs(tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text(
        "contig-a\t1\t0\ncontig-a\t2\t0\ncontig-a\t3\t2\ncontig-a\t4\t4\n",
        encoding="utf-8",
    )
    output = tmp_path / "coverage.tsv"

    main(
        ["--sample-id", "sample-1", "--depth-tsv", str(depth), "--output", str(output)],
    )

    rows = read_rows(output)
    assert rows == [
        {
            "sample_id": "sample-1",
            "qseqid": "contig-a",
            "contig_length": "4",
            "covered_bases_1x": "2",
            "breadth_1x": "0.5",
            "raw_aligned_bases": "6",
            "mean_depth_full": "1.5",
            "median_depth_full": "1.0",
            "median_depth_positive": "3.0",
            "depth_p95": "4",
            "depth_p99": "4",
            "max_depth": "4",
            "crumbs_p95": "6",
            "crumbs_p99": "6",
        },
    ]


def test_nearest_rank_percentiles_include_zero_depth_positions(tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    zero_rows = [f"contig-a\t{position}\t0" for position in range(1, 100)]
    depth.write_text(
        "\n".join([*zero_rows, "contig-a\t100\t100"]),
        encoding="utf-8",
    )
    output = tmp_path / "coverage.tsv"

    main(
        ["--sample-id", "sample-1", "--depth-tsv", str(depth), "--output", str(output)],
    )

    [row] = read_rows(output)
    assert row["depth_p95"] == "0"
    assert row["depth_p99"] == "0"
    assert row["crumbs_p95"] == "0"
    assert row["crumbs_p99"] == "0"
    assert row["max_depth"] == "100"


def test_missing_positions_are_treated_as_zero_depth(tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text(
        "contig-a\t1\t2\ncontig-a\t3\t4\n",
        encoding="utf-8",
    )
    output = tmp_path / "coverage.tsv"

    main(
        ["--sample-id", "sample-1", "--depth-tsv", str(depth), "--output", str(output)],
    )

    [row] = read_rows(output)
    assert row["contig_length"] == "3"
    assert row["covered_bases_1x"] == "2"
    assert row["breadth_1x"] == "0.6666666666666666"
    assert row["median_depth_full"] == "2.0"


def test_all_zero_depth_contig_is_reported_with_zero_crumbs(tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text(
        "contig-a\t1\t0\ncontig-a\t2\t0\n",
        encoding="utf-8",
    )
    output = tmp_path / "coverage.tsv"

    main(
        ["--sample-id", "sample-1", "--depth-tsv", str(depth), "--output", str(output)],
    )

    [row] = read_rows(output)
    assert row["covered_bases_1x"] == "0"
    assert row["breadth_1x"] == "0.0"
    assert row["median_depth_positive"] == "0.0"
    assert row["crumbs_p95"] == "0"
    assert row["crumbs_p99"] == "0"


def test_depth_rows_can_stream_from_stdin(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sys, "stdin", StringIO("contig-a\t1\t3\n"))
    output = tmp_path / "coverage.tsv"

    main(["--sample-id", "sample-1", "--depth-tsv", "-", "--output", str(output)])

    [row] = read_rows(output)
    assert row["qseqid"] == "contig-a"
    assert row["crumbs_p99"] == "3"


def test_repeated_contig_blocks_are_rejected(tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text(
        "contig-a\t1\t1\ncontig-b\t1\t1\ncontig-a\t2\t1\n",
        encoding="utf-8",
    )
    output = tmp_path / "coverage.tsv"

    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "--sample-id",
                "sample-1",
                "--depth-tsv",
                str(depth),
                "--output",
                str(output),
            ],
        )

    assert exc_info.value.code == 1
