"""Tests for BLAST query-filtering summaries."""

from __future__ import annotations

import csv
import sys
from typing import TYPE_CHECKING
from unittest.mock import patch

import pytest
from summarize_blast_filtering import main, query_ids

if TYPE_CHECKING:
    from pathlib import Path


def test_query_ids_accepts_zero_byte_and_header_only_inputs(tmp_path: Path) -> None:
    empty = tmp_path / "empty.tsv"
    empty.touch()
    header_only = tmp_path / "header_only.tsv"
    header_only.write_text("sample\tqseqid\ttask\n", encoding="utf-8")

    assert query_ids(empty) == set()
    assert query_ids(header_only) == set()


def test_query_ids_resolves_reordered_qseqid_and_deduplicates_hits(
    tmp_path: Path,
) -> None:
    annotated = tmp_path / "annotated.tsv"
    annotated.write_text(
        "sample\ttask\tqseqid\n"
        "sample_A\tmegablast\tquery_001\n"
        "sample_A\tmegablast\tquery_001\n"
        "sample_A\tmegablast\tquery_002\n",
        encoding="utf-8",
    )

    assert query_ids(annotated) == {"query_001", "query_002"}


def test_query_ids_rejects_header_without_qseqid(tmp_path: Path) -> None:
    annotated = tmp_path / "annotated.tsv"
    annotated.write_text(
        "sample\ttask\nsample_A\tmegablast\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="qseqid"):
        query_ids(annotated)


def test_summary_counts_queries_by_named_qseqid(tmp_path: Path) -> None:
    annotated = tmp_path / "annotated.tsv"
    annotated.write_text(
        "sample\tqseqid\ttask\n"
        "sample_A\tquery_001\tmegablast\n"
        "sample_A\tquery_001\tmegablast\n"
        "sample_A\tquery_002\tmegablast\n",
        encoding="utf-8",
    )
    retained = tmp_path / "retained.tsv"
    retained.write_text(
        "task\tsample\tqseqid\nmegablast\tsample_A\tquery_001\n",
        encoding="utf-8",
    )
    summary = tmp_path / "summary.tsv"

    with patch.object(
        sys,
        "argv",
        [
            "summarize_blast_filtering.py",
            "--sample-id",
            "sample_A",
            "--query-class",
            "single_read",
            "--stage",
            "megablast_virus_filter",
            "--input",
            str(annotated),
            "--output",
            str(retained),
            "--summary",
            str(summary),
        ],
    ):
        main()

    with summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert rows == [
        {
            "sample_id": "sample_A",
            "query_class": "single_read",
            "stage": "megablast_virus_filter",
            "queries_in": "2",
            "queries_retained": "1",
            "queries_removed": "1",
        },
    ]
