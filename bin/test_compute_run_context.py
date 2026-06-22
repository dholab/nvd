"""Tests for state-free run-context computation."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

import pytest
from compute_run_context import read_sample_ids
from py_nvd.read_inputs import ResolutionError
from py_nvd.run_context import compute_sample_set_id

if TYPE_CHECKING:
    from pathlib import Path


def touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")
    return path


def write_samplesheet(path: Path, rows: list[dict[str, str]]) -> Path:
    fieldnames = [
        "sample_id",
        "srr",
        "platform",
        "fastq1",
        "fastq2",
        "fastq1_glob",
        "fastq2_glob",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_read_sample_ids_uses_resolved_samplesheet_semantics(tmp_path: Path) -> None:
    r1 = touch(tmp_path / "sample_R1.fastq.gz")
    r2 = touch(tmp_path / "sample_R2.fastq.gz")
    samplesheet = write_samplesheet(
        tmp_path / "samples.csv",
        [
            {
                "sample_id": "S1",
                "platform": "illumina",
                "fastq1": str(r1),
                "fastq2": str(r2),
            },
            {
                "sample_id": "S2",
                "platform": "ont",
                "srr": "SRR123456",
            },
        ],
    )

    assert read_sample_ids(samplesheet) == ["S1", "S2"]


def test_read_sample_ids_rejects_duplicate_samples_like_read_resolution(
    tmp_path: Path,
) -> None:
    samplesheet = write_samplesheet(
        tmp_path / "samples.csv",
        [
            {"sample_id": "S1", "platform": "ont", "srr": "SRR1"},
            {"sample_id": "S1", "platform": "ont", "srr": "SRR2"},
        ],
    )

    with pytest.raises(ResolutionError, match='sample_id "S1" appears more than once'):
        read_sample_ids(samplesheet)


def test_read_sample_ids_rejects_missing_fastq_before_context_hash(
    tmp_path: Path,
) -> None:
    samplesheet = write_samplesheet(
        tmp_path / "samples.csv",
        [
            {
                "sample_id": "S1",
                "platform": "illumina",
                "fastq1": str(tmp_path / "missing_R1.fastq.gz"),
            },
        ],
    )

    with pytest.raises(ResolutionError, match="This file was not found"):
        read_sample_ids(samplesheet)


def test_context_hash_is_derived_from_resolved_sample_ids(tmp_path: Path) -> None:
    samplesheet = write_samplesheet(
        tmp_path / "samples.csv",
        [
            {"sample_id": "S2", "platform": "ont", "srr": "SRR2"},
            {"sample_id": "S1", "platform": "ont", "srr": "SRR1"},
        ],
    )

    assert compute_sample_set_id(read_sample_ids(samplesheet)) == compute_sample_set_id(
        ["S1", "S2"],
    )
