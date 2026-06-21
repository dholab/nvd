"""Tests for samplesheet generation helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from py_nvd.cli.commands.samplesheet import (
    SamplesheetGenerationError,
    _extract_stem,
    scan_fastq_directory,
)

if TYPE_CHECKING:
    from pathlib import Path


def touch(path: Path) -> None:
    """Create an empty FASTQ fixture."""
    path.write_text("", encoding="utf-8")


def test_extract_stem_keeps_casava_suffix_without_sanitize() -> None:
    """Default sample IDs preserve existing CASAVA sample/lane suffix behavior."""
    assert _extract_stem("patient-001_S7_L003_R1_001.fastq.gz") == "patient-001_S7_L003"


def test_extract_stem_strips_casava_suffix_with_sanitize() -> None:
    """Sanitized sample IDs strip Illumina sample and lane suffixes."""
    assert (
        _extract_stem("patient-001_S7_L003_R1_001.fastq.gz", sanitize=True)
        == "patient-001"
    )


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("tumor_R1_control_S7_L003_R1_001.fastq.gz", "tumor_R1_control"),
        ("sample_R1_S7_L003_R2_001.fastq.gz", "sample_R1"),
        ("sample_S7like_L003_R1_001.fastq.gz", "sample_S7like_L003"),
        ("sample_S7_L003_extra_R1_001.fastq.gz", "sample_S7_L003_extra"),
    ],
)
def test_extract_stem_sanitize_only_strips_terminal_casava_suffix(
    filename: str,
    expected: str,
) -> None:
    """Sanitization preserves read-like tokens inside the sample ID."""
    assert _extract_stem(filename, sanitize=True) == expected


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("tumor_R1_control_S7_L003_R1_001.fastq.gz", "tumor_R1_control_S7_L003"),
        ("sample_R1_S7_L003_R2_001.fastq.gz", "sample_R1_S7_L003"),
        ("sample_R1.fastq.gz", "sample"),
        ("sample_R1_001.fastq.gz", "sample"),
        ("sample.1.fq.gz", "sample"),
    ],
)
def test_extract_stem_documents_read_indicator_ambiguity(
    filename: str,
    expected: str,
) -> None:
    """Only final read-indicator tokens are treated as read structure."""
    assert _extract_stem(filename) == expected


def test_scan_fastq_directory_sanitizes_sample_ids_after_pairing(
    tmp_path: Path,
) -> None:
    """Sanitization changes the sample ID while preserving paired FASTQ paths."""
    read1 = tmp_path / "patient-001_S7_L003_R1_001.fastq.gz"
    read2 = tmp_path / "patient-001_S7_L003_R2_001.fastq.gz"
    touch(read1)
    touch(read2)

    samples, warnings = scan_fastq_directory(tmp_path, sanitize=True)

    assert warnings == []
    assert len(samples) == 1
    assert samples[0].sample_id == "patient-001"
    assert samples[0].fastq1 == str(read1.resolve())
    assert samples[0].fastq2 == str(read2.resolve())


def test_scan_fastq_directory_does_not_sanitize_by_default(tmp_path: Path) -> None:
    """Existing samplesheet generation behavior is preserved by default."""
    read1 = tmp_path / "patient-001_S7_L003_R1_001.fastq.gz"
    read2 = tmp_path / "patient-001_S7_L003_R2_001.fastq.gz"
    touch(read1)
    touch(read2)

    samples, warnings = scan_fastq_directory(tmp_path)

    assert warnings == []
    assert len(samples) == 1
    assert samples[0].sample_id == "patient-001_S7_L003"


def test_scan_fastq_directory_sanitizes_ids_without_collapsing_lanes(
    tmp_path: Path,
) -> None:
    """Sanitization does not use the sanitized ID as the pairing key."""
    for lane in ("L001", "L002"):
        touch(tmp_path / f"patient-001_S7_{lane}_R1_001.fastq.gz")
        touch(tmp_path / f"patient-001_S7_{lane}_R2_001.fastq.gz")

    samples, warnings = scan_fastq_directory(tmp_path, sanitize=True)

    assert warnings == []
    assert [sample.sample_id for sample in samples] == ["patient-001", "patient-001"]
    assert len(samples) == 2


def test_scan_fastq_directory_groups_casava_lanes_into_glob_columns(
    tmp_path: Path,
) -> None:
    """Lane grouping emits glob columns, not glob patterns in exact path columns."""
    for lane in ("L001", "L002"):
        touch(tmp_path / f"patient-001_S7_{lane}_R1_001.fastq.gz")
        touch(tmp_path / f"patient-001_S7_{lane}_R2_001.fastq.gz")

    samples, warnings = scan_fastq_directory(
        tmp_path,
        sanitize=True,
        group_lanes=True,
    )

    assert warnings == []
    assert len(samples) == 1
    assert samples[0].sample_id == "patient-001"
    assert samples[0].fastq1 == ""
    assert samples[0].fastq2 == ""
    assert samples[0].fastq1_glob == str(
        tmp_path.resolve()
        / "patient-001_S7_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz",
    )
    assert samples[0].fastq2_glob == str(
        tmp_path.resolve()
        / "patient-001_S7_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9].fastq.gz",
    )


def test_scan_fastq_directory_group_lanes_preserves_read_like_sample_tokens(
    tmp_path: Path,
) -> None:
    """CASAVA parsing does not mistake read-like sample tokens for read markers."""
    for lane in ("L001", "L002"):
        touch(tmp_path / f"tumor_R1_control_S7_{lane}_R1_001.fastq.gz")
        touch(tmp_path / f"tumor_R1_control_S7_{lane}_R2_001.fastq.gz")

    samples, warnings = scan_fastq_directory(
        tmp_path,
        sanitize=True,
        group_lanes=True,
    )

    assert warnings == []
    assert [sample.sample_id for sample in samples] == ["tumor_R1_control"]


def test_scan_fastq_directory_group_lanes_rejects_missing_lane_mate(
    tmp_path: Path,
) -> None:
    """Generation fails before an expensive run when lane pairs are incomplete."""
    touch(tmp_path / "sample_S1_L001_R1_001.fastq.gz")
    touch(tmp_path / "sample_S1_L002_R1_001.fastq.gz")
    touch(tmp_path / "sample_S1_L002_R2_001.fastq.gz")

    with pytest.raises(SamplesheetGenerationError, match="missing matching R2"):
        scan_fastq_directory(tmp_path, group_lanes=True)


def test_scan_fastq_directory_group_lanes_rejects_non_casava_names(
    tmp_path: Path,
) -> None:
    """Lane grouping is CASAVA-only rather than guessing arbitrary names."""
    touch(tmp_path / "sample.R1.fastq.gz")
    touch(tmp_path / "sample.R2.fastq.gz")

    with pytest.raises(SamplesheetGenerationError, match="CASAVA-style"):
        scan_fastq_directory(tmp_path, group_lanes=True)
