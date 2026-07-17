"""Tests for samplesheet generation helpers."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

import pytest
from typer.testing import CliRunner

from py_nvd.cli.app import app
from py_nvd.cli.commands.samplesheet import (
    FastqGrouping,
    SamplesheetGenerationError,
    _extract_stem,
    scan_fastq_directory,
)
from py_nvd.samplesheet_validation import validate_samplesheet

if TYPE_CHECKING:
    from pathlib import Path

runner = CliRunner()


def touch(path: Path) -> Path:
    """Create an empty FASTQ fixture."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")
    return path


def write_samplesheet(path: Path, rows: list[dict[str, str]]) -> Path:
    """Write a samplesheet fixture with optional glob columns."""
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
    assert samples[0].fastq1 == str(read1.absolute())
    assert samples[0].fastq2 == str(read2.absolute())


def test_scan_fastq_directory_preserves_symlink_directory_namespace(
    tmp_path: Path,
) -> None:
    """Generated paths keep the user-provided stable symlink namespace."""
    real_dir = tmp_path / "real-fastqs"
    read1 = touch(real_dir / "patient-001_S7_L003_R1_001.fastq.gz")
    read2 = touch(real_dir / "patient-001_S7_L003_R2_001.fastq.gz")
    symlink_dir = tmp_path / "stable-fastqs"
    symlink_dir.symlink_to(real_dir, target_is_directory=True)

    samples, warnings = scan_fastq_directory(symlink_dir, sanitize=True)

    assert warnings == []
    assert len(samples) == 1
    assert samples[0].fastq1 == str(symlink_dir.absolute() / read1.name)
    assert samples[0].fastq2 == str(symlink_dir.absolute() / read2.name)


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


def test_samplesheet_generate_rejects_duplicate_sanitized_ids_before_write(
    tmp_path: Path,
) -> None:
    """CLI generation fails safely when sanitization would duplicate IDs."""
    fastq_dir = tmp_path / "fastqs"
    for lane in ("L001", "L002"):
        touch(fastq_dir / f"patient-001_S7_{lane}_R1_001.fastq.gz")
        touch(fastq_dir / f"patient-001_S7_{lane}_R2_001.fastq.gz")
    output = tmp_path / "samplesheet.csv"

    result = runner.invoke(
        app,
        [
            "samplesheet",
            "generate",
            "--from-dir",
            str(fastq_dir),
            "--platform",
            "illumina",
            "--sanitize",
            "--output",
            str(output),
            "--force",
        ],
    )

    assert result.exit_code != 0
    assert not output.exists()
    assert "--group-by illumina-lanes" in result.output
    assert "omit" in result.output
    assert "--sanitize" in result.output


def test_scan_fastq_directory_groups_casava_lanes_into_glob_columns(
    tmp_path: Path,
) -> None:
    """Illumina lane grouping emits glob columns, not exact path columns."""
    for lane in ("L001", "L002"):
        touch(tmp_path / f"patient-001_S7_{lane}_R1_001.fastq.gz")
        touch(tmp_path / f"patient-001_S7_{lane}_R2_001.fastq.gz")

    samples, warnings = scan_fastq_directory(
        tmp_path,
        sanitize=True,
        group_by=FastqGrouping.ILLUMINA_LANES,
    )

    assert warnings == []
    assert len(samples) == 1
    assert samples[0].sample_id == "patient-001"
    assert samples[0].fastq1 == ""
    assert samples[0].fastq2 == ""
    assert samples[0].fastq1_glob == str(
        tmp_path.absolute()
        / "patient-001_S7_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq.gz",
    )
    assert samples[0].fastq2_glob == str(
        tmp_path.absolute()
        / "patient-001_S7_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9].fastq.gz",
    )


def test_scan_fastq_directory_group_by_illumina_lanes_preserves_read_like_sample_tokens(
    tmp_path: Path,
) -> None:
    """CASAVA parsing does not mistake read-like sample tokens for read markers."""
    for lane in ("L001", "L002"):
        touch(tmp_path / f"tumor_R1_control_S7_{lane}_R1_001.fastq.gz")
        touch(tmp_path / f"tumor_R1_control_S7_{lane}_R2_001.fastq.gz")

    samples, warnings = scan_fastq_directory(
        tmp_path,
        sanitize=True,
        group_by=FastqGrouping.ILLUMINA_LANES,
    )

    assert warnings == []
    assert [sample.sample_id for sample in samples] == ["tumor_R1_control"]


def test_scan_fastq_directory_group_by_illumina_lanes_rejects_missing_lane_mate(
    tmp_path: Path,
) -> None:
    """Generation fails before an expensive run when lane pairs are incomplete."""
    touch(tmp_path / "sample_S1_L001_R1_001.fastq.gz")
    touch(tmp_path / "sample_S1_L002_R1_001.fastq.gz")
    touch(tmp_path / "sample_S1_L002_R2_001.fastq.gz")

    with pytest.raises(SamplesheetGenerationError, match="missing matching R2"):
        scan_fastq_directory(tmp_path, group_by=FastqGrouping.ILLUMINA_LANES)


def test_scan_fastq_directory_group_by_illumina_lanes_rejects_non_casava_names(
    tmp_path: Path,
) -> None:
    """Lane grouping is CASAVA-only rather than guessing arbitrary names."""
    touch(tmp_path / "sample.R1.fastq.gz")
    touch(tmp_path / "sample.R2.fastq.gz")

    with pytest.raises(SamplesheetGenerationError, match="CASAVA-style"):
        scan_fastq_directory(tmp_path, group_by=FastqGrouping.ILLUMINA_LANES)


def test_scan_fastq_directory_groups_ont_barcodes_into_single_glob_columns(
    tmp_path: Path,
) -> None:
    touch(tmp_path / "barcode02" / "chunk2.fastq.gz")
    touch(tmp_path / "barcode02" / "chunk10.fastq.gz")
    touch(tmp_path / "barcode01" / "chunk1.fastq.gz")

    samples, warnings = scan_fastq_directory(
        tmp_path,
        group_by=FastqGrouping.ONT_BARCODES,
    )

    assert warnings == []
    assert [sample.sample_id for sample in samples] == ["barcode01", "barcode02"]
    assert samples[0].fastq1_glob == str(
        tmp_path.absolute() / "barcode01" / "*.fastq.gz",
    )
    assert samples[0].fastq2_glob == ""
    assert samples[1].fastq1_glob == str(
        tmp_path.absolute() / "barcode02" / "*.fastq.gz",
    )
    assert samples[1].fastq2_glob == ""


def test_scan_fastq_directory_group_by_ont_barcodes_skips_non_barcode_dirs(
    tmp_path: Path,
) -> None:
    touch(tmp_path / "barcode01" / "chunk1.fastq.gz")
    touch(tmp_path / "unclassified" / "chunk1.fastq.gz")

    samples, warnings = scan_fastq_directory(
        tmp_path,
        group_by=FastqGrouping.ONT_BARCODES,
    )

    assert [sample.sample_id for sample in samples] == ["barcode01"]
    assert warnings == ["Skipping non-barcode directory: unclassified"]


def test_scan_fastq_directory_group_by_ont_barcodes_rejects_mixed_extensions(
    tmp_path: Path,
) -> None:
    touch(tmp_path / "barcode01" / "chunk1.fastq.gz")
    touch(tmp_path / "barcode01" / "chunk2.fastq")

    with pytest.raises(SamplesheetGenerationError, match="mixed FASTQ extensions"):
        scan_fastq_directory(tmp_path, group_by=FastqGrouping.ONT_BARCODES)


def test_samplesheet_generate_requires_group_by_value(tmp_path: Path) -> None:
    fastq_dir = tmp_path / "fastqs"
    touch(fastq_dir / "reads.fastq.gz")

    result = runner.invoke(
        app,
        [
            "samplesheet",
            "generate",
            "--from-dir",
            str(fastq_dir),
            "--platform",
            "ont",
            "--group-by",
        ],
    )

    assert result.exit_code != 0


def test_samplesheet_generate_rejects_cross_platform_grouping(tmp_path: Path) -> None:
    fastq_dir = tmp_path / "fastqs"
    touch(fastq_dir / "reads.fastq.gz")

    result = runner.invoke(
        app,
        [
            "samplesheet",
            "generate",
            "--from-dir",
            str(fastq_dir),
            "--platform",
            "ont",
            "--group-by",
            "illumina-lanes",
        ],
    )

    assert result.exit_code != 0
    assert "--group-by illumina-lanes requires --platform illumina" in result.output


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (
            {"sample_id": "S1", "platform": "illumina", "fastq1": "reads.fastq.gz"},
            "FASTQ paths and glob patterns must be absolute",
        ),
        (
            {"sample_id": "S1", "platform": "illumina", "fastq1": "*.fastq.gz"},
            "fastq1 column contains glob characters",
        ),
        (
            {"sample_id": "S1", "platform": "illumina", "fastq2": "/reads/r2.fastq.gz"},
            "fastq2 was provided without fastq1",
        ),
        (
            {
                "sample_id": "S1",
                "platform": "illumina",
                "fastq2_glob": "*_R2_*.fastq.gz",
            },
            "fastq2_glob was provided without fastq1_glob",
        ),
        (
            {"sample_id": "S1", "platform": "illumina"},
            "must provide exact FASTQ path",
        ),
        (
            {"sample_id": "S1", "platform": "pacbio", "srr": "SRR123"},
            "Invalid platform",
        ),
        (
            {"sample_id": "S1", "platform": "illumina", "fastq1_glob": "*.fastq.gz"},
            "FASTQ paths and glob patterns must be absolute",
        ),
    ],
)
def test_validate_samplesheet_reuses_read_resolution_errors(
    tmp_path: Path,
    row: dict[str, str],
    expected: str,
) -> None:
    """Validation catches the same read declaration errors as a run preflight."""
    result = validate_samplesheet(write_samplesheet(tmp_path / "samples.csv", [row]))

    assert not result.valid
    assert any(expected in error for error in result.errors)


def test_validate_samplesheet_rejects_missing_exact_file(tmp_path: Path) -> None:
    """Exact FASTQ paths must resolve during validation."""
    missing = tmp_path / "missing_R1.fastq.gz"

    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1": str(missing),
                },
            ],
        ),
    )

    assert not result.valid
    assert any("was not found" in error for error in result.errors)


def test_validate_samplesheet_rejects_unmatched_glob(tmp_path: Path) -> None:
    """Glob patterns must match at least one visible FASTQ."""
    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1_glob": str(tmp_path / "*_R1_*.fastq.gz"),
                },
            ],
        ),
    )

    assert not result.valid
    assert any("did not match any files" in error for error in result.errors)


def test_validate_samplesheet_rejects_mismatched_paired_glob_lanes(
    tmp_path: Path,
) -> None:
    """Paired globs must describe the same CASAVA lane set."""
    lane_dir = tmp_path / "lanes"
    touch(lane_dir / "Sample_S1_L001_R1_001.fastq.gz")
    touch(lane_dir / "Sample_S1_L002_R1_001.fastq.gz")
    touch(lane_dir / "Sample_S1_L001_R2_001.fastq.gz")

    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1_glob": str(lane_dir / "*_R1_*.fastq.gz"),
                    "fastq2_glob": str(lane_dir / "*_R2_*.fastq.gz"),
                },
            ],
        ),
    )

    assert not result.valid
    assert any("do not describe the same lane set" in error for error in result.errors)


def test_validate_samplesheet_accepts_terminal_marker_paired_globs(
    tmp_path: Path,
) -> None:
    """Paired glob validation accepts safe terminal read markers."""
    lane_dir = tmp_path / "lanes"
    touch(lane_dir / "Sample.R1.fastq.gz")
    touch(lane_dir / "Sample.R2.fastq.gz")

    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1_glob": str(lane_dir / "*.R1.fastq.gz"),
                    "fastq2_glob": str(lane_dir / "*.R2.fastq.gz"),
                },
            ],
        ),
    )

    assert result.valid
    assert result.samples == ["S1"]


def test_validate_samplesheet_rejects_nonterminal_paired_globs(tmp_path: Path) -> None:
    """Paired glob validation does not pair arbitrary sorted glob results."""
    lane_dir = tmp_path / "lanes"
    touch(lane_dir / "Sample_R1_extra.fastq.gz")
    touch(lane_dir / "Sample_R2_extra.fastq.gz")

    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1_glob": str(lane_dir / "*_R1_*.fastq.gz"),
                    "fastq2_glob": str(lane_dir / "*_R2_*.fastq.gz"),
                },
            ],
        ),
    )

    assert not result.valid
    assert any("simple terminal read markers" in error for error in result.errors)


def test_validate_samplesheet_rejects_mixed_glob_compressions(tmp_path: Path) -> None:
    """Glob-expanded reads for one sample must use one compression kind."""
    reads_dir = tmp_path / "reads"
    touch(reads_dir / "chunk1.fastq.gz")
    touch(reads_dir / "chunk2.fastq.xz")

    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "ont",
                    "fastq1_glob": str(reads_dir / "*.fastq*"),
                },
            ],
        ),
    )

    assert not result.valid
    assert any("must use one compression type" in error for error in result.errors)


def test_validate_samplesheet_rejects_duplicate_sample_ids(tmp_path: Path) -> None:
    """NVD expects one row per emitted sample ID."""
    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {"sample_id": "S1", "platform": "illumina", "srr": "SRR1"},
                {"sample_id": "S1", "platform": "illumina", "srr": "SRR2"},
            ],
        ),
    )

    assert not result.valid
    assert any(
        'sample_id "S1" appears more than once' in error for error in result.errors
    )


@pytest.mark.parametrize(
    "sample_id",
    ["patient 1", "../patient-1", "patient'1", "patient;1", "patient@1"],
)
def test_validate_samplesheet_rejects_unsafe_sample_ids(
    tmp_path: Path,
    sample_id: str,
) -> None:
    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [{"sample_id": sample_id, "platform": "illumina", "srr": "SRR1"}],
        ),
    )

    assert not result.valid
    assert any(
        "may contain only letters, numbers, dots, underscores, and hyphens" in error
        for error in result.errors
    )


def test_validate_samplesheet_accepts_safe_sample_id_punctuation(
    tmp_path: Path,
) -> None:
    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "patient-1_v2.3",
                    "platform": "illumina",
                    "srr": "SRR1",
                },
            ],
        ),
    )

    assert result.valid


def test_validate_samplesheet_reports_shadowed_sources_as_warnings(
    tmp_path: Path,
) -> None:
    """Validation surfaces resolver warnings without failing the samplesheet."""
    r1 = touch(tmp_path / "Sample_S1_L001_R1_001.fastq.gz")
    r2 = touch(tmp_path / "Sample_S1_L001_R2_001.fastq.gz")

    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1": str(r1),
                    "fastq2": str(r2),
                    "srr": "SRR123",
                },
            ],
        ),
    )

    assert result.valid
    assert result.warnings == ['sample "S1": using fastq1/fastq2; ignoring srr']


def test_validate_samplesheet_warns_on_suspicious_manual_sra_accession(
    tmp_path: Path,
) -> None:
    """Manual SRA rows use the same warning-only accession shape check as generation."""
    result = validate_samplesheet(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "srr": "not-an-accession",
                },
            ],
        ),
    )

    assert result.valid
    assert result.warnings == [
        "sample \"S1\": SRA accession 'not-an-accession' does not look like an SRR/ERR/DRR accession",
    ]
