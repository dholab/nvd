"""Tests for samplesheet read-input resolution."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING

import pytest
from resolve_read_inputs import (
    ResolutionError,
    main,
    read_rows,
    resolve_rows,
    write_jsonl,
)

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


def read_jsonl(path: Path) -> list[dict[str, object]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def test_exact_paired_files_take_precedence_with_warning(tmp_path: Path) -> None:
    r1 = touch(tmp_path / "exact" / "sample_R1.fastq.gz")
    r2 = touch(tmp_path / "exact" / "sample_R2.fastq.gz")

    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "srr": "SRR123",
                    "platform": "illumina",
                    "fastq1": str(r1),
                    "fastq2": str(r2),
                    "fastq1_glob": str(tmp_path / "lanes" / "*_R1_*.fastq.gz"),
                    "fastq2_glob": str(tmp_path / "lanes" / "*_R2_*.fastq.gz"),
                },
            ],
        ),
    )

    resolved = resolve_rows(rows, cwd=tmp_path)

    assert resolved.samples[0].source == "paired_files"
    assert resolved.samples[0].r1 == (r1,)
    assert resolved.samples[0].r2 == (r2,)
    assert resolved.samples[0].warnings == (
        "using fastq1/fastq2; ignoring fastq1_glob, fastq2_glob, and srr",
    )
    assert resolved.warnings == (
        'sample "S1": using fastq1/fastq2; ignoring fastq1_glob, fastq2_glob, and srr',
    )


def test_single_file_accepts_nanopore_alias(tmp_path: Path) -> None:
    reads = touch(tmp_path / "reads.fastq.xz")
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "ONT1",
                    "platform": "nanopore",
                    "fastq1": str(reads),
                },
            ],
        ),
    )

    resolved = resolve_rows(rows, cwd=tmp_path)

    assert resolved.samples[0].source == "single_file"
    assert resolved.samples[0].platform == "ont"
    assert resolved.samples[0].reads == (reads,)
    assert resolved.samples[0].warnings == (
        "platform alias 'nanopore' was normalized to 'ont'",
    )


def test_exact_paths_are_absolute_but_symlink_preserving(tmp_path: Path) -> None:
    target = touch(tmp_path / "archive" / "sample.fastq.gz")
    link = tmp_path / "stable" / "sample.fastq.gz"
    link.parent.mkdir(parents=True)
    link.symlink_to(target)
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "ONT1",
                    "platform": "ont",
                    "fastq1": str(link),
                },
            ],
        ),
    )

    resolved = resolve_rows(rows, cwd=tmp_path)

    assert resolved.samples[0].reads == (link,)
    assert resolved.samples[0].reads[0] != target


def test_paired_globs_pair_by_casava_key_not_path_order(tmp_path: Path) -> None:
    lane_dir = tmp_path / "lanes"
    r1_l002 = touch(lane_dir / "Sample_S1_L002_R1_001.fastq.gz")
    r2_l002 = touch(lane_dir / "Sample_S1_L002_R2_001.fastq.gz")
    r1_l001 = touch(lane_dir / "Sample_S1_L001_R1_001.fastq.gz")
    r2_l001 = touch(lane_dir / "Sample_S1_L001_R2_001.fastq.gz")
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1_glob": str(lane_dir / "*_R1_*.fastq.gz"),
                    "fastq2_glob": str(lane_dir / "*_R2_*.fastq.gz"),
                    "srr": "SRR123",
                },
            ],
        ),
    )

    resolved = resolve_rows(rows, cwd=tmp_path)

    assert resolved.samples[0].source == "paired_globs"
    assert resolved.samples[0].r1 == (r1_l001, r1_l002)
    assert resolved.samples[0].r2 == (r2_l001, r2_l002)
    assert resolved.samples[0].warnings == ("using fastq1_glob/fastq2_glob; ignoring srr",)


def test_single_glob_sorts_naturally_without_casava_requirement(tmp_path: Path) -> None:
    reads10 = touch(tmp_path / "ont" / "chunk10.fastq.gz")
    reads2 = touch(tmp_path / "ont" / "chunk2.fastq.gz")
    reads1 = touch(tmp_path / "ont" / "chunk1.fastq.gz")
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "ONT1",
                    "platform": "ont",
                    "fastq1_glob": str(tmp_path / "ont" / "*.fastq.gz"),
                },
            ],
        ),
    )

    resolved = resolve_rows(rows, cwd=tmp_path)

    assert resolved.samples[0].source == "single_glob"
    assert resolved.samples[0].reads == (reads1, reads2, reads10)


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        (
            {"sample_id": "S1", "platform": "illumina", "fastq1": "*.fastq.gz"},
            "fastq1 column contains glob characters",
        ),
        (
            {"sample_id": "S1", "platform": "pacbio", "srr": "SRR123"},
            'Invalid platform for sample "S1"',
        ),
        (
            {"sample_id": "S1", "platform": "illumina", "fastq2_glob": "*_R2_*.fastq.gz"},
            "fastq2_glob was provided without fastq1_glob",
        ),
    ],
)
def test_friendly_validation_errors(
    tmp_path: Path,
    row: dict[str, str],
    expected: str,
) -> None:
    rows = read_rows(write_samplesheet(tmp_path / "samples.csv", [row]))

    with pytest.raises(ResolutionError, match=expected):
        resolve_rows(rows, cwd=tmp_path)


def test_duplicate_sample_id_errors(tmp_path: Path) -> None:
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {"sample_id": "S1", "platform": "illumina", "srr": "SRR1"},
                {"sample_id": "S1", "platform": "illumina", "srr": "SRR2"},
            ],
        ),
    )

    with pytest.raises(ResolutionError, match='sample_id "S1" appears more than once'):
        resolve_rows(rows, cwd=tmp_path)


def test_paired_glob_rejects_missing_pair(tmp_path: Path) -> None:
    lane_dir = tmp_path / "lanes"
    touch(lane_dir / "Sample_S1_L001_R1_001.fastq.gz")
    touch(lane_dir / "Sample_S1_L002_R1_001.fastq.gz")
    touch(lane_dir / "Sample_S1_L001_R2_001.fastq.gz")
    rows = read_rows(
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

    with pytest.raises(ResolutionError, match="do not describe the same lane set"):
        resolve_rows(rows, cwd=tmp_path)


def test_paired_glob_rejects_non_casava_names(tmp_path: Path) -> None:
    lane_dir = tmp_path / "lanes"
    touch(lane_dir / "Sample.R1.fastq.gz")
    touch(lane_dir / "Sample.R2.fastq.gz")
    rows = read_rows(
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

    with pytest.raises(ResolutionError, match="CASAVA-style Illumina lane name"):
        resolve_rows(rows, cwd=tmp_path)


def test_glob_rejects_mixed_compression_within_sample(tmp_path: Path) -> None:
    reads_dir = tmp_path / "ont"
    touch(reads_dir / "chunk1.fastq.gz")
    touch(reads_dir / "chunk2.fastq.xz")
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "ONT1",
                    "platform": "ont",
                    "fastq1_glob": str(reads_dir / "*.fastq*"),
                },
            ],
        ),
    )

    with pytest.raises(ResolutionError, match="must use one compression type"):
        resolve_rows(rows, cwd=tmp_path)


def test_zst_suffix_is_supported_for_deacon_extension_based_decoding(tmp_path: Path) -> None:
    reads = touch(tmp_path / "reads.fastq.zst")
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "ONT1",
                    "platform": "ont",
                    "fastq1": str(reads),
                },
            ],
        ),
    )

    resolved = resolve_rows(rows, cwd=tmp_path)

    assert resolved.samples[0].source == "single_file"
    assert resolved.samples[0].reads == (reads,)


def test_relative_exact_paths_are_rejected(tmp_path: Path) -> None:
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "ONT1",
                    "platform": "ont",
                    "fastq1": "reads.fastq.gz",
                },
            ],
        ),
    )

    with pytest.raises(ResolutionError, match="must be absolute"):
        resolve_rows(rows, cwd=tmp_path)


def test_relative_glob_patterns_are_rejected(tmp_path: Path) -> None:
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "ONT1",
                    "platform": "ont",
                    "fastq1_glob": "*.fastq.gz",
                },
            ],
        ),
    )

    with pytest.raises(ResolutionError, match="must be absolute"):
        resolve_rows(rows, cwd=tmp_path)


def test_write_jsonl_uses_expected_shape(tmp_path: Path) -> None:
    r1 = touch(tmp_path / "Sample_S1_L001_R1_001.fastq.gz")
    r2 = touch(tmp_path / "Sample_S1_L001_R2_001.fastq.gz")
    rows = read_rows(
        write_samplesheet(
            tmp_path / "samples.csv",
            [
                {
                    "sample_id": "S1",
                    "platform": "illumina",
                    "fastq1": str(r1),
                    "fastq2": str(r2),
                },
            ],
        ),
    )
    resolved = resolve_rows(rows, cwd=tmp_path)
    output = tmp_path / "resolved.jsonl"

    write_jsonl(resolved.samples, output)

    assert read_jsonl(output) == [
        {
            "sample_id": "S1",
            "platform": "illumina",
            "source": "paired_files",
            "r1": [str(r1)],
            "r2": [str(r2)],
            "warnings": [],
        },
    ]


def test_cli_writes_jsonl_and_warnings_to_stderr(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    r1 = touch(tmp_path / "Sample_S1_L001_R1_001.fastq.gz")
    r2 = touch(tmp_path / "Sample_S1_L001_R2_001.fastq.gz")
    samplesheet = write_samplesheet(
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
    )
    output = tmp_path / "resolved.jsonl"
    monkeypatch.setattr(
        "sys.argv",
        [
            "resolve_read_inputs.py",
            "--samplesheet",
            str(samplesheet),
            "--output-jsonl",
            str(output),
        ],
    )

    main()

    captured = capsys.readouterr()
    assert 'Warning: sample "S1": using fastq1/fastq2; ignoring srr' in captured.err
    assert read_jsonl(output)[0]["warnings"] == ["using fastq1/fastq2; ignoring srr"]
