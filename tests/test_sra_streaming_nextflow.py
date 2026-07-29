from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    ("mode_args", "expected_messages"),
    [
        (
            ["--experimental", "false", "--skip_fastqc", "true"],
            [
                "SRA Toolkit streaming is available only when experimental features are enabled",
                "experimental = true",
                "stream_sra = false",
            ],
        ),
        (
            ["--experimental", "true", "--skip_fastqc", "false"],
            [
                "SRA Toolkit streaming cannot run with raw FastQC enabled",
                "skip_fastqc = true",
                "stream_sra = false",
            ],
        ),
    ],
)
def test_sra_streaming_rejects_incompatible_modes_with_corrective_message(
    tmp_path: Path,
    mode_args: list[str],
    expected_messages: list[str],
) -> None:
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text(
        "sample_id,srr,platform,fastq1,fastq2\n",
        encoding="utf-8",
    )

    completed = subprocess.run(  # noqa: S603
        [  # noqa: S607
            "nextflow",
            "run",
            ".",
            "-profile",
            "test",
            "--samplesheet",
            str(samplesheet),
            "--stream_sra",
            "true",
            *mode_args,
            "--skip_blast",
            "true",
            "--skip_assembly",
            "true",
        ],
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=True,
        timeout=30,
    )

    diagnostics = completed.stdout + completed.stderr
    assert completed.returncode != 0
    for message in expected_messages:
        assert message in diagnostics
