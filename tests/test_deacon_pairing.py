"""Focused behavioral tests for pair-atomic Deacon read depletion."""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DEACON_MODULE = ROOT / "modules" / "deacon"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None or shutil.which("deacon") is None,
    reason="DEACON_DEPLETE behavior requires Nextflow and locked Deacon",
)

CONTAMINANT = (
    "ACGTTGCAGATCGTAGCTAGGCTAACGATCGTACGGCATTCGATGCTAGCATCGTACGTATCGATGGCATACGT"
)
NONCONTAMINANT = (
    "TTAGCCGTATGACCTGAGTCCATGGTACCACTGATCGTGCATAGGTCACCTAGTGCAGATCCTGACGTAGCA"
)
MATE_ONE = 1
MATE_TWO = 2


def write_interleaved_fastq(
    path: Path,
    *,
    contaminant_mate: int,
) -> None:
    """Write one SRA-style interleaved pair with one contaminant mate."""
    sequences = {
        MATE_ONE: CONTAMINANT if contaminant_mate == MATE_ONE else NONCONTAMINANT,
        MATE_TWO: CONTAMINANT if contaminant_mate == MATE_TWO else NONCONTAMINANT,
    }
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for mate in (1, 2):
            sequence = sequences[mate]
            handle.write(f"@spot.{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def build_deacon_index(tmp_path: Path) -> Path:
    """Build a small exact-match Deacon index for the contaminant fixture."""
    reference = tmp_path / "contaminant.fasta"
    reference.write_text(f">contaminant\n{CONTAMINANT}\n", encoding="utf-8")
    index = tmp_path / "contaminant.idx"
    completed = subprocess.run(  # noqa: S603
        [  # noqa: S607
            "deacon",
            "index",
            "build",
            "-k",
            "31",
            "-w",
            "1",
            "-o",
            str(index),
            str(reference),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=45,
    )
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    return index


def run_depletion(
    tmp_path: Path,
    *,
    read_structure: str,
    contaminant_mate: int = 1,
) -> str:
    """Run the public DEACON_DEPLETE process and return its FASTQ text."""
    reads = tmp_path / "reads.fastq.gz"
    if read_structure == "interleaved":
        write_interleaved_fastq(reads, contaminant_mate=contaminant_mate)
    else:
        with gzip.open(reads, "wt", encoding="utf-8") as handle:
            handle.write(
                f"@single\n{NONCONTAMINANT}\n+\n{'I' * len(NONCONTAMINANT)}\n",
            )
    index = build_deacon_index(tmp_path)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.host_abs_threshold = 1
params.host_rel_threshold = 0.0

include {{ DEACON_DEPLETE }} from '{DEACON_MODULE}'

workflow {{
    DEACON_DEPLETE(Channel.of(tuple(
        'sample_A',
        'illumina',
        '{read_structure}',
        'single_read',
        file('{reads}'),
        file('{index}'),
    )))
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "-C", "/dev/null", "run", str(workflow)],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=45,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    outputs = list(
        (tmp_path / "work").glob("**/sample_A.single_read.depleted.fastq.gz")
    )
    assert len(outputs) == 1, outputs
    with gzip.open(outputs[0], "rt", encoding="utf-8") as handle:
        return handle.read()


@pytest.mark.parametrize("contaminant_mate", [MATE_ONE, MATE_TWO])
def test_depletion_removes_complete_interleaved_pair_when_either_mate_matches(
    tmp_path: Path,
    contaminant_mate: int,
) -> None:
    """A contaminant match on either mate removes the complete pair."""
    output = run_depletion(
        tmp_path,
        read_structure="interleaved",
        contaminant_mate=contaminant_mate,
    )

    assert output == ""


def test_depletion_keeps_direct_file_filtering_for_single_reads(tmp_path: Path) -> None:
    """Single reads retain Deacon's native direct gzip input path."""
    output = run_depletion(tmp_path, read_structure="single")

    assert "@single\n" in output
