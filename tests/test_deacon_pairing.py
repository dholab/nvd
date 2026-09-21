"""Focused behavioral tests for pair-atomic Deacon read depletion."""

from __future__ import annotations

import gzip
import json
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
    canonical_names: bool = False,
) -> None:
    """Write one SRA-style interleaved pair with one contaminant mate."""
    sequences = {
        MATE_ONE: CONTAMINANT if contaminant_mate == MATE_ONE else NONCONTAMINANT,
        MATE_TWO: CONTAMINANT if contaminant_mate == MATE_TWO else NONCONTAMINANT,
    }
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for mate in (1, 2):
            sequence = sequences[mate]
            separator = "/" if canonical_names else "."
            handle.write(
                f"@spot{separator}{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n",
            )


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
    check_pairs: bool = False,
    canonical_names: bool = False,
) -> tuple[str, str]:
    """Run DEACON_DEPLETE and return its FASTQ text and rendered task command."""
    reads = tmp_path / "reads.fastq.gz"
    if read_structure == "interleaved":
        write_interleaved_fastq(
            reads,
            contaminant_mate=contaminant_mate,
            canonical_names=canonical_names,
        )
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
params.check_pairs = {str(check_pairs).lower()}

include {{ DEACON_DEPLETE }} from '{DEACON_MODULE}'

workflow {{
    meta = [id: 'sample_A', platform: 'illumina', read_structure: '{read_structure}', query_class: 'single_read']
    DEACON_DEPLETE(Channel.of(tuple(
        meta,
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
        output = handle.read()
    commands = list((tmp_path / "work").glob("**/.command.sh"))
    assert len(commands) == 1, commands
    return output, commands[0].read_text(encoding="utf-8")


@pytest.mark.parametrize("contaminant_mate", [MATE_ONE, MATE_TWO])
def test_depletion_removes_complete_interleaved_pair_when_either_mate_matches(
    tmp_path: Path,
    contaminant_mate: int,
) -> None:
    """A contaminant match on either mate removes the complete pair."""
    output, _command = run_depletion(
        tmp_path,
        read_structure="interleaved",
        contaminant_mate=contaminant_mate,
    )

    assert output == ""


def test_depletion_keeps_direct_file_filtering_for_single_reads(tmp_path: Path) -> None:
    """Single reads retain Deacon's native direct gzip input path."""
    output, command = run_depletion(tmp_path, read_structure="single")

    assert "@single\n" in output
    assert "reads.fastq.gz" in command
    assert "--interleaved" not in command
    assert "--check-pairs" not in command
    assert "gzip -dc" not in command


def test_interleaved_depletion_uses_explicit_deacon_pair_mode(tmp_path: Path) -> None:
    """Interleaved files use Deacon's native compressed-input interface."""
    _output, command = run_depletion(tmp_path, read_structure="interleaved")

    assert "--interleaved" in command
    assert "gzip -dc" not in command
    assert " - -" not in command


def test_interleaved_depletion_forwards_opt_in_pair_validation(
    tmp_path: Path,
) -> None:
    """NVD only asks Deacon to validate names when explicitly enabled."""
    unchecked_dir = tmp_path / "unchecked"
    checked_dir = tmp_path / "checked"
    unchecked_dir.mkdir()
    checked_dir.mkdir()

    _unchecked_output, unchecked_command = run_depletion(
        unchecked_dir,
        read_structure="interleaved",
    )
    _checked_output, checked_command = run_depletion(
        checked_dir,
        read_structure="interleaved",
        check_pairs=True,
        canonical_names=True,
    )

    assert "--check-pairs" not in unchecked_command
    assert "--check-pairs" in checked_command


def test_single_input_enrichment_filters_records_independently(
    tmp_path: Path,
) -> None:
    """A single-input target pass can discard one former mate by itself."""
    reads = tmp_path / "postmerge.fastq.gz"
    with gzip.open(reads, "wt", encoding="utf-8") as handle:
        handle.write(f"@target\n{CONTAMINANT}\n+\n{'I' * len(CONTAMINANT)}\n")
        handle.write(
            f"@nontarget\n{NONCONTAMINANT}\n+\n{'I' * len(NONCONTAMINANT)}\n",
        )

    index = build_deacon_index(tmp_path)
    output = tmp_path / "enriched.fastq.gz"
    summary = tmp_path / "summary.json"
    completed = subprocess.run(  # noqa: S603
        [
            "deacon",
            "filter",
            "--threads",
            "1",
            "--abs-threshold",
            "1",
            "--rel-threshold",
            "0.0",
            "--summary",
            str(summary),
            "--output",
            str(output),
            str(index),
            str(reads),
        ],  # noqa: S607
        check=False,
        capture_output=True,
        text=True,
        timeout=45,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        enriched_reads = handle.read()
    assert "@target\n" in enriched_reads
    assert "@nontarget\n" not in enriched_reads

    summary_data = json.loads(summary.read_text(encoding="utf-8"))
    assert summary_data["seqs_in"] == 2
    assert summary_data["seqs_out"] == 1
