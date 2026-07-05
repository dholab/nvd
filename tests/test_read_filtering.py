"""Focused behavioral tests for read quality, length, and entropy filtering."""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BBMAP_MODULE = ROOT / "modules" / "bbmap"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None or shutil.which("bbduk.sh") is None,
    reason="FILTER_READS behavior requires Nextflow and locked BBTools",
)

HIGH_COMPLEXITY = "ACGTTGCAGATCGTAGCTAGGCTAACGATCGTACGGCATTCGATGCTAGCATCGTACGTA"
LOW_COMPLEXITY = "T" * len(HIGH_COMPLEXITY)
DINUCLEOTIDE_REPEAT = "AC" * (len(HIGH_COMPLEXITY) // 2)


def write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    """Write named sequence and quality triples as gzip-compressed FASTQ."""
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for name, sequence, quality in records:
            handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


def run_filter(
    tmp_path: Path,
    records: list[tuple[str, str, str]],
    *,
    filter_reads: bool,
    filter_low_complexity_reads: bool,
    read_structure: str = "single",
    min_read_length: int = 50,
) -> str:
    """Run the public FILTER_READS Nextflow process and return its FASTQ text."""
    reads = tmp_path / "reads.fastq.gz"
    write_fastq(reads, records)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.filter_reads = {str(filter_reads).lower()}
params.filter_low_complexity_reads = {str(filter_low_complexity_reads).lower()}
params.min_read_entropy = 0.9
params.min_read_length = {min_read_length}
params.max_read_length = null

include {{ FILTER_READS }} from '{BBMAP_MODULE}'

workflow {{
    FILTER_READS(Channel.of(tuple('sample_A', 'illumina', '{read_structure}', 'single_read', file('{reads}'), 20)))
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
        (tmp_path / "work").glob("**/sample_A.single_read.filtered.fastq.gz")
    )
    assert len(outputs) == 1, outputs
    with gzip.open(outputs[0], "rt", encoding="utf-8") as handle:
        return handle.read()


def test_entropy_filter_does_not_enable_quality_filtering(tmp_path: Path) -> None:
    """Entropy-only filtering keeps a complex read despite low base quality."""
    output = run_filter(
        tmp_path,
        [
            ("complex-low-quality", HIGH_COMPLEXITY, "!" * len(HIGH_COMPLEXITY)),
            ("homopolymer-high-quality", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
            (
                "dinucleotide-high-quality",
                DINUCLEOTIDE_REPEAT,
                "I" * len(DINUCLEOTIDE_REPEAT),
            ),
        ],
        filter_reads=False,
        filter_low_complexity_reads=True,
        min_read_length=len(HIGH_COMPLEXITY) + 1,
    )

    assert "@complex-low-quality\n" in output
    assert "@homopolymer-high-quality\n" not in output
    assert "@dinucleotide-high-quality\n" not in output


def test_quality_filter_does_not_enable_entropy_filtering(tmp_path: Path) -> None:
    """Quality-only filtering keeps a high-quality low-complexity read."""
    output = run_filter(
        tmp_path,
        [
            ("complex-low-quality", HIGH_COMPLEXITY, "!" * len(HIGH_COMPLEXITY)),
            ("homopolymer-high-quality", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
        ],
        filter_reads=True,
        filter_low_complexity_reads=False,
    )

    assert "@complex-low-quality\n" not in output
    assert "@homopolymer-high-quality\n" in output


def test_combined_filters_can_emit_valid_empty_fastq(tmp_path: Path) -> None:
    """An all-rejected batch remains a valid empty gzip output."""
    output = run_filter(
        tmp_path,
        [
            ("complex-low-quality", HIGH_COMPLEXITY, "!" * len(HIGH_COMPLEXITY)),
            ("homopolymer-high-quality", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
        ],
        filter_reads=True,
        filter_low_complexity_reads=True,
    )

    assert output == ""


def test_entropy_filter_removes_complete_interleaved_pair(tmp_path: Path) -> None:
    """An entropy failure removes both mates from an intact interleaved pair."""
    output = run_filter(
        tmp_path,
        [
            ("pair/1", HIGH_COMPLEXITY, "I" * len(HIGH_COMPLEXITY)),
            ("pair/2", LOW_COMPLEXITY, "I" * len(LOW_COMPLEXITY)),
        ],
        filter_reads=False,
        filter_low_complexity_reads=True,
        read_structure="interleaved",
    )

    assert output == ""


def test_default_entropy_filters_observed_repeat_rich_reads(tmp_path: Path) -> None:
    """The default rejects seven examples but preserves a mostly complex read."""
    observed = [
        (
            "000110",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        ),
        (
            "000113",
            "TCCTCTCTCTCTTCCTCTCTCTCTTCCTCTTCTACTTTTGTTGACGACGAAGAAATCGAGAGACAGATACGAGAGGAGCTGAGGTGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA",
        ),
        (
            "000138",
            "TAATCCGCCCACACACACACATACACATACACACACCCACACCCACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
        ),
        (
            "000139",
            "TGGCTGATGGCGTCATGCGTTCGTTTGTTATTTTGTTGTTGTTTTTCGTTGTTTTCTTTTGTTGGTTCAAGGATACGATGGTGGTCTTCCTTGTGTTTTTCTCTTTTTTTCTTTTTTTCTTTTTGTTTTCATCTTTTTTTTTTTTTT",
        ),
        (
            "000141",
            "TTGAAAGCAAGACGAAGCAAAACAGAATGACGGGAAAAGAGGCTCTTCCAAGAAGAAACCAGGAAAAACGCACGTGGTTTGTTGCTTTTTTCGTTTTGGCGCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        ),
        (
            "000144",
            "TGTGTGATGTGTGCGTACGTCTGTGCGTCTGACATCCCTAGGGTCGGCATCAGCCTGTGCTTGCCCGCCGCGTGCGCCCAGGACATCGAGGAGGGCTACGAGTTCCTCCTCCCTTACCTCGCGTACGCACACACACACACACACACA",
        ),
        (
            "000145",
            "TCCCGTCCATTATCACTCATTTTTTGCTCATTGTTGTTGTTGTTGTTGCTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG",
        ),
        (
            "000146",
            "TCGGAGAAATTTATTGGCTAAATGCCAGTATAAAATGGGCGAAGTGTCAAACCTTAAAAAGATAGAAACGAATGAAAAAAGGAAAGAACTGACAGTTGAGTGAAAGAATGAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        ),
        (
            "control",
            "ACGTTGCAGATCGTAGCTAGGCTAACGATCGTACGGCATTCGATGCTAGCATCGTACGTATCGATGGCATACGTAGCTGATCGTACGATCGTAGCTAGCATGCTACGATCGTACGTAGCATCGATGCTAGCATCGTAGCTAGCATGCA",
        ),
    ]
    output = run_filter(
        tmp_path,
        [(name, sequence, "I" * len(sequence)) for name, sequence in observed],
        filter_reads=False,
        filter_low_complexity_reads=True,
    )
    retained_ids = {
        line.removeprefix("@")
        for index, line in enumerate(output.splitlines())
        if index % 4 == 0
    }

    assert retained_ids == {"000144", "control"}


def test_either_filter_family_schedules_filter_reads() -> None:
    """Either read-filter family should schedule the shared FILTER_READS process."""
    preprocessing = (ROOT / "subworkflows" / "preprocess_reads.nf").read_text(
        encoding="utf-8",
    )

    assert (
        "params.filter_reads || params.filter_low_complexity_reads" in preprocessing
    )
