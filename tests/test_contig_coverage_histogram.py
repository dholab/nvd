"""Nextflow contract test for the plain-text contig coverage histogram."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
SAMTOOLS_MODULE = ROOT / "modules" / "samtools"
CONTIG_ID = "nvdContigQuery_water_000001"


def write_indexed_bam(path: Path) -> Path:
    """Write two separated, fully covered 20-base intervals on a 100-base contig."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG_ID, "LN": 100}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for index, start in enumerate((0, 50), start=1):
            read = pysam.AlignedSegment()
            read.query_name = f"read-{index}"
            read.query_sequence = "A" * 20
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = ((0, 20),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 20)
            bam.write(read)

    pysam.index(str(path))
    return Path(f"{path}.bai")


def test_process_renders_unicode_coverage_for_stable_contig_id(
    tmp_path: Path,
) -> None:
    """The report should preserve NVD IDs and the fixture's separated coverage pattern."""
    nextflow = shutil.which("nextflow")
    assert nextflow is not None

    bam = tmp_path / "water.filtered.bam"
    bai = write_indexed_bam(bam)
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ RENDER_CONTIG_COVERAGE_HISTOGRAM }} from '{SAMTOOLS_MODULE}'

workflow {{
    RENDER_CONTIG_COVERAGE_HISTOGRAM(Channel.of(tuple(
        'water',
        file('{bam}'),
        file('{bai}'),
    )))
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    completed = subprocess.run(  # noqa: S603
        [nextflow, "-C", "/dev/null", "run", str(workflow)],
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics

    reports = list((tmp_path / "work").rglob("*.contig_coverage.histogram.txt"))
    assert len(reports) == 1, diagnostics
    report = reports[0].read_text(encoding="utf-8")

    assert CONTIG_ID in report
    assert "Number of reads: 2" in report
    assert "Covered bases:   40bp" in report
    assert "Percent covered: 40%" in report
    assert "Mean coverage:   0.4x" in report
    assert "Histo bin width: 1bp" in report
    assert "Histo max bin:   100%" in report

    # At one base per bin, the two reads form separate 20-base plateaus with
    # 30-base uncovered intervals between and after them.
    expected_bins = "█" * 20 + " " * 30 + "█" * 20 + " " * 30
    histogram_rows = {
        line.split("│")[1] for line in report.splitlines() if line.startswith(">")
    }
    assert histogram_rows == {expected_bins}
