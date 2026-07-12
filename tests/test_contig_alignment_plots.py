"""Nextflow contract test for interactive per-contig alignment reports."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
REPORTING_MODULE = ROOT / "modules" / "reporting"
MAPPED_CONTIGS = (
    "nvdContigQuery_water_000001",
    "nvdContigQuery_water_000002",
)
UNMAPPED_CONTIG = "nvdContigQuery_water_000003"
CONTIGS = (*MAPPED_CONTIGS, UNMAPPED_CONTIG)


def write_reference(path: Path) -> None:
    """Write three stable NVD contigs of equal length."""
    with path.open("w", encoding="utf-8") as handle:
        for contig in CONTIGS:
            handle.write(f">{contig}\n{'A' * 100}\n")


def write_indexed_bam(path: Path) -> Path:
    """Write one mapped read on each of the first two reference contigs."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": 100} for contig in CONTIGS],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for reference_id, contig in enumerate(MAPPED_CONTIGS):
            read = pysam.AlignedSegment()
            read.query_name = f"read-for-{contig}"
            read.query_sequence = "A" * 20
            read.flag = 0
            read.reference_id = reference_id
            read.reference_start = 10
            read.mapping_quality = 60
            read.cigar = ((0, 20),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 20)
            bam.write(read)

    pysam.index(str(path))
    return Path(f"{path}.bai")


def test_process_renders_self_contained_html_for_covered_contigs_only(
    tmp_path: Path,
) -> None:
    """Mapped stable IDs should get offline HTML reports; empty references should not."""
    nextflow = shutil.which("nextflow")
    assert nextflow is not None

    reference = tmp_path / "water.target_filtered.fasta"
    write_reference(reference)
    bam = tmp_path / "water.filtered.bam"
    bai = write_indexed_bam(bam)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ RENDER_CONTIG_ALIGNMENT_PLOTS }} from '{REPORTING_MODULE}'

workflow {{
    RENDER_CONTIG_ALIGNMENT_PLOTS(Channel.of(tuple(
        'water',
        file('{reference}'),
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

    reports = list((tmp_path / "work").rglob("*.alignoth.html"))
    assert {report.name for report in reports} == {
        f"{contig}.alignoth.html" for contig in MAPPED_CONTIGS
    }
    assert {report.parent.name for report in reports} == {"water.alignoth"}
    assert all(UNMAPPED_CONTIG not in report.name for report in reports)

    for report in reports:
        html = report.read_text(encoding="utf-8")
        assert "<html" in html.lower()
        assert "vegaEmbed" in html
        assert "decompressFromUTF16" in html
        assert "cdn.jsdelivr.net" not in html
        assert 'src="http' not in html
