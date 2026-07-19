"""Focused contract tests for single-read quality filtering."""

from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BBMAP_MODULE = ROOT / "modules" / "bbmap"


def test_single_read_filter_retains_passing_former_mate(tmp_path: Path) -> None:
    reads = tmp_path / "former_mates.fastq.gz"
    with gzip.open(reads, "wt", encoding="utf-8") as handle:
        handle.write("@pair/1\nACGT\n+\nIIII\n")
        handle.write("@pair/2\nACGT\n+\n!!!!\n")

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.min_read_length = 1
params.max_read_length = null
params.filter_reads = true
params.filter_low_complexity_reads = false
params.min_read_quality_illumina = 20
params.min_read_quality_nanopore = 20

include {{ FILTER_READS }} from '{BBMAP_MODULE}'

workflow {{
    meta = [id: 'sample_A', platform: 'illumina', read_structure: 'single', query_class: 'single_read']
    FILTER_READS(Channel.of(tuple(meta, file('{reads}'))))
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
        timeout=30,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    outputs = list(
        (tmp_path / "work").glob("**/sample_A.single_read.filtered.fastq.gz"),
    )
    assert len(outputs) == 1, outputs
    with gzip.open(outputs[0], "rt", encoding="utf-8") as handle:
        assert handle.read() == "@pair/1\nACGT\n+\nIIII\n"
