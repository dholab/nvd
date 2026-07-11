"""Focused behavioral tests for target-enrichment routing and publishing."""

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
    shutil.which("nextflow") is None,
    reason="target-enrichment behavior requires Nextflow",
)


def write_fake_deacon(bin_dir: Path) -> None:
    """Install a deterministic Deacon stand-in for process-boundary tests."""
    executable = bin_dir / "deacon"
    executable.write_text(
        """#!/usr/bin/env python3
import gzip
import json
import sys

arguments = sys.argv[1:]
summary = arguments[arguments.index("--summary") + 1]
output = arguments[arguments.index("--output") + 1]
with gzip.open(output, "wt") as stream:
    stream.write("@retained\\nACGT\\n+\\nIIII\\n")
with open(summary, "w") as stream:
    json.dump({"seqs_in": 1, "seqs_out": 1, "seqs_removed": 0}, stream)
""",
        encoding="utf-8",
    )
    executable.chmod(0o755)


def test_target_enrichment_outputs_are_published(tmp_path: Path) -> None:
    """Configured target enrichment publishes retained reads and its summary."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_deacon(bin_dir)

    reads = tmp_path / "reads.fastq.gz"
    with gzip.open(reads, "wt", encoding="utf-8") as handle:
        handle.write("@source\nACGT\n+\nIIII\n")
    target_index = tmp_path / "target.idx"
    target_index.touch()

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

params.no_enrichment = false
params.virus_index = '{target_index}'
params.virus_index_url = null
params.virus_reference_fasta = null
params.virus_abs_threshold = 1
params.virus_rel_threshold = 0.0

include {{ DEACON_ENRICH_TARGET_READS }} from '{DEACON_MODULE}'

workflow {{
    DEACON_ENRICH_TARGET_READS(Channel.of(tuple(
        [
            id: 'sample_A',
            platform: 'illumina',
            read_mode: 'single',
            r1_count: 1,
            deacon_read_structure: 'single',
        ],
        [file('{reads}')],
        file('{target_index}'),
        true,
    )))
}}
""",
        encoding="utf-8",
    )

    results = tmp_path / "results"
    config_text = (
        (ROOT / "nextflow.config")
        .read_text(encoding="utf-8")
        .replace(
            "${projectDir}/conf/results.config",
            str(ROOT / "conf" / "results.config"),
        )
        .replace(
            'results                   = "${launchDir}/results"',
            f"results                   = '{results}'",
        )
        .replace(
            "virus_index               = null",
            f"virus_index               = '{target_index}'",
        )
    )
    (tmp_path / "nextflow.config").write_text(
        config_text,
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["PATH"] = f"{bin_dir}:{environment['PATH']}"
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "run", str(workflow)],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=45,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    published = results / "nvd" / "01_virus_enriched_reads"
    assert (published / "sample_A.human_virus.fastq.gz").is_file()
    assert (published / "sample_A.deacon_filter.json").is_file()
