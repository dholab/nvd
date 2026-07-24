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
PREPROCESS_READS = ROOT / "subworkflows" / "preprocess_reads"

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None,
    reason="target-enrichment behavior requires Nextflow",
)


def write_fake_deacon(bin_dir: Path, *, seqs_out: int = 1) -> None:
    """Install a deterministic Deacon stand-in for process-boundary tests."""
    executable = bin_dir / "deacon"
    executable.write_text(
        f"""#!/usr/bin/env python3
import gzip
import json
import sys

arguments = sys.argv[1:]
summary = arguments[arguments.index("--summary") + 1]
output = arguments[arguments.index("--output") + 1]
with gzip.open(output, "wt") as stream:
    if {seqs_out}:
        stream.write("@retained\\nACGT\\n+\\nIIII\\n")
with open(summary, "w") as stream:
    json.dump({{"seqs_in": 1, "seqs_out": {seqs_out}, "seqs_removed": {1 - seqs_out}}}, stream)
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
        ["nextflow", "run", "-profile", "test", str(workflow)],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=45,
    )

    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    target_results = results / "nvd" / "01_target_enrichment"
    assert (target_results / "reads" / "sample_A.target_enriched.fastq.gz").is_file()
    assert (target_results / "summaries" / "sample_A.deacon_filter.json").is_file()


def test_zero_target_reads_stop_before_preprocessing(tmp_path: Path) -> None:
    """A sample with no retained target reads emits no preprocessing batch."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_deacon(bin_dir, seqs_out=0)

    lib_dir = tmp_path / "lib"
    lib_dir.mkdir()
    shutil.copy2(ROOT / "lib" / "NvdUtils.groovy", lib_dir)
    shutil.copy2(ROOT / "lib" / "NvdReporting.groovy", lib_dir)
    assets_dir = tmp_path / "assets"
    assets_dir.mkdir()
    (assets_dir / "README.md").write_text("test placeholder", encoding="utf-8")

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
params.skip_fastqc = true
params.merge_pairs = false
params.virus_index = '{target_index}'
params.virus_index_url = null
params.virus_reference_fasta = null
params.virus_abs_threshold = 1
params.virus_rel_threshold = 0.0
params.dedup = false
params.dedup_seq = false
params.trim_adapters = false
params.host_index = null
params.host_index_url = null
params.host_contaminants_fasta = null
params.filter_reads = false
params.filter_low_complexity_reads = false
params.min_read_quality_illumina = 20
params.min_read_quality_nanopore = 12
params.min_read_length = 50
params.min_consecutive_bases = 200

include {{ PREPROCESS_READS }} from '{PREPROCESS_READS}'

workflow {{
    PREPROCESS_READS(Channel.of(tuple(
        [
            id: 'sample_A',
            platform: 'illumina',
            read_mode: 'single',
            r1_count: 1,
            deacon_read_structure: 'single',
        ],
        [file('{reads}')],
    )))

    PREPROCESS_READS.out.reads.view {{ sample_id, _platform, _structure, output ->
        "READ_OUTPUT:${{sample_id}}:${{output.size()}}"
    }}
    PREPROCESS_READS.out.read_counts.view {{ sample_id, count ->
        "INPUT_READS:${{sample_id}}:${{count}}"
    }}
    PREPROCESS_READS.out.complete_empty_samples.view {{ sample_id, platform ->
        "COMPLETE_EMPTY:${{sample_id}}:${{platform}}"
    }}
}}
""",
        encoding="utf-8",
    )

    environment = os.environ.copy()
    environment["PATH"] = f"{bin_dir}:{environment['PATH']}"
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
    assert "INPUT_READS:sample_A:1" in completed.stdout
    assert "COMPLETE_EMPTY:sample_A:illumina" in completed.stdout
    assert "READ_OUTPUT:sample_A" not in completed.stdout
