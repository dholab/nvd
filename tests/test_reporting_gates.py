"""Pin which reporting stages run by default and which --skip_big_tables removes."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REPORTING = ROOT / "subworkflows" / "reporting.nf"
RESULTS_CONFIG = ROOT / "conf" / "results.config"

STILL_EXPERIMENTAL = {
    "SOURMASH_SKETCH_QUERY_METAGENOME",
    "SOURMASH_COLLECT_QUERY_SKETCHES",
    "SOURMASH_COMPARE_QUERY_SKETCHES",
    "COMPUTE_SAMPLE_SKETCH_DISTANCES",
    "COMPUTE_SAMPLE_ORDINATION",
    "REPORT_POSSIBLE_SAMPLE_MIXUPS",
    "SUMMARIZE_SAMPLE_SIMILARITY",
    "PLOT_SAMPLE_ORDINATION",
    "ASSEMBLE_WITH_MYLOASM",
    "ASSEMBLE_WITH_METAMDBG",
    "ASSEMBLE_WITH_METAFLYE",
    "REPORT_LONG_READ_ASSEMBLY_ELIGIBILITY",
    "EXTRACT_UNIQUE_CONTIGS",
}


def publish_blocks(text: str) -> dict[str, str]:
    """Map each withName selector to the text of its publishDir block."""
    blocks: dict[str, str] = {}
    for match in re.finditer(
        r"withName:\s*'([^']+)'\s*\{(.*?)\n    \}", text, re.DOTALL
    ):
        blocks[match.group(1)] = match.group(2)
    return blocks


def test_reporting_subworkflow_only_forwards_experimental_to_multiqc() -> None:
    """The invitation flag is the one surviving reference; nothing is gated on it."""
    text = REPORTING.read_text(encoding="utf-8")
    assert text.count("params.experimental") == 1
    assert "channel.value(params.experimental == true)," in text


def test_only_similarity_qc_and_long_read_assembly_publish_behind_experimental() -> (
    None
):
    blocks = publish_blocks(RESULTS_CONFIG.read_text(encoding="utf-8"))
    gated = {name for name, body in blocks.items() if "params.experimental" in body}
    assert gated == STILL_EXPERIMENTAL


def test_big_table_publishing_follows_skip_big_tables() -> None:
    blocks = publish_blocks(RESULTS_CONFIG.read_text(encoding="utf-8"))
    for name in (
        "BUILD_QUERY_BIG_TABLE",
        "BUILD_TAXON_BIG_TABLE",
        "CONCATENATE_QUERY_BIG_TABLE",
        "CONCATENATE_TAXON_BIG_TABLE",
    ):
        assert "enabled: !params.skip_big_tables" in blocks[name], name
    assert (
        "enabled: !params.skip_big_tables && !params.skip_blast"
        in blocks["EMIT_BEST_HIT_SEQUENCE_EVIDENCE"]
    )
    for name in (
        "BUILD_SEQUENCE_FLOW",
        "ESTIMATE_CRUMBS_PROFILE",
        "RENDER_CRUMBS_TAXBURST",
        "RENDER_MERGED_CRUMBS_TAXBURST",
    ):
        assert "enabled: true" in blocks[name], name


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="needs Nextflow")
@pytest.mark.parametrize(
    ("skip_big_tables", "expected"),
    [("null", "enabled:true]"), ("true", "enabled:false]")],
)
def test_rendered_config_resolves_big_table_publishing(
    tmp_path: Path, skip_big_tables: str, expected: str
) -> None:
    """The null default publishes the big tables; skip_big_tables = true turns it off.

    Publish rules evaluate their `enabled:` expressions when results.config is
    included, so the params must be set before the include, the way the
    pipeline's own nextflow.config orders them.
    """
    config = tmp_path / "render.config"
    config.write_text(
        f"""\
params.results = '{tmp_path / "results"}'
params.experimental = false
params.skip_unassembled_read_queries = null
params.skip_blast = null
params.skip_big_tables = {skip_big_tables}
params.no_enrichment = false
params.background_index = null
params.virus_index = null
params.virus_index_url = null
params.virus_reference_fasta = null
includeConfig '{RESULTS_CONFIG}'
""",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "-C", str(config), "config", "-flat"],  # noqa: S607
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    line = next(
        line
        for line in completed.stdout.splitlines()
        if line.startswith("process.'withName:BUILD_QUERY_BIG_TABLE'.publishDir")
    )
    assert line.strip().endswith(expected), line


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="needs Nextflow")
@pytest.mark.parametrize(
    ("virus_index", "background_index", "expected"),
    [
        ("'/refs/virus.idx'", "null", "enabled:true"),
        ("null", "'/refs/background.idx'", "enabled:true"),
        ("null", "null", "enabled:false"),
    ],
)
def test_rendered_config_publishes_step_one_outputs_for_either_filter(
    tmp_path: Path, virus_index: str, background_index: str, expected: str
) -> None:
    """Step-one reads, summaries, and the run-level report publish whenever step 1 is a real filter."""
    config = tmp_path / "render.config"
    config.write_text(
        f"""\
params.results = '{tmp_path / "results"}'
params.experimental = false
params.skip_unassembled_read_queries = null
params.skip_blast = null
params.skip_big_tables = null
params.no_enrichment = false
params.virus_index = {virus_index}
params.virus_index_url = null
params.virus_reference_fasta = null
params.background_index = {background_index}
includeConfig '{RESULTS_CONFIG}'
""",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    completed = subprocess.run(  # noqa: S603
        ["nextflow", "-C", str(config), "config", "-flat"],  # noqa: S607
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    for selector in (
        "process.'withName:^(DEACON_ENRICH_TARGET_READS|DEACON_ENRICH_SRA_READS)$'.publishDir",
        "process.'withName:TARGET_ENRICHMENT_REPORT'.publishDir",
    ):
        line = next(
            line for line in completed.stdout.splitlines() if line.startswith(selector)
        )
        assert line.count(expected) == 2, line
