"""Step-one Deacon filter policy: enrichment, background depletion, or passthrough."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None,
    reason="policy helpers run inside Nextflow's Groovy runtime",
)

SCRIPT = """\
nextflow.enable.dsl = 2

workflow {
    if (params.validate) {
        NvdUtils.validateStepOneFilter(params)
    }
    def policy = NvdUtils.stepOneFilterPolicy(params)
    println "POLICY enrich=${policy.target_enrichment_enabled} background=${policy.background_depletion_enabled} abs=${policy.abs_threshold} rel=${policy.rel_threshold}"
}
"""

BASE_CONFIG = """\
params.virus_index = null
params.virus_index_url = null
params.virus_reference_fasta = null
params.virus_abs_threshold = 1
params.virus_rel_threshold = 0.0
params.background_index = null
params.background_abs_threshold = 1
params.background_rel_threshold = 0.0
params.host_abs_threshold = 2
params.host_rel_threshold = 0.01
params.no_enrichment = false
params.validate = false
"""


def write_harness(tmp_path: Path, script: str) -> tuple[Path, Path]:
    """Stage NvdUtils, the script, and a config carrying the pipeline defaults."""
    lib = tmp_path / "lib"
    lib.mkdir(exist_ok=True)
    shutil.copy2(ROOT / "lib" / "NvdUtils.groovy", lib / "NvdUtils.groovy")
    workflow = tmp_path / "main.nf"
    workflow.write_text(script, encoding="utf-8")
    config = tmp_path / "harness.config"
    config.write_text(BASE_CONFIG, encoding="utf-8")
    return workflow, config


def run_policy(tmp_path: Path, *parameters: str) -> subprocess.CompletedProcess[str]:
    """Nextflow parses a CLI value of `null` as the string 'null', so defaults live in the config."""
    workflow, config = write_harness(tmp_path, SCRIPT)
    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    return subprocess.run(  # noqa: S603
        ["nextflow", "-C", str(config), "run", str(workflow), *parameters],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=90,
    )


def test_virus_index_selects_enrichment(tmp_path: Path) -> None:
    virus = tmp_path / "virus.idx"
    virus.touch()
    completed = run_policy(tmp_path, "--virus_index", str(virus))
    assert completed.returncode == 0, completed.stderr
    assert "POLICY enrich=true background=false abs=1 rel=0.0" in completed.stdout


def test_background_index_selects_depletion_with_its_thresholds(tmp_path: Path) -> None:
    background = tmp_path / "background.idx"
    background.touch()
    completed = run_policy(tmp_path, "--background_index", str(background))
    assert completed.returncode == 0, completed.stderr
    assert "POLICY enrich=false background=true abs=1 rel=0.0" in completed.stdout


def test_neither_index_is_passthrough(tmp_path: Path) -> None:
    completed = run_policy(tmp_path)
    assert completed.returncode == 0, completed.stderr
    assert "POLICY enrich=false background=false" in completed.stdout


def test_no_enrichment_lets_a_preset_virus_index_coexist(tmp_path: Path) -> None:
    virus = tmp_path / "virus.idx"
    background = tmp_path / "background.idx"
    virus.touch()
    background.touch()
    completed = run_policy(
        tmp_path,
        "--virus_index",
        str(virus),
        "--background_index",
        str(background),
        "--no_enrichment",
        "true",
        "--validate",
        "true",
    )
    assert completed.returncode == 0, completed.stderr
    assert "POLICY enrich=false background=true abs=1 rel=0.0" in completed.stdout


def test_validation_stops_on_conflict(tmp_path: Path) -> None:
    virus = tmp_path / "virus.idx"
    background = tmp_path / "background.idx"
    virus.touch()
    background.touch()
    completed = run_policy(
        tmp_path,
        "--virus_index",
        str(virus),
        "--background_index",
        str(background),
        "--validate",
        "true",
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "background_index" in output
    assert "virus_index" in output
    assert "--no_enrichment" in output


def test_validation_rejects_bare_flag(tmp_path: Path) -> None:
    """`--background_index` with no value becomes boolean true in Nextflow."""
    completed = run_policy(tmp_path, "--background_index", "--validate", "true")
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "background_index must be a path" in output


def test_validation_rejects_missing_file(tmp_path: Path) -> None:
    missing = tmp_path / "absent.idx"
    completed = run_policy(
        tmp_path, "--background_index", str(missing), "--validate", "true"
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "does not exist" in output
    assert str(missing) in output


CONTIG_SCRIPT = """\
nextflow.enable.dsl = 2

workflow {
    def policy = NvdUtils.contigFilterPolicy(params, params.use_depletion)
    println "CONTIG enrich=${policy.target_enrichment_enabled} abs=${policy.target_abs_threshold} rel=${policy.target_rel_threshold} dep=${policy.depletion_enabled} dabs=${policy.depletion_abs_threshold}"
}
"""


def run_contig_policy(
    tmp_path: Path, *parameters: str
) -> subprocess.CompletedProcess[str]:
    workflow, config = write_harness(tmp_path, CONTIG_SCRIPT)
    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    return subprocess.run(  # noqa: S603
        ["nextflow", "-C", str(config), "run", str(workflow), *parameters],  # noqa: S607
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
        timeout=90,
    )


def test_contig_filter_uses_background_thresholds_in_background_mode(
    tmp_path: Path,
) -> None:
    background = tmp_path / "background.idx"
    background.touch()
    completed = run_contig_policy(
        tmp_path,
        "--background_index",
        str(background),
        "--background_abs_threshold",
        "4",
        "--use_depletion",
        "false",
    )
    assert completed.returncode == 0, completed.stderr
    assert "CONTIG enrich=false abs=4 rel=0.0 dep=false dabs=null" in completed.stdout


def test_contig_filter_keeps_virus_thresholds_and_host_depletion_in_enrichment_mode(
    tmp_path: Path,
) -> None:
    virus = tmp_path / "virus.idx"
    virus.touch()
    completed = run_contig_policy(
        tmp_path,
        "--virus_index",
        str(virus),
        "--use_depletion",
        "true",
    )
    assert completed.returncode == 0, completed.stderr
    assert "CONTIG enrich=true abs=1 rel=0.0 dep=true dabs=2" in completed.stdout
