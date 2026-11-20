"""skip_contig_filter sends every contig past the step-one screen."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None,
    reason="index selection runs inside Nextflow's Groovy runtime",
)

SCRIPT = """\
nextflow.enable.dsl = 2

workflow {
    def chosen = NvdUtils.contigScreeningIndex(params, file(params.step_one_index), params.project_dir)
    println "INDEX ${chosen.name}"
}
"""


def run_index_choice(
    tmp_path: Path, *parameters: str
) -> subprocess.CompletedProcess[str]:
    """Defaults live in a config: a CLI value of `null` would arrive as the string 'null'."""
    lib = tmp_path / "lib"
    lib.mkdir(exist_ok=True)
    shutil.copy2(ROOT / "lib" / "NvdUtils.groovy", lib / "NvdUtils.groovy")
    workflow = tmp_path / "main.nf"
    workflow.write_text(SCRIPT, encoding="utf-8")
    step_one = tmp_path / "step_one.idx"
    step_one.touch()
    config = tmp_path / "harness.config"
    config.write_text(
        f"""\
params.skip_contig_filter = null
params.step_one_index = '{step_one}'
params.project_dir = '{ROOT}'
""",
        encoding="utf-8",
    )
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


def test_default_keeps_the_step_one_index(tmp_path: Path) -> None:
    completed = run_index_choice(tmp_path)
    assert completed.returncode == 0, completed.stderr
    assert "INDEX step_one.idx" in completed.stdout


def test_skip_substitutes_the_empty_index(tmp_path: Path) -> None:
    completed = run_index_choice(tmp_path, "--skip_contig_filter")
    assert completed.returncode == 0, completed.stderr
    assert "INDEX empty_deacon.k31w1.idx" in completed.stdout
