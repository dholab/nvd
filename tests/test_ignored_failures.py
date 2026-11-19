"""Making ignored task failures audible at the end of a run.

Nearly every NVD process ends its retry budget at ``errorStrategy 'ignore'``.
That is the right policy, but Nextflow then reports the run as a success, so the
failure is invisible. These tests cover the run-end warning that reads
Nextflow's own count of ignored tasks.

This is deliberately separate from the per-sample completeness check: a failure
can lose a run-level product without losing a sample, which is exactly what
silently cost the MultiQC report on the scheduled runs.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None,
    reason="the helper runs inside Nextflow's Groovy runtime",
)

# NONE stands for "the helper returned null", so a quiet run is distinguishable
# from a run that produced an empty message.
SCRIPT = """\
nextflow.enable.dsl = 2

workflow {
    def message = NvdUtils.describeIgnoredFailures(
        params.ignored as int,
        params.succeeded as int,
        params.report,
    )
    println "RESULT=${message == null ? 'NONE' : message}"
}
"""


def describe_ignored_failures(ignored: int, succeeded: int) -> str:
    """Run the helper inside Nextflow and return what it produced."""
    nextflow = shutil.which("nextflow")
    assert nextflow is not None

    workflow_dir = Path(os.environ["PYTEST_CURRENT_TEST_TMPDIR"])
    workflow = workflow_dir / "main.nf"
    workflow.write_text(SCRIPT, encoding="utf-8")
    lib = workflow_dir / "lib"
    lib.mkdir(exist_ok=True)
    shutil.copy2(ROOT / "lib" / "NvdUtils.groovy", lib / "NvdUtils.groovy")

    completed = subprocess.run(  # noqa: S603
        [
            nextflow,
            "-C",
            "/dev/null",
            "run",
            str(workflow),
            "--ignored",
            str(ignored),
            "--succeeded",
            str(succeeded),
            "--report",
            "results/report.html",
        ],
        cwd=workflow_dir,
        env={**os.environ, "NXF_ANSI_LOG": "false"},
        text=True,
        capture_output=True,
        check=False,
    )
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics

    results = [line for line in completed.stdout.splitlines() if line.startswith("RESULT=")]
    assert len(results) == 1, diagnostics
    return results[0].removeprefix("RESULT=")


@pytest.fixture(autouse=True)
def _workflow_tmpdir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Give the helper a writable directory for the generated workflow."""
    monkeypatch.setenv("PYTEST_CURRENT_TEST_TMPDIR", str(tmp_path))


def test_a_clean_run_says_nothing() -> None:
    """A healthy run must stay quiet, or the warning becomes noise."""
    assert describe_ignored_failures(0, 120) == "NONE"


def test_one_ignored_failure_is_reported_and_points_at_the_report() -> None:
    """This is the case that silently cost the MultiQC report in production."""
    message = describe_ignored_failures(1, 119)

    assert message != "NONE"
    assert "1 task failed" in message
    # The reader has to be able to act on the log line alone.
    assert "results/report.html" in message


def test_the_counts_are_both_reported() -> None:
    """A denominator tells the reader whether the run mostly worked."""
    message = describe_ignored_failures(3, 117)

    assert "3 tasks failed" in message
    assert "117 succeeded" in message


def test_the_task_noun_agrees_with_the_count() -> None:
    """Small thing, but the warning is read by people under time pressure."""
    assert "1 task failed" in describe_ignored_failures(1, 10)
    assert "2 tasks failed" in describe_ignored_failures(2, 10)


def test_a_run_where_everything_was_ignored_still_reports() -> None:
    """The worst case must not be mistaken for a clean run."""
    message = describe_ignored_failures(14, 0)

    assert "14 tasks failed" in message
    assert "0 succeeded" in message
