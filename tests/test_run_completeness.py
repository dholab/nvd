"""Reconciling the samplesheet against the run's results.

Several NVD processes use ``errorStrategy 'ignore'``, and an ignored Nextflow
task emits nothing into its output channel. A sample can therefore disappear
between the samplesheet and the final outputs while the run still reports
success. That has happened twice: a crashed assembler (since fixed) and a failed
``FETCH_FASTQ`` download (still possible). These tests cover the backstop that
names such samples instead of letting them vanish quietly.
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

# Mirrors the wiring in workflows/nvd_main.nf, including the one-element list
# wrapping. Without that wrapping `combine` flattens the two rosters into a
# single list and the diff silently compares garbage, so the shape is part of
# what is being tested here rather than an incidental detail.
SCRIPT = """\
nextflow.enable.dsl = 2

workflow {
    ch_expected = channel.fromList(params.expected.tokenize(','))
        .collect()
        .map { sample_ids -> [sample_ids] }
    // NONE stands for "no samples reached the results". An empty --observed ""
    // cannot be used: Nextflow turns an empty-valued param into boolean true,
    // which then has no tokenize method.
    def observed_ids = params.observed == 'NONE' ? [] : params.observed.tokenize(',')
    ch_observed = (observed_ids ? channel.fromList(observed_ids) : channel.empty())
        .collect()
        .ifEmpty([])
        .map { sample_ids -> [sample_ids] }

    ch_expected.combine(ch_observed).subscribe { expected, observed ->
        def missing = NvdUtils.missingSampleIds(expected, observed)
        println "EXPECTED=${expected.size()} OBSERVED=${observed.size()} MISSING=${missing.join('|')}"
    }
}
"""


def run_completeness(expected: str, observed: str) -> str:
    """Run the helper inside Nextflow and return the printed summary line."""
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
            "--expected",
            expected,
            "--observed",
            observed,
        ],
        cwd=workflow_dir,
        env={**os.environ, "NXF_ANSI_LOG": "false"},
        text=True,
        capture_output=True,
        check=False,
    )
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics

    summary = [line for line in completed.stdout.splitlines() if line.startswith("EXPECTED=")]
    assert len(summary) == 1, diagnostics
    return summary[0]


@pytest.fixture(autouse=True)
def _workflow_tmpdir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Give run_completeness a writable directory for the generated workflow."""
    monkeypatch.setenv("PYTEST_CURRENT_TEST_TMPDIR", str(tmp_path))


def test_a_complete_run_reports_nothing_missing() -> None:
    """The normal case must stay quiet, or the warning becomes noise."""
    assert run_completeness("a,b,c", "a,b,c") == "EXPECTED=3 OBSERVED=3 MISSING="


def test_a_dropped_sample_is_named() -> None:
    """This is the case that cost a sample silently in production."""
    assert run_completeness("a,b,c", "a,c") == "EXPECTED=3 OBSERVED=2 MISSING=b"


def test_several_dropped_samples_are_named_in_a_stable_order() -> None:
    """Sorted output keeps the warning text comparable between runs."""
    assert run_completeness("c,b,a", "b") == "EXPECTED=3 OBSERVED=1 MISSING=a|c"


def test_a_run_that_lost_every_sample_still_reports() -> None:
    """The worst case must not collapse to an empty channel and warn nothing.

    An empty observed roster is exactly what a run whose every sample failed
    looks like, and `ifEmpty` is what keeps the check alive for it.
    """
    assert run_completeness("a,b", "NONE") == "EXPECTED=2 OBSERVED=0 MISSING=a|b"


def test_duplicate_input_rows_are_counted_once() -> None:
    """A sample id repeated on the samplesheet must not be reported twice."""
    assert run_completeness("a,a,b", "a") == "EXPECTED=3 OBSERVED=1 MISSING=b"


def test_an_unexpected_extra_result_is_not_an_error() -> None:
    """Observing more than was asked for is not this check's business.

    Only disappearance matters here; an extra id should not be reported as
    missing, and must not make the check throw.
    """
    assert run_completeness("a", "a,b") == "EXPECTED=1 OBSERVED=2 MISSING="
