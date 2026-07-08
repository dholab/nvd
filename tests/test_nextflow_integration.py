"""Slow end-to-end tests for the checked-in mini viral fixtures."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
MANIFEST = DATA / "reference.manifest.json"
SAMPLESHEET = DATA / "integration_sra_samplesheet.csv"
DEACON_INDEX = DATA / "mini_virus_deacon.k31w1.idx"
BLAST_DB_PREFIX = "mini_virus_blast"
BLAST_DB = DATA
E2E_OUTPUT_DIR = ROOT / ".e2e"


def write_mini_taxdump(taxonomy_dir: Path) -> None:
    """Write a small taxonomy tree covering the fixture virus taxids."""
    taxonomy_dir.mkdir(parents=True, exist_ok=True)
    (taxonomy_dir / "nodes.dmp").write_text(
        """\
1\t|\t1\t|\tno rank\t|\t
10239\t|\t1\t|\tsuperkingdom\t|\t
10240\t|\t10239\t|\tfamily\t|\t
10242\t|\t10240\t|\tgenus\t|\t
10244\t|\t10242\t|\tspecies\t|\t
10255\t|\t10240\t|\tgenus\t|\t
10258\t|\t10255\t|\tspecies\t|\t
""",
        encoding="utf-8",
    )
    (taxonomy_dir / "names.dmp").write_text(
        """\
1\t|\troot\t|\t\t|\tscientific name\t|
10239\t|\tViruses\t|\t\t|\tscientific name\t|
10240\t|\tPoxviridae\t|\t\t|\tscientific name\t|
10242\t|\tOrthopoxvirus\t|\t\t|\tscientific name\t|
10244\t|\tMonkeypox virus\t|\t\t|\tscientific name\t|
10255\t|\tParapoxvirus\t|\t\t|\tscientific name\t|
10258\t|\tOrf virus\t|\t\t|\tscientific name\t|
""",
        encoding="utf-8",
    )
    (taxonomy_dir / "merged.dmp").write_text("", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest() -> dict[str, Any]:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def verify_fixture_checksums(manifest: dict[str, Any]) -> None:
    missing: list[str] = []
    mismatched: list[str] = []
    for filename, metadata in manifest["files"].items():
        path = DATA / filename
        if not path.exists():
            missing.append(filename)
            continue
        if sha256(path) != metadata["sha256"]:
            mismatched.append(filename)

    assert not missing, "Missing integration fixture files: " + ", ".join(missing)
    assert not mismatched, "Integration fixture checksum mismatch: " + ", ".join(
        mismatched,
    )


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_e2e_run_dir() -> Path:
    output_root = Path(os.environ.get("NVD_E2E_OUTPUT_DIR", E2E_OUTPUT_DIR))
    run_id = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ") + f"-pid{os.getpid()}"
    run_dir = output_root / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    (output_root / "latest.txt").write_text(str(run_dir), encoding="utf-8")
    return run_dir


def run_nextflow() -> tuple[subprocess.CompletedProcess[str], Path]:
    profile = os.environ.get("NVD_INTEGRATION_PROFILE", "test")
    show_progress = os.environ.get("NVD_E2E_SHOW_PROGRESS") == "1"
    run_dir = make_e2e_run_dir()
    results_dir = run_dir / "results"
    work_dir = run_dir / "work"
    taxonomy_dir = run_dir / "taxonomy"
    write_mini_taxdump(taxonomy_dir)
    command = [
        "nextflow",
        "run",
        ".",
        "-profile",
        profile,
        "--samplesheet",
        str(SAMPLESHEET),
        "--virus_index",
        str(DEACON_INDEX),
        "--blast_db",
        str(BLAST_DB),
        "--blast_db_prefix",
        BLAST_DB_PREFIX,
        "--taxonomy_dir",
        str(taxonomy_dir),
        "--experiment_id",
        "1",
        "--results",
        str(results_dir),
        "--work_dir",
        str(work_dir),
    ]
    completed = subprocess.run(  # noqa: S603
        command,
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=not show_progress,
        timeout=60 * 60 * 2,
    )
    (run_dir / "nextflow.stdout.log").write_text(
        completed.stdout or "",
        encoding="utf-8",
    )
    (run_dir / "nextflow.stderr.log").write_text(
        completed.stderr or "",
        encoding="utf-8",
    )
    return completed, run_dir


@pytest.mark.slow
@pytest.mark.network
def test_mini_sra_viral_pipeline_completes() -> None:
    """Tiny SRA runs should complete through enrichment, assembly, and BLAST."""
    manifest = load_manifest()
    verify_fixture_checksums(manifest)

    completed, run_dir = run_nextflow()
    assert completed.returncode == 0, (
        f"Nextflow integration run failed.\n\nRun directory: {run_dir}\n\nSTDOUT:\n"
        + (completed.stdout or "<not captured; check terminal output or Nextflow logs>")
        + "\n\nSTDERR:\n"
        + (completed.stderr or "<not captured; check terminal output or Nextflow logs>")
    )

    results_root = run_dir / "results" / "nvd"
    final_dir = results_root / "07_merged_blast_results" / "final"
    final_blast_files = sorted(final_dir.glob("*_blast.final.tsv"))
    assert final_blast_files, f"No final BLAST TSVs found in {final_dir}"

    experiment_blast = results_root / "08_experiment_summary" / "experiment_blast_results.tsv"
    assert experiment_blast.is_file(), f"Missing experiment BLAST summary: {experiment_blast}"

    final_text = "\n".join(path.read_text(encoding="utf-8") for path in final_blast_files)
    experiment_rows = read_tsv_rows(experiment_blast)
    assert experiment_rows, f"No experiment BLAST rows found in {experiment_blast}"

    for run_info in manifest["sra_runs"]:
        sample_rows = [
            row for row in experiment_rows if row.get("sample") == run_info["sample_id"]
        ]
        assert sample_rows, f"No experiment BLAST rows found for {run_info['sample_id']}"
        assert any(row.get("staxids") == str(run_info["taxid"]) for row in sample_rows), (
            f"No {run_info['taxid']} BLAST taxid found for {run_info['sample_id']}"
        )
        assert any(run_info["expected_organism"] in row.get("rank", "") for row in sample_rows), (
            f"No {run_info['expected_organism']} lineage found for "
            f"{run_info['sample_id']}"
        )
        expected_tasks = run_info.get("expected_tasks", [])
        observed_tasks = {row.get("task") for row in sample_rows}
        assert set(expected_tasks) <= observed_tasks, (
            f"Missing expected BLAST tasks for {run_info['sample_id']}: "
            f"expected {sorted(expected_tasks)}, observed {sorted(observed_tasks)}"
        )

    assert "Orf virus" in final_text
    assert "Monkeypox virus" in final_text
