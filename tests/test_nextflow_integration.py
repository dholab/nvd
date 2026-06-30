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
SOURMASH_REF_FASTA = DATA / "mini_virus_reference.fasta"
SOURMASH_LINEAGES = DATA / "mini_sourmash_lineages.csv"
BLAST_DB_PREFIX = "mini_virus_blast"
BLAST_DB = DATA
E2E_OUTPUT_DIR = ROOT / ".e2e"
LOCAL_FASTQ_FIXTURES = {
    "hits_r1": DATA / "hits_only_R1.fastq.gz",
    "hits_r2": DATA / "hits_only_R2.fastq.gz",
    "water_plus_hits_r1": DATA / "water_plus_hits_R1.fastq.gz",
    "water_plus_hits_r2": DATA / "water_plus_hits_R2.fastq.gz",
}
LOCAL_E2E_SAMPLES = (
    {
        "sample_id": "local_hits_exact",
        "expected_organism": "Orf virus",
        "taxid": 10258,
        "source": "paired_files",
        "fastq1": LOCAL_FASTQ_FIXTURES["hits_r1"],
        "fastq2": LOCAL_FASTQ_FIXTURES["hits_r2"],
    },
    {
        "sample_id": "local_hits_glob",
        "expected_organism": "Orf virus",
        "taxid": 10258,
        "source": "paired_globs",
        "fastq1_glob": "local_hits_glob_L*_R1_001.fastq.gz",
        "fastq2_glob": "local_hits_glob_L*_R2_001.fastq.gz",
    },
)


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


def read_delimited_rows(path: Path, *, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    return read_delimited_rows(path, delimiter=",")


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    return read_delimited_rows(path, delimiter="\t")


def make_e2e_run_dir() -> Path:
    output_root = Path(os.environ.get("NVD_E2E_OUTPUT_DIR", E2E_OUTPUT_DIR))
    run_id = datetime.now(tz=UTC).strftime("%Y%m%dT%H%M%SZ") + f"-pid{os.getpid()}"
    run_dir = output_root / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    (output_root / "latest.txt").write_text(str(run_dir), encoding="utf-8")
    return run_dir


def integration_experimental_enabled() -> bool:
    # The pipeline already has an `experimental` parameter; this environment
    # variable only tells the test harness whether to pass that parameter and
    # assert the extra outputs. A pytest option would make that relationship more
    # explicit, but would add plumbing for little gain compared with the existing
    # just recipes and CI entrypoints.
    return os.environ.get("NVD_INTEGRATION_EXPERIMENTAL") == "1"


def local_sample_row(sample: dict[str, Any], local_fastq_dir: Path) -> dict[str, str]:
    fastq1 = sample.get("fastq1")
    fastq2 = sample.get("fastq2")
    fastq1_glob = sample.get("fastq1_glob")
    fastq2_glob = sample.get("fastq2_glob")
    return {
        "sample_id": str(sample["sample_id"]),
        "srr": "",
        "platform": "illumina",
        "fastq1": str(fastq1.resolve()) if isinstance(fastq1, Path) else "",
        "fastq2": str(fastq2.resolve()) if isinstance(fastq2, Path) else "",
        "fastq1_glob": str(local_fastq_dir / fastq1_glob)
        if isinstance(fastq1_glob, str)
        else "",
        "fastq2_glob": str(local_fastq_dir / fastq2_glob)
        if isinstance(fastq2_glob, str)
        else "",
    }


def write_augmented_samplesheet(run_dir: Path) -> Path:
    """Write the e2e samplesheet with portable SRA rows plus absolute local FASTQs."""
    local_fastq_dir = run_dir / "local_fastqs"
    local_fastq_dir.mkdir(parents=True, exist_ok=True)

    glob_links = {
        "local_hits_glob_L001_R1_001.fastq.gz": LOCAL_FASTQ_FIXTURES["hits_r1"],
        "local_hits_glob_L001_R2_001.fastq.gz": LOCAL_FASTQ_FIXTURES["hits_r2"],
        "local_hits_glob_L002_R1_001.fastq.gz": LOCAL_FASTQ_FIXTURES["water_plus_hits_r1"],
        "local_hits_glob_L002_R2_001.fastq.gz": LOCAL_FASTQ_FIXTURES["water_plus_hits_r2"],
    }
    for name, target in glob_links.items():
        (local_fastq_dir / name).symlink_to(target.resolve())

    fieldnames = [
        "sample_id",
        "srr",
        "platform",
        "fastq1",
        "fastq2",
        "fastq1_glob",
        "fastq2_glob",
    ]
    rows: list[dict[str, str]] = []
    with SAMPLESHEET.open(newline="", encoding="utf-8") as handle:
        rows.extend(
            {column: row.get(column, "") for column in fieldnames}
            for row in csv.DictReader(handle)
        )

    rows.extend(local_sample_row(sample, local_fastq_dir) for sample in LOCAL_E2E_SAMPLES)

    samplesheet = run_dir / "integration_samplesheet.csv"
    with samplesheet.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return samplesheet


def run_nextflow() -> tuple[subprocess.CompletedProcess[str], Path]:
    profile = os.environ.get("NVD_INTEGRATION_PROFILE", "test")
    show_progress = os.environ.get("NVD_E2E_SHOW_PROGRESS") == "1"
    experimental = integration_experimental_enabled()
    run_dir = make_e2e_run_dir()
    results_dir = run_dir / "results"
    work_dir = run_dir / "work"
    taxonomy_dir = run_dir / "taxonomy"
    print(f"NVD e2e run directory: {run_dir}", flush=True)
    print(f"NVD e2e results directory: {results_dir}", flush=True)
    print(f"NVD e2e work directory: {work_dir}", flush=True)
    samplesheet = write_augmented_samplesheet(run_dir)
    write_mini_taxdump(taxonomy_dir)
    command = [
        "nextflow",
        "run",
        ".",
        "-profile",
        profile,
        "--samplesheet",
        str(samplesheet),
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
    if experimental:
        command.extend(
            [
                "--experimental",
                "true",
                "--sourmash_ref_fasta",
                str(SOURMASH_REF_FASTA),
                "--sourmash_lineages_path",
                str(SOURMASH_LINEAGES),
            ],
        )
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
    experimental = integration_experimental_enabled()

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

    experiment_blast = (
        results_root / "08_experiment_summary" / "experiment_blast_results.tsv"
    )
    assert experiment_blast.is_file(), (
        f"Missing experiment BLAST summary: {experiment_blast}"
    )

    final_text = "\n".join(
        path.read_text(encoding="utf-8") for path in final_blast_files
    )
    experiment_rows = read_tsv_rows(experiment_blast)
    assert experiment_rows, f"No experiment BLAST rows found in {experiment_blast}"

    for run_info in manifest["sra_runs"]:
        sample_rows = [
            row for row in experiment_rows if row.get("sample") == run_info["sample_id"]
        ]
        assert sample_rows, (
            f"No experiment BLAST rows found for {run_info['sample_id']}"
        )
        assert any(
            row.get("staxids") == str(run_info["taxid"]) for row in sample_rows
        ), f"No {run_info['taxid']} BLAST taxid found for {run_info['sample_id']}"
        assert any(
            run_info["expected_organism"] in row.get("rank", "") for row in sample_rows
        ), (
            f"No {run_info['expected_organism']} lineage found for "
            f"{run_info['sample_id']}"
        )
        expected_tasks = run_info.get("expected_tasks", [])
        observed_tasks = {row.get("task") for row in sample_rows}
        assert set(expected_tasks) <= observed_tasks, (
            f"Missing expected BLAST tasks for {run_info['sample_id']}: "
            f"expected {sorted(expected_tasks)}, observed {sorted(observed_tasks)}"
        )

    resolved_manifest = results_root / "00_input_resolution" / "resolved_reads.jsonl"
    resolved_records = [
        json.loads(line)
        for line in resolved_manifest.read_text(encoding="utf-8").splitlines()
    ]
    resolved_by_sample = {record["sample_id"]: record for record in resolved_records}
    for run_info in LOCAL_E2E_SAMPLES:
        record = resolved_by_sample[run_info["sample_id"]]
        assert record["source"] == run_info["source"]

    glob_record = resolved_by_sample["local_hits_glob"]
    assert [Path(path).name for path in glob_record["r1"]] == [
        "local_hits_glob_L001_R1_001.fastq.gz",
        "local_hits_glob_L002_R1_001.fastq.gz",
    ]
    assert [Path(path).name for path in glob_record["r2"]] == [
        "local_hits_glob_L001_R2_001.fastq.gz",
        "local_hits_glob_L002_R2_001.fastq.gz",
    ]

    assert "Orf virus" in final_text
    assert "Monkeypox virus" in final_text

    if experimental:
        sourmash_root = results_root / "experimental_sourmash"
        ref_dir = sourmash_root / "reference_profiling" / "reference"
        gather_dir = sourmash_root / "reference_profiling" / "gather"
        merged_taxburst_dir = sourmash_root / "reference_profiling" / "reports" / "taxburst"
        taxburst_dir = merged_taxburst_dir / "per_sample"
        sankey_dir = sourmash_root / "reference_profiling" / "reports" / "sankey"

        ref_sketches = sorted(ref_dir.glob("sourmash_reference.k31.scaled50.sig.zip"))
        assert ref_sketches, f"Missing sourmash reference sketch in {ref_dir}"

        expected_species_by_sample = {
            run_info["sample_id"]: run_info["expected_organism"]
            for run_info in manifest["sra_runs"]
        }
        merged_taxburst_html = merged_taxburst_dir / "sourmash.taxburst.html"
        assert merged_taxburst_html.is_file(), (
            f"Missing merged sourmash taxburst report: {merged_taxburst_html}"
        )
        assert merged_taxburst_html.stat().st_size > 0, (
            f"Empty merged sourmash taxburst report: {merged_taxburst_html}"
        )

        for sample_id, expected_species in expected_species_by_sample.items():
            gather_csv = gather_dir / f"{sample_id}.sourmash.gather.csv"
            assert gather_csv.is_file(), f"Missing sourmash gather CSV: {gather_csv}"
            gather_rows = read_csv_rows(gather_csv)
            assert gather_rows, f"No sourmash gather rows found for {sample_id}"
            assert any(
                expected_species in (row.get("name") or row.get("match_name", ""))
                for row in gather_rows
            ), f"No {expected_species} sourmash gather hit found for {sample_id}"

            taxburst_html = taxburst_dir / f"{sample_id}.sourmash.taxburst.html"
            taxburst_json = taxburst_dir / f"{sample_id}.sourmash.taxburst.json"
            for report in (taxburst_html, taxburst_json):
                assert report.is_file(), f"Missing sourmash taxburst report: {report}"
                assert report.stat().st_size > 0, (
                    f"Empty sourmash taxburst report: {report}"
                )

            sankey_html = sankey_dir / f"{sample_id}.sourmash.sankey.html"
            assert sankey_html.is_file(), (
                f"Missing sourmash Sankey report: {sankey_html}"
            )
            assert sankey_html.stat().st_size > 0, (
                f"Empty sourmash Sankey report: {sankey_html}"
            )
