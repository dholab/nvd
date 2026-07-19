"""Slow end-to-end tests for the checked-in mini viral fixtures."""

from __future__ import annotations

import csv
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
10239\t|\t1\t|\tacellular root\t|\t
2732408\t|\t10239\t|\tphylum\t|\t
2732506\t|\t2732408\t|\tclass\t|\t
2732544\t|\t2732506\t|\torder\t|\t
10240\t|\t2732544\t|\tfamily\t|\t
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
2732408\t|\tNucleocytoviricota\t|\t\t|\tscientific name\t|
2732506\t|\tPokkesviricetes\t|\t\t|\tscientific name\t|
2732544\t|\tChitovirales\t|\t\t|\tscientific name\t|
10240\t|\tPoxviridae\t|\t\t|\tscientific name\t|
10242\t|\tOrthopoxvirus\t|\t\t|\tscientific name\t|
10244\t|\tMonkeypox virus\t|\t\t|\tscientific name\t|
10255\t|\tParapoxvirus\t|\t\t|\tscientific name\t|
10258\t|\tOrf virus\t|\t\t|\tscientific name\t|
""",
        encoding="utf-8",
    )
    (taxonomy_dir / "merged.dmp").write_text("", encoding="utf-8")


def load_manifest() -> dict[str, Any]:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def verify_fixture_files(manifest: dict[str, Any]) -> None:
    missing: list[str] = []
    for filename in manifest["files"]:
        path = DATA / filename
        if not path.exists():
            missing.append(filename)

    assert not missing, "Missing integration fixture files: " + ", ".join(missing)


def selected_manifest_sra_runs(
    manifest: dict[str, Any],
    *,
    samplesheet: Path = SAMPLESHEET,
) -> list[dict[str, Any]]:
    """Return manifest metadata only for SRA rows selected by the samplesheet."""
    manifest_by_sample = {
        str(run_info["sample_id"]): run_info for run_info in manifest["sra_runs"]
    }

    with samplesheet.open(newline="", encoding="utf-8") as handle:
        selected_sample_ids = [
            row.get("sample_id", "") for row in csv.DictReader(handle)
        ]

    return [
        manifest_by_sample[sample_id]
        for sample_id in selected_sample_ids
        if sample_id in manifest_by_sample
    ]


def test_selected_manifest_sra_runs_follow_samplesheet_rows(tmp_path: Path) -> None:
    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text(
        "sample_id,srr,platform,fastq1,fastq2\n"
        "monkeypox_pt1020_2026,ERR17356125,illumina,,\n",
        encoding="utf-8",
    )

    selected = selected_manifest_sra_runs(load_manifest(), samplesheet=samplesheet)

    assert [run_info["sample_id"] for run_info in selected] == [
        "monkeypox_pt1020_2026",
    ]


def test_mini_sourmash_lineages_support_bioboxes() -> None:
    """The mini sourmash taxonomy fixture must include BioBoxes taxpaths."""
    with SOURMASH_LINEAGES.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    ranks = (
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
    )
    assert rows, f"No sourmash lineage rows found in {SOURMASH_LINEAGES}"
    for row in rows:
        taxpath = row.get("taxpath", "")
        taxids = taxpath.split("|") if taxpath else []
        assert len(taxids) == len(ranks), row
        for rank, taxid in zip(ranks, taxids, strict=True):
            assert bool(row.get(rank)) == bool(taxid), row


def run_sourmash(*args: str | Path) -> subprocess.CompletedProcess[str]:
    """Run sourmash from the test environment."""
    return subprocess.run(  # noqa: S603
        ["sourmash", *(str(arg) for arg in args)],  # noqa: S607
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=True,
        timeout=120,
    )


def test_sourmash_tax_metagenome_writes_all_formats_with_strain_taxids(
    tmp_path: Path,
) -> None:
    """Sourmash 4.9.4 consumes the complete positional taxonomy contract."""
    query = tmp_path / "query.sig.zip"
    reference = tmp_path / "reference.sig.zip"
    gather = tmp_path / "gather.csv"
    output_base = tmp_path / "profile"
    commands = [
        (
            "sketch",
            "dna",
            SOURMASH_REF_FASTA,
            "-p",
            "dna,k=31,scaled=50,abund",
            "-o",
            query,
        ),
        (
            "sketch",
            "dna",
            SOURMASH_REF_FASTA,
            "--singleton",
            "-p",
            "dna,k=31,scaled=50",
            "-o",
            reference,
        ),
        (
            "gather",
            query,
            reference,
            "-k",
            "31",
            "--scaled",
            "50",
            "-o",
            gather,
        ),
        (
            "tax",
            "metagenome",
            "--gather-csv",
            gather,
            "--taxonomy-csv",
            SOURMASH_LINEAGES,
            "--keep-identifier-versions",
            "--use-abundances",
            "--output-format",
            "csv_summary",
            "lineage_summary",
            "krona",
            "kreport",
            "bioboxes",
            "--rank",
            "species",
            "--output-base",
            output_base,
        ),
    ]

    for command in commands:
        result = run_sourmash(*command)
        assert result.returncode == 0, result.stderr

    outputs = (
        tmp_path / "profile.summarized.csv",
        tmp_path / "profile.lineage_summary.tsv",
        tmp_path / "profile.krona.tsv",
        tmp_path / "profile.kreport.txt",
        tmp_path / "profile.bioboxes.profile",
    )
    for output in outputs:
        assert output.stat().st_size > 0, output

    bioboxes_lines = outputs[-1].read_text(encoding="utf-8").splitlines()
    assert not any("None" in line for line in bioboxes_lines)
    profile_rows = [
        line.split("\t")
        for line in bioboxes_lines
        if line and not line.startswith(("#", "@"))
    ]
    assert profile_rows
    rank_order = (
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "strain",
    )
    for taxid, rank, taxpath, taxpath_names, percentage in profile_rows:
        taxids = taxpath.split("|")
        names = taxpath_names.split("|")
        assert taxid == taxids[-1]
        assert rank == rank_order[len(taxids) - 1]
        assert len(taxids) == len(names)
        assert all(taxids)
        assert 0 <= float(percentage) <= 100  # noqa: PLR2004


def test_mini_nvd_taxdump_preserves_noncanonical_virus_root_rank(
    tmp_path: Path,
) -> None:
    """The mini NVD taxonomy should exercise modern NCBI virus ranks."""
    write_mini_taxdump(tmp_path)

    nodes_text = (tmp_path / "nodes.dmp").read_text(encoding="utf-8")

    assert "10239\t|\t1\t|\tacellular root\t|" in nodes_text
    assert "10239\t|\t1\t|\tsuperkingdom\t|" not in nodes_text


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


def integration_skip_assembly_enabled() -> bool:
    return os.environ.get("NVD_INTEGRATION_SKIP_ASSEMBLY") == "1"


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
        "local_hits_glob_L002_R1_001.fastq.gz": LOCAL_FASTQ_FIXTURES[
            "water_plus_hits_r1"
        ],
        "local_hits_glob_L002_R2_001.fastq.gz": LOCAL_FASTQ_FIXTURES[
            "water_plus_hits_r2"
        ],
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

    rows.extend(
        local_sample_row(sample, local_fastq_dir) for sample in LOCAL_E2E_SAMPLES
    )

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
    skip_assembly = integration_skip_assembly_enabled()
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
    if skip_assembly:
        command.extend(["--skip_assembly", "true"])
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
    verify_fixture_files(manifest)
    selected_sra_runs = selected_manifest_sra_runs(manifest)
    experimental = integration_experimental_enabled()
    skip_assembly = integration_skip_assembly_enabled()

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

    experiment_blast = (
        results_root / "08_experiment_summary" / "experiment_blast_results.tsv"
    )

    if skip_assembly:
        assert not final_blast_files, (
            f"Skip-assembly run unexpectedly produced final BLAST TSVs: {final_blast_files}"
        )
        if experiment_blast.exists():
            assert experiment_blast.read_text(encoding="utf-8") == ""
    else:
        assert experiment_blast.is_file(), (
            f"Missing experiment BLAST summary: {experiment_blast}"
        )
        assert final_blast_files, f"No final BLAST TSVs found in {final_dir}"

        final_text = "\n".join(
            path.read_text(encoding="utf-8") for path in final_blast_files
        )
        experiment_rows = read_tsv_rows(experiment_blast)
        assert experiment_rows, f"No experiment BLAST rows found in {experiment_blast}"

        # Assert per-sample biological expectations for the rows actually under
        # test. Coupling this loop to every manifest row makes a deleted
        # samplesheet row fail as a missing output, even though the pipeline did
        # exactly what the samplesheet requested.
        for run_info in selected_sra_runs:
            sample_rows = [
                row
                for row in experiment_rows
                if row.get("sample") == run_info["sample_id"]
            ]
            assert sample_rows, (
                f"No experiment BLAST rows found for {run_info['sample_id']}"
            )
            assert any(
                row.get("staxids") == str(run_info["taxid"]) for row in sample_rows
            ), f"No {run_info['taxid']} BLAST taxid found for {run_info['sample_id']}"
            assert any(
                run_info["expected_organism"] in row.get("rank", "")
                for row in sample_rows
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

        resolved_manifest = (
            results_root / "00_input_resolution" / "resolved_reads.jsonl"
        )
        resolved_records = [
            json.loads(line)
            for line in resolved_manifest.read_text(encoding="utf-8").splitlines()
        ]
        resolved_by_sample = {
            record["sample_id"]: record for record in resolved_records
        }
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

        for organism in {
            str(run_info["expected_organism"]) for run_info in selected_sra_runs
        }:
            assert organism in final_text

    if experimental:
        rapid_screening_root = results_root / "experimental_rapid_screening"
        sourmash_root = rapid_screening_root / "engines" / "sourmash"
        ref_dir = sourmash_root / "reference_profiling" / "reference"
        gather_dir = sourmash_root / "reference_profiling" / "gather"
        merged_taxburst_dir = (
            sourmash_root / "reference_profiling" / "reports" / "taxburst"
        )
        taxburst_dir = merged_taxburst_dir / "per_sample"
        sankey_dir = sourmash_root / "reference_profiling" / "reports" / "sankey"

        ref_sketches = sorted(ref_dir.glob("sourmash_reference.k31.scaled50.sig.zip"))
        assert ref_sketches, f"Missing sourmash reference sketch in {ref_dir}"

        expected_species_by_sample = {
            run_info["sample_id"]: run_info["expected_organism"]
            for run_info in selected_sra_runs
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

        eval_root = rapid_screening_root / "eval"
        for eval_artifact in (
            eval_root / "database" / "rapid_screening_eval.duckdb",
            eval_root / "exports" / "screening_signal_followup_by_sample_rank.tsv",
            eval_root / "exports" / "screening_signals_without_same_rank_followup.tsv",
            eval_root / "reports" / "rapid_screening_eval.html",
        ):
            assert eval_artifact.is_file(), (
                f"Missing rapid-screening eval artifact: {eval_artifact}"
            )
            assert eval_artifact.stat().st_size > 0, (
                f"Empty rapid-screening eval artifact: {eval_artifact}"
            )
