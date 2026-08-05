"""Slow end-to-end tests for the checked-in mini viral fixtures."""

from __future__ import annotations

import csv
import json
import os
import secrets
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


def read_tsv_document(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None, f"Missing TSV header: {path}"
        return reader.fieldnames, list(reader)


def assert_concatenated_tsv(output: Path, inputs: list[Path]) -> None:
    ordered_inputs = sorted(inputs, key=lambda path: (path.name, str(path)))
    if not ordered_inputs:
        assert output.read_text(encoding="utf-8") == ""
        return

    expected_header: list[str] | None = None
    expected_rows: list[dict[str, str]] = []
    for path in ordered_inputs:
        header, rows = read_tsv_document(path)
        if expected_header is None:
            expected_header = header
        else:
            assert header == expected_header, f"Unexpected per-sample header: {path}"
        expected_rows.extend(rows)

    output_header, output_rows = read_tsv_document(output)
    assert output_header == expected_header
    assert output_rows == expected_rows


def assert_long_read_assembly_outputs(
    results_root: Path,
    selected_sra_runs: list[dict[str, Any]],
) -> None:
    assembly_root = results_root / "03_assembled_contigs"
    eligibility_report = (
        assembly_root / "decisions" / "long_read_assembly_eligibility.tsv"
    )
    assert eligibility_report.is_file(), (
        f"Missing long-read assembly eligibility report: {eligibility_report}"
    )
    eligibility_rows = read_tsv_rows(eligibility_report)
    expected_long_read_samples = {
        run_info["sample_id"]
        for run_info in selected_sra_runs
        if run_info["platform"] != "illumina"
    }
    assert {row["sample_id"] for row in eligibility_rows} == (
        expected_long_read_samples
    )
    for row in eligibility_rows:
        for assembler in ("metamdbg", "myloasm", "metaflye"):
            decision = row[f"{assembler}_decision"]
            qualifying_reads = int(row[f"{assembler}_qualifying_reads"])
            assert decision in {"run", "skip"}
            assert (decision == "run") == (qualifying_reads > 0)

    assert not (assembly_root / "long_read_assemblers").exists()
    for assembler in ("metamdbg", "myloasm", "metaflye"):
        assembler_dir = assembly_root / assembler
        assert assembler_dir.is_dir(), (
            f"Missing {assembler} output directory: {assembler_dir}"
        )


def assert_read_profiles_respect_length_filter(results_root: Path) -> None:
    """Assert final read profiles satisfy their declared minimum length."""
    profile_dir = results_root / "02_preprocessed_reads" / "profiles"
    profiles = sorted(profile_dir.glob("*.fastx_profile.json"))
    assert profiles, f"No final read profiles found in {profile_dir}"

    for path in profiles:
        profile = json.loads(path.read_text(encoding="utf-8"))
        thresholds = {
            threshold["name"]: float(threshold["value"])
            for threshold in profile["thresholds"]
        }
        minimum = profile["length"]["min"]
        expected_minimum = thresholds["min_read_length"]
        assert minimum is not None, f"{path} does not report a minimum read length"
        assert float(minimum) >= expected_minimum, (
            f"{path} reports minimum read length {minimum}, below the declared "
            f"filter threshold {expected_minimum}"
        )


def assert_target_enrichment_outputs(
    results_root: Path,
    expected_sample_ids: set[str],
) -> None:
    """Assert per-sample and experiment-wide enrichment results are published."""
    enriched_reads_dir = results_root / "01_target_enrichment" / "reads"
    for sample_id in expected_sample_ids:
        enriched_reads = enriched_reads_dir / f"{sample_id}.target_enriched.fastq.gz"
        assert enriched_reads.is_file(), (
            f"Missing target enrichment result: {enriched_reads}"
        )
        assert enriched_reads.stat().st_size > 0, (
            f"Empty target enrichment result: {enriched_reads}"
        )

    summary_dir = results_root / "01_target_enrichment" / "summaries"
    for sample_id in expected_sample_ids:
        summary = summary_dir / f"{sample_id}.deacon_filter.json"
        assert summary.is_file(), f"Missing target enrichment summary: {summary}"
        assert summary.stat().st_size > 0, f"Empty target enrichment summary: {summary}"
    experiment_summary = summary_dir / "target_enrichment_summary.tsv"
    assert experiment_summary.is_file(), (
        f"Missing target enrichment summary: {experiment_summary}"
    )
    assert experiment_summary.stat().st_size > 0, (
        f"Empty target enrichment summary: {experiment_summary}"
    )

    plot_dir = results_root / "01_target_enrichment" / "plots"
    expected_plots = (
        "target_enriched_bases_ranked.html",
        "target_retained_vs_filtered_stacked.html",
        "target_reads_vs_bases_scatter.html",
    )
    for filename in expected_plots:
        plot = plot_dir / filename
        assert plot.is_file(), f"Missing target enrichment plot: {plot}"
        assert plot.stat().st_size > 0, f"Empty target enrichment plot: {plot}"


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
        "--filter_reads",
        "true",
    ]
    if experimental:
        command.extend(
            [
                "--experimental",
                "true",
                "--merge_pairs",
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


def write_local_only_samplesheet(run_dir: Path) -> Path:
    """Write a samplesheet with only the local FASTQ fixtures (no SRA rows).

    Dropping the SRA rows keeps the LIMS end-to-end test off the network so it
    can be marked ``slow`` alone. Mirrors the local portion of
    ``write_augmented_samplesheet``.
    """
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
    rows = [local_sample_row(sample, local_fastq_dir) for sample in LOCAL_E2E_SAMPLES]

    samplesheet = run_dir / "integration_samplesheet.csv"
    with samplesheet.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return samplesheet


LABKEY_SECRET_NAME = "LABKEY_API_KEY"


def get_labkey_secret(env: dict[str, str]) -> str | None:
    """Return the current ``LABKEY_API_KEY`` secret value, or ``None`` if unset.

    ``nextflow secrets get`` prints the plaintext value, or the literal ``null``
    when the secret does not exist.
    """
    result = subprocess.run(  # noqa: S603
        ["nextflow", "secrets", "get", LABKEY_SECRET_NAME],  # noqa: S607
        cwd=ROOT,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )
    value = (result.stdout or "").strip()
    if result.returncode != 0 or value in {"", "null"}:
        return None
    return value


def set_labkey_secret(env: dict[str, str], value: str = "mock-api-key") -> None:
    """Register a LabKey API-key secret so the LIMS processes can launch.

    The LabKey processes declare ``secret 'LABKEY_API_KEY'``; Nextflow refuses
    to launch them unless the secret exists. The dummy value is never validated
    by the mock endpoint, so any placeholder works. Stored in the default
    Nextflow secrets store (``$NXF_HOME``).
    """
    subprocess.run(  # noqa: S603
        ["nextflow", "secrets", "set", LABKEY_SECRET_NAME, value],  # noqa: S607
        cwd=ROOT,
        env=env,
        check=True,
        text=True,
        capture_output=True,
    )


def delete_labkey_secret(env: dict[str, str]) -> None:
    subprocess.run(  # noqa: S603
        ["nextflow", "secrets", "delete", LABKEY_SECRET_NAME],  # noqa: S607
        cwd=ROOT,
        env=env,
        check=False,
        text=True,
        capture_output=True,
    )


def restore_labkey_secret(env: dict[str, str], previous: str | None) -> None:
    """Restore the pre-existing secret, or delete the dummy if there was none.

    Kept exception-safe so a restore hiccup never masks a test failure raised in
    the surrounding ``finally``.
    """
    try:
        if previous is not None:
            set_labkey_secret(env, previous)
        else:
            delete_labkey_secret(env)
    except Exception:  # noqa: BLE001, S110 - never mask the test's own failure
        pass


def run_labkey_nextflow(  # noqa: PLR0913
    *,
    server: str,
    webdav: str,
    project: str,
    schema: str,
    hits_list: str,
    fasta_list: str,
    experiment_id: int,
    experimental: bool = True,
    cert_file: Path | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    """Run the local-only pipeline with LabKey enabled against the given endpoint.

    Shared by the mocked and real-LabKey e2e tests. ``experimental`` (default on)
    turns on the read-derived query classes; set it off for a contig-only run.
    ``cert_file`` is set only for the mock (its self-signed TLS cert); a real
    LabKey server presents a publicly-trusted cert and needs no override.
    """
    profile = os.environ.get("NVD_INTEGRATION_PROFILE", "test")
    show_progress = os.environ.get("NVD_E2E_SHOW_PROGRESS") == "1"
    run_dir = make_e2e_run_dir()
    results_dir = run_dir / "results"
    work_dir = run_dir / "work"
    taxonomy_dir = run_dir / "taxonomy"
    print(f"NVD LIMS e2e run directory: {run_dir}", flush=True)
    samplesheet = write_local_only_samplesheet(run_dir)
    write_mini_taxdump(taxonomy_dir)

    env = os.environ.copy()
    if cert_file is not None:
        # Make both the requests-based LabKey client and the urllib WebDAV client
        # trust the mock's self-signed TLS certificate.
        env["SSL_CERT_FILE"] = str(cert_file)
        env["REQUESTS_CA_BUNDLE"] = str(cert_file)

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
        str(experiment_id),
        "--results",
        str(results_dir),
        "--work_dir",
        str(work_dir),
        "--filter_reads",
        "true",
        # Disable deacon target enrichment. The mini deacon index shares no
        # k-mers with the local fixtures, so enrichment would drop every read
        # (0/N enriched) and nothing would reach BLAST. With enrichment off all
        # reads are retained, assemble into contigs, and hit the mini BLAST db,
        # producing the per-batch hits the LIMS upload consumes -- fully offline.
        "--no_enrichment",
        "true",
        "--labkey",
        "true",
        "--labkey_server",
        server,
        "--labkey_webdav",
        webdav,
        "--labkey_project_name",
        project,
        "--labkey_schema",
        schema,
        "--labkey_blast_meta_hits_list",
        hits_list,
        "--labkey_blast_fasta_list",
        fasta_list,
    ]
    if experimental:
        # Read-derived query classes (overlap_merged_pair, single_read) alongside
        # short_assembly_contig, so the run exercises the full per-query-class
        # split. Needs the sourmash rapid-screening fixtures experimental turns on.
        command += [
            "--experimental",
            "true",
            "--merge_pairs",
            "true",
            "--sourmash_ref_fasta",
            str(SOURMASH_REF_FASTA),
            "--sourmash_lineages_path",
            str(SOURMASH_LINEAGES),
        ]
    completed = subprocess.run(  # noqa: S603
        command,
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=not show_progress,
        env=env,
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


def data_hits_inserts(mock: Any) -> list[dict[str, Any]]:
    """Eager per-batch BLAST hit inserts (rows carrying a ``query_class``).

    Distinguishes the real per-(sample, query_class) hit uploads from the
    schema-probe row ``validate_labkey.py`` inserts-then-deletes, which has no
    ``query_class`` column.
    """
    return [
        op
        for op in mock.inserts
        if op["table"] == "hits"
        and any(row.get("query_class") for row in op["rows"])
    ]


# Marked `labkey` (with the real test) so it runs only via `just e2e-labkey` and
# is excluded from `just e2e`/CI. Even though it uses a mock endpoint and a dummy
# secret, CI has no LabKey credentials/secret store, so LIMS e2e stays opt-in.
@pytest.mark.slow
@pytest.mark.network
@pytest.mark.labkey
def test_lims_enabled_pipeline_uploads_eagerly_and_dedups() -> None:
    """LabKey-enabled run uploads hits per batch; a rerun skips present combos."""
    from tests.scripts.mock_labkey_server import mock_labkey_server

    with mock_labkey_server() as mock:
        secret_env = os.environ.copy()
        # Capture any developer-set secret so we restore (never destroy) it.
        previous_secret = get_labkey_secret(secret_env)
        set_labkey_secret(secret_env)
        mock_coords = {
            "server": mock.base_url,
            "webdav": mock.webdav_url,
            "project": "test",
            "schema": "lists",
            "hits_list": "hits",
            "fasta_list": "fasta",
        }
        try:
            first, run_dir = run_labkey_nextflow(
                **mock_coords, experiment_id=1, cert_file=mock.cert_file
            )
            assert first.returncode == 0, (
                f"LabKey e2e run failed.\n\nRun directory: {run_dir}\n\nSTDOUT:\n"
                + (first.stdout or "<not captured; see nextflow logs in run dir>")
                + "\n\nSTDERR:\n"
                + (first.stderr or "<not captured; see nextflow logs in run dir>")
            )

            results_root = run_dir / "results" / "nvd"
            assert (results_root / "13_labkey_uploads").exists(), (
                f"Missing LabKey uploads results dir under {results_root}"
            )

            first_hits = data_hits_inserts(mock)
            assert first_hits, "no eager per-batch hits insert carrying query_class"
            for op in first_hits:
                for row in op["rows"]:
                    assert row.get("query_class"), f"hits row missing query_class: {row}"

            observed_classes = {
                row.get("query_class") for op in first_hits for row in op["rows"]
            }
            print(f"query_class values uploaded: {sorted(observed_classes)}", flush=True)
            assert "short_assembly_contig" in observed_classes, (
                f"expected the contig query class; observed {sorted(observed_classes)}"
            )
            assert observed_classes & {"overlap_merged_pair", "single_read"}, (
                "expected read-derived query classes too (experimental mode splits "
                f"queries by read type); observed {sorted(observed_classes)}"
            )

            before = len(first_hits)
            second, second_dir = run_labkey_nextflow(
                **mock_coords, experiment_id=1, cert_file=mock.cert_file
            )
            assert second.returncode == 0, (
                f"Second LabKey e2e run failed.\n\nRun directory: {second_dir}\n\n"
                "STDOUT:\n"
                + (second.stdout or "<not captured; see nextflow logs in run dir>")
                + "\n\nSTDERR:\n"
                + (second.stderr or "<not captured; see nextflow logs in run dir>")
            )
            assert len(data_hits_inserts(mock)) == before, (
                "second run must skip already-present combos (dedup on the "
                "destination hits list), but new hits inserts were recorded"
            )
        finally:
            restore_labkey_secret(secret_env, previous_secret)


# Opt-in real-LabKey e2e target. The LabKey coordinates come from an `nvd`
# preset (a validated bundle of pipeline params); the API key stays a Nextflow
# secret (a params preset must not hold secrets). The preset's existence is the
# gate: absent -> skip; present-but-incomplete or missing secret -> fail loudly.
DEFAULT_LABKEY_PRESET = "labkey-e2e"
# local field -> the pipeline param the preset stores it under.
LABKEY_PRESET_PARAMS = {
    "server": "labkey_server",
    "project": "labkey_project_name",
    "schema": "labkey_schema",
    "hits_list": "labkey_blast_meta_hits_list",
    "fasta_list": "labkey_blast_fasta_list",
    "webdav": "labkey_webdav",
}


def labkey_test_preset_name() -> str:
    """Preset name for the real-LabKey e2e: ``NVD_TEST_LABKEY_PRESET`` or default."""
    return os.environ.get("NVD_TEST_LABKEY_PRESET", DEFAULT_LABKEY_PRESET)


def load_labkey_test_target() -> dict[str, str] | None:
    """Resolve the real-LabKey coordinates from an ``nvd`` preset plus the secret.

    Returns ``None`` when the preset does not exist (the test skips). When the
    preset exists it must carry every LabKey coordinate param, and the
    ``LABKEY_API_KEY`` Nextflow secret must be set; otherwise the test fails
    naming exactly what is missing. The returned dict adds ``api_key`` (read from
    the secret store) to the preset's coordinates.
    """
    from py_nvd.presets import get_preset_store

    name = labkey_test_preset_name()
    preset = get_preset_store().get(name)
    if preset is None:
        return None
    params = preset.params
    missing = [
        param
        for param in LABKEY_PRESET_PARAMS.values()
        if not str(params.get(param, "")).strip()
    ]
    if missing:
        pytest.fail(
            f"preset {name!r} is missing LabKey coordinates: {', '.join(missing)} "
            "(register them with `nvd preset register`)",
        )
    api_key = get_labkey_secret(os.environ.copy())
    if not api_key:
        pytest.fail(
            f"preset {name!r} found, but the LABKEY_API_KEY secret is not set; "
            'set it with `nvd secrets set LABKEY_API_KEY "<key>"` to run',
        )
    target = {
        field: str(params[param]).strip()
        for field, param in LABKEY_PRESET_PARAMS.items()
    }
    target["api_key"] = api_key
    return target


def real_labkey_api(cfg: dict[str, str]) -> Any:
    """Build a LabKey APIWrapper for verifying and cleaning up the real list."""
    from labkey.api_wrapper import APIWrapper

    return APIWrapper(cfg["server"], cfg["project"], api_key=cfg["api_key"])


def experiment_rows(
    api: Any,
    schema: str,
    table: str,
    experiment_id: int,
) -> list[dict[str, Any]]:
    """Return every row on ``table`` belonging to ``experiment_id``."""
    from labkey.query import QueryFilter

    result = api.query.select_rows(
        schema_name=schema,
        query_name=table,
        filter_array=[QueryFilter("experiment", experiment_id, "eq")],
    )
    return result.get("rows", []) if result else []


def delete_experiment_rows(
    api: Any,
    schema: str,
    table: str,
    experiment_id: int,
) -> None:
    """Delete every row on ``table`` for ``experiment_id`` (by its list Key)."""
    rows = experiment_rows(api, schema, table, experiment_id)
    keys = [{"Key": row["Key"]} for row in rows if "Key" in row]
    if keys:
        api.query.delete_rows(schema_name=schema, query_name=table, rows=keys)


def delete_webdav_experiment(webdav: str, api_key: str, experiment_id: int) -> None:
    """Recursively delete the experiment's WebDAV upload directory.

    The upload processes write to ``{webdav}/{experiment_id}/...`` (LabKey WebDAV
    Basic auth is ``apikey:<key>``). A collection DELETE removes it and every file
    under it, so a run leaves no WebDAV residue. A 404 (nothing was uploaded) is
    treated as already-clean.
    """
    import base64  # noqa: PLC0415
    import urllib.error  # noqa: PLC0415
    import urllib.request  # noqa: PLC0415

    url = webdav.rstrip("/") + f"/{experiment_id}/"
    token = base64.b64encode(f"apikey:{api_key}".encode()).decode()
    request = urllib.request.Request(  # noqa: S310 - fixed LabKey WebDAV https URL
        url,
        method="DELETE",
        headers={"Authorization": f"Basic {token}"},
    )
    try:
        urllib.request.urlopen(request, timeout=60)  # noqa: S310
    except urllib.error.HTTPError as exc:
        if exc.code != 404:
            raise


# Marked `labkey` so it only runs via the dedicated `just e2e-labkey` target and
# is excluded from `just e2e`/CI (`-m '... and not labkey'`). Otherwise, once the
# nvd preset is registered, a plain `just e2e` would hit the real LabKey server.
@pytest.mark.slow
@pytest.mark.network
@pytest.mark.labkey
def test_lims_enabled_real_labkey_uploads_and_dedups() -> None:
    """Against a real LabKey list: first run inserts per batch; a rerun is deduped.

    Opt-in via an `nvd` preset (coordinates) plus the LABKEY_API_KEY secret; see
    labkey_test_preset_name(). Uses a random experiment_id for isolation and
    removes all of that experiment's rows from both lists at the end, so the test
    lists are left as they were found.
    """
    cfg = load_labkey_test_target()
    if cfg is None:
        pytest.skip(
            f"register an `nvd` preset named {labkey_test_preset_name()!r} with the "
            "LabKey coordinates (and set the LABKEY_API_KEY secret) to run the real "
            "LabKey e2e",
        )

    coords = {
        "server": cfg["server"],
        "webdav": cfg["webdav"],
        "project": cfg["project"],
        "schema": cfg["schema"],
        "hits_list": cfg["hits_list"],
        "fasta_list": cfg["fasta_list"],
    }
    # Pin the experiment id with NVD_TEST_LABKEY_EXPERIMENT_ID to rerun against
    # an id that already has rows (verifies the check-ahead skips duplicates
    # across invocations, not just within one). Otherwise use a random,
    # collision-avoiding id in the int32 range LabKey's integer column accepts.
    pinned = os.environ.get("NVD_TEST_LABKEY_EXPERIMENT_ID")
    experiment_id = int(pinned) if pinned else 2_000_000_000 + secrets.randbelow(147_000_000)

    # Experimental on (default) emits the read-derived query classes too; set
    # NVD_TEST_LABKEY_EXPERIMENTAL=0 for a faster contig-only run.
    experimental = os.environ.get("NVD_TEST_LABKEY_EXPERIMENTAL", "1") != "0"
    print(
        f"NVD real LabKey e2e experiment_id={experiment_id} "
        f"(pinned={bool(pinned)}, experimental={experimental})",
        flush=True,
    )

    # Keep the run's rows and WebDAV files by default so they can be inspected in
    # LabKey. Set NVD_TEST_LABKEY_KEEP=0 to have the test remove everything it
    # created and leave the lists and file store as it found them.
    keep = os.environ.get("NVD_TEST_LABKEY_KEEP", "1") != "0"

    # The API key already lives in the Nextflow secret store (the pipeline reads
    # it directly); we read it only to drive verification/cleanup, and never
    # mutate the store here.
    api = real_labkey_api(cfg)
    try:
        # For a random (ephemeral) id, start clean so the first run genuinely
        # inserts. For a pinned id, keep any existing rows so a rerun exercises
        # the cross-invocation check-ahead dedup instead of re-inserting.
        if not pinned:
            delete_experiment_rows(api, cfg["schema"], cfg["hits_list"], experiment_id)
            delete_experiment_rows(api, cfg["schema"], cfg["fasta_list"], experiment_id)

        first, run_dir = run_labkey_nextflow(
            **coords, experiment_id=experiment_id, experimental=experimental
        )
        assert first.returncode == 0, (
            f"Real LabKey e2e run failed.\n\nRun directory: {run_dir}\n\nSTDOUT:\n"
            + (first.stdout or "<not captured; see nextflow logs in run dir>")
            + "\n\nSTDERR:\n"
            + (first.stderr or "<not captured; see nextflow logs in run dir>")
        )

        hits_after_first = experiment_rows(
            api, cfg["schema"], cfg["hits_list"], experiment_id
        )
        assert hits_after_first, "no hits present for the test experiment after the first run"
        for row in hits_after_first:
            assert row.get("query_class"), f"hits row missing query_class: {row}"
        observed_classes = {row.get("query_class") for row in hits_after_first}
        print(f"query_class values present: {sorted(observed_classes)}", flush=True)
        assert "short_assembly_contig" in observed_classes, (
            f"expected the contig query class; observed {sorted(observed_classes)}"
        )
        if experimental:
            assert observed_classes & {"overlap_merged_pair", "single_read"}, (
                "expected read-derived query classes too (experimental mode splits "
                f"queries by read type); observed {sorted(observed_classes)}"
            )
        first_hits_count = len(hits_after_first)
        first_fasta_count = len(
            experiment_rows(api, cfg["schema"], cfg["fasta_list"], experiment_id),
        )

        # Insert-again is blocked: the rerun dedups against the destination lists.
        second, second_dir = run_labkey_nextflow(
            **coords, experiment_id=experiment_id, experimental=experimental
        )
        assert second.returncode == 0, (
            f"Second real LabKey e2e run failed.\n\nRun directory: {second_dir}\n\n"
            "STDOUT:\n"
            + (second.stdout or "<not captured; see nextflow logs in run dir>")
            + "\n\nSTDERR:\n"
            + (second.stderr or "<not captured; see nextflow logs in run dir>")
        )
        assert (
            len(experiment_rows(api, cfg["schema"], cfg["hits_list"], experiment_id))
            == first_hits_count
        ), "second run must not add hits rows (dedup on the destination hits list)"
        assert (
            len(experiment_rows(api, cfg["schema"], cfg["fasta_list"], experiment_id))
            == first_fasta_count
        ), "second run must not add contig FASTA rows (dedup on the FASTA list)"
    finally:
        # Leave the lists and WebDAV store as we found them (unless KEEP is set).
        # Best-effort so a cleanup hiccup never masks a test failure raised above.
        if keep:
            print(
                f"Leaving experiment {experiment_id} in place for inspection "
                f"(hits list {cfg['hits_list']!r}, fasta list {cfg['fasta_list']!r}, "
                f"WebDAV dir {experiment_id}/). Set NVD_TEST_LABKEY_KEEP=0 to auto-clean.",
                flush=True,
            )
        else:
            try:
                delete_experiment_rows(api, cfg["schema"], cfg["hits_list"], experiment_id)
                delete_experiment_rows(api, cfg["schema"], cfg["fasta_list"], experiment_id)
                delete_webdav_experiment(cfg["webdav"], cfg["api_key"], experiment_id)
            except Exception:  # noqa: BLE001, S110 - cleanup is best-effort; never mask a failure
                pass


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
    assert not (results_root / "13_labkey_uploads").exists()
    expected_sample_ids = {
        str(run_info["sample_id"])
        for run_info in (*LOCAL_E2E_SAMPLES, *selected_sra_runs)
    }
    assert_target_enrichment_outputs(results_root, expected_sample_ids)
    assert_read_profiles_respect_length_filter(results_root)
    merged_blast_dir = results_root / "07_merged_blast_results"
    final_dir = merged_blast_dir / "final"
    final_blast_files = sorted(final_dir.glob("*_blast.final.tsv"))

    experiment_blast = (
        results_root / "12_experiment_summary" / "experiment_blast_results.tsv"
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
        assert_concatenated_tsv(experiment_blast, final_blast_files)

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
        if not skip_assembly:
            assert_long_read_assembly_outputs(results_root, selected_sra_runs)

        big_tables_dir = results_root / "11_big_tables"
        table_inputs = {
            "query_big_table.tsv": big_tables_dir / "query" / "per_sample",
            "taxon_big_table.tsv": big_tables_dir / "taxon" / "per_sample",
        }
        for filename, per_sample_dir in table_inputs.items():
            featured_table = results_root / filename
            grouped_table = big_tables_dir / filename
            assert featured_table.is_file(), (
                f"Missing featured Big Table: {featured_table}"
            )
            assert grouped_table.is_file(), (
                f"Missing grouped Big Table: {grouped_table}"
            )
            assert featured_table.read_bytes() == grouped_table.read_bytes()
            assert_concatenated_tsv(
                grouped_table,
                list(per_sample_dir.glob("*.tsv")),
            )

        sourmash_root = (
            results_root
            / "08_metagenomic_profiles"
            / "rapid_screening"
            / "engines"
            / "sourmash"
        )
        ref_dir = sourmash_root / "reference"
        gather_dir = sourmash_root / "gather"
        merged_taxburst_dir = sourmash_root / "plots" / "taxburst"
        taxburst_dir = merged_taxburst_dir / "per_sample"
        sankey_dir = sourmash_root / "plots" / "sankey"

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

        eval_root = results_root / "10_rapid_screening_eval"
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
