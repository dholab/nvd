"""Slow end-to-end tests for the checked-in mini viral fixtures."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import secrets
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest

from tests.e2e_parity import find_parity_violations, load_parity, summarize_detections
from tests.scripts.build_synthetic_metagenomes import REPRODUCIBLE_CPUS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
MANIFEST = DATA / "reference.manifest.json"
SAMPLESHEET = DATA / "integration_sra_samplesheet.csv"
# Provenance for the committed downsamples that stand in for the SRA rows when
# the e2e runs offline, which is the default. Keyed on the same sample ids as
# SAMPLESHEET, so switching between the two changes only where reads come from.
OFFLINE_SRA_MANIFEST_PATH = DATA / "e2e_sra_fixtures.manifest.json"
DEACON_INDEX = DATA / "mini_virus_deacon.k31w1.idx"
BLAST_DB_PREFIX = "mini_virus_blast"
BLAST_DB = DATA
E2E_OUTPUT_DIR = ROOT / ".e2e"
LOCAL_FASTQ_FIXTURES = {
    "hits_r1": DATA / "hits_only_R1.fastq.gz",
    "hits_r2": DATA / "hits_only_R2.fastq.gz",
    "water_plus_hits_r1": DATA / "water_plus_hits_R1.fastq.gz",
    "water_plus_hits_r2": DATA / "water_plus_hits_R2.fastq.gz",
}
# WHO risk groups the mini BLAST fixtures can report. Kept here as a named
# constant so a new expected_organism fails with a clear lookup error instead
# of a bare KeyError on a dict literal buried inside an assertion loop.
WHO_RISK_GROUP_BY_ORGANISM = {
    "Orf virus": "RG2",
    "Monkeypox virus": "RG3",
}


def assert_detections_agree_with_the_full_size_run(
    experiment_rows: list[dict[str, str]],
    selected_sra_runs: list[dict[str, Any]],
) -> None:
    """Check the reduced samplesheet still detects what all five runs detected.

    The samplesheet keeps one SRA run per distinct code path rather than the
    original five. The per-sample assertions above prove each remaining run
    finds the organism it must; this proves the reduction did not *add* a
    detection, which would mean the pipeline behaves differently on the smaller
    input. Detecting less is expected and allowed -- see tests/e2e_parity.py.

    Only the SRA rows are compared. The synthetic metagenomes have no
    counterpart in the recording (they were added after it) and carry their own
    assertions, so including them would be an unmapped-sample error.

    This runs only in the full tier, because only there are the reads the same
    ones the recording was made from. The offline tier reads *downsampled*
    copies, and fewer reads mean less assembly, which pushes more sequences onto
    the read-level BLAST route where short spurious hits surface. A measured
    example: the offline ONT fixture reports taxid 10244 on six rows out of 716,
    every one a 31 bp alignment against a conserved poxvirus region, which the
    full run does not report. Raising that fixture's fraction 2.5x did not
    remove it, so it is a property of comparing unlike inputs rather than
    something a bigger fixture fixes. Comparing a downsampled run against a
    full-read recording would therefore fail for a reason that is not a
    behaviour change, so the comparison belongs where it is like-for-like.
    """
    if not integration_sra_downloads_enabled():
        print(  # noqa: T201
            "PARITY SKIPPED: the offline tier reads downsampled fixtures, which "
            "are not the reads the recording was made from. Run `just e2e-full` "
            "to compare against it.",
            flush=True,
        )
        return

    recorded_detections = load_parity()["detections"]
    sra_sample_ids = {str(run_info["sample_id"]) for run_info in selected_sra_runs}
    observed_detections = {
        sample_id: detection
        for sample_id, detection in summarize_detections(experiment_rows).items()
        if sample_id in sra_sample_ids
    }
    # The reduction kept the samplesheet's sample ids, so the recording is keyed
    # by the same names and the mapping is the identity.
    sample_map = {sample_id: sample_id for sample_id in sra_sample_ids}

    errors, warnings = find_parity_violations(
        observed_detections,
        recorded_detections,
        sample_map,
    )
    for warning in warnings:
        print(f"PARITY WARNING: {warning}")  # noqa: T201
    assert not errors, "Detections disagree with the recorded full-size run:\n" + "\n".join(
        errors,
    )

# The synthetic wastewater metagenomes. Unlike the SRA rows these are built
# locally rather than downloaded, so they cost no network transfer, and unlike
# the water* fixtures they genuinely contain target-virus reads and therefore
# carry biological expectations. They are the only shotgun metagenomic samples
# in the suite, which is the pipeline's real target case: the SRA rows are all
# hybrid-capture or cell-culture libraries.
#
# Both were simulated from the exact accessions in mini_virus_reference.fasta,
# which holds complete genomes, so reads from anywhere in a genome can be
# detected rather than only those overlapping a slice.
# Committed to the repository rather than downloaded or generated on demand, so
# the end-to-end test runs offline with no network and no build step. They are
# ~11.2 MB of the repository, which buys the only shotgun metagenomic coverage
# in the suite. The builder stays available for regeneration, and the pinned
# checksums below now verify the committed bytes rather than a fresh build.
SYNTHETIC_FIXTURE_DIR = Path(
    os.environ.get(
        "NVD_SYNTHETIC_FIXTURE_DIR",
        DATA,
    ),
)
SYNTHETIC_BUILDER = ROOT / "tests" / "scripts" / "build_synthetic_metagenomes.py"
# Both of these are part of what the checksums below pin. The worker count
# matters because InSilicoSeq seeds one worker per CPU, so the same seed at a
# different worker count yields different reads.
# One run per distinct code path fits in this; three paired-Illumina Orf runs
# did not. Asserted by test_samplesheet_keeps_one_run_per_distinct_code_path so
# a future row cannot quietly restore the old 420 MB download.
E2E_DOWNLOAD_BUDGET_MB = 150
SYNTHETIC_FIXTURE_SEED = 42
SYNTHETIC_FIXTURE_CPUS = 4
# Pinned so a drifting build fails loudly here rather than quietly changing
# what the end-to-end test is asserting against. The builder is deterministic:
# same seed, same reads, same bytes.
SYNTHETIC_FIXTURE_SHA256 = {
    "synthetic_ww_lowcx_R1.fastq.gz": "f5266f2a847596c748f24756aac5044103088859a2da6c9f578600a8d42118f8",
    "synthetic_ww_lowcx_R2.fastq.gz": "dc33caffcd6597898dc98fb96a9a0ff8caa3142d06c749b5b64d5a6eb64b5fc3",
    "synthetic_ww_highcx_R1.fastq.gz": "ecef42aada5f84328f0069b91f3e01d5248d2c03e3a4ee7a46dade09b3222156",
    "synthetic_ww_highcx_R2.fastq.gz": "656e61b7b1d8041807a2cc9ec15ce530c40f7898e23691822a85c2a9d336eacf",
}
SYNTHETIC_METAGENOME_SAMPLES = (
    {
        "sample_id": "synthetic_ww_lowcx",
        # Orf at about 6x depth, so this sample assembles and exercises the
        # contig BLAST route.
        "expected_taxa": ({"expected_organism": "Orf virus", "taxid": 10258},),
        "fastq1": SYNTHETIC_FIXTURE_DIR / "synthetic_ww_lowcx_R1.fastq.gz",
        "fastq2": SYNTHETIC_FIXTURE_DIR / "synthetic_ww_lowcx_R2.fastq.gz",
    },
    {
        "sample_id": "synthetic_ww_highcx",
        # Orf at 0.46% and Mpox at 0.37% of reads, below the assembly floor, so
        # this sample exercises detection from unassembled reads instead. That
        # is the abundance real wastewater shows: the modelled run found
        # Poxviridae on 21 reads across 18 samples.
        #
        # Read-level BLAST is therefore this sample's ONLY detection route, so
        # the skip_unassembled_read_queries variant necessarily finds nothing
        # here. That is the variant working, not a regression, and the flag is
        # what keeps the two facts coupled in one place.
        "detects_only_from_unassembled_reads": True,
        "expected_taxa": (
            {"expected_organism": "Orf virus", "taxid": 10258},
            {"expected_organism": "Monkeypox virus", "taxid": 10244},
        ),
        "fastq1": SYNTHETIC_FIXTURE_DIR / "synthetic_ww_highcx_R1.fastq.gz",
        "fastq2": SYNTHETIC_FIXTURE_DIR / "synthetic_ww_highcx_R2.fastq.gz",
    },
)

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
3431483\t|\t10242\t|\tspecies\t|\t
10244\t|\t3431483\t|\tno rank\t|\t
10255\t|\t10240\t|\tgenus\t|\t
3431389\t|\t10255\t|\tspecies\t|\t
10258\t|\t3431389\t|\tno rank\t|\t
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
3431483\t|\tOrthopoxvirus monkeypox\t|\t\t|\tscientific name\t|
10244\t|\tMonkeypox virus\t|\t\t|\tscientific name\t|
10255\t|\tParapoxvirus\t|\t\t|\tscientific name\t|
3431389\t|\tParapoxvirus orf\t|\t\t|\tscientific name\t|
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


def read_fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current_id: str | None = None
    sequence_lines = 0
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line.startswith(">"):
            current_id = raw_line[1:].split(maxsplit=1)[0]
            assert current_id not in records, (
                f"Duplicate FASTA ID {current_id} in {path}"
            )
            records[current_id] = ""
            sequence_lines = 0
        else:
            assert current_id is not None, f"Sequence before FASTA header in {path}"
            sequence_lines += 1
            assert sequence_lines == 1, (
                f"Wrapped FASTA sequence for {current_id} in {path}"
            )
            records[current_id] += raw_line
    return records


def assert_best_hit_sequence_evidence_outputs(results_root: Path) -> None:
    best_hit_sequences_dir = results_root / "10_best_hit_sequences"
    expected_artifacts = {
        "query_sequences.fasta",
        "selected_references.fasta",
        "best_hit_placements.bed",
    }
    assert best_hit_sequences_dir.is_dir(), (
        f"Missing best-hit sequence evidence directory: {best_hit_sequences_dir}"
    )
    assert {
        path.name
        for path in best_hit_sequences_dir.iterdir()
        if path.is_file() and path.name != "versions.yml"
    } == expected_artifacts
    qbt_rows = read_tsv_rows(results_root / "query_big_table.tsv")
    assert qbt_rows, "Missing Query Big Table rows for best-hit sequence evidence"
    query_fasta_ids = set(
        read_fasta_records(best_hit_sequences_dir / "query_sequences.fasta"),
    )
    reference_fasta_ids = set(
        read_fasta_records(best_hit_sequences_dir / "selected_references.fasta"),
    )
    expected_query_ids = {f"{row['sample_id']}|{row['qseqid']}" for row in qbt_rows}
    expected_reference_ids = {row["best_hit_reference_accession"] for row in qbt_rows}
    assert query_fasta_ids == expected_query_ids
    assert reference_fasta_ids == expected_reference_ids

    strand_map = {"plus": "+", "minus": "-", "+": "+", "-": "-"}
    expected_bed_rows = sorted(
        [
            [
                row["best_hit_reference_accession"],
                str(
                    min(
                        int(row["best_hit_reference_start_1based"]),
                        int(row["best_hit_reference_end_1based"]),
                    )
                    - 1,
                ),
                str(
                    max(
                        int(row["best_hit_reference_start_1based"]),
                        int(row["best_hit_reference_end_1based"]),
                    ),
                ),
                f"{row['sample_id']}|{row['qseqid']}",
                "0",
                strand_map[row["best_hit_reference_strand"]],
            ]
            for row in qbt_rows
        ],
        key=lambda row: (row[0], int(row[1]), int(row[2]), row[3], row[5]),
    )
    observed_bed_rows = [
        line.split("\t")
        for line in (best_hit_sequences_dir / "best_hit_placements.bed")
        .read_text(encoding="utf-8")
        .splitlines()
        if line
    ]
    assert observed_bed_rows == expected_bed_rows


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


def expected_organism_names(
    selected_sra_runs: list[dict[str, Any]],
) -> set[str]:
    """Every organism the run should name, across SRA and synthetic samples."""
    return {str(run_info["expected_organism"]) for run_info in selected_sra_runs} | {
        str(taxon["expected_organism"])
        for sample in SYNTHETIC_METAGENOME_SAMPLES
        for taxon in sample["expected_taxa"]
    }


def assert_synthetic_metagenome_taxa_detected(
    experiment_rows: list[dict[str, str]],
) -> None:
    """Assert each synthetic metagenome reports the taxa buried in it.

    These are the suite's only shotgun samples, and the point of including them
    is that a target virus inside a real community is still found. Checked per
    expected taxon, because a metagenome carries several rather than the one an
    SRA row declares.

    Every tier asserts every sample's taxa in full. The single exception is the
    opt-in skip_unassembled_read_queries variant, which turns off read-level
    querying: a sample whose only detection route is read-level BLAST is then
    checked the opposite way, because the variant removed the route rather than
    the pipeline failing to use it.
    """
    for sample in SYNTHETIC_METAGENOME_SAMPLES:
        sample_id = str(sample["sample_id"])
        reads_only = bool(sample.get("detects_only_from_unassembled_reads"))
        if reads_only and integration_skip_unassembled_read_queries_enabled():
            assert_sample_reports_nothing_without_read_queries(
                experiment_rows,
                sample_id=sample_id,
            )
            continue
        assert_expected_taxa_detected(
            experiment_rows,
            sample_id=sample_id,
            expected_taxa=sample["expected_taxa"],
        )


def assert_sample_reports_nothing_without_read_queries(
    experiment_rows: list[dict[str, str]],
    *,
    sample_id: str,
) -> None:
    """Assert a reads-only sample finds nothing once read querying is turned off.

    Two different things in this suite both talk about unassembled reads, and
    confusing them is easy. This sample *detects from* unassembled reads, which
    is its only route. The --skip_unassembled_read_queries flag *turns off*
    querying unassembled reads. So that one opt-in variant removes this sample's
    only route, and finding nothing is the flag working.

    This runs ONLY under that variant. Every other tier, the default e2e
    included, asserts this sample's taxa in full: the normal run finds it on
    1,279 single-read queries carrying both expected taxids.

    Deliberately an assertion rather than a skip, so it proves the flag took
    effect instead of letting the sample pass for some unrelated reason.
    """
    unexpected_rows = [row for row in experiment_rows if row.get("sample") == sample_id]

    assert not unexpected_rows, (
        f"{sample_id} detects only from unassembled reads, but the "
        f"skip_unassembled_read_queries run still reported {len(unexpected_rows)} "
        "rows for it, so the flag did not take effect."
    )
    print(
        "READ-LEVEL QUERYING IS OFF IN THIS VARIANT, SO "
        f"{sample_id} CORRECTLY REPORTED NO ROWS (its taxa are asserted in full "
        "by every other tier)",
    )


def assert_expected_taxa_detected(
    experiment_rows: list[dict[str, str]],
    *,
    sample_id: str,
    expected_taxa: tuple[dict[str, Any], ...],
) -> None:
    """Assert a sample reports each taxon it is expected to contain.

    Shaped for samples with more than one expected taxon, which is the normal
    case for a metagenome and why this does not reuse the single-organism
    checks the SRA rows run.
    """
    sample_rows = [row for row in experiment_rows if row.get("sample") == sample_id]
    assert sample_rows, f"No experiment BLAST rows found for {sample_id}"

    for taxon in expected_taxa:
        organism = str(taxon["expected_organism"])
        taxid = str(taxon["taxid"])
        expected_risk_group = WHO_RISK_GROUP_BY_ORGANISM[organism]

        assert any(row.get("staxids") == taxid for row in sample_rows), (
            f"No {taxid} ({organism}) BLAST taxid found for {sample_id}"
        )
        assert any(
            row.get("adjusted_taxid") == taxid
            and row.get("adjusted_taxid_name") == organism
            for row in sample_rows
        ), f"No {organism} adjusted taxon found for {sample_id}"
        assert any(
            row.get("adjusted_taxid") == taxid
            and row.get("who_risk_group") == expected_risk_group
            for row in sample_rows
        ), (
            f"No {expected_risk_group} WHO risk group found for {taxid} "
            f"({organism}) in {sample_id}"
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


def assert_successful_nvd_multiqc_outputs(
    results_root: Path,
    *,
    experimental: bool,
    skip_assembly: bool,
) -> None:
    """Assert the healthy fixture completed ancillary reporting and publication."""
    report = results_root / "multiqc_report.html"
    data = results_root / "12_experiment_summary" / "multiqc_data"
    nvd_inputs = data / "nvd_inputs"
    manifest_path = nvd_inputs / "nvd_report_manifest.json"
    raw_fastqc = results_root / "00_input_preparation" / "raw_fastq_qc" / "fastqc"

    assert report.is_file(), f"Missing NVD MultiQC report: {report}"
    assert report.stat().st_size > 0, f"Empty NVD MultiQC report: {report}"
    assert data.is_dir(), f"Missing NVD MultiQC data directory: {data}"
    assert any(data.iterdir()), f"Empty NVD MultiQC data directory: {data}"
    assert manifest_path.is_file(), (
        f"Missing retained NVD MultiQC manifest: {manifest_path}"
    )
    assert (nvd_inputs / "nvd_sample_roster_mqc.yaml").is_file()
    expected_domain_inputs = [
        "nvd_target_enrichment_mqc.yaml",
        "nvd_depletion_mqc.yaml",
        "nvd_fastx_profiles_mqc.yaml",
        "nvd_fastx_single_read_length_distribution_mqc.yaml",
        "nvd_fastx_quality_distribution_mqc.yaml",
        "nvd_assembly_mqc.yaml",
        "nvd_prepared_blast_query_batches_mqc.yaml",
    ]
    if experimental:
        expected_domain_inputs.append(
            "nvd_fastx_overlap_merged_pair_length_distribution_mqc.yaml",
        )
    if not skip_assembly:
        expected_domain_inputs.append(
            "nvd_fastx_filtered_contigs_length_distribution_mqc.yaml",
        )
    for filename in expected_domain_inputs:
        assert (nvd_inputs / filename).is_file(), (
            f"Missing retained NVD MultiQC input: {filename}"
        )
    assert raw_fastqc.is_dir(), f"Missing retained raw FastQC directory: {raw_fastqc}"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    raw_units = manifest.get("raw_fastqc")
    assert raw_units, (
        f"No observed raw FastQC units in retained manifest: {manifest_path}"
    )
    expected_zips = sorted(str(unit["zip_alias"]) for unit in raw_units)
    expected_htmls = sorted(str(unit["html_alias"]) for unit in raw_units)
    observed_zips = sorted(path.name for path in raw_fastqc.glob("*.raw.*_fastqc.zip"))
    observed_htmls = sorted(
        path.name for path in raw_fastqc.glob("*.raw.*_fastqc.html")
    )
    assert observed_zips == expected_zips
    assert observed_htmls == expected_htmls


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


def integration_sra_downloads_enabled() -> bool:
    """Whether to fetch the virus samples from SRA instead of using the fixtures.

    Off by default, which is what makes the e2e offline. The full tier turns it
    on so the FETCH_FASTQ/sracha download path stays covered; nothing else
    exercises it, because every other row is a local file.
    """
    return os.environ.get("NVD_E2E_USE_SRA") == "1"


def integration_skip_assembly_enabled() -> bool:
    return os.environ.get("NVD_INTEGRATION_SKIP_ASSEMBLY") == "1"


def integration_skip_unassembled_read_queries_enabled() -> bool:
    # Read querying is on by default from v3.4.0, so the contig-only path needs
    # its own switch to stay covered.
    return os.environ.get("NVD_INTEGRATION_SKIP_UNASSEMBLED_READ_QUERIES") == "1"


def synthetic_builder_command() -> list[str]:
    """Build the builder invocation, pinning everything the reads depend on.

    Seed and worker count are both passed explicitly rather than left to the
    builder's defaults. InSilicoSeq seeds one worker per CPU from the base
    seed, so the worker count changes the generated reads: a default that moved
    would silently produce different fixtures that no longer match
    SYNTHETIC_FIXTURE_SHA256.
    """
    return [
        sys.executable,
        str(SYNTHETIC_BUILDER),
        "--out-dir",
        str(SYNTHETIC_FIXTURE_DIR),
        "--cache-dir",
        str(SYNTHETIC_FIXTURE_DIR.parent / "refs"),
        "--work-dir",
        str(SYNTHETIC_FIXTURE_DIR.parent / "work"),
        "--seed",
        str(SYNTHETIC_FIXTURE_SEED),
        "--cpus",
        str(SYNTHETIC_FIXTURE_CPUS),
    ]


def sha256_of_file(path: Path) -> str:
    """Return the hex sha256 of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def synthetic_fixtures_are_current() -> bool:
    """Report whether every synthetic FASTQ is present with its pinned sha256."""
    return all(
        (SYNTHETIC_FIXTURE_DIR / name).is_file()
        and sha256_of_file(SYNTHETIC_FIXTURE_DIR / name) == expected_sha256
        for name, expected_sha256 in SYNTHETIC_FIXTURE_SHA256.items()
    )


def ensure_synthetic_metagenomes() -> None:
    """Verify the committed synthetic fixtures, rebuilding only if they are absent.

    Normally this returns immediately: the fixtures are committed to
    ``tests/data`` and match their pinned checksums, so the end-to-end test
    needs neither the network nor a build step. The rebuild path remains for a
    regeneration or a custom ``NVD_SYNTHETIC_FIXTURE_DIR``, and it does need
    network access to fetch the RefSeq genomes it simulates from. The builder is
    deterministic, so a rebuild reproduces the pinned checksums byte for byte.
    """
    if synthetic_fixtures_are_current():
        return

    SYNTHETIC_FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    print(
        f"Building synthetic metagenome fixtures in {SYNTHETIC_FIXTURE_DIR} "
        "(about 40s, cached for later runs)",
        flush=True,
    )
    completed = subprocess.run(  # noqa: S603
        synthetic_builder_command(),
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, (
        "Failed to build the synthetic metagenome fixtures.\n\nSTDOUT:\n"
        + completed.stdout
        + "\n\nSTDERR:\n"
        + completed.stderr
    )

    mismatched = {
        name: sha256_of_file(SYNTHETIC_FIXTURE_DIR / name)
        for name, expected_sha256 in SYNTHETIC_FIXTURE_SHA256.items()
        if sha256_of_file(SYNTHETIC_FIXTURE_DIR / name) != expected_sha256
    }
    assert not mismatched, (
        "The synthetic metagenome build did not reproduce its pinned "
        f"checksums: {mismatched}. Either the builder or a tool version "
        "changed; re-pin SYNTHETIC_FIXTURE_SHA256 only after confirming why."
    )


def local_sample_row(sample: dict[str, Any], local_fastq_dir: Path) -> dict[str, str]:
    fastq1 = sample.get("fastq1")
    fastq2 = sample.get("fastq2")
    fastq1_glob = sample.get("fastq1_glob")
    fastq2_glob = sample.get("fastq2_glob")
    return {
        "sample_id": str(sample["sample_id"]),
        "srr": "",
        # Defaults to illumina for the local fixtures that predate this
        # argument. The offline SRA rows set it explicitly, and the ONT one must
        # keep its platform or the long-read route would never be exercised.
        "platform": str(sample.get("platform", "illumina")),
        "fastq1": str(fastq1.resolve()) if isinstance(fastq1, Path) else "",
        "fastq2": str(fastq2.resolve()) if isinstance(fastq2, Path) else "",
        "fastq1_glob": str(local_fastq_dir / fastq1_glob)
        if isinstance(fastq1_glob, str)
        else "",
        "fastq2_glob": str(local_fastq_dir / fastq2_glob)
        if isinstance(fastq2_glob, str)
        else "",
    }


def load_offline_sra_manifest() -> dict[str, Any]:
    """Read the provenance manifest for the committed SRA downsamples."""
    with OFFLINE_SRA_MANIFEST_PATH.open(encoding="utf-8") as handle:
        return json.load(handle)


def offline_sra_sample_rows() -> list[dict[str, Any]]:
    """Build local-FASTQ rows from the committed downsamples of the SRA runs.

    Sample ids match the full runs', so the parity comparison stays an identity
    map and every downstream assertion keeps working; only the file names carry
    the _s42 provenance marker. See the manifest's maintainer_note.
    """
    rows: list[dict[str, Any]] = []
    for fixture in load_offline_sra_manifest()["fixtures"]:
        file_names = [file_entry["name"] for file_entry in fixture["files"]]
        row: dict[str, Any] = {
            "sample_id": fixture["sample_id"],
            "platform": fixture["platform"],
            "fastq1": DATA / file_names[0],
        }
        if len(file_names) > 1:
            row["fastq2"] = DATA / file_names[1]
        rows.append(row)
    return rows


def write_augmented_samplesheet(run_dir: Path) -> Path:
    """Write the e2e samplesheet: local FASTQ rows, or SRA accessions on request.

    By default every row points at a committed fixture, so the run needs no
    network. With NVD_E2E_USE_SRA=1 the three virus samples become SRA
    accession rows instead, which is what keeps the FETCH_FASTQ/sracha download
    path under test in the full tier.
    """
    ensure_synthetic_metagenomes()
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
    if integration_sra_downloads_enabled():
        with SAMPLESHEET.open(newline="", encoding="utf-8") as handle:
            rows.extend(
                {column: row.get(column, "") for column in fieldnames}
                for row in csv.DictReader(handle)
            )
    else:
        rows.extend(
            local_sample_row(sample, local_fastq_dir)
            for sample in offline_sra_sample_rows()
        )

    rows.extend(
        local_sample_row(sample, local_fastq_dir) for sample in LOCAL_E2E_SAMPLES
    )
    rows.extend(
        local_sample_row(sample, local_fastq_dir)
        for sample in SYNTHETIC_METAGENOME_SAMPLES
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
            ],
        )
    if skip_assembly:
        command.extend(["--skip_assembly", "true"])
    if integration_skip_unassembled_read_queries_enabled():
        command.extend(["--skip_unassembled_read_queries", "true"])
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
        command += [
            "--experimental",
            "true",
            "--merge_pairs",
            "true",
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
        if op["table"] == "hits" and any(row.get("query_class") for row in op["rows"])
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
            assert (
                results_root
                / "13_labkey_uploads"
                / "upload_logs"
                / "final_labkey_upload.log"
            ).exists(), f"Missing final LabKey upload log under {results_root}"

            first_hits = data_hits_inserts(mock)
            assert first_hits, "no eager per-batch hits insert carrying query_class"
            for op in first_hits:
                for row in op["rows"]:
                    assert row.get("query_class"), (
                        f"hits row missing query_class: {row}"
                    )

            observed_classes = {
                row.get("query_class") for op in first_hits for row in op["rows"]
            }
            print(
                f"query_class values uploaded: {sorted(observed_classes)}", flush=True
            )
            assert "short_assembly_contig" in observed_classes, (
                f"expected the contig query class; observed {sorted(observed_classes)}"
            )
            assert observed_classes & {"overlap_merged_pair", "single_read"}, (
                "expected read-derived query classes too (queries are split by "
                f"read type by default); observed {sorted(observed_classes)}"
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
    experiment_id = (
        int(pinned) if pinned else 2_000_000_000 + secrets.randbelow(147_000_000)
    )

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
        assert hits_after_first, (
            "no hits present for the test experiment after the first run"
        )
        for row in hits_after_first:
            assert row.get("query_class"), f"hits row missing query_class: {row}"
        observed_classes = {row.get("query_class") for row in hits_after_first}
        print(f"query_class values present: {sorted(observed_classes)}", flush=True)
        assert "short_assembly_contig" in observed_classes, (
            f"expected the contig query class; observed {sorted(observed_classes)}"
        )
        if experimental:
            assert observed_classes & {"overlap_merged_pair", "single_read"}, (
                "expected read-derived query classes too (queries are split by "
                f"read type by default); observed {sorted(observed_classes)}"
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
                delete_experiment_rows(
                    api, cfg["schema"], cfg["hits_list"], experiment_id
                )
                delete_experiment_rows(
                    api, cfg["schema"], cfg["fasta_list"], experiment_id
                )
                delete_webdav_experiment(cfg["webdav"], cfg["api_key"], experiment_id)
            except Exception:  # noqa: BLE001, S110 - cleanup is best-effort; never mask a failure
                pass


@pytest.mark.slow
@pytest.mark.e2e
def test_viral_pipeline_completes_end_to_end() -> None:
    """The whole pipeline should run through enrichment, assembly and BLAST.

    Not marked `network`: by default every read comes from a committed fixture
    and the run needs none, which `unshare -rn` verifies. `just e2e-full` sets
    NVD_E2E_USE_SRA=1 to fetch the source accessions instead, and that tier is
    the only thing covering FETCH_FASTQ/sracha.
    """
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
        for run_info in (
            *LOCAL_E2E_SAMPLES,
            *SYNTHETIC_METAGENOME_SAMPLES,
            *selected_sra_runs,
        )
    }
    assert_successful_nvd_multiqc_outputs(
        results_root,
        experimental=experimental,
        skip_assembly=skip_assembly,
    )
    assert_target_enrichment_outputs(results_root, expected_sample_ids)
    assert_read_profiles_respect_length_filter(results_root)
    resolved_manifest = (
        results_root
        / "00_input_preparation"
        / "input_resolution"
        / "resolved_reads.jsonl"
    )
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
    merged_blast_dir = results_root / "07_merged_blast_results"
    final_dir = merged_blast_dir / "final"
    final_blast_files = sorted(final_dir.glob("*_blast.final.tsv"))

    experiment_blast = (
        results_root / "12_experiment_summary" / "experiment_blast_results.tsv"
    )

    if skip_assembly and integration_skip_unassembled_read_queries_enabled():
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
        assert "who_risk_group" in experiment_rows[0]

        # Assert per-sample biological expectations for the rows actually under
        # test. Coupling this loop to every manifest row makes a deleted
        # samplesheet row fail as a missing output, even though the pipeline did
        # exactly what the samplesheet requested.
        for run_info in selected_sra_runs:
            expected_risk_group = WHO_RISK_GROUP_BY_ORGANISM[
                run_info["expected_organism"]
            ]
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
                row.get("adjusted_taxid") == str(run_info["taxid"])
                and row.get("who_risk_group") == expected_risk_group
                for row in sample_rows
            ), (
                f"No {expected_risk_group} WHO risk group found for "
                f"{run_info['taxid']} in {run_info['sample_id']}"
            )
            assert any(
                row.get("adjusted_taxid") == str(run_info["taxid"])
                and row.get("adjusted_taxid_name") == run_info["expected_organism"]
                for row in sample_rows
            ), (
                f"No {run_info['expected_organism']} adjusted taxon found for "
                f"{run_info['sample_id']}"
            )
            expected_tasks = run_info.get("expected_tasks", [])
            observed_tasks = {row.get("task") for row in sample_rows}
            assert set(expected_tasks) <= observed_tasks, (
                f"Missing expected BLAST tasks for {run_info['sample_id']}: "
                f"expected {sorted(expected_tasks)}, observed {sorted(observed_tasks)}"
            )

        assert_synthetic_metagenome_taxa_detected(experiment_rows)

        assert_detections_agree_with_the_full_size_run(
            experiment_rows,
            selected_sra_runs,
        )

        for organism in expected_organism_names(selected_sra_runs):
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
            rows = read_tsv_rows(featured_table)
            assert rows, f"No rows in featured Big Table: {featured_table}"
            columns = list(rows[0])
            taxid_column = (
                "assigned_taxid" if filename == "query_big_table.tsv" else "taxid"
            )
            assert columns.index("who_risk_group") == columns.index(taxid_column) + 1
            if filename == "query_big_table.tsv":
                placement_columns = {
                    "best_hit_reference_accession",
                    "best_hit_reference_title",
                    "best_hit_alignment_length",
                    "best_hit_query_start_1based",
                    "best_hit_query_end_1based",
                    "best_hit_reference_length",
                    "best_hit_reference_start_1based",
                    "best_hit_reference_end_1based",
                    "best_hit_reference_strand",
                }
                assert placement_columns <= set(columns)
                for row in rows:
                    assert row["best_hit_reference_accession"]
                    assert int(row["best_hit_query_start_1based"]) <= int(
                        row["best_hit_query_end_1based"],
                    )
                    assert int(row["best_hit_reference_start_1based"]) <= int(
                        row["best_hit_reference_end_1based"],
                    )

        assert_best_hit_sequence_evidence_outputs(results_root)


def test_synthetic_metagenomes_target_accessions_are_in_the_mini_blast_db() -> None:
    """A synthetic target is only findable if the mini BLAST DB carries it.

    The mini reference holds complete genomes, so this also means reads from
    anywhere in a target genome can be detected, not only those overlapping a
    slice of it.
    """
    taxid_map_rows = (
        (DATA / "mini_virus_blast_taxid_map.tsv")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    mini_db_taxids = {row.split("\t")[1] for row in taxid_map_rows if row}

    for sample in SYNTHETIC_METAGENOME_SAMPLES:
        for taxon in sample["expected_taxa"]:
            assert str(taxon["taxid"]) in mini_db_taxids, (
                f"{sample['sample_id']} expects taxid {taxon['taxid']}, which "
                "the mini BLAST database cannot report"
            )


def test_every_expected_organism_has_a_declared_risk_group() -> None:
    """Guards the lookup that used to be a dict literal inside the assertions."""
    declared_organisms = {
        str(run_info["expected_organism"]) for run_info in load_manifest()["sra_runs"]
    } | {
        str(taxon["expected_organism"])
        for sample in SYNTHETIC_METAGENOME_SAMPLES
        for taxon in sample["expected_taxa"]
    }
    missing = declared_organisms - set(WHO_RISK_GROUP_BY_ORGANISM)
    assert not missing, f"No WHO risk group declared for {missing}"


def test_samplesheet_keeps_one_run_per_distinct_code_path() -> None:
    """The full e2e should not pay to download the same code path three times.

    It previously ran three paired-Illumina Orf runs, 420 MB in total. One run
    per (platform, organism) is enough for the SRA rows, because the shotgun
    metagenomic case is covered by the locally built synthetic samples.
    """
    manifest_by_sample = {
        str(run_info["sample_id"]): run_info for run_info in load_manifest()["sra_runs"]
    }
    with SAMPLESHEET.open(newline="", encoding="utf-8") as handle:
        selected = [
            manifest_by_sample[row["sample_id"]]
            for row in csv.DictReader(handle)
            if row["sample_id"] in manifest_by_sample
        ]

    code_paths = [
        (run_info["platform"], run_info["expected_organism"]) for run_info in selected
    ]
    assert len(code_paths) == len(set(code_paths)), (
        f"Duplicated code paths in the e2e samplesheet: {code_paths}"
    )
    assert sum(run_info["size_mb"] for run_info in selected) <= E2E_DOWNLOAD_BUDGET_MB


def test_synthetic_samples_are_the_only_shotgun_metagenomic_rows() -> None:
    """At least one metagenomic sample must be present; the SRA rows are not.

    Every SRA row is hybrid-capture, pre-extracted reads, or cell culture, so
    without these the suite would never test the pipeline's actual target case.
    """
    assert len(SYNTHETIC_METAGENOME_SAMPLES) >= 1
    sample_ids = [str(sample["sample_id"]) for sample in SYNTHETIC_METAGENOME_SAMPLES]
    assert len(sample_ids) == len(set(sample_ids))
    # One sample must sit above the assembly floor and one below it, so both
    # the contig route and the unassembled-read route are exercised.
    assert {
        len(sample["expected_taxa"]) for sample in SYNTHETIC_METAGENOME_SAMPLES
    } == {
        1,
        2,
    }


def test_synthetic_fixture_checksums_cover_every_declared_file() -> None:
    """Every FASTQ a synthetic row points at must have a pinned checksum."""
    referenced = {
        Path(str(sample[mate])).name
        for sample in SYNTHETIC_METAGENOME_SAMPLES
        for mate in ("fastq1", "fastq2")
    }
    assert referenced == set(SYNTHETIC_FIXTURE_SHA256)


def test_synthetic_sample_rows_carry_absolute_paired_paths(tmp_path: Path) -> None:
    """The samplesheet needs absolute paths, since Nextflow runs elsewhere."""
    for sample in SYNTHETIC_METAGENOME_SAMPLES:
        row = local_sample_row(sample, tmp_path)
        assert row["platform"] == "illumina"
        assert row["srr"] == ""
        assert Path(row["fastq1"]).is_absolute()
        assert Path(row["fastq2"]).is_absolute()
        assert row["fastq1_glob"] == ""


def test_exactly_the_trace_level_metagenome_is_flagged_as_reads_only() -> None:
    """The reads-only flag must stay pinned to the sample whose biology earns it.

    The flag is what lets the skip_unassembled_read_queries variant pass, so a
    fixture edit that spread it to an assembling sample would silently stop
    checking that sample's taxa. Pin both directions.
    """
    flagged = {
        str(sample["sample_id"])
        for sample in SYNTHETIC_METAGENOME_SAMPLES
        if sample.get("detects_only_from_unassembled_reads")
    }

    assert flagged == {"synthetic_ww_highcx"}, (
        "Only the trace-level metagenome detects solely from unassembled reads. "
        f"Flagged: {sorted(flagged)}"
    )
    # The low-coverage sample assembles, so it must keep its taxa assertions in
    # every variant. This is the half that a careless edit would break.
    assert not any(
        sample.get("detects_only_from_unassembled_reads")
        for sample in SYNTHETIC_METAGENOME_SAMPLES
        if sample["sample_id"] == "synthetic_ww_lowcx"
    )


def test_assert_expected_taxa_detected_requires_each_declared_taxon() -> None:
    rows = [
        {
            "sample": "synthetic_ww_highcx",
            "staxids": "10258",
            "adjusted_taxid": "10258",
            "adjusted_taxid_name": "Orf virus",
            "who_risk_group": "RG2",
        },
    ]
    orf_only = ({"expected_organism": "Orf virus", "taxid": 10258},)
    assert_expected_taxa_detected(
        rows,
        sample_id="synthetic_ww_highcx",
        expected_taxa=orf_only,
    )

    # Mpox is declared but absent from the rows, so this must fail.
    with pytest.raises(AssertionError, match="10244"):
        assert_expected_taxa_detected(
            rows,
            sample_id="synthetic_ww_highcx",
            expected_taxa=(
                *orf_only,
                {"expected_organism": "Monkeypox virus", "taxid": 10244},
            ),
        )

    with pytest.raises(AssertionError, match="No experiment BLAST rows"):
        assert_expected_taxa_detected(
            rows,
            sample_id="absent_sample",
            expected_taxa=orf_only,
        )


def test_offline_sra_rows_point_at_committed_files_and_keep_their_platform() -> None:
    """The offline rows must be local files, and the ONT row must stay ONT.

    local_sample_row used to hardcode illumina, which would have processed the
    ONT fixture as a short-read sample and quietly skipped the long-read
    assembly route the fixture exists to cover.
    """
    rows = offline_sra_sample_rows()
    assert len(rows) == len(load_offline_sra_manifest()["fixtures"])

    platforms = {}
    for sample in rows:
        row = local_sample_row(sample, Path("/unused"))
        assert not row["srr"], (
            f"{row['sample_id']} still carries an SRA accession, so the offline "
            "run would try to download it"
        )
        assert Path(row["fastq1"]).is_file(), (
            f"{row['sample_id']} points at a missing file: {row['fastq1']}"
        )
        platforms[row["sample_id"]] = row["platform"]

    assert "ont" in platforms.values(), (
        f"No offline row is ONT, so the long-read route is uncovered: {platforms}"
    )


def test_the_offline_samplesheet_matches_the_sra_samplesheet_sample_ids() -> None:
    """Switching tiers must change only where reads come from, not which samples.

    The parity check compares against a recording keyed on these ids, so a
    mismatch here would silently turn that comparison into a no-op.
    """
    with SAMPLESHEET.open(newline="", encoding="utf-8") as handle:
        sra_sample_ids = {row["sample_id"] for row in csv.DictReader(handle)}
    offline_sample_ids = {str(row["sample_id"]) for row in offline_sra_sample_rows()}
    assert offline_sample_ids == sra_sample_ids


def test_parity_runs_in_the_full_tier_and_is_skipped_offline(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Parity must still run somewhere, and say so when it does not.

    The skip exists because the offline tier reads downsampled fixtures rather
    than the reads the recording was made from. If NVD_E2E_USE_SRA ever stopped
    enabling the comparison, the check would be dead rather than deferred, and
    nothing else would notice.
    """
    # A detection the recording cannot contain, so parity must object when it runs.
    impossible_rows = [
        {"sample": "orf_virus_ov_pt001_2024", "adjusted_taxid": "999999", "saccver": "X"},
    ]
    selected = [{"sample_id": "orf_virus_ov_pt001_2024"}]

    monkeypatch.delenv("NVD_E2E_USE_SRA", raising=False)
    assert_detections_agree_with_the_full_size_run(impossible_rows, selected)
    assert "PARITY SKIPPED" in capsys.readouterr().out

    monkeypatch.setenv("NVD_E2E_USE_SRA", "1")
    with pytest.raises(AssertionError, match="Detections disagree"):
        assert_detections_agree_with_the_full_size_run(impossible_rows, selected)


def test_sra_downloads_are_off_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The e2e must be offline unless explicitly asked for the full tier."""
    monkeypatch.delenv("NVD_E2E_USE_SRA", raising=False)
    assert not integration_sra_downloads_enabled()
    monkeypatch.setenv("NVD_E2E_USE_SRA", "1")
    assert integration_sra_downloads_enabled()


def test_synthetic_fixtures_are_committed_so_the_e2e_needs_no_build() -> None:
    """The fixtures must ship in the repository, so the e2e can run offline.

    This is the test that keeps the end-to-end test offline. If the fixtures
    ever stop being committed, or one drifts from its pinned checksum, the e2e
    would quietly fall back to the builder, which needs the network to fetch
    RefSeq genomes.
    """
    assert SYNTHETIC_FIXTURE_DIR == DATA, (
        "The synthetic fixtures must default to the committed copies in tests/data."
    )
    for name in SYNTHETIC_FIXTURE_SHA256:
        assert (DATA / name).is_file(), f"{name} is not committed to tests/data"
    assert synthetic_fixtures_are_current(), (
        "A committed synthetic fixture does not match its pinned sha256."
    )


def test_synthetic_fixtures_are_current_detects_a_corrupt_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A stale or truncated cache must trigger a rebuild rather than be used."""
    monkeypatch.setattr(
        "tests.test_nextflow_integration.SYNTHETIC_FIXTURE_DIR",
        tmp_path,
    )
    assert not synthetic_fixtures_are_current()

    for name in SYNTHETIC_FIXTURE_SHA256:
        (tmp_path / name).write_bytes(b"not the real reads")
    assert not synthetic_fixtures_are_current()


def test_synthetic_builder_command_pins_the_reproducibility_parameters() -> None:
    """Seed and worker count must be explicit, not inherited from defaults.

    Both change the generated reads, so a moved default would produce fixtures
    that silently no longer match SYNTHETIC_FIXTURE_SHA256.
    """
    command = synthetic_builder_command()

    assert command[command.index("--seed") + 1] == str(SYNTHETIC_FIXTURE_SEED)
    assert command[command.index("--cpus") + 1] == str(SYNTHETIC_FIXTURE_CPUS)
    for flag in ("--out-dir", "--cache-dir", "--work-dir"):
        assert Path(command[command.index(flag) + 1]).is_absolute()


def test_synthetic_builder_cpus_matches_the_builder_constant() -> None:
    """The e2e and the builder must agree on the worker count."""
    assert SYNTHETIC_FIXTURE_CPUS == REPRODUCIBLE_CPUS
