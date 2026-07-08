"""Tests for CRUMBS profile estimation."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING

import polars as pl
import pytest
from estimate_crumbs_profile import CrumbsProfileError, build_taxon_table, main
from hypothesis import given, settings
from hypothesis import strategies as st

if TYPE_CHECKING:
    from pathlib import Path

EXPECTED_CONTIG_A_CRUMBS = 50.0
EXPECTED_HOMO_TOTAL_CRUMBS = 200.0
EXPECTED_HOMO_TAXON_CRUMBS = 0.5
EXPECTED_PAN_TAXON_CRUMBS = 1.0
EXPECTED_ZERO_TAXON_CRUMBS = 0.0
EXPECTED_CONTIG_COUNT = 4
EXPECTED_BIOBOXES_TAXA = 2
EXPECTED_HOMO_PERCENTAGE = 33.33333333333333
EXPECTED_PAN_PERCENTAGE = 66.66666666666666
EXPECTED_FRAGMENTED_TAXON_COUNT = 2
BLAST_FIELDS = [
    "sample",
    "qseqid",
    "adjusted_taxid",
    "adjusted_taxid_name",
    "adjusted_taxid_rank",
    "adjustment_method",
]
COVERAGE_FIELDS = [
    "sample_id",
    "qseqid",
    "contig_length",
    "covered_bases_1x",
    "breadth_1x",
    "raw_aligned_bases",
    "mean_depth_full",
    "median_depth_full",
    "median_depth_positive",
    "depth_p95",
    "depth_p99",
    "max_depth",
    "crumbs_p95",
    "crumbs_p99",
]
TAXONOMY_FIELDS = ["taxon_id", "taxon_name", "rank", "taxpath", "taxpathsn", "rankpath"]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bioboxes_data_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        data_lines = [
            line for line in handle if not line.startswith("@") or line.startswith("@@")
        ]
    reader = csv.DictReader(data_lines, delimiter="\t")
    return list(reader)


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def write_minimal_inputs(
    tmp_path: Path,
    *,
    blast_rows: list[dict[str, object]] | None = None,
    coverage_rows: list[dict[str, object]] | None = None,
    taxonomy_rows: list[dict[str, object]] | None = None,
) -> tuple[Path, Path, Path]:
    blast = write_tsv(
        tmp_path / "sample-1.blast.final.tsv",
        BLAST_FIELDS,
        blast_rows
        if blast_rows is not None
        else [
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
        ],
    )
    coverage = write_tsv(
        tmp_path / "sample-1.crumbs.coverage.tsv",
        COVERAGE_FIELDS,
        coverage_rows
        if coverage_rows is not None
        else [
            {
                "sample_id": "sample-1",
                "qseqid": "contig-a",
                "contig_length": 100,
                "covered_bases_1x": 100,
                "breadth_1x": 1.0,
                "raw_aligned_bases": 100,
                "mean_depth_full": 1.0,
                "median_depth_full": 1,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 100,
                "crumbs_p99": 100,
            },
        ],
    )
    taxonomy = write_tsv(
        tmp_path / "profile_taxonomy.tsv",
        TAXONOMY_FIELDS,
        taxonomy_rows
        if taxonomy_rows is not None
        else [
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo|Homo sapiens",
                "rankpath": "no rank|superkingdom|genus|species",
            },
        ],
    )
    return blast, coverage, taxonomy


def run_estimator(tmp_path: Path, blast: Path, coverage: Path, taxonomy: Path) -> None:
    main(
        [
            "--sample-id",
            "sample-1",
            "--blast-tsv",
            str(blast),
            "--coverage-tsv",
            str(coverage),
            "--profile-taxonomy-tsv",
            str(taxonomy),
            "--output-dir",
            str(tmp_path),
        ],
    )


def taxon_metadata(taxon_ids: set[str]) -> pl.DataFrame:
    return pl.DataFrame(
        [
            {
                "taxon_id": taxon_id,
                "taxon_name": taxon_id,
                "rank": "species",
                "taxpath": f"root|{taxon_id}",
                "taxpathsn": f"root|{taxon_id}",
                "rankpath": "no rank|species",
            }
            for taxon_id in sorted(taxon_ids)
        ],
    )


def contig_evidence_row(
    *,
    qseqid: str,
    taxon_id: str,
    contig_length: int,
    crumbs_score: float,
    breadth_1x: float = 1.0,
) -> dict[str, object]:
    return {
        "sample_id": "sample-1",
        "qseqid": qseqid,
        "adjusted_taxid": taxon_id,
        "adjusted_taxid_name": taxon_id,
        "adjusted_taxid_rank": "species",
        "adjustment_method": "dominant",
        "contig_length": contig_length,
        "covered_bases_1x": round(contig_length * breadth_1x),
        "breadth_1x": breadth_1x,
        "raw_aligned_bases": crumbs_score,
        "mean_depth_full": crumbs_score / contig_length,
        "median_depth_full": 0.0,
        "median_depth_positive": 0.0,
        "depth_p95": 0.0,
        "depth_p99": 0.0,
        "max_depth": 0.0,
        "crumbs_p95": crumbs_score,
        "crumbs_p99": crumbs_score,
        "winsorization_percentile": "p99",
        "winsorization_cap": 0.0,
        "crumbs_score": crumbs_score,
    }


def test_estimates_taxon_crumbs_and_bioboxes_profile(tmp_path: Path) -> None:
    blast = write_tsv(
        tmp_path / "sample-1.blast.final.tsv",
        [
            "sample",
            "qseqid",
            "sseqid",
            "adjusted_taxid",
            "adjusted_taxid_name",
            "adjusted_taxid_rank",
            "adjustment_method",
        ],
        [
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "sseqid": "ref-a1",
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "sseqid": "ref-a2",
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
            {
                "sample": "sample-1",
                "qseqid": "contig-b",
                "sseqid": "ref-b1",
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
            {
                "sample": "sample-1",
                "qseqid": "contig-c",
                "sseqid": "ref-c1",
                "adjusted_taxid": 9598,
                "adjusted_taxid_name": "Pan troglodytes",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "lca",
            },
            {
                "sample": "sample-1",
                "qseqid": "contig-zero",
                "sseqid": "ref-z1",
                "adjusted_taxid": 9605,
                "adjusted_taxid_name": "Homo",
                "adjusted_taxid_rank": "genus",
                "adjustment_method": "lca",
            },
        ],
    )
    coverage = write_tsv(
        tmp_path / "sample-1.crumbs.coverage.tsv",
        [
            "sample_id",
            "qseqid",
            "contig_length",
            "covered_bases_1x",
            "breadth_1x",
            "raw_aligned_bases",
            "mean_depth_full",
            "median_depth_full",
            "median_depth_positive",
            "depth_p95",
            "depth_p99",
            "max_depth",
            "crumbs_p95",
            "crumbs_p99",
        ],
        [
            {
                "sample_id": "sample-1",
                "qseqid": "contig-a",
                "contig_length": 100,
                "covered_bases_1x": 50,
                "breadth_1x": 0.5,
                "raw_aligned_bases": 50,
                "mean_depth_full": 0.5,
                "median_depth_full": 0,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 50,
                "crumbs_p99": 50,
            },
            {
                "sample_id": "sample-1",
                "qseqid": "contig-b",
                "contig_length": 300,
                "covered_bases_1x": 150,
                "breadth_1x": 0.5,
                "raw_aligned_bases": 150,
                "mean_depth_full": 0.5,
                "median_depth_full": 0.5,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 150,
                "crumbs_p99": 150,
            },
            {
                "sample_id": "sample-1",
                "qseqid": "contig-c",
                "contig_length": 100,
                "covered_bases_1x": 100,
                "breadth_1x": 1.0,
                "raw_aligned_bases": 100,
                "mean_depth_full": 1.0,
                "median_depth_full": 1,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 100,
                "crumbs_p99": 100,
            },
            {
                "sample_id": "sample-1",
                "qseqid": "contig-zero",
                "contig_length": 100,
                "covered_bases_1x": 0,
                "breadth_1x": 0.0,
                "raw_aligned_bases": 0,
                "mean_depth_full": 0.0,
                "median_depth_full": 0,
                "median_depth_positive": 0,
                "depth_p95": 0,
                "depth_p99": 0,
                "max_depth": 0,
                "crumbs_p95": 0,
                "crumbs_p99": 0,
            },
        ],
    )
    taxonomy = write_tsv(
        tmp_path / "profile_taxonomy.tsv",
        TAXONOMY_FIELDS,
        [
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo|Homo sapiens",
                "rankpath": "no rank|superkingdom|genus|species",
            },
            {
                "taxon_id": 9598,
                "taxon_name": "Pan troglodytes",
                "rank": "species",
                "taxpath": "131567|2759|9596|9598",
                "taxpathsn": "cellular organisms|Eukaryota|Pan|Pan troglodytes",
                "rankpath": "no rank|superkingdom|genus|species",
            },
            {
                "taxon_id": 9605,
                "taxon_name": "Homo",
                "rank": "genus",
                "taxpath": "131567|2759|9605",
                "taxpathsn": "cellular organisms|Eukaryota|Homo",
                "rankpath": "no rank|superkingdom|genus",
            },
        ],
    )

    main(
        [
            "--sample-id",
            "sample-1",
            "--blast-tsv",
            str(blast),
            "--coverage-tsv",
            str(coverage),
            "--profile-taxonomy-tsv",
            str(taxonomy),
            "--output-dir",
            str(tmp_path),
        ],
    )

    contigs = read_tsv(tmp_path / "sample-1.crumbs.contigs.tsv")
    assert [row["qseqid"] for row in contigs] == [
        "contig-a",
        "contig-b",
        "contig-c",
        "contig-zero",
    ]
    assert float(contigs[0]["crumbs_score"]) == EXPECTED_CONTIG_A_CRUMBS

    taxa = read_tsv(tmp_path / "sample-1.crumbs.taxa.tsv")
    taxa_by_id = {row["taxon_id"]: row for row in taxa}
    assert taxa_by_id["9606"]["n_contigs"] == "2"
    assert taxa_by_id["9606"]["total_contig_length"] == "400"
    assert float(taxa_by_id["9606"]["total_crumbs_score"]) == EXPECTED_HOMO_TOTAL_CRUMBS
    assert float(taxa_by_id["9606"]["taxon_crumbs"]) == EXPECTED_HOMO_TAXON_CRUMBS
    assert float(taxa_by_id["9598"]["taxon_crumbs"]) == EXPECTED_PAN_TAXON_CRUMBS
    assert float(taxa_by_id["9605"]["taxon_crumbs"]) == EXPECTED_ZERO_TAXON_CRUMBS
    assert taxa_by_id["9605"]["n_zero_crumbs_contigs"] == "1"
    assert taxa_by_id["9606"]["rankpath"] == "no rank|superkingdom|genus|species"
    assert taxa_by_id["9605"]["rankpath"] == "no rank|superkingdom|genus"

    profile = bioboxes_data_rows(tmp_path / "sample-1.crumbs.bioboxes.profile.tsv")
    profile_by_id = {row["@@TAXID"]: row for row in profile}
    assert set(profile_by_id) == {"9606", "9598"}
    assert float(profile_by_id["9606"]["PERCENTAGE"]) == pytest.approx(
        EXPECTED_HOMO_PERCENTAGE,
    )
    assert float(profile_by_id["9598"]["PERCENTAGE"]) == pytest.approx(
        EXPECTED_PAN_PERCENTAGE,
    )
    profile_lines = read_lines(tmp_path / "sample-1.crumbs.bioboxes.profile.tsv")
    assert profile_lines[:3] == [
        "@SampleID:sample-1",
        "@Version:0.9.1",
        "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE",
    ]

    qc = json.loads((tmp_path / "sample-1.crumbs.qc.json").read_text(encoding="utf-8"))
    assert qc["sample_id"] == "sample-1"
    assert qc["n_contigs"] == EXPECTED_CONTIG_COUNT
    assert qc["n_bioboxes_taxa"] == EXPECTED_BIOBOXES_TAXA
    assert qc["crumbs_score_column"] == "crumbs_p99"


def test_missing_profile_taxonomy_row_fails_loudly(
    tmp_path: Path,
) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(tmp_path, taxonomy_rows=[])

    with pytest.raises(
        CrumbsProfileError,
        match="profile taxonomy is missing adjusted taxids: 9606",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_assigned_contig_missing_coverage_fails_loudly(
    tmp_path: Path,
) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        coverage_rows=[
            {
                "sample_id": "sample-1",
                "qseqid": "unassigned-contig",
                "contig_length": 100,
                "covered_bases_1x": 100,
                "breadth_1x": 1.0,
                "raw_aligned_bases": 100,
                "mean_depth_full": 1.0,
                "median_depth_full": 1,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 100,
                "crumbs_p99": 100,
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="coverage TSV is missing assigned contigs: contig-a",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_read_query_assignments_use_support_count_without_coverage_row(
    tmp_path: Path,
) -> None:
    blast = write_tsv(
        tmp_path / "sample-1.blast.final.tsv",
        [
            "sample",
            "qseqid",
            "qlen",
            "adjusted_taxid",
            "adjusted_taxid_name",
            "adjusted_taxid_rank",
            "adjustment_method",
            "evidence_class",
            "support_record_count",
        ],
        [
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "qlen": 100,
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
                "evidence_class": "short_assembly_contig",
                "support_record_count": 1,
            },
            {
                "sample": "sample-1",
                "qseqid": "nvdMergeReadQuery_sample-1_000001",
                "qlen": 50,
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
                "evidence_class": "overlap_merged_pair",
                "support_record_count": 7,
            },
        ],
    )
    coverage = write_tsv(
        tmp_path / "sample-1.crumbs.coverage.tsv",
        COVERAGE_FIELDS,
        [
            {
                "sample_id": "sample-1",
                "qseqid": "contig-a",
                "contig_length": 100,
                "covered_bases_1x": 100,
                "breadth_1x": 1.0,
                "raw_aligned_bases": 100,
                "mean_depth_full": 1.0,
                "median_depth_full": 1,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 100,
                "crumbs_p99": 100,
            },
        ],
    )
    taxonomy = write_tsv(
        tmp_path / "profile_taxonomy.tsv",
        TAXONOMY_FIELDS,
        [
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo|Homo sapiens",
                "rankpath": "no rank|superkingdom|genus|species",
            },
        ],
    )

    run_estimator(tmp_path, blast, coverage, taxonomy)

    contigs = read_tsv(tmp_path / "sample-1.crumbs.contigs.tsv")
    read_query = next(
        row for row in contigs if row["qseqid"] == "nvdMergeReadQuery_sample-1_000001"
    )
    assert read_query["contig_length"] == "50"
    assert read_query["covered_bases_1x"] == "0"
    assert read_query["breadth_1x"] == "0.0"
    assert read_query["raw_aligned_bases"] == "7.0"
    assert read_query["crumbs_score"] == "7.0"
    assert read_query["winsorization_percentile"] == "read_query_support"

    [taxon] = read_tsv(tmp_path / "sample-1.crumbs.taxa.tsv")
    assert taxon["n_contigs"] == "2"
    assert float(taxon["total_crumbs_score"]) == 107.0


def test_conflicting_repeated_blast_assignments_fail_loudly(
    tmp_path: Path,
) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        blast_rows=[
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo",
                "adjusted_taxid_rank": "genus",
                "adjustment_method": "dominant",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="conflicting assignments for contigs: contig-a",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_profile_assignment_uses_strongest_blast_row_when_tasks_disagree(
    tmp_path: Path,
) -> None:
    blast = write_tsv(
        tmp_path / "sample-1.blast.final.tsv",
        [
            "sample",
            "qseqid",
            "evalue",
            "bitscore",
            "adjusted_taxid",
            "adjusted_taxid_name",
            "adjusted_taxid_rank",
            "adjustment_method",
        ],
        [
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "evalue": "1.98e-9",
                "bitscore": "58.4",
                "adjusted_taxid": 10244,
                "adjusted_taxid_name": "Monkeypox virus",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "evalue": "8.09e-44",
                "bitscore": "172.0",
                "adjusted_taxid": 10258,
                "adjusted_taxid_name": "Orf virus",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "lca",
            },
        ],
    )
    coverage = write_tsv(
        tmp_path / "sample-1.crumbs.coverage.tsv",
        COVERAGE_FIELDS,
        [
            {
                "sample_id": "sample-1",
                "qseqid": "contig-a",
                "contig_length": 100,
                "covered_bases_1x": 100,
                "breadth_1x": 1.0,
                "raw_aligned_bases": 100,
                "mean_depth_full": 1.0,
                "median_depth_full": 1,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 100,
                "crumbs_p99": 100,
            },
        ],
    )
    taxonomy = write_tsv(
        tmp_path / "profile_taxonomy.tsv",
        TAXONOMY_FIELDS,
        [
            {
                "taxon_id": 10244,
                "taxon_name": "Monkeypox virus",
                "rank": "species",
                "taxpath": "10239|10244",
                "taxpathsn": "Viruses|Monkeypox virus",
                "rankpath": "superkingdom|species",
            },
            {
                "taxon_id": 10258,
                "taxon_name": "Orf virus",
                "rank": "species",
                "taxpath": "10239|10258",
                "taxpathsn": "Viruses|Orf virus",
                "rankpath": "superkingdom|species",
            },
        ],
    )

    run_estimator(tmp_path, blast, coverage, taxonomy)

    [contig] = read_tsv(tmp_path / "sample-1.crumbs.contigs.tsv")
    assert contig["adjusted_taxid"] == "10258"
    assert contig["adjusted_taxid_name"] == "Orf virus"
    assert contig["adjustment_method"] == "lca"


def test_equally_strong_conflicting_profile_assignments_fail_loudly(
    tmp_path: Path,
) -> None:
    _blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        taxonomy_rows=[
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo|Homo sapiens",
                "rankpath": "no rank|superkingdom|genus|species",
            },
            {
                "taxon_id": 9598,
                "taxon_name": "Pan troglodytes",
                "rank": "species",
                "taxpath": "131567|2759|9596|9598",
                "taxpathsn": "cellular organisms|Eukaryota|Pan|Pan troglodytes",
                "rankpath": "no rank|superkingdom|genus|species",
            },
        ],
    )
    blast = write_tsv(
        tmp_path / "sample-1.blast.final.tsv",
        [
            "sample",
            "qseqid",
            "evalue",
            "bitscore",
            "adjusted_taxid",
            "adjusted_taxid_name",
            "adjusted_taxid_rank",
            "adjustment_method",
        ],
        [
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "evalue": "1e-20",
                "bitscore": "200.0",
                "adjusted_taxid": 9606,
                "adjusted_taxid_name": "Homo sapiens",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
            {
                "sample": "sample-1",
                "qseqid": "contig-a",
                "evalue": "1e-20",
                "bitscore": "200.0",
                "adjusted_taxid": 9598,
                "adjusted_taxid_name": "Pan troglodytes",
                "adjusted_taxid_rank": "species",
                "adjustment_method": "dominant",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="conflicting assignments for contigs: contig-a",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_duplicate_coverage_rows_fail_loudly(
    tmp_path: Path,
) -> None:
    duplicate_coverage = [
        {
            "sample_id": "sample-1",
            "qseqid": "contig-a",
            "contig_length": 100,
            "covered_bases_1x": 100,
            "breadth_1x": 1.0,
            "raw_aligned_bases": 100,
            "mean_depth_full": 1.0,
            "median_depth_full": 1,
            "median_depth_positive": 1,
            "depth_p95": 1,
            "depth_p99": 1,
            "max_depth": 1,
            "crumbs_p95": 100,
            "crumbs_p99": 100,
        },
        {
            "sample_id": "sample-1",
            "qseqid": "contig-a",
            "contig_length": 100,
            "covered_bases_1x": 100,
            "breadth_1x": 1.0,
            "raw_aligned_bases": 100,
            "mean_depth_full": 1.0,
            "median_depth_full": 1,
            "median_depth_positive": 1,
            "depth_p95": 1,
            "depth_p99": 1,
            "max_depth": 1,
            "crumbs_p95": 100,
            "crumbs_p99": 100,
        },
    ]
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        coverage_rows=duplicate_coverage,
    )

    with pytest.raises(
        CrumbsProfileError,
        match="duplicate qseqid rows for sample sample-1: contig-a",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_null_numeric_coverage_values_fail_loudly(tmp_path: Path) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        coverage_rows=[
            {
                "sample_id": "sample-1",
                "qseqid": "contig-a",
                "contig_length": "",
                "covered_bases_1x": 100,
                "breadth_1x": 1.0,
                "raw_aligned_bases": 100,
                "mean_depth_full": 1.0,
                "median_depth_full": 1,
                "median_depth_positive": 1,
                "depth_p95": 1,
                "depth_p99": 1,
                "max_depth": 1,
                "crumbs_p95": 100,
                "crumbs_p99": 100,
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="invalid coverage values for contigs: contig-a",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


@pytest.mark.parametrize(
    "bad_updates",
    [
        {"covered_bases_1x": 101},
        {"breadth_1x": 0.25},
        {"mean_depth_full": 0.25},
        {"depth_p95": 2},
        {"depth_p99": 2},
        {"crumbs_p99": 101},
        {"crumbs_p99": 100, "depth_p99": 0.5},
    ],
)
def test_impossible_coverage_relationships_fail_loudly(
    tmp_path: Path,
    bad_updates: dict[str, object],
) -> None:
    coverage_row = {
        "sample_id": "sample-1",
        "qseqid": "contig-a",
        "contig_length": 100,
        "covered_bases_1x": 100,
        "breadth_1x": 1.0,
        "raw_aligned_bases": 100,
        "mean_depth_full": 1.0,
        "median_depth_full": 1,
        "median_depth_positive": 1,
        "depth_p95": 1,
        "depth_p99": 1,
        "max_depth": 1,
        "crumbs_p95": 100,
        "crumbs_p99": 100,
    } | bad_updates
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        coverage_rows=[coverage_row],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="invalid coverage values for contigs: contig-a",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_profile_taxonomy_path_contract_is_validated(
    tmp_path: Path,
) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        taxonomy_rows=[
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605",
                "taxpathsn": "cellular organisms|Eukaryota|Homo|Homo sapiens",
                "rankpath": "no rank|superkingdom|genus|species",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="path length mismatch for taxon_id 9606",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_profile_taxonomy_taxpath_must_end_with_taxon_id(tmp_path: Path) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        taxonomy_rows=[
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605",
                "taxpathsn": "cellular organisms|Eukaryota|Homo sapiens",
                "rankpath": "no rank|superkingdom|species",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="taxpath for taxon_id 9606 does not end with taxon_id",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_profile_taxonomy_taxpathsn_must_end_with_taxon_name(tmp_path: Path) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        taxonomy_rows=[
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo",
                "rankpath": "no rank|superkingdom|species",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="taxpathsn for taxon_id 9606 does not end with taxon_name",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_profile_taxonomy_rankpath_length_must_match_taxpath(tmp_path: Path) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        taxonomy_rows=[
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo|Homo sapiens",
                "rankpath": "no rank|superkingdom|species",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="path length mismatch for taxon_id 9606",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_profile_taxonomy_rankpath_must_end_with_terminal_rank(
    tmp_path: Path,
) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        taxonomy_rows=[
            {
                "taxon_id": 9606,
                "taxon_name": "Homo sapiens",
                "rank": "species",
                "taxpath": "131567|2759|9605|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo|Homo sapiens",
                "rankpath": "no rank|superkingdom|genus|genus",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="rankpath for taxon_id 9606 does not end with rank",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_profile_taxonomy_blank_required_values_fail(tmp_path: Path) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        taxonomy_rows=[
            {
                "taxon_id": 9606,
                "taxon_name": "",
                "rank": "species",
                "taxpath": "131567|2759|9606",
                "taxpathsn": "cellular organisms|Eukaryota|Homo sapiens",
                "rankpath": "no rank|superkingdom|species",
            },
        ],
    )

    with pytest.raises(
        CrumbsProfileError,
        match="profile taxonomy contains null or blank required values",
    ):
        run_estimator(tmp_path, blast, coverage, taxonomy)


def test_all_zero_crumbs_sample_keeps_sidecars_without_bioboxes_taxa(
    tmp_path: Path,
) -> None:
    blast, coverage, taxonomy = write_minimal_inputs(
        tmp_path,
        coverage_rows=[
            {
                "sample_id": "sample-1",
                "qseqid": "contig-a",
                "contig_length": 100,
                "covered_bases_1x": 0,
                "breadth_1x": 0.0,
                "raw_aligned_bases": 0,
                "mean_depth_full": 0.0,
                "median_depth_full": 0,
                "median_depth_positive": 0,
                "depth_p95": 0,
                "depth_p99": 0,
                "max_depth": 0,
                "crumbs_p95": 0,
                "crumbs_p99": 0,
            },
        ],
    )

    run_estimator(tmp_path, blast, coverage, taxonomy)

    taxa = read_tsv(tmp_path / "sample-1.crumbs.taxa.tsv")
    assert taxa[0]["taxon_id"] == "9606"
    assert float(taxa[0]["taxon_crumbs"]) == 0.0
    assert bioboxes_data_rows(tmp_path / "sample-1.crumbs.bioboxes.profile.tsv") == []
    qc = json.loads((tmp_path / "sample-1.crumbs.qc.json").read_text(encoding="utf-8"))
    assert qc["n_bioboxes_taxa"] == 0


@settings(max_examples=100)
@given(
    rows=st.lists(
        st.tuples(
            st.sampled_from(["tax-a", "tax-b", "tax-c"]),
            st.integers(min_value=1, max_value=10_000),
            st.floats(
                min_value=0,
                max_value=1_000_000,
                allow_nan=False,
                allow_infinity=False,
            ),
        ),
        min_size=1,
        max_size=30,
    ),
)
def test_taxon_table_obeys_independent_crumbs_formula(
    rows: list[tuple[str, int, float]],
) -> None:
    contigs = pl.DataFrame(
        [
            contig_evidence_row(
                qseqid=f"contig-{index}",
                taxon_id=taxon_id,
                contig_length=contig_length,
                crumbs_score=crumbs_score,
            )
            for index, (taxon_id, contig_length, crumbs_score) in enumerate(rows)
        ],
    )
    taxonomy = taxon_metadata({taxon_id for taxon_id, _length, _crumbs in rows})

    taxa = build_taxon_table(
        sample_id="sample-1",
        contigs=contigs,
        profile_taxonomy=taxonomy,
    )

    expected: dict[str, dict[str, float]] = {}
    for taxon_id, contig_length, crumbs_score in rows:
        bucket = expected.setdefault(taxon_id, {"length": 0.0, "crumbs": 0.0})
        bucket["length"] += contig_length
        bucket["crumbs"] += crumbs_score

    expected_taxon_crumbs = {
        taxon_id: values["crumbs"] / values["length"]
        for taxon_id, values in expected.items()
    }
    percentage_denominator = sum(expected_taxon_crumbs.values())

    for row in taxa.iter_rows(named=True):
        taxon_id = row["taxon_id"]
        assert row["taxon_crumbs"] == pytest.approx(expected_taxon_crumbs[taxon_id])
        expected_percentage = (
            expected_taxon_crumbs[taxon_id] / percentage_denominator * 100
            if percentage_denominator > 0
            else 0.0
        )
        assert row["percentage_emitted"] == pytest.approx(expected_percentage)

    expected_sum = 100.0 if percentage_denominator > 0 else 0.0
    assert taxa.get_column("percentage_emitted").sum() == pytest.approx(expected_sum)


@settings(max_examples=100)
@given(
    depth_a=st.floats(
        min_value=0.001,
        max_value=1_000,
        allow_nan=False,
        allow_infinity=False,
    ),
    depth_b=st.floats(
        min_value=0.001,
        max_value=1_000,
        allow_nan=False,
        allow_infinity=False,
    ),
    taxon_a_lengths=st.lists(
        st.integers(min_value=1, max_value=10_000),
        min_size=1,
        max_size=20,
    ),
    taxon_b_length=st.integers(min_value=1, max_value=10_000),
)
def test_taxon_percentage_is_invariant_to_same_depth_fragmentation(
    depth_a: float,
    depth_b: float,
    taxon_a_lengths: list[int],
    taxon_b_length: int,
) -> None:
    contig_rows = [
        contig_evidence_row(
            qseqid=f"tax-a-contig-{index}",
            taxon_id="tax-a",
            contig_length=contig_length,
            crumbs_score=depth_a * contig_length,
        )
        for index, contig_length in enumerate(taxon_a_lengths)
    ]
    contig_rows.append(
        contig_evidence_row(
            qseqid="tax-b-contig-0",
            taxon_id="tax-b",
            contig_length=taxon_b_length,
            crumbs_score=depth_b * taxon_b_length,
        ),
    )
    taxa = build_taxon_table(
        sample_id="sample-1",
        contigs=pl.DataFrame(contig_rows),
        profile_taxonomy=taxon_metadata({"tax-a", "tax-b"}),
    )

    assert taxa.height == EXPECTED_FRAGMENTED_TAXON_COUNT
    percentages = {
        row["taxon_id"]: row["percentage_emitted"] for row in taxa.iter_rows(named=True)
    }
    assert percentages["tax-a"] == pytest.approx(depth_a / (depth_a + depth_b) * 100)
    assert percentages["tax-b"] == pytest.approx(depth_b / (depth_a + depth_b) * 100)
