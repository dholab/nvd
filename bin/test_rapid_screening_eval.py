"""Tests for the rapid-screening eval."""

from __future__ import annotations

import csv
from pathlib import Path

import duckdb
import pytest
from evaluate_rapid_screening import build_outputs
from render_rapid_screening_eval import render_html

ROOT = Path(__file__).resolve().parents[1]


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_build_outputs_classifies_followup_reference_only_and_no_followup_signals(
    tmp_path: Path,
) -> None:
    sample_counts = tmp_path / "sample_read_counts.tsv"
    sample_counts.write_text("sample-1\t1000\n", encoding="utf-8")
    lineages = write_csv(
        tmp_path / "lineages.csv",
        [
            "ident",
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
        [
            {
                "ident": "NC_1",
                "superkingdom": "Viruses",
                "phylum": "Pisuviricota",
                "class": "Pisoniviricetes",
                "order": "Picornavirales",
                "family": "Picornaviridae",
                "genus": "Enterovirus",
                "species": "Enterovirus C",
            },
            {
                "ident": "NC_2",
                "superkingdom": "Viruses",
                "phylum": "Negarnaviricota",
                "class": "Insthoviricetes",
                "order": "Articulavirales",
                "family": "Orthomyxoviridae",
                "genus": "Alphainfluenzavirus",
                "species": "Influenza A virus",
            },
            {
                "ident": "WVDB|contig-a",
                "superkingdom": "Viruses",
                "phylum": "WVDBphylum",
                "class": "WVDBclass",
                "order": "WVDBorder",
                "family": "WVDBfamily",
                "genus": "WVDBgenus",
                "species": "WVDBspecies",
            },
        ],
    )
    tax_summary = write_csv(
        tmp_path / "sample-1.sourmash.tax_metagenome.summarized.csv",
        ["rank", "lineage", "fraction", "f_weighted_at_rank"],
        [
            {
                "rank": "family",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae",
                "fraction": "0.2",
                "f_weighted_at_rank": "0.2",
            },
            {
                "rank": "species",
                "lineage": "Viruses;Negarnaviricota;Insthoviricetes;Articulavirales;Orthomyxoviridae;Alphainfluenzavirus;Influenza A virus",
                "fraction": "0.0000000001",
                "f_weighted_at_rank": "0.0000000001",
            },
            {
                "rank": "species",
                "lineage": "Viruses;WVDBphylum;WVDBclass;WVDBorder;WVDBfamily;WVDBgenus;WVDBspecies",
                "fraction": "0.3",
                "f_weighted_at_rank": "0.3",
            },
            {
                "rank": "species",
                "lineage": "Viruses;UnknownPhylum;UnknownClass;UnknownOrder;UnknownFamily;UnknownGenus;Unknown species",
                "fraction": "0.1",
                "f_weighted_at_rank": "0.1",
            },
        ],
    )
    gather = write_csv(
        tmp_path / "sample-1.sourmash.gather.csv",
        ["name", "sum_weighted_found", "total_weighted_hashes"],
        [{"name": "NC_1", "sum_weighted_found": "80", "total_weighted_hashes": "100"}],
    )
    crumbs_taxa = write_tsv(
        tmp_path / "sample-1.crumbs.taxa.tsv",
        [
            "sample_id",
            "taxon_id",
            "taxon_name",
            "rank",
            "taxpath",
            "taxpathsn",
            "rankpath",
            "taxon_crumbs",
            "percentage_emitted",
        ],
        [
            {
                "sample_id": "sample-1",
                "taxon_id": "12059",
                "taxon_name": "Enterovirus C",
                "rank": "species",
                "taxpath": "1|2|3|4|5|6|12059",
                "taxpathsn": "Viruses|Pisuviricota|Pisoniviricetes|Picornavirales|Picornaviridae|Enterovirus|Enterovirus C",
                "rankpath": "superkingdom|phylum|class|order|family|genus|species",
                "taxon_crumbs": "2.5",
                "percentage_emitted": "100",
            },
        ],
    )

    build_outputs(
        sample_read_counts=sample_counts,
        lineages_csv=lineages,
        sourmash_tax_summaries=[tax_summary],
        sourmash_gather_csvs=[gather],
        crumbs_taxa_tsvs=[crumbs_taxa],
        output_database=tmp_path / "rapid_screening_eval.duckdb",
        output_dir=tmp_path,
    )

    followup = read_tsv(
        tmp_path / "exports" / "screening_signal_followup_by_sample_rank.tsv",
    )
    assert {(row["comparison_rank"], row["followup_bin"]) for row in followup} == {
        ("family", "same_rank_followup_evidence"),
        ("species", "screening_reference_only"),
        ("species", "no_same_rank_followup_evidence"),
        ("species", "not_comparable"),
    }
    assert {row["screening_engine"] for row in followup} == {"sourmash"}
    assert list(followup[0]) == [
        "sample_id",
        "screening_engine",
        "comparison_rank",
        "followup_bin",
        "n_signals",
        "total_primary_screening_mass",
    ]

    without_followup = read_tsv(
        tmp_path / "exports" / "screening_signals_without_same_rank_followup.tsv",
    )
    assert len(without_followup) == 1
    assert without_followup[0]["screening_engine"] == "sourmash"
    assert without_followup[0]["screening_taxon_name"] == "Influenza A virus"
    assert float(without_followup[0]["primary_screening_mass"]) == pytest.approx(1e-10)
    assert float(without_followup[0]["sourmash_fraction"]) == pytest.approx(1e-10)
    assert float(without_followup[0]["sourmash_f_weighted_at_rank"]) == pytest.approx(
        1e-10,
    )
    assert list(without_followup[0]) == [
        "sample_id",
        "screening_engine",
        "comparison_rank",
        "screening_taxon_name",
        "screening_lineage",
        "primary_screening_mass",
        "sourmash_fraction",
        "sourmash_f_weighted_at_rank",
        "lineage_source",
        "source_path",
    ]

    con = duckdb.connect(str(tmp_path / "rapid_screening_eval.duckdb"), read_only=True)
    try:
        sample = con.execute("select * from samples").fetchone()
        assert sample[0] == "sample-1"
        assert sample[4] == pytest.approx(0.8)
        assert sample[5] == pytest.approx(0.2)
        tables = {row[0] for row in con.execute("show tables").fetchall()}
        assert {
            "screening_detections",
            "screening_signal_followup_by_rank",
            "screening_signals_without_same_rank_followup",
        } <= tables
        assert {
            "sourmash_detections",
            "sourmash_followup_by_rank",
            "unsupported_sourmash_signals",
        }.isdisjoint(tables)
        detection_columns = [
            row[1]
            for row in con.execute(
                "pragma table_info('screening_detections')",
            ).fetchall()
        ]
        assert detection_columns == [
            "sample_id",
            "screening_engine",
            "source_path",
            "rank",
            "screening_taxon_name",
            "screening_lineage",
            "primary_screening_mass",
            "sourmash_fraction",
            "sourmash_f_weighted_at_rank",
            "lineage_source",
            "followup_bin",
        ]
        screening_detection = con.execute(
            """
            select screening_engine,
                   screening_taxon_name,
                   sourmash_fraction,
                   sourmash_f_weighted_at_rank
            from screening_detections
            where screening_taxon_name = 'Picornaviridae'
            """,
        ).fetchone()
        assert screening_detection == (
            "sourmash",
            "Picornaviridae",
            pytest.approx(0.2),
            pytest.approx(0.2),
        )
    finally:
        con.close()

    report = render_html(
        database=tmp_path / "rapid_screening_eval.duckdb",
        followup_tsv=tmp_path
        / "exports"
        / "screening_signal_followup_by_sample_rank.tsv",
        without_followup_tsv=tmp_path
        / "exports"
        / "screening_signals_without_same_rank_followup.tsv",
    )
    assert "Rapid screening eval" in report
    assert "Influenza A virus" in report


def test_shared_wvdb_and_non_wvdb_lineage_is_not_reference_only(tmp_path: Path) -> None:
    sample_counts = tmp_path / "sample_read_counts.tsv"
    sample_counts.write_text("sample-1\t1000\n", encoding="utf-8")
    lineages = write_csv(
        tmp_path / "lineages.csv",
        [
            "ident",
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
        [
            {
                "ident": "WVDB|contig-a",
                "superkingdom": "Viruses",
                "phylum": "Pisuviricota",
                "class": "Pisoniviricetes",
                "order": "Picornavirales",
                "family": "Sharedfamily",
                "genus": "Sharedgenus",
                "species": "Shared species",
            },
            {
                "ident": "NC_1",
                "superkingdom": "Viruses",
                "phylum": "Pisuviricota",
                "class": "Pisoniviricetes",
                "order": "Picornavirales",
                "family": "Sharedfamily",
                "genus": "Sharedgenus",
                "species": "Shared species",
            },
        ],
    )
    tax_summary = write_csv(
        tmp_path / "sample-1.sourmash.tax_metagenome.summarized.csv",
        ["rank", "lineage", "fraction", "f_weighted_at_rank"],
        [
            {
                "rank": "species",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;Sharedfamily;Sharedgenus;Shared species",
                "fraction": "0.25",
                "f_weighted_at_rank": "0.25",
            },
        ],
    )

    build_outputs(
        sample_read_counts=sample_counts,
        lineages_csv=lineages,
        sourmash_tax_summaries=[tax_summary],
        sourmash_gather_csvs=[],
        crumbs_taxa_tsvs=[],
        output_database=tmp_path / "rapid_screening_eval.duckdb",
        output_dir=tmp_path,
    )

    without_followup = read_tsv(
        tmp_path / "exports" / "screening_signals_without_same_rank_followup.tsv",
    )
    assert len(without_followup) == 1
    assert without_followup[0]["lineage_source"] == "non_wvdb+wvdb"


def test_rapid_screening_eval_wiring_and_publish_hierarchy() -> None:
    """Rapid screening should own the stage and eval while sourmash stays the engine."""
    source_text = {
        "workflows/nvd_main.nf": (ROOT / "workflows" / "nvd_main.nf").read_text(
            encoding="utf-8",
        ),
        "subworkflows/rapid_screening.nf": (
            ROOT / "subworkflows" / "rapid_screening.nf"
        ).read_text(encoding="utf-8"),
        "subworkflows/rapid_screening_eval.nf": (
            ROOT / "subworkflows" / "rapid_screening_eval.nf"
        ).read_text(encoding="utf-8"),
        "modules/rapid_screening_eval.nf": (
            ROOT / "modules" / "rapid_screening_eval.nf"
        ).read_text(encoding="utf-8"),
        "conf/results.config": (ROOT / "conf" / "results.config").read_text(
            encoding="utf-8",
        ),
    }
    expected_by_source = {
        "workflows/nvd_main.nf": (
            'include { RAPID_SCREENING         } from "../subworkflows/rapid_screening"',
            'include { RAPID_SCREENING_EVAL    } from "../subworkflows/rapid_screening_eval"',
            "rapid_screening = RAPID_SCREENING(PREPROCESS_READS.out.reads)",
            "RAPID_SCREENING_EVAL(",
        ),
        "subworkflows/rapid_screening.nf": ("workflow RAPID_SCREENING",),
        "subworkflows/rapid_screening_eval.nf": (
            "workflow RAPID_SCREENING_EVAL",
            "BUILD_RAPID_SCREENING_EVAL_DB(",
            "RENDER_RAPID_SCREENING_EVAL(",
        ),
        "modules/rapid_screening_eval.nf": (
            "process BUILD_RAPID_SCREENING_EVAL_DB",
            "process RENDER_RAPID_SCREENING_EVAL",
            'path "rapid_screening_eval.duckdb"',
            'path "exports/screening_signal_followup_by_sample_rank.tsv"',
            'path "exports/screening_signals_without_same_rank_followup.tsv"',
            'path "rapid_screening_eval.html"',
        ),
        "conf/results.config": (
            'experimental_rapid_screening = params.nvd_results + "/experimental_rapid_screening"',
            'sourmash_screening_engine = params.experimental_rapid_screening + "/engines/sourmash"',
            'path: { params.experimental_rapid_screening + "/eval" }',
            'path: { params.experimental_rapid_screening + "/eval/reports" }',
        ),
    }
    missing = [
        f"{source}: missing {fragment}"
        for source, fragments in expected_by_source.items()
        for fragment in fragments
        if fragment not in source_text[source]
    ]

    assert not missing, "\n".join(missing)
    for retired_path in (
        ROOT / "subworkflows" / "metagenome_profiling.nf",
        ROOT / "subworkflows" / "sourmash_evaluation_harness.nf",
        ROOT / "modules" / "sourmash_evaluation.nf",
    ):
        assert not retired_path.exists()
