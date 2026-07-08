"""Tests for taxon-level Big Table construction."""

# ruff: noqa: COM812

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

from build_taxon_big_table import main

if TYPE_CHECKING:
    from pathlib import Path


QUERY_COLUMNS = [
    "sample_id",
    "assigned_taxid_name",
    "assigned_taxid_rank",
    "support_tier",
    "query_class",
    "crumbs_score",
    "qlen",
    "assigned_taxid",
    "assignment_method",
    "qseqid",
    "best_hit_qcov",
    "best_hit_pident",
    "best_hit_evalue",
    "best_hit_bitscore",
    "retained_reference_count",
    "assignment_reference_count",
    "assignment_taxid_count",
    "support_tier_rule",
    "support_note",
    "support_record_count",
    "mapped_reads",
    "producer",
    "source_id",
    "best_hit_alignment_length",
    "blast_db_version",
    "virus_index_version",
    "nextflow_run_id",
]
CRUMBS_TAXA_COLUMNS = [
    "sample_id",
    "taxon_id",
    "taxon_name",
    "rank",
    "total_contig_length",
    "total_crumbs_score",
    "taxon_crumbs",
    "percentage_emitted",
]
EXPECTED_LEFT_TO_RIGHT_COLUMNS = [
    "sample_id",
    "taxon_name",
    "taxon_rank",
    "support_tier",
    "taxon_crumbs",
    "relative_crumbs_percent",
    "supporting_query_count",
    "taxid",
    "total_query_span",
    "total_crumbs_score",
    "strong_query_count",
    "moderate_query_count",
    "weak_query_count",
    "review_query_count",
    "redacted_query_count",
    "supporting_genome_like_contig_count",
    "supporting_long_contig_count",
    "supporting_short_contig_count",
    "supporting_merged_pair_count",
    "supporting_single_read_count",
    "support_tier_rule",
    "support_note",
]


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def query_row(**overrides: object) -> dict[str, object]:
    row = {
        "sample_id": "sample-1",
        "assigned_taxid_name": "Alpha virus",
        "assigned_taxid_rank": "species",
        "support_tier": "strong",
        "query_class": "long_assembly_contig",
        "crumbs_score": "100",
        "qlen": "1000",
        "assigned_taxid": "111",
        "assignment_method": "dominant",
        "qseqid": "q1",
        "best_hit_qcov": "0.95",
        "best_hit_pident": "99",
        "best_hit_evalue": "1e-50",
        "best_hit_bitscore": "100",
        "retained_reference_count": "1",
        "assignment_reference_count": "1",
        "assignment_taxid_count": "1",
        "support_tier_rule": "long_contig_dominant_high_qcov",
        "support_note": "strong query",
        "support_record_count": "1",
        "mapped_reads": "10",
        "producer": "spades",
        "source_id": "NODE_1",
        "best_hit_alignment_length": "950",
        "blast_db_version": "nt-test",
        "virus_index_version": "virus-test",
        "nextflow_run_id": "run-test",
    }
    row.update(overrides)
    return row


def test_taxon_big_table_aggregates_query_support_and_crumbs(tmp_path: Path) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    crumbs_taxa = tmp_path / "sample-1.crumbs.taxa.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        query_big_table,
        [
            query_row(qseqid="q1", crumbs_score="100", qlen="1000"),
            query_row(qseqid="q2", crumbs_score="200", qlen="2000"),
        ],
        QUERY_COLUMNS,
    )
    write_tsv(
        crumbs_taxa,
        [
            {
                "sample_id": "sample-1",
                "taxon_id": "111",
                "taxon_name": "Alpha virus",
                "rank": "species",
                "total_contig_length": "3000",
                "total_crumbs_score": "300",
                "taxon_crumbs": "0.1",
                "percentage_emitted": "100",
            },
        ],
        CRUMBS_TAXA_COLUMNS,
    )

    main(
        [
            "--query-big-table",
            str(query_big_table),
            "--crumbs-taxa-tsv",
            str(crumbs_taxa),
            "--output",
            str(output),
        ]
    )

    [row] = read_tsv(output)
    assert list(row) == EXPECTED_LEFT_TO_RIGHT_COLUMNS
    assert row["taxid"] == "111"
    assert row["taxon_name"] == "Alpha virus"
    assert row["support_tier"] == "strong"
    assert row["support_tier_rule"] == "multi_query_strong_support"
    assert row["supporting_query_count"] == "2"
    assert row["supporting_long_contig_count"] == "2"
    assert row["total_query_span"] == "3000"
    assert row["total_crumbs_score"] == "300"
    assert row["taxon_crumbs"] == "0.1"
    assert row["relative_crumbs_percent"] == "100"
    assert row["support_note"] == (
        "2 strong query assignments support Alpha virus (species)."
    )


def test_single_strong_single_read_does_not_become_strong_taxon_support(
    tmp_path: Path,
) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        query_big_table,
        [
            query_row(
                qseqid="read-q1",
                query_class="single_read",
                qlen="150",
                crumbs_score="150",
                support_tier="strong",
                support_tier_rule="single_read_dominant_high_qcov",
            ),
        ],
        QUERY_COLUMNS,
    )

    main(["--query-big-table", str(query_big_table), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["support_tier"] == "moderate"
    assert row["support_tier_rule"] == "single_read_strong_query_lacks_corroboration"
    assert row["supporting_query_count"] == "1"
    assert row["supporting_single_read_count"] == "1"
    assert row["support_note"] == (
        "1 strong read-derived query assignment supports Alpha virus (species); "
        "no additional query is assigned to this taxon."
    )


def test_single_strong_short_contig_does_not_become_strong_taxon_support(
    tmp_path: Path,
) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        query_big_table,
        [
            query_row(
                qseqid="short-contig-q1",
                query_class="short_assembly_contig",
                qlen="450",
                crumbs_score="450",
                support_tier="strong",
                support_tier_rule="short_contig_dominant_high_qcov",
            ),
        ],
        QUERY_COLUMNS,
    )

    main(["--query-big-table", str(query_big_table), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["support_tier"] == "moderate"
    assert (
        row["support_tier_rule"]
        == "single_short_contig_strong_query_lacks_corroboration"
    )
    assert row["supporting_query_count"] == "1"
    assert row["supporting_short_contig_count"] == "1"
    assert row["support_note"] == (
        "1 strong short-contig query assignment supports Alpha virus (species); "
        "no additional query is assigned to this taxon."
    )


def test_single_strong_long_contig_has_rank_safe_note(tmp_path: Path) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(query_big_table, [query_row()], QUERY_COLUMNS)

    main(["--query-big-table", str(query_big_table), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["support_tier"] == "strong"
    assert row["support_note"] == (
        "1 strong long-contig query assignment supports Alpha virus (species)."
    )


def test_weak_query_does_not_downgrade_strong_long_contig_support(
    tmp_path: Path,
) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        query_big_table,
        [
            query_row(qseqid="strong-long"),
            query_row(
                qseqid="weak-read",
                query_class="single_read",
                support_tier="weak",
                qlen="150",
            ),
        ],
        QUERY_COLUMNS,
    )

    main(["--query-big-table", str(query_big_table), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["support_tier"] == "strong"
    assert row["support_tier_rule"] == "strong_long_or_genome_like_query"
    assert row["support_note"] == (
        "1 strong long-contig query assignment supports Alpha virus (species); "
        "1 additional query is assigned to this taxon."
    )


def test_all_weak_queries_have_rank_safe_note(tmp_path: Path) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        query_big_table,
        [
            query_row(qseqid="weak-1", support_tier="weak"),
            query_row(qseqid="weak-2", support_tier="weak"),
        ],
        QUERY_COLUMNS,
    )

    main(["--query-big-table", str(query_big_table), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["support_tier"] == "weak"
    assert row["support_tier_rule"] == "all_weak_queries"
    assert row["support_note"] == (
        "All 2 query assignments to Alpha virus (species) have low best-hit "
        "query coverage."
    )


def test_mixed_support_note_reports_strong_and_additional_queries(
    tmp_path: Path,
) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        query_big_table,
        [
            query_row(
                qseqid="strong-short",
                query_class="short_assembly_contig",
                support_tier="strong",
            ),
            query_row(
                qseqid="weak-read",
                query_class="single_read",
                support_tier="weak",
            ),
        ],
        QUERY_COLUMNS,
    )

    main(["--query-big-table", str(query_big_table), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["support_tier"] == "moderate"
    assert row["support_tier_rule"] == "mixed_support_with_strong_query"
    assert row["support_note"] == (
        "1 strong and 1 additional query assignment support Alpha virus (species)."
    )


def test_lca_queries_aggregate_under_the_assigned_ancestor(tmp_path: Path) -> None:
    query_big_table = tmp_path / "query_big_table.tsv"
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        query_big_table,
        [
            query_row(
                qseqid="lca-1",
                assigned_taxid="10239",
                assigned_taxid_name="Viruses",
                assigned_taxid_rank="superkingdom",
                assignment_method="lca",
            ),
            query_row(
                qseqid="lca-2",
                assigned_taxid="10239",
                assigned_taxid_name="Viruses",
                assigned_taxid_rank="superkingdom",
                assignment_method="lca",
            ),
        ],
        QUERY_COLUMNS,
    )

    main(["--query-big-table", str(query_big_table), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["taxid"] == "10239"
    assert row["taxon_name"] == "Viruses"
    assert row["taxon_rank"] == "superkingdom"
