"""Tests for query-level Big Table construction."""

# ruff: noqa: COM812

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

import pytest
from build_query_big_table import main
from py_nvd.blast_consensus import ASSIGNMENT_BITSCORE_FRACTION

if TYPE_CHECKING:
    from pathlib import Path


BLAST_COLUMNS = [
    "task",
    "sample",
    "qseqid",
    "qlen",
    "sseqid",
    "stitle",
    "length",
    "pident",
    "evalue",
    "bitscore",
    "sscinames",
    "staxids",
    "rank",
    "adjusted_taxid",
    "adjusted_taxid_name",
    "adjusted_taxid_rank",
    "adjustment_method",
    "query_class",
    "producer",
    "source_id",
    "support_record_count",
    "mapped_reads",
    "total_reads",
    "blast_db_version",
    "virus_index_version",
    "nextflow_run_id",
]
EXPECTED_LEFT_TO_RIGHT_COLUMNS = [
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


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def base_hit(**overrides: str) -> dict[str, str]:
    row = {
        "task": "megablast",
        "sample": "sample-1",
        "qseqid": "nvdContigQuery_sample-1_000001",
        "qlen": "200",
        "sseqid": "ref-alpha-1",
        "stitle": "Alpha virus reference",
        "length": "180",
        "pident": "99.0",
        "evalue": "1e-50",
        "bitscore": "100",
        "sscinames": "Alpha virus",
        "staxids": "111",
        "rank": "acellular root:Viruses; family:Alpha; genus:Alphavirus; species:Alpha virus",
        "adjusted_taxid": "111",
        "adjusted_taxid_name": "Alpha virus",
        "adjusted_taxid_rank": "species",
        "adjustment_method": "dominant",
        "query_class": "long_assembly_contig",
        "producer": "spades",
        "source_id": "NODE_1",
        "support_record_count": "1",
        "mapped_reads": "12",
        "total_reads": "1000",
        "blast_db_version": "nt-test",
        "virus_index_version": "virus-test",
        "nextflow_run_id": "run-test",
    }
    row.update(overrides)
    return row


def test_query_big_table_collapses_retained_hits_to_one_auditable_assignment_row(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    crumbs_tsv = tmp_path / "crumbs.tsv"
    output = tmp_path / "query_big_table.tsv"
    rows = [
        base_hit(
            sseqid="ref-alpha-2",
            stitle="Second Alpha reference",
            evalue="1e-20",
            bitscore="70",
            staxids="112",
            rank="acellular root:Viruses; family:Alpha; genus:Betavirus; species:Beta virus",
        ),
        base_hit(),
    ]
    write_tsv(blast_tsv, rows, BLAST_COLUMNS)
    write_tsv(
        crumbs_tsv,
        [{"sample_id": "sample-1", "qseqid": rows[0]["qseqid"], "crumbs_score": "500"}],
        ["sample_id", "qseqid", "crumbs_score"],
    )

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
            "--crumbs-tsv",
            str(crumbs_tsv),
            "--output",
            str(output),
        ]
    )

    [row] = read_tsv(output)
    assert list(row) == EXPECTED_LEFT_TO_RIGHT_COLUMNS
    assert row["sample_id"] == "sample-1"
    assert row["qseqid"] == "nvdContigQuery_sample-1_000001"
    assert row["assigned_taxid"] == "111"
    assert row["assigned_taxid_name"] == "Alpha virus"
    assert row["assignment_method"] == "dominant"
    assert row["best_hit_qcov"] == "0.9"
    assert row["best_hit_evalue"] == "1e-50"
    assert row["best_hit_alignment_length"] == "180"
    assert row["retained_reference_count"] == "2"
    assert row["assignment_reference_count"] == "1"
    assert row["assignment_taxid_count"] == "1"
    assert row["crumbs_score"] == "500"
    assert row["support_tier"] == "strong"
    assert row["support_tier_rule"] == "long_contig_dominant_high_qcov"
    threshold = f"{ASSIGNMENT_BITSCORE_FRACTION:.0%}"
    assert row["support_note"] == (
        "The best alignment spans 180 of 200 query bases (90.0%). Of 2 retained "
        f"references, 1 has a bitscore at least {threshold} of the best bitscore; the "
        "available taxid from those references resolves to Alpha virus (species)."
    )
    assert "adjusted_taxid" not in row
    assert not any(column.startswith("second_") for column in row)
    assert "confidence_score" not in row


def test_query_big_table_sorts_rows_by_human_salience(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    crumbs_tsv = tmp_path / "crumbs.tsv"
    output = tmp_path / "query_big_table.tsv"
    weak_high_mass = base_hit(
        qseqid="nvdContigQuery_sample-1_weak",
        query_class="short_assembly_contig",
        qlen="200",
        length="50",
        bitscore="100",
        evalue="1e-20",
    )
    strong_low_mass = base_hit(
        qseqid="nvdContigQuery_sample-1_strong_low",
        query_class="long_assembly_contig",
        qlen="200",
        length="190",
        bitscore="120",
        evalue="1e-60",
    )
    strong_high_mass = base_hit(
        qseqid="nvdContigQuery_sample-1_strong_high",
        query_class="long_assembly_contig",
        qlen="200",
        length="190",
        bitscore="120",
        evalue="1e-60",
    )
    write_tsv(
        blast_tsv, [weak_high_mass, strong_low_mass, strong_high_mass], BLAST_COLUMNS
    )
    write_tsv(
        crumbs_tsv,
        [
            {
                "sample_id": "sample-1",
                "qseqid": weak_high_mass["qseqid"],
                "crumbs_score": "10000",
            },
            {
                "sample_id": "sample-1",
                "qseqid": strong_low_mass["qseqid"],
                "crumbs_score": "10",
            },
            {
                "sample_id": "sample-1",
                "qseqid": strong_high_mass["qseqid"],
                "crumbs_score": "100",
            },
        ],
        ["sample_id", "qseqid", "crumbs_score"],
    )

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
            "--crumbs-tsv",
            str(crumbs_tsv),
            "--output",
            str(output),
        ]
    )

    rows = read_tsv(output)
    assert [row["qseqid"] for row in rows] == [
        "nvdContigQuery_sample-1_strong_high",
        "nvdContigQuery_sample-1_strong_low",
        "nvdContigQuery_sample-1_weak",
    ]


def test_query_big_table_describes_low_coverage_lca_assignment(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    output = tmp_path / "query_big_table.tsv"
    rows = [
        base_hit(
            qseqid="nvdContigQuery_sample-1_000002",
            query_class="short_assembly_contig",
            qlen="200",
            length="60",
            bitscore="100",
            evalue="1e-15",
            rank="acellular root:Viruses; family:Alpha; genus:Alphavirus; species:Alpha virus",
            adjusted_taxid="10239",
            adjusted_taxid_name="Viruses",
            adjusted_taxid_rank="superkingdom",
            adjustment_method="lca",
        ),
        base_hit(
            qseqid="nvdContigQuery_sample-1_000002",
            query_class="short_assembly_contig",
            sseqid="ref-host",
            stitle="Host reference",
            length="55",
            bitscore="95",
            evalue="1e-14",
            staxids="9606",
            rank="cellular root:cellular organisms; domain:Eukaryota; family:Hominidae; genus:Homo; species:Homo sapiens",
            adjusted_taxid="10239",
            adjusted_taxid_name="Viruses",
            adjusted_taxid_rank="superkingdom",
            adjustment_method="lca",
        ),
    ]
    write_tsv(blast_tsv, rows, BLAST_COLUMNS)

    main(["--blast-tsv", str(blast_tsv), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["best_hit_qcov"] == "0.3"
    assert row["retained_reference_count"] == "2"
    assert row["assignment_reference_count"] == "2"
    assert row["assignment_taxid_count"] == "2"
    assert row["support_tier"] == "weak"
    assert row["support_tier_rule"] == "short_contig_lca_low_qcov"
    assert row["support_note"] == (
        "The best alignment spans 60 of 200 query bases (30.0%). Of 2 retained "
        "references, 2 have bitscores at least 95% of the best bitscore and map "
        "to 2 taxids; their LCA is Viruses (superkingdom)."
    )


def test_query_big_table_handles_single_retained_reference(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    output = tmp_path / "query_big_table.tsv"
    write_tsv(blast_tsv, [base_hit()], BLAST_COLUMNS)

    main(["--blast-tsv", str(blast_tsv), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["retained_reference_count"] == "1"
    assert row["assignment_reference_count"] == "1"
    assert row["assignment_taxid_count"] == "1"
    assert row["support_tier"] == "strong"


def test_query_big_table_counts_references_before_expanded_canonical_taxids(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    output = tmp_path / "query_big_table.tsv"
    query = "nvdContigQuery_sample-1_multi_taxid"
    rows = [
        base_hit(
            qseqid=query,
            sseqid="ref-multi",
            staxids="111",
            bitscore="100",
            adjusted_taxid="10239",
            adjusted_taxid_name="Viruses",
            adjusted_taxid_rank="superkingdom",
            adjustment_method="lca",
        ),
        base_hit(
            qseqid=query,
            sseqid="ref-multi",
            staxids="112",
            bitscore="100",
            adjusted_taxid="10239",
            adjusted_taxid_name="Viruses",
            adjusted_taxid_rank="superkingdom",
            adjustment_method="lca",
        ),
        base_hit(
            qseqid=query,
            sseqid="ref-boundary",
            staxids="113",
            bitscore="95",
            adjusted_taxid="10239",
            adjusted_taxid_name="Viruses",
            adjusted_taxid_rank="superkingdom",
            adjustment_method="lca",
        ),
        base_hit(
            qseqid=query,
            sseqid="ref-below",
            staxids="114",
            bitscore="94.9",
            adjusted_taxid="10239",
            adjusted_taxid_name="Viruses",
            adjusted_taxid_rank="superkingdom",
            adjustment_method="lca",
        ),
    ]
    write_tsv(blast_tsv, rows, BLAST_COLUMNS)

    main(["--blast-tsv", str(blast_tsv), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["retained_reference_count"] == "3"
    assert row["assignment_reference_count"] == "2"
    assert row["assignment_taxid_count"] == "3"


def test_query_big_table_excludes_null_normalized_taxids_from_assignment_count(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    output = tmp_path / "query_big_table.tsv"
    query = "nvdContigQuery_sample-1_missing_taxid"
    rows = [
        base_hit(qseqid=query, sseqid="ref-valid", staxids="111"),
        base_hit(qseqid=query, sseqid="ref-unavailable", staxids=""),
    ]
    write_tsv(blast_tsv, rows, BLAST_COLUMNS)

    main(["--blast-tsv", str(blast_tsv), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["retained_reference_count"] == "2"
    assert row["assignment_reference_count"] == "2"
    assert row["assignment_taxid_count"] == "1"
    assert row["support_note"] == (
        "The best alignment spans 180 of 200 query bases (90.0%). Of 2 retained "
        "references, 2 have bitscores at least 95% of the best bitscore; the "
        "available taxid from those references resolves to Alpha virus (species)."
    )


def test_query_big_table_uses_moderate_tier_for_partial_coverage(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    output = tmp_path / "query_big_table.tsv"
    write_tsv(blast_tsv, [base_hit(length="120")], BLAST_COLUMNS)

    main(["--blast-tsv", str(blast_tsv), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["best_hit_qcov"] == "0.6"
    assert row["support_tier"] == "moderate"
    assert row["support_tier_rule"] == "long_contig_dominant_partial_qcov"


def test_query_big_table_selects_best_hit_with_all_tie_breakers(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    output = tmp_path / "query_big_table.tsv"
    rows = [
        base_hit(sseqid="ref-worse-evalue", evalue="1e-40", pident="100", length="200"),
        base_hit(sseqid="ref-worse-pident", evalue="1e-50", pident="98", length="200"),
        base_hit(sseqid="ref-shorter", evalue="1e-50", pident="99", length="170"),
        base_hit(sseqid="ref-best", evalue="1e-50", pident="99", length="180"),
    ]
    write_tsv(blast_tsv, rows, BLAST_COLUMNS)

    main(["--blast-tsv", str(blast_tsv), "--output", str(output)])

    [row] = read_tsv(output)
    assert row["best_hit_evalue"] == "1e-50"
    assert row["best_hit_pident"] == "99"
    assert row["best_hit_alignment_length"] == "180"


def test_query_big_table_collapses_conflicting_crumbs_scores_to_absence(
    tmp_path: Path,
) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    crumbs_tsv = tmp_path / "crumbs.tsv"
    output = tmp_path / "query_big_table.tsv"
    row = base_hit()
    write_tsv(blast_tsv, [row], BLAST_COLUMNS)
    write_tsv(
        crumbs_tsv,
        [
            {"sample_id": "sample-1", "qseqid": row["qseqid"], "crumbs_score": "10"},
            {"sample_id": "sample-1", "qseqid": row["qseqid"], "crumbs_score": "20"},
        ],
        ["sample_id", "qseqid", "crumbs_score"],
    )

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
            "--crumbs-tsv",
            str(crumbs_tsv),
            "--output",
            str(output),
        ]
    )

    [result] = read_tsv(output)
    assert result["crumbs_score"] == ""


def test_query_big_table_reports_missing_consumed_columns(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "final_blast.tsv"
    output = tmp_path / "query_big_table.tsv"
    columns_without_producer = [column for column in BLAST_COLUMNS if column != "producer"]
    row = base_hit()
    write_tsv(
        blast_tsv,
        [{key: value for key, value in row.items() if key != "producer"}],
        columns_without_producer,
    )

    with pytest.raises(ValueError, match="producer"):
        main(["--blast-tsv", str(blast_tsv), "--output", str(output)])
