"""Tests for higher-risk taxonomic findings report inputs."""

from __future__ import annotations

from typing import TYPE_CHECKING

from py_nvd.multiqc_domains import higher_risk_taxonomic_findings_section
from py_nvd.multiqc_packages import ReportPackage, TaxonBigTableReceipt
from py_nvd.multiqc_taxonomy import (
    collect_taxon_big_tables,
    higher_risk_taxon_rows,
)

if TYPE_CHECKING:
    from pathlib import Path


TAXON_HEADER = (
    "sample_id\ttaxon_name\ttaxon_rank\tsupport_tier\ttaxon_crumbs\t"
    "relative_crumbs_percent\tsupporting_query_count\ttaxid\twho_risk_group\t"
    "total_query_span\ttotal_crumbs_score\tstrong_query_count\t"
    "moderate_query_count\tweak_query_count\treview_query_count\t"
    "redacted_query_count\tsupporting_genome_like_contig_count\t"
    "supporting_long_contig_count\tsupporting_short_contig_count\t"
    "supporting_merged_pair_count\tsupporting_single_read_count\t"
    "support_tier_rule\tsupport_note"
)


def write_taxon_package(
    root: Path,
    *,
    sample_id: str,
    rows: tuple[tuple[str, ...], ...],
) -> ReportPackage:
    root.mkdir()
    (root / "table.payload").write_text(
        TAXON_HEADER + "\n" + "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return ReportPackage(
        root=root,
        receipt=TaxonBigTableReceipt(sample_id=sample_id),
    )


def test_higher_risk_taxon_rows_preserve_big_table_information_and_sort_by_risk(
    tmp_path: Path,
) -> None:
    sample_a = write_taxon_package(
        tmp_path / "sample_a",
        sample_id="sample_A",
        rows=(
            (
                "sample_A",
                "Taxon alpha",
                "species",
                "strong",
                "8.5",
                "12.5",
                "4",
                "101",
                "RG2",
                "4200",
                "35700",
                "3",
                "1",
                "0",
                "0",
                "0",
                "1",
                "1",
                "2",
                "0",
                "0",
                "multi_query_strong_support",
                "3 strong query assignments support Taxon alpha (species).",
            ),
        ),
    )
    sample_b = write_taxon_package(
        tmp_path / "sample_b",
        sample_id="sample_B",
        rows=(
            (
                "sample_B",
                "Taxon beta",
                "genus",
                "weak",
                "1.25",
                "2.0",
                "2",
                "202",
                "RG3",
                "250",
                "312.5",
                "0",
                "0",
                "2",
                "0",
                "0",
                "0",
                "0",
                "0",
                "0",
                "2",
                "all_weak_queries",
                "All 2 query assignments have low best-hit query coverage.",
            ),
            (
                "sample_B",
                "Taxon background",
                "family",
                "strong",
                "50.0",
                "80.0",
                "100",
                "303",
                "RG1",
                "100000",
                "5000000",
                "100",
                "0",
                "0",
                "0",
                "0",
                "0",
                "2",
                "98",
                "0",
                "0",
                "multi_query_strong_support",
                "100 strong query assignments support Taxon background (family).",
            ),
        ),
    )

    tables = collect_taxon_big_tables((sample_a, sample_b))
    rows, eligible_count = higher_risk_taxon_rows(tables, max_rows=100)

    assert eligible_count == 2
    assert [row.who_risk_group for row in rows] == ["RG3", "RG2"]
    assert rows[0].model_dump(exclude={"taxon_key"}) == {
        "sample_id": "sample_B",
        "taxon_name": "Taxon beta",
        "taxon_rank": "genus",
        "support_tier": "weak",
        "taxon_crumbs": 1.25,
        "relative_crumbs_percent": 2.0,
        "supporting_query_count": 2,
        "taxid": 202,
        "who_risk_group": "RG3",
        "total_query_span": 250,
        "total_crumbs_score": 312.5,
        "strong_query_count": 0,
        "moderate_query_count": 0,
        "weak_query_count": 2,
        "review_query_count": 0,
        "redacted_query_count": 0,
        "supporting_genome_like_contig_count": 0,
        "supporting_long_contig_count": 0,
        "supporting_short_contig_count": 0,
        "supporting_merged_pair_count": 0,
        "supporting_single_read_count": 2,
        "support_tier_rule": "all_weak_queries",
        "support_note": "All 2 query assignments have low best-hit query coverage.",
        "problem": None,
        "validation_details": None,
    }


def test_higher_risk_section_exposes_the_taxon_big_table_columns(
    tmp_path: Path,
) -> None:
    package = write_taxon_package(
        tmp_path / "sample_a",
        sample_id="sample_A",
        rows=(
            (
                "sample_A",
                "Taxon alpha",
                "species",
                "strong",
                "8.5",
                "12.5",
                "4",
                "101",
                "RG4",
                "4200",
                "35700",
                "3",
                "1",
                "0",
                "0",
                "0",
                "1",
                "1",
                "2",
                "0",
                "0",
                "multi_query_strong_support",
                "3 strong query assignments support Taxon alpha (species).",
            ),
        ),
    )

    section = higher_risk_taxonomic_findings_section((package,))

    assert section is not None
    assert section.section_name == "Higher-Risk Taxonomic Findings"
    assert section.description == (
        "Showing 1 of 1 WHO risk group 2-4 findings. "
        "Complete results: taxon_big_table.tsv."
    )
    assert tuple(section.headers) == (
        "taxon_name",
        "taxid",
        "taxon_rank",
        "who_risk_group",
        "support_tier",
        "support_note",
        "taxon_crumbs",
        "relative_crumbs_percent",
        "total_crumbs_score",
        "supporting_query_count",
        "total_query_span",
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
        "problem",
        "validation_details",
    )
