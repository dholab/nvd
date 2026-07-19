"""Tests for higher-risk taxonomic findings report inputs."""

from __future__ import annotations

from typing import TYPE_CHECKING

from py_nvd.multiqc_domains import (
    higher_risk_taxonomic_findings_section,
    taxonomy_sections,
)
from py_nvd.multiqc_packages import ReportPackage, TaxonBigTableReceipt
from py_nvd.multiqc_taxonomy import (
    ParsedTaxonBigTable,
    collect_taxon_big_tables,
    faintest_taxon_rows,
    higher_risk_taxon_rows,
    strongest_taxon_rows,
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


def taxon_row(
    sample_id: str,
    taxon_name: str,
    *,
    taxid: int,
    taxon_crumbs: str,
    risk_group: str = "",
) -> tuple[str, ...]:
    return (
        sample_id,
        taxon_name,
        "species",
        "strong",
        taxon_crumbs,
        taxon_crumbs,
        "2",
        str(taxid),
        risk_group,
        "200",
        "400",
        "2",
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        "2",
        "0",
        "0",
        "multi_query_strong_support",
        f"2 strong query assignments support {taxon_name} (species).",
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
        "Rows may also appear in the signal-strength views. "
        "All taxon signals are detailed in the Taxon Big Table. "
        "Queries supporting each taxon are detailed in the Query Big Table."
    )
    assert tuple(section.headers) == (
        "taxon_name",
        "taxon_rank",
        "support_tier",
        "taxon_crumbs",
        "relative_crumbs_percent",
        "supporting_query_count",
        "taxid",
        "who_risk_group",
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
        "problem",
        "validation_details",
    )


def test_blank_optional_risk_group_does_not_invalidate_higher_risk_rows(
    tmp_path: Path,
) -> None:
    package = write_taxon_package(
        tmp_path / "sample_a",
        sample_id="sample_A",
        rows=(
            taxon_row(
                "sample_A",
                "Higher-risk taxon",
                taxid=101,
                taxon_crumbs="12.0",
                risk_group="RG2",
            ),
            taxon_row(
                "sample_A",
                "Unannotated family",
                taxid=202,
                taxon_crumbs="3.0",
            ),
        ),
    )

    tables = collect_taxon_big_tables((package,))
    rows, eligible_count = higher_risk_taxon_rows(tables)

    assert isinstance(tables[0], ParsedTaxonBigTable)
    assert eligible_count == 1
    assert [row.taxon_name for row in rows] == ["Higher-risk taxon"]


def test_strongest_and_faintest_taxon_rows_are_bounded_per_sample(
    tmp_path: Path,
) -> None:
    sample_a = write_taxon_package(
        tmp_path / "sample_a",
        sample_id="sample_A",
        rows=(
            taxon_row(
                "sample_A",
                "A strongest",
                taxid=101,
                taxon_crumbs="20.0",
            ),
            taxon_row(
                "sample_A",
                "A faintest",
                taxid=102,
                taxon_crumbs="1.0",
            ),
            taxon_row(
                "sample_A",
                "A zero",
                taxid=103,
                taxon_crumbs="0.0",
            ),
            taxon_row(
                "sample_A",
                "A unknown",
                taxid=104,
                taxon_crumbs="",
            ),
        ),
    )
    sample_b = write_taxon_package(
        tmp_path / "sample_b",
        sample_id="sample_B",
        rows=(
            taxon_row(
                "sample_B",
                "B strongest",
                taxid=201,
                taxon_crumbs="8.0",
            ),
            taxon_row(
                "sample_B",
                "B faintest",
                taxid=202,
                taxon_crumbs="0.5",
            ),
        ),
    )
    tables = collect_taxon_big_tables((sample_a, sample_b))

    strongest, strongest_eligible = strongest_taxon_rows(
        tables,
        max_rows_per_sample=1,
    )
    faintest, faintest_eligible = faintest_taxon_rows(
        tables,
        max_rows_per_sample=1,
    )

    assert strongest_eligible == 4
    assert [row.taxon_name for row in strongest] == ["A strongest", "B strongest"]
    assert faintest_eligible == 4
    assert [row.taxon_name for row in faintest] == ["A faintest", "B faintest"]


def test_taxonomy_sections_lead_with_risk_strongest_and_faintest_views(
    tmp_path: Path,
) -> None:
    package = write_taxon_package(
        tmp_path / "sample_a",
        sample_id="sample_A",
        rows=(
            taxon_row(
                "sample_A",
                "Taxon alpha",
                taxid=101,
                taxon_crumbs="8.5",
                risk_group="RG3",
            ),
        ),
    )

    sections = taxonomy_sections((package,))

    assert tuple(sections) == (
        "nvd_higher_risk_taxonomic_findings",
        "nvd_strongest_taxonomic_signals",
        "nvd_faintest_taxonomic_signals",
    )
    assert sections["nvd_strongest_taxonomic_signals"].section_name == (
        "Strongest Taxonomic Signals"
    )
    assert sections["nvd_faintest_taxonomic_signals"].section_name == (
        "Faintest Taxonomic Signals"
    )
    for section in sections.values():
        assert "taxon_big_table.tsv" not in section.description
        assert "All taxon signals are detailed in the Taxon Big Table." in (
            section.description
        )
        assert (
            "Queries supporting each taxon are detailed in the Query Big Table."
            in section.description
        )
    for section_id in (
        "nvd_strongest_taxonomic_signals",
        "nvd_faintest_taxonomic_signals",
    ):
        assert "Taxon CRUMBS greater than zero" in sections[section_id].description
        assert "positive taxon signals" not in sections[section_id].description


def test_default_signal_budgets_keep_ten_strongest_and_five_faintest(
    tmp_path: Path,
) -> None:
    package = write_taxon_package(
        tmp_path / "sample_a",
        sample_id="sample_A",
        rows=tuple(
            taxon_row(
                "sample_A",
                f"Signal {signal}",
                taxid=100 + signal,
                taxon_crumbs=str(signal),
            )
            for signal in range(1, 13)
        ),
    )
    tables = collect_taxon_big_tables((package,))

    strongest, _strongest_eligible = strongest_taxon_rows(tables)
    faintest, _faintest_eligible = faintest_taxon_rows(tables)

    assert [row.taxon_name for row in strongest] == [
        "Signal 12",
        "Signal 11",
        "Signal 10",
        "Signal 9",
        "Signal 8",
        "Signal 7",
        "Signal 6",
        "Signal 5",
        "Signal 4",
        "Signal 3",
    ]
    assert [row.taxon_name for row in faintest] == [
        "Signal 1",
        "Signal 2",
        "Signal 3",
        "Signal 4",
        "Signal 5",
    ]
