"""Tests for merged taxonomic profile sunburst rendering."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING

from render_multisample_taxburst import main, render_html

if TYPE_CHECKING:
    from pathlib import Path


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "lineage", "fraction", "f_weighted_at_rank"],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_crumbs_taxa(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "taxon_id",
                "taxon_name",
                "rank",
                "taxpath",
                "taxpathsn",
                "rankpath",
                "percentage_emitted",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def test_per_sample_crumbs_html_preserves_full_ncbi_paths(tmp_path: Path) -> None:
    taxa = tmp_path / "madison.crumbs.taxa.tsv"
    write_crumbs_taxa(
        taxa,
        [
            {
                "sample_id": "Madison",
                "taxon_id": "11308",
                "taxon_name": "Orthomyxoviridae",
                "rank": "family",
                "taxpath": "10239|2559587|2732396|11308",
                "taxpathsn": "Viruses|Riboviria|Negarnaviricota|Orthomyxoviridae",
                "rankpath": "acellular root|realm|phylum|family",
                "percentage_emitted": "60",
            },
            {
                "sample_id": "Madison",
                "taxon_id": "10912",
                "taxon_name": "Rotavirus",
                "rank": "genus",
                "taxpath": "10239|2559587|2732396|10912",
                "taxpathsn": "Viruses|Riboviria|Negarnaviricota|Rotavirus",
                "rankpath": "acellular root|realm|phylum|genus",
                "percentage_emitted": "30",
            },
            {
                "sample_id": "Madison",
                "taxon_id": "1",
                "taxon_name": "root",
                "rank": "no rank",
                "taxpath": "1",
                "taxpathsn": "root",
                "rankpath": "no rank",
                "percentage_emitted": "10",
            },
        ],
    )

    html = render_html([("Madison", taxa)], input_format="crumbs")

    assert "<dataset>Madison</dataset>" in html
    assert '<node name="Viruses">' in html
    assert "<rank><val>acellular root</val></rank>" in html
    assert '<node name="Riboviria">' in html
    assert '<node name="Orthomyxoviridae">' in html
    assert '<node name="Rotavirus">' in html
    assert '<node name="root">' in html
    assert "<count><val>1000.0</val></count>" in html


def test_per_sample_crumbs_cli_writes_taxburst_json(tmp_path: Path) -> None:
    taxa = tmp_path / "sample.crumbs.taxa.tsv"
    output_html = tmp_path / "sample.crumbs.taxburst.html"
    output_json = tmp_path / "sample.crumbs.taxburst.json"
    write_crumbs_taxa(
        taxa,
        [
            {
                "sample_id": "sample",
                "taxon_id": "11308",
                "taxon_name": "Orthomyxoviridae",
                "rank": "family",
                "taxpath": "10239|11308",
                "taxpathsn": "Viruses|Orthomyxoviridae",
                "rankpath": "acellular root|family",
                "percentage_emitted": "100",
            },
        ],
    )

    main(
        [
            "--input-format",
            "crumbs",
            "--summary",
            f"sample={taxa}",
            "--output",
            str(output_html),
            "--output-json",
            str(output_json),
        ],
    )

    assert output_html.is_file()
    assert json.loads(output_json.read_text(encoding="utf-8")) == [
        {
            "name": "Viruses",
            "rank": "acellular root",
            "count": 1000.0,
            "children": [
                {
                    "name": "Orthomyxoviridae",
                    "rank": "family",
                    "count": 1000.0,
                },
            ],
        },
    ]


def test_crumbs_rendering_accumulates_direct_and_duplicate_mass(tmp_path: Path) -> None:
    taxa = tmp_path / "sample.crumbs.taxa.tsv"
    write_crumbs_taxa(
        taxa,
        [
            {
                "sample_id": "sample",
                "taxon_id": "10239",
                "taxon_name": "Viruses",
                "rank": "acellular root",
                "taxpath": "10239",
                "taxpathsn": "Viruses",
                "rankpath": "acellular root",
                "percentage_emitted": "10",
            },
            {
                "sample_id": "sample",
                "taxon_id": "10780",
                "taxon_name": "Parvoviridae",
                "rank": "family",
                "taxpath": "10239|10780",
                "taxpathsn": "Viruses|Parvoviridae",
                "rankpath": "acellular root|family",
                "percentage_emitted": "40",
            },
            {
                "sample_id": "sample",
                "taxon_id": "10780",
                "taxon_name": "Parvoviridae",
                "rank": "family",
                "taxpath": "10239|10780",
                "taxpathsn": "Viruses|Parvoviridae",
                "rankpath": "acellular root|family",
                "percentage_emitted": "50",
            },
        ],
    )

    html = render_html([("sample", taxa)], input_format="crumbs")

    assert '<node name="Viruses">' in html
    assert '<node name="Parvoviridae">' in html
    expected_full_mass_nodes = 2
    assert html.count("<count><val>1000.0</val></count>") == expected_full_mass_nodes
    assert "<count><val>900.0</val></count>" in html


def test_merged_crumbs_html_conserves_each_samples_mass(tmp_path: Path) -> None:
    water = tmp_path / "water.crumbs.taxa.tsv"
    madison = tmp_path / "madison.crumbs.taxa.tsv"
    write_crumbs_taxa(
        water,
        [
            {
                "sample_id": "Water",
                "taxon_id": "10780",
                "taxon_name": "Parvoviridae",
                "rank": "family",
                "taxpath": "10239|10780",
                "taxpathsn": "Viruses|Parvoviridae",
                "rankpath": "acellular root|family",
                "percentage_emitted": "40",
            },
            {
                "sample_id": "Water",
                "taxon_id": "10780",
                "taxon_name": "Parvoviridae",
                "rank": "family",
                "taxpath": "10239|10780",
                "taxpathsn": "Viruses|Parvoviridae",
                "rankpath": "acellular root|family",
                "percentage_emitted": "60",
            },
        ],
    )
    write_crumbs_taxa(
        madison,
        [
            {
                "sample_id": "Madison",
                "taxon_id": "11308",
                "taxon_name": "Orthomyxoviridae",
                "rank": "family",
                "taxpath": "10239|11308",
                "taxpathsn": "Viruses|Orthomyxoviridae",
                "rankpath": "acellular root|family",
                "percentage_emitted": "25",
            },
            {
                "sample_id": "Madison",
                "taxon_id": "10912",
                "taxon_name": "Rotavirus",
                "rank": "genus",
                "taxpath": "10239|10912",
                "taxpathsn": "Viruses|Rotavirus",
                "rankpath": "acellular root|genus",
                "percentage_emitted": "75",
            },
        ],
    )

    html = render_html(
        [("Water", water), ("Madison", madison)],
        input_format="crumbs",
    )

    assert "<count><val>1000.0</val><val>1000.0</val></count>" in html
    assert '<node name="Parvoviridae">' in html
    assert "<count><val>1000.0</val><val>0.0</val></count>" in html
    assert '<node name="Orthomyxoviridae">' in html
    assert "<count><val>0.0</val><val>250.0</val></count>" in html
    assert '<node name="Rotavirus">' in html


def test_crumbs_rendering_keeps_duplicate_names_with_distinct_taxon_ids(
    tmp_path: Path,
) -> None:
    taxa = tmp_path / "sample.crumbs.taxa.tsv"
    write_crumbs_taxa(
        taxa,
        [
            {
                "sample_id": "sample",
                "taxon_id": "1001",
                "taxon_name": "Shared name",
                "rank": "species",
                "taxpath": "10239|1001",
                "taxpathsn": "Viruses|Shared name",
                "rankpath": "acellular root|species",
                "percentage_emitted": "40",
            },
            {
                "sample_id": "sample",
                "taxon_id": "1002",
                "taxon_name": "Shared name",
                "rank": "species",
                "taxpath": "10239|1002",
                "taxpathsn": "Viruses|Shared name",
                "rankpath": "acellular root|species",
                "percentage_emitted": "60",
            },
        ],
    )

    html = render_html([("sample", taxa)], input_format="crumbs")

    expected_nodes = 2
    assert html.count('<node name="Shared name">') == expected_nodes
    assert "<count><val>400.0</val></count>" in html
    assert "<count><val>600.0</val></count>" in html


def test_merged_html_contains_one_dataset_per_sample(tmp_path: Path) -> None:
    sample_a = tmp_path / "sample_a.csv"
    sample_b = tmp_path / "sample_b.csv"
    write_summary(
        sample_a,
        [
            {
                "rank": "superkingdom",
                "lineage": "Viruses",
                "fraction": "0.2",
                "f_weighted_at_rank": "0.2",
            },
        ],
    )
    write_summary(
        sample_b,
        [
            {
                "rank": "superkingdom",
                "lineage": "Viruses",
                "fraction": "0.7",
                "f_weighted_at_rank": "0.7",
            },
        ],
    )

    html = render_html([("S25", sample_a), ("S26", sample_b)])

    assert "<dataset>S25</dataset>" in html
    assert "<dataset>S26</dataset>" in html
    assert "#lastDataset { font: 11px sans-serif !important; }" in html
    assert '<node name="Viruses">' in html
    assert "<count><val>200.0</val><val>700.0</val></count>" in html


def test_merged_html_uses_zero_for_absent_taxa(tmp_path: Path) -> None:
    sample_a = tmp_path / "sample_a.csv"
    sample_b = tmp_path / "sample_b.csv"
    write_summary(
        sample_a,
        [
            {
                "rank": "superkingdom",
                "lineage": "Viruses",
                "fraction": "1.0",
                "f_weighted_at_rank": "1.0",
            },
            {
                "rank": "phylum",
                "lineage": "Viruses;Pisuviricota",
                "fraction": "1.0",
                "f_weighted_at_rank": "1.0",
            },
        ],
    )
    write_summary(
        sample_b,
        [
            {
                "rank": "superkingdom",
                "lineage": "Viruses",
                "fraction": "1.0",
                "f_weighted_at_rank": "1.0",
            },
        ],
    )

    html = render_html([("S25", sample_a), ("S26", sample_b)])

    assert '<node name="Pisuviricota">' in html
    assert "<count><val>1000.0</val><val>0.0</val></count>" in html
