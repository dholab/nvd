"""Tests for merged taxonomic profile sunburst rendering."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

from render_multisample_taxburst import render_html

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
