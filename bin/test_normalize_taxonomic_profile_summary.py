"""Tests for taxonomic profile summary display-name normalization."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

from normalize_taxonomic_profile_summary import (
    apply_lineage_replacements,
    normalize_taxonomic_profile_summary,
    read_normalization_map,
    write_normalization_map,
)

if TYPE_CHECKING:
    from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rank", "lineage", "fraction"])
        writer.writeheader()
        writer.writerows(rows)


def retained_lineages(rows: list[dict[str, str]]) -> list[str]:
    return [
        row["lineage"]
        for row in rows
        if not (row["lineage"] == "unclassified" and row["rank"] != "superkingdom")
    ]


def node_names(lineages: list[str]) -> list[str]:
    prefixes: set[tuple[str, ...]] = set()
    for lineage in lineages:
        parts = lineage.split(";")
        prefixes.update(tuple(parts[:depth]) for depth in range(1, len(parts) + 1))
    return [prefix[-1] for prefix in prefixes]


def assert_globally_unique_node_names(rows: list[dict[str, str]]) -> None:
    names = node_names(retained_lineages(rows))
    assert len(names) == len(set(names))


def assert_no_missing_parent_prefix(rows: list[dict[str, str]]) -> None:
    prefixes: set[tuple[str, ...]] = set()
    for lineage in retained_lineages(rows):
        parts = lineage.split(";")
        prefixes.update(tuple(parts[:depth]) for depth in range(1, len(parts) + 1))

    for prefix in prefixes:
        if len(prefix) > 1:
            assert prefix[:-1] in prefixes


def test_duplicate_unclassified_nodes_are_disambiguated_consistently(
    tmp_path: Path,
) -> None:
    """Repeated placeholder family names should not break tree structure."""
    source = tmp_path / "summary.csv"
    normalized = tmp_path / "normalized.csv"
    write_summary(
        source,
        [
            {
                "rank": "family",
                "lineage": "Viruses;unclassified phylum;unclassified class;Unclassified",
                "fraction": "0.1",
            },
            {
                "rank": "family",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;Unclassified",
                "fraction": "0.2",
            },
            {
                "rank": "genus",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;Unclassified;unclassified genus",
                "fraction": "0.3",
            },
            {
                "rank": "family",
                "lineage": "Viruses;Kitrinoviricota;Alsuviricetes;Hepelivirales;Unclassified",
                "fraction": "0.4",
            },
        ],
    )

    normalize_taxonomic_profile_summary(input_path=source, output_path=normalized)

    rows = read_rows(normalized)
    lineages = [row["lineage"] for row in rows]
    assert_globally_unique_node_names(rows)
    assert_no_missing_parent_prefix(rows)
    assert lineages[1] == (
        "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;"
        "Unclassified [family under Picornavirales]"
    )
    assert lineages[2] == (
        "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;"
        "Unclassified [family under Picornavirales];unclassified genus"
    )
    assert [row["fraction"] for row in rows] == ["0.1", "0.2", "0.3", "0.4"]


def test_duplicate_unclassified_genus_uses_original_parent_context(
    tmp_path: Path,
) -> None:
    """Child labels should not recursively embed normalized parent labels."""
    source = tmp_path / "summary.csv"
    normalized = tmp_path / "normalized.csv"
    write_summary(
        source,
        [
            {
                "rank": "genus",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;Unclassified;unclassified genus",
                "fraction": "0.1",
            },
            {
                "rank": "genus",
                "lineage": "Viruses;unclassified phylum;unclassified class;Unclassified;Unclassified;unclassified genus",
                "fraction": "0.2",
            },
            {
                "rank": "genus",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;Picornaviridae;unclassified genus",
                "fraction": "0.3",
            },
        ],
    )

    normalize_taxonomic_profile_summary(input_path=source, output_path=normalized)

    rows = read_rows(normalized)
    input_depths = [len(row["lineage"].split(";")) for row in read_rows(source)]
    output_depths = [len(row["lineage"].split(";")) for row in rows]
    joined_lineages = "\n".join(row["lineage"] for row in rows)
    assert_globally_unique_node_names(rows)
    assert_no_missing_parent_prefix(rows)
    assert output_depths == input_depths
    assert (
        "unclassified genus [genus under Picornavirales > Unclassified]"
        in joined_lineages
    )
    assert "under Unclassified [family" not in joined_lineages


def test_unique_summary_is_written_unchanged(tmp_path: Path) -> None:
    """Already unique taxonomic summaries should be a textual no-op."""
    source = tmp_path / "summary.csv"
    normalized = tmp_path / "normalized.csv"
    source.write_text(
        "rank,lineage,fraction\n"
        "superkingdom,Viruses,1.0\n"
        "phylum,Viruses;Pisuviricota,0.5\n"
        "class,Viruses;Pisuviricota;Pisoniviricetes,0.5\n",
        encoding="utf-8",
    )

    normalize_taxonomic_profile_summary(input_path=source, output_path=normalized)

    assert normalized.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")


def test_non_superkingdom_bare_unclassified_rows_are_left_alone(
    tmp_path: Path,
) -> None:
    """Rows ignored by taxburst's tree parser should not be renamed."""
    source = tmp_path / "summary.csv"
    normalized = tmp_path / "normalized.csv"
    write_summary(
        source,
        [
            {"rank": "superkingdom", "lineage": "unclassified", "fraction": "0.7"},
            {"rank": "phylum", "lineage": "unclassified", "fraction": "0.2"},
            {"rank": "class", "lineage": "unclassified", "fraction": "0.1"},
        ],
    )

    normalize_taxonomic_profile_summary(input_path=source, output_path=normalized)

    assert read_rows(normalized) == read_rows(source)


def test_global_map_makes_sample_normalization_consistent(tmp_path: Path) -> None:
    """Run-level maps should avoid sample-specific display names."""
    sample_a = tmp_path / "sample_a.csv"
    sample_b = tmp_path / "sample_b.csv"
    normalized_a = tmp_path / "sample_a.normalized.csv"
    normalized_b = tmp_path / "sample_b.normalized.csv"
    normalization_map = tmp_path / "normalization_map.json"
    shared_lineage = "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;Unclassified"
    write_summary(
        sample_a,
        [
            {"rank": "superkingdom", "lineage": "Viruses", "fraction": "1.0"},
            {
                "rank": "phylum",
                "lineage": "Viruses;Pisuviricota",
                "fraction": "1.0",
            },
            {
                "rank": "class",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes",
                "fraction": "1.0",
            },
            {
                "rank": "order",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales",
                "fraction": "1.0",
            },
            {"rank": "family", "lineage": shared_lineage, "fraction": "1.0"},
        ],
    )
    write_summary(
        sample_b,
        [
            {"rank": "superkingdom", "lineage": "Viruses", "fraction": "1.0"},
            {
                "rank": "phylum",
                "lineage": "Viruses;Pisuviricota",
                "fraction": "0.5",
            },
            {
                "rank": "class",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes",
                "fraction": "0.5",
            },
            {
                "rank": "order",
                "lineage": "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales",
                "fraction": "0.5",
            },
            {"rank": "family", "lineage": shared_lineage, "fraction": "0.5"},
            {
                "rank": "phylum",
                "lineage": "Viruses;Kitrinoviricota",
                "fraction": "0.5",
            },
            {
                "rank": "class",
                "lineage": "Viruses;Kitrinoviricota;Alsuviricetes",
                "fraction": "0.5",
            },
            {
                "rank": "order",
                "lineage": "Viruses;Kitrinoviricota;Alsuviricetes;Hepelivirales",
                "fraction": "0.5",
            },
            {
                "rank": "family",
                "lineage": "Viruses;Kitrinoviricota;Alsuviricetes;Hepelivirales;Unclassified",
                "fraction": "0.5",
            },
        ],
    )

    write_normalization_map(
        input_paths=[sample_a, sample_b],
        output_path=normalization_map,
    )
    replacements = read_normalization_map(normalization_map)
    apply_lineage_replacements(
        input_path=sample_a,
        output_path=normalized_a,
        replacements=replacements,
    )
    apply_lineage_replacements(
        input_path=sample_b,
        output_path=normalized_b,
        replacements=replacements,
    )

    expected = (
        "Viruses;Pisuviricota;Pisoniviricetes;Picornavirales;"
        "Unclassified [family under Picornavirales]"
    )
    assert read_rows(normalized_a)[-1]["lineage"] == expected
    assert read_rows(normalized_b)[4]["lineage"] == expected
