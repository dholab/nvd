"""Behavior tests for WHO risk-group annotation."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

from annotate_risk_groups import main
from py_nvd import taxonomy

if TYPE_CHECKING:
    from pathlib import Path


def write_table(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
    *,
    delimiter: str,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def read_table(path: Path, *, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_lookup(path: Path, rows: list[dict[str, str]]) -> None:
    normalized = [
        {
            "canonical_taxid": "",
            "lookup_name": "",
            "human_classification": "",
            **row,
        }
        for row in rows
    ]
    write_table(
        path,
        normalized,
        ["canonical_taxid", "lookup_name", "human_classification"],
        delimiter="\t",
    )


def run_annotation(
    command: str,
    input_path: Path,
    lookup_path: Path,
    output_path: Path,
    extra_args: list[str],
) -> None:
    args = [
        command,
        "--input",
        str(input_path),
        "--lookup",
        str(lookup_path),
        "--output",
        str(output_path),
        *extra_args,
    ]
    main(args)


def run_sourmash_annotation(
    input_path: Path,
    bioboxes_path: Path,
    lookup_path: Path,
    output_path: Path,
) -> None:
    run_annotation(
        "sourmash-summary",
        input_path,
        lookup_path,
        output_path,
        ["--bioboxes", str(bioboxes_path)],
    )


def run_blast_annotation(
    input_path: Path,
    lookup_path: Path,
    taxonomy_dir: Path,
    output_path: Path,
) -> None:
    run_annotation(
        "blast",
        input_path,
        lookup_path,
        output_path,
        ["--taxonomy-dir", str(taxonomy_dir)],
    )


def write_bioboxes(path: Path) -> None:
    path.write_text(
        "# Taxonomic Profiling Output\n"
        "@SampleID:sample-1\n"
        "@Version:0.10.0\n"
        "@Ranks:superkingdom|species|strain\n"
        "@__program__:sourmash\n"
        "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n"
        "123\tspecies\t1|123\tViruses|Example virus\t25.0\n"
        "124\tstrain\t1|123|124\tViruses|Example virus|Example strain\t15.0\n"
        "125\tstrain\t1|123|125\tViruses|Example virus|Exceptional strain\t10.0\n",
        encoding="utf-8",
    )


def write_test_taxonomy(path: Path) -> Path:
    path.mkdir()
    (path / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "122\t|\t1\t|\tfamily\t|\n"
        "123\t|\t122\t|\tspecies\t|\n"
        "124\t|\t123\t|\tno rank\t|\n"
        "125\t|\t123\t|\tno rank\t|\n"
        "999\t|\t1\t|\tspecies\t|\n",
        encoding="utf-8",
    )
    (path / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "122\t|\tExample family\t|\t\t|\tscientific name\t|\n"
        "123\t|\tExample virus\t|\t\t|\tscientific name\t|\n"
        "124\t|\tExample strain\t|\t\t|\tscientific name\t|\n"
        "125\t|\tExceptional strain\t|\t\t|\tscientific name\t|\n"
        "999\t|\tOther virus\t|\t\t|\tscientific name\t|\n",
        encoding="utf-8",
    )
    (path / "merged.dmp").write_text("", encoding="utf-8")
    taxonomy.ensure_taxonomy_available(path)
    return path


def test_sourmash_summary_adds_risk_group_without_dropping_unmatched_rows(
    tmp_path: Path,
) -> None:
    summary = tmp_path / "sample.summarized.csv"
    lookup = tmp_path / "risk-groups.tsv"
    output = tmp_path / "sample.with-risk-groups.csv"
    bioboxes = tmp_path / "sample.bioboxes.profile"
    write_bioboxes(bioboxes)
    write_table(
        summary,
        [
            {
                "query_name": "sample-1",
                "rank": "species",
                "fraction": "0.25",
                "lineage": "Viruses;Example virus",
                "query_md5": "abc123",
            },
            {
                "query_name": "sample-1",
                "rank": "strain",
                "fraction": "0.15",
                "lineage": "Viruses;Example virus;Example strain",
                "query_md5": "abc123",
            },
            {
                "query_name": "sample-1",
                "rank": "strain",
                "fraction": "0.10",
                "lineage": "Viruses;Example virus;Exceptional strain",
                "query_md5": "abc123",
            },
            {
                "query_name": "sample-1",
                "rank": "species",
                "fraction": "0.75",
                "lineage": "unclassified",
                "query_md5": "abc123",
            },
        ],
        ["query_name", "rank", "fraction", "lineage", "query_md5"],
        delimiter=",",
    )
    write_lookup(
        lookup,
        [
            {"human_classification": "RG2", "canonical_taxid": "123"},
            {"human_classification": "RG3", "canonical_taxid": "125"},
        ],
    )

    run_sourmash_annotation(summary, bioboxes, lookup, output)

    rows = read_table(output, delimiter=",")
    assert list(rows[0]) == [
        "query_name",
        "rank",
        "fraction",
        "lineage",
        "who_risk_group",
        "query_md5",
    ]
    assert [row["who_risk_group"] for row in rows] == ["RG2", "RG2", "RG3", ""]


def test_sourmash_descendant_does_not_classify_ancestor_or_malformed_path(
    tmp_path: Path,
) -> None:
    summary = tmp_path / "sample.summarized.csv"
    lookup = tmp_path / "risk-groups.tsv"
    output = tmp_path / "sample.with-risk-groups.csv"
    bioboxes = tmp_path / "sample.bioboxes.profile"
    input_rows = [
        {"lineage": "Viruses;Example virus", "fraction": "0.5"},
        {
            "lineage": "Viruses;Example virus;Exceptional strain",
            "fraction": "0.25",
        },
        {"lineage": "Viruses;Malformed virus", "fraction": "0.25"},
    ]
    write_table(
        summary,
        input_rows,
        ["lineage", "fraction"],
        delimiter=",",
    )
    bioboxes.write_text(
        "# Taxonomic Profiling Output\n"
        "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n"
        "123\tspecies\t1|123\tViruses|Example virus\t50.0\n"
        "125\tstrain\t1|123|125\tViruses|Example virus|Exceptional strain\t25.0\n"
        "999\tspecies\t1|123\tViruses|Malformed virus\t25.0\n",
        encoding="utf-8",
    )
    write_lookup(
        lookup,
        [{"canonical_taxid": "125", "human_classification": "RG3"}],
    )

    run_sourmash_annotation(summary, bioboxes, lookup, output)

    rows = read_table(output, delimiter=",")
    assert [row["who_risk_group"] for row in rows] == ["", "RG3", ""]
    assert [
        {key: value for key, value in row.items() if key != "who_risk_group"}
        for row in rows
    ] == input_rows


def test_blast_results_add_risk_group_by_consensus_taxid(tmp_path: Path) -> None:
    blast = tmp_path / "sample.blast.tsv"
    lookup = tmp_path / "risk-groups.tsv"
    output = tmp_path / "sample.blast.with-risk-groups.tsv"
    taxonomy_dir = write_test_taxonomy(tmp_path / "taxonomy")
    input_rows = [
        {
            "task": "megablast",
            "sample": "sample-1",
            "adjusted_taxid": "123",
            "adjusted_taxid_name": "Example virus",
            "adjusted_taxid_rank": "species",
            "adjustment_method": "dominant",
        },
        {
            "task": "blastn",
            "sample": "sample-1",
            "adjusted_taxid": "not-a-taxid",
            "adjusted_taxid_name": "Malformed",
            "adjusted_taxid_rank": "no rank",
            "adjustment_method": "unresolved",
        },
        {
            "task": "blastn",
            "sample": "sample-1",
            "adjusted_taxid": "",
            "adjusted_taxid_name": "Missing",
            "adjusted_taxid_rank": "no rank",
            "adjustment_method": "unresolved",
        },
        {
            "task": "megablast",
            "sample": "sample-1",
            "adjusted_taxid": "999",
            "adjusted_taxid_name": "Unlisted virus",
            "adjusted_taxid_rank": "species",
            "adjustment_method": "dominant",
        },
    ]
    write_table(
        blast,
        input_rows,
        [
            "task",
            "sample",
            "adjusted_taxid",
            "adjusted_taxid_name",
            "adjusted_taxid_rank",
            "adjustment_method",
        ],
        delimiter="\t",
    )
    write_lookup(
        lookup,
        [
            {
                "human_classification": "RG2",
                "canonical_taxid": "123",
            },
        ],
    )

    run_blast_annotation(blast, lookup, taxonomy_dir, output)

    rows = read_table(output, delimiter="\t")
    assert list(rows[0]) == [
        "task",
        "sample",
        "adjusted_taxid",
        "adjusted_taxid_name",
        "adjusted_taxid_rank",
        "who_risk_group",
        "adjustment_method",
    ]
    assert [row["adjusted_taxid"] for row in rows] == [
        "123",
        "not-a-taxid",
        "",
        "999",
    ]
    assert [row["adjustment_method"] for row in rows] == [
        "dominant",
        "unresolved",
        "unresolved",
        "dominant",
    ]
    assert [row["who_risk_group"] for row in rows] == ["RG2", "", "", ""]
    assert [
        {key: value for key, value in row.items() if key != "who_risk_group"}
        for row in rows
    ] == input_rows


def test_conflicting_risk_groups_choose_highest_and_ignore_invalid_labels(
    tmp_path: Path,
) -> None:
    blast = tmp_path / "sample.blast.tsv"
    lookup = tmp_path / "risk-groups.tsv"
    output = tmp_path / "sample.blast.with-risk-groups.tsv"
    taxonomy_dir = write_test_taxonomy(tmp_path / "taxonomy")
    write_table(
        blast,
        [
            {
                "adjusted_taxid": "123",
                "adjusted_taxid_rank": "species",
            },
            {
                "adjusted_taxid": "999",
                "adjusted_taxid_rank": "species",
            },
        ],
        ["adjusted_taxid", "adjusted_taxid_rank"],
        delimiter="\t",
    )
    write_lookup(
        lookup,
        [
            {
                "canonical_taxid": "123",
                "human_classification": "RG1",
            },
            {
                "canonical_taxid": "123",
                "human_classification": "RG2",
            },
            {
                "canonical_taxid": "123",
                "human_classification": "RG3",
            },
            {
                "canonical_taxid": "123",
                "human_classification": "RG4",
            },
            {
                "canonical_taxid": "123",
                "human_classification": "RG5",
            },
            {
                "canonical_taxid": "999",
                "human_classification": "RG5",
            },
        ],
    )

    run_blast_annotation(blast, lookup, taxonomy_dir, output)

    rows = read_table(output, delimiter="\t")
    assert [row["who_risk_group"] for row in rows] == ["RG4", ""]


def test_nearest_classified_ancestor_is_inherited_and_descendant_overrides(
    tmp_path: Path,
) -> None:
    blast = tmp_path / "sample.blast.tsv"
    lookup = tmp_path / "risk-groups.tsv"
    output = tmp_path / "sample.blast.with-risk-groups.tsv"
    taxonomy_dir = write_test_taxonomy(tmp_path / "taxonomy")
    write_table(
        blast,
        [
            {"adjusted_taxid": "124", "adjusted_taxid_rank": "no rank"},
            {"adjusted_taxid": "125", "adjusted_taxid_rank": "no rank"},
        ],
        ["adjusted_taxid", "adjusted_taxid_rank"],
        delimiter="\t",
    )
    write_lookup(
        lookup,
        [
            {"canonical_taxid": "122", "human_classification": "RG1"},
            {"canonical_taxid": "123", "human_classification": "RG2"},
            {"canonical_taxid": "125", "human_classification": "RG3"},
        ],
    )

    run_blast_annotation(blast, lookup, taxonomy_dir, output)

    rows = read_table(output, delimiter="\t")
    assert [row["who_risk_group"] for row in rows] == ["RG2", "RG3"]


def test_descendant_classification_does_not_flow_to_ancestor(tmp_path: Path) -> None:
    blast = tmp_path / "sample.blast.tsv"
    lookup = tmp_path / "risk-groups.tsv"
    output = tmp_path / "sample.blast.with-risk-groups.tsv"
    taxonomy_dir = write_test_taxonomy(tmp_path / "taxonomy")
    write_table(
        blast,
        [{"adjusted_taxid": "123", "adjusted_taxid_rank": "species"}],
        ["adjusted_taxid", "adjusted_taxid_rank"],
        delimiter="\t",
    )
    write_lookup(
        lookup,
        [{"canonical_taxid": "125", "human_classification": "RG3"}],
    )

    run_blast_annotation(blast, lookup, taxonomy_dir, output)

    [row] = read_table(output, delimiter="\t")
    assert row["who_risk_group"] == ""


def test_unused_lookup_conflict_does_not_block_unique_match(tmp_path: Path) -> None:
    blast = tmp_path / "sample.blast.tsv"
    lookup = tmp_path / "risk-groups.tsv"
    output = tmp_path / "sample.blast.with-risk-groups.tsv"
    taxonomy_dir = write_test_taxonomy(tmp_path / "taxonomy")
    write_table(
        blast,
        [{"adjusted_taxid": "123", "adjusted_taxid_rank": "species"}],
        ["adjusted_taxid", "adjusted_taxid_rank"],
        delimiter="\t",
    )
    write_lookup(
        lookup,
        [
            {"canonical_taxid": "123", "human_classification": "RG2"},
            {"canonical_taxid": "999", "human_classification": "RG2"},
            {"canonical_taxid": "999", "human_classification": "RG3"},
        ],
    )

    run_blast_annotation(blast, lookup, taxonomy_dir, output)

    [row] = read_table(output, delimiter="\t")
    assert row["who_risk_group"] == "RG2"
