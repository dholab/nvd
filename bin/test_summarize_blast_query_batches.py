from __future__ import annotations

import csv
from pathlib import Path

import pytest
from summarize_blast_query_batches import QueryBatchSummaryError, main


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_summarizes_present_and_absent_query_classes(tmp_path: Path) -> None:
    short_fasta = tmp_path / "sample-1.short_assembly_contig.blast_queries.fasta"
    short_fasta.write_text(
        ">q1 query_class=short_assembly_contig\nACGT\n"
        ">q2 query_class=short_assembly_contig\nTGCA\n",
        encoding="utf-8",
    )
    short_lookup = tmp_path / "sample-1.query_sequences.sqlite"
    short_lookup.write_bytes(b"sqlite placeholder")
    single_fasta = tmp_path / "sample-1.single_read.blast_queries.fasta"
    single_fasta.write_text(">r1 query_class=single_read\nAAAA\n", encoding="utf-8")
    single_lookup = tmp_path / "sample-1.single_read.query_sequences.sqlite"
    single_lookup.write_bytes(b"sqlite placeholder")
    output = tmp_path / "sample-1.blast_query_batches.tsv"

    main(
        [
            "--sample-id",
            "sample-1",
            "--platform",
            "illumina",
            "--batch",
            "short_assembly_contig",
            str(short_fasta),
            str(short_lookup),
            "--batch",
            "single_read",
            str(single_fasta),
            str(single_lookup),
            "--output",
            str(output),
        ],
    )

    rows = read_tsv(output)
    by_class = {row["query_class"]: row for row in rows}
    assert [row["query_class"] for row in rows] == [
        "short_assembly_contig",
        "long_assembly_contig",
        "overlap_merged_pair",
        "single_read",
    ]
    assert by_class["short_assembly_contig"]["n_query_sequences"] == "2"
    assert by_class["short_assembly_contig"]["query_fasta_present"] == "true"
    assert by_class["short_assembly_contig"]["query_lookup_present"] == "true"
    assert by_class["long_assembly_contig"]["n_query_sequences"] == "0"
    assert by_class["long_assembly_contig"]["query_fasta_present"] == "false"
    assert by_class["long_assembly_contig"]["query_lookup"] == ""
    assert by_class["overlap_merged_pair"]["query_source"] == "read_query"
    assert by_class["single_read"]["n_query_sequences"] == "1"


def test_rejects_duplicate_query_class_batches(tmp_path: Path) -> None:
    fasta = tmp_path / "sample-1.single_read.blast_queries.fasta"
    fasta.write_text(">r1\nACGT\n", encoding="utf-8")
    lookup = tmp_path / "sample-1.single_read.query_sequences.sqlite"
    lookup.write_bytes(b"sqlite placeholder")

    with pytest.raises(QueryBatchSummaryError, match="duplicate BLAST query batch"):
        main(
            [
                "--sample-id",
                "sample-1",
                "--platform",
                "illumina",
                "--batch",
                "single_read",
                str(fasta),
                str(lookup),
                "--batch",
                "single_read",
                str(fasta),
                str(lookup),
                "--output",
                str(tmp_path / "summary.tsv"),
            ],
        )


def test_rejects_unknown_query_class_batches(tmp_path: Path) -> None:
    fasta = tmp_path / "sample-1.unknown.blast_queries.fasta"
    fasta.write_text(">r1\nACGT\n", encoding="utf-8")
    lookup = tmp_path / "sample-1.unknown.query_sequences.sqlite"
    lookup.write_bytes(b"sqlite placeholder")

    with pytest.raises(QueryBatchSummaryError, match="unsupported BLAST query classes"):
        main(
            [
                "--sample-id",
                "sample-1",
                "--platform",
                "illumina",
                "--batch",
                "unknown",
                str(fasta),
                str(lookup),
                "--output",
                str(tmp_path / "summary.tsv"),
            ],
        )


def test_no_batches_emits_explicit_zero_rows_for_every_query_class(
    tmp_path: Path,
) -> None:
    output = tmp_path / "sample-1.blast_query_batches.tsv"

    main(
        [
            "--sample-id",
            "sample-1",
            "--platform",
            "illumina",
            "--output",
            str(output),
        ],
    )

    rows = read_tsv(output)
    assert [row["query_class"] for row in rows] == [
        "short_assembly_contig",
        "long_assembly_contig",
        "overlap_merged_pair",
        "single_read",
    ]
    assert all(row["n_query_sequences"] == "0" for row in rows)
    assert all(row["query_fasta_present"] == "false" for row in rows)
    assert all(row["query_lookup_present"] == "false" for row in rows)
