"""Tests for final BLAST result enrichment."""

from __future__ import annotations

import csv
import sqlite3
from typing import TYPE_CHECKING

from finalize_blast_results import main

if TYPE_CHECKING:
    from pathlib import Path


BLAST_HEADER = [
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
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_blast_tsv(path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BLAST_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "task": "megablast",
                "sample": "sample-1",
                "qseqid": "nvdContigQuery_sample-1_000001",
                "qlen": "4",
                "sseqid": "ref-1",
                "stitle": "reference",
                "length": "4",
                "pident": "100.0",
                "evalue": "1e-10",
                "bitscore": "42.0",
                "sscinames": "Virus example",
                "staxids": "10239",
                "rank": "superkingdom",
                "adjusted_taxid": "10239",
                "adjusted_taxid_name": "Viruses",
                "adjusted_taxid_rank": "superkingdom",
                "adjustment_method": "dominant",
            },
        )


def write_mixed_blast_tsv(path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BLAST_HEADER, delimiter="\t")
        writer.writeheader()
        for qseqid in (
            "nvdContigQuery_sample-1_000001",
            "nvdReadQuery_sample-1_000001",
        ):
            writer.writerow(
                {
                    "task": "megablast",
                    "sample": "sample-1",
                    "qseqid": qseqid,
                    "qlen": "4",
                    "sseqid": "ref-1",
                    "stitle": "reference",
                    "length": "4",
                    "pident": "100.0",
                    "evalue": "1e-10",
                    "bitscore": "42.0",
                    "sscinames": "Virus example",
                    "staxids": "10239",
                    "rank": "superkingdom",
                    "adjusted_taxid": "10239",
                    "adjusted_taxid_name": "Viruses",
                    "adjusted_taxid_rank": "superkingdom",
                    "adjustment_method": "dominant",
                },
            )


def write_counts(path: Path) -> None:
    path.write_text("nvdContigQuery_sample-1_000001\t7\n", encoding="utf-8")


def write_query_lookup(path: Path) -> None:
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            create table query_sequences (
                qseqid text primary key,
                sample_id text not null,
                evidence_class text not null,
                producer text not null,
                source_id text not null,
                support_record_count integer not null,
                length integer not null,
                sha256 text not null
            )
            """,
        )
        connection.execute(
            """
            insert into query_sequences (
                qseqid,
                sample_id,
                evidence_class,
                producer,
                source_id,
                support_record_count,
                length,
                sha256
            ) values (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "nvdContigQuery_sample-1_000001",
                "sample-1",
                "short_assembly_contig",
                "spades",
                "NODE_1_length_4_cov_1.0",
                1,
                4,
                "sha",
            ),
        )


def write_read_query_lookup(path: Path) -> None:
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            create table query_sequences (
                qseqid text primary key,
                sample_id text not null,
                evidence_class text not null,
                producer text not null,
                source_id text not null,
                support_record_count integer not null,
                length integer not null,
                sha256 text not null
            )
            """,
        )
        connection.execute(
            """
            insert into query_sequences (
                qseqid,
                sample_id,
                evidence_class,
                producer,
                source_id,
                support_record_count,
                length,
                sha256
            ) values (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "nvdReadQuery_sample-1_000001",
                "sample-1",
                "single_read",
                "source_read",
                "readA/1",
                3,
                4,
                "sha-read",
            ),
        )
def test_final_blast_rows_include_collected_contig_metadata(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "blast.tsv"
    counts = tmp_path / "counts.tsv"
    query_lookup = tmp_path / "sample.query_sequences.sqlite"
    output = tmp_path / "final.tsv"
    write_blast_tsv(blast_tsv)
    write_counts(counts)
    write_query_lookup(query_lookup)

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
            "--contig-counts",
            str(counts),
            "--query-lookup",
            str(query_lookup),
            "--output",
            str(output),
            "--total-reads",
            "10",
            "--blast-db-version",
            "nt-test",
            "--virus-index-version",
            "virus-test",
            "--run-id",
            "run-test",
        ],
    )

    [row] = read_tsv(output)
    assert row["qseqid"] == "nvdContigQuery_sample-1_000001"
    assert row["evidence_class"] == "short_assembly_contig"
    assert row["producer"] == "spades"
    assert row["source_id"] == "NODE_1_length_4_cov_1.0"
    assert row["support_record_count"] == "1"
    assert row["mapped_reads"] == "7"
    assert "evidence_length" not in row


def test_finalizer_can_run_without_query_lookup(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "blast.tsv"
    counts = tmp_path / "counts.tsv"
    output = tmp_path / "final.tsv"
    write_blast_tsv(blast_tsv)
    write_counts(counts)

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
            "--contig-counts",
            str(counts),
            "--output",
            str(output),
            "--total-reads",
            "10",
            "--blast-db-version",
            "nt-test",
            "--virus-index-version",
            "virus-test",
            "--run-id",
            "run-test",
        ],
    )

    [row] = read_tsv(output)
    assert row["evidence_class"] == ""
    assert row["producer"] == ""
    assert row["source_id"] == ""
    assert row["support_record_count"] == ""


def test_finalizer_loads_metadata_from_multiple_query_lookups(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "blast.tsv"
    counts = tmp_path / "counts.tsv"
    contig_lookup = tmp_path / "sample.contig.query_sequences.sqlite"
    read_lookup = tmp_path / "sample.single_read.query_sequences.sqlite"
    output = tmp_path / "final.tsv"
    write_mixed_blast_tsv(blast_tsv)
    write_counts(counts)
    write_query_lookup(contig_lookup)
    write_read_query_lookup(read_lookup)

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
            "--contig-counts",
            str(counts),
            "--query-lookup",
            str(contig_lookup),
            "--query-lookup",
            str(read_lookup),
            "--output",
            str(output),
            "--total-reads",
            "10",
            "--blast-db-version",
            "nt-test",
            "--virus-index-version",
            "virus-test",
            "--run-id",
            "run-test",
        ],
    )

    rows = read_tsv(output)
    read_row = next(row for row in rows if row["qseqid"].startswith("nvdReadQuery_"))
    assert read_row["evidence_class"] == "single_read"
    assert read_row["producer"] == "source_read"
    assert read_row["source_id"] == "readA/1"
    assert read_row["support_record_count"] == "3"
    assert read_row["mapped_reads"] == "0"
