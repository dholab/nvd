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
                "qseqid": "nvdContig1_sample-1_000001",
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
    path.write_text("nvdContig1_sample-1_000001\t7\n", encoding="utf-8")


def write_contig_lookup(path: Path) -> None:
    with sqlite3.connect(path) as connection:
        connection.execute(
            """
            create table contigs (
                qseqid text primary key,
                sample_id text not null,
                evidence_class text not null,
                producer text not null,
                contig_id text not null,
                length integer not null,
                sha256 text not null
            )
            """,
        )
        connection.execute(
            """
            insert into contigs (
                qseqid,
                sample_id,
                evidence_class,
                producer,
                contig_id,
                length,
                sha256
            ) values (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "nvdContig1_sample-1_000001",
                "sample-1",
                "short_assembly_contig",
                "spades",
                "NODE_1_length_4_cov_1.0",
                4,
                "sha",
            ),
        )


def test_final_blast_rows_include_collected_contig_metadata(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "blast.tsv"
    counts = tmp_path / "counts.tsv"
    contig_lookup = tmp_path / "sample.contigs.sqlite"
    output = tmp_path / "final.tsv"
    write_blast_tsv(blast_tsv)
    write_counts(counts)
    write_contig_lookup(contig_lookup)

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
            "--contig-counts",
            str(counts),
            "--contig-lookup",
            str(contig_lookup),
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
    assert row["qseqid"] == "nvdContig1_sample-1_000001"
    assert row["evidence_class"] == "short_assembly_contig"
    assert row["producer"] == "spades"
    assert row["contig_id"] == "NODE_1_length_4_cov_1.0"
    assert row["mapped_reads"] == "7"
    assert "evidence_length" not in row


def test_finalizer_can_run_without_contig_lookup(tmp_path: Path) -> None:
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
    assert row["contig_id"] == ""
