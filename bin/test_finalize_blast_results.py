"""Tests for final BLAST result enrichment."""

from __future__ import annotations

import csv
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest
from finalize_blast_results import main

BIN = Path(__file__).resolve().parent
SCRIPT = BIN / "finalize_blast_results.py"

MERGED_HEADER = "task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids\trank"


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


def write_query_blast_tsv(path: Path, qseqids: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BLAST_HEADER, delimiter="\t")
        writer.writeheader()
        for qseqid in qseqids:
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
                query_class text not null,
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
                query_class,
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
                query_class text not null,
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
                query_class,
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
    assert row["query_class"] == "short_assembly_contig"
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
    assert row["query_class"] == ""
    assert row["producer"] == ""
    assert row["source_id"] == ""
    assert row["support_record_count"] == ""


def test_finalizer_loads_metadata_from_multiple_query_lookups(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "blast.tsv"
    counts = tmp_path / "counts.tsv"
    contig_lookup = tmp_path / "sample.contig.query_sequences.sqlite"
    read_lookup = tmp_path / "sample.single_read.query_sequences.sqlite"
    output = tmp_path / "final.tsv"
    write_query_blast_tsv(
        blast_tsv,
        (
            "nvdContigQuery_sample-1_000001",
            "nvdReadQuery_sample-1_000001",
        ),
    )
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
    assert read_row["query_class"] == "single_read"
    assert read_row["producer"] == "source_read"
    assert read_row["source_id"] == "readA/1"
    assert read_row["support_record_count"] == "3"
    assert read_row["mapped_reads"] == "0"


def test_rows_default_to_zero_without_contig_counts(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "blast.tsv"
    read_lookup = tmp_path / "sample.read.query_sequences.sqlite"
    output = tmp_path / "final.tsv"
    write_query_blast_tsv(blast_tsv, ("nvdReadQuery_sample-1_000001",))
    write_read_query_lookup(read_lookup)

    main(
        [
            "--blast-tsv",
            str(blast_tsv),
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

    assert {row["mapped_reads"] for row in read_tsv(output)} == {"0"}


def test_contig_rows_require_contig_counts(tmp_path: Path) -> None:
    blast_tsv = tmp_path / "blast.tsv"
    query_lookup = tmp_path / "sample.query_sequences.sqlite"
    output = tmp_path / "final.tsv"
    write_blast_tsv(blast_tsv)
    write_query_lookup(query_lookup)

    with pytest.raises(
        ValueError,
        match="contig query requires mapped-read counts: nvdContigQuery_sample-1_000001",
    ):
        main(
            [
                "--blast-tsv",
                str(blast_tsv),
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


def _lookup(path: Path, rows: list[tuple[str, str, str, str, int]]) -> None:
    conn = sqlite3.connect(path)
    conn.execute(
        "CREATE TABLE query_sequences (qseqid TEXT, query_class TEXT, producer TEXT, source_id TEXT, support_record_count INTEGER)",
    )
    conn.executemany("INSERT INTO query_sequences VALUES (?,?,?,?,?)", rows)
    conn.commit()
    conn.close()


def _run(blast_tsv: Path, out: Path, lookups: list[Path], contig_counts: Path | None) -> None:
    cmd = [sys.executable, str(SCRIPT), "--blast-tsv", str(blast_tsv), "--output", str(out),
           "--total-reads", "1000", "--blast-db-version", "v1", "--virus-index-version", "vx", "--run-id", "R1"]
    for lk in lookups:
        cmd += ["--query-lookup", str(lk)]
    if contig_counts is not None:
        cmd += ["--contig-counts", str(contig_counts)]
    subprocess.run(cmd, check=True)


def test_per_batch_concat_equals_per_sample(tmp_path: Path) -> None:
    """Per-batch enrichment must equal per-sample enrichment, row-for-row."""
    # Two batches for one sample: a contig batch and a single_read batch.
    contig_tsv = tmp_path / "contig.tsv"
    contig_tsv.write_text(MERGED_HEADER + "\nmegablast\ts1\tc1\t300\trefA\ttA\t250\t99.0\t1e-9\t400\tHomo\t9606\tspecies:x\n")
    read_tsv = tmp_path / "read.tsv"
    read_tsv.write_text(MERGED_HEADER + "\nmegablast\ts1\tnvdReadQuery_1\t150\trefB\ttB\t140\t98.0\t1e-8\t200\tHomo\t9606\tspecies:x\n")

    contig_lk = tmp_path / "contig.sqlite"
    _lookup(contig_lk, [("c1", "short_assembly_contig", "spades", "c1", 1)])
    read_lk = tmp_path / "read.sqlite"
    _lookup(read_lk, [("nvdReadQuery_1", "single_read", "source_read", "r1", 3)])

    counts = tmp_path / "counts.tsv"
    counts.write_text("c1\t42\n")

    # Per-batch: enrich each batch with its own lookup, concat.
    out_contig = tmp_path / "contig.final.tsv"
    out_read = tmp_path / "read.final.tsv"
    _run(contig_tsv, out_contig, [contig_lk, read_lk], counts)   # sample's grouped lookups; unused rows harmless
    _run(read_tsv, out_read, [contig_lk, read_lk], None)
    per_batch_lines = out_contig.read_text().splitlines() + out_read.read_text().splitlines()[1:]

    # Per-sample: concat batches first, then enrich once.
    merged = tmp_path / "merged.tsv"
    merged.write_text(contig_tsv.read_text() + "\n".join(read_tsv.read_text().splitlines()[1:]) + "\n")
    out_sample = tmp_path / "sample.final.tsv"
    _run(merged, out_sample, [contig_lk, read_lk], counts)
    per_sample_lines = out_sample.read_text().splitlines()

    assert sorted(per_batch_lines) == sorted(per_sample_lines)
