"""Tests for selecting BLAST query FASTA records."""

from __future__ import annotations

import sqlite3
from typing import TYPE_CHECKING

import pytest
from select_blast_queries import BlastQuerySelectionError, main, parse_fasta

if TYPE_CHECKING:
    from pathlib import Path


def write_lookup(path: Path, rows: list[tuple[str, str]]) -> None:
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
        connection.executemany(
            """
            insert into contigs (
                qseqid,
                sample_id,
                evidence_class,
                producer,
                contig_id,
                length,
                sha256
            ) values (?, 'sample-1', ?, 'spades', 'NODE_1', 4, 'unused')
            """,
            rows,
        )


def test_selects_current_assembly_contig_classes_without_mutating_lookup(
    tmp_path: Path,
) -> None:
    input_fasta = tmp_path / "screened.fasta"
    input_fasta.write_text(
        ">nvdContig1_sample-1_000001 evidence_class=short_assembly_contig\nACGT\n"
        ">nvdContig1_sample-1_000002 evidence_class=long_assembly_contig\nAACCGG\n",
        encoding="utf-8",
    )
    lookup = tmp_path / "sample.contigs.sqlite"
    write_lookup(
        lookup,
        [
            ("nvdContig1_sample-1_000001", "short_assembly_contig"),
            ("nvdContig1_sample-1_000002", "long_assembly_contig"),
        ],
    )
    before = lookup.read_bytes()
    output_fasta = tmp_path / "selected.fasta"

    main(
        [
            "--input-fasta",
            str(input_fasta),
            "--contig-lookup",
            str(lookup),
            "--output-fasta",
            str(output_fasta),
        ],
    )

    assert output_fasta.read_text(encoding="utf-8") == (
        ">nvdContig1_sample-1_000001 evidence_class=short_assembly_contig\nACGT\n"
        ">nvdContig1_sample-1_000002 evidence_class=long_assembly_contig\nAACCGG\n"
    )
    assert lookup.read_bytes() == before


def test_empty_screened_fasta_produces_empty_selected_fasta(tmp_path: Path) -> None:
    input_fasta = tmp_path / "screened.fasta"
    input_fasta.write_text("", encoding="utf-8")
    lookup = tmp_path / "sample.contigs.sqlite"
    write_lookup(lookup, [])
    output_fasta = tmp_path / "selected.fasta"

    main(
        [
            "--input-fasta",
            str(input_fasta),
            "--contig-lookup",
            str(lookup),
            "--output-fasta",
            str(output_fasta),
        ],
    )

    assert output_fasta.read_text(encoding="utf-8") == ""


def test_missing_lookup_row_fails_loudly(tmp_path: Path) -> None:
    input_fasta = tmp_path / "screened.fasta"
    input_fasta.write_text(">nvdContig1_sample-1_000001\nACGT\n", encoding="utf-8")
    lookup = tmp_path / "sample.contigs.sqlite"
    write_lookup(lookup, [])
    output_fasta = tmp_path / "selected.fasta"

    with pytest.raises(
        BlastQuerySelectionError,
        match="missing from the contig lookup",
    ):
        main(
            [
                "--input-fasta",
                str(input_fasta),
                "--contig-lookup",
                str(lookup),
                "--output-fasta",
                str(output_fasta),
            ],
        )


def test_unadmitted_class_fails_loudly(tmp_path: Path) -> None:
    input_fasta = tmp_path / "screened.fasta"
    input_fasta.write_text(">nvdContig1_sample-1_000001\nACGT\n", encoding="utf-8")
    lookup = tmp_path / "sample.contigs.sqlite"
    write_lookup(lookup, [("nvdContig1_sample-1_000001", "overlap_merged_pair")])
    output_fasta = tmp_path / "selected.fasta"

    with pytest.raises(BlastQuerySelectionError, match="admits only"):
        main(
            [
                "--input-fasta",
                str(input_fasta),
                "--contig-lookup",
                str(lookup),
                "--output-fasta",
                str(output_fasta),
            ],
        )


def test_duplicate_qseqid_fails_loudly(tmp_path: Path) -> None:
    input_fasta = tmp_path / "screened.fasta"
    input_fasta.write_text(
        ">nvdContig1_sample-1_000001\nACGT\n>nvdContig1_sample-1_000001 duplicate\nTGCA\n",
        encoding="utf-8",
    )
    lookup = tmp_path / "sample.contigs.sqlite"
    write_lookup(lookup, [("nvdContig1_sample-1_000001", "short_assembly_contig")])
    output_fasta = tmp_path / "selected.fasta"

    with pytest.raises(BlastQuerySelectionError, match="Duplicate FASTA qseqid"):
        main(
            [
                "--input-fasta",
                str(input_fasta),
                "--contig-lookup",
                str(lookup),
                "--output-fasta",
                str(output_fasta),
            ],
        )


def test_parse_fasta_rejects_sequence_before_header() -> None:
    with pytest.raises(BlastQuerySelectionError, match="before any header"):
        parse_fasta(["ACGT\n"])
