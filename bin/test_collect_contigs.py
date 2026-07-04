"""Tests for collecting contigs under stable NVD query IDs."""

from __future__ import annotations

import hashlib
import sqlite3
from typing import TYPE_CHECKING

import pytest
from collect_contigs import ContigCollectionError, main, parse_fasta

if TYPE_CHECKING:
    from pathlib import Path


def lookup_rows(path: Path) -> list[dict[str, object]]:
    with sqlite3.connect(path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(
            """
            select
                qseqid,
                sample_id,
                evidence_class,
                producer,
                contig_id,
                length,
                sha256
            from contigs
            order by qseqid
            """,
        ).fetchall()
    return [dict(row) for row in rows]


def fasta_headers(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]


def test_collects_assembly_contigs_with_stable_ids_and_lookup(
    tmp_path: Path,
) -> None:
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text(
        ">NODE_1_length_4_cov_1.0\nACGT\n"
        ">NODE_2_length_6_cov_2.0 extra assembler text\nAACCGG\n",
        encoding="utf-8",
    )
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    main(
        [
            "--sample-id",
            "sample-1",
            "--input-fasta",
            str(input_fasta),
            "--output-fasta",
            str(output_fasta),
            "--contig-lookup",
            str(contig_lookup),
        ],
    )

    assert fasta_headers(output_fasta) == [
        ">nvdContig1_sample-1_000001 evidence_class=assembly_contig producer=spades contig_id=NODE_1_length_4_cov_1.0",
        ">nvdContig1_sample-1_000002 evidence_class=assembly_contig producer=spades contig_id=NODE_2_length_6_cov_2.0",
    ]
    assert lookup_rows(contig_lookup) == [
        {
            "qseqid": "nvdContig1_sample-1_000001",
            "sample_id": "sample-1",
            "evidence_class": "assembly_contig",
            "producer": "spades",
            "contig_id": "NODE_1_length_4_cov_1.0",
            "length": 4,
            "sha256": hashlib.sha256(b"ACGT").hexdigest(),
        },
        {
            "qseqid": "nvdContig1_sample-1_000002",
            "sample_id": "sample-1",
            "evidence_class": "assembly_contig",
            "producer": "spades",
            "contig_id": "NODE_2_length_6_cov_2.0",
            "length": 6,
            "sha256": hashlib.sha256(b"AACCGG").hexdigest(),
        },
    ]


def test_qseqid_does_not_include_mutable_class(tmp_path: Path) -> None:
    expected_qseqid = "nvdContig1_sample-1_000001"
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text(">NODE_1\nACGT\n", encoding="utf-8")
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    main(
        [
            "--sample-id",
            "sample-1",
            "--input-fasta",
            str(input_fasta),
            "--output-fasta",
            str(output_fasta),
            "--contig-lookup",
            str(contig_lookup),
            "--evidence-class",
            "genome_like_assembly_contig",
        ],
    )

    [header] = fasta_headers(output_fasta)
    first_token = header[1:].split(maxsplit=1)[0]
    assert first_token == expected_qseqid
    assert "genome_like_assembly_contig" not in first_token
    assert "evidence_class=genome_like_assembly_contig" in header


def test_rejects_duplicate_contig_ids(tmp_path: Path) -> None:
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text(">NODE_1\nACGT\n>NODE_1\nTGCA\n", encoding="utf-8")
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    with pytest.raises(ContigCollectionError, match="Duplicate FASTA identifier"):
        main(
            [
                "--sample-id",
                "sample-1",
                "--input-fasta",
                str(input_fasta),
                "--output-fasta",
                str(output_fasta),
                "--contig-lookup",
                str(contig_lookup),
            ],
        )


def test_rejects_sequence_before_header(tmp_path: Path) -> None:
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text("ACGT\n", encoding="utf-8")

    with (
        input_fasta.open(encoding="utf-8") as handle,
        pytest.raises(
            ContigCollectionError,
            match="before any header",
        ),
    ):
        parse_fasta(handle)


def test_lookup_is_queryable_by_qseqid(tmp_path: Path) -> None:
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text(">NODE_1\nACGT\n", encoding="utf-8")
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    main(
        [
            "--sample-id",
            "sample-1",
            "--input-fasta",
            str(input_fasta),
            "--output-fasta",
            str(output_fasta),
            "--contig-lookup",
            str(contig_lookup),
        ],
    )

    with sqlite3.connect(contig_lookup) as connection:
        row = connection.execute(
            "select evidence_class, producer, contig_id from contigs where qseqid = ?",
            ("nvdContig1_sample-1_000001",),
        ).fetchone()
    assert row == ("assembly_contig", "spades", "NODE_1")
