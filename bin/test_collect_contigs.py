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


def producer_run_rows(path: Path) -> list[dict[str, object]]:
    with sqlite3.connect(path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(
            """
            select
                sample_id,
                producer,
                n_contigs,
                n_short_assembly_contigs,
                n_long_assembly_contigs,
                long_contig_min_length
            from contig_producer_runs
            order by sample_id, producer
            """,
        ).fetchall()
    return [dict(row) for row in rows]


def fasta_headers(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]


def collect_args(
    *,
    input_fasta: Path,
    output_fasta: Path,
    contig_lookup: Path,
    producer: str = "spades",
    long_contig_min_length: int = 10_000,
) -> list[str]:
    return [
        "--sample-id",
        "sample-1",
        "--input-fasta",
        str(input_fasta),
        "--output-fasta",
        str(output_fasta),
        "--contig-lookup",
        str(contig_lookup),
        "--producer",
        producer,
        "--long-contig-min-length",
        str(long_contig_min_length),
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
        collect_args(
            input_fasta=input_fasta,
            output_fasta=output_fasta,
            contig_lookup=contig_lookup,
        ),
    )

    assert fasta_headers(output_fasta) == [
        ">nvdContig1_sample-1_000001 evidence_class=short_assembly_contig producer=spades contig_id=NODE_1_length_4_cov_1.0",
        ">nvdContig1_sample-1_000002 evidence_class=short_assembly_contig producer=spades contig_id=NODE_2_length_6_cov_2.0",
    ]
    assert lookup_rows(contig_lookup) == [
        {
            "qseqid": "nvdContig1_sample-1_000001",
            "sample_id": "sample-1",
            "evidence_class": "short_assembly_contig",
            "producer": "spades",
            "contig_id": "NODE_1_length_4_cov_1.0",
            "length": 4,
            "sha256": hashlib.sha256(b"ACGT").hexdigest(),
        },
        {
            "qseqid": "nvdContig1_sample-1_000002",
            "sample_id": "sample-1",
            "evidence_class": "short_assembly_contig",
            "producer": "spades",
            "contig_id": "NODE_2_length_6_cov_2.0",
            "length": 6,
            "sha256": hashlib.sha256(b"AACCGG").hexdigest(),
        },
    ]


def test_qseqid_does_not_include_length_based_class(tmp_path: Path) -> None:
    expected_qseqid = "nvdContig1_sample-1_000001"
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text(">NODE_1\nACGT\n", encoding="utf-8")
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    main(
        collect_args(
            input_fasta=input_fasta,
            output_fasta=output_fasta,
            contig_lookup=contig_lookup,
            long_contig_min_length=4,
        ),
    )

    [header] = fasta_headers(output_fasta)
    first_token = header[1:].split(maxsplit=1)[0]
    assert first_token == expected_qseqid
    assert "long_assembly_contig" not in first_token
    assert "evidence_class=long_assembly_contig" in header


def test_splits_short_and_long_assembly_contigs_at_threshold(tmp_path: Path) -> None:
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text(
        ">short\nAAA\n>long\nAAAA\n",
        encoding="utf-8",
    )
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    main(
        collect_args(
            input_fasta=input_fasta,
            output_fasta=output_fasta,
            contig_lookup=contig_lookup,
            producer="flye",
            long_contig_min_length=4,
        ),
    )

    rows = lookup_rows(contig_lookup)
    assert [row["evidence_class"] for row in rows] == [
        "short_assembly_contig",
        "long_assembly_contig",
    ]
    assert {row["producer"] for row in rows} == {"flye"}
    assert producer_run_rows(contig_lookup) == [
        {
            "sample_id": "sample-1",
            "producer": "flye",
            "n_contigs": 2,
            "n_short_assembly_contigs": 1,
            "n_long_assembly_contigs": 1,
            "long_contig_min_length": 4,
        },
    ]


def test_records_when_producer_emits_no_contigs(tmp_path: Path) -> None:
    input_fasta = tmp_path / "empty.fasta"
    input_fasta.write_text("", encoding="utf-8")
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    main(
        collect_args(
            input_fasta=input_fasta,
            output_fasta=output_fasta,
            contig_lookup=contig_lookup,
            long_contig_min_length=10_000,
        ),
    )

    assert output_fasta.read_text(encoding="utf-8") == ""
    assert lookup_rows(contig_lookup) == []
    assert producer_run_rows(contig_lookup) == [
        {
            "sample_id": "sample-1",
            "producer": "spades",
            "n_contigs": 0,
            "n_short_assembly_contigs": 0,
            "n_long_assembly_contigs": 0,
            "long_contig_min_length": 10_000,
        },
    ]


def test_rejects_duplicate_contig_ids(tmp_path: Path) -> None:
    input_fasta = tmp_path / "contigs.fasta"
    input_fasta.write_text(">NODE_1\nACGT\n>NODE_1\nTGCA\n", encoding="utf-8")
    output_fasta = tmp_path / "sample.contigs.fasta"
    contig_lookup = tmp_path / "sample.contigs.sqlite"

    with pytest.raises(ContigCollectionError, match="Duplicate FASTA identifier"):
        main(
            collect_args(
                input_fasta=input_fasta,
                output_fasta=output_fasta,
                contig_lookup=contig_lookup,
            ),
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
        collect_args(
            input_fasta=input_fasta,
            output_fasta=output_fasta,
            contig_lookup=contig_lookup,
        ),
    )

    with sqlite3.connect(contig_lookup) as connection:
        row = connection.execute(
            "select evidence_class, producer, contig_id from contigs where qseqid = ?",
            ("nvdContig1_sample-1_000001",),
        ).fetchone()
    assert row == ("short_assembly_contig", "spades", "NODE_1")
