#!/usr/bin/env python3
"""Collect SPAdes contigs under stable NVD query IDs."""

from __future__ import annotations

import argparse
import hashlib
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
WRAP_WIDTH = 80


class ContigCollectionError(Exception):
    """Raised when FASTA records cannot be collected as contigs."""


@dataclass(frozen=True)
class FastaRecord:
    """A parsed FASTA record with its original identifier and sequence."""

    identifier: str
    sequence: str


@dataclass(frozen=True)
class CollectedContig:
    """A contig collected under a stable NVD query ID."""

    qseqid: str
    sample_id: str
    evidence_class: str
    producer: str
    contig_id: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode()).hexdigest()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect SPAdes contigs under stable NVD query IDs.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--contig-lookup", required=True, type=Path)
    parser.add_argument("--producer", default="spades")
    parser.add_argument("--evidence-class", default="assembly_contig")
    return parser.parse_args(argv)


def parse_fasta(handle: TextIO) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_identifier: str | None = None
    current_sequence: list[str] = []

    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_identifier is not None:
                records.append(
                    FastaRecord(current_identifier, "".join(current_sequence)),
                )
            header = line[1:].strip()
            if not header:
                msg = f"FASTA header on line {line_number} is empty."
                raise ContigCollectionError(msg)
            current_identifier = header.split(maxsplit=1)[0]
            current_sequence = []
        elif current_identifier is None:
            msg = f"FASTA sequence line {line_number} appears before any header."
            raise ContigCollectionError(msg)
        else:
            current_sequence.append(line)

    if current_identifier is not None:
        records.append(FastaRecord(current_identifier, "".join(current_sequence)))

    return records


def validate_inputs(
    *,
    sample_id: str,
    evidence_class: str,
    producer: str,
) -> None:
    for label, value in (
        ("sample_id", sample_id),
        ("evidence_class", evidence_class),
        ("producer", producer),
    ):
        if not value:
            msg = f"{label} must not be empty."
            raise ContigCollectionError(msg)
        if not SAMPLE_ID_RE.match(value):
            msg = f"{label} contains characters unsafe for FASTA key-value labels: {value!r}."
            raise ContigCollectionError(msg)


def collect_contigs(
    records: list[FastaRecord],
    *,
    sample_id: str,
    evidence_class: str,
    producer: str,
) -> list[CollectedContig]:
    validate_inputs(
        sample_id=sample_id,
        evidence_class=evidence_class,
        producer=producer,
    )

    seen_contigs: set[str] = set()
    contigs: list[CollectedContig] = []
    for index, record in enumerate(records, start=1):
        if not record.sequence:
            msg = f"FASTA record {record.identifier!r} has no sequence."
            raise ContigCollectionError(msg)
        if record.identifier in seen_contigs:
            msg = f"Duplicate FASTA identifier: {record.identifier!r}."
            raise ContigCollectionError(msg)
        seen_contigs.add(record.identifier)
        contigs.append(
            CollectedContig(
                qseqid=f"nvdContig1_{sample_id}_{index:06d}",
                sample_id=sample_id,
                evidence_class=evidence_class,
                producer=producer,
                contig_id=record.identifier,
                sequence=record.sequence,
            ),
        )
    return contigs


def wrapped(sequence: str) -> str:
    return "\n".join(
        sequence[start : start + WRAP_WIDTH]
        for start in range(0, len(sequence), WRAP_WIDTH)
    )


def write_fasta(path: Path, contigs: list[CollectedContig]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in contigs:
            handle.write(
                f">{record.qseqid} "
                f"evidence_class={record.evidence_class} "
                f"producer={record.producer} "
                f"contig_id={record.contig_id}\n",
            )
            handle.write(f"{wrapped(record.sequence)}\n")


def write_contig_lookup(path: Path, contigs: list[CollectedContig]) -> None:
    if path.exists():
        path.unlink()

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
            ) values (?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    record.qseqid,
                    record.sample_id,
                    record.evidence_class,
                    record.producer,
                    record.contig_id,
                    record.length,
                    record.sha256,
                )
                for record in contigs
            ],
        )
        connection.execute("create index contigs_sample_id on contigs(sample_id)")
        connection.execute(
            "create unique index contigs_contig_id on contigs(contig_id)",
        )


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    with args.input_fasta.open(encoding="utf-8") as handle:
        fasta_records = parse_fasta(handle)
    contigs = collect_contigs(
        fasta_records,
        sample_id=args.sample_id,
        evidence_class=args.evidence_class,
        producer=args.producer,
    )
    write_fasta(args.output_fasta, contigs)
    write_contig_lookup(args.contig_lookup, contigs)


if __name__ == "__main__":
    main()
