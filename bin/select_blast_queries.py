#!/usr/bin/env python3
"""Select FASTA records admitted to BLAST by class-aware rules."""

from __future__ import annotations

import argparse
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO
from urllib.parse import quote

WRAP_WIDTH = 80
ADMITTED_CLASSES = frozenset({"short_assembly_contig", "long_assembly_contig"})


class BlastQuerySelectionError(Exception):
    """Raised when BLAST query selection cannot validate candidate records."""


@dataclass(frozen=True)
class FastaRecord:
    """A FASTA record preserving the full header text."""

    header: str
    sequence: str

    @property
    def qseqid(self) -> str:
        return self.header.split(maxsplit=1)[0]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select FASTA records admitted to BLAST by class-aware rules.",
    )
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--query-lookup", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    return parser.parse_args(argv)


def parse_fasta(handle: TextIO) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    current_header: str | None = None
    current_sequence: list[str] = []

    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append(FastaRecord(current_header, "".join(current_sequence)))
            current_header = line[1:].strip()
            if not current_header:
                msg = f"FASTA header on line {line_number} is empty."
                raise BlastQuerySelectionError(msg)
            current_sequence = []
        elif current_header is None:
            msg = f"FASTA sequence line {line_number} appears before any header."
            raise BlastQuerySelectionError(msg)
        else:
            current_sequence.append(line)

    if current_header is not None:
        records.append(FastaRecord(current_header, "".join(current_sequence)))

    return records


def sqlite_readonly_uri(path: Path) -> str:
    quoted_path = quote(str(path.resolve()), safe="/")
    return f"file:{quoted_path}?mode=ro"


def load_classes_by_qseqid(path: Path) -> dict[str, str]:
    with sqlite3.connect(sqlite_readonly_uri(path), uri=True) as connection:
        rows = connection.execute(
            """
            select qseqid, evidence_class
            from query_sequences
            """,
        ).fetchall()
    return dict(rows)


def select_blast_queries(
    records: list[FastaRecord],
    *,
    classes_by_qseqid: dict[str, str],
) -> list[FastaRecord]:
    selected: list[FastaRecord] = []
    seen_qseqids: set[str] = set()
    for record in records:
        if not record.sequence:
            msg = f"FASTA record {record.qseqid!r} has no sequence."
            raise BlastQuerySelectionError(msg)
        if record.qseqid in seen_qseqids:
            msg = f"Duplicate FASTA qseqid: {record.qseqid!r}."
            raise BlastQuerySelectionError(msg)
        seen_qseqids.add(record.qseqid)

        evidence_class = classes_by_qseqid.get(record.qseqid)
        if evidence_class is None:
            msg = f"FASTA qseqid {record.qseqid!r} is missing from the query lookup."
            raise BlastQuerySelectionError(msg)
        if evidence_class not in ADMITTED_CLASSES:
            admitted = ", ".join(sorted(ADMITTED_CLASSES))
            msg = (
                f"FASTA qseqid {record.qseqid!r} has class {evidence_class!r}; "
                f"this selector currently admits only: {admitted}."
            )
            raise BlastQuerySelectionError(msg)
        selected.append(record)
    return selected


def wrapped(sequence: str) -> str:
    return "\n".join(
        sequence[start : start + WRAP_WIDTH]
        for start in range(0, len(sequence), WRAP_WIDTH)
    )


def write_fasta(path: Path, records: list[FastaRecord]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            handle.write(f"{wrapped(record.sequence)}\n")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    with args.input_fasta.open(encoding="utf-8") as handle:
        records = parse_fasta(handle)
    classes_by_qseqid = load_classes_by_qseqid(args.query_lookup)
    selected = select_blast_queries(records, classes_by_qseqid=classes_by_qseqid)
    write_fasta(args.output_fasta, selected)


if __name__ == "__main__":
    main()
