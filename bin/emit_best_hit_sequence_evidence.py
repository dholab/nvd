#!/usr/bin/env python3
"""Emit sequence evidence for Query Big Table representative best hits."""

# ruff: noqa: EM101, EM102, TRY003

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator


OUTPUT_FILENAMES = (
    "query_sequences.fasta",
    "selected_references.fasta",
    "best_hit_placements.bed",
)
BLASTDBCMD_FIELD_COUNT = 2
REQUIRED_QBT_COLUMNS = {
    "sample_id",
    "qseqid",
    "qlen",
    "best_hit_reference_accession",
    "best_hit_reference_title",
    "best_hit_query_start_1based",
    "best_hit_query_end_1based",
    "best_hit_reference_length",
    "best_hit_reference_start_1based",
    "best_hit_reference_end_1based",
    "best_hit_reference_strand",
}
STRANDS = {"plus": "+", "minus": "-", "+": "+", "-": "-"}


class BestHitSequenceEvidenceError(Exception):
    """Raised when the evidence artifacts would be incomplete or inconsistent."""


@dataclass(frozen=True)
class Placement:
    """One Query Big Table representative best-HSP envelope."""

    sample_id: str
    qseqid: str
    qlen: int
    reference_accession: str
    reference_title: str
    reference_length: int
    reference_start: int
    reference_end: int
    strand: str

    @property
    def query_id(self) -> str:
        return f"{self.sample_id}|{self.qseqid}"

    @property
    def bed_start(self) -> int:
        return min(self.reference_start, self.reference_end) - 1

    @property
    def bed_end(self) -> int:
        return max(self.reference_start, self.reference_end)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Emit FASTA and BED evidence for Query Big Table best hits.",
    )
    parser.add_argument("--query-big-table", required=True, type=Path)
    parser.add_argument(
        "--query-fasta",
        action="append",
        nargs=2,
        metavar=("SAMPLE_ID", "FASTA"),
        default=[],
    )
    parser.add_argument("--blast-db", required=True, type=Path)
    parser.add_argument("--blast-db-prefix", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--blastdbcmd-bin", default="blastdbcmd")
    return parser.parse_args(argv)


def fasta_records(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    """Yield the first header token and unwrapped sequence from FASTA text."""
    identifier: str | None = None
    sequence: list[str] = []
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                yield identifier, "".join(sequence)
            header = line[1:].strip()
            if not header:
                raise BestHitSequenceEvidenceError(
                    f"FASTA header on line {line_number} is empty.",
                )
            identifier = header.split(maxsplit=1)[0]
            sequence = []
        elif identifier is None:
            raise BestHitSequenceEvidenceError(
                f"FASTA sequence line {line_number} appears before any header.",
            )
        else:
            sequence.append(line)
    if identifier is not None:
        yield identifier, "".join(sequence)


def positive_int(value: str, column: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise BestHitSequenceEvidenceError(
            f"QBT column {column!r} must be an integer, got {value!r}.",
        ) from error
    if result < 1:
        raise BestHitSequenceEvidenceError(
            f"QBT column {column!r} must be positive, got {value!r}.",
        )
    return result


def load_placements(path: Path) -> list[Placement]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(REQUIRED_QBT_COLUMNS - set(reader.fieldnames or []))
        if missing:
            raise BestHitSequenceEvidenceError(
                f"Query Big Table missing required columns: {', '.join(missing)}",
            )

        placements: list[Placement] = []
        queries: set[tuple[str, str]] = set()
        references: dict[str, tuple[str, int]] = {}
        for line_number, row in enumerate(reader, start=2):
            query = (row["sample_id"], row["qseqid"])
            if query in queries:
                raise BestHitSequenceEvidenceError(
                    "Query Big Table must contain one row per query; duplicate "
                    f"row for {query[0]}|{query[1]}.",
                )
            queries.add(query)

            qlen = positive_int(row["qlen"], "qlen")
            query_start = positive_int(
                row["best_hit_query_start_1based"],
                "best_hit_query_start_1based",
            )
            query_end = positive_int(
                row["best_hit_query_end_1based"],
                "best_hit_query_end_1based",
            )
            reference_length = positive_int(
                row["best_hit_reference_length"],
                "best_hit_reference_length",
            )
            reference_start = positive_int(
                row["best_hit_reference_start_1based"],
                "best_hit_reference_start_1based",
            )
            reference_end = positive_int(
                row["best_hit_reference_end_1based"],
                "best_hit_reference_end_1based",
            )
            if max(query_start, query_end) > qlen:
                raise BestHitSequenceEvidenceError(
                    f"Query placement on line {line_number} exceeds qlen {qlen}.",
                )
            if max(reference_start, reference_end) > reference_length:
                raise BestHitSequenceEvidenceError(
                    f"Reference placement on line {line_number} exceeds "
                    f"best_hit_reference_length {reference_length}.",
                )

            accession = row["best_hit_reference_accession"]
            title = row["best_hit_reference_title"].strip()
            metadata = (title, reference_length)
            previous = references.get(accession)
            if previous is not None and previous != metadata:
                raise BestHitSequenceEvidenceError(
                    f"Conflicting metadata for reference {accession!r}.",
                )
            references[accession] = metadata

            strand = STRANDS.get(row["best_hit_reference_strand"])
            if strand is None:
                raise BestHitSequenceEvidenceError(
                    f"Invalid reference strand on line {line_number}.",
                )
            placements.append(
                Placement(
                    sample_id=query[0],
                    qseqid=query[1],
                    qlen=qlen,
                    reference_accession=accession,
                    reference_title=title,
                    reference_length=reference_length,
                    reference_start=reference_start,
                    reference_end=reference_end,
                    strand=strand,
                ),
            )
    return placements


def load_query_sequences(
    query_fastas: list[list[str]],
    placements: list[Placement],
) -> dict[tuple[str, str], str]:
    expected = {
        (placement.sample_id, placement.qseqid): placement.qlen
        for placement in placements
    }
    sequences: dict[tuple[str, str], str] = {}
    for sample_id, fasta_text in query_fastas:
        fasta_path = Path(fasta_text)
        with fasta_path.open(encoding="utf-8") as handle:
            for qseqid, sequence in fasta_records(handle):
                key = (sample_id, qseqid)
                if key not in expected:
                    continue
                if key in sequences:
                    raise BestHitSequenceEvidenceError(
                        f"Duplicate query sequence for {sample_id}|{qseqid}.",
                    )
                if len(sequence) != expected[key]:
                    raise BestHitSequenceEvidenceError(
                        f"Query sequence length for {sample_id}|{qseqid} is "
                        f"{len(sequence)}, expected qlen {expected[key]}.",
                    )
                sequences[key] = sequence

    missing = sorted(set(expected) - set(sequences))
    if missing:
        sample_id, qseqid = missing[0]
        raise BestHitSequenceEvidenceError(
            f"Missing query sequence for {sample_id}|{qseqid}.",
        )
    return sequences


def write_query_fasta(
    path: Path,
    sequences: dict[tuple[str, str], str],
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for sample_id, qseqid in sorted(sequences):
            handle.write(f">{sample_id}|{qseqid}\n{sequences[(sample_id, qseqid)]}\n")


def write_reference_fasta(
    path: Path,
    placements: list[Placement],
    args: argparse.Namespace,
) -> None:
    references = {
        placement.reference_accession: (
            placement.reference_title,
            placement.reference_length,
        )
        for placement in placements
    }
    accessions = sorted(references)
    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary = Path(temporary_directory)
        batch_path = temporary / "accessions.txt"
        raw_path = temporary / "references.tsv"
        batch_path.write_text("\n".join(accessions) + "\n", encoding="utf-8")
        command = [
            args.blastdbcmd_bin,
            "-db",
            str(args.blast_db / args.blast_db_prefix),
            "-entry_batch",
            str(batch_path),
            "-target_only",
            "-outfmt",
            "%a\t%s",
        ]
        try:
            with raw_path.open("w", encoding="utf-8") as stdout:
                completed = subprocess.run(  # noqa: S603
                    command,
                    check=False,
                    stdout=stdout,
                    stderr=subprocess.PIPE,
                    text=True,
                )
        except OSError as error:
            raise BestHitSequenceEvidenceError(
                f"Failed to launch blastdbcmd: {error}",
            ) from error
        if completed.returncode != 0:
            raise BestHitSequenceEvidenceError(
                f"blastdbcmd failed with exit status {completed.returncode}: {(completed.stderr or '').strip()}",
            )

        seen: set[str] = set()
        with (
            path.open("w", encoding="utf-8") as output,
            raw_path.open(
                encoding="utf-8",
            ) as raw,
        ):
            for index, line in enumerate(raw):
                if index >= len(accessions):
                    raise BestHitSequenceEvidenceError(
                        "blastdbcmd returned more references than requested.",
                    )
                expected_accession = accessions[index]
                fields = line.rstrip("\n").split("\t")
                if (
                    len(fields) != BLASTDBCMD_FIELD_COUNT
                    or fields[0] != expected_accession
                    or not fields[1]
                ):
                    raise BestHitSequenceEvidenceError(
                        f"Invalid blastdbcmd result for {expected_accession!r}.",
                    )
                accession, sequence = fields
                title, expected_length = references[accession]
                if len(sequence) != expected_length:
                    raise BestHitSequenceEvidenceError(
                        f"Reference sequence length for {accession} is "
                        f"{len(sequence)}, expected best_hit_reference_length "
                        f"{expected_length}.",
                    )
                header = f"{accession} {title}" if title else accession
                output.write(f">{header}\n{sequence}\n")
                seen.add(accession)

        missing = sorted(set(accessions) - seen)
        if missing:
            raise BestHitSequenceEvidenceError(
                f"Missing reference sequence for {', '.join(missing)}.",
            )


def write_bed(path: Path, placements: list[Placement]) -> None:
    ordered = sorted(
        placements,
        key=lambda placement: (
            placement.reference_accession,
            placement.bed_start,
            placement.bed_end,
            placement.query_id,
        ),
    )
    with path.open("w", encoding="utf-8") as handle:
        for placement in ordered:
            handle.write(
                f"{placement.reference_accession}\t{placement.bed_start}\t"
                f"{placement.bed_end}\t{placement.query_id}\t0\t{placement.strand}\n",
            )


def emit_best_hit_sequence_evidence(args: argparse.Namespace) -> None:
    placements = load_placements(args.query_big_table)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not placements:
        for filename in OUTPUT_FILENAMES:
            (args.output_dir / filename).write_text("", encoding="utf-8")
        return

    sequences = load_query_sequences(args.query_fasta, placements)
    write_query_fasta(args.output_dir / "query_sequences.fasta", sequences)
    write_reference_fasta(
        args.output_dir / "selected_references.fasta",
        placements,
        args,
    )
    write_bed(args.output_dir / "best_hit_placements.bed", placements)


def main(argv: list[str] | None = None) -> None:
    emit_best_hit_sequence_evidence(parse_args(argv))


if __name__ == "__main__":
    main()
