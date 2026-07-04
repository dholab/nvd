#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///

"""Build a conservative union FASTA from assembler contigs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

DNA_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)
VALID_SEQUENCE = frozenset("ACGTRYSWKMBDHVN")
SAFE_ID_RE = re.compile(r"[^A-Za-z0-9_.-]+")
PROVENANCE_COLUMNS = (
    "sample_id",
    "assembler",
    "original_header",
    "contig_id",
    "length",
    "sha256",
    "representative_id",
    "disposition",
    "reason",
    "container_id",
    "orientation",
    "start_1based",
    "end_1based",
)
SCHEMA = """
create table contig (
  contig_id text primary key,
  sample_id text not null,
  assembler text not null,
  original_header text not null,
  source_fasta text not null,
  length integer not null,
  sha256 text not null,
  canonical_sha256 text not null,
  emitted integer not null,
  representative_id text
);

create table duplicate (
  duplicate_id text not null,
  representative_id text not null,
  orientation text not null,
  proof text not null,
  primary key (duplicate_id, representative_id)
);

create table containment (
  contained_id text not null,
  container_id text not null,
  orientation text not null,
  start_1based integer not null,
  end_1based integer not null,
  proof text not null,
  primary key (contained_id, container_id, orientation, start_1based, end_1based)
);
"""


class Orientation(StrEnum):
    SAME = "same"
    REVERSE_COMPLEMENT = "reverse_complement"


class Proof(StrEnum):
    EXACT_DUPLICATE = "canonical_sha256_exact_duplicate"
    EXACT_CONTAINMENT = "python_verified_exact_containment"


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


@dataclass(frozen=True)
class SourceContig:
    contig_id: str
    representative_id: str
    sample_id: str
    assembler: str
    original_header: str
    source_fasta: str
    sequence: str
    sha256: str
    canonical_sha256: str

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class FinalizedUnionPaths:
    fasta: Path
    provenance_tsv: Path
    summary_json: Path


def safe_id(value: str) -> str:
    cleaned = SAFE_ID_RE.sub("_", value.strip())
    return cleaned.strip("_") or "sample"


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1].upper()


def canonical_hash(sequence: str) -> str:
    forward = sequence_hash(sequence)
    reverse = sequence_hash(reverse_complement(sequence))
    return min(forward, reverse)


@dataclass(frozen=True)
class InputFasta:
    index: int
    assembler: str
    path: Path


def duplicate_orientation(sequence: str, representative: str) -> Orientation:
    if sequence.upper() == representative.upper():
        return Orientation.SAME
    if reverse_complement(sequence) == representative.upper():
        return Orientation.REVERSE_COMPLEMENT
    msg = "duplicate group member is not identical to representative on either strand"
    raise ValueError(msg)


def containment_interval(
    query: str,
    target: str,
) -> tuple[Orientation, int, int] | None:
    query_upper = query.upper()
    target_upper = target.upper()
    index = target_upper.find(query_upper)
    if index >= 0:
        return (Orientation.SAME, index + 1, index + len(query_upper))

    rc_query = reverse_complement(query_upper)
    index = target_upper.find(rc_query)
    if index >= 0:
        return (Orientation.REVERSE_COMPLEMENT, index + 1, index + len(rc_query))
    return None


def validate_sequence(sequence: str, *, path: Path, header: str) -> str:
    normalized = sequence.upper()
    invalid = sorted(set(normalized) - VALID_SEQUENCE)
    if invalid:
        msg = (
            f"Unsupported FASTA symbols in {path} record {header!r}: {''.join(invalid)}"
        )
        raise ValueError(msg)
    return normalized


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    sequence_parts: list[str] = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if header is not None:
                    if not sequence_parts:
                        msg = f"FASTA record {header!r} in {path} has no sequence"
                        raise ValueError(msg)
                    sequence = validate_sequence(
                        "".join(sequence_parts),
                        path=path,
                        header=header,
                    )
                    records.append(FastaRecord(header, sequence))
                header = stripped[1:].strip()
                if not header:
                    msg = f"FASTA record in {path} has an empty header"
                    raise ValueError(msg)
                sequence_parts = []
            else:
                if header is None:
                    msg = f"FASTA sequence appears before first header in {path}"
                    raise ValueError(msg)
                sequence_parts.append(stripped)

    if header is not None:
        if not sequence_parts:
            msg = f"FASTA record {header!r} in {path} has no sequence"
            raise ValueError(msg)
        sequence = validate_sequence("".join(sequence_parts), path=path, header=header)
        records.append(FastaRecord(header, sequence))
    return records


def write_fasta(records: list[FastaRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            for start in range(0, len(record.sequence), 80):
                handle.write(record.sequence[start : start + 80] + "\n")


def read_inputs_manifest(path: Path) -> list[InputFasta]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"assembler", "source_fasta"}
        if not required <= set(reader.fieldnames or []):
            msg = "inputs manifest must contain assembler and source_fasta columns"
            raise ValueError(msg)
        inputs: list[InputFasta] = []
        for index, row in enumerate(reader):
            assembler = row["assembler"].strip()
            source_fasta = row["source_fasta"].strip()
            if not assembler or not source_fasta:
                msg = f"inputs manifest row {index + 2} has empty assembler or source_fasta"
                raise ValueError(msg)
            inputs.append(
                InputFasta(index=index, assembler=assembler, path=Path(source_fasta)),
            )
        return inputs


def initialize_database(path: Path) -> sqlite3.Connection:
    if path.exists():
        path.unlink()
    connection = sqlite3.connect(path)
    connection.executescript(SCHEMA)
    return connection


def connect_existing_database(path: Path) -> sqlite3.Connection:
    if not path.is_file():
        msg = f"SQLite sidecar does not exist: {path}"
        raise FileNotFoundError(msg)
    return sqlite3.connect(path)


def source_contig_sort_key(contig: SourceContig) -> tuple[str, str, str, str, str]:
    return (
        contig.canonical_sha256,
        contig.sha256,
        contig.original_header,
        contig.source_fasta,
        contig.assembler,
    )


def build_source_contigs(
    *,
    sample_id: str,
    inputs: list[InputFasta],
) -> list[SourceContig]:
    sample_prefix = safe_id(sample_id)
    provisional: list[tuple[str, str, str, str, str, str, str]] = []
    for input_fasta in inputs:
        for record_index, record in enumerate(read_fasta(input_fasta.path)):
            sha = sequence_hash(record.sequence)
            canonical = canonical_hash(record.sequence)
            provenance_hash = hashlib.sha256(
                "\t".join(
                    [
                        str(input_fasta.index),
                        str(record_index),
                        input_fasta.assembler,
                        record.header,
                        sha,
                    ],
                ).encode("utf-8"),
            ).hexdigest()
            provisional.append(
                (
                    canonical,
                    sha,
                    record.header,
                    str(input_fasta.path),
                    input_fasta.assembler,
                    record.sequence,
                    provenance_hash,
                ),
            )

    provisional.sort()
    contigs: list[SourceContig] = []
    for (
        canonical,
        sha,
        header,
        source_fasta,
        assembler,
        sequence,
        provenance_hash,
    ) in provisional:
        contig_id = f"{sample_prefix}__source_{provenance_hash}"
        representative_id = f"{sample_prefix}__union_{canonical}"
        contigs.append(
            SourceContig(
                contig_id=contig_id,
                representative_id=representative_id,
                sample_id=sample_id,
                assembler=assembler,
                original_header=header,
                source_fasta=source_fasta,
                sequence=sequence,
                sha256=sha,
                canonical_sha256=canonical,
            ),
        )
    return contigs


def write_contigs(connection: sqlite3.Connection, contigs: list[SourceContig]) -> None:
    connection.executemany(
        """
        insert into contig (
          contig_id, sample_id, assembler, original_header, source_fasta,
          length, sha256, canonical_sha256, emitted, representative_id
        ) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                contig.contig_id,
                contig.sample_id,
                contig.assembler,
                contig.original_header,
                contig.source_fasta,
                contig.length,
                contig.sha256,
                contig.canonical_sha256,
                0,
                contig.representative_id,
            )
            for contig in contigs
        ],
    )


def prepare_union(
    *,
    sample_id: str,
    inputs_manifest: Path,
    representatives_fasta: Path,
    sqlite_path: Path,
    meta_json: Path,
) -> None:
    inputs = read_inputs_manifest(inputs_manifest)
    contigs = build_source_contigs(sample_id=sample_id, inputs=inputs)
    connection = initialize_database(sqlite_path)
    write_contigs(connection, contigs)

    representatives: list[FastaRecord] = []
    for _canonical, group in group_by_canonical(contigs):
        representative = min(group, key=source_contig_sort_key)
        representatives.append(
            FastaRecord(representative.representative_id, representative.sequence),
        )
        connection.execute(
            "update contig set emitted = 1 where contig_id = ?",
            (representative.contig_id,),
        )
        for duplicate in group:
            if duplicate.contig_id == representative.contig_id:
                continue
            orientation = duplicate_orientation(
                duplicate.sequence,
                representative.sequence,
            )
            connection.execute(
                """
                insert into duplicate (
                  duplicate_id, representative_id, orientation, proof
                ) values (?, ?, ?, ?)
                """,
                (
                    duplicate.contig_id,
                    representative.representative_id,
                    orientation.value,
                    Proof.EXACT_DUPLICATE.value,
                ),
            )

    representatives.sort(key=lambda record: record.header)
    write_fasta(representatives, representatives_fasta)
    meta = {
        "sample_id": sample_id,
        "representative_count": len(representatives),
        "max_contig_length": max(
            (len(record.sequence) for record in representatives),
            default=0,
        ),
    }
    meta_json.write_text(
        json.dumps(meta, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    connection.commit()
    connection.close()


def group_by_canonical(
    contigs: list[SourceContig],
) -> list[tuple[str, list[SourceContig]]]:
    groups: dict[str, list[SourceContig]] = {}
    for contig in contigs:
        groups.setdefault(contig.canonical_sha256, []).append(contig)
    return sorted(groups.items())


def read_representatives(path: Path) -> dict[str, str]:
    representatives: dict[str, str] = {}
    for record in read_fasta(path):
        if record.header in representatives:
            msg = f"Duplicate representative FASTA header: {record.header}"
            raise ValueError(msg)
        representatives[record.header] = record.sequence
    return representatives


def parse_candidate_rows(path: Path) -> list[dict[str, str]]:
    fields = [
        "query",
        "target",
        "pident",
        "alnlen",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "qlen",
        "tlen",
        "qcov",
        "tcov",
        "qframe",
        "tframe",
    ]
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            if parts == fields:
                continue
            if len(parts) != len(fields):
                msg = f"MMseqs2 candidate row has {len(parts)} fields; expected {len(fields)}"
                raise ValueError(msg)
            rows.append(dict(zip(fields, parts, strict=True)))
    return rows


def canonical_containments(
    connection: sqlite3.Connection,
) -> dict[str, sqlite3.Row]:
    containments: dict[str, sqlite3.Row] = {}
    rows = connection.execute(
        """
        select
          containment.contained_id,
          containment.container_id,
          containment.orientation,
          containment.start_1based,
          containment.end_1based,
          containment.proof,
          max(contig.length) as container_length
        from containment
        join contig on contig.representative_id = containment.container_id
        group by
          containment.contained_id,
          containment.container_id,
          containment.orientation,
          containment.start_1based,
          containment.end_1based,
          containment.proof
        order by
          containment.contained_id,
          container_length,
          containment.container_id,
          containment.orientation,
          containment.start_1based,
          containment.end_1based
        """,
    )
    for row in rows:
        containments.setdefault(row["contained_id"], row)
    return containments


def write_provenance(
    connection: sqlite3.Connection,
    *,
    path: Path,
) -> None:
    duplicates = {
        row["duplicate_id"]: row
        for row in connection.execute(
            "select duplicate_id, orientation, proof from duplicate",
        )
    }
    containments = canonical_containments(connection)
    contigs = connection.execute(
        """
        select
          sample_id,
          assembler,
          original_header,
          contig_id,
          length,
          sha256,
          representative_id,
          emitted
        from contig
        order by assembler, original_header, contig_id
        """,
    )

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROVENANCE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for contig in contigs:
            duplicate = duplicates.get(contig["contig_id"])
            containment = containments.get(contig["representative_id"])
            if duplicate is not None:
                disposition = "duplicate"
                reason = duplicate["proof"]
                orientation = duplicate["orientation"]
                container_id = ""
                start_1based = ""
                end_1based = ""
            elif containment is not None:
                disposition = "contained"
                reason = containment["proof"]
                orientation = containment["orientation"]
                container_id = containment["container_id"]
                start_1based = containment["start_1based"]
                end_1based = containment["end_1based"]
            elif contig["emitted"]:
                disposition = "emitted"
                reason = "selected_union_representative"
                orientation = ""
                container_id = ""
                start_1based = ""
                end_1based = ""
            else:
                msg = f"Contig {contig['contig_id']} has no final union disposition"
                raise ValueError(msg)

            writer.writerow(
                {
                    "sample_id": contig["sample_id"],
                    "assembler": contig["assembler"],
                    "original_header": contig["original_header"],
                    "contig_id": contig["contig_id"],
                    "length": contig["length"],
                    "sha256": contig["sha256"],
                    "representative_id": contig["representative_id"],
                    "disposition": disposition,
                    "reason": reason,
                    "container_id": container_id,
                    "orientation": orientation,
                    "start_1based": start_1based,
                    "end_1based": end_1based,
                },
            )


def write_summary(
    connection: sqlite3.Connection,
    *,
    sample_id: str,
    path: Path,
) -> None:
    def count(query: str) -> int:
        return int(connection.execute(query).fetchone()[0])

    summary = {
        "sample_id": sample_id,
        "source_contig_count": count("select count(*) from contig"),
        "representative_count": count(
            "select count(distinct representative_id) from contig",
        ),
        "emitted_contig_count": count("select count(*) from contig where emitted = 1"),
        "exact_duplicate_count": count("select count(*) from duplicate"),
        "exact_containment_count": count(
            "select count(distinct contained_id) from containment",
        ),
    }
    path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def finalize_union(
    *,
    sample_id: str,
    representatives_fasta: Path,
    candidates_tsv: Path,
    sqlite_path: Path,
    outputs: FinalizedUnionPaths,
) -> None:
    representatives = read_representatives(representatives_fasta)
    connection = connect_existing_database(sqlite_path)
    connection.row_factory = sqlite3.Row
    contained: set[str] = set()

    for row in parse_candidate_rows(candidates_tsv):
        query_id = row["query"]
        target_id = row["target"]
        if query_id == target_id:
            continue
        query = representatives.get(query_id)
        target = representatives.get(target_id)
        if query is None or target is None:
            sys.stderr.write(
                "MMseqs2 candidate referenced a contig absent from representatives FASTA: "
                f"query={query_id} target={target_id}\n",
            )
            continue
        if len(query) == len(target):
            if containment_interval(query, target) is not None:
                sys.stderr.write(
                    "Equal-length exact duplicate candidate appeared after prepare duplicate collapse: "
                    f"query={query_id} target={target_id}\n",
                )
            continue
        if len(query) > len(target):
            continue

        interval = containment_interval(query, target)
        if interval is None:
            continue
        orientation, start_1based, end_1based = interval
        contained.add(query_id)
        connection.execute(
            """
            insert or ignore into containment (
              contained_id, container_id, orientation, start_1based, end_1based, proof
            ) values (?, ?, ?, ?, ?, ?)
            """,
            (
                query_id,
                target_id,
                orientation.value,
                start_1based,
                end_1based,
                Proof.EXACT_CONTAINMENT.value,
            ),
        )

    for representative_id in contained:
        connection.execute(
            "update contig set emitted = 0 where representative_id = ?",
            (representative_id,),
        )

    emitted_records = [
        FastaRecord(header, sequence)
        for header, sequence in sorted(representatives.items())
        if header not in contained
    ]
    write_fasta(emitted_records, outputs.fasta)
    write_provenance(connection, path=outputs.provenance_tsv)
    write_summary(connection, sample_id=sample_id, path=outputs.summary_json)
    connection.commit()
    connection.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--sample-id", required=True)
    prepare.add_argument("--inputs-manifest", type=Path, required=True)
    prepare.add_argument("--representatives-fasta", type=Path, required=True)
    prepare.add_argument("--sqlite", type=Path, required=True)
    prepare.add_argument("--meta-json", type=Path, required=True)

    finalize = subparsers.add_parser("finalize")
    finalize.add_argument("--sample-id", required=True)
    finalize.add_argument("--representatives-fasta", type=Path, required=True)
    finalize.add_argument("--candidates-tsv", type=Path, required=True)
    finalize.add_argument("--sqlite", type=Path, required=True)
    finalize.add_argument("--output-fasta", type=Path, required=True)
    finalize.add_argument("--provenance-tsv", type=Path, required=True)
    finalize.add_argument("--summary-json", type=Path, required=True)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            prepare_union(
                sample_id=args.sample_id,
                inputs_manifest=args.inputs_manifest,
                representatives_fasta=args.representatives_fasta,
                sqlite_path=args.sqlite,
                meta_json=args.meta_json,
            )
        elif args.command == "finalize":
            finalize_union(
                sample_id=args.sample_id,
                representatives_fasta=args.representatives_fasta,
                candidates_tsv=args.candidates_tsv,
                sqlite_path=args.sqlite,
                outputs=FinalizedUnionPaths(
                    fasta=args.output_fasta,
                    provenance_tsv=args.provenance_tsv,
                    summary_json=args.summary_json,
                ),
            )
        else:
            parser.error(f"unknown command: {args.command}")
    except (OSError, ValueError, sqlite3.Error) as exc:
        sys.stderr.write(f"{exc}\n")
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
