#!/usr/bin/env python3
"""Normalize mapback-unmapped read batches into BLAST query FASTA and lookup rows."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import shutil
import sqlite3
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, TextIO

WRAP_WIDTH = 80
SAFE_TOKEN = re.compile(r"^[A-Za-z0-9_.-]+$")
READ_QUERY_PREFIXES = {
    "overlap_merged_pair": "nvdMergeReadQuery",
    "single_read": "nvdReadQuery",
}
SUPPORT_COUNT_POLICIES = frozenset({"abundance", "one"})
CASAVA_MATE_COMMENT = re.compile(r"(?P<mate>[12]):[YN]:\d+:\S*")
CASAVA_CLUSTER_ID = re.compile(
    r"[^:\s]+:\d+:[^:\s]+:\d+:\d+:\d+:\d+(?:/(?P<mate>[12]))?",
)


class ReadQueryNormalizationError(Exception):
    """Raised when read-query normalization cannot complete safely."""


@dataclass(frozen=True)
class VsearchUnique:
    """One exact unique sequence emitted by vsearch."""

    source_id: str
    abundance: int
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode()).hexdigest()


@dataclass(frozen=True)
class NormalizedReadQuery:
    """A read-derived BLAST query with NVD identity and lookup metadata."""

    qseqid: str
    sample_id: str
    query_class: str
    producer: str
    source_id: str
    support_record_count: int
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode()).hexdigest()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize mapback-unmapped read batches into BLAST query FASTA and lookup rows.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument(
        "--query-class",
        required=True,
        choices=sorted(READ_QUERY_PREFIXES),
    )
    parser.add_argument("--producer", required=True)
    parser.add_argument(
        "--fastq",
        action="append",
        type=Path,
        required=True,
        help="Compressed or plain FASTQ input. Repeat to normalize multiple files as one batch.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--support-count-policy",
        required=True,
        choices=sorted(SUPPORT_COUNT_POLICIES),
        help="Use vsearch abundance or one support unit per unique query.",
    )
    parser.add_argument("--vsearch-bin", default="vsearch")
    return parser.parse_args(argv)


def validate_safe_token(name: str, value: str) -> str:
    if not value or not SAFE_TOKEN.fullmatch(value):
        message = f"{name} must contain only letters, numbers, dots, underscores, or hyphens: {value!r}"
        raise ReadQueryNormalizationError(message)
    return value


def open_fastq_bytes(path: Path) -> BinaryIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def run_vsearch_fastx_uniques(
    fastqs: list[Path],
    *,
    output_fasta: Path,
    vsearch_bin: str,
) -> None:
    missing = [path for path in fastqs if not path.is_file()]
    if missing:
        message = "FASTQ inputs contain missing files:\n" + "\n".join(
            f"  {path}" for path in missing
        )
        raise ReadQueryNormalizationError(message)

    command = [
        vsearch_bin,
        "--fastx_uniques",
        "-",
        "--fastaout",
        str(output_fasta),
        "--sizeout",
        "--fasta_width",
        "0",
        "--minseqlength",
        "1",
        "--notrunclabels",
        "--quiet",
    ]
    process = subprocess.Popen(command, stdin=subprocess.PIPE)
    assert process.stdin is not None
    try:
        for fastq in fastqs:
            with open_fastq_bytes(fastq) as handle:
                shutil.copyfileobj(handle, process.stdin)
    finally:
        process.stdin.close()
    return_code = process.wait()
    if return_code != 0:
        message = f"vsearch failed with exit code {return_code}: {' '.join(command)}"
        raise ReadQueryNormalizationError(message)


def canonical_source_id(source_label: str) -> str:
    source_id, *comment = source_label.split(maxsplit=1)
    if not comment:
        return source_id

    casava_match = CASAVA_MATE_COMMENT.fullmatch(comment[0])
    cluster_match = CASAVA_CLUSTER_ID.fullmatch(source_id)
    if casava_match is None or cluster_match is None:
        return source_id

    casava_mate = casava_match.group("mate")
    existing_mate = cluster_match.group("mate")
    if existing_mate is None:
        return f"{source_id}/{casava_mate}"
    if existing_mate != casava_mate:
        message = f"conflicting mate annotations in FASTQ label: {source_label!r}"
        raise ReadQueryNormalizationError(message)
    return source_id


def parse_vsearch_header(header: str) -> tuple[str, int]:
    try:
        source_id, size_text = header.rsplit(";size=", maxsplit=1)
    except ValueError as error:
        message = f"vsearch FASTA header lacks ';size=' abundance: {header!r}"
        raise ReadQueryNormalizationError(message) from error
    try:
        abundance = int(size_text)
    except ValueError as error:
        message = f"vsearch FASTA header has non-integer abundance: {header!r}"
        raise ReadQueryNormalizationError(message) from error
    if abundance < 1:
        message = f"vsearch FASTA abundance must be positive: {header!r}"
        raise ReadQueryNormalizationError(message)
    return canonical_source_id(source_id), abundance


def parse_vsearch_fasta(handle: TextIO) -> list[VsearchUnique]:
    uniques: list[VsearchUnique] = []
    current_header: str | None = None
    current_sequence: list[str] = []

    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                uniques.append(unique_from_parts(current_header, current_sequence))
            current_header = line[1:].strip()
            if not current_header:
                message = f"FASTA header on line {line_number} is empty."
                raise ReadQueryNormalizationError(message)
            current_sequence = []
        elif current_header is None:
            message = f"FASTA sequence line {line_number} appears before any header."
            raise ReadQueryNormalizationError(message)
        else:
            current_sequence.append(line)

    if current_header is not None:
        uniques.append(unique_from_parts(current_header, current_sequence))
    return uniques


def unique_from_parts(header: str, sequence_parts: list[str]) -> VsearchUnique:
    source_id, abundance = parse_vsearch_header(header)
    sequence = "".join(sequence_parts)
    if not sequence:
        message = f"vsearch FASTA record {source_id!r} has no sequence."
        raise ReadQueryNormalizationError(message)
    return VsearchUnique(source_id=source_id, abundance=abundance, sequence=sequence)


def normalize_unique_reads(
    uniques: list[VsearchUnique],
    *,
    sample_id: str,
    query_class: str,
    producer: str,
    support_count_policy: str,
) -> list[NormalizedReadQuery]:
    validate_safe_token("sample_id", sample_id)
    validate_safe_token("producer", producer)
    if query_class not in READ_QUERY_PREFIXES:
        message = f"unsupported query_class: {query_class!r}"
        raise ReadQueryNormalizationError(message)
    if support_count_policy not in SUPPORT_COUNT_POLICIES:
        message = f"unsupported support_count_policy: {support_count_policy!r}"
        raise ReadQueryNormalizationError(message)

    prefix = READ_QUERY_PREFIXES[query_class]
    normalized: list[NormalizedReadQuery] = []
    seen_source_ids: set[str] = set()
    for index, unique in enumerate(uniques, start=1):
        if unique.source_id in seen_source_ids:
            message = (
                f"duplicate vsearch representative source_id: {unique.source_id!r}"
            )
            raise ReadQueryNormalizationError(message)
        seen_source_ids.add(unique.source_id)
        support_record_count = (
            unique.abundance if support_count_policy == "abundance" else 1
        )
        normalized.append(
            NormalizedReadQuery(
                qseqid=f"{prefix}_{sample_id}_{index:06d}",
                sample_id=sample_id,
                query_class=query_class,
                producer=producer,
                source_id=unique.source_id,
                support_record_count=support_record_count,
                sequence=unique.sequence,
            ),
        )
    return normalized


def wrapped(sequence: str) -> str:
    return "\n".join(
        sequence[start : start + WRAP_WIDTH]
        for start in range(0, len(sequence), WRAP_WIDTH)
    )


def write_fasta(path: Path, queries: list[NormalizedReadQuery]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for query in queries:
            handle.write(
                f">{query.qseqid} "
                f"query_class={query.query_class} "
                f"producer={query.producer} "
                f"source_id={query.source_id}\n",
            )
            handle.write(f"{wrapped(query.sequence)}\n")


def write_query_lookup(path: Path, queries: list[NormalizedReadQuery]) -> None:
    if path.exists():
        path.unlink()
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
        connection.executemany(
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
            [
                (
                    query.qseqid,
                    query.sample_id,
                    query.query_class,
                    query.producer,
                    query.source_id,
                    query.support_record_count,
                    query.length,
                    query.sha256,
                )
                for query in queries
            ],
        )
        connection.execute(
            "create index query_sequences_sample_id on query_sequences(sample_id)",
        )


def output_paths(
    output_dir: Path,
    *,
    sample_id: str,
    query_class: str,
) -> tuple[Path, Path]:
    return (
        output_dir / f"{sample_id}.{query_class}.blast_queries.fasta",
        output_dir / f"{sample_id}.{query_class}.query_sequences.sqlite",
    )


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    validate_safe_token("sample_id", args.sample_id)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    vsearch_fasta = (
        args.output_dir / f"{args.sample_id}.{args.query_class}.vsearch.uniques.fasta"
    )
    run_vsearch_fastx_uniques(
        args.fastq,
        output_fasta=vsearch_fasta,
        vsearch_bin=args.vsearch_bin,
    )
    with vsearch_fasta.open(encoding="utf-8") as handle:
        uniques = parse_vsearch_fasta(handle)
    queries = normalize_unique_reads(
        uniques,
        sample_id=args.sample_id,
        query_class=args.query_class,
        producer=args.producer,
        support_count_policy=args.support_count_policy,
    )
    if not queries:
        return
    fasta_path, lookup_path = output_paths(
        args.output_dir,
        sample_id=args.sample_id,
        query_class=args.query_class,
    )
    write_fasta(fasta_path, queries)
    write_query_lookup(lookup_path, queries)


if __name__ == "__main__":
    main()
