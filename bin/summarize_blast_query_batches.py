#!/usr/bin/env python3
"""Summarize per-sample BLAST query batches by query class."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import polars as pl

EXPECTED_QUERY_CLASSES = (
    "short_assembly_contig",
    "long_assembly_contig",
    "overlap_merged_pair",
    "single_read",
)
QUERY_SOURCES = {
    "short_assembly_contig": "contig",
    "long_assembly_contig": "contig",
    "overlap_merged_pair": "read_query",
    "single_read": "read_query",
}
OUTPUT_COLUMNS = [
    "sample_id",
    "platform",
    "query_class",
    "query_source",
    "n_query_sequences",
    "query_fasta_present",
    "query_lookup_present",
    "query_fasta",
    "query_lookup",
]
ACTUAL_SCHEMA = {
    "query_class": pl.String,
    "n_query_sequences": pl.Int64,
    "query_fasta_present": pl.Boolean,
    "query_lookup_present": pl.Boolean,
    "query_fasta": pl.String,
    "query_lookup": pl.String,
}


class QueryBatchSummaryError(ValueError):
    """Raised when BLAST query batch summaries cannot be produced safely."""


@dataclass(frozen=True)
class QueryBatch:
    query_class: str
    fasta: Path
    lookup: Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize per-sample BLAST query batches by query class.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--batch",
        action="append",
        nargs=3,
        metavar=("QUERY_CLASS", "QUERY_FASTA", "QUERY_LOOKUP"),
        default=[],
        help="One BLAST query batch. Repeat for each emitted query class.",
    )
    return parser.parse_args(argv)


def count_fasta_records(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def batches_from_args(raw_batches: list[list[str]]) -> list[QueryBatch]:
    return [
        QueryBatch(query_class=query_class, fasta=Path(fasta), lookup=Path(lookup))
        for query_class, fasta, lookup in raw_batches
    ]


def expected_query_classes() -> pl.LazyFrame:
    return pl.DataFrame(
        {
            "query_class": list(EXPECTED_QUERY_CLASSES),
            "query_source": [QUERY_SOURCES[query_class] for query_class in EXPECTED_QUERY_CLASSES],
        },
    ).lazy()


def validate_unique_query_classes(batches: list[QueryBatch]) -> None:
    if not batches:
        return
    duplicates = (
        pl.DataFrame({"query_class": [batch.query_class for batch in batches]})
        .lazy()
        .group_by("query_class")
        .len()
        .filter(pl.col("len") > 1)
        .select("query_class")
        .collect()
    )
    if duplicates.height > 0:
        query_class = duplicates.get_column("query_class").item(0)
        message = f"duplicate BLAST query batch for query_class {query_class!r}"
        raise QueryBatchSummaryError(message)


def actual_batches(batches: list[QueryBatch]) -> pl.LazyFrame:
    if not batches:
        return pl.DataFrame(schema=ACTUAL_SCHEMA).lazy()
    return pl.DataFrame(
        {
            "query_class": [batch.query_class for batch in batches],
            "n_query_sequences": [count_fasta_records(batch.fasta) for batch in batches],
            "query_fasta_present": [batch.fasta.is_file() for batch in batches],
            "query_lookup_present": [batch.lookup.is_file() for batch in batches],
            "query_fasta": [str(batch.fasta) for batch in batches],
            "query_lookup": [str(batch.lookup) for batch in batches],
        },
        schema=ACTUAL_SCHEMA,
    ).lazy()


def summarize_batches(
    *,
    sample_id: str,
    platform: str,
    batches: list[QueryBatch],
) -> pl.LazyFrame:
    validate_unique_query_classes(batches)
    return (
        expected_query_classes()
        .join(actual_batches(batches), on="query_class", how="left")
        .with_columns(
            pl.lit(sample_id).alias("sample_id"),
            pl.lit(platform).alias("platform"),
            pl.col("n_query_sequences").fill_null(0),
            pl.col("query_fasta_present").fill_null(False),
            pl.col("query_lookup_present").fill_null(False),
            pl.col("query_fasta").fill_null(""),
            pl.col("query_lookup").fill_null(""),
        )
        .select(OUTPUT_COLUMNS)
    )


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    summary = summarize_batches(
        sample_id=args.sample_id,
        platform=args.platform,
        batches=batches_from_args(args.batch),
    )
    summary.collect().write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
