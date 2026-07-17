#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "polars",
# ]
# ///

"""
Select top-scoring BLAST hits per query sequence.

This module retains the top K BLAST reference positions by bitscore for each
query sequence, including every reference tied at the Kth position. It handles
BLAST-specific quirks like semicolon-delimited taxids after reference ranking so
one multi-taxid reference cannot consume several retention positions.

Key features:
- Memory-efficient lazy evaluation with Polars
- Handles semicolon-delimited staxids correctly
- Configurable retention count (default: top 5 hits per query)
- Per-query ranking within each staged BLAST result file
- Explicit ten-field parsing with malformed-row failures

Usage:
    python select_top_blast_hits.py -i blast_results.txt -o top_hits.txt \\
        --blast-retention-count 5

Algorithm:
    1. Parse BLAST results as TSV with lazy evaluation
    2. Collapse repeated rows for the same query/reference to their best row
    3. Rank references by descending bitscore within each query group
    4. Retain the top K positions and every boundary tie
    5. Split semicolon-delimited taxids into separate rows
    6. Materialize the retained result and write it to the output file
"""

import argparse
import csv
import os
from pathlib import Path

import polars as pl

EXPECTED_BLAST_FIELDS = 10


def require_blast_row_width(path: str | Path) -> None:
    """Require every headerless BLAST row to contain exactly ten fields."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(
            csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE),
            start=1,
        ):
            if len(row) != EXPECTED_BLAST_FIELDS:
                message = (
                    f"BLAST row {line_number} has {len(row)} fields; "
                    f"expected {EXPECTED_BLAST_FIELDS}"
                )
                raise ValueError(message)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for BLAST hit selection.

    Returns:
        Namespace with parsed arguments including I/O paths and retention count

    Side effects:
        - Reads sys.argv
        - May exit program if arguments invalid or --help requested
    """
    parser = argparse.ArgumentParser(
        description="Select top BLAST hits per query sequence.",
    )

    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="Path to the input BLAST results text file.",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        required=True,
        help="Path to write the retained BLAST references.",
    )
    parser.add_argument(
        "--blast-retention-count",
        type=int,
        default=5,
        help="Number of top hits to retain per query (default: 5).",
    )

    return parser.parse_args()


def select_top_hits(blast_txt: str | Path, top_k: int = 5) -> pl.LazyFrame:
    """
    Select the top K highest-scoring BLAST hits per query sequence.

    Processes BLAST results using lazy evaluation to handle large files efficiently.
    Correctly handles BLAST's semicolon-delimited staxids by ranking reference
    sequences before exploding them into separate taxid rows. Uses minimum
    ranking so every bitscore tie crossing the retention boundary survives.

    Args:
        blast_txt: Path to BLAST results file (TSV format)
        top_k: Number of top hits to retain per query (default: 5)

    Returns:
        LazyFrame with retained references per query, sorted by query, score,
        and reference ID. Temporary columns are removed from output.

    Side effects: None (pure function, returns lazy computation graph)

    Minimum expected input columns:
        - staxids: Taxids (may be semicolon-delimited)
        - bitscore: Alignment score (used for ranking)
        - qseqid: Query sequence ID

    Algorithm:
        1. Collapse repeated query/reference rows to the best alignment row
        2. Rank references by descending bitscore within each query
        3. Retain the first K positions and every tie at the Kth position
        4. Split semicolon-delimited taxids into lists
        5. Explode lists into separate rows (one row per taxid)
        6. Convert taxids back to integers
        7. Remove duplicates and temporary columns
        8. Sort for readability

    Note: Uses lazy evaluation - no data is loaded until .collect() or .sink_*() called.
    """
    return (
        # Parse the exact headerless BLAST -outfmt 6 contract. Numeric or ragged
        # rows fail here rather than being silently discarded or truncated.
        pl.scan_csv(
            blast_txt,
            separator="\t",
            quote_char=None,
            has_header=False,
            schema={
                "qseqid": pl.String,
                "qlen": pl.Int64,
                "sseqid": pl.String,
                "stitle": pl.String,
                "length": pl.Int64,
                "pident": pl.Float64,
                "evalue": pl.Float64,
                "bitscore": pl.Float64,
                "sscinames": pl.String,
                "staxids": pl.String,
            },
        )
        # One reference can have multiple HSP rows. Represent it by the best
        # row before assigning retention positions.
        .sort(
            ["qseqid", "bitscore", "evalue", "pident", "length"],
            descending=[False, True, False, True, True],
        )
        .unique(subset=["qseqid", "sseqid"], keep="first", maintain_order=True)
        # Minimum ranking gives every tied reference the first ordinal position
        # occupied by its score, preserving ties across the Kth boundary.
        .with_columns(
            pl.col("bitscore")
            .rank(method="min", descending=True)
            .over(["qseqid"])
            .alias("_rank"),
        )
        .filter(pl.col("_rank") <= top_k)
        # Expand multi-taxid references only after reference retention.
        .with_columns(pl.col("staxids").cast(pl.Utf8).str.split(by=";").alias("_taxid"))
        .explode("_taxid")
        .with_columns(pl.col("_taxid").str.to_integer(strict=False).alias("staxids"))
        .unique(maintain_order=True)
        .drop("_rank", "_taxid")
        # sort by the aforementioned grouping for readability; these files are written
        # in plain human-readable text
        .sort(["qseqid", "bitscore", "sseqid"], descending=[False, True, False])
    )


def main() -> None:
    """
    Main entry point for selecting top BLAST hits.

    Orchestrates the pipeline:
    1. Parse command-line arguments
    2. Select top K hits per query using lazy evaluation
    3. Materialize the retention-reduced result and write it to the output file

    Side effects:
        - Reads input BLAST file
        - Writes output file
        - May exit on invalid arguments or I/O errors

    Note: The scan and transformations remain lazy until the retained result is
          collected for CSV output.
    """
    args = parse_args()

    # Handle empty input file
    if os.path.getsize(args.input_file) == 0:
        Path(args.output_file).touch()
        return

    require_blast_row_width(args.input_file)
    top_k_lf = select_top_hits(args.input_file, args.blast_retention_count)
    top_k_lf.collect().write_csv(args.output_file, separator="\t")


if __name__ == "__main__":
    main()
