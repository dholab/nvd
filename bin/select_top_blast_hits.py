#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "polars",
# ]
# ///

"""
Select top-scoring BLAST hits per query sequence.

This module filters BLAST results to retain only the top K hits (by bitscore)
for each query sequence. It handles BLAST-specific quirks like semicolon-delimited
taxids (when a single hit matches multiple taxa) by exploding them into separate
rows before ranking.

Key features:
- Memory-efficient lazy evaluation with Polars
- Handles semicolon-delimited staxids correctly
- Configurable retention count (default: top 5 hits per query)
- Grouping by task/sample/qseqid for multi-sample BLAST outputs
- Streaming output to avoid loading full results into memory

Usage:
    python select_top_blast_hits.py -i blast_results.txt -o top_hits.txt \\
        --blast-retention-count 5

Algorithm:
    1. Parse BLAST results as TSV with lazy evaluation
    2. Split semicolon-delimited taxids into separate rows
    3. Rank hits by descending bitscore within each query group
    4. Filter to top K ranks
    5. Stream directly to output file
"""

import argparse
from pathlib import Path

import polars as pl
import polars.selectors as cs


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for BLAST hit selection.

    Returns:
        Namespace with parsed arguments including I/O paths and retention count

    Side effects:
        - Reads sys.argv
        - May exit program if arguments invalid or --help requested
    """
    parser = argparse.ArgumentParser(description="Select top BLAST hits per query sequence.")

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
        help="Path to write the filtered BLAST results (top hits).",
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
    Correctly handles BLAST's semicolon-delimited staxids by exploding them into
    separate rows before ranking. Uses ordinal ranking on bitscores within each
    task/sample/qseqid group.

    Args:
        blast_txt: Path to BLAST results file (TSV format)
        top_k: Number of top hits to retain per query (default: 5)

    Returns:
        LazyFrame with top K hits per query, sorted by task/sample/qseqid.
        Temporary columns (_tmp1, _rank) are removed from output.

    Side effects: None (pure function, returns lazy computation graph)

    Minimum expected input columns:
        - staxids: Taxids (may be semicolon-delimited)
        - bitscore: Alignment score (used for ranking)
        - task: BLAST task identifier (e.g., "megablast", "blastn")
        - sample: Sample identifier
        - qseqid: Query sequence ID

    Algorithm:
        1. Cast staxids to strings (handles mixed types)
        2. Split semicolon-delimited taxids into lists
        3. Explode lists into separate rows (one row per taxid)
        4. Convert taxids back to integers
        5. Remove duplicates
        6. Rank by descending bitscore within each group
        7. Filter to top K ranks
        8. Clean up temporary columns
        9. Sort for readability

    Note: Uses lazy evaluation - no data is loaded until .collect() or .sink_*() called.
    """
    return (
        # open a lazy file scanner than uses tab as its separator and refrains from inferring
        # the schema. This is a concession to the fact that BLAST taxids are sometimes semicolon
        # delimited and therefore not always numbers.
        pl.scan_csv(
            blast_txt,
            separator="\t",
            infer_schema_length=10000,
        )
        # cast staxids into a column of strings and cast bitscores as floats. Bitscores should always
        # be floats, but they often get truncated because CSV is an untyped data encoding.
        .with_columns(pl.col("staxids").cast(pl.Utf8), pl.col("bitscore").cast(pl.Float64))
        # split on any instances of a semicolon in staxids, which will replace the semicolon-delimited strings
        # with lists/arrays of taxid strings. Then, explode each item in those arrays into their own rows
        # and the convert the staxids column back into integers now that the semicolons have been handled
        .with_columns(pl.col("staxids").str.split(by=";").alias("_tmp1"))
        .explode("_tmp1")
        .with_columns(pl.col("_tmp1").str.to_integer().alias("staxids"))
        # get rid of duplicate entries just in case
        .unique()
        # use ranking by descending-order bitscores to skim <= 5 hits off the top of
        # each task-sample-qseqid grouping
        .with_columns(
            pl.col("bitscore")
            .rank(method="ordinal", descending=True)  # highest score -> rank 1
            .over(["task", "sample", "qseqid"])
            .alias("_rank"),
        )
        .filter(pl.col("_rank") <= top_k)
        # get rid of temporary columns, which by convention I prefix with "_"
        .drop(cs.starts_with("_"))
        # sort by the aforementioned grouping for readability; these files are written
        # in plain human-readable text
        .sort(["task", "sample", "qseqid"])
    )


def main() -> None:
    """
    Main entry point for selecting top BLAST hits.

    Orchestrates the pipeline:
    1. Parse command-line arguments
    2. Select top K hits per query using lazy evaluation
    3. Stream results directly to output file

    Side effects:
        - Reads input BLAST file
        - Writes output file
        - May exit on invalid arguments or I/O errors

    Note: Uses streaming I/O to handle arbitrarily large BLAST result files
          without loading them entirely into memory.
    """
    args = parse_args()
    top_k_lf = select_top_hits(args.input_file, args.blast_retention_count)
    top_k_lf.sink_csv(args.output_file, separator="\t")


if __name__ == "__main__":
    main()
