#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "polars",
# ]
# ///

import argparse
from pathlib import Path

import polars as pl
import polars.selectors as cs


def parse_args() -> argparse.Namespace:
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
    args = parse_args()
    top_k_lf = select_top_hits(Path(args.input_file), args.blast_retention_count)
    top_k_lf.sink_csv(args.output_file, separator="\t")


if __name__ == "__main__":
    main()
