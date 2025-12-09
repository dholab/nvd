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
    parser = argparse.ArgumentParser(
        description="Concatenated two potentially large BLAST results tables, one from MEGABLAST and one from BLASTN",
    )

    parser.add_argument(
        "--megablast-hits",
        required=True,
        help="Path to the MEGABLAST results text file.",
    )
    parser.add_argument(
        "--blastn-hits",
        required=True,
        help="Path to the BLASTN results text file.",
    )
    parser.add_argument(
        "--output-file",
        required=True,
        help="Path to write the concatenated BLAST results.",
    )

    return parser.parse_args()


def concat_blast_tables(
    megablast_hits: str | Path, blastn_hits: str | Path
) -> pl.LazyFrame:
    mb_hits = (
        # open a lazy file scanner than uses tab as its separator and refrains from inferring
        # the schema. This is a concession to the fact that BLAST taxids are sometimes semicolon
        # delimited and therefore not always numbers.
        pl.scan_csv(
            megablast_hits,
            separator="\t",
            infer_schema_length=10000,
            schema_overrides={"staxids": pl.Utf8, "bitscore": pl.Float64},
        )
        # split on any instances of a semicolon in staxids, which will replace the semicolon-delimited strings
        # with lists/arrays of taxid strings. Then, explode each item in those arrays into their own rows
        # and the convert the staxids column back into integers now that the semicolons have been handled
        .with_columns(pl.col("staxids").str.split(by=";").alias("_tmp1"))
        .explode("_tmp1")
        .with_columns(pl.col("_tmp1").str.to_integer().alias("staxids"))
        # get rid of duplicate entries just in case
        .unique()
        # get rid of temporary columns, which by convention I prefix with "_"
        .drop(cs.starts_with("_"))
    )
    bn_hits = (
        # do all the same things but for the blastn hits
        pl.scan_csv(
            blastn_hits,
            separator="\t",
            infer_schema_length=10000,
            schema_overrides={"staxids": pl.Utf8, "bitscore": pl.Float64},
        )
        .with_columns(pl.col("staxids").str.split(by=";").alias("_tmp1"))
        .explode("_tmp1")
        .with_columns(pl.col("_tmp1").str.to_integer().alias("staxids"))
        .unique()
        .drop(cs.starts_with("_"))
    )

    # return the concatenated lazyframes, sorted usefully by sample, query seq id, and blast task type
    return pl.concat([mb_hits, bn_hits]).sort(["sample", "qseqid", "task"])


def main() -> None:
    args = parse_args()
    concat_lf = concat_blast_tables(args.megablast_hits, args.blastn_hits)
    concat_lf.sink_csv(args.output_file, separator="\t")


if __name__ == "__main__":
    main()
