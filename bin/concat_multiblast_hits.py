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
        description="Concatenate two potentially large BLAST results tables, one from MEGABLAST and one from BLASTN",
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


def has_data(filepath: str | Path) -> bool:
    """Check if file exists and has more than just a header."""
    path = Path(filepath)
    if not path.exists() or path.stat().st_size == 0:
        return False
    with open(path) as f:
        return sum(1 for _ in f) > 1


def process_blast_file(filepath: str | Path) -> pl.LazyFrame:
    """Load and process a BLAST results file."""
    return (
        pl.scan_csv(
            filepath,
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


def concat_blast_tables(
    megablast_hits: str | Path,
    blastn_hits: str | Path,
) -> pl.LazyFrame | None:
    """Concatenate BLAST tables, handling empty files gracefully."""
    mb_has_data = has_data(megablast_hits)
    bn_has_data = has_data(blastn_hits)

    if mb_has_data and bn_has_data:
        # Both have data - concat them
        return pl.concat(
            [
                process_blast_file(megablast_hits),
                process_blast_file(blastn_hits),
            ],
        ).sort(["sample", "qseqid", "task"])
    if mb_has_data:
        # Only megablast has data
        return process_blast_file(megablast_hits).sort(["sample", "qseqid", "task"])
    if bn_has_data:
        # Only blastn has data
        return process_blast_file(blastn_hits).sort(["sample", "qseqid", "task"])
    # Neither has data
    return None


def main() -> None:
    args = parse_args()
    concat_lf = concat_blast_tables(args.megablast_hits, args.blastn_hits)

    if concat_lf is not None:
        concat_lf.sink_csv(args.output_file, separator="\t")
    else:
        # Create empty file with header
        Path(args.output_file).write_text(
            "task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids\trank\n",
        )


if __name__ == "__main__":
    main()
