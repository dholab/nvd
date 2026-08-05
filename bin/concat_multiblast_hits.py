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
        description="Concatenate one or more potentially large BLAST results tables.",
    )

    parser.add_argument(
        "--blast-hits",
        action="append",
        required=True,
        help="Path to a filtered BLAST results text file. Repeat for multiple files.",
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


def first_header(blast_hits: list[str | Path]) -> str | None:
    """Return the first header line from existing non-empty input files."""
    for filepath in blast_hits:
        path = Path(filepath)
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open(encoding="utf-8") as handle:
            header = handle.readline()
        if header:
            return header.rstrip("\n")
    return None


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
    blast_hits: list[str | Path],
) -> pl.LazyFrame | None:
    """Concatenate BLAST tables, handling empty files gracefully."""
    frames = [process_blast_file(path) for path in blast_hits if has_data(path)]
    if not frames:
        return None
    if len(frames) == 1:
        return frames[0].sort(["sample", "qseqid", "task"])
    return pl.concat(frames).sort(["sample", "qseqid", "task"])


def main() -> None:
    args = parse_args()
    concat_lf = concat_blast_tables(args.blast_hits)

    if concat_lf is not None:
        concat_lf.sink_csv(args.output_file, separator="\t")
    else:
        header = first_header(args.blast_hits) or (
            "task\tsample\tqseqid\tqlen\tsseqid\tstitle\tlength\tpident\tevalue\tbitscore\tsscinames\tstaxids\trank"
        )
        Path(args.output_file).write_text(f"{header}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
