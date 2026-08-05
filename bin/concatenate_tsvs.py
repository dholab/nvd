#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#   "polars>=1.27.1",
# ]
# ///
"""Validate and concatenate typed TSV result artifacts with bounded memory."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence


NUMERIC_DTYPES = {
    "qlen": pl.Int64,
    "length": pl.Int64,
    "pident": pl.Float64,
    "evalue": pl.Float64,
    "bitscore": pl.Float64,
    "staxids": pl.Int64,
    "adjusted_taxid": pl.Int64,
    "support_record_count": pl.Int64,
    "mapped_reads": pl.Int64,
    "total_reads": pl.Int64,
    "crumbs_score": pl.Float64,
    "assigned_taxid": pl.Int64,
    "best_hit_qcov": pl.Float64,
    "best_hit_pident": pl.Float64,
    "best_hit_evalue": pl.Float64,
    "best_hit_bitscore": pl.Float64,
    "retained_reference_count": pl.Int64,
    "assignment_reference_count": pl.Int64,
    "assignment_taxid_count": pl.Int64,
    "best_hit_alignment_length": pl.Int64,
    "taxon_crumbs": pl.Float64,
    "relative_crumbs_percent": pl.Float64,
    "supporting_query_count": pl.Int64,
    "taxid": pl.Int64,
    "total_query_span": pl.Int64,
    "total_crumbs_score": pl.Float64,
    "strong_query_count": pl.Int64,
    "moderate_query_count": pl.Int64,
    "weak_query_count": pl.Int64,
    "review_query_count": pl.Int64,
    "redacted_query_count": pl.Int64,
    "supporting_genome_like_contig_count": pl.Int64,
    "supporting_long_contig_count": pl.Int64,
    "supporting_short_contig_count": pl.Int64,
    "supporting_merged_pair_count": pl.Int64,
    "supporting_single_read_count": pl.Int64,
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate and concatenate typed TSVs in filename order.",
    )
    parser.add_argument("--input", type=Path, nargs="*", default=[])
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def scan_tsvs(paths: Sequence[Path]) -> pl.LazyFrame:
    try:
        columns = (
            pl.scan_csv(
                paths[0],
                separator="\t",
                infer_schema=False,
            )
            .collect_schema()
            .names()
        )
        schema_overrides = {
            column: NUMERIC_DTYPES[column]
            for column in columns
            if column in NUMERIC_DTYPES
        }
        frame = pl.scan_csv(
            paths,
            separator="\t",
            infer_schema=False,
            schema_overrides=schema_overrides,
        )
        frame.collect_schema()
    except pl.exceptions.PolarsError as error:
        inputs = ", ".join(str(path) for path in paths)
        msg = f"could not read matching input TSV schemas for {inputs}: {error}"
        raise ValueError(msg) from error
    return frame


def concatenate_tsvs(
    paths: Sequence[Path],
    output: Path,
) -> None:
    if not paths:
        output.write_text("", encoding="utf-8")
        return

    ordered_paths = sorted(
        paths,
        key=lambda candidate: (candidate.name, str(candidate)),
    )
    missing_inputs = [path for path in ordered_paths if not path.is_file()]
    if missing_inputs:
        msg = f"input TSV does not exist: {missing_inputs[0]}"
        raise FileNotFoundError(msg)

    frame = scan_tsvs(ordered_paths)
    try:
        frame.sink_csv(
            output,
            separator="\t",
            maintain_order=True,
        )
    except pl.exceptions.PolarsError as error:
        inputs = ", ".join(str(path) for path in ordered_paths)
        msg = f"could not concatenate input TSVs {inputs}: {error}"
        raise ValueError(msg) from error


def input_paths(explicit_inputs: Sequence[Path], input_dir: Path | None) -> list[Path]:
    if input_dir is not None and not input_dir.is_dir():
        msg = f"{input_dir}: input directory does not exist"
        raise FileNotFoundError(msg)
    directory_inputs = sorted(input_dir.glob("*.tsv")) if input_dir else []
    return [*explicit_inputs, *directory_inputs]


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    paths = input_paths(args.input, args.input_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    concatenate_tsvs(paths, args.output)


if __name__ == "__main__":
    main()
