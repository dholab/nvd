#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Concatenate TSVs that already share a header, keeping the header once."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Concatenate TSVs sharing a header, writing the header once.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the concatenated TSV.",
    )
    parser.add_argument(
        "tsv_files",
        nargs="*",
        help="Input TSV paths to concatenate, in order.",
    )
    return parser.parse_args(argv)


def has_data(path: Path) -> bool:
    """Check if file exists and has more than just a header."""
    return path.exists() and path.stat().st_size > 0


def concat_tsv_with_header(tsv_files: list[Path], output: Path) -> None:
    """Write the first non-empty file's header once, then all data rows."""
    inputs = [path for path in tsv_files if has_data(path)]

    if not inputs:
        output.write_text("")
        return

    with output.open("w", encoding="utf-8") as out:
        with inputs[0].open(encoding="utf-8") as first:
            header = first.readline()
            out.write(header)
            for line in first:
                out.write(line)
        for tsv in inputs[1:]:
            with tsv.open(encoding="utf-8") as f:
                f.readline()  # skip header
                for line in f:
                    out.write(line)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    concat_tsv_with_header(
        [Path(path) for path in args.tsv_files],
        Path(args.output),
    )


if __name__ == "__main__":
    main()
