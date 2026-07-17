#!/usr/bin/env python3
"""Resolve NVD samplesheet read declarations into a canonical JSONL manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from py_nvd.read_inputs import ResolutionError, read_rows, resolve_rows, write_jsonl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve NVD samplesheet read declarations into JSONL.",
    )
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--output-jsonl", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        resolution = resolve_rows(read_rows(args.samplesheet))
    except ResolutionError as error:
        sys.stderr.write(f"{error}\n")
        raise SystemExit(1) from error

    write_jsonl(resolution.samples, args.output_jsonl)
    for warning in resolution.warnings:
        sys.stderr.write(f"Warning: {warning}\n")


if __name__ == "__main__":
    main()
