#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
# ]
# ///

"""Emit the stateless NVD run context for a samplesheet."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from loguru import logger
from py_nvd.run_context import compute_sample_set_id


def configure_logging(verbose: bool) -> None:  # noqa: FBT001
    logger.remove()
    logger.add(sys.stderr, level="INFO" if verbose else "WARNING", format="{message}")


def read_sample_ids(samplesheet: Path) -> list[str]:
    """Read non-comment sample IDs from an NVD samplesheet."""
    sample_ids: list[str] = []
    with samplesheet.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "sample_id" not in reader.fieldnames:
            msg = "Samplesheet must contain a sample_id column"
            raise ValueError(msg)

        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id or sample_id.startswith("#"):
                continue
            sample_ids.append(sample_id)

    if not sample_ids:
        msg = "Samplesheet contains no sample IDs"
        raise ValueError(msg)
    return sorted(set(sample_ids))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute stateless NVD run context")
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    try:
        sample_ids = read_sample_ids(args.samplesheet)
    except (OSError, ValueError) as exc:
        logger.error(str(exc))
        raise SystemExit(1) from exc

    sample_set_id = compute_sample_set_id(sample_ids)
    logger.info(f"Computed sample_set_id {sample_set_id} for {len(sample_ids)} samples")
    print(sample_set_id, end="")  # noqa: T201 - stdout is the Nextflow value channel


if __name__ == "__main__":
    main()
