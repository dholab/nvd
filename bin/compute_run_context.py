#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
# ]
# ///

"""Emit the state-free NVD run context for a samplesheet."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from loguru import logger
from py_nvd.read_inputs import ResolutionError, read_rows, resolve_rows
from py_nvd.run_context import compute_sample_set_id


def configure_logging(verbose: bool) -> None:  # noqa: FBT001
    logger.remove()
    logger.add(sys.stderr, level="INFO" if verbose else "WARNING", format="{message}")


def read_sample_ids(samplesheet: Path) -> list[str]:
    """Read sample IDs from the shared NVD samplesheet resolver."""
    resolution = resolve_rows(read_rows(samplesheet))
    return [sample.sample_id for sample in resolution.samples]


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute state-free NVD run context")
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    configure_logging(args.verbose)

    try:
        sample_ids = read_sample_ids(args.samplesheet)
    except (OSError, ResolutionError, ValueError) as exc:
        logger.error(str(exc))
        raise SystemExit(1) from exc

    sample_set_id = compute_sample_set_id(sample_ids)
    logger.info(f"Computed sample_set_id {sample_set_id} for {len(sample_ids)} samples")
    print(sample_set_id, end="")  # noqa: T201 - stdout is the Nextflow value channel


if __name__ == "__main__":
    main()
