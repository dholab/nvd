#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
# ]
# ///

"""Ensure the NCBI taxonomy cache is prepared for distributed workers."""

import argparse
import sys
import tarfile
import urllib.error
from pathlib import Path

from loguru import logger
from py_nvd.paths import get_taxdump_dir
from py_nvd.taxonomy import (
    TaxonomyBuildError,
    TaxonomyUnavailableError,
    ensure_taxonomy_available,
    format_taxonomy_warning,
)

_TAXONOMY_PREFLIGHT_ERRORS = (
    OSError,
    tarfile.TarError,
    urllib.error.URLError,
    TaxonomyBuildError,
    TaxonomyUnavailableError,
)


def configure_logging(verbose: bool = False) -> None:  # noqa: FBT001, FBT002
    """Configure loguru to write to stderr only."""
    logger.remove()
    logger.add(
        sys.stderr,
        format="<level>{level: <8}</level> | {message}",
        level="DEBUG" if verbose else "INFO",
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Ensure the NCBI taxonomy cache is ready.",
    )
    parser.add_argument(
        "--taxonomy-dir",
        type=Path,
        default=None,
        help="Explicit taxonomy directory; overrides state directory fallback",
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="Fail if taxonomy cannot be prepared during preflight",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug-level logging",
    )
    return parser.parse_args()


def main() -> None:
    """Prepare taxonomy cache if needed."""
    args = parse_args()
    configure_logging(verbose=args.verbose)

    taxdump_dir = get_taxdump_dir(
        taxonomy_dir=args.taxonomy_dir,
    )
    sqlite_path = taxdump_dir / "taxonomy.sqlite"

    logger.info("Ensuring taxonomy database is available...")
    if sqlite_path.exists():
        logger.info(f"Taxonomy database found: {taxdump_dir}")
        print(taxdump_dir, end="")  # noqa: T201 - stdout is the Nextflow value channel
        return

    logger.info("Taxonomy database not found, downloading from NCBI...")
    logger.info("This may take a few minutes on first run.")
    try:
        ensure_taxonomy_available(
            taxonomy_dir=args.taxonomy_dir,
        )
        logger.info(f"Taxonomy database ready: {taxdump_dir}")
        print(taxdump_dir, end="")  # noqa: T201 - stdout is the Nextflow value channel
    except _TAXONOMY_PREFLIGHT_ERRORS as error:
        if args.sync:
            logger.error(f"Taxonomy preflight failed and --sync was specified: {error}")
            sys.exit(1)

        warning = format_taxonomy_warning(
            operation="Pre-downloading taxonomy database",
            context="Preparing taxonomy for distributed workers",
            error=error,
            taxdump_dir=taxdump_dir,
            will_download=True,
        )
        logger.warning(warning)
        logger.warning(
            "Workers will attempt to download taxonomy when needed. "
            "This may cause delays and potential version inconsistencies.",
        )
        print(taxdump_dir, end="")  # noqa: T201 - stdout is the Nextflow value channel


if __name__ == "__main__":
    main()
