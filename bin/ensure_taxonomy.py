#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
# ]
# ///

"""Ensure the NCBI taxonomy cache is prepared for distributed workers."""

import argparse
import os
import sys
import tarfile
import urllib.error
from pathlib import Path

from loguru import logger
from py_nvd.paths import get_taxdump_dir
from py_nvd.taxonomy import (
    DEFAULT_TAXONOMY_MAX_AGE_DAYS,
    TaxonomyAction,
    TaxonomyBuildError,
    TaxonomyMode,
    TaxonomyOfflineError,
    TaxonomyPolicy,
    TaxonomyRefresh,
    TaxonomyUnavailableError,
    ensure_taxonomy_available,
    execute_taxonomy_plan,
    format_taxonomy_warning,
    inspect_taxonomy,
    plan_pipeline_taxonomy,
)

_TAXONOMY_PREFLIGHT_ERRORS = (
    OSError,
    tarfile.TarError,
    urllib.error.URLError,
    TaxonomyBuildError,
    TaxonomyOfflineError,
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
        help="Explicit taxonomy directory; overrides NVD_TAXONOMY_DB and config-dir fallback",
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="Fail if taxonomy cannot be prepared during preflight",
    )
    parser.add_argument(
        "--taxonomy-refresh",
        choices=[refresh.value for refresh in TaxonomyRefresh],
        default=TaxonomyRefresh.MISSING.value,
        help="Admin refresh policy: missing, stale, or force",
    )
    parser.add_argument(
        "--taxonomy-mode",
        choices=[mode.value for mode in TaxonomyMode],
        default=None,
        help="Pipeline taxonomy mode: read_only or missing",
    )
    parser.add_argument(
        "--taxonomy-max-age-days",
        type=int,
        default=DEFAULT_TAXONOMY_MAX_AGE_DAYS,
        help="Freshness threshold for --taxonomy-refresh stale",
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
    env_offline = os.environ.get("NVD_TAXONOMY_OFFLINE", "").lower() in ("1", "true")
    requires_prepared_taxonomy = (
        args.taxonomy_mode == TaxonomyMode.READ_ONLY.value
        or args.taxonomy_mode is None
        and env_offline
    )

    taxdump_dir = get_taxdump_dir(
        taxonomy_dir=args.taxonomy_dir,
    )

    logger.info("Ensuring taxonomy database is available...")
    try:
        if args.taxonomy_mode is not None:
            policy = TaxonomyPolicy.from_env_and_params(
                mode=args.taxonomy_mode,
                max_age_days=args.taxonomy_max_age_days,
                offline=False,
            )
            status = inspect_taxonomy(taxonomy_dir=args.taxonomy_dir)
            plan = plan_pipeline_taxonomy(status, policy)
            if plan.action is TaxonomyAction.FAIL:
                raise TaxonomyOfflineError(plan.reason)
            execute_taxonomy_plan(plan)
        elif env_offline:
            policy = TaxonomyPolicy.from_env_and_params(
                mode=None,
                max_age_days=args.taxonomy_max_age_days,
                offline=True,
            )
            status = inspect_taxonomy(taxonomy_dir=args.taxonomy_dir)
            plan = plan_pipeline_taxonomy(status, policy)
            if plan.action is TaxonomyAction.FAIL:
                raise TaxonomyOfflineError(plan.reason)
            execute_taxonomy_plan(plan)
        else:
            logger.info(f"Taxonomy refresh policy: {args.taxonomy_refresh}")
            ensure_taxonomy_available(
                taxonomy_dir=args.taxonomy_dir,
                refresh=args.taxonomy_refresh,
                max_age_days=args.taxonomy_max_age_days,
            )
        logger.info(f"Taxonomy database ready: {taxdump_dir}")
        print(taxdump_dir, end="")  # noqa: T201 - stdout is the Nextflow value channel
    except _TAXONOMY_PREFLIGHT_ERRORS as error:
        if args.sync or requires_prepared_taxonomy:
            logger.error(f"Taxonomy preflight failed: {error}")
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
