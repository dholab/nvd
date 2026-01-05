#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
# ]
# ///

"""
Mark a pipeline run as completed or failed.

This script is called at the end of a Nextflow workflow to update the
run status in the state database. It should be gated on completion of
all sample processing (and uploads, if LabKey is enabled).

Usage:
    complete_run.py --run-id run_001 --state-dir /path/to/state
    complete_run.py --run-id run_001 --state-dir /path/to/state --status failed

Exit codes:
    0: Success (run marked as completed/failed, or run not found - non-fatal)
    1: Error (invalid arguments, database error, etc.)
"""

import argparse
import sys
from pathlib import Path

from loguru import logger
from py_nvd.state import complete_run


def configure_logging(verbosity: int) -> None:
    """Configure loguru logging based on verbosity level."""
    logger.remove()

    if verbosity == 0:
        level = "WARNING"
    elif verbosity == 1:
        level = "INFO"
    else:
        level = "DEBUG"

    logger.add(
        sys.stderr,
        level=level,
        format="<level>{level: <8}</level> | {message}",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mark a pipeline run as completed or failed",
    )
    parser.add_argument(
        "--run-id",
        type=str,
        required=True,
        help="Run identifier (workflow.runName)",
    )
    parser.add_argument(
        "--state-dir",
        type=Path,
        required=True,
        help="Path to state directory containing state.sqlite",
    )
    parser.add_argument(
        "--status",
        type=str,
        choices=["completed", "failed"],
        default="completed",
        help="Status to set (default: completed)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (-v for INFO, -vv for DEBUG)",
    )

    args = parser.parse_args()

    configure_logging(args.verbose)

    logger.info(f"Marking run {args.run_id} as {args.status}")

    try:
        success = complete_run(
            run_id=args.run_id,
            status=args.status,
            state_dir=args.state_dir,
        )

        if success:
            logger.info(f"Run {args.run_id} marked as {args.status}")
            print(f"Run {args.run_id} marked as {args.status}")
        else:
            # Run not found - this is non-fatal
            # Could happen if CHECK_RUN_STATE failed to register the run
            logger.warning(f"Run {args.run_id} not found in state database")
            print(f"Warning: Run {args.run_id} not found in state database")

    except Exception as e:
        logger.error(f"Failed to complete run: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
