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

Graceful Degradation:
    By default, if the state database is unavailable (permission errors,
    corruption, etc.), this script will warn and exit successfully.
    Use --sync to require the state database and fail if unavailable.

Usage:
    complete_run.py --run-id run_001 --state-dir /path/to/state
    complete_run.py --run-id run_001 --state-dir /path/to/state --status failed

    # Require state database (fail if unavailable)
    complete_run.py --run-id run_001 --state-dir /path/to/state --sync

Exit codes:
    0: Success (run marked as completed/failed, or state DB unavailable without --sync)
    1: Error (invalid arguments, or state DB unavailable with --sync)
"""

import argparse
import sys
from pathlib import Path

from loguru import logger
from py_nvd.db import (
    create_backup,
    format_state_warning,
    get_state_db_path,
    StateUnavailableError,
)
from py_nvd.state import complete_run, release_all_locks_for_run


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
        "--backup",
        action="store_true",
        default=False,
        help="Create a timestamped backup of the state database after completion",
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="Require state database synchronization (fail if unavailable)",
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

    # Track whether state operations succeeded
    state_available = True

    # Try to mark run as completed
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
        db_path = get_state_db_path(args.state_dir)
        if args.sync:
            raise StateUnavailableError(
                db_path=db_path,
                operation="Marking run as completed",
                reason=str(e),
                original_error=e,
            ) from e

        # Graceful degradation: warn and continue
        state_available = False
        warning = format_state_warning(
            operation="Marking run as completed",
            context=f"Run '{args.run_id}' with status '{args.status}'",
            error=e,
            db_path=db_path,
            consequences=[
                "Run completion status NOT recorded",
                "Run may appear as 'running' in state queries",
            ],
        )
        logger.warning(warning)
        print(f"Warning: Could not mark run as {args.status} (state DB unavailable)")

    # Try to release locks (only if state was available for completion)
    if state_available:
        try:
            locks_released = release_all_locks_for_run(
                run_id=args.run_id,
                state_dir=args.state_dir,
            )
            if locks_released > 0:
                logger.info(f"Released {locks_released} remaining sample lock(s)")
                print(f"Released {locks_released} remaining sample lock(s)")
        except Exception as e:
            # Lock release failure is non-fatal
            logger.warning(f"Failed to release locks (non-fatal): {e}")

    # Create backup if requested (typically on successful completion)
    if args.backup and state_available:
        try:
            backup_path = create_backup(state_dir=args.state_dir)
            logger.info(f"Created backup: {backup_path}")
            print(f"Created backup: {backup_path}")
        except FileNotFoundError as e:
            logger.warning(f"Could not create backup: {e}")
        except Exception as e:
            # Backup failure is non-fatal - the run completed successfully
            logger.warning(f"Backup failed (non-fatal): {e}")


if __name__ == "__main__":
    main()
