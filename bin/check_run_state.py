#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
# ]
# ///

"""
Check and register pipeline run state.

This script prevents duplicate processing of the same sample set by:
1. Reading sample IDs from a samplesheet CSV
2. Computing a deterministic sample_set_id from the sorted sample list
3. Checking if this sample set was already processed
4. Registering the run and samples with provenance metadata

The sample_set_id is output to stdout for Nextflow to capture.
All logging goes to stderr to avoid polluting the stdout output.

Usage:
    check_run_state.py --samplesheet samples.csv --run-id run_001 \\
        --state-dir /path/to/state

Exit codes:
    0: Success, run registered
    1: Error (sample set already processed, no samples found, etc.)
"""

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

from loguru import logger
from py_nvd.state import (
    compute_sample_set_id,
    get_run_by_sample_set,
    register_processed_sample,
    register_run,
)


@dataclass
class RunRegistration:
    """
    Bundles run registration data and provides registration logic.

    Attributes:
        run_id: Unique run identifier (e.g., workflow.runName)
        sample_ids: List of sample IDs in this run
        state_dir: Path to state directory (None for default)
        experiment_id: Optional experiment ID for tracking
        blast_db_version: Optional BLAST database version for provenance
        stat_db_version: Optional STAT database version for provenance
    """

    run_id: str
    sample_ids: list[str]
    state_dir: str | None = None
    experiment_id: int | None = None
    blast_db_version: str | None = None
    stat_db_version: str | None = None

    @property
    def sample_set_id(self) -> str:
        """Compute deterministic sample set ID from sample list."""
        return compute_sample_set_id(self.sample_ids)

    def check_existing(self) -> None:
        """
        Check if this sample set was already processed.

        Raises:
            SystemExit: If sample set was already processed
        """
        logger.info("Checking for existing runs with this sample set...")
        logger.debug(f"Querying state database for sample_set_id: {self.sample_set_id}")

        existing = get_run_by_sample_set(self.sample_set_id, state_dir=self.state_dir)

        if existing:
            logger.error("=" * 60)
            logger.error("DUPLICATE SAMPLE SET DETECTED")
            logger.error("=" * 60)
            logger.error("This sample set was already processed in a previous run.")
            logger.error(f"  Previous run ID: {existing.run_id}")
            logger.error(f"  Sample set ID: {self.sample_set_id}")
            logger.error(f"  Status: {existing.status}")
            if existing.experiment_id:
                logger.error(f"  Experiment ID: {existing.experiment_id}")
            logger.error("=" * 60)
            logger.error("Options:")
            logger.error("  - Use a different sample set")
            logger.error("  - Clear the state database if re-processing is intended")
            logger.error("=" * 60)
            sys.exit(1)

        logger.info("No existing run found - safe to proceed")

    def register(self) -> None:
        """
        Register the run and all samples in the state database.

        Raises:
            SystemExit: If run registration fails (e.g., race condition)
        """
        logger.info(f"Registering run: {self.run_id}")
        logger.debug(f"Sample set ID: {self.sample_set_id}")

        run = register_run(
            self.run_id,
            self.sample_set_id,
            self.experiment_id,
            state_dir=self.state_dir,
        )

        if not run:
            logger.error("=" * 60)
            logger.error("RUN REGISTRATION FAILED")
            logger.error("=" * 60)
            logger.error(
                f"Sample set {self.sample_set_id} may have been claimed by a concurrent run.",
            )
            logger.error("This can happen when multiple pipelines start simultaneously")
            logger.error("with overlapping sample sets.")
            logger.error("=" * 60)
            sys.exit(1)

        logger.debug(f"Run registered successfully: {run}")

        # Register each sample with provenance
        logger.info(
            f"Registering {len(self.sample_ids)} samples with provenance metadata..."
        )
        for i, sid in enumerate(self.sample_ids, 1):
            logger.debug(f"  [{i}/{len(self.sample_ids)}] Registering sample: {sid}")
            register_processed_sample(
                sample_id=sid,
                sample_set_id=self.sample_set_id,
                run_id=self.run_id,
                blast_db_version=self.blast_db_version,
                stat_db_version=self.stat_db_version,
                state_dir=self.state_dir,
            )

        logger.info(f"All {len(self.sample_ids)} samples registered")
        self._log_summary()

    def _log_summary(self) -> None:
        """Log registration summary details."""
        logger.debug("Registration summary:")
        logger.debug(f"  Run ID: {self.run_id}")
        logger.debug(f"  Sample set ID: {self.sample_set_id}")
        logger.debug(f"  State dir: {self.state_dir or '(default)'}")
        logger.debug(f"  Sample count: {len(self.sample_ids)}")
        logger.debug(f"  BLAST DB: {self.blast_db_version or '(not specified)'}")
        logger.debug(f"  STAT DB: {self.stat_db_version or '(not specified)'}")
        if self.experiment_id:
            logger.debug(f"  Experiment ID: {self.experiment_id}")


def configure_logging(verbose: bool = False) -> None:  # noqa: FBT001, FBT002
    """Configure loguru to write to stderr only."""
    logger.remove()
    logger.add(
        sys.stderr,
        format="<level>{level: <8}</level> | {message}",
        level="DEBUG" if verbose else "INFO",
    )


def parse_samplesheet(samplesheet_path: Path) -> list[str]:
    """
    Parse sample IDs from a samplesheet CSV.

    Args:
        samplesheet_path: Path to the samplesheet CSV file

    Returns:
        Sorted, deduplicated list of sample IDs

    Raises:
        ValueError: If no valid sample IDs are found
    """
    sample_ids = []

    with open(samplesheet_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row.get("sample_id", "").strip()
            if sample_id and not sample_id.startswith("#"):
                sample_ids.append(sample_id)

    if not sample_ids:
        msg = f"No sample IDs found in samplesheet: {samplesheet_path}"
        raise ValueError(msg)

    # Deduplicate and sort for deterministic sample_set_id
    return sorted(set(sample_ids))


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Check and register pipeline run state.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--samplesheet",
        type=Path,
        required=True,
        help="Path to samplesheet CSV file",
    )
    parser.add_argument(
        "--run-id",
        required=True,
        help="Unique run identifier (e.g., workflow.runName)",
    )
    parser.add_argument(
        "--state-dir",
        help="State directory for SQLite database (default: ~/.nvd/ or NVD_STATE_DIR)",
    )
    parser.add_argument(
        "--experiment-id",
        type=int,
        help="Experiment ID for tracking",
    )
    parser.add_argument(
        "--blast-db-version",
        help="BLAST database version for provenance",
    )
    parser.add_argument(
        "--stat-db-version",
        help="STAT database version for provenance",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug-level logging",
    )

    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    configure_logging(verbose=args.verbose)

    # Log header
    logger.info("=" * 60)
    logger.info("RUN STATE CHECK AND REGISTRATION")
    logger.info("=" * 60)
    logger.debug(f"Samplesheet: {args.samplesheet}")
    logger.debug(f"Run ID: {args.run_id}")
    logger.debug(f"State dir: {args.state_dir or '(default)'}")
    if args.experiment_id:
        logger.debug(f"Experiment ID: {args.experiment_id}")
    if args.blast_db_version:
        logger.debug(f"BLAST DB version: {args.blast_db_version}")
    if args.stat_db_version:
        logger.debug(f"STAT DB version: {args.stat_db_version}")
    logger.info("=" * 60)

    # Parse samplesheet
    logger.info(f"Parsing samplesheet: {args.samplesheet}")
    try:
        sample_ids = parse_samplesheet(args.samplesheet)
    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)
    except FileNotFoundError:
        logger.error(f"Samplesheet not found: {args.samplesheet}")
        sys.exit(1)

    logger.info(f"Found {len(sample_ids)} unique samples")
    logger.debug(f"Sample IDs: {sample_ids}")

    # Create registration object
    registration = RunRegistration(
        run_id=args.run_id,
        sample_ids=sample_ids,
        state_dir=args.state_dir,
        experiment_id=args.experiment_id,
        blast_db_version=args.blast_db_version,
        stat_db_version=args.stat_db_version,
    )

    logger.debug(f"Computed sample_set_id: {registration.sample_set_id}")

    # Check for existing run and register
    registration.check_existing()
    registration.register()

    logger.info("=" * 60)
    logger.info("REGISTRATION SUCCESSFUL")
    logger.info("=" * 60)

    # Output sample_set_id to stdout for Nextflow to capture
    # IMPORTANT: end="" prevents trailing newline - stdout emit captures raw text
    print(registration.sample_set_id, end="")  # noqa: T201


if __name__ == "__main__":
    main()
