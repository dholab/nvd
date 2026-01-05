#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
# ]
# ///

"""
Check and register pipeline run state.

This script enables resumable runs by checking per-sample upload status:
1. Reading sample IDs from a samplesheet CSV
2. Computing a deterministic sample_set_id from the sorted sample list
3. Checking which samples have already been uploaded to LabKey
4. Registering the run for provenance tracking

Key behavior:
- If ALL samples are already uploaded: exit with error (nothing to do)
- If SOME samples are uploaded: warn and continue (Nextflow handles caching)
- If NO samples are uploaded: proceed normally

Individual samples are NOT registered here - they are marked 'completed'
by REGISTER_HITS after processing succeeds, and 'uploaded' after LabKey
upload succeeds.

The sample_set_id is output to stdout for Nextflow to capture.
All logging goes to stderr to avoid polluting the stdout output.

Usage:
    check_run_state.py --samplesheet samples.csv --run-id run_001 \\
        --state-dir /path/to/state

Exit codes:
    0: Success, run registered (or some samples need processing)
    1: Error (all samples already uploaded, no samples found, etc.)
"""

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

from loguru import logger
from py_nvd.state import (
    compute_sample_set_id,
    get_uploaded_sample_ids,
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
    """

    run_id: str
    sample_ids: list[str]
    state_dir: str | None = None
    experiment_id: int | None = None

    @property
    def sample_set_id(self) -> str:
        """Compute deterministic sample set ID from sample list."""
        return compute_sample_set_id(self.sample_ids)

    def check_existing(self) -> list[str]:
        """
        Check sample upload status and determine which samples need processing.

        Returns:
            List of sample IDs that need processing (not yet uploaded)

        Raises:
            SystemExit: If all samples are already uploaded
        """
        logger.info("Checking for previously uploaded samples...")

        uploaded = get_uploaded_sample_ids(self.sample_ids, state_dir=self.state_dir)

        if uploaded:
            not_uploaded = [s for s in self.sample_ids if s not in uploaded]

            if not not_uploaded:
                # ALL samples already uploaded
                logger.error("=" * 60)
                logger.error("ALL SAMPLES ALREADY UPLOADED")
                logger.error("=" * 60)
                logger.error(f"All {len(self.sample_ids)} samples in this set were")
                logger.error("previously uploaded to LabKey.")
                logger.error("=" * 60)
                logger.error("Options:")
                logger.error("  - Use a different sample set")
                logger.error("  - Use 'nvd state samples' to review upload status")
                logger.error("=" * 60)
                sys.exit(1)

            # SOME samples uploaded - warn and continue
            logger.warning("=" * 60)
            logger.warning(
                f"{len(uploaded)} of {len(self.sample_ids)} samples already uploaded:",
            )
            for sid in sorted(uploaded)[:10]:
                logger.warning(f"  - {sid}")
            if len(uploaded) > 10:
                logger.warning(f"  ... and {len(uploaded) - 10} more")
            logger.warning("These samples will be processed but NOT re-uploaded.")
            logger.warning("=" * 60)

            return not_uploaded

        logger.info("No previously uploaded samples - all samples will be processed")
        return self.sample_ids

    def register(self) -> None:
        """
        Register the run in the state database.

        Note: Individual samples are NOT registered here. They are marked
        'completed' by REGISTER_HITS after processing succeeds.

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
            # Run registration failed - likely sample_set_id already exists
            # This is now a warning, not an error, since we check per-sample status
            logger.warning(
                f"Run registration returned None for sample_set_id {self.sample_set_id}",
            )
            logger.warning(
                "This may indicate a concurrent run or previous incomplete run.",
            )
            logger.warning("Proceeding anyway - per-sample status will gate uploads.")
        else:
            logger.debug(f"Run registered successfully: {run}")

        self._log_summary()

    def _log_summary(self) -> None:
        """Log registration summary details."""
        logger.debug("Run registration summary:")
        logger.debug(f"  Run ID: {self.run_id}")
        logger.debug(f"  Sample set ID: {self.sample_set_id}")
        logger.debug(f"  State dir: {self.state_dir or '(default)'}")
        logger.debug(f"  Total samples in set: {len(self.sample_ids)}")
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
        help="BLAST database version (ignored - now passed to register_hits.py)",
    )
    parser.add_argument(
        "--stat-db-version",
        help="STAT database version (ignored - now passed to register_hits.py)",
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
    logger.info("RUN STATE CHECK")
    logger.info("=" * 60)
    logger.debug(f"Samplesheet: {args.samplesheet}")
    logger.debug(f"Run ID: {args.run_id}")
    logger.debug(f"State dir: {args.state_dir or '(default)'}")
    if args.experiment_id:
        logger.debug(f"Experiment ID: {args.experiment_id}")
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
    )

    logger.debug(f"Computed sample_set_id: {registration.sample_set_id}")

    # Check for already-uploaded samples and register run
    samples_to_process = registration.check_existing()
    registration.register()

    logger.info("=" * 60)
    logger.info("RUN STATE CHECK COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Samples to process: {len(samples_to_process)}")
    if len(samples_to_process) < len(sample_ids):
        logger.info(
            f"Samples already uploaded (will skip upload): "
            f"{len(sample_ids) - len(samples_to_process)}",
        )

    # Output sample_set_id to stdout for Nextflow to capture
    # IMPORTANT: end="" prevents trailing newline - stdout emit captures raw text
    print(registration.sample_set_id, end="")  # noqa: T201


if __name__ == "__main__":
    main()
