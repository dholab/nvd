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
    acquire_sample_locks,
    compute_sample_set_id,
    DEFAULT_LOCK_TTL_HOURS,
    get_uploaded_sample_ids,
    register_run,
    update_run_id,
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
        lock_ttl_hours: Lock TTL in hours (default: 72)
    """

    run_id: str
    sample_ids: list[str]
    state_dir: str | None = None
    experiment_id: int | None = None
    lock_ttl_hours: int = DEFAULT_LOCK_TTL_HOURS

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

    def register(self) -> bool:
        """
        Register the run in the state database.

        Note: Individual samples are NOT registered here. They are marked
        'completed' by REGISTER_HITS after processing succeeds.

        Returns:
            True if a new run was registered, False if sample_set_id already exists
            (indicating a resume scenario).
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
            # Run registration failed - sample_set_id already exists
            # This is expected on resume; we'll update the run_id after acquiring locks
            logger.info(
                f"Sample set {self.sample_set_id} already has a run registered.",
            )
            logger.info("Will update run_id after acquiring sample locks.")
            self._log_summary()
            return False

        logger.debug(f"Run registered successfully: {run}")
        self._log_summary()
        return True

    def _log_summary(self) -> None:
        """Log registration summary details."""
        logger.debug("Run registration summary:")
        logger.debug(f"  Run ID: {self.run_id}")
        logger.debug(f"  Sample set ID: {self.sample_set_id}")
        logger.debug(f"  State dir: {self.state_dir or '(default)'}")
        logger.debug(f"  Total samples in set: {len(self.sample_ids)}")
        if self.experiment_id:
            logger.debug(f"  Experiment ID: {self.experiment_id}")

    def acquire_locks(self) -> list[str]:
        """
        Acquire sample locks to prevent duplicate processing.

        Returns:
            List of sample IDs that were successfully locked

        Raises:
            SystemExit: If ALL samples are locked by other runs
        """
        logger.info("Acquiring sample locks...")

        acquired, conflicts = acquire_sample_locks(
            self.sample_ids,
            self.run_id,
            ttl_hours=self.lock_ttl_hours,
            state_dir=self.state_dir,
        )

        if conflicts:
            # Group conflicts by run_id for cleaner logging
            by_run: dict[str, list[str]] = {}
            for lock in conflicts:
                by_run.setdefault(lock.run_id, []).append(lock.sample_id)

            # Check if this is a fingerprint conflict (same run_id, different machine)
            same_run_conflicts = [c for c in conflicts if c.run_id == self.run_id]
            if same_run_conflicts:
                lock = same_run_conflicts[0]
                logger.error("=" * 60)
                logger.error("RUN CONFLICT DETECTED")
                logger.error("=" * 60)
                logger.error(
                    f"Run '{self.run_id}' is already active on another machine!"
                )
                logger.error(f"  Hostname: {lock.hostname}")
                logger.error(f"  Username: {lock.username}")
                logger.error(f"  Locked at: {lock.locked_at}")
                logger.error(f"  Expires at: {lock.expires_at}")
                logger.error("=" * 60)
                logger.error("This run_id cannot be resumed from this machine.")
                logger.error("If the other run crashed, wait for lock expiry or use:")
                logger.error(f"  nvd state unlock --run-id {self.run_id} --force")
                logger.error("=" * 60)
                sys.exit(1)

            # Different runs hold the locks
            logger.warning("=" * 60)
            logger.warning(f"{len(conflicts)} samples locked by other runs:")
            for run_id, samples in by_run.items():
                logger.warning(f"  Run '{run_id}': {len(samples)} samples")
                for sid in samples[:5]:
                    logger.warning(f"    - {sid}")
                if len(samples) > 5:
                    logger.warning(f"    ... and {len(samples) - 5} more")
            logger.warning("These samples will be SKIPPED in this run.")
            logger.warning("=" * 60)

        if not acquired:
            logger.error("=" * 60)
            logger.error("ALL SAMPLES LOCKED")
            logger.error("=" * 60)
            logger.error("All samples are currently being processed by other runs.")
            logger.error("No samples available to process.")
            logger.error("=" * 60)
            sys.exit(1)

        logger.info(f"Acquired locks for {len(acquired)} samples")
        if conflicts:
            logger.info(f"Skipping {len(conflicts)} samples (locked by other runs)")

        return acquired


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
        "--lock-ttl",
        type=int,
        default=DEFAULT_LOCK_TTL_HOURS,
        help=f"Lock TTL in hours (default: {DEFAULT_LOCK_TTL_HOURS})",
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
        lock_ttl_hours=args.lock_ttl,
    )

    logger.debug(f"Computed sample_set_id: {registration.sample_set_id}")

    # Check for already-uploaded samples and register run
    samples_to_process = registration.check_existing()
    is_new_run = registration.register()

    # Acquire locks for samples we'll process
    locked_samples = registration.acquire_locks()

    # If this is a resume (run already existed), update the run_id now that
    # we've successfully acquired locks. This keeps the runs table in sync
    # with the sample_locks table, which updates run_id on re-acquisition.
    if not is_new_run and locked_samples:
        logger.info(f"Updating run_id to '{registration.run_id}' for resumed run")
        updated_run = update_run_id(
            registration.sample_set_id,
            registration.run_id,
            state_dir=registration.state_dir,
        )
        if updated_run:
            logger.debug(f"Run updated: {updated_run}")
        else:
            # This shouldn't happen if register() returned False
            logger.warning("Failed to update run_id - run may not exist")

    # The actual samples to process are those that are:
    # 1. Not already uploaded (samples_to_process)
    # 2. Successfully locked (locked_samples)
    final_samples = [s for s in samples_to_process if s in locked_samples]

    logger.info("=" * 60)
    logger.info("RUN STATE CHECK COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Samples to process: {len(final_samples)}")
    if len(final_samples) < len(sample_ids):
        skipped_uploaded = len(sample_ids) - len(samples_to_process)
        skipped_locked = len(samples_to_process) - len(final_samples)
        if skipped_uploaded > 0:
            logger.info(
                f"Samples already uploaded (will skip upload): {skipped_uploaded}"
            )
        if skipped_locked > 0:
            logger.info(f"Samples locked by other runs (skipped): {skipped_locked}")

    # Output sample_set_id to stdout for Nextflow to capture
    # IMPORTANT: end="" prevents trailing newline - stdout emit captures raw text
    print(registration.sample_set_id, end="")  # noqa: T201


if __name__ == "__main__":
    main()
