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

Graceful Degradation:
    By default, if the state database is unavailable (permission errors,
    corruption, etc.), this script will warn and continue without coordination.
    Use --sync to require the state database and fail if unavailable.

Usage:
    check_run_state.py --samplesheet samples.csv --run-id run_001 \\
        --state-dir /path/to/state

    # Require state database (fail if unavailable)
    check_run_state.py --samplesheet samples.csv --run-id run_001 --sync

Exit codes:
    0: Success, run registered (or some samples need processing)
    1: Error (all samples already uploaded, no samples found, etc.)
"""

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path

from loguru import logger
from py_nvd.db import (
    StateUnavailableError,
    format_state_warning,
    get_state_db_path,
    get_taxdump_dir,
)
from py_nvd.state import (
    DEFAULT_LOCK_TTL_HOURS,
    acquire_sample_locks,
    compute_sample_set_id,
    get_uploaded_sample_ids,
    register_run,
    update_run_id,
)
from py_nvd.taxonomy import ensure_taxonomy_available, format_taxonomy_warning


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
        sync: If True, require state DB and fail if unavailable
        state_available: Tracks whether state DB operations succeeded
    """

    run_id: str
    sample_ids: list[str]
    state_dir: str | None = None
    experiment_id: str | None = None
    upload_types: list[str] | None = None
    lock_ttl_hours: int = DEFAULT_LOCK_TTL_HOURS
    sync: bool = False
    state_available: bool = field(default=True, init=False)

    @property
    def sample_set_id(self) -> str:
        """Compute deterministic sample set ID from sample list."""
        return compute_sample_set_id(self.sample_ids)

    def _handle_state_error(
        self,
        error: Exception,
        operation: str,
        context: str,
        consequences: list[str],
    ) -> None:
        """
        Handle a state database error with graceful degradation or failure.

        Args:
            error: The exception that occurred
            operation: Short description of the operation
            context: Detailed context about what was being done
            consequences: List of consequences of continuing without state

        Raises:
            StateUnavailableError: If sync=True
        """
        db_path = get_state_db_path(self.state_dir)

        if self.sync:
            raise StateUnavailableError(
                db_path=db_path,
                operation=operation,
                reason=str(error),
                original_error=error,
            ) from error

        # Graceful degradation: warn and continue
        self.state_available = False
        warning = format_state_warning(
            operation=operation,
            context=context,
            error=error,
            db_path=db_path,
            consequences=consequences,
        )
        logger.warning(warning)

    def check_existing(self) -> list[str]:
        """
        Check sample upload status and determine which samples need processing.

        Returns:
            List of sample IDs that need processing (not yet uploaded)

        Raises:
            SystemExit: If all samples are already uploaded
            StateUnavailableError: If sync=True and state DB unavailable
        """
        logger.info("Checking for previously uploaded samples...")

        try:
            uploaded = get_uploaded_sample_ids(
                self.sample_ids,
                upload_types=self.upload_types,
                state_dir=self.state_dir,
            )
        except Exception as e:
            self._handle_state_error(
                error=e,
                operation="Checking uploaded samples",
                context=f"Run '{self.run_id}' with {len(self.sample_ids)} samples",
                consequences=[
                    "Duplicate detection DISABLED (samples may have been uploaded before)",
                    "All samples will be processed regardless of prior upload status",
                ],
            )
            # Continue without duplicate detection
            return self.sample_ids

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
            (indicating a resume scenario), or None if state DB unavailable.

        Raises:
            StateUnavailableError: If sync=True and state DB unavailable
        """
        if not self.state_available:
            logger.info("Skipping run registration (state database unavailable)")
            return True  # Treat as new run since we can't check

        logger.info(f"Registering run: {self.run_id}")
        logger.debug(f"Sample set ID: {self.sample_set_id}")

        try:
            run = register_run(
                self.run_id,
                self.sample_set_id,
                self.experiment_id,
                state_dir=self.state_dir,
            )
        except Exception as e:
            self._handle_state_error(
                error=e,
                operation="Registering run",
                context=f"Run '{self.run_id}' with sample_set_id '{self.sample_set_id}'",
                consequences=[
                    "Run provenance tracking DISABLED",
                    "Resume detection DISABLED (run may be re-registered on resume)",
                ],
            )
            return True  # Treat as new run

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
        logger.debug(f"  State available: {self.state_available}")
        if self.experiment_id:
            logger.debug(f"  Experiment ID: {self.experiment_id}")

    def acquire_locks(self) -> list[str]:
        """
        Acquire sample locks to prevent duplicate processing.

        Returns:
            List of sample IDs that were successfully locked

        Raises:
            SystemExit: If ALL samples are locked by other runs
            StateUnavailableError: If sync=True and state DB unavailable
        """
        if not self.state_available:
            logger.info("Skipping lock acquisition (state database unavailable)")
            logger.warning("Sample locking DISABLED - concurrent runs may conflict")
            return self.sample_ids  # Return all samples as "locked"

        logger.info("Acquiring sample locks...")

        try:
            acquired, conflicts = acquire_sample_locks(
                self.sample_ids,
                self.run_id,
                ttl_hours=self.lock_ttl_hours,
                state_dir=self.state_dir,
            )
        except Exception as e:
            self._handle_state_error(
                error=e,
                operation="Acquiring sample locks",
                context=f"Run '{self.run_id}' attempting to lock {len(self.sample_ids)} samples",
                consequences=[
                    "Sample locking DISABLED (concurrent runs may process same samples)",
                    "All samples will be processed without coordination",
                ],
            )
            return self.sample_ids  # Return all samples as "locked"

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
                    f"Run '{self.run_id}' is already active on another machine!",
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

    def update_run_id_on_resume(self, locked_samples: list[str]) -> None:
        """
        Update run_id for a resumed run after acquiring locks.

        Args:
            locked_samples: List of samples that were successfully locked
        """
        if not self.state_available:
            logger.debug("Skipping run_id update (state database unavailable)")
            return

        if not locked_samples:
            return

        logger.info(f"Updating run_id to '{self.run_id}' for resumed run")
        try:
            updated_run = update_run_id(
                self.sample_set_id,
                self.run_id,
                state_dir=self.state_dir,
            )
            if updated_run:
                logger.debug(f"Run updated: {updated_run}")
            else:
                # This shouldn't happen if register() returned False
                logger.warning("Failed to update run_id - run may not exist")
        except Exception as e:
            # Non-fatal - log warning but continue
            logger.warning(f"Failed to update run_id (non-fatal): {e}")


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
        "--taxonomy-dir",
        help="Explicit taxonomy directory (takes precedence over --state-dir for taxonomy)",
    )
    parser.add_argument(
        "--experiment-id",
        type=str,
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
        "--upload-types",
        help="Comma-separated upload types to check for duplicate detection "
        "(e.g., 'blast,blast_fasta' or 'gottcha2,gottcha2_fasta'). "
        "When set, only uploads matching these types count as 'already uploaded'.",
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="Require state database synchronization (fail if unavailable)",
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
    logger.debug(f"Sync mode: {args.sync}")
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

    # Ensure taxonomy database is available before workers start
    # This downloads from NCBI if needed (once, on submit node) so workers
    # don't all race to download simultaneously
    logger.info("Ensuring taxonomy database is available...")
    taxdump_dir = get_taxdump_dir(args.state_dir, args.taxonomy_dir)
    if (taxdump_dir / "taxonomy.sqlite").exists():
        logger.info(f"Taxonomy database found: {taxdump_dir}")
    else:
        logger.info("Taxonomy database not found, downloading from NCBI...")
        logger.info("(This may take a few minutes on first run)")
        try:
            ensure_taxonomy_available(args.state_dir, args.taxonomy_dir)
            logger.info(f"Taxonomy database ready: {taxdump_dir}")
        except Exception as e:
            # Taxonomy download failed
            if args.sync:
                logger.error(f"Taxonomy download failed and --sync was specified: {e}")
                sys.exit(1)
            else:
                # Graceful degradation: warn but continue
                # Workers will attempt lazy download when they need taxonomy
                warning = format_taxonomy_warning(
                    operation="Pre-downloading taxonomy database",
                    context="Preparing taxonomy for distributed workers",
                    error=e,
                    taxdump_dir=taxdump_dir,
                    will_download=True,
                )
                logger.warning(warning)
                logger.warning(
                    "Workers will attempt to download taxonomy when needed. "
                    "This may cause delays and potential version inconsistencies.",
                )

    # Create registration object
    upload_types = (
        [t.strip() for t in args.upload_types.split(",") if t.strip()]
        if args.upload_types
        else None
    )

    registration = RunRegistration(
        run_id=args.run_id,
        sample_ids=sample_ids,
        state_dir=args.state_dir,
        experiment_id=args.experiment_id,
        upload_types=upload_types,
        lock_ttl_hours=args.lock_ttl,
        sync=args.sync,
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
    if not is_new_run:
        registration.update_run_id_on_resume(locked_samples)

    # The actual samples to process are those that are:
    # 1. Not already uploaded (samples_to_process)
    # 2. Successfully locked (locked_samples)
    final_samples = [s for s in samples_to_process if s in locked_samples]

    logger.info("=" * 60)
    logger.info("RUN STATE CHECK COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Samples to process: {len(final_samples)}")
    if not registration.state_available:
        logger.warning("State coordination DISABLED - running without state database")
    if len(final_samples) < len(sample_ids):
        skipped_uploaded = len(sample_ids) - len(samples_to_process)
        skipped_locked = len(samples_to_process) - len(final_samples)
        if skipped_uploaded > 0:
            logger.info(
                f"Samples already uploaded (will skip upload): {skipped_uploaded}",
            )
        if skipped_locked > 0:
            logger.info(f"Samples locked by other runs (skipped): {skipped_locked}")

    # Output sample_set_id to stdout for Nextflow to capture
    # IMPORTANT: end="" prevents trailing newline - stdout emit captures raw text
    print(registration.sample_set_id, end="")  # noqa: T201


if __name__ == "__main__":
    main()
