#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "blake3",
#     "loguru",
# ]
# ///

"""
Register BLAST hits with idempotent keys in the state database.

This script is a thin wrapper around py_nvd.hits that:
1. Parses CLI arguments
2. Reads FASTA and BLAST result files (I/O)
3. Calls py_nvd.hits API for hit registration

The actual hit key computation, sequence compression, and database operations
are handled by the py_nvd.hits module.

Usage:
    register_hits.py --contigs contigs.fasta --blast-results merged.tsv \\
        --state-dir /path/to/state --sample-set-id abc123 --sample-id sample_a \\
        --run-id my_run_2024

Exit codes:
    0: Success
    1: Error (missing files, database errors, etc.)
"""

import argparse
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import polars as pl
from Bio import SeqIO
from loguru import logger
from py_nvd.hits import record_hit_observation, register_hit
from py_nvd.state import mark_sample_completed


@dataclass(frozen=True)
class HitRegistrationContext:
    """
    Context for registering hits from a single sample.

    Bundles all the parameters needed for hit registration, providing
    validation and type safety over loose string arguments.
    """

    state_dir: Path
    sample_set_id: str
    sample_id: str
    run_date: str

    def __post_init__(self) -> None:
        """Validate the context after initialization."""
        if not self.sample_set_id:
            msg = "sample_set_id cannot be empty"
            raise ValueError(msg)
        if not self.sample_id:
            msg = "sample_id cannot be empty"
            raise ValueError(msg)
        if not self.run_date:
            msg = "run_date cannot be empty"
            raise ValueError(msg)


def parse_fasta(fasta_path: Path) -> dict[str, str]:
    """
    Parse a FASTA file into a dict of {contig_id: sequence}.

    Uses Biopython's SeqIO for robust handling of edge cases.
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}


QSEQID_COLUMN = "qseqid"


def parse_blast_contig_ids(blast_results_path: Path) -> set[str]:
    """
    Extract contig IDs (qseqid) from merged BLAST results TSV.

    Uses Polars for robust TSV parsing with column name validation.
    The qseqid column contains the original contig identifiers from SPAdes.

    Returns set of unique contig IDs that had BLAST hits.

    Raises:
        ValueError: If the qseqid column is not found in the TSV.
    """
    df = pl.read_csv(blast_results_path, separator="\t")

    if QSEQID_COLUMN not in df.columns:
        msg = (
            f"Expected column '{QSEQID_COLUMN}' not found in BLAST results. "
            f"Available columns: {df.columns}"
        )
        raise ValueError(msg)

    return set(df[QSEQID_COLUMN].unique().to_list())


def register_hits_from_contigs(
    contigs: dict[str, str],
    contig_ids_with_hits: set[str],
    context: HitRegistrationContext,
) -> tuple[int, int]:
    """
    Register hits for contigs that had BLAST results.

    For each contig with BLAST hits, registers the hit (if new) and
    records the observation linking it to this sample/run.

    Args:
        contigs: Dict mapping contig IDs to sequences
        contig_ids_with_hits: Set of contig IDs that had BLAST hits
        context: Registration context with sample/run metadata

    Returns:
        Tuple of (new_hits_count, total_observations_count).
    """
    new_hits = 0
    observations = 0

    for contig_id in contig_ids_with_hits:
        if contig_id not in contigs:
            logger.warning(f"Contig {contig_id} in BLAST results but not in FASTA")
            continue

        seq = contigs[contig_id]
        if not seq:
            logger.warning(f"Empty sequence for contig {contig_id}")
            continue

        # Register the hit (idempotent - returns existing if already present)
        hit, is_new = register_hit(
            seq=seq,
            first_seen_date=context.run_date,
            state_dir=context.state_dir,
        )

        if is_new:
            new_hits += 1

        # Record the observation (idempotent - no-op if duplicate)
        observation = record_hit_observation(
            hit_key=hit.hit_key,
            sample_set_id=context.sample_set_id,
            sample_id=context.sample_id,
            run_date=context.run_date,
            contig_id=contig_id,
            state_dir=context.state_dir,
        )

        if observation is not None:
            observations += 1

    return new_hits, observations


def configure_logging(verbosity: int) -> None:
    """
    Configure loguru logging based on verbosity level.

    Args:
        verbosity: 0 = WARNING, 1 = INFO, 2+ = DEBUG
    """
    logger.remove()  # Remove default handler

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
        description="Register BLAST hits with idempotent keys in state database",
    )
    parser.add_argument(
        "--contigs",
        type=Path,
        required=True,
        help="Path to contigs FASTA file",
    )
    parser.add_argument(
        "--blast-results",
        type=Path,
        required=True,
        help="Path to merged BLAST results TSV",
    )
    parser.add_argument(
        "--state-dir",
        type=Path,
        required=True,
        help="Path to state directory containing state.sqlite",
    )
    parser.add_argument(
        "--sample-set-id",
        type=str,
        required=True,
        help="Sample set identifier for this run",
    )
    parser.add_argument(
        "--sample-id",
        type=str,
        required=True,
        help="Sample identifier",
    )
    parser.add_argument(
        "--run-date",
        type=str,
        default=None,
        help="Run date in ISO8601 format (default: now)",
    )
    parser.add_argument(
        "--run-id",
        type=str,
        required=True,
        help="Run identifier (workflow.runName) for state tracking",
    )
    parser.add_argument(
        "--blast-db-version",
        type=str,
        default=None,
        help="BLAST database version for provenance tracking",
    )
    parser.add_argument(
        "--stat-db-version",
        type=str,
        default=None,
        help="STAT database version for provenance tracking",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity (-v for INFO, -vv for DEBUG)",
    )

    args = parser.parse_args()

    # Configure logging
    configure_logging(args.verbose)

    # Validate inputs
    if not args.contigs.exists():
        logger.error(f"Contigs file not found: {args.contigs}")
        sys.exit(1)

    if not args.blast_results.exists():
        logger.error(f"BLAST results file not found: {args.blast_results}")
        sys.exit(1)

    # Build registration context
    run_date = args.run_date or datetime.now(UTC).isoformat()

    try:
        context = HitRegistrationContext(
            state_dir=args.state_dir,
            sample_set_id=args.sample_set_id,
            sample_id=args.sample_id,
            run_date=run_date,
        )
    except ValueError as e:
        logger.error(f"Invalid registration context: {e}")
        sys.exit(1)

    # Parse inputs
    logger.info(f"Parsing contigs from {args.contigs}")
    contigs = parse_fasta(args.contigs)
    logger.info(f"Found {len(contigs)} contigs")

    logger.info(f"Parsing BLAST results from {args.blast_results}")
    contig_ids_with_hits = parse_blast_contig_ids(args.blast_results)
    logger.info(f"Found {len(contig_ids_with_hits)} contigs with hits")

    if not contig_ids_with_hits:
        logger.info("No hits to register")
        new_hits, observations = 0, 0
    else:
        # Register hits
        logger.info(f"Registering hits in {context.state_dir}")
        new_hits, observations = register_hits_from_contigs(
            contigs=contigs,
            contig_ids_with_hits=contig_ids_with_hits,
            context=context,
        )

    logger.info(
        f"Registered {new_hits} new hits, {observations} observations for {context.sample_id}",
    )

    # Mark sample as completed in state database
    # This enables per-sample resume and upload gating
    try:
        mark_sample_completed(
            sample_id=context.sample_id,
            sample_set_id=context.sample_set_id,
            run_id=args.run_id,
            blast_db_version=args.blast_db_version,
            stat_db_version=args.stat_db_version,
            state_dir=str(context.state_dir),
        )
        logger.info(f"Marked sample {context.sample_id} as completed in state database")
    except Exception as e:
        # Log but don't fail - hits are registered, that's the critical part
        # State tracking is for resume/upload gating, not data integrity
        logger.warning(f"Failed to mark sample completed (non-fatal): {e}")


if __name__ == "__main__":
    main()
