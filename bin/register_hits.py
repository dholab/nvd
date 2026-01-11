#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "biopython",
#     "blake3",
#     "loguru",
#     "polars",
# ]
# ///

"""
Register BLAST hits to parquet files for crash-safe storage.

This script is a thin wrapper around py_nvd.hits that:
1. Parses CLI arguments
2. Reads FASTA and BLAST result files (I/O)
3. Computes hit keys and metadata
4. Writes hits atomically to a per-sample parquet file

Each sample's hits are stored in a Hive-partitioned structure:
    {state_dir}/hits/month=NULL/{sample_set_id}/{sample_id}/data.parquet

The month=NULL partition holds uncompacted data. After running `nvd hits compact`,
data moves to month=YYYY-MM partitions for better query performance.

Usage:
    register_hits.py --contigs contigs.fasta --blast-results merged.tsv \\
        --state-dir /path/to/state --sample-set-id abc123 --sample-id sample_a \\
        --run-id my_run_2024

Exit codes:
    0: Success
    1: Error (missing files, write errors, etc.)
"""

import argparse
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import polars as pl
from Bio import SeqIO
from loguru import logger
from py_nvd.hits import (
    HitRecord,
    calculate_gc_content,
    compress_sequence,
    compute_hit_key,
    write_hits_parquet,
)
from py_nvd.state import mark_sample_completed, release_sample_lock


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


def build_hit_records(
    contigs: dict[str, str],
    contig_ids_with_hits: set[str],
    context: HitRegistrationContext,
) -> list[HitRecord]:
    """
    Build HitRecord objects for contigs that had BLAST results.

    Args:
        contigs: Dict mapping contig IDs to sequences
        contig_ids_with_hits: Set of contig IDs that had BLAST hits
        context: Registration context with sample/run metadata

    Returns:
        List of HitRecord objects ready to be written.
    """
    hit_records: list[HitRecord] = []

    for contig_id in contig_ids_with_hits:
        if contig_id not in contigs:
            logger.warning(f"Contig {contig_id} in BLAST results but not in FASTA")
            continue

        seq = contigs[contig_id]
        if not seq:
            logger.warning(f"Empty sequence for contig {contig_id}")
            continue

        hit_key = compute_hit_key(seq)
        compressed = compress_sequence(seq)
        gc_content = calculate_gc_content(seq)

        hit_records.append(
            HitRecord(
                hit_key=hit_key,
                sequence_length=len(seq),
                sequence_compressed=compressed,
                gc_content=gc_content,
                sample_set_id=context.sample_set_id,
                sample_id=context.sample_id,
                run_date=context.run_date,
                contig_id=contig_id,
            )
        )

    return hit_records


def write_hits_to_path(hits: list[HitRecord], output_path: Path) -> Path:
    """
    Write hit records to a specific parquet file path.

    Uses atomic write (temp file + rename) for crash safety.

    Args:
        hits: List of HitRecord objects to write
        output_path: Destination path for the parquet file

    Returns:
        Path to the written parquet file
    """
    # Convert to polars DataFrame
    if not hits:
        # Create empty DataFrame with correct schema
        df = pl.DataFrame(
            schema={
                "hit_key": pl.Utf8,
                "sequence_length": pl.Int64,
                "sequence_compressed": pl.Binary,
                "gc_content": pl.Float64,
                "sample_set_id": pl.Utf8,
                "sample_id": pl.Utf8,
                "run_date": pl.Utf8,
                "contig_id": pl.Utf8,
            }
        )
    else:
        df = pl.DataFrame(
            [
                {
                    "hit_key": h.hit_key,
                    "sequence_length": h.sequence_length,
                    "sequence_compressed": h.sequence_compressed,
                    "gc_content": h.gc_content,
                    "sample_set_id": h.sample_set_id,
                    "sample_id": h.sample_id,
                    "run_date": h.run_date,
                    "contig_id": h.contig_id,
                }
                for h in hits
            ]
        )

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Atomic write: write to temp file, then rename
    tmp_path = output_path.with_suffix(".parquet.tmp")
    df.write_parquet(tmp_path)
    tmp_path.rename(output_path)

    return output_path


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
        "--labkey",
        action="store_true",
        default=False,
        help="LabKey integration is enabled (lock released after upload, not here)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Write parquet to this path instead of state directory (for Nextflow publishDir)",
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

    # Build hit records
    hit_records = build_hit_records(
        contigs=contigs,
        contig_ids_with_hits=contig_ids_with_hits,
        context=context,
    )
    logger.info(f"Computed {len(hit_records)} hit records")

    # Write parquet file
    if args.output:
        # Write to explicit output path (for Nextflow publishDir)
        parquet_path = write_hits_to_path(hit_records, args.output)
        logger.info(f"Wrote {len(hit_records)} hits to {parquet_path}")
    else:
        # Write to state directory (legacy behavior)
        parquet_path = write_hits_parquet(
            hits=hit_records,
            sample_id=context.sample_id,
            sample_set_id=context.sample_set_id,
            state_dir=context.state_dir,
        )
        logger.info(f"Wrote {len(hit_records)} hits to {parquet_path}")

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

    # Release sample lock if LabKey is disabled
    # When LabKey is enabled, locks are released after successful upload
    if not args.labkey:
        try:
            released = release_sample_lock(
                sample_id=context.sample_id,
                run_id=args.run_id,
                state_dir=str(context.state_dir),
            )
            if released:
                logger.info(f"Released lock for sample {context.sample_id}")
            else:
                logger.debug(f"No lock to release for sample {context.sample_id}")
        except Exception as e:
            logger.warning(f"Failed to release sample lock (non-fatal): {e}")


if __name__ == "__main__":
    main()
