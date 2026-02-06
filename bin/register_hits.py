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

Graceful Degradation:
    By default, if the state database is unavailable for marking samples as
    completed or releasing locks, this script will warn and continue. The
    parquet file (the critical output) is still written. Use --sync to require
    state database operations and fail if unavailable.

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
from py_nvd.db import StateUnavailableError, format_state_warning, get_state_db_path
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

    When state_dir is None (stateless mode), state operations are skipped
    and --output must be specified for parquet output.
    """

    state_dir: Path | None
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


@dataclass(frozen=True)
class ContigClassification:
    """
    Classification data for a single contig from LCA-annotated BLAST results.

    Each contig gets a consensus taxonomic assignment (either dominant taxid
    or LCA of near-ties). This data is extracted from the merged BLAST results
    after LCA annotation.
    """

    contig_id: str
    adjusted_taxid: int | None
    adjusted_taxid_name: str | None
    adjusted_taxid_rank: str | None

    def __post_init__(self) -> None:
        assert self.contig_id, "contig_id cannot be empty"
        assert isinstance(self.contig_id, str), (
            f"contig_id must be str, got {type(self.contig_id)}"
        )


def parse_blast_classifications(
    blast_results_path: Path,
) -> dict[str, ContigClassification]:
    """
    Extract classification data per contig from LCA-annotated BLAST results TSV.

    Uses Polars for robust TSV parsing. The input file has multiple rows per contig
    (one per BLAST hit), but they share the same adjusted_* values from LCA analysis.
    We take the first row per contig to get the classification.

    Args:
        blast_results_path: Path to LCA-annotated merged BLAST results TSV

    Returns:
        Dict mapping contig_id to ContigClassification

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

    # Get one row per contig (they all have the same adjusted_* values)
    # Use first() aggregation to pick one representative row
    unique_contigs = df.group_by(QSEQID_COLUMN).agg(
        pl.col("adjusted_taxid").first(),
        pl.col("adjusted_taxid_name").first(),
        pl.col("adjusted_taxid_rank").first(),
    )

    classifications: dict[str, ContigClassification] = {}
    for row in unique_contigs.iter_rows(named=True):
        contig_id = row[QSEQID_COLUMN]

        # Handle type conversion for adjusted_taxid (may be float from CSV parsing)
        raw_taxid = row["adjusted_taxid"]
        if raw_taxid is not None:
            try:
                adjusted_taxid = int(raw_taxid)
            except (ValueError, TypeError):
                adjusted_taxid = None
        else:
            adjusted_taxid = None

        classifications[contig_id] = ContigClassification(
            contig_id=contig_id,
            adjusted_taxid=adjusted_taxid,
            adjusted_taxid_name=row["adjusted_taxid_name"],
            adjusted_taxid_rank=row["adjusted_taxid_rank"],
        )

    return classifications


def build_hit_records(
    contigs: dict[str, str],
    classifications: dict[str, ContigClassification],
    context: HitRegistrationContext,
) -> list[HitRecord]:
    """
    Build HitRecord objects for contigs that had BLAST results.

    Args:
        contigs: Dict mapping contig IDs to sequences
        classifications: Dict mapping contig IDs to their classification data
        context: Registration context with sample/run metadata

    Returns:
        List of HitRecord objects ready to be written.
    """
    hit_records: list[HitRecord] = []

    for contig_id, classification in classifications.items():
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
                adjusted_taxid=classification.adjusted_taxid,
                adjusted_taxid_name=classification.adjusted_taxid_name,
                adjusted_taxid_rank=classification.adjusted_taxid_rank,
            ),
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
                "adjusted_taxid": pl.Int64,
                "adjusted_taxid_name": pl.Utf8,
                "adjusted_taxid_rank": pl.Utf8,
            },
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
                    "adjusted_taxid": h.adjusted_taxid,
                    "adjusted_taxid_name": h.adjusted_taxid_name,
                    "adjusted_taxid_rank": h.adjusted_taxid_rank,
                }
                for h in hits
            ],
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
        required=False,
        default=None,
        help="Path to state directory containing state.sqlite (optional in stateless mode)",
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

    # Configure logging
    configure_logging(args.verbose)

    # Validate inputs
    if not args.contigs.exists():
        logger.error(f"Contigs file not found: {args.contigs}")
        sys.exit(1)

    if not args.blast_results.exists():
        logger.error(f"BLAST results file not found: {args.blast_results}")
        sys.exit(1)

    # In stateless mode (no state_dir), --output is required
    if args.state_dir is None and args.output is None:
        logger.error(
            "--output is required when --state-dir is not provided (stateless mode)",
        )
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
    classifications = parse_blast_classifications(args.blast_results)
    logger.info(f"Found {len(classifications)} contigs with hits")

    # Build hit records
    hit_records = build_hit_records(
        contigs=contigs,
        classifications=classifications,
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
    # Skip in stateless mode (no state_dir)
    if context.state_dir is not None:
        try:
            mark_sample_completed(
                sample_id=context.sample_id,
                sample_set_id=context.sample_set_id,
                run_id=args.run_id,
                blast_db_version=args.blast_db_version,
                stat_db_version=args.stat_db_version,
                state_dir=str(context.state_dir),
            )
            logger.info(
                f"Marked sample {context.sample_id} as completed in state database",
            )
        except Exception as e:
            db_path = get_state_db_path(context.state_dir)
            if args.sync:
                raise StateUnavailableError(
                    db_path=db_path,
                    operation="Marking sample as completed",
                    reason=str(e),
                    original_error=e,
                ) from e

            # Graceful degradation: warn but don't fail
            # Hits are registered (parquet written), that's the critical part
            # State tracking is for resume/upload gating, not data integrity
            warning = format_state_warning(
                operation="Marking sample as completed",
                context=f"Sample '{context.sample_id}' in run '{args.run_id}'",
                error=e,
                db_path=db_path,
                consequences=[
                    "Sample completion status NOT recorded",
                    "Sample may be re-processed on pipeline resume",
                    "Upload gating may not work correctly",
                ],
            )
            logger.warning(warning)
    else:
        logger.debug("Skipping state operations (stateless mode)")

    # Release sample lock if LabKey is disabled
    # When LabKey is enabled, locks are released after successful upload
    # Skip in stateless mode (no state_dir)
    if not args.labkey and context.state_dir is not None:
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
            db_path = get_state_db_path(context.state_dir)
            if args.sync:
                raise StateUnavailableError(
                    db_path=db_path,
                    operation="Releasing sample lock",
                    reason=str(e),
                    original_error=e,
                ) from e

            # Graceful degradation: warn but don't fail
            warning = format_state_warning(
                operation="Releasing sample lock",
                context=f"Sample '{context.sample_id}' in run '{args.run_id}'",
                error=e,
                db_path=db_path,
                consequences=[
                    "Sample lock NOT released",
                    "Lock will expire after TTL (default: 72 hours)",
                    "Other runs may be blocked from processing this sample",
                ],
            )
            logger.warning(warning)


if __name__ == "__main__":
    main()
