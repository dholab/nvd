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
Register GOTTCHA2 hits to parquet files for crash-safe storage.

This script is the GOTTCHA2 counterpart to register_hits.py (BLAST). It:
1. Parses CLI arguments
2. Reads the extracted FASTA from GOTTCHA2 (post-SANITIZE_EXTRACTED_FASTA)
3. Extracts taxonomy from FASTA headers (LEVEL, NAME, TAXID)
4. Computes hit keys and metadata
5. Writes hits atomically to a per-sample parquet file

GOTTCHA2 FASTA headers encode taxonomy directly:
    >{read_name}{mate}|{ref}:{start}..{end} LEVEL={level} NAME={name} TAXID={taxid}

Names use underscores in place of spaces (restored during parsing).

Each sample's hits are stored in a Hive-partitioned structure:
    {state_dir}/hits/month=NULL/{sample_set_id}/{sample_id}/data.parquet

Graceful Degradation:
    By default, if the state database is unavailable for marking samples as
    completed or releasing locks, this script will warn and continue. The
    parquet file (the critical output) is still written. Use --sync to require
    state database operations and fail if unavailable.

Usage:
    register_gottcha2_hits.py --fasta extracted.fasta --full-tsv profile.tsv \\
        --state-dir /path/to/state --sample-set-id abc123 --sample-id sample_a \\
        --run-id my_run_2024

Exit codes:
    0: Success
    1: Error (missing files, write errors, etc.)
"""

import argparse
import re
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

# Regex for parsing GOTTCHA2 FASTA header metadata.
# Matches: LEVEL={level} NAME={name} TAXID={taxid}
# NAME may contain underscores (spaces replaced by GOTTCHA2).
HEADER_PATTERN = re.compile(
    r"LEVEL=(?P<level>\S+)\s+NAME=(?P<name>\S+)\s+TAXID=(?P<taxid>\d+)",
)


@dataclass(frozen=True)
class Gottcha2HitRegistrationContext:
    """
    Context for registering GOTTCHA2 hits from a single sample.

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


@dataclass(frozen=True)
class Gottcha2Classification:
    """
    Classification data for a single sequence from GOTTCHA2 extracted FASTA.

    Taxonomy is parsed directly from the FASTA header, unlike BLAST where
    it comes from a separate LCA-annotated TSV.
    """

    sequence_id: str
    taxid: int
    name: str
    rank: str

    def __post_init__(self) -> None:
        assert self.sequence_id, "sequence_id cannot be empty"
        assert isinstance(self.sequence_id, str), (
            f"sequence_id must be str, got {type(self.sequence_id)}"
        )


def parse_gottcha2_fasta(
    fasta_path: Path,
) -> list[tuple[str, str, Gottcha2Classification]]:
    """
    Parse a GOTTCHA2 extracted FASTA file into sequences with classifications.

    Each FASTA entry has taxonomy encoded in its header:
        >{read_name}{mate}|{ref}:{start}..{end} LEVEL={level} NAME={name} TAXID={taxid}

    Names have underscores replacing spaces (restored here).

    Args:
        fasta_path: Path to the GOTTCHA2 extracted FASTA file

    Returns:
        List of (sequence_id, sequence, classification) tuples.
        Entries with unparseable headers are logged and skipped.
    """
    results: list[tuple[str, str, Gottcha2Classification]] = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        if not seq:
            logger.warning(f"Empty sequence for {record.id}, skipping")
            continue

        # The full description line (after >) includes the ID and metadata
        description = record.description
        match = HEADER_PATTERN.search(description)

        if match is None:
            logger.warning(
                f"Could not parse GOTTCHA2 header for {record.id}, skipping: "
                f"{description}",
            )
            continue

        taxid_str = match.group("taxid")
        name_raw = match.group("name")
        level = match.group("level")

        # Restore spaces from underscores in taxon name
        name = name_raw.replace("_", " ")

        classification = Gottcha2Classification(
            sequence_id=record.id,
            taxid=int(taxid_str),
            name=name,
            rank=level,
        )

        results.append((record.id, seq, classification))

    return results


# Priority ranking for taxonomic levels: lower value = more specific.
# GOTTCHA2's -ef extraction duplicates each read across all levels; we keep
# only the most specific classification per read (strain > species).
_LEVEL_PRIORITY: dict[str, int] = {"strain": 0, "species": 1}


def dedupe_by_taxonomic_level(
    parsed_entries: list[tuple[str, str, Gottcha2Classification]],
) -> list[tuple[str, str, Gottcha2Classification]]:
    """
    Collapse parsed FASTA entries to one taxonomic level per read.

    GOTTCHA2's ``-ef`` extraction emits each aligned read once per taxonomic
    level (strain, species, genus, family, order, class, phylum, superkingdom),
    inflating the entry count up to 8x. This function keeps only the most
    specific classification per read:

        1. STRAIN entries present  -> keep only STRAIN
        2. No strain, has SPECIES  -> keep only SPECIES
        3. Neither                 -> discard (no informative sub-species classification)

    Read identity is the token before the first ``|`` in the sequence ID::

        {read_name}{mate}|{ref}:{start}..{end}

    The same read appears at multiple taxonomic levels with an identical key,
    so this token serves as the grouping key for deduplication.

    This mirrors the selection logic in ``labkey_upload_gottcha2_fasta.py``
    (``select_records_for_upload``) so that hit registration and LabKey uploads
    agree on which entries represent each read.

    Args:
        parsed_entries: List of ``(sequence_id, sequence, classification)``
            tuples from :func:`parse_gottcha2_fasta`.

    Returns:
        Filtered list with at most one taxonomic level per read.
    """
    per_read: dict[str, dict[str, list[tuple[str, str, Gottcha2Classification]]]] = {}

    for entry in parsed_entries:
        sequence_id, _seq, classification = entry
        read_key = sequence_id.split("|", 1)[0]
        rank = classification.rank.lower()

        if rank not in _LEVEL_PRIORITY:
            continue

        bucket = per_read.setdefault(read_key, {})
        bucket.setdefault(rank, []).append(entry)

    selected: list[tuple[str, str, Gottcha2Classification]] = []
    for levels in per_read.values():
        best_rank = min(levels, key=lambda r: _LEVEL_PRIORITY[r])
        selected.extend(levels[best_rank])

    return selected


def build_gottcha2_hit_records(
    parsed_entries: list[tuple[str, str, Gottcha2Classification]],
    context: Gottcha2HitRegistrationContext,
) -> list[HitRecord]:
    """
    Build HitRecord objects from parsed GOTTCHA2 FASTA entries.

    Args:
        parsed_entries: List of (sequence_id, sequence, classification) tuples
            from parse_gottcha2_fasta()
        context: Registration context with sample/run metadata

    Returns:
        List of HitRecord objects ready to be written.
    """
    hit_records: list[HitRecord] = []

    for sequence_id, seq, classification in parsed_entries:
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
                sequence_id=sequence_id,
                source="gottcha2",
                adjusted_taxid=classification.taxid,
                adjusted_taxid_name=classification.name,
                adjusted_taxid_rank=classification.rank,
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
        # Create empty DataFrame with correct schema.
        # Includes contig_id=NULL for backward compatibility (see COALESCE strategy).
        df = pl.DataFrame(
            schema={
                "hit_key": pl.Utf8,
                "sequence_length": pl.Int64,
                "sequence_compressed": pl.Binary,
                "gc_content": pl.Float64,
                "sample_set_id": pl.Utf8,
                "sample_id": pl.Utf8,
                "run_date": pl.Utf8,
                "sequence_id": pl.Utf8,
                "contig_id": pl.Utf8,
                "source": pl.Utf8,
                "adjusted_taxid": pl.Int64,
                "adjusted_taxid_name": pl.Utf8,
                "adjusted_taxid_rank": pl.Utf8,
            },
        )
    else:
        n = len(hits)
        df = pl.DataFrame(
            {
                "hit_key": [h.hit_key for h in hits],
                "sequence_length": [h.sequence_length for h in hits],
                "sequence_compressed": [h.sequence_compressed for h in hits],
                "gc_content": [h.gc_content for h in hits],
                "sample_set_id": [h.sample_set_id for h in hits],
                "sample_id": [h.sample_id for h in hits],
                "run_date": [h.run_date for h in hits],
                "sequence_id": pl.Series(
                    "sequence_id",
                    [h.sequence_id for h in hits],
                    dtype=pl.Utf8,
                ),
                "contig_id": pl.Series("contig_id", [None] * n, dtype=pl.Utf8),
                "source": [h.source for h in hits],
                "adjusted_taxid": [h.adjusted_taxid for h in hits],
                "adjusted_taxid_name": [h.adjusted_taxid_name for h in hits],
                "adjusted_taxid_rank": [h.adjusted_taxid_rank for h in hits],
            },
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
        description="Register GOTTCHA2 hits with idempotent keys in state database",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to GOTTCHA2 extracted FASTA file (post-SANITIZE_EXTRACTED_FASTA)",
    )
    parser.add_argument(
        "--full-tsv",
        type=Path,
        required=True,
        help="Path to GOTTCHA2 full profile TSV (retained for provenance)",
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
    if not args.fasta.exists():
        logger.error(f"FASTA file not found: {args.fasta}")
        sys.exit(1)

    if not args.full_tsv.exists():
        logger.error(f"Full TSV file not found: {args.full_tsv}")
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
        context = Gottcha2HitRegistrationContext(
            state_dir=args.state_dir,
            sample_set_id=args.sample_set_id,
            sample_id=args.sample_id,
            run_date=run_date,
        )
    except ValueError as e:
        logger.error(f"Invalid registration context: {e}")
        sys.exit(1)

    # Parse GOTTCHA2 FASTA with embedded taxonomy
    logger.info(f"Parsing GOTTCHA2 FASTA from {args.fasta}")
    parsed_entries = parse_gottcha2_fasta(args.fasta)
    logger.info(f"Found {len(parsed_entries)} sequences with taxonomy")

    # Collapse to one taxonomic level per read (strain > species, discard higher).
    # GOTTCHA2 -ef duplicates each read across up to 8 levels; this brings parity
    # with what labkey_upload_gottcha2_fasta.py sends to LabKey.
    deduped_entries = dedupe_by_taxonomic_level(parsed_entries)
    removed = len(parsed_entries) - len(deduped_entries)
    logger.info(
        f"After taxonomic dedup: {len(deduped_entries)} entries "
        f"({removed} duplicates removed)"
    )

    # Log the full TSV path for provenance
    logger.info(f"GOTTCHA2 profile TSV: {args.full_tsv}")

    # Build hit records
    hit_records = build_gottcha2_hit_records(
        parsed_entries=deduped_entries,
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
    # Skip in stateless mode (no state_dir)
    if context.state_dir is not None:
        try:
            mark_sample_completed(
                sample_id=context.sample_id,
                sample_set_id=context.sample_set_id,
                run_id=args.run_id,
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
