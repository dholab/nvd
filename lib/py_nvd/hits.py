# ruff: noqa: S608
"""
Hit management for metagenomic sequence identification.

Provides idempotent hit registration and querying. A "hit" is a sequence
identified by BLAST or GOTTCHA2. Hit keys are computed from canonical sequences
(lexicographically smaller of forward/reverse-complement), enabling:

- Cross-run deduplication: same sequence in different runs → same key
- Cross-sample deduplication: same sequence in different samples → same key
- Strand-agnostic identity: forward and reverse-complement → same key

Sequence Compression:
    Sequences are stored using 2-bit encoding (A=00, C=01, G=10, T=11) with
    zlib compression. Non-ACGT bases are coerced to N and their positions
    tracked separately for lossless recovery.

Storage (Parquet-based with Hive Partitioning):
    This module uses crash-proof Parquet files with Hive-style partitioning.
    Each sample's hits are stored in a separate parquet file, enabling:

    - Atomic writes: temp file + rename pattern prevents partial/corrupt files
    - Isolation: one sample's failure cannot affect another sample's data
    - Concurrent reads: multiple processes can query simultaneously
    - Efficient queries: DuckDB provides fast analytical queries across all files
    - Partition pruning: time-based queries skip irrelevant files

    Storage Layout (Hive-partitioned by month):
        Uncompacted (fresh from pipeline):
            {state_dir}/hits/month=NULL/{sample_set_id}/{sample_id}/data.parquet

        Compacted (after running `nvd hits compact`):
            {state_dir}/hits/month=2026-01/data.parquet

    The month=NULL partition holds uncompacted data. After compaction, data
    moves to month=YYYY-MM partitions sorted by hit_key for better compression
    and query performance.

    Write Pattern:
        1. Write to .data.parquet.tmp
        2. Atomic rename to data.parquet
        3. Readers never see partial files (glob ignores .tmp files)

    Query Pattern:
        DuckDB reads all parquet files via glob with hive_partitioning=true.
        The `month` column is extracted from the directory structure.
        The `first_seen_date` is computed at query time as MIN(run_date).
        The `effective_month` column normalizes month for both compacted and
        uncompacted data via COALESCE(month, strftime(run_date, '%Y-%m')).
"""

from __future__ import annotations

import fcntl
import re
import shutil
import zlib
from contextlib import contextmanager, suppress
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from collections.abc import Generator

import blake3
import duckdb
import polars as pl

from py_nvd.db import get_hits_dir
from py_nvd.models import (
    CategorySummary,
    CompactionResult,
    ContigQuality,
    DeleteResult,
    Hit,
    HitObservation,
    HitStats,
    MonthCompactionInfo,
    NovelTaxon,
    RareTaxon,
    RecurringHit,
    RunComparison,
    RunReport,
    SampleSummary,
    TaxonChange,
    TaxonHistory,
    TaxonRunPresence,
    TaxonSummary,
    TimelineBucket,
)

HIT_KEY_LENGTH = 32
BASES_PER_BYTE = 4
INT_BYTE_SIZE = 4
MONTH_FORMAT_LENGTH = 7  # "YYYY-MM" is 7 characters
SHORT_PATTERN_MAX_LENGTH = 4  # Word boundary threshold for regex patterns
MAX_PERCENTAGE = 100  # Upper bound for percentage values

COMPLEMENT = str.maketrans(
    "ACGTacgtNnRYSWKMBDHVryswkmbdhv",
    "TGCAtgcaNnYRSWMKVHDByrswmkvhdb",
)

# Canonical parquet schema for hit records. This is the target state after
# compaction normalizes legacy files. Column order matters for readability
# but not for correctness (DuckDB uses column names, not positions).
CANONICAL_HIT_SCHEMA: dict[str, pl.DataType] = {
    "hit_key": pl.Utf8,
    "sequence_length": pl.Int64,
    "sequence_compressed": pl.Binary,
    "gc_content": pl.Float64,
    "sample_set_id": pl.Utf8,
    "sample_id": pl.Utf8,
    "run_date": pl.Utf8,
    "sequence_id": pl.Utf8,
    "source": pl.Utf8,
    "adjusted_taxid": pl.Int64,
    "adjusted_taxid_name": pl.Utf8,
    "adjusted_taxid_rank": pl.Utf8,
}

# Column renames applied during normalization: old_name -> new_name.
SCHEMA_RENAMES: dict[str, str] = {
    "contig_id": "sequence_id",
}

# Default values for columns missing from legacy parquet files.
SCHEMA_DEFAULTS: dict[str, object] = {
    "source": "blast",
    "sequence_id": None,
}


def normalize_schema(df: pl.DataFrame) -> pl.DataFrame:
    """
    Normalize a DataFrame to the canonical hit schema.

    Applies column renames (e.g. contig_id -> sequence_id), adds missing
    columns with default values, and retains a NULL contig_id column for
    backward compatibility with the DuckDB view's COALESCE/EXCLUDE pattern.
    This is idempotent — normalizing an already-normalized DataFrame is a
    no-op.

    Used by compaction to ensure compacted output uses a consistent schema
    regardless of what schema the input files had.
    """
    for old_name, new_name in SCHEMA_RENAMES.items():
        if old_name in df.columns and new_name not in df.columns:
            df = df.rename({old_name: new_name})

    for col_name, default_value in SCHEMA_DEFAULTS.items():
        if col_name not in df.columns:
            df = df.with_columns(pl.lit(default_value).alias(col_name))

    # Retain contig_id as a NULL column so the DuckDB view's
    # EXCLUDE (contig_id, sequence_id) + COALESCE pattern always works,
    # even when reading only compacted files.
    if "contig_id" not in df.columns:
        df = df.with_columns(pl.lit(None).cast(pl.Utf8).alias("contig_id"))

    # Select canonical columns plus the transition contig_id stub.
    output_cols = [c for c in CANONICAL_HIT_SCHEMA if c in df.columns]
    if "contig_id" not in output_cols:
        output_cols.append("contig_id")
    return df.select(output_cols)


def is_valid_hit_key(key: str) -> bool:
    """
    Check if a string is a valid hit key (32 lowercase hex characters).

    Args:
        key: The string to validate

    Returns:
        True if key is exactly 32 lowercase hex characters, False otherwise.
    """
    return len(key) == HIT_KEY_LENGTH and all(c in "0123456789abcdef" for c in key)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def canonical_sequence(seq: str) -> str:
    """
    Return the canonical form of a sequence.

    The canonical form is the lexicographically smaller of the sequence
    and its reverse complement. This ensures strand-agnostic identity.
    """
    assert seq, "seq cannot be empty"
    assert isinstance(seq, str), f"seq must be str, got {type(seq)}"

    seq = seq.upper()
    revcomp = reverse_complement(seq)
    return min(seq, revcomp)


def compute_hit_key(seq: str) -> str:
    """
    Compute a deterministic hit key for a sequence.

    Returns a 32-character hex string (128 bits) from BLAKE3 hash
    of the canonical sequence.
    """
    assert seq, "seq cannot be empty"

    canonical = canonical_sequence(seq)
    digest = blake3.blake3(canonical.encode()).hexdigest()
    result = digest[:HIT_KEY_LENGTH]

    assert len(result) == HIT_KEY_LENGTH, (
        f"expected {HIT_KEY_LENGTH} chars, got {len(result)}"
    )
    return result


def compress_sequence(seq: str) -> bytes:
    """
    Compress a DNA sequence using 2-bit encoding with N position tracking.

    Encoding scheme:
    - A=00, C=01, G=10, T=11
    - Non-ACGT characters are coerced to N, stored as 00 (A)
    - N positions are appended as 4-byte little-endian integers
    - Final result is zlib compressed

    Format: zlib([2-bit-data][n_count:u32][n_pos_1:u32][n_pos_2:u32]...)
    """
    assert seq, "seq cannot be empty"
    assert isinstance(seq, str), f"seq must be str, got {type(seq)}"

    seq = seq.upper()

    n_positions = [i for i, base in enumerate(seq) if base not in "ACGT"]

    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    bits = [mapping.get(base, 0) for base in seq]

    result = bytearray()
    for i in range(0, len(bits), BASES_PER_BYTE):
        chunk = bits[i : i + BASES_PER_BYTE]
        while len(chunk) < BASES_PER_BYTE:
            chunk.append(0)
        byte = (chunk[0] << 6) | (chunk[1] << 4) | (chunk[2] << 2) | chunk[3]
        result.append(byte)

    result.extend(len(n_positions).to_bytes(INT_BYTE_SIZE, "little"))
    for pos in n_positions:
        result.extend(pos.to_bytes(INT_BYTE_SIZE, "little"))

    compressed = zlib.compress(bytes(result), level=9)
    assert len(compressed) > 0, "compression should produce non-empty result"
    return compressed


def decompress_sequence(data: bytes, length: int) -> str:
    """
    Decompress a sequence from 2-bit encoding with N position tracking.

    Args:
        data: Compressed sequence blob from compress_sequence()
        length: Original sequence length in bases

    Returns:
        Original DNA sequence
    """
    assert data, "data cannot be empty"
    assert length > 0, f"length must be positive, got {length}"

    raw = zlib.decompress(data)

    two_bit_len = (length + BASES_PER_BYTE - 1) // BASES_PER_BYTE
    two_bit_data = raw[:two_bit_len]

    n_count = int.from_bytes(raw[two_bit_len : two_bit_len + INT_BYTE_SIZE], "little")
    n_positions = []
    offset = two_bit_len + INT_BYTE_SIZE
    for _ in range(n_count):
        n_positions.append(
            int.from_bytes(raw[offset : offset + INT_BYTE_SIZE], "little"),
        )
        offset += INT_BYTE_SIZE

    mapping = ["A", "C", "G", "T"]
    bases = []
    for byte in two_bit_data:
        bases.append(mapping[(byte >> 6) & 0x3])
        bases.append(mapping[(byte >> 4) & 0x3])
        bases.append(mapping[(byte >> 2) & 0x3])
        bases.append(mapping[byte & 0x3])
    bases = bases[:length]

    for pos in n_positions:
        bases[pos] = "N"

    result = "".join(bases)
    assert len(result) == length, f"expected {length} bases, got {len(result)}"
    return result


def calculate_gc_content(seq: str) -> float:
    """Calculate GC content as a fraction (0.0 to 1.0)."""
    assert isinstance(seq, str), f"seq must be str, got {type(seq)}"

    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    total = len(seq)
    result = gc_count / total if total > 0 else 0.0

    assert 0.0 <= result <= 1.0, f"GC content must be 0.0-1.0, got {result}"
    return result


def _hits_glob_pattern(state_dir: Path | str | None = None) -> str:
    """
    Get the glob pattern for reading all parquet hit files.

    Returns a pattern like: /path/to/state/hits/**/*.parquet
    """
    hits_dir = get_hits_dir(state_dir)
    return str(hits_dir / "**" / "*.parquet")


def _get_hits_view_sql(glob_pattern: str) -> str:
    """
    Generate SQL to create a view over all parquet files.

    Uses hive_partitioning=true to extract `month` from directory structure.
    Uses union_by_name=true for schema evolution compatibility.

    The view adds an `effective_month` column that:
    - Returns the Hive-extracted `month` for compacted files (e.g., month=2026-01)
    - Computes month from `run_date` for uncompacted files (month=NULL)
    """
    assert glob_pattern, "glob_pattern cannot be empty"
    assert glob_pattern.endswith(".parquet"), "glob_pattern must end with .parquet"

    return f"""
        CREATE OR REPLACE VIEW hits AS
        SELECT
            * EXCLUDE (contig_id, sequence_id),
            COALESCE(sequence_id, contig_id) AS sequence_id,
            COALESCE(month::VARCHAR, strftime(run_date::DATE, '%Y-%m')) AS effective_month
        FROM read_parquet('{glob_pattern}', hive_partitioning=true, union_by_name=true)
    """  # SQL uses module constants, not user input


def query_hits(state_dir: Path | str | None = None) -> duckdb.DuckDBPyConnection:
    """
    Get a DuckDB connection configured to query hits.

    Creates an in-memory DuckDB connection with a 'hits' view that reads
    from all parquet files in the hits directory. The view is recreated
    on each call to pick up any new files.

    Args:
        state_dir: Optional state directory override

    Returns:
        DuckDB connection with 'hits' view available for querying

    Example:
        con = query_hits()
        result = con.execute("SELECT COUNT(DISTINCT hit_key) FROM hits").fetchone()
    """
    glob_pattern = _hits_glob_pattern(state_dir)
    con = duckdb.connect()  # In-memory, no persistence needed

    # Only create the view if there are parquet files to read
    # Otherwise DuckDB will error on empty glob
    hits_dir = get_hits_dir(state_dir)
    if any(hits_dir.glob("**/*.parquet")):
        con.execute(_get_hits_view_sql(glob_pattern))

    return con


def has_any_hits(state_dir: Path | str | None = None) -> bool:
    """
    Check if any parquet hit files exist.

    Useful for early-exit in query functions when the database is empty.
    """
    hits_dir = get_hits_dir(state_dir)
    return any(hits_dir.glob("**/*.parquet"))


@dataclass(frozen=True)
class HitRecord:
    """
    A single hit record ready for parquet serialization.

    This is the denormalized form that combines hit data with observation metadata.
    Each row in the parquet file represents one observation of a hit in a sample.

    Classification fields (adjusted_taxid, adjusted_taxid_name, adjusted_taxid_rank)
    contain the consensus taxonomic assignment. For BLAST hits this is the LCA-adjusted
    classification; for GOTTCHA2 hits this is the direct taxonomic assignment from
    profiling. These may be None for hits from older runs that predate classification
    tracking.
    """

    hit_key: str
    sequence_length: int
    sequence_compressed: bytes
    gc_content: float
    sample_set_id: str
    sample_id: str
    run_date: str
    sequence_id: str | None
    source: Literal["blast", "gottcha2"]
    adjusted_taxid: int | None = None
    adjusted_taxid_name: str | None = None
    adjusted_taxid_rank: str | None = None

    def __post_init__(self) -> None:
        """Validate invariants after initialization."""
        assert len(self.hit_key) == HIT_KEY_LENGTH, (
            f"hit_key must be {HIT_KEY_LENGTH} hex chars, got {len(self.hit_key)}"
        )
        assert all(c in "0123456789abcdef" for c in self.hit_key), (
            "hit_key must be lowercase hex"
        )
        assert self.sequence_length > 0, (
            f"sequence_length must be positive, got {self.sequence_length}"
        )
        assert len(self.sequence_compressed) > 0, "sequence_compressed cannot be empty"
        assert 0.0 <= self.gc_content <= 1.0, (
            f"gc_content must be 0.0-1.0, got {self.gc_content}"
        )
        assert self.sample_set_id, "sample_set_id cannot be empty"
        assert self.sample_id, "sample_id cannot be empty"
        assert self.run_date, "run_date cannot be empty"
        assert self.source in ("blast", "gottcha2"), (
            f"source must be 'blast' or 'gottcha2', got {self.source!r}"
        )


def write_hits_parquet(
    hits: list[HitRecord],
    sample_id: str,
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> Path:
    """
    Write hits for a sample to a parquet file atomically.

    Uses the temp file + rename pattern to ensure readers never see partial files.
    If a file already exists for this sample, it is overwritten (last writer wins).

    Files are written to a Hive-partitioned structure:
        hits/month=NULL/{sample_set_id}/{sample_id}/data.parquet

    The month=NULL partition allows DuckDB to read uncompacted and compacted files
    together with hive_partitioning=true. After compaction, data moves to
    month={YYYY-MM}/ partitions.

    Args:
        hits: List of HitRecord objects to write
        sample_id: The sample identifier (used for subdirectory)
        sample_set_id: The sample set identifier (used for subdirectory)
        state_dir: Optional state directory override

    Returns:
        Path to the written parquet file

    Raises:
        OSError: If the atomic rename fails
    """
    assert sample_id, "sample_id cannot be empty"
    assert sample_set_id, "sample_set_id cannot be empty"

    hits_dir = get_hits_dir(state_dir)
    sample_dir = hits_dir / "month=NULL" / sample_set_id / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    final_path = sample_dir / "data.parquet"
    temp_path = sample_dir / ".data.parquet.tmp"

    if not hits:
        # Write empty parquet with correct schema.
        # Includes contig_id=NULL for backward compatibility with union_by_name
        # reads against legacy parquet files (see COALESCE strategy in plan).
        hits_df = pl.DataFrame(
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
        # Build DataFrame from hit records.
        # Writes contig_id=NULL alongside sequence_id so that
        # union_by_name=true always surfaces both columns when reading
        # a mix of old and new parquet files.
        n = len(hits)
        hits_df = pl.DataFrame(
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

    # Write to temp file first
    hits_df.write_parquet(temp_path)

    # Atomic rename (POSIX guarantees this is atomic on same filesystem)
    temp_path.rename(final_path)

    return final_path


def get_sample_parquet_path(
    sample_id: str,
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> Path:
    """
    Get the path where a sample's uncompacted parquet file would be stored.

    Returns the path in the month=NULL partition:
        hits/month=NULL/{sample_set_id}/{sample_id}/data.parquet

    Does not check if the file exists.
    """
    assert sample_id, "sample_id cannot be empty"
    assert sample_set_id, "sample_set_id cannot be empty"

    hits_dir = get_hits_dir(state_dir)
    return hits_dir / "month=NULL" / sample_set_id / sample_id / "data.parquet"


def sample_hits_exist(
    sample_id: str,
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> bool:
    """
    Check if a parquet file exists for a specific sample.

    Useful for idempotency checks before writing.
    """
    assert sample_id, "sample_id cannot be empty"
    assert sample_set_id, "sample_set_id cannot be empty"

    return get_sample_parquet_path(sample_id, sample_set_id, state_dir).exists()


def count_hits(state_dir: Path | str | None = None) -> int:
    """Return the total number of unique hits across all parquet files."""
    if not has_any_hits(state_dir):
        return 0

    con = query_hits(state_dir)
    result = con.execute("SELECT COUNT(DISTINCT hit_key) FROM hits").fetchone()
    assert result is not None, "COUNT query should always return a row"
    count = result[0]
    assert isinstance(count, int), f"count must be int, got {type(count)}"
    assert count >= 0, f"count must be non-negative, got {count}"
    return count


def count_hit_observations(state_dir: Path | str | None = None) -> int:
    """Return the total number of hit observations across all parquet files."""
    if not has_any_hits(state_dir):
        return 0

    con = query_hits(state_dir)
    result = con.execute("SELECT COUNT(*) FROM hits").fetchone()
    assert result is not None, "COUNT query should always return a row"
    count = result[0]
    assert isinstance(count, int), f"count must be int, got {type(count)}"
    assert count >= 0, f"count must be non-negative, got {count}"
    return count


def get_hit(hit_key: str, state_dir: Path | str | None = None) -> Hit | None:
    """
    Get a hit by its key.

    Since parquet stores denormalized data, this aggregates across all
    observations to compute first_seen_date as MIN(run_date). Classification
    is taken from the most recent observation with non-null adjusted_taxid.
    """
    assert hit_key, "hit_key cannot be empty"
    assert len(hit_key) == HIT_KEY_LENGTH, (
        f"hit_key must be {HIT_KEY_LENGTH} hex chars, got {len(hit_key)}"
    )

    if not has_any_hits(state_dir):
        return None

    con = query_hits(state_dir)
    # Get basic hit info with first_seen_date
    # Then join with most recent classification (where adjusted_taxid is not null)
    result = con.execute(
        """
        WITH hit_base AS (
            SELECT
                hit_key,
                ANY_VALUE(sequence_length) as sequence_length,
                ANY_VALUE(sequence_compressed) as sequence_compressed,
                ANY_VALUE(gc_content) as gc_content,
                MIN(run_date) as first_seen_date
            FROM hits
            WHERE hit_key = ?
            GROUP BY hit_key
        ),
        most_recent_classification AS (
            SELECT
                hit_key,
                adjusted_taxid,
                adjusted_taxid_name,
                adjusted_taxid_rank
            FROM hits
            WHERE hit_key = ? AND adjusted_taxid IS NOT NULL
            ORDER BY run_date DESC
            LIMIT 1
        )
        SELECT
            h.hit_key,
            h.sequence_length,
            h.sequence_compressed,
            h.gc_content,
            h.first_seen_date,
            c.adjusted_taxid,
            c.adjusted_taxid_name,
            c.adjusted_taxid_rank
        FROM hit_base h
        LEFT JOIN most_recent_classification c ON h.hit_key = c.hit_key
        """,
        [hit_key.lower(), hit_key.lower()],
    ).fetchone()

    if result is None:
        return None

    return Hit(
        hit_key=result[0],
        sequence_length=result[1],
        sequence_compressed=result[2],
        gc_content=result[3],
        first_seen_date=result[4],
        top_taxid=result[5],
        top_taxid_name=result[6],
        top_taxid_rank=result[7],
    )


def list_hits(
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[Hit]:
    """
    List all unique hits, ordered by first_seen_date descending.

    Classification is taken from the most recent observation with non-null
    adjusted_taxid for each hit.

    Args:
        limit: Maximum number of hits to return
        state_dir: Optional state directory override

    Returns:
        List of Hit objects
    """
    assert limit is None or limit > 0, "limit must be positive if provided"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Get basic hit info, then join with most recent classification per hit
    sql = """
        WITH hit_base AS (
            SELECT
                hit_key,
                ANY_VALUE(sequence_length) as sequence_length,
                ANY_VALUE(sequence_compressed) as sequence_compressed,
                ANY_VALUE(gc_content) as gc_content,
                MIN(run_date) as first_seen_date
            FROM hits
            GROUP BY hit_key
        ),
        most_recent_classification AS (
            SELECT DISTINCT ON (hit_key)
                hit_key,
                adjusted_taxid,
                adjusted_taxid_name,
                adjusted_taxid_rank
            FROM hits
            WHERE adjusted_taxid IS NOT NULL
            ORDER BY hit_key, run_date DESC
        )
        SELECT
            h.hit_key,
            h.sequence_length,
            h.sequence_compressed,
            h.gc_content,
            h.first_seen_date,
            c.adjusted_taxid,
            c.adjusted_taxid_name,
            c.adjusted_taxid_rank
        FROM hit_base h
        LEFT JOIN most_recent_classification c ON h.hit_key = c.hit_key
        ORDER BY h.first_seen_date DESC
    """

    if limit is not None:
        # f-string is safe here because limit is validated as positive int above
        sql += f" LIMIT {limit}"

    rows = con.execute(sql).fetchall()

    return [
        Hit(
            hit_key=row[0],
            sequence_length=row[1],
            sequence_compressed=row[2],
            gc_content=row[3],
            first_seen_date=row[4],
            top_taxid=row[5],
            top_taxid_name=row[6],
            top_taxid_rank=row[7],
        )
        for row in rows
    ]


def list_hit_observations(
    hit_key: str | None = None,
    sample_id: str | None = None,
    sample_set_id: str | None = None,
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[HitObservation]:
    """
    List hit observations with optional filters.

    Args:
        hit_key: Filter to observations of a specific hit
        sample_id: Filter to observations from a specific sample
        sample_set_id: Filter to observations from a specific sample set
        limit: Maximum number of observations to return
        state_dir: Optional state directory override

    Returns:
        List of HitObservation objects, ordered by run_date descending
    """
    assert hit_key is None or hit_key, (
        "hit_key cannot be empty string (use None for no filter)"
    )
    assert sample_id is None or sample_id, (
        "sample_id cannot be empty string (use None for no filter)"
    )
    assert sample_set_id is None or sample_set_id, (
        "sample_set_id cannot be empty string"
    )
    assert limit is None or limit > 0, "limit must be positive if provided"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    sql = "SELECT hit_key, sample_set_id, sample_id, run_date, sequence_id, source FROM hits"
    conditions: list[str] = []
    params: list[str] = []

    if hit_key is not None:
        conditions.append("hit_key = ?")
        params.append(hit_key.lower())

    if sample_id is not None:
        conditions.append("sample_id = ?")
        params.append(sample_id)

    if sample_set_id is not None:
        conditions.append("sample_set_id = ?")
        params.append(sample_set_id)

    if conditions:
        sql += " WHERE " + " AND ".join(conditions)

    sql += " ORDER BY run_date DESC"

    if limit is not None:
        # f-string is safe here because limit is validated as positive int above
        sql += f" LIMIT {limit}"

    rows = con.execute(sql, params).fetchall()

    return [
        HitObservation(
            hit_key=row[0],
            sample_set_id=row[1],
            sample_id=row[2],
            run_date=row[3],
            sequence_id=row[4],
            source=row[5],
        )
        for row in rows
    ]


def lookup_hit(
    query: str,
    state_dir: Path | str | None = None,
) -> Hit | None:
    """
    Find a hit by key or sequence.

    Auto-detects whether the query is a hit key (32 hex characters) or a
    DNA sequence. For sequences, computes the canonical hash and looks up.

    Args:
        query: Either a 32-character hex hit key, or a DNA sequence
        state_dir: Optional state directory override

    Returns:
        The Hit if found, None otherwise.
    """
    assert query, "query cannot be empty"
    assert isinstance(query, str), f"query must be str, got {type(query)}"

    query = query.strip()

    # Check if it looks like a hex key (32 chars, all hex digits)
    if is_valid_hit_key(query.lower()):
        return get_hit(query.lower(), state_dir)

    # Otherwise treat as sequence, compute key and look up
    hit_key = compute_hit_key(query)
    return get_hit(hit_key, state_dir)


def get_hit_sequence(hit: Hit) -> str:
    """
    Decompress and return the sequence for a hit.

    Args:
        hit: The Hit object (must have sequence_compressed and sequence_length)

    Returns:
        The decompressed DNA sequence
    """
    assert hit is not None, "hit cannot be None"
    assert hit.sequence_compressed, "hit must have sequence_compressed"
    assert hit.sequence_length > 0, "hit must have positive sequence_length"

    return decompress_sequence(hit.sequence_compressed, hit.sequence_length)


@dataclass(frozen=True)
class SampleSetStats:
    """
    Summary statistics for hits within a specific sample set (run).

    Used for run-specific reporting (e.g., Slack notifications) as opposed
    to cumulative stats across all runs.
    """

    sample_set_id: str
    total_observations: int  # Number of hit observations in this sample set
    unique_hits: int  # Distinct hit_keys in this sample set
    unique_samples: int  # Distinct sample_ids in this sample set
    date_first: str | None  # Earliest run_date in this sample set
    date_last: str | None  # Latest run_date in this sample set


def _empty_hit_stats() -> HitStats:
    """Return HitStats with zero counts and None distributions."""
    return HitStats(
        total_hits=0,
        total_observations=0,
        unique_samples=0,
        unique_runs=0,
        length_min=None,
        length_max=None,
        length_median=None,
        gc_min=None,
        gc_max=None,
        gc_median=None,
        date_first=None,
        date_last=None,
    )


def get_stats_for_sample_set(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> SampleSetStats | None:
    """
    Compute statistics for hits within a specific sample set (run).

    Returns run-specific metrics useful for notifications and reporting.
    Unlike get_hit_stats() which aggregates across all data, this focuses
    on a single sample_set_id.

    Args:
        sample_set_id: The sample set identifier to filter by
        state_dir: Optional state directory override

    Returns:
        SampleSetStats for the sample set, or None if no hits exist for it.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return None

    con = query_hits(state_dir)

    # Check if this sample set has any data
    count_result = con.execute(
        "SELECT COUNT(*) FROM hits WHERE sample_set_id = ?",
        [sample_set_id],
    ).fetchone()
    assert count_result is not None, "COUNT query should always return a row"

    if count_result[0] == 0:
        return None

    # Get all stats in one query
    stats_result = con.execute(
        """
        SELECT
            COUNT(*) as total_observations,
            COUNT(DISTINCT hit_key) as unique_hits,
            COUNT(DISTINCT sample_id) as unique_samples,
            MIN(run_date) as date_first,
            MAX(run_date) as date_last
        FROM hits
        WHERE sample_set_id = ?
        """,
        [sample_set_id],
    ).fetchone()
    assert stats_result is not None, "aggregate query should always return a row"

    return SampleSetStats(
        sample_set_id=sample_set_id,
        total_observations=stats_result[0],
        unique_hits=stats_result[1],
        unique_samples=stats_result[2],
        date_first=stats_result[3],
        date_last=stats_result[4],
    )


def get_hit_stats(state_dir: Path | str | None = None) -> HitStats:
    """
    Compute summary statistics for all hits.

    Returns aggregate metrics including counts, length/GC distributions,
    and date range. Useful for dashboard displays and health checks.

    Since parquet stores denormalized data, first_seen_date is computed
    as MIN(run_date) for each unique hit.

    Args:
        state_dir: Optional state directory override

    Returns:
        HitStats with all computed metrics. Counts are 0 and distributions
        are None when there are no hits.
    """
    if not has_any_hits(state_dir):
        return _empty_hit_stats()

    con = query_hits(state_dir)

    # Basic counts - use explicit assertions before subscripting
    hits_result = con.execute("SELECT COUNT(DISTINCT hit_key) FROM hits").fetchone()
    assert hits_result is not None, "COUNT query should always return a row"
    total_hits = hits_result[0]

    obs_result = con.execute("SELECT COUNT(*) FROM hits").fetchone()
    assert obs_result is not None, "COUNT query should always return a row"
    total_observations = obs_result[0]

    samples_result = con.execute(
        "SELECT COUNT(DISTINCT sample_id) FROM hits",
    ).fetchone()
    assert samples_result is not None, "COUNT query should always return a row"
    unique_samples = samples_result[0]

    runs_result = con.execute(
        "SELECT COUNT(DISTINCT sample_set_id) FROM hits",
    ).fetchone()
    assert runs_result is not None, "COUNT query should always return a row"
    unique_runs = runs_result[0]

    assert isinstance(total_hits, int), "total_hits must be int"
    assert total_hits >= 0, "total_hits must be non-negative"
    assert isinstance(total_observations, int), "total_observations must be int"
    assert total_observations >= 0, "total_observations must be non-negative"

    # If no hits, return early with empty stats
    if total_hits == 0:
        return _empty_hit_stats()

    # Length stats (min, max) - need to aggregate by hit_key first to get unique hits
    length_row = con.execute("""
        SELECT MIN(sequence_length), MAX(sequence_length)
        FROM (SELECT DISTINCT hit_key, sequence_length FROM hits)
    """).fetchone()
    assert length_row is not None, "MIN/MAX query should always return a row"
    length_min, length_max = length_row[0], length_row[1]

    # GC stats (min, max)
    gc_row = con.execute("""
        SELECT MIN(gc_content), MAX(gc_content)
        FROM (SELECT DISTINCT hit_key, gc_content FROM hits)
    """).fetchone()
    assert gc_row is not None, "MIN/MAX query should always return a row"
    gc_min, gc_max = gc_row[0], gc_row[1]

    # Date range (first_seen_date = MIN(run_date) per hit)
    date_row = con.execute("""
        SELECT MIN(first_seen), MAX(first_seen)
        FROM (
            SELECT hit_key, MIN(run_date) as first_seen
            FROM hits
            GROUP BY hit_key
        )
    """).fetchone()
    assert date_row is not None, "MIN/MAX query should always return a row"
    date_first, date_last = date_row[0], date_row[1]

    # Median calculations - DuckDB has MEDIAN function!
    length_median_row = con.execute("""
        SELECT MEDIAN(sequence_length)
        FROM (SELECT DISTINCT hit_key, sequence_length FROM hits)
    """).fetchone()
    assert length_median_row is not None, "MEDIAN query should always return a row"
    length_median = (
        float(length_median_row[0]) if length_median_row[0] is not None else None
    )

    gc_median_row = con.execute("""
        SELECT MEDIAN(gc_content)
        FROM (SELECT DISTINCT hit_key, gc_content FROM hits)
    """).fetchone()
    assert gc_median_row is not None, "MEDIAN query should always return a row"
    gc_median = float(gc_median_row[0]) if gc_median_row[0] is not None else None

    return HitStats(
        total_hits=total_hits,
        total_observations=total_observations,
        unique_samples=unique_samples,
        unique_runs=unique_runs,
        length_min=length_min,
        length_max=length_max,
        length_median=length_median,
        gc_min=gc_min,
        gc_max=gc_max,
        gc_median=gc_median,
        date_first=date_first,
        date_last=date_last,
    )


# Taxa that are filtered out by default in get_top_taxa() because they
# represent unclassified hits or host contamination, not viral findings.
DEFAULT_NOISE_TAXA = frozenset({"root", "Homo sapiens"})

# Mapping from virus family to search patterns for adjusted_taxid_name.
# Used by get_taxa_by_category() to group findings by virus family for triage.
#
# Patterns are matched with word boundaries (e.g., "hiv" matches "HIV-1" but
# not "archivirus"). Matching is case-insensitive. First matching family wins.
#
# Based on human_virus_families from nextflow.config (the families we enrich
# for during pipeline processing). Update this mapping when adding new
# families to the enrichment list.
#
# This is a starting point—patterns may need refinement as we encounter
# taxa with unexpected naming conventions in NCBI taxonomy.
FAMILY_SEARCH_PATTERNS: dict[str, tuple[str, ...]] = {
    "Adenoviridae": ("adenovirus",),
    "Anelloviridae": ("anellovirus", "torque teno", "ttv"),
    "Arenaviridae": ("arenavirus", "lassa", "junin", "machupo"),
    "Arteriviridae": ("arterivirus",),
    "Astroviridae": ("astrovirus",),
    "Bornaviridae": ("bornavirus",),
    "Peribunyaviridae": ("bunyavirus", "orthobunyavirus", "la crosse"),
    "Caliciviridae": ("calicivirus", "norovirus", "sapovirus"),
    "Coronaviridae": ("coronavirus", "sars-cov", "mers-cov"),
    "Filoviridae": ("filovirus", "ebola", "marburg"),
    "Flaviviridae": (
        "flavivirus",
        "dengue",
        "zika",
        "hepatitis c",
        "yellow fever",
        "west nile",
    ),
    "Hepadnaviridae": ("hepadnavirus", "hepatitis b"),
    "Hepeviridae": ("hepevirus", "hepatitis e"),
    "Orthoherpesviridae": (
        "herpesvirus",
        "herpes simplex",
        "varicella",
        "epstein-barr",
        "cytomegalovirus",
    ),
    "Orthomyxoviridae": ("influenza", "orthomyxovirus"),
    "Papillomaviridae": ("papillomavirus", "hpv"),
    "Paramyxoviridae": (
        "paramyxovirus",
        "measles",
        "mumps",
        "parainfluenza",
        "henipavirus",
        "nipah",
    ),
    "Parvoviridae": ("parvovirus", "bocavirus", "erythrovirus"),
    "Picobirnaviridae": ("picobirnavirus",),
    "Picornaviridae": (
        "picornavirus",
        "rhinovirus",
        "enterovirus",
        "coxsackie",
        "poliovirus",
        "hepatitis a",
        "parechovirus",
    ),
    "Pneumoviridae": ("pneumovirus", "respiratory syncytial", "rsv", "metapneumovirus"),
    "Polyomaviridae": ("polyomavirus", "merkel cell", "jc virus", "bk virus"),
    "Poxviridae": ("poxvirus", "variola", "monkeypox", "mpox", "vaccinia", "molluscum"),
    "Sedoreoviridae": ("rotavirus", "reovirus", "orbivirus"),
    "Retroviridae": ("retrovirus", "hiv", "htlv", "lentivirus", "t-lymphotropic"),
    "Rhabdoviridae": ("rhabdovirus", "rabies", "lyssavirus", "vesiculovirus"),
    "Togaviridae": ("togavirus", "alphavirus", "rubivirus", "rubella", "chikungunya"),
    "Kolmioviridae": ("deltavirus", "hepatitis d", "hepatitis delta"),
}


# Compiled regex patterns for family classification.
# Short patterns (<=4 chars like "hiv", "rsv") use word boundaries to prevent
# false positives (e.g., "hiv" should not match "archivirus").
# Longer patterns use substring matching to handle compound names
# (e.g., "herpesvirus" should match "gammaherpesvirus").
# Compiled once at module load for efficiency.
def _build_family_pattern(patterns: tuple[str, ...]) -> re.Pattern[str]:
    """Build a regex pattern for a family's search patterns."""
    parts = []
    for p in patterns:
        escaped = re.escape(p)
        if len(p) <= SHORT_PATTERN_MAX_LENGTH:
            # Short patterns need word boundaries to avoid false positives
            parts.append(r"\b" + escaped + r"\b")
        else:
            # Longer patterns can use substring matching
            parts.append(escaped)
    return re.compile("(" + "|".join(parts) + ")", re.IGNORECASE)


_FAMILY_PATTERNS_COMPILED: list[tuple[str, re.Pattern[str]]] = [
    (family, _build_family_pattern(patterns))
    for family, patterns in FAMILY_SEARCH_PATTERNS.items()
]


def get_top_taxa(
    sample_set_id: str | None = None,
    limit: int = 10,
    exclude_noise: bool = True,  # noqa: FBT001, FBT002  # existing public API
    state_dir: Path | str | None = None,
) -> list[TaxonSummary]:
    """
    Get the most prevalent taxa by sample count.

    Returns taxa ordered by the number of samples they appear in (descending),
    then by hit count. This ordering prioritizes breadth of detection over
    depth, which is more meaningful for interpreting metagenomic findings.

    Note: Results require expert interpretation. Common false positives include:
    - Endogenous retroviruses (ERVs) matching HIV/HTLV
    - Host sequences matching viral databases at conserved regions
    - Database artifacts and chimeric reference sequences

    Args:
        sample_set_id: Filter to a specific run. If None, queries all runs.
        limit: Maximum number of taxa to return (default: 10)
        exclude_noise: If True, excludes 'root' and 'Homo sapiens' (default: True)
        state_dir: Optional state directory override

    Returns:
        List of TaxonSummary ordered by sample_count DESC, hit_count DESC.
        Empty list if no hits exist or no taxa match the filters.
    """
    assert limit > 0, "limit must be positive"
    assert sample_set_id is None or sample_set_id, (
        "sample_set_id cannot be empty string (use None for all runs)"
    )

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Build WHERE clause
    conditions: list[str] = ["adjusted_taxid_name IS NOT NULL"]
    params: list[str | int] = []

    if sample_set_id is not None:
        conditions.append("sample_set_id = ?")
        params.append(sample_set_id)

    if exclude_noise:
        # Use placeholders for the noise taxa to avoid SQL injection
        placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)
        conditions.append(f"adjusted_taxid_name NOT IN ({placeholders})")
        params.extend(DEFAULT_NOISE_TAXA)

    where_clause = " AND ".join(conditions)

    sql = f"""
        SELECT
            adjusted_taxid as taxid,
            adjusted_taxid_name as name,
            adjusted_taxid_rank as rank,
            COUNT(DISTINCT sample_id) as sample_count,
            COUNT(DISTINCT hit_key) as hit_count,
            ROUND(AVG(sequence_length))::INTEGER as avg_contig_length
        FROM hits
        WHERE {where_clause}
        GROUP BY adjusted_taxid, adjusted_taxid_name, adjusted_taxid_rank
        ORDER BY sample_count DESC, hit_count DESC
        LIMIT {limit}
    """  # SQL uses module constants, not user input

    rows = con.execute(sql, params).fetchall()

    return [
        TaxonSummary(
            taxid=row[0],
            name=row[1],
            rank=row[2],
            sample_count=row[3],
            hit_count=row[4],
            avg_contig_length=row[5] if row[5] is not None else 0,
        )
        for row in rows
    ]


def get_highlights_string(
    sample_set_id: str | None = None,
    limit: int = 5,
    state_dir: Path | str | None = None,
) -> str:
    """
    Get a one-liner summary of top taxa for display (e.g., Slack messages).

    Returns a formatted string like:
    "EBV (12 samples), Rhinovirus A (8), HPV (7), JC polyomavirus (5)"

    Uses get_top_taxa() internally, so the same filtering and ordering applies.

    Args:
        sample_set_id: Filter to a specific run. If None, queries all runs.
        limit: Maximum number of taxa to include (default: 5)
        state_dir: Optional state directory override

    Returns:
        Formatted string of top taxa with sample counts.
        Empty string if no taxa found.
    """
    taxa = get_top_taxa(
        sample_set_id=sample_set_id,
        limit=limit,
        exclude_noise=True,
        state_dir=state_dir,
    )

    if not taxa:
        return ""

    return ", ".join(f"{t.name} ({t.sample_count})" for t in taxa)


def get_sample_summaries(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> list[SampleSummary]:
    """
    Get per-sample breakdown of hits for a run.

    Returns a summary for each sample showing total hits and taxa detected.
    Ordered by viral_taxa_count descending (most interesting samples first),
    then by total_hits descending.

    Noise taxa ('root', 'Homo sapiens') are excluded from viral_taxa_count
    and taxa_detected.

    Args:
        sample_set_id: The run to query (required)
        state_dir: Optional state directory override

    Returns:
        List of SampleSummary ordered by viral_taxa_count DESC, total_hits DESC.
        Empty list if no hits exist for the sample set.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Build noise filter placeholders
    # Note: We cast adjusted_taxid_name to VARCHAR to handle edge cases where
    # parquet schema inference produces INTEGER when all values are NULL.
    noise_placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)

    sql = f"""
        SELECT
            sample_id,
            COUNT(DISTINCT hit_key) as total_hits,
            COUNT(DISTINCT adjusted_taxid_name) FILTER (
                WHERE adjusted_taxid_name IS NOT NULL
                  AND CAST(adjusted_taxid_name AS VARCHAR) NOT IN ({noise_placeholders})
            ) as viral_taxa_count,
            list(DISTINCT adjusted_taxid_name ORDER BY adjusted_taxid_name) FILTER (
                WHERE adjusted_taxid_name IS NOT NULL
                  AND CAST(adjusted_taxid_name AS VARCHAR) NOT IN ({noise_placeholders})
            ) as taxa_detected
        FROM hits
        WHERE sample_set_id = ?
        GROUP BY sample_id
        ORDER BY viral_taxa_count DESC, total_hits DESC
    """  # SQL uses module constants, not user input

    # Parameters: noise taxa twice (for COUNT and list filters), then sample_set_id
    params = [*DEFAULT_NOISE_TAXA, *DEFAULT_NOISE_TAXA, sample_set_id]

    rows = con.execute(sql, params).fetchall()

    return [
        SampleSummary(
            sample_id=row[0],
            total_hits=row[1],
            viral_taxa_count=row[2],
            taxa_detected=tuple(row[3]) if row[3] else (),
        )
        for row in rows
    ]


def get_negative_samples(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> list[str]:
    """
    Get samples with no viral hits (only noise taxa like 'root' or 'Homo sapiens').

    Useful for identifying samples that may have failed processing, had
    insufficient sequencing depth, or genuinely have no detectable viral content.

    Args:
        sample_set_id: The run to query (required)
        state_dir: Optional state directory override

    Returns:
        List of sample_id strings for samples with no viral taxa.
        Sorted alphabetically. Empty list if all samples have viral hits.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Convert frozenset to list for consistent parameter ordering
    noise_list = list(DEFAULT_NOISE_TAXA)
    noise_placeholders = ", ".join("?" for _ in noise_list)

    # Note: We cast adjusted_taxid_name to VARCHAR to handle edge cases where
    # parquet schema inference produces INTEGER when all values are NULL.
    sql = f"""
        SELECT sample_id
        FROM hits
        WHERE sample_set_id = ?
        GROUP BY sample_id
        HAVING COUNT(DISTINCT adjusted_taxid_name) FILTER (
            WHERE adjusted_taxid_name IS NOT NULL
              AND CAST(adjusted_taxid_name AS VARCHAR) NOT IN ({noise_placeholders})
        ) = 0
        ORDER BY sample_id
    """  # SQL uses module constants, not user input

    params = [sample_set_id, *noise_list]

    rows = con.execute(sql, params).fetchall()

    return [row[0] for row in rows]


def get_contig_quality(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> ContigQuality | None:
    """
    Get assembly quality metrics for contigs in a run.

    Computes statistics about contig lengths and GC content, useful for
    QC and troubleshooting. Each unique hit_key is counted once (not per
    observation).

    Args:
        sample_set_id: The run to query (required)
        state_dir: Optional state directory override

    Returns:
        ContigQuality with computed metrics, or None if no hits exist
        for the sample set.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return None

    con = query_hits(state_dir)

    sql = """
        WITH unique_contigs AS (
            SELECT DISTINCT
                hit_key,
                sequence_length,
                gc_content
            FROM hits
            WHERE sample_set_id = ?
        )
        SELECT
            COUNT(*) as total_contigs,
            COALESCE(SUM(sequence_length), 0) as total_bases,
            COALESCE(AVG(sequence_length), 0) as mean_length,
            COALESCE(MEDIAN(sequence_length), 0)::INTEGER as median_length,
            COALESCE(MAX(sequence_length), 0) as max_length,
            COUNT(*) FILTER (WHERE sequence_length >= 500) as contigs_500bp_plus,
            COUNT(*) FILTER (WHERE sequence_length >= 1000) as contigs_1kb_plus,
            COALESCE(AVG(gc_content), 0) as mean_gc
        FROM unique_contigs
    """

    row = con.execute(sql, [sample_set_id]).fetchone()

    if row is None or row[0] == 0:
        return None

    return ContigQuality(
        total_contigs=row[0],
        total_bases=row[1],
        mean_length=row[2],
        median_length=row[3],
        max_length=row[4],
        contigs_500bp_plus=row[5],
        contigs_1kb_plus=row[6],
        mean_gc=row[7],
    )


def _classify_taxon_to_family(taxon_name: str) -> str:
    """
    Classify a taxon name to a virus family using pattern matching.

    Searches FAMILY_SEARCH_PATTERNS for a case-insensitive word-boundary match.
    First matching family wins. Returns "Other" if no pattern matches.

    Word boundaries prevent false positives like "hiv" matching "archivirus".

    Args:
        taxon_name: The adjusted_taxid_name to classify

    Returns:
        Family name (e.g., "Orthoherpesviridae") or "Other"
    """
    for family, pattern in _FAMILY_PATTERNS_COMPILED:
        if pattern.search(taxon_name):
            return family
    return "Other"


def get_taxa_by_category(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> list[CategorySummary]:
    """
    Group findings by virus family for high-level triage.

    Returns a summary for each virus family detected, showing how many
    samples and contigs belong to that family. Families are determined
    by pattern matching taxon names against FAMILY_SEARCH_PATTERNS.

    Taxa that don't match any known family are grouped under "Other".
    Noise taxa ('root', 'Homo sapiens') are excluded.

    For specific taxon-level details, use get_top_taxa() instead.

    Args:
        sample_set_id: The run to query (required)
        state_dir: Optional state directory override

    Returns:
        List of CategorySummary ordered by sample_count DESC, hit_count DESC.
        Empty list if no hits exist for the sample set.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Build noise filter
    noise_placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)

    # Get all unique (taxon_name, sample_id, hit_key) combinations
    # We need to fetch taxon names to classify them in Python
    sql = f"""
        SELECT DISTINCT
            adjusted_taxid_name as taxon_name,
            sample_id,
            hit_key
        FROM hits
        WHERE sample_set_id = ?
          AND adjusted_taxid_name IS NOT NULL
          AND adjusted_taxid_name NOT IN ({noise_placeholders})
    """  # SQL uses module constants, not user input

    params = [sample_set_id, *DEFAULT_NOISE_TAXA]
    rows = con.execute(sql, params).fetchall()

    if not rows:
        return []

    # Classify each taxon and aggregate
    # Structure: {family: {"samples": set, "hits": set, "taxa": set}}  # noqa: ERA001
    family_data: dict[str, dict[str, set[str]]] = {}

    for taxon_name, sample_id, hit_key in rows:
        family = _classify_taxon_to_family(taxon_name)

        if family not in family_data:
            family_data[family] = {"samples": set(), "hits": set(), "taxa": set()}

        family_data[family]["samples"].add(sample_id)
        family_data[family]["hits"].add(hit_key)
        family_data[family]["taxa"].add(taxon_name)

    # Convert to CategorySummary objects
    results = [
        CategorySummary(
            category=family,
            sample_count=len(data["samples"]),
            hit_count=len(data["hits"]),
            taxa=tuple(sorted(data["taxa"])),
        )
        for family, data in family_data.items()
    ]

    # Sort by sample_count DESC, hit_count DESC, then category name for stability
    results.sort(key=lambda x: (-x.sample_count, -x.hit_count, x.category))

    return results


def get_run_report(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> RunReport | None:
    """
    Generate a comprehensive report for a pipeline run.

    Combines counts, quality metrics, top findings, and per-sample summaries
    into a single RunReport object. This is the primary function for run-level
    reporting, suitable for Slack notifications, CLI display, or JSON export.

    Internally calls get_top_taxa(), get_sample_summaries(), get_negative_samples(),
    and get_contig_quality() to gather data.

    Args:
        sample_set_id: The run to report on (required)
        state_dir: Optional state directory override

    Returns:
        RunReport with all computed metrics, or None if no hits exist
        for the sample set.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return None

    con = query_hits(state_dir)

    # Check if this sample set has any data
    count_result = con.execute(
        "SELECT COUNT(*) FROM hits WHERE sample_set_id = ?",
        [sample_set_id],
    ).fetchone()
    assert count_result is not None, "COUNT query should always return a row"

    if count_result[0] == 0:
        return None

    # Get basic counts and run_date in one query
    # Note: We cast adjusted_taxid_name to VARCHAR to handle edge cases where
    # parquet schema inference produces INTEGER when all values are NULL.
    noise_placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)

    summary_sql = f"""
        SELECT
            COUNT(DISTINCT sample_id) as samples_analyzed,
            COUNT(DISTINCT hit_key) as unique_hits,
            MIN(run_date) as run_date,
            COUNT(DISTINCT adjusted_taxid_name) FILTER (
                WHERE adjusted_taxid_name IS NOT NULL
                  AND CAST(adjusted_taxid_name AS VARCHAR) NOT IN ({noise_placeholders})
            ) as viral_taxa_found
        FROM hits
        WHERE sample_set_id = ?
    """  # SQL uses module constants, not user input

    summary_params = [*DEFAULT_NOISE_TAXA, sample_set_id]
    summary_row = con.execute(summary_sql, summary_params).fetchone()
    assert summary_row is not None, "summary query should always return a row"

    samples_analyzed = summary_row[0]
    unique_hits = summary_row[1]
    run_date = summary_row[2]
    viral_taxa_found = summary_row[3]

    # Get quality metrics
    quality = get_contig_quality(sample_set_id, state_dir)
    median_contig_length = quality.median_length if quality else 0
    contigs_over_500bp = quality.contigs_500bp_plus if quality else 0

    # Get top findings (top 5 taxa)
    top_findings = get_top_taxa(
        sample_set_id=sample_set_id,
        limit=5,
        exclude_noise=True,
        state_dir=state_dir,
    )

    # Get per-sample summaries
    sample_summaries = get_sample_summaries(sample_set_id, state_dir)

    # Get negative samples
    negative_samples = get_negative_samples(sample_set_id, state_dir)

    # Calculate samples with viral hits
    samples_negative = len(negative_samples)
    samples_with_viral_hits = samples_analyzed - samples_negative

    return RunReport(
        sample_set_id=sample_set_id,
        run_date=run_date,
        samples_analyzed=samples_analyzed,
        samples_with_viral_hits=samples_with_viral_hits,
        samples_negative=samples_negative,
        unique_hits=unique_hits,
        viral_taxa_found=viral_taxa_found,
        median_contig_length=median_contig_length,
        contigs_over_500bp=contigs_over_500bp,
        top_findings=tuple(top_findings),
        sample_summaries=tuple(sample_summaries),
    )


def get_novel_taxa(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> list[NovelTaxon]:
    """
    Get taxa in this run that have never appeared in any previous run.

    Identifies first-time sightings by comparing taxa in the specified run
    against all taxa seen in prior runs. Useful for detecting emerging
    pathogens, new contamination sources, or simply tracking discovery.

    Noise taxa ('root', 'Homo sapiens') are excluded.

    Args:
        sample_set_id: The run to analyze (required)
        state_dir: Optional state directory override

    Returns:
        List of NovelTaxon ordered by sample_count DESC, hit_count DESC.
        Empty list if no novel taxa or no hits exist.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Build noise filter
    noise_placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)

    # Find taxa in current run that don't exist in any other run.
    # Optimized: compute counts in first CTE, use NOT EXISTS instead of
    # NOT IN + re-join pattern. This avoids scanning hits twice for counts.
    sql = f"""
        WITH current_taxa AS (
            SELECT
                adjusted_taxid as taxid,
                adjusted_taxid_name as name,
                adjusted_taxid_rank as rank,
                COUNT(DISTINCT sample_id) as sample_count,
                COUNT(DISTINCT hit_key) as hit_count
            FROM hits
            WHERE sample_set_id = ?
              AND adjusted_taxid_name IS NOT NULL
              AND adjusted_taxid_name NOT IN ({noise_placeholders})
            GROUP BY adjusted_taxid, adjusted_taxid_name, adjusted_taxid_rank
        ),
        historical_taxa AS (
            SELECT DISTINCT adjusted_taxid_name as name
            FROM hits
            WHERE sample_set_id != ?
              AND adjusted_taxid_name IS NOT NULL
              AND adjusted_taxid_name NOT IN ({noise_placeholders})
        )
        SELECT
            c.taxid,
            c.name,
            c.rank,
            c.sample_count,
            c.hit_count
        FROM current_taxa c
        WHERE NOT EXISTS (SELECT 1 FROM historical_taxa h WHERE h.name = c.name)
        ORDER BY c.sample_count DESC, c.hit_count DESC
    """  # SQL uses module constants, not user input

    # Parameters: sample_set_id, noise taxa, sample_set_id, noise taxa
    params = [
        sample_set_id,
        *DEFAULT_NOISE_TAXA,
        sample_set_id,
        *DEFAULT_NOISE_TAXA,
    ]

    rows = con.execute(sql, params).fetchall()

    result = [
        NovelTaxon(
            taxid=row[0],
            name=row[1],
            rank=row[2],
            sample_count=row[3],
            hit_count=row[4],
        )
        for row in rows
    ]

    assert all(t.name not in DEFAULT_NOISE_TAXA for t in result), (
        "noise taxa should be filtered from results"
    )
    return result


def get_top_movers(
    sample_set_id: str,
    limit: int = 10,
    state_dir: Path | str | None = None,
) -> list[TaxonChange]:
    """
    Get taxa with the largest prevalence change vs historical baseline.

    Compares each taxon's prevalence in the specified run against its
    historical average prevalence across all prior runs. Returns taxa
    sorted by absolute change (largest movers first).

    Prevalence is calculated as: (samples with taxon / total samples) * 100

    Noise taxa ('root', 'Homo sapiens') are excluded.

    Args:
        sample_set_id: The run to analyze (required)
        limit: Maximum number of taxa to return (default: 10)
        state_dir: Optional state directory override

    Returns:
        List of TaxonChange ordered by |change_pct| DESC.
        Empty list if no hits exist or this is the first run.
    """
    assert sample_set_id, "sample_set_id cannot be empty"
    assert limit > 0, "limit must be positive"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Check if there are any prior runs
    prior_runs = con.execute(
        "SELECT COUNT(DISTINCT sample_set_id) FROM hits WHERE sample_set_id != ?",
        [sample_set_id],
    ).fetchone()

    if prior_runs is None or prior_runs[0] == 0:
        return []

    baseline_run_count = prior_runs[0]

    # Build noise filter
    noise_placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)

    sql = f"""
        WITH
        -- Sample counts per run (for prevalence calculation)
        run_sample_counts AS (
            SELECT sample_set_id, COUNT(DISTINCT sample_id) as sample_count
            FROM hits
            GROUP BY sample_set_id
        ),
        -- Current run's taxa with prevalence
        current_prevalence AS (
            SELECT
                adjusted_taxid as taxid,
                adjusted_taxid_name as name,
                adjusted_taxid_rank as rank,
                COUNT(DISTINCT sample_id) as sample_count,
                COUNT(DISTINCT sample_id) * 100.0 /
                    (SELECT sample_count FROM run_sample_counts WHERE sample_set_id = ?) as prevalence_pct
            FROM hits
            WHERE sample_set_id = ?
              AND adjusted_taxid_name IS NOT NULL
              AND adjusted_taxid_name NOT IN ({noise_placeholders})
            GROUP BY adjusted_taxid, adjusted_taxid_name, adjusted_taxid_rank
        ),
        -- Historical prevalence per taxon (average across prior runs)
        historical_prevalence AS (
            SELECT
                adjusted_taxid_name as name,
                AVG(prevalence_pct) as avg_prevalence_pct
            FROM (
                SELECT
                    sample_set_id,
                    adjusted_taxid_name,
                    COUNT(DISTINCT sample_id) * 100.0 / r.sample_count as prevalence_pct
                FROM hits h
                JOIN run_sample_counts r USING (sample_set_id)
                WHERE sample_set_id != ?
                  AND adjusted_taxid_name IS NOT NULL
                  AND adjusted_taxid_name NOT IN ({noise_placeholders})
                GROUP BY sample_set_id, adjusted_taxid_name, r.sample_count
            )
            GROUP BY adjusted_taxid_name
        )
        SELECT
            c.taxid,
            c.name,
            c.rank,
            c.prevalence_pct as current_prevalence_pct,
            COALESCE(h.avg_prevalence_pct, 0) as baseline_prevalence_pct,
            c.prevalence_pct - COALESCE(h.avg_prevalence_pct, 0) as change_pct,
            c.sample_count as current_sample_count
        FROM current_prevalence c
        LEFT JOIN historical_prevalence h ON c.name = h.name
        ORDER BY ABS(change_pct) DESC
        LIMIT {limit}
    """  # SQL uses module constants, not user input

    # Parameters for the query
    params = [
        sample_set_id,  # run_sample_counts subquery
        sample_set_id,  # current_prevalence WHERE
        *DEFAULT_NOISE_TAXA,  # current_prevalence noise filter
        sample_set_id,  # historical_prevalence WHERE
        *DEFAULT_NOISE_TAXA,  # historical_prevalence noise filter
    ]

    rows = con.execute(sql, params).fetchall()

    result = [
        TaxonChange(
            taxid=row[0],
            name=row[1],
            rank=row[2],
            current_prevalence_pct=row[3],
            baseline_prevalence_pct=row[4],
            change_pct=row[5],
            current_sample_count=row[6],
            baseline_run_count=baseline_run_count,
        )
        for row in rows
    ]

    assert all(t.name not in DEFAULT_NOISE_TAXA for t in result), (
        "noise taxa should be filtered from results"
    )
    return result


def get_rare_taxa(
    sample_set_id: str,
    threshold_pct: float = 5.0,
    state_dir: Path | str | None = None,
) -> list[RareTaxon]:
    """
    Get taxa in this run that appear in few historical runs.

    Identifies taxa present in the current run that have historically
    appeared in less than threshold_pct of all runs. Useful for flagging
    unusual findings that warrant closer inspection.

    Noise taxa ('root', 'Homo sapiens') are excluded.

    Args:
        sample_set_id: The run to analyze (required)
        threshold_pct: Maximum historical run percentage to be considered rare
                       (default: 5.0, meaning taxa in <5% of runs)
        state_dir: Optional state directory override

    Returns:
        List of RareTaxon ordered by historical_run_pct ASC (rarest first),
        then by current_sample_count DESC.
        Empty list if no rare taxa or no hits exist.
    """
    assert sample_set_id, "sample_set_id cannot be empty"
    assert 0 < threshold_pct <= MAX_PERCENTAGE, (
        "threshold_pct must be between 0 and 100"
    )

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Build noise filter
    noise_placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)

    sql = f"""
        WITH
        -- Total number of runs in the system
        total_runs AS (
            SELECT COUNT(DISTINCT sample_set_id) as count FROM hits
        ),
        -- Count of runs where each taxon appears (across all runs)
        taxa_run_counts AS (
            SELECT
                adjusted_taxid_name as name,
                COUNT(DISTINCT sample_set_id) as run_count
            FROM hits
            WHERE adjusted_taxid_name IS NOT NULL
              AND adjusted_taxid_name NOT IN ({noise_placeholders})
            GROUP BY adjusted_taxid_name
        ),
        -- Taxa in the current run with their details
        current_taxa AS (
            SELECT
                adjusted_taxid as taxid,
                adjusted_taxid_name as name,
                adjusted_taxid_rank as rank,
                COUNT(DISTINCT sample_id) as sample_count
            FROM hits
            WHERE sample_set_id = ?
              AND adjusted_taxid_name IS NOT NULL
              AND adjusted_taxid_name NOT IN ({noise_placeholders})
            GROUP BY adjusted_taxid, adjusted_taxid_name, adjusted_taxid_rank
        )
        SELECT
            c.taxid,
            c.name,
            c.rank,
            t.run_count as historical_run_count,
            t.run_count * 100.0 / r.count as historical_run_pct,
            c.sample_count as current_sample_count
        FROM current_taxa c
        JOIN taxa_run_counts t ON c.name = t.name
        CROSS JOIN total_runs r
        WHERE t.run_count * 100.0 / r.count < ?
        ORDER BY historical_run_pct ASC, current_sample_count DESC
    """  # SQL uses module constants, not user input

    # Parameters
    params = [
        *DEFAULT_NOISE_TAXA,  # taxa_run_counts noise filter
        sample_set_id,  # current_taxa WHERE
        *DEFAULT_NOISE_TAXA,  # current_taxa noise filter
        threshold_pct,  # threshold filter
    ]

    rows = con.execute(sql, params).fetchall()

    result = [
        RareTaxon(
            taxid=row[0],
            name=row[1],
            rank=row[2],
            historical_run_count=row[3],
            historical_run_pct=row[4],
            current_sample_count=row[5],
        )
        for row in rows
    ]

    assert all(t.name not in DEFAULT_NOISE_TAXA for t in result), (
        "noise taxa should be filtered from results"
    )
    assert all(t.historical_run_pct < threshold_pct for t in result), (
        "all results should be below threshold"
    )
    return result


def get_taxon_history(
    taxon_name: str,
    state_dir: Path | str | None = None,
) -> TaxonHistory | None:
    """
    Get the complete history of a taxon across all runs.

    Provides a detailed track record showing when and how often a taxon
    has been detected. Useful for understanding whether a finding is
    typical or unusual for this dataset.

    Args:
        taxon_name: The taxon name to look up (case-sensitive)
        state_dir: Optional state directory override

    Returns:
        TaxonHistory with per-run details, or None if the taxon has
        never been seen.
    """
    assert taxon_name, "taxon_name cannot be empty"

    if not has_any_hits(state_dir):
        return None

    con = query_hits(state_dir)

    # Check if this taxon exists
    exists = con.execute(
        "SELECT COUNT(*) FROM hits WHERE adjusted_taxid_name = ?",
        [taxon_name],
    ).fetchone()

    if exists is None or exists[0] == 0:
        return None

    # Get taxon metadata (taxid, rank) from most recent observation
    taxon_info = con.execute(
        """
        SELECT adjusted_taxid, adjusted_taxid_rank
        FROM hits
        WHERE adjusted_taxid_name = ?
        ORDER BY run_date DESC
        LIMIT 1
        """,
        [taxon_name],
    ).fetchone()

    taxid = taxon_info[0] if taxon_info else None
    rank = taxon_info[1] if taxon_info else None

    # Get total runs in system
    total_runs_result = con.execute(
        "SELECT COUNT(DISTINCT sample_set_id) FROM hits",
    ).fetchone()
    total_runs_in_system = total_runs_result[0] if total_runs_result else 0

    # Get per-run presence with prevalence
    run_history_sql = """
        WITH run_sample_counts AS (
            SELECT sample_set_id, COUNT(DISTINCT sample_id) as total_samples
            FROM hits
            GROUP BY sample_set_id
        )
        SELECT
            h.sample_set_id,
            MIN(h.run_date) as run_date,
            COUNT(DISTINCT h.sample_id) as sample_count,
            COUNT(DISTINCT h.sample_id) * 100.0 / r.total_samples as prevalence_pct
        FROM hits h
        JOIN run_sample_counts r ON h.sample_set_id = r.sample_set_id
        WHERE h.adjusted_taxid_name = ?
        GROUP BY h.sample_set_id, r.total_samples
        ORDER BY run_date ASC
    """

    run_rows = con.execute(run_history_sql, [taxon_name]).fetchall()

    if not run_rows:
        return None

    run_history = tuple(
        TaxonRunPresence(
            sample_set_id=row[0],
            run_date=row[1],
            sample_count=row[2],
            prevalence_pct=row[3],
        )
        for row in run_rows
    )

    total_runs_seen = len(run_history)
    first_seen_date = run_history[0].run_date
    last_seen_date = run_history[-1].run_date
    overall_prevalence_pct = (
        (total_runs_seen / total_runs_in_system * 100)
        if total_runs_in_system > 0
        else 0.0
    )

    result = TaxonHistory(
        taxid=taxid,
        name=taxon_name,
        rank=rank,
        total_runs_seen=total_runs_seen,
        total_runs_in_system=total_runs_in_system,
        overall_prevalence_pct=overall_prevalence_pct,
        first_seen_date=first_seen_date,
        last_seen_date=last_seen_date,
        run_history=run_history,
    )

    assert result.total_runs_seen == len(result.run_history), (
        "total_runs_seen must match run_history length"
    )
    assert result.total_runs_seen <= result.total_runs_in_system, (
        "cannot have seen more runs than exist"
    )
    return result


def get_run_comparison(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> RunComparison | None:
    """
    Compare a run's metrics to historical norms.

    Provides context for interpreting a run by comparing its key metrics
    (sample count, hit count, taxa count, contig length) against the
    historical average across all prior runs.

    Args:
        sample_set_id: The run to analyze (required)
        state_dir: Optional state directory override

    Returns:
        RunComparison with current metrics and historical baselines,
        or None if the run doesn't exist.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return None

    con = query_hits(state_dir)

    # Check if this run exists
    exists = con.execute(
        "SELECT COUNT(*) FROM hits WHERE sample_set_id = ?",
        [sample_set_id],
    ).fetchone()

    if exists is None or exists[0] == 0:
        return None

    # Build noise filter for taxa counting
    noise_placeholders = ", ".join("?" for _ in DEFAULT_NOISE_TAXA)

    # Get current run metrics
    current_sql = f"""
        SELECT
            COUNT(DISTINCT sample_id) as samples_analyzed,
            COUNT(DISTINCT hit_key) as unique_hits,
            COUNT(DISTINCT adjusted_taxid_name) FILTER (
                WHERE adjusted_taxid_name IS NOT NULL
                  AND adjusted_taxid_name NOT IN ({noise_placeholders})
            ) as taxa_found,
            COALESCE(MEDIAN(sequence_length), 0)::INTEGER as median_contig_length,
            MIN(run_date) as run_date
        FROM hits
        WHERE sample_set_id = ?
    """  # SQL uses module constants, not user input

    current_params = [*DEFAULT_NOISE_TAXA, sample_set_id]
    current_row = con.execute(current_sql, current_params).fetchone()

    if current_row is None:
        return None

    samples_analyzed = current_row[0]
    unique_hits = current_row[1]
    taxa_found = current_row[2]
    median_contig_length = current_row[3]
    run_date = current_row[4]

    # Get historical metrics (from all runs except current)
    historical_sql = f"""
        WITH run_metrics AS (
            SELECT
                sample_set_id,
                COUNT(DISTINCT sample_id) as samples,
                COUNT(DISTINCT hit_key) as hits,
                COUNT(DISTINCT adjusted_taxid_name) FILTER (
                    WHERE adjusted_taxid_name IS NOT NULL
                      AND adjusted_taxid_name NOT IN ({noise_placeholders})
                ) as taxa,
                COALESCE(MEDIAN(sequence_length), 0) as median_length
            FROM hits
            WHERE sample_set_id != ?
            GROUP BY sample_set_id
        )
        SELECT
            COUNT(*) as run_count,
            COALESCE(AVG(samples), 0) as avg_samples,
            COALESCE(AVG(hits), 0) as avg_hits,
            COALESCE(AVG(taxa), 0) as avg_taxa,
            COALESCE(AVG(median_length), 0) as avg_median_length
        FROM run_metrics
    """  # SQL uses module constants, not user input

    historical_params = [*DEFAULT_NOISE_TAXA, sample_set_id]
    historical_row = con.execute(historical_sql, historical_params).fetchone()

    historical_run_count = historical_row[0] if historical_row else 0
    avg_samples = historical_row[1] if historical_row else 0.0
    avg_hits = historical_row[2] if historical_row else 0.0
    avg_taxa = historical_row[3] if historical_row else 0.0
    avg_median_length = historical_row[4] if historical_row else 0.0

    result = RunComparison(
        sample_set_id=sample_set_id,
        run_date=run_date,
        samples_analyzed=samples_analyzed,
        unique_hits=unique_hits,
        taxa_found=taxa_found,
        median_contig_length=median_contig_length,
        historical_run_count=historical_run_count,
        avg_samples=avg_samples,
        avg_hits=avg_hits,
        avg_taxa=avg_taxa,
        avg_median_length=avg_median_length,
    )

    assert result.samples_analyzed >= 0, "samples_analyzed must be non-negative"
    assert result.historical_run_count >= 0, "historical_run_count must be non-negative"
    return result


def get_recurring_hits(
    min_samples: int = 2,
    min_runs: int | None = None,
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[RecurringHit]:
    """
    Find hits that appear in multiple samples or runs.

    Recurring hits may indicate persistent infections, contamination,
    or lab reagent artifacts. This query helps identify sequences that
    warrant further investigation.

    Classification is the most common (mode) adjusted_taxid across all
    observations for each hit.

    Args:
        min_samples: Minimum number of distinct samples (default: 2)
        min_runs: Minimum number of distinct runs (sample_set_ids), or None for no filter
        limit: Maximum number of results to return, or None for all
        state_dir: Optional state directory override

    Returns:
        List of RecurringHit objects, ordered by sample_count DESC, then
        observation_count DESC. Empty list if no hits meet the criteria.
    """
    assert min_samples >= 1, "min_samples must be at least 1"
    assert min_runs is None or min_runs >= 1, "min_runs must be at least 1 if provided"
    assert limit is None or limit > 0, "limit must be positive if provided"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # Build query with most common classification per hit
    # MODE() returns the most frequent value; we use it on adjusted_taxid
    # and then get the corresponding name/rank from any row with that taxid
    sql = """
        WITH hit_stats AS (
            SELECT
                hit_key,
                ANY_VALUE(sequence_length) as sequence_length,
                ANY_VALUE(gc_content) as gc_content,
                MIN(run_date) as first_seen_date,
                MAX(run_date) as last_seen_date,
                COUNT(DISTINCT sample_id) as sample_count,
                COUNT(DISTINCT sample_set_id) as run_count,
                COUNT(*) as observation_count,
                MODE(adjusted_taxid) as top_taxid
            FROM hits
            GROUP BY hit_key
            HAVING sample_count >= ?
        ),
        taxid_names AS (
            -- Get one (name, rank) pair per (hit_key, adjusted_taxid) combination.
            -- ORDER BY ensures deterministic selection when multiple rows exist.
            SELECT DISTINCT ON (hit_key, adjusted_taxid)
                hit_key,
                adjusted_taxid,
                adjusted_taxid_name,
                adjusted_taxid_rank
            FROM hits
            WHERE adjusted_taxid IS NOT NULL
            ORDER BY hit_key, adjusted_taxid
        )
        SELECT
            s.hit_key,
            s.sequence_length,
            s.gc_content,
            s.first_seen_date,
            s.last_seen_date,
            s.sample_count,
            s.run_count,
            s.observation_count,
            s.top_taxid,
            t.adjusted_taxid_name,
            t.adjusted_taxid_rank
        FROM hit_stats s
        LEFT JOIN taxid_names t ON s.hit_key = t.hit_key AND s.top_taxid = t.adjusted_taxid
    """
    params: list[int] = [min_samples]

    if min_runs is not None:
        sql += " WHERE s.run_count >= ?"
        params.append(min_runs)

    sql += " ORDER BY s.sample_count DESC, s.observation_count DESC"

    if limit is not None:
        # f-string is safe here because limit is validated as positive int above
        sql += f" LIMIT {limit}"

    rows = con.execute(sql, params).fetchall()

    result = [
        RecurringHit(
            hit_key=row[0],
            sequence_length=row[1],
            gc_content=row[2],
            first_seen_date=row[3],
            last_seen_date=row[4],
            sample_count=row[5],
            run_count=row[6],
            observation_count=row[7],
            top_taxid=row[8],
            top_taxid_name=row[9],
            top_taxid_rank=row[10],
        )
        for row in rows
    ]
    assert all(isinstance(r, RecurringHit) for r in result), (
        "result must contain only RecurringHit objects"
    )
    return result


def get_discovery_timeline(
    granularity: Literal["day", "week", "month", "year"] = "month",
    state_dir: Path | str | None = None,
) -> list[TimelineBucket]:
    """
    Count new hits discovered per time period.

    Groups hits by their first_seen_date (computed as MIN(run_date) per hit)
    and counts how many were discovered in each period.

    Args:
        granularity: Time period grouping - "day", "week", "month", or "year"
        state_dir: Optional state directory override

    Returns:
        List of TimelineBucket objects ordered by period ascending.
        Empty list if no hits exist.
    """
    valid_granularities = {"day", "week", "month", "year"}
    assert granularity in valid_granularities, (
        f"granularity must be one of {valid_granularities}"
    )

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # DuckDB strftime format strings for each granularity
    # Note: DuckDB strftime takes (TIMESTAMP, format) not (format, TIMESTAMP) like SQLite
    format_map = {
        "day": "%Y-%m-%d",
        "week": "%Y-W%W",  # ISO week number
        "month": "%Y-%m",
        "year": "%Y",
    }
    date_format = format_map[granularity]

    # First compute first_seen_date per hit, then group by period
    # Cast the string date to TIMESTAMP for strftime
    sql = f"""
        SELECT
            strftime(CAST(first_seen_date AS TIMESTAMP), '{date_format}') as period,
            COUNT(*) as new_hits
        FROM (
            SELECT hit_key, MIN(run_date) as first_seen_date
            FROM hits
            GROUP BY hit_key
        )
        GROUP BY period
        ORDER BY period ASC
    """  # SQL uses module constants, not user input

    rows = con.execute(sql).fetchall()

    result = [
        TimelineBucket(
            period=row[0],
            new_hits=row[1],
        )
        for row in rows
    ]
    assert all(isinstance(b, TimelineBucket) for b in result), (
        "result must contain only TimelineBucket objects"
    )
    return result


def list_hits_with_observations(
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[tuple[Hit, HitObservation]]:
    """
    List hits joined with their observations.

    Returns one tuple per observation, so a hit observed in multiple samples
    will appear multiple times (once per observation).

    Note:
        Unlike get_hit() and list_hits() which return the *most recent*
        classification, this function returns the classification from each
        *specific observation*. This is the per-observation semantic, not
        aggregated. Use this when you need to see how classification may
        have changed over time for a given hit.

    Args:
        limit: Maximum number of rows to return
        state_dir: Optional state directory override

    Returns:
        List of (Hit, HitObservation) tuples, ordered by first_seen_date DESC
    """
    assert limit is None or limit > 0, "limit must be positive if provided"

    if not has_any_hits(state_dir):
        return []

    con = query_hits(state_dir)

    # We need to compute first_seen_date per hit and join back
    # Include classification from each observation
    sql = """
        SELECT
            h.hit_key,
            h.sequence_length,
            h.sequence_compressed,
            h.gc_content,
            first_seen.first_seen_date,
            h.sample_set_id,
            h.sample_id,
            h.run_date,
            h.sequence_id,
            h.adjusted_taxid,
            h.adjusted_taxid_name,
            h.adjusted_taxid_rank,
            h.source
        FROM hits h
        INNER JOIN (
            SELECT hit_key, MIN(run_date) as first_seen_date
            FROM hits
            GROUP BY hit_key
        ) first_seen ON h.hit_key = first_seen.hit_key
        ORDER BY first_seen.first_seen_date DESC, h.sample_id
    """

    if limit is not None:
        # f-string is safe here because limit is validated as positive int above
        sql += f" LIMIT {limit}"

    rows = con.execute(sql).fetchall()

    result = []
    for row in rows:
        hit = Hit(
            hit_key=row[0],
            sequence_length=row[1],
            sequence_compressed=row[2],
            gc_content=row[3],
            first_seen_date=row[4],
            top_taxid=row[9],
            top_taxid_name=row[10],
            top_taxid_rank=row[11],
        )
        observation = HitObservation(
            hit_key=row[0],
            sample_set_id=row[5],
            sample_id=row[6],
            run_date=row[7],
            sequence_id=row[8],
            source=row[12],
        )
        result.append((hit, observation))

    assert all(isinstance(t, tuple) and len(t) == 2 for t in result), (  # noqa: PLR2004
        "result must be list of 2-tuples"
    )
    return result


def count_sample_set_observations(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> int:
    """
    Count hit observations for a specific sample set.

    Searches both uncompacted (month=NULL/{sample_set_id}/) and compacted
    (month=*/) partitions for observations matching the sample_set_id.

    Args:
        sample_set_id: The sample set identifier
        state_dir: Optional state directory override

    Returns:
        Number of observations (rows) in parquet files for this sample set
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    if not has_any_hits(state_dir):
        return 0

    con = query_hits(state_dir)
    result = con.execute(
        "SELECT COUNT(*) FROM hits WHERE sample_set_id = ?",
        [sample_set_id],
    ).fetchone()

    assert result is not None, "COUNT query should always return a row"
    count = result[0]
    assert isinstance(count, int), f"count must be int, got {type(count)}"
    assert count >= 0, f"count must be non-negative, got {count}"
    return count


def delete_sample_hits(
    sample_id: str,
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> bool:
    """
    Delete the parquet file for a specific sample (uncompacted data only).

    Deletes the entire sample directory:
        hits/month=NULL/{sample_set_id}/{sample_id}/

    Args:
        sample_id: The sample identifier
        sample_set_id: The sample set identifier
        state_dir: Optional state directory override

    Returns:
        True if the directory existed and was deleted, False if it didn't exist
    """
    assert sample_id, "sample_id cannot be empty"
    assert sample_set_id, "sample_set_id cannot be empty"

    hits_dir = get_hits_dir(state_dir)
    sample_dir = hits_dir / "month=NULL" / sample_set_id / sample_id

    if not sample_dir.exists():
        return False

    shutil.rmtree(sample_dir)
    return True


def _delete_sample_set_hits_uncompacted(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> int:
    """
    Delete uncompacted parquet files for a sample set.

    Removes the sample set subdirectory under month=NULL/:
        hits/month=NULL/{sample_set_id}/

    Args:
        sample_set_id: The sample set identifier
        state_dir: Optional state directory override

    Returns:
        Number of parquet files that were deleted
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    hits_dir = get_hits_dir(state_dir)
    sample_set_dir = hits_dir / "month=NULL" / sample_set_id

    if not sample_set_dir.exists():
        return 0

    parquet_files = list(sample_set_dir.rglob("*.parquet"))
    count = len(parquet_files)

    assert sample_set_dir.is_dir(), f"expected {sample_set_dir} to be a directory"

    shutil.rmtree(sample_set_dir)

    return count


def _delete_sample_set_hits_compacted(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> int:
    """
    Delete compacted hits for a sample set by rewriting affected month files.

    Scans month=YYYY-MM directories (skipping month=NULL) and rewrites any
    that contain data for the given sample_set_id, excluding those rows.

    Args:
        sample_set_id: The sample set identifier
        state_dir: Optional state directory override

    Returns:
        Number of month files that were rewritten
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    hits_dir = get_hits_dir(state_dir)
    assert hits_dir.exists(), f"hits directory must exist: {hits_dir}"

    months_rewritten = 0

    for month_dir in hits_dir.glob("month=*"):
        if month_dir.name == "month=NULL":
            continue

        data_file = month_dir / "data.parquet"
        if not data_file.exists():
            continue

        con = duckdb.connect()

        # Check if this sample_set has data in this month
        # Note: file path must be interpolated (DuckDB limitation), but
        # sample_set_id is parameterized to avoid injection
        has_data_row = con.execute(
            f"""
            SELECT COUNT(*) > 0
            FROM read_parquet('{data_file}')
            WHERE sample_set_id = $1
            """,  # SQL uses module constants, not user input
            [sample_set_id],
        ).fetchone()
        assert has_data_row is not None, "COUNT query should always return a row"
        has_data = has_data_row[0]

        if not has_data:
            con.close()
            continue

        temp_file = month_dir / ".data.parquet.tmp"
        con.execute(
            f"""
            COPY (
                SELECT * FROM read_parquet('{data_file}')
                WHERE sample_set_id != $1
                ORDER BY hit_key, run_date
            )
            TO '{temp_file}'
            (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000)
            """,  # SQL uses module constants, not user input
            [sample_set_id],
        )

        row_count_row = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{temp_file}')",  # SQL uses module constants, not user input
        ).fetchone()
        assert row_count_row is not None, "COUNT query should always return a row"
        row_count = row_count_row[0]

        con.close()

        if row_count == 0:
            temp_file.unlink()
            data_file.unlink()
            month_dir.rmdir()
        else:
            temp_file.rename(data_file)

        months_rewritten += 1

    return months_rewritten


def delete_sample_set_hits(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> DeleteResult:
    """
    Delete all hits for a sample set (both uncompacted and compacted).

    For uncompacted data: removes the sample set subdirectory under month=NULL/.
    For compacted data: rewrites affected month files to exclude the sample_set_id.

    Args:
        sample_set_id: The sample set identifier
        state_dir: Optional state directory override

    Returns:
        DeleteResult with counts of uncompacted files deleted and compacted
        months rewritten.
    """
    assert sample_set_id, "sample_set_id cannot be empty"

    uncompacted = _delete_sample_set_hits_uncompacted(sample_set_id, state_dir)
    compacted = _delete_sample_set_hits_compacted(sample_set_id, state_dir)

    assert uncompacted >= 0, "uncompacted count cannot be negative"
    assert compacted >= 0, "compacted count cannot be negative"

    return DeleteResult(
        uncompacted_files_deleted=uncompacted,
        compacted_months_rewritten=compacted,
    )


@dataclass(frozen=True)
class _FileSnapshot:
    """A frozen snapshot of an uncompacted parquet file at inventory time.

    Used to ensure compaction operates on a fixed set of files and only
    deletes files that haven't been modified since the snapshot was taken.
    """

    path: Path
    mtime_ns: int
    size: int


def _snapshot_uncompacted_files(
    state_dir: Path | str | None = None,
) -> list[_FileSnapshot]:
    """Snapshot all uncompacted parquet files with their stat metadata.

    Returns a frozen list of (path, mtime_ns, size) tuples representing
    the exact files present at call time. Subsequent compaction steps
    should operate only on these files, ignoring any that arrive later.
    """
    hits_dir = get_hits_dir(state_dir)
    uncompacted_dir = hits_dir / "month=NULL"

    if not uncompacted_dir.exists():
        return []

    snapshots = []
    for p in uncompacted_dir.rglob("*.parquet"):
        st = p.stat()
        snapshots.append(
            _FileSnapshot(path=p, mtime_ns=st.st_mtime_ns, size=st.st_size),
        )
    return snapshots


@contextmanager
def _compaction_lock(hits_dir: Path) -> Generator[None]:
    """Acquire an exclusive advisory lock for compaction. Fails fast if held.

    This prevents two concurrent ``compact_hits`` invocations from
    clobbering each other's temp files or double-deleting source files.
    It does not affect Nextflow writer processes.
    """
    lock_path = hits_dir / ".compact.lock"
    lock_file = lock_path.open("w")
    try:
        fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as e:
        lock_file.close()
        msg = f"Another compaction is already running (lock: {lock_path})"
        raise RuntimeError(msg) from e
    try:
        yield
    finally:
        fcntl.flock(lock_file, fcntl.LOCK_UN)
        lock_file.close()
        with suppress(OSError):
            lock_path.unlink()


def _get_uncompacted_glob(state_dir: Path | str | None = None) -> str:
    """Get glob pattern for uncompacted parquet files."""
    hits_dir = get_hits_dir(state_dir)
    result = str(hits_dir / "month=NULL" / "**" / "*.parquet")
    assert "month=NULL" in result, "glob pattern must target uncompacted partition"
    return result


def _inventory_uncompacted(
    file_paths: list[Path],
    month: str | None = None,
) -> list[MonthCompactionInfo]:
    """
    Inventory uncompacted data grouped by month.

    Operates on an explicit list of file paths rather than globbing, so
    that the caller can freeze the file set once and use it consistently
    across inventory, compaction, and deletion.

    Args:
        file_paths: Explicit list of uncompacted parquet file paths to inventory.
        month: Optional filter to a specific month (YYYY-MM format)

    Returns:
        List of MonthCompactionInfo for each month with uncompacted data
    """
    assert month is None or (len(month) == MONTH_FORMAT_LENGTH and month[4] == "-"), (
        f"month must be YYYY-MM format, got {month!r}"
    )

    if not file_paths:
        return []

    file_list_sql = "[" + ",".join(f"'{p}'" for p in file_paths) + "]"

    con = duckdb.connect()

    month_filter = ""
    if month is not None:
        month_filter = f"WHERE strftime(run_date::DATE, '%Y-%m') = '{month}'"

    result = con.execute(
        f"""
        SELECT
            strftime(run_date::DATE, '%Y-%m') AS target_month,
            COUNT(*) AS observation_count,
            COUNT(DISTINCT hit_key) AS unique_hits,
            COUNT(DISTINCT sample_set_id) AS sample_sets
        FROM read_parquet({file_list_sql})
        {month_filter}
        GROUP BY target_month
        ORDER BY target_month DESC
    """,  # SQL uses module constants, not user input
    ).fetchall()

    inventory = []
    for row in result:
        target_month, obs_count, unique_hits, sample_sets = row

        file_count = con.execute(
            f"""
            SELECT COUNT(DISTINCT filename)
            FROM read_parquet({file_list_sql}, filename=true)
            WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'
        """,  # SQL uses module constants, not user input
        ).fetchone()

        inventory.append(
            MonthCompactionInfo(
                month=target_month,
                observation_count=obs_count,
                unique_hits=unique_hits,
                sample_sets=sample_sets,
                source_files=file_count[0] if file_count else 0,
            ),
        )

    assert all(isinstance(m, MonthCompactionInfo) for m in inventory), (
        "inventory must contain only MonthCompactionInfo objects"
    )
    return inventory


def compact_hits(  # noqa: C901, PLR0912, PLR0915  # complexity is inherent to this function
    state_dir: Path | str | None = None,
    month: str | None = None,
    dry_run: bool = False,  # noqa: FBT001, FBT002  # existing public API
    keep_source: bool = False,  # noqa: FBT001, FBT002  # existing public API
) -> CompactionResult:
    """
    Compact uncompacted hits into monthly partitions.

    Reads data from hits/month=NULL/**, groups by month (from run_date),
    sorts by hit_key, and writes to hits/month={YYYY-MM}/data.parquet.

    If a compacted file already exists for a month, the new data is merged
    with the existing data.

    Args:
        state_dir: Optional state directory override
        month: Only compact data for this month (YYYY-MM format).
               Default: compact all uncompacted data.
        dry_run: If True, show what would be compacted without doing it.
        keep_source: If True, keep uncompacted files after compaction.
                     Default: delete after successful compaction.

    Returns:
        CompactionResult with details of what was (or would be) compacted.
    """
    assert month is None or (len(month) == MONTH_FORMAT_LENGTH and month[4] == "-"), (
        f"month must be YYYY-MM format, got {month!r}"
    )

    hits_dir = get_hits_dir(state_dir)

    # Freeze the set of uncompacted files. Every subsequent step operates
    # on this exact list — new files arriving mid-compaction are invisible.
    snapshot = _snapshot_uncompacted_files(state_dir)
    file_paths = [s.path for s in snapshot]

    # Build the stat lookup for conditional deletion
    snapshot_stats = {s.path: (s.mtime_ns, s.size) for s in snapshot}

    # Inventory from the frozen file list
    inventory = _inventory_uncompacted(file_paths, month)

    if not inventory:
        return CompactionResult(
            dry_run=dry_run,
            months=[],
            total_observations=0,
            total_source_files=0,
        )

    if dry_run:
        return CompactionResult(
            dry_run=True,
            months=inventory,
            total_observations=sum(m.observation_count for m in inventory),
            total_source_files=sum(m.source_files for m in inventory),
        )

    # Acquire exclusive lock so two concurrent compactions cannot clobber
    # each other's temp files or double-delete source files.
    with _compaction_lock(hits_dir):
        errors: list[str] = []
        compacted_months: list[MonthCompactionInfo] = []
        files_to_delete: list[Path] = []

        # Build the DuckDB file list once from the frozen snapshot
        file_list_sql = "[" + ",".join(f"'{p}'" for p in file_paths) + "]"

        con = duckdb.connect()

        for info in inventory:
            target_month = info.month
            month_dir = hits_dir / f"month={target_month}"
            month_dir.mkdir(parents=True, exist_ok=True)

            existing_file = month_dir / "data.parquet"
            temp_file = month_dir / ".data.parquet.tmp"

            try:
                # Build source list: uncompacted data for this month
                # Plus existing compacted file if it exists.
                # Uses the frozen file list (not a glob) so concurrent writers
                # cannot inject files into the read set.
                if existing_file.exists():
                    # Use SELECT * with UNION ALL BY NAME for schema evolution:
                    # - Old compacted files may lack new columns (e.g., adjusted_taxid)
                    # - New uncompacted files may have columns old compacted files lack
                    # UNION ALL BY NAME fills missing columns with NULL
                    source_sql = f"""
                        SELECT *
                        FROM read_parquet({file_list_sql}, union_by_name=true)
                        WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'

                        UNION ALL BY NAME

                        SELECT *
                        FROM read_parquet('{existing_file}')
                    """  # SQL uses module constants, not user input
                else:
                    source_sql = f"""
                        SELECT *
                        FROM read_parquet({file_list_sql}, union_by_name=true)
                        WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'
                    """  # SQL uses module constants, not user input

                # Read into Polars, normalize schema, then write.
                # normalize_schema() handles column renames (contig_id ->
                # sequence_id) and adds missing columns (source) with defaults,
                # ensuring compacted output always uses the canonical schema.
                raw_df = con.execute(
                    f"{source_sql} ORDER BY hit_key, run_date",
                ).pl()
                normalized_df = normalize_schema(raw_df)
                normalized_df.write_parquet(
                    temp_file,
                    compression="zstd",
                    row_group_size=100_000,
                )

                # Atomic rename
                temp_file.rename(existing_file)

                # Verify compaction succeeded
                assert existing_file.exists(), (
                    f"compacted file should exist: {existing_file}"
                )
                assert not temp_file.exists(), (
                    f"temp file should be cleaned up: {temp_file}"
                )

                compacted_months.append(info)

                # Track files to delete (if not keeping source).
                # Queries the frozen file list, not a fresh glob.
                if not keep_source:
                    file_result = con.execute(
                        f"""
                        SELECT DISTINCT filename
                        FROM read_parquet({file_list_sql}, filename=true)
                        WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'
                    """,  # SQL uses module constants, not user input
                    ).fetchall()

                    for (filepath,) in file_result:
                        files_to_delete.append(Path(filepath))

            except Exception as e:  # noqa: BLE001  # intentional catch-all to collect errors per month
                errors.append(
                    f"Failed to compact month {target_month}: {type(e).__name__}: {e}",
                )
                # Clean up temp file if it exists
                if temp_file.exists():
                    temp_file.unlink()

        # Delete source files after all compaction is done.
        # Only delete files whose (mtime_ns, size) still match the snapshot —
        # if a concurrent writer replaced a file at the same path, skip it.
        if not keep_source and not errors:
            deleted_dirs: set[Path] = set()
            for filepath in files_to_delete:
                expected = snapshot_stats.get(filepath)
                if expected is None:
                    continue
                try:
                    st = filepath.stat()
                except FileNotFoundError:
                    continue
                if (st.st_mtime_ns, st.st_size) == expected:
                    filepath.unlink()
                    deleted_dirs.add(filepath.parent)

            # Clean up empty directories (sample_id dirs, then sample_set_id dirs)
            for dir_path in sorted(
                deleted_dirs,
                key=lambda p: len(p.parts),
                reverse=True,
            ):
                try:
                    if dir_path.exists() and not any(dir_path.iterdir()):
                        dir_path.rmdir()
                        # Also try to remove parent if empty
                        parent = dir_path.parent
                        if parent.exists() and not any(parent.iterdir()):
                            parent.rmdir()
                except OSError:
                    pass  # Directory not empty or other issue, skip

        return CompactionResult(
            dry_run=False,
            months=compacted_months,
            total_observations=sum(m.observation_count for m in compacted_months),
            total_source_files=sum(m.source_files for m in compacted_months),
            errors=errors,
        )
