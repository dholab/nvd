"""
Hit management for metagenomic BLAST results.

Provides idempotent hit registration and querying. A "hit" is a contig sequence
that had BLAST results. Hit keys are computed from canonical sequences
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

import shutil
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import blake3
import duckdb
import polars as pl

from py_nvd.db import get_hits_dir
from py_nvd.models import (
    CompactionResult,
    DeleteResult,
    Hit,
    HitObservation,
    HitStats,
    MonthCompactionInfo,
    RecurringHit,
    TimelineBucket,
)

HIT_KEY_LENGTH = 32
BASES_PER_BYTE = 4
INT_BYTE_SIZE = 4

COMPLEMENT = str.maketrans(
    "ACGTacgtNnRYSWKMBDHVryswkmbdhv",
    "TGCAtgcaNnYRSWMKVHDByrswmkvhdb",
)


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
            *,
            COALESCE(month, strftime(run_date::DATE, '%Y-%m')) AS effective_month
        FROM read_parquet('{glob_pattern}', hive_partitioning=true, union_by_name=true)
    """


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
    """

    hit_key: str
    sequence_length: int
    sequence_compressed: bytes
    gc_content: float
    sample_set_id: str
    sample_id: str
    run_date: str
    contig_id: str | None

    def __post_init__(self) -> None:
        """Validate invariants after initialization."""
        assert len(self.hit_key) == 32, (
            f"hit_key must be 32 hex chars, got {len(self.hit_key)}"
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
        # Write empty parquet with correct schema
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
        # Build DataFrame from hit records
        df = pl.DataFrame(
            {
                "hit_key": [h.hit_key for h in hits],
                "sequence_length": [h.sequence_length for h in hits],
                "sequence_compressed": [h.sequence_compressed for h in hits],
                "gc_content": [h.gc_content for h in hits],
                "sample_set_id": [h.sample_set_id for h in hits],
                "sample_id": [h.sample_id for h in hits],
                "run_date": [h.run_date for h in hits],
                "contig_id": [h.contig_id for h in hits],
            }
        )

    # Write to temp file first
    df.write_parquet(temp_path)

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
    assert isinstance(count, int) and count >= 0, (
        f"count must be non-negative int, got {count}"
    )
    return count


def count_hit_observations(state_dir: Path | str | None = None) -> int:
    """Return the total number of hit observations across all parquet files."""
    if not has_any_hits(state_dir):
        return 0

    con = query_hits(state_dir)
    result = con.execute("SELECT COUNT(*) FROM hits").fetchone()
    assert result is not None, "COUNT query should always return a row"
    count = result[0]
    assert isinstance(count, int) and count >= 0, (
        f"count must be non-negative int, got {count}"
    )
    return count


def get_hit(hit_key: str, state_dir: Path | str | None = None) -> Hit | None:
    """
    Get a hit by its key.

    Since parquet stores denormalized data, this aggregates across all
    observations to compute first_seen_date as MIN(run_date).
    """
    assert hit_key, "hit_key cannot be empty"
    assert len(hit_key) == 32, f"hit_key must be 32 hex chars, got {len(hit_key)}"

    if not has_any_hits(state_dir):
        return None

    con = query_hits(state_dir)
    result = con.execute(
        """
        SELECT
            hit_key,
            sequence_length,
            sequence_compressed,
            gc_content,
            MIN(run_date) as first_seen_date
        FROM hits
        WHERE hit_key = ?
        GROUP BY hit_key, sequence_length, sequence_compressed, gc_content
        """,
        [hit_key.lower()],
    ).fetchone()

    if result is None:
        return None

    return Hit(
        hit_key=result[0],
        sequence_length=result[1],
        sequence_compressed=result[2],
        gc_content=result[3],
        first_seen_date=result[4],
    )


def list_hits(
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[Hit]:
    """
    List all unique hits, ordered by first_seen_date descending.

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

    sql = """
        SELECT
            hit_key,
            sequence_length,
            sequence_compressed,
            gc_content,
            MIN(run_date) as first_seen_date
        FROM hits
        GROUP BY hit_key, sequence_length, sequence_compressed, gc_content
        ORDER BY first_seen_date DESC
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

    sql = "SELECT hit_key, sample_set_id, sample_id, run_date, contig_id FROM hits"
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
            contig_id=row[4],
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
        "SELECT COUNT(DISTINCT sample_id) FROM hits"
    ).fetchone()
    assert samples_result is not None, "COUNT query should always return a row"
    unique_samples = samples_result[0]

    runs_result = con.execute(
        "SELECT COUNT(DISTINCT sample_set_id) FROM hits"
    ).fetchone()
    assert runs_result is not None, "COUNT query should always return a row"
    unique_runs = runs_result[0]

    assert isinstance(total_hits, int) and total_hits >= 0, (
        "total_hits must be non-negative"
    )
    assert isinstance(total_observations, int) and total_observations >= 0, (
        "total_observations must be non-negative"
    )

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

    # Build query - note we compute first_seen_date as MIN(run_date)
    sql = """
        SELECT
            hit_key,
            ANY_VALUE(sequence_length) as sequence_length,
            ANY_VALUE(gc_content) as gc_content,
            MIN(run_date) as first_seen_date,
            MAX(run_date) as last_seen_date,
            COUNT(DISTINCT sample_id) as sample_count,
            COUNT(DISTINCT sample_set_id) as run_count,
            COUNT(*) as observation_count
        FROM hits
        GROUP BY hit_key
        HAVING sample_count >= ?
    """
    params: list[int] = [min_samples]

    if min_runs is not None:
        sql += " AND run_count >= ?"
        params.append(min_runs)

    sql += " ORDER BY sample_count DESC, observation_count DESC"

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
    """

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
            h.contig_id
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
        )
        observation = HitObservation(
            hit_key=row[0],
            sample_set_id=row[5],
            sample_id=row[6],
            run_date=row[7],
            contig_id=row[8],
        )
        result.append((hit, observation))

    assert all(isinstance(t, tuple) and len(t) == 2 for t in result), (
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
    assert isinstance(count, int) and count >= 0, (
        f"count must be non-negative int, got {count}"
    )
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
        has_data = con.execute(
            f"""
            SELECT COUNT(*) > 0
            FROM read_parquet('{data_file}')
            WHERE sample_set_id = $1
            """,
            [sample_set_id],
        ).fetchone()[0]

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
            """,
            [sample_set_id],
        )

        row_count = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{temp_file}')"
        ).fetchone()[0]

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


def _get_uncompacted_glob(state_dir: Path | str | None = None) -> str:
    """Get glob pattern for uncompacted parquet files."""
    hits_dir = get_hits_dir(state_dir)
    result = str(hits_dir / "month=NULL" / "**" / "*.parquet")
    assert "month=NULL" in result, "glob pattern must target uncompacted partition"
    return result


def _inventory_uncompacted(
    state_dir: Path | str | None = None,
    month: str | None = None,
) -> list[MonthCompactionInfo]:
    """
    Inventory uncompacted data grouped by month.

    Args:
        state_dir: Optional state directory override
        month: Optional filter to a specific month (YYYY-MM format)

    Returns:
        List of MonthCompactionInfo for each month with uncompacted data
    """
    assert month is None or (len(month) == 7 and month[4] == "-"), (
        f"month must be YYYY-MM format, got {month!r}"
    )

    uncompacted_glob = _get_uncompacted_glob(state_dir)
    hits_dir = get_hits_dir(state_dir)
    uncompacted_dir = hits_dir / "month=NULL"

    if not uncompacted_dir.exists():
        return []

    parquet_files = list(uncompacted_dir.rglob("*.parquet"))
    if not parquet_files:
        return []

    con = duckdb.connect()

    month_filter = ""
    if month is not None:
        month_filter = f"WHERE strftime(run_date::DATE, '%Y-%m') = '{month}'"

    result = con.execute(f"""
        SELECT
            strftime(run_date::DATE, '%Y-%m') AS target_month,
            COUNT(*) AS observation_count,
            COUNT(DISTINCT hit_key) AS unique_hits,
            COUNT(DISTINCT sample_set_id) AS sample_sets
        FROM read_parquet('{uncompacted_glob}')
        {month_filter}
        GROUP BY target_month
        ORDER BY target_month DESC
    """).fetchall()

    # Count source files per month by examining the data
    inventory = []
    for row in result:
        target_month, obs_count, unique_hits, sample_sets = row

        # Count files that have data for this month
        file_count = con.execute(f"""
            SELECT COUNT(DISTINCT filename)
            FROM read_parquet('{uncompacted_glob}', filename=true)
            WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'
        """).fetchone()

        inventory.append(
            MonthCompactionInfo(
                month=target_month,
                observation_count=obs_count,
                unique_hits=unique_hits,
                sample_sets=sample_sets,
                source_files=file_count[0] if file_count else 0,
            )
        )

    assert all(isinstance(m, MonthCompactionInfo) for m in inventory), (
        "inventory must contain only MonthCompactionInfo objects"
    )
    return inventory


def compact_hits(
    state_dir: Path | str | None = None,
    month: str | None = None,
    dry_run: bool = False,
    keep_source: bool = False,
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
    assert month is None or len(month) == 7, (
        f"month must be YYYY-MM format, got {month!r}"
    )

    hits_dir = get_hits_dir(state_dir)
    uncompacted_glob = _get_uncompacted_glob(state_dir)

    # Inventory what we have to compact
    inventory = _inventory_uncompacted(state_dir, month)

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

    # Perform compaction for each month
    errors: list[str] = []
    compacted_months: list[MonthCompactionInfo] = []
    files_to_delete: list[Path] = []

    con = duckdb.connect()

    for info in inventory:
        target_month = info.month
        month_dir = hits_dir / f"month={target_month}"
        month_dir.mkdir(parents=True, exist_ok=True)

        existing_file = month_dir / "data.parquet"
        temp_file = month_dir / ".data.parquet.tmp"

        try:
            # Build source list: uncompacted data for this month
            # Plus existing compacted file if it exists
            if existing_file.exists():
                source_sql = f"""
                    SELECT
                        hit_key,
                        sequence_length,
                        sequence_compressed,
                        gc_content,
                        sample_set_id,
                        sample_id,
                        run_date,
                        contig_id
                    FROM read_parquet('{uncompacted_glob}')
                    WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'

                    UNION ALL

                    SELECT
                        hit_key,
                        sequence_length,
                        sequence_compressed,
                        gc_content,
                        sample_set_id,
                        sample_id,
                        run_date,
                        contig_id
                    FROM read_parquet('{existing_file}')
                """
            else:
                source_sql = f"""
                    SELECT
                        hit_key,
                        sequence_length,
                        sequence_compressed,
                        gc_content,
                        sample_set_id,
                        sample_id,
                        run_date,
                        contig_id
                    FROM read_parquet('{uncompacted_glob}')
                    WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'
                """

            # Write sorted, compressed parquet
            con.execute(f"""
                COPY (
                    {source_sql}
                    ORDER BY hit_key, run_date
                )
                TO '{temp_file}'
                (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000)
            """)

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

            # Track files to delete (if not keeping source)
            if not keep_source:
                # Find all uncompacted files that contributed to this month
                file_result = con.execute(f"""
                    SELECT DISTINCT filename
                    FROM read_parquet('{uncompacted_glob}', filename=true)
                    WHERE strftime(run_date::DATE, '%Y-%m') = '{target_month}'
                """).fetchall()

                for (filepath,) in file_result:
                    files_to_delete.append(Path(filepath))

        except Exception as e:
            errors.append(
                f"Failed to compact month {target_month}: {type(e).__name__}: {e}"
            )
            # Clean up temp file if it exists
            if temp_file.exists():
                temp_file.unlink()

    # Delete source files after all compaction is done
    if not keep_source and not errors:
        deleted_dirs: set[Path] = set()
        for filepath in files_to_delete:
            if filepath.exists():
                filepath.unlink()
                # Track parent directories for cleanup
                deleted_dirs.add(filepath.parent)

        # Clean up empty directories (sample_id dirs, then sample_set_id dirs)
        for dir_path in sorted(deleted_dirs, key=lambda p: len(p.parts), reverse=True):
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
