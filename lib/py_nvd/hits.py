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
"""

from __future__ import annotations

import zlib
from pathlib import Path
from typing import Literal

import blake3

from py_nvd.db import connect
from py_nvd.models import Hit, HitObservation, HitStats, RecurringHit, TimelineBucket

# Hit key length in hex characters (128 bits = 32 hex chars)
HIT_KEY_LENGTH = 32

# Number of bases packed per byte in 2-bit encoding
BASES_PER_BYTE = 4

# Byte size for storing integer counts/positions in compressed format
INT_BYTE_SIZE = 4

# Complement table for reverse complement
# Handles standard bases, N, and IUPAC ambiguity codes
COMPLEMENT = str.maketrans(
    "ACGTacgtNnRYSWKMBDHVryswkmbdhv",
    "TGCAtgcaNnYRSWMKVHDByrswmkvhdb",
)


# --- Sequence Operations ---


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def canonical_sequence(seq: str) -> str:
    """
    Return the canonical form of a sequence.

    The canonical form is the lexicographically smaller of the sequence
    and its reverse complement. This ensures strand-agnostic identity.
    """
    seq = seq.upper()
    revcomp = reverse_complement(seq)
    return min(seq, revcomp)


def compute_hit_key(seq: str) -> str:
    """
    Compute a deterministic hit key for a sequence.

    Returns a 32-character hex string (128 bits) from BLAKE3 hash
    of the canonical sequence.
    """
    canonical = canonical_sequence(seq)
    digest = blake3.blake3(canonical.encode()).hexdigest()
    return digest[:HIT_KEY_LENGTH]


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
    seq = seq.upper()

    # Identify N positions (including any non-ACGT)
    n_positions = [i for i, base in enumerate(seq) if base not in "ACGT"]

    # 2-bit encode (N -> A as placeholder)
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    bits = [mapping.get(base, 0) for base in seq]  # non-ACGT -> 0

    # Pack 4 bases per byte
    result = bytearray()
    for i in range(0, len(bits), BASES_PER_BYTE):
        chunk = bits[i : i + BASES_PER_BYTE]
        while len(chunk) < BASES_PER_BYTE:
            chunk.append(0)  # pad final byte
        byte = (chunk[0] << 6) | (chunk[1] << 4) | (chunk[2] << 2) | chunk[3]
        result.append(byte)

    # Append N count and positions
    result.extend(len(n_positions).to_bytes(INT_BYTE_SIZE, "little"))
    for pos in n_positions:
        result.extend(pos.to_bytes(INT_BYTE_SIZE, "little"))

    return zlib.compress(bytes(result), level=9)


def decompress_sequence(data: bytes, length: int) -> str:
    """
    Decompress a sequence from 2-bit encoding with N position tracking.

    Args:
        data: Compressed sequence blob from compress_sequence()
        length: Original sequence length in bases

    Returns:
        Original DNA sequence
    """
    raw = zlib.decompress(data)

    # Calculate where 2-bit data ends
    two_bit_len = (length + BASES_PER_BYTE - 1) // BASES_PER_BYTE
    two_bit_data = raw[:two_bit_len]

    # Read N positions
    n_count = int.from_bytes(raw[two_bit_len : two_bit_len + INT_BYTE_SIZE], "little")
    n_positions = []
    offset = two_bit_len + INT_BYTE_SIZE
    for _ in range(n_count):
        n_positions.append(
            int.from_bytes(raw[offset : offset + INT_BYTE_SIZE], "little"),
        )
        offset += INT_BYTE_SIZE

    # Decode 2-bit to bases
    mapping = ["A", "C", "G", "T"]
    result = []
    for byte in two_bit_data:
        result.append(mapping[(byte >> 6) & 0x3])
        result.append(mapping[(byte >> 4) & 0x3])
        result.append(mapping[(byte >> 2) & 0x3])
        result.append(mapping[byte & 0x3])
    result = list("".join(result[:length]))

    # Restore N positions
    for pos in n_positions:
        result[pos] = "N"

    return "".join(result)


def calculate_gc_content(seq: str) -> float:
    """Calculate GC content as a fraction (0.0 to 1.0)."""
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    total = len(seq)
    return gc_count / total if total > 0 else 0.0


# --- Hit Registration ---


def register_hit(
    seq: str,
    first_seen_date: str,
    state_dir: Path | str | None = None,
) -> tuple[Hit, bool]:
    """
    Register a hit sequence in the database.

    Computes the canonical sequence and hit key, then inserts into the hits
    table if not already present. This is idempotent: registering the same
    sequence multiple times returns the existing hit.

    Args:
        seq: The DNA sequence to register
        first_seen_date: ISO8601 date string for when this hit was first seen
        state_dir: Optional state directory override

    Returns:
        Tuple of (Hit, is_new) where is_new is True if this was a new hit
    """
    hit_key = compute_hit_key(seq)

    with connect(state_dir) as conn:
        # Check if hit already exists
        row = conn.execute(
            "SELECT * FROM hits WHERE hit_key = ?",
            (hit_key,),
        ).fetchone()

        if row is not None:
            # Existing hit
            return (
                Hit(
                    hit_key=row["hit_key"],
                    sequence_length=row["sequence_length"],
                    sequence_compressed=row["sequence_compressed"],
                    gc_content=row["gc_content"],
                    first_seen_date=row["first_seen_date"],
                ),
                False,
            )

        # New hit - insert
        compressed = compress_sequence(seq)
        gc_content = calculate_gc_content(seq)

        conn.execute(
            """
            INSERT INTO hits (hit_key, sequence_length, sequence_compressed, gc_content, first_seen_date)
            VALUES (?, ?, ?, ?, ?)
            """,
            (hit_key, len(seq), compressed, gc_content, first_seen_date),
        )
        conn.commit()

        return (
            Hit(
                hit_key=hit_key,
                sequence_length=len(seq),
                sequence_compressed=compressed,
                gc_content=gc_content,
                first_seen_date=first_seen_date,
            ),
            True,
        )


def record_hit_observation(
    hit_key: str,
    sample_set_id: str,
    sample_id: str,
    run_date: str,
    contig_id: str | None = None,
    state_dir: Path | str | None = None,
) -> HitObservation | None:
    """
    Record an observation of a hit in a specific sample/run.

    This is idempotent: recording the same observation multiple times
    (same hit_key, sample_set_id, sample_id) is a no-op.

    Args:
        hit_key: The hit key (from compute_hit_key or register_hit)
        sample_set_id: The sample set this observation belongs to
        sample_id: The sample where this hit was observed
        run_date: ISO8601 date string for when this observation was made
        contig_id: Optional original contig ID from SPAdes (for traceability)
        state_dir: Optional state directory override

    Returns:
        The HitObservation if newly recorded, None if it was a duplicate
    """
    with connect(state_dir) as conn:
        cursor = conn.execute(
            """
            INSERT OR IGNORE INTO hit_observations 
                (hit_key, sample_set_id, sample_id, run_date, contig_id)
            VALUES (?, ?, ?, ?, ?)
            """,
            (hit_key, sample_set_id, sample_id, run_date, contig_id),
        )
        conn.commit()

        if cursor.rowcount == 0:
            # Duplicate - already existed
            return None

        # Return the newly created observation
        return HitObservation(
            hit_key=hit_key,
            sample_set_id=sample_set_id,
            sample_id=sample_id,
            run_date=run_date,
            contig_id=contig_id,
        )


# --- Hit Queries ---


def get_hit(hit_key: str, state_dir: Path | str | None = None) -> Hit | None:
    """Get a hit by its key."""
    with connect(state_dir) as conn:
        row = conn.execute(
            "SELECT * FROM hits WHERE hit_key = ?",
            (hit_key,),
        ).fetchone()

        if row is None:
            return None

        return Hit(
            hit_key=row["hit_key"],
            sequence_length=row["sequence_length"],
            sequence_compressed=row["sequence_compressed"],
            gc_content=row["gc_content"],
            first_seen_date=row["first_seen_date"],
        )


def list_hits(
    limit: int | None = None,
    state_dir: Path | str | None = None,
) -> list[Hit]:
    """
    List all hits, ordered by first_seen_date descending.

    Args:
        limit: Maximum number of hits to return
        state_dir: Optional state directory override

    Returns:
        List of Hit objects
    """
    with connect(state_dir) as conn:
        query = "SELECT * FROM hits ORDER BY first_seen_date DESC"
        params: list[int] = []

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = conn.execute(query, params).fetchall()

        return [
            Hit(
                hit_key=row["hit_key"],
                sequence_length=row["sequence_length"],
                sequence_compressed=row["sequence_compressed"],
                gc_content=row["gc_content"],
                first_seen_date=row["first_seen_date"],
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
    with connect(state_dir) as conn:
        query = "SELECT * FROM hit_observations"
        conditions: list[str] = []
        params: list[str | int] = []

        if hit_key is not None:
            conditions.append("hit_key = ?")
            params.append(hit_key)

        if sample_id is not None:
            conditions.append("sample_id = ?")
            params.append(sample_id)

        if sample_set_id is not None:
            conditions.append("sample_set_id = ?")
            params.append(sample_set_id)

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        query += " ORDER BY run_date DESC"

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = conn.execute(query, params).fetchall()

        return [
            HitObservation(
                hit_key=row["hit_key"],
                sample_set_id=row["sample_set_id"],
                sample_id=row["sample_id"],
                run_date=row["run_date"],
                contig_id=row["contig_id"],
            )
            for row in rows
        ]


def get_hit_sequence(hit: Hit) -> str:
    """
    Decompress and return the sequence for a hit.

    This is a convenience function that handles decompression.

    Args:
        hit: The Hit object (must have sequence_compressed and sequence_length)

    Returns:
        The decompressed DNA sequence
    """
    return decompress_sequence(hit.sequence_compressed, hit.sequence_length)


def count_hits(state_dir: Path | str | None = None) -> int:
    """Return the total number of unique hits in the database."""
    with connect(state_dir) as conn:
        return conn.execute("SELECT COUNT(*) FROM hits").fetchone()[0]


def count_hit_observations(state_dir: Path | str | None = None) -> int:
    """Return the total number of hit observations in the database."""
    with connect(state_dir) as conn:
        return conn.execute("SELECT COUNT(*) FROM hit_observations").fetchone()[0]


# --- Joined Queries ---


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
    with connect(state_dir) as conn:
        query = """
            SELECT 
                h.hit_key, h.sequence_length, h.sequence_compressed, 
                h.gc_content, h.first_seen_date,
                o.sample_set_id, o.sample_id, o.run_date, o.contig_id
            FROM hits h
            INNER JOIN hit_observations o ON h.hit_key = o.hit_key
            ORDER BY h.first_seen_date DESC, o.sample_id
        """
        params: list[int] = []

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = conn.execute(query, params).fetchall()

        return [
            (
                Hit(
                    hit_key=row["hit_key"],
                    sequence_length=row["sequence_length"],
                    sequence_compressed=row["sequence_compressed"],
                    gc_content=row["gc_content"],
                    first_seen_date=row["first_seen_date"],
                ),
                HitObservation(
                    hit_key=row["hit_key"],
                    sample_set_id=row["sample_set_id"],
                    sample_id=row["sample_id"],
                    run_date=row["run_date"],
                    contig_id=row["contig_id"],
                ),
            )
            for row in rows
        ]


# --- Statistics / Aggregations ---


def get_hit_stats(state_dir: Path | str | None = None) -> HitStats:
    """
    Compute summary statistics for all hits.

    Returns aggregate metrics including counts, length/GC distributions,
    and date range. Useful for dashboard displays and health checks.

    Args:
        state_dir: Optional state directory override

    Returns:
        HitStats with all computed metrics. Counts are 0 and distributions
        are None when the database is empty.
    """
    with connect(state_dir) as conn:
        # Basic counts
        total_hits = conn.execute("SELECT COUNT(*) FROM hits").fetchone()[0]
        total_observations = conn.execute(
            "SELECT COUNT(*) FROM hit_observations"
        ).fetchone()[0]
        unique_samples = conn.execute(
            "SELECT COUNT(DISTINCT sample_id) FROM hit_observations"
        ).fetchone()[0]
        unique_runs = conn.execute(
            "SELECT COUNT(DISTINCT sample_set_id) FROM hit_observations"
        ).fetchone()[0]

        # If no hits, return early with empty stats
        if total_hits == 0:
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

        # Length stats (min, max)
        length_row = conn.execute(
            "SELECT MIN(sequence_length), MAX(sequence_length) FROM hits"
        ).fetchone()
        length_min, length_max = length_row[0], length_row[1]

        # GC stats (min, max)
        gc_row = conn.execute(
            "SELECT MIN(gc_content), MAX(gc_content) FROM hits"
        ).fetchone()
        gc_min, gc_max = gc_row[0], gc_row[1]

        # Date range
        date_row = conn.execute(
            "SELECT MIN(first_seen_date), MAX(first_seen_date) FROM hits"
        ).fetchone()
        date_first, date_last = date_row[0], date_row[1]

        # Median calculations require fetching sorted values
        # SQLite doesn't have a built-in MEDIAN function
        lengths = [
            row[0]
            for row in conn.execute(
                "SELECT sequence_length FROM hits ORDER BY sequence_length"
            ).fetchall()
        ]
        length_median = _compute_median(lengths)

        gc_values = [
            row[0]
            for row in conn.execute(
                "SELECT gc_content FROM hits ORDER BY gc_content"
            ).fetchall()
        ]
        gc_median = _compute_median(gc_values)

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


def _compute_median(values: list[int | float]) -> float | None:
    """Compute median of a sorted list of values."""
    if not values:
        return None
    n = len(values)
    mid = n // 2
    if n % 2 == 0:
        return (values[mid - 1] + values[mid]) / 2
    return float(values[mid])


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
    with connect(state_dir) as conn:
        # Build query with dynamic HAVING clause
        query = """
            SELECT
                h.hit_key,
                h.sequence_length,
                h.gc_content,
                h.first_seen_date,
                MAX(o.run_date) as last_seen_date,
                COUNT(DISTINCT o.sample_id) as sample_count,
                COUNT(DISTINCT o.sample_set_id) as run_count,
                COUNT(*) as observation_count
            FROM hits h
            INNER JOIN hit_observations o ON h.hit_key = o.hit_key
            GROUP BY h.hit_key
            HAVING sample_count >= ?
        """
        params: list[int] = [min_samples]

        if min_runs is not None:
            query += " AND run_count >= ?"
            params.append(min_runs)

        query += " ORDER BY sample_count DESC, observation_count DESC"

        if limit is not None:
            query += " LIMIT ?"
            params.append(limit)

        rows = conn.execute(query, params).fetchall()

        return [
            RecurringHit(
                hit_key=row["hit_key"],
                sequence_length=row["sequence_length"],
                gc_content=row["gc_content"],
                first_seen_date=row["first_seen_date"],
                last_seen_date=row["last_seen_date"],
                sample_count=row["sample_count"],
                run_count=row["run_count"],
                observation_count=row["observation_count"],
            )
            for row in rows
        ]


def get_discovery_timeline(
    granularity: Literal["day", "week", "month", "year"] = "month",
    state_dir: Path | str | None = None,
) -> list[TimelineBucket]:
    """
    Count new hits discovered per time period.

    Groups hits by their first_seen_date and counts how many were discovered
    in each period. Useful for visualizing discovery rate over time.

    Args:
        granularity: Time period grouping - "day", "week", "month", or "year"
        state_dir: Optional state directory override

    Returns:
        List of TimelineBucket objects ordered by period ascending.
        Empty list if no hits exist.
    """
    # SQLite strftime format strings for each granularity
    format_map = {
        "day": "%Y-%m-%d",
        "week": "%Y-W%W",  # ISO week number
        "month": "%Y-%m",
        "year": "%Y",
    }
    date_format = format_map[granularity]

    with connect(state_dir) as conn:
        query = f"""
            SELECT
                strftime('{date_format}', first_seen_date) as period,
                COUNT(*) as new_hits
            FROM hits
            GROUP BY period
            ORDER BY period ASC
        """

        rows = conn.execute(query).fetchall()

        return [
            TimelineBucket(
                period=row["period"],
                new_hits=row["new_hits"],
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

    This is a convenience function for CLI and interactive use.

    Args:
        query: Either a 32-character hex hit key, or a DNA sequence
        state_dir: Optional state directory override

    Returns:
        The Hit if found, None otherwise.
    """
    query = query.strip()

    # Check if it looks like a hex key (32 chars, all hex digits)
    if len(query) == HIT_KEY_LENGTH and all(
        c in "0123456789abcdefABCDEF" for c in query
    ):
        return get_hit(query.lower(), state_dir)

    # Otherwise treat as sequence, compute key and look up
    hit_key = compute_hit_key(query)
    return get_hit(hit_key, state_dir)


# --- Cleanup Operations ---


def delete_observations_for_sample_set(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> int:
    """
    Delete all hit observations for a sample set.

    Use this when pruning a run to remove its associated observations.
    After deleting observations, call delete_orphaned_hits() to clean up
    any hits that no longer have observations.

    Args:
        sample_set_id: The sample set whose observations to delete
        state_dir: Optional state directory override

    Returns:
        Number of observations deleted
    """
    with connect(state_dir) as conn:
        count = conn.execute(
            "DELETE FROM hit_observations WHERE sample_set_id = ?",
            (sample_set_id,),
        ).rowcount
        conn.commit()
        return count


def delete_orphaned_hits(state_dir: Path | str | None = None) -> int:
    """
    Delete hits that have no observations.

    A hit without observations has no provenance (no sample, no run context)
    and is not useful. Call this after deleting observations to garbage
    collect orphaned hits.

    Args:
        state_dir: Optional state directory override

    Returns:
        Number of orphaned hits deleted
    """
    with connect(state_dir) as conn:
        count = conn.execute(
            """
            DELETE FROM hits
            WHERE hit_key NOT IN (SELECT DISTINCT hit_key FROM hit_observations)
            """
        ).rowcount
        conn.commit()
        return count


def prune_hits_for_sample_set(
    sample_set_id: str,
    state_dir: Path | str | None = None,
) -> tuple[int, int]:
    """
    Delete observations for a sample set and garbage collect orphaned hits.

    This is the recommended way to clean up hits when pruning a run. It
    performs both operations atomically in a single transaction.

    Args:
        sample_set_id: The sample set whose observations to delete
        state_dir: Optional state directory override

    Returns:
        Tuple of (observations_deleted, orphaned_hits_deleted)
    """
    with connect(state_dir) as conn:
        # Delete observations for this sample set
        observations_deleted = conn.execute(
            "DELETE FROM hit_observations WHERE sample_set_id = ?",
            (sample_set_id,),
        ).rowcount

        # Garbage collect orphaned hits
        orphaned_hits_deleted = conn.execute(
            """
            DELETE FROM hits
            WHERE hit_key NOT IN (SELECT DISTINCT hit_key FROM hit_observations)
            """
        ).rowcount

        conn.commit()
        return observations_deleted, orphaned_hits_deleted
