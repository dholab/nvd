"""Tests for py_nvd.hits module.

These tests verify:
- Sequence operations (reverse complement, canonical form)
- Hit key computation (deterministic, strand-agnostic)
- Sequence compression/decompression (round-trip, N handling)
- GC content calculation
- Parquet-based hit queries
"""

import pytest

from py_nvd.hits import (
    HitRecord,
    SampleSetStats,
    calculate_gc_content,
    canonical_sequence,
    compact_hits,
    compress_sequence,
    compute_hit_key,
    count_hit_observations,
    count_hits,
    count_sample_set_observations,
    decompress_sequence,
    delete_sample_hits,
    delete_sample_set_hits,
    get_discovery_timeline,
    get_hit,
    get_hit_sequence,
    get_hit_stats,
    get_recurring_hits,
    get_sample_parquet_path,
    get_stats_for_sample_set,
    is_valid_hit_key,
    list_hit_observations,
    list_hits,
    list_hits_with_observations,
    lookup_hit,
    query_hits,
    reverse_complement,
    sample_hits_exist,
    write_hits_parquet,
)
from py_nvd.models import DeleteResult


def make_hit_record(
    seq: str,
    sample_set_id: str,
    sample_id: str,
    run_date: str,
    contig_id: str | None = None,
) -> HitRecord:
    """Helper to create a HitRecord from a sequence."""
    return HitRecord(
        hit_key=compute_hit_key(seq),
        sequence_length=len(seq),
        sequence_compressed=compress_sequence(seq),
        gc_content=calculate_gc_content(seq),
        sample_set_id=sample_set_id,
        sample_id=sample_id,
        run_date=run_date,
        contig_id=contig_id,
    )


class TestReverseComplement:
    def test_palindrome(self):
        assert reverse_complement("ACGT") == "ACGT"

    def test_asymmetric(self):
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("CCCC") == "GGGG"
        assert reverse_complement("AAACCC") == "GGGTTT"

    def test_palindrome_longer(self):
        assert reverse_complement("GCGC") == "GCGC"

    def test_with_n(self):
        assert reverse_complement("ACNGT") == "ACNGT"
        assert reverse_complement("NNNN") == "NNNN"

    def test_lowercase(self):
        assert reverse_complement("acgt") == "acgt"
        assert reverse_complement("aaaa") == "tttt"

    def test_mixed_case(self):
        assert reverse_complement("AcGt") == "aCgT"


class TestCanonicalSequence:
    def test_already_canonical(self):
        assert canonical_sequence("AAAA") == "AAAA"

    def test_needs_revcomp(self):
        assert canonical_sequence("TTTT") == "AAAA"

    def test_palindrome(self):
        assert canonical_sequence("ACGT") == "ACGT"

    def test_lowercase_normalized(self):
        assert canonical_sequence("aaaa") == "AAAA"
        assert canonical_sequence("tttt") == "AAAA"

    def test_complex_sequence(self):
        seq = "AAACCCGGG"
        assert canonical_sequence(seq) == "AAACCCGGG"

    def test_complex_sequence_reverse(self):
        seq = "CCCGGGTTT"
        assert canonical_sequence(seq) == "AAACCCGGG"


class TestComputeHitKey:
    def test_deterministic(self):
        seq = "ACGTACGTACGT"
        key1 = compute_hit_key(seq)
        key2 = compute_hit_key(seq)
        assert key1 == key2

    def test_length(self):
        key = compute_hit_key("ACGT")
        assert len(key) == 32

    def test_hex_format(self):
        key = compute_hit_key("ACGT")
        int(key, 16)

    def test_strand_agnostic(self):
        seq = "AAACCCGGG"
        revcomp = reverse_complement(seq)
        assert compute_hit_key(seq) == compute_hit_key(revcomp)

    def test_different_sequences_different_keys(self):
        key1 = compute_hit_key("AAAA")
        key2 = compute_hit_key("CCCC")
        assert key1 != key2

    def test_case_insensitive(self):
        key1 = compute_hit_key("ACGT")
        key2 = compute_hit_key("acgt")
        assert key1 == key2


class TestCompression:
    def test_round_trip_simple(self):
        seq = "ACGTACGTACGT"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == seq

    def test_round_trip_all_bases(self):
        seq = "AAAACCCCGGGGTTTT"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == seq

    def test_round_trip_with_n(self):
        seq = "ACGTNNNACGT"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == seq

    def test_round_trip_n_at_start(self):
        seq = "NNNACGT"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == seq

    def test_round_trip_n_at_end(self):
        seq = "ACGTNNN"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == seq

    def test_round_trip_all_n(self):
        seq = "NNNNNNNN"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == seq

    def test_round_trip_long_sequence(self):
        seq = "ACGT" * 1000 + "NNN" + "TGCA" * 500
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == seq

    def test_ambiguous_bases_become_n(self):
        seq = "ACRYSWGT"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == "ACNNNNGT"

    def test_all_iupac_ambiguous(self):
        seq = "RYSWKMBDHV"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == "NNNNNNNNNN"

    def test_compression_ratio(self):
        seq = "ACGT" * 1000
        compressed = compress_sequence(seq)
        assert len(compressed) < len(seq) * 0.30

    def test_lowercase_normalized(self):
        seq = "acgtacgt"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == "ACGTACGT"

    def test_odd_length(self):
        for length in [1, 2, 3, 5, 7, 13, 17]:
            seq = (
                "ACGT"[:length]
                if length <= 4
                else "ACGT" * (length // 4) + "ACG"[: length % 4]
            )
            seq = seq[:length]
            compressed = compress_sequence(seq)
            decompressed = decompress_sequence(compressed, len(seq))
            assert decompressed == seq.upper()


class TestGCContent:
    def test_all_gc(self):
        assert calculate_gc_content("GCGCGC") == 1.0

    def test_all_at(self):
        assert calculate_gc_content("ATATAT") == 0.0

    def test_balanced(self):
        assert calculate_gc_content("ACGT") == 0.5

    def test_empty(self):
        assert calculate_gc_content("") == 0.0

    def test_with_n(self):
        assert calculate_gc_content("GCNN") == 0.5

    def test_lowercase(self):
        assert calculate_gc_content("gcgc") == 1.0

    def test_typical_viral(self):
        seq = "AATTAATTGC"
        assert calculate_gc_content(seq) == 0.2


class TestGetHit:
    """Tests for get_hit()."""

    def test_get_existing_hit(self, temp_state_dir):
        """Getting an existing hit returns it."""
        seq = "ACGTACGTACGT"
        hit_record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([hit_record], "sample_a", "set_001", temp_state_dir)

        retrieved = get_hit(hit_record.hit_key, temp_state_dir)

        assert retrieved is not None
        assert retrieved.hit_key == hit_record.hit_key
        assert retrieved.sequence_length == len(seq)

    def test_get_nonexistent_hit(self, temp_state_dir):
        """Getting a nonexistent hit returns None."""
        result = get_hit("a" * 32, temp_state_dir)
        assert result is None


class TestListHits:
    """Tests for list_hits()."""

    def test_list_empty(self, temp_state_dir):
        """Listing hits from empty database returns empty list."""
        hits = list_hits(state_dir=temp_state_dir)
        assert hits == []

    def test_list_multiple_hits(self, temp_state_dir):
        """Listing hits returns all registered hits."""
        records = [
            make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "set_001", "sample_a", "2024-01-02T00:00:00Z"),
            make_hit_record("AAATTT", "set_001", "sample_a", "2024-01-03T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        hits = list_hits(state_dir=temp_state_dir)

        assert len(hits) == 3

    def test_list_with_limit(self, temp_state_dir):
        """Listing hits respects limit parameter."""
        records = [
            make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "set_001", "sample_a", "2024-01-02T00:00:00Z"),
            make_hit_record("AAATTT", "set_001", "sample_a", "2024-01-03T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        hits = list_hits(limit=2, state_dir=temp_state_dir)

        assert len(hits) == 2

    def test_list_ordered_by_date_desc(self, temp_state_dir):
        """Hits are ordered by first_seen_date descending."""
        records = [
            make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "set_001", "sample_a", "2024-01-03T00:00:00Z"),
            make_hit_record("AAATTT", "set_001", "sample_a", "2024-01-02T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        hits = list_hits(state_dir=temp_state_dir)

        assert hits[0].first_seen_date == "2024-01-03T00:00:00Z"
        assert hits[1].first_seen_date == "2024-01-02T00:00:00Z"
        assert hits[2].first_seen_date == "2024-01-01T00:00:00Z"


class TestListHitObservations:
    """Tests for list_hit_observations()."""

    def test_list_empty(self, temp_state_dir):
        """Listing observations from empty database returns empty list."""
        obs = list_hit_observations(state_dir=temp_state_dir)
        assert obs == []

    def test_list_all_observations(self, temp_state_dir):
        """Listing without filters returns all observations."""
        records = [
            make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAA", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
            make_hit_record("CCCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records[:2], "sample_a", "set_001", temp_state_dir)
        write_hits_parquet([records[2]], "sample_b", "set_001", temp_state_dir)

        obs = list_hit_observations(state_dir=temp_state_dir)

        assert len(obs) == 3

    def test_filter_by_hit_key(self, temp_state_dir):
        """Filtering by hit_key returns only matching observations."""
        hit_key_aaaa = compute_hit_key("AAAA")
        records = [
            make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAA", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
            make_hit_record("CCCC", "set_001", "sample_c", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        obs = list_hit_observations(hit_key=hit_key_aaaa, state_dir=temp_state_dir)

        assert len(obs) == 2
        assert all(o.hit_key == hit_key_aaaa for o in obs)

    def test_filter_by_sample_id(self, temp_state_dir):
        """Filtering by sample_id returns only matching observations."""
        records = [
            make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("CCCC", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        obs = list_hit_observations(sample_id="sample_a", state_dir=temp_state_dir)

        assert len(obs) == 1
        assert obs[0].sample_id == "sample_a"

    def test_filter_by_sample_set_id(self, temp_state_dir):
        """Filtering by sample_set_id returns only matching observations."""
        write_hits_parquet(
            [make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAAA", "set_002", "sample_a", "2024-01-02T00:00:00Z")],
            "sample_a",
            "set_002",
            temp_state_dir,
        )

        obs = list_hit_observations(sample_set_id="set_001", state_dir=temp_state_dir)

        assert len(obs) == 1
        assert obs[0].sample_set_id == "set_001"


class TestCountFunctions:
    """Tests for count_hits() and count_hit_observations()."""

    def test_count_hits_empty(self, temp_state_dir):
        """Counting hits in empty database returns 0."""
        assert count_hits(temp_state_dir) == 0

    def test_count_hits(self, temp_state_dir):
        """Counting hits returns correct count."""
        records = [
            make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAATTT", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        assert count_hits(temp_state_dir) == 3

    def test_count_observations_empty(self, temp_state_dir):
        """Counting observations in empty database returns 0."""
        assert count_hit_observations(temp_state_dir) == 0

    def test_count_observations(self, temp_state_dir):
        """Counting observations returns correct count."""
        records = [
            make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAA", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        assert count_hit_observations(temp_state_dir) == 2


class TestHivePartitioning:
    """Tests for Hive partition column extraction."""

    def test_month_is_null_for_uncompacted(self, temp_state_dir):
        """Uncompacted files have month=NULL extracted from path."""
        import duckdb

        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        hits_glob = temp_state_dir / "hits" / "**" / "*.parquet"
        con = duckdb.connect()
        result = con.execute(f"""
            SELECT month FROM read_parquet('{hits_glob}', hive_partitioning=true)
        """).fetchone()

        assert result is not None
        assert result[0] is None  # month=NULL in path becomes SQL NULL

    def test_effective_month_computed_from_run_date(self, temp_state_dir):
        """effective_month is computed from run_date when month is NULL."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-03-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        con = query_hits(temp_state_dir)
        result = con.execute("SELECT month, effective_month FROM hits").fetchone()

        assert result is not None
        assert result[0] is None  # month is NULL
        assert result[1] == "2024-03"  # effective_month computed from run_date

    def test_directory_structure_is_hive_partitioned(self, temp_state_dir):
        """Written files follow month=NULL/{sample_set_id}/{sample_id}/ structure."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        path = write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        assert path.name == "data.parquet"
        assert path.parent.name == "sample_a"
        assert path.parent.parent.name == "set_001"
        assert path.parent.parent.parent.name == "month=NULL"


class TestGetHitSequence:
    """Tests for get_hit_sequence()."""

    def test_get_sequence_simple(self, temp_state_dir):
        """Getting sequence from hit returns original sequence."""
        seq = "ACGTACGTACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        hit = get_hit(record.hit_key, temp_state_dir)
        assert hit is not None
        recovered = get_hit_sequence(hit)

        assert recovered == seq

    def test_get_sequence_with_n(self, temp_state_dir):
        """Getting sequence preserves N positions."""
        seq = "ACGTNNNACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        hit = get_hit(record.hit_key, temp_state_dir)
        assert hit is not None
        recovered = get_hit_sequence(hit)

        assert recovered == seq

    def test_get_sequence_long(self, temp_state_dir):
        """Getting long sequence works correctly."""
        seq = "ACGT" * 1000 + "NNN" + "TGCA" * 500
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        hit = get_hit(record.hit_key, temp_state_dir)
        assert hit is not None
        recovered = get_hit_sequence(hit)

        assert recovered == seq


class TestListHitsWithObservations:
    """Tests for list_hits_with_observations()."""

    def test_empty_database(self, temp_state_dir):
        """Empty database returns empty list."""
        result = list_hits_with_observations(state_dir=temp_state_dir)
        assert result == []

    def test_returns_joined_data(self, temp_state_dir):
        """Returns hits joined with their observations."""
        record = make_hit_record(
            "AAACCC",
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = list_hits_with_observations(state_dir=temp_state_dir)

        assert len(result) == 1
        returned_hit, returned_obs = result[0]
        assert returned_hit.hit_key == record.hit_key
        assert returned_obs.sample_id == "sample_a"

    def test_one_hit_multiple_observations(self, temp_state_dir):
        """Hit observed in multiple samples appears multiple times."""
        records = [
            make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        result = list_hits_with_observations(state_dir=temp_state_dir)

        assert len(result) == 2
        sample_ids = {obs.sample_id for _, obs in result}
        assert sample_ids == {"sample_a", "sample_b"}

    def test_respects_limit(self, temp_state_dir):
        """Limit parameter restricts number of rows."""
        records = [
            make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "set_001", "sample_c", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        result = list_hits_with_observations(limit=2, state_dir=temp_state_dir)

        assert len(result) == 2


class TestDeleteSampleHits:
    """Tests for delete_sample_hits()."""

    def test_deletes_sample_file(self, temp_state_dir):
        """Deletes the parquet file for a sample."""
        record = make_hit_record(
            "AAACCC",
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        assert count_hits(temp_state_dir) == 1

        deleted = delete_sample_hits("sample_a", "set_001", temp_state_dir)

        assert deleted is True
        assert count_hits(temp_state_dir) == 0

    def test_returns_false_for_nonexistent(self, temp_state_dir):
        """Returns False when file doesn't exist."""
        deleted = delete_sample_hits("nonexistent", "set_001", temp_state_dir)
        assert deleted is False


class TestDeleteSampleSetHits:
    """Tests for delete_sample_set_hits()."""

    def test_deletes_all_files_in_set(self, temp_state_dir):
        """Deletes all parquet files for a sample set."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAAGGG", "set_001", "sample_b", "2024-01-01T00:00:00Z")],
            "sample_b",
            "set_001",
            temp_state_dir,
        )

        assert count_hits(temp_state_dir) == 2

        result = delete_sample_set_hits("set_001", temp_state_dir)

        assert result == DeleteResult(
            uncompacted_files_deleted=2,
            compacted_months_rewritten=0,
        )
        assert count_hits(temp_state_dir) == 0

    def test_only_deletes_matching_set(self, temp_state_dir):
        """Only deletes files for the specified sample set."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-01T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAAGGG", "set_002", "sample_b", "2024-01-01T00:00:00Z")],
            "sample_b",
            "set_002",
            temp_state_dir,
        )

        result = delete_sample_set_hits("set_001", temp_state_dir)

        assert result == DeleteResult(
            uncompacted_files_deleted=1,
            compacted_months_rewritten=0,
        )
        assert count_hits(temp_state_dir) == 1

    def test_returns_zero_for_nonexistent(self, temp_state_dir):
        """Returns zeros when sample set doesn't exist."""
        result = delete_sample_set_hits("nonexistent", temp_state_dir)
        assert result == DeleteResult(
            uncompacted_files_deleted=0,
            compacted_months_rewritten=0,
        )

    def test_deletes_compacted_data(self, temp_state_dir):
        """Rewrites compacted month files to exclude deleted sample set."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-15T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAAGGG", "set_002", "sample_b", "2024-01-20T00:00:00Z")],
            "sample_b",
            "set_002",
            temp_state_dir,
        )

        compact_hits(temp_state_dir)
        assert count_hits(temp_state_dir) == 2

        result = delete_sample_set_hits("set_001", temp_state_dir)

        assert result == DeleteResult(
            uncompacted_files_deleted=0,
            compacted_months_rewritten=1,
        )
        assert count_hits(temp_state_dir) == 1

        observations = list_hit_observations(state_dir=temp_state_dir)
        assert len(observations) == 1
        assert observations[0].sample_set_id == "set_002"

    def test_deletes_compacted_across_multiple_months(self, temp_state_dir):
        """Rewrites multiple month files when sample set spans months."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-15T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAAGGG", "set_001", "sample_b", "2024-02-15T00:00:00Z")],
            "sample_b",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAATTT", "set_002", "sample_c", "2024-01-20T00:00:00Z")],
            "sample_c",
            "set_002",
            temp_state_dir,
        )

        compact_hits(temp_state_dir)
        assert count_hits(temp_state_dir) == 3

        result = delete_sample_set_hits("set_001", temp_state_dir)

        assert result == DeleteResult(
            uncompacted_files_deleted=0,
            compacted_months_rewritten=2,
        )
        assert count_hits(temp_state_dir) == 1

    def test_deletes_both_uncompacted_and_compacted(self, temp_state_dir):
        """Deletes from both uncompacted and compacted storage."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-15T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )

        compact_hits(temp_state_dir)

        write_hits_parquet(
            [make_hit_record("AAAGGG", "set_001", "sample_b", "2024-01-20T00:00:00Z")],
            "sample_b",
            "set_001",
            temp_state_dir,
        )

        assert count_hits(temp_state_dir) == 2

        result = delete_sample_set_hits("set_001", temp_state_dir)

        assert result == DeleteResult(
            uncompacted_files_deleted=1,
            compacted_months_rewritten=1,
        )
        assert count_hits(temp_state_dir) == 0

    def test_removes_empty_month_directory(self, temp_state_dir):
        """Removes month directory when all data is deleted."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-15T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )

        compact_hits(temp_state_dir)

        month_dir = temp_state_dir / "hits" / "month=2024-01"
        assert month_dir.exists()

        delete_sample_set_hits("set_001", temp_state_dir)

        assert not month_dir.exists()

    def test_raises_on_empty_sample_set_id(self, temp_state_dir):
        """Raises AssertionError for empty sample_set_id."""
        with pytest.raises(AssertionError, match="cannot be empty"):
            delete_sample_set_hits("", temp_state_dir)


class TestGetHitStats:
    """Tests for get_hit_stats()."""

    def test_empty_database(self, temp_state_dir):
        """Empty database returns zero counts and None for distributions."""
        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 0
        assert stats.total_observations == 0
        assert stats.unique_samples == 0
        assert stats.unique_runs == 0
        assert stats.length_min is None
        assert stats.length_max is None
        assert stats.length_median is None
        assert stats.gc_min is None
        assert stats.gc_max is None
        assert stats.gc_median is None
        assert stats.date_first is None
        assert stats.date_last is None

    def test_single_hit(self, temp_state_dir):
        """Single hit returns correct stats."""
        record = make_hit_record(
            "ACGTACGT",
            "run_001",
            "sample_a",
            "2024-01-15T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "run_001", temp_state_dir)

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 1
        assert stats.total_observations == 1
        assert stats.unique_samples == 1
        assert stats.unique_runs == 1
        assert stats.length_min == 8
        assert stats.length_max == 8
        assert stats.length_median == 8.0
        assert stats.gc_min == 0.5
        assert stats.gc_max == 0.5
        assert stats.gc_median == 0.5

    def test_multiple_hits_length_distribution(self, temp_state_dir):
        """Multiple hits with varying lengths."""
        records = [
            make_hit_record("ACGTACGT", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record(
                "ACGTACGTACGT",
                "run_001",
                "sample_a",
                "2024-01-02T00:00:00Z",
            ),
            make_hit_record(
                "ACGTACGTACGTACGT",
                "run_001",
                "sample_a",
                "2024-01-03T00:00:00Z",
            ),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 3
        assert stats.length_min == 8
        assert stats.length_max == 16
        assert stats.length_median == 12.0

    def test_multiple_hits_gc_distribution(self, temp_state_dir):
        """Multiple hits with varying GC content."""
        records = [
            make_hit_record("AAAAAAAA", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("ACGTACGT", "run_001", "sample_a", "2024-01-02T00:00:00Z"),
            make_hit_record("GCGCGCGC", "run_001", "sample_a", "2024-01-03T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 3
        assert stats.gc_min == 0.0
        assert stats.gc_max == 1.0
        assert stats.gc_median == 0.5

    def test_date_range(self, temp_state_dir):
        """Date range spans first and last seen dates."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-06-20T00:00:00Z"),
            make_hit_record("AAATTT", "run_001", "sample_a", "2024-03-10T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        stats = get_hit_stats(temp_state_dir)

        assert stats.date_first == "2024-01-15T00:00:00Z"
        assert stats.date_last == "2024-06-20T00:00:00Z"

    def test_multiple_runs_and_samples(self, temp_state_dir):
        """Correctly counts unique runs and samples."""
        # hit1 in run_001: sample_a, sample_b
        # hit1 in run_002: sample_a (same sample, different run)
        # hit2 in run_001: sample_c
        write_hits_parquet(
            [
                make_hit_record(
                    "AAACCC",
                    "run_001",
                    "sample_a",
                    "2024-01-01T00:00:00Z",
                ),
                make_hit_record(
                    "AAACCC",
                    "run_001",
                    "sample_b",
                    "2024-01-01T00:00:00Z",
                ),
            ],
            "sample_a",
            "run_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAACCC", "run_002", "sample_a", "2024-01-02T00:00:00Z")],
            "sample_a",
            "run_002",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAAGGG", "run_001", "sample_c", "2024-01-01T00:00:00Z")],
            "sample_c",
            "run_001",
            temp_state_dir,
        )

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 2
        assert stats.total_observations == 4
        assert stats.unique_samples == 3  # sample_a, sample_b, sample_c
        assert stats.unique_runs == 2  # run_001, run_002

    def test_median_even_count(self, temp_state_dir):
        """Median calculation for even number of values."""
        # Lengths: 8, 12, 16, 20 -> median = (12 + 16) / 2 = 14
        records = [
            make_hit_record("ACGTACGT", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record(
                "ACGTACGTACGT",
                "run_001",
                "sample_a",
                "2024-01-02T00:00:00Z",
            ),
            make_hit_record(
                "ACGTACGTACGTACGT",
                "run_001",
                "sample_a",
                "2024-01-03T00:00:00Z",
            ),
            make_hit_record(
                "ACGTACGTACGTACGTACGT",
                "run_001",
                "sample_a",
                "2024-01-04T00:00:00Z",
            ),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        stats = get_hit_stats(temp_state_dir)

        assert stats.length_median == 14.0


class TestGetRecurringHits:
    """Tests for get_recurring_hits()."""

    def test_empty_database(self, temp_state_dir):
        """Empty database returns empty list."""
        result = get_recurring_hits(state_dir=temp_state_dir)
        assert result == []

    def test_no_recurring_hits(self, temp_state_dir):
        """Hits with only one sample each don't appear."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-01T00:00:00Z")],
            "sample_a",
            "run_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAAGGG", "run_001", "sample_b", "2024-01-01T00:00:00Z")],
            "sample_b",
            "run_001",
            temp_state_dir,
        )

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)
        assert result == []

    def test_finds_recurring_hit(self, temp_state_dir):
        """Hit appearing in multiple samples is returned."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "run_001", "sample_b", "2024-01-02T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].sample_count == 2
        assert result[0].observation_count == 2

    def test_min_samples_filter(self, temp_state_dir):
        """min_samples parameter filters correctly."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "run_001", "sample_b", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "run_001", "sample_c", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        assert len(get_recurring_hits(min_samples=2, state_dir=temp_state_dir)) == 1
        assert len(get_recurring_hits(min_samples=3, state_dir=temp_state_dir)) == 1
        assert len(get_recurring_hits(min_samples=4, state_dir=temp_state_dir)) == 0

    def test_min_runs_filter(self, temp_state_dir):
        """min_runs parameter filters correctly."""
        # Observed in 2 samples across 2 runs
        write_hits_parquet(
            [make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-01T00:00:00Z")],
            "sample_a",
            "run_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [make_hit_record("AAACCC", "run_002", "sample_b", "2024-01-02T00:00:00Z")],
            "sample_b",
            "run_002",
            temp_state_dir,
        )

        # Should appear with min_runs=1 or 2, not 3
        assert (
            len(get_recurring_hits(min_samples=2, min_runs=1, state_dir=temp_state_dir))
            == 1
        )
        assert (
            len(get_recurring_hits(min_samples=2, min_runs=2, state_dir=temp_state_dir))
            == 1
        )
        assert (
            len(get_recurring_hits(min_samples=2, min_runs=3, state_dir=temp_state_dir))
            == 0
        )

    def test_limit_parameter(self, temp_state_dir):
        """limit parameter restricts results."""
        # Both hits in 2 samples
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "run_001", "sample_b", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_b", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        assert len(get_recurring_hits(min_samples=2, state_dir=temp_state_dir)) == 2
        assert (
            len(get_recurring_hits(min_samples=2, limit=1, state_dir=temp_state_dir))
            == 1
        )

    def test_ordered_by_sample_count(self, temp_state_dir):
        """Results ordered by sample_count descending."""
        # hit1 (AAACCC) in 2 samples, hit2 (AAAGGG) in 3 samples
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAACCC", "run_001", "sample_b", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_b", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_c", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 2
        assert result[0].sample_count == 3  # AAAGGG first (more samples)
        assert result[1].sample_count == 2  # AAACCC second

    def test_last_seen_date(self, temp_state_dir):
        """last_seen_date is the most recent observation."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("AAACCC", "run_001", "sample_b", "2024-06-20T00:00:00Z"),
            make_hit_record("AAACCC", "run_002", "sample_c", "2024-03-10T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].first_seen_date == "2024-01-15T00:00:00Z"
        assert result[0].last_seen_date == "2024-06-20T00:00:00Z"

    def test_returns_correct_metadata(self, temp_state_dir):
        """Returned RecurringHit has correct hit metadata."""
        records = [
            make_hit_record(
                "GCGCGCGC",
                "run_001",
                "sample_a",
                "2024-01-01T00:00:00Z",
            ),  # 100% GC, 8 bp
            make_hit_record("GCGCGCGC", "run_001", "sample_b", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].sequence_length == 8
        assert result[0].gc_content == 1.0


class TestGetDiscoveryTimeline:
    """Tests for get_discovery_timeline()."""

    def test_empty_database(self, temp_state_dir):
        """Empty database returns empty list."""
        result = get_discovery_timeline(state_dir=temp_state_dir)
        assert result == []

    def test_single_hit(self, temp_state_dir):
        """Single hit returns one bucket."""
        record = make_hit_record(
            "AAACCC",
            "run_001",
            "sample_a",
            "2024-01-15T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "run_001", temp_state_dir)

        result = get_discovery_timeline(granularity="month", state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].period == "2024-01"
        assert result[0].new_hits == 1

    def test_multiple_hits_same_month(self, temp_state_dir):
        """Multiple hits in same month are grouped."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-10T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-01-20T00:00:00Z"),
            make_hit_record("AAATTT", "run_001", "sample_a", "2024-01-25T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_discovery_timeline(granularity="month", state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].period == "2024-01"
        assert result[0].new_hits == 3

    def test_multiple_months(self, temp_state_dir):
        """Hits across months create separate buckets."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-01-20T00:00:00Z"),
            make_hit_record("AAATTT", "run_001", "sample_a", "2024-03-10T00:00:00Z"),
            make_hit_record("CCCGGG", "run_001", "sample_a", "2024-06-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_discovery_timeline(granularity="month", state_dir=temp_state_dir)

        assert len(result) == 3
        assert result[0].period == "2024-01"
        assert result[0].new_hits == 2
        assert result[1].period == "2024-03"
        assert result[1].new_hits == 1
        assert result[2].period == "2024-06"
        assert result[2].new_hits == 1

    def test_ordered_ascending(self, temp_state_dir):
        """Results are ordered by period ascending."""
        # Register in non-chronological order
        records = [
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-06-01T00:00:00Z"),
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("AAATTT", "run_001", "sample_a", "2024-03-10T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_discovery_timeline(granularity="month", state_dir=temp_state_dir)

        periods = [b.period for b in result]
        assert periods == ["2024-01", "2024-03", "2024-06"]

    def test_granularity_day(self, temp_state_dir):
        """Day granularity groups by exact date."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-01-15T12:00:00Z"),
            make_hit_record("AAATTT", "run_001", "sample_a", "2024-01-16T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_discovery_timeline(granularity="day", state_dir=temp_state_dir)

        assert len(result) == 2
        assert result[0].period == "2024-01-15"
        assert result[0].new_hits == 2
        assert result[1].period == "2024-01-16"
        assert result[1].new_hits == 1

    def test_granularity_year(self, temp_state_dir):
        """Year granularity groups by year."""
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2023-06-15T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("AAATTT", "run_001", "sample_a", "2024-12-31T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_discovery_timeline(granularity="year", state_dir=temp_state_dir)

        assert len(result) == 2
        assert result[0].period == "2023"
        assert result[0].new_hits == 1
        assert result[1].period == "2024"
        assert result[1].new_hits == 2

    def test_granularity_week(self, temp_state_dir):
        """Week granularity groups by ISO week."""
        # 2024-01-01 is Monday of week 01
        # 2024-01-08 is Monday of week 02
        records = [
            make_hit_record("AAACCC", "run_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("AAAGGG", "run_001", "sample_a", "2024-01-07T00:00:00Z"),
            make_hit_record("AAATTT", "run_001", "sample_a", "2024-01-08T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "run_001", temp_state_dir)

        result = get_discovery_timeline(granularity="week", state_dir=temp_state_dir)

        assert len(result) == 2
        # First bucket should have 2 hits, second should have 1
        assert result[0].new_hits == 2
        assert result[1].new_hits == 1


class TestLookupHit:
    """Tests for lookup_hit()."""

    def test_lookup_by_key(self, temp_state_dir):
        """Lookup by hit key returns the hit."""
        seq = "ACGTACGTACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = lookup_hit(record.hit_key, state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == record.hit_key

    def test_lookup_by_sequence(self, temp_state_dir):
        """Lookup by sequence returns the hit."""
        seq = "ACGTACGTACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = lookup_hit(seq, state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == record.hit_key

    def test_lookup_by_reverse_complement(self, temp_state_dir):
        """Lookup by reverse complement finds the same hit."""
        seq = "AAACCCGGG"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        revcomp = reverse_complement(seq)
        result = lookup_hit(revcomp, state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == record.hit_key

    def test_lookup_not_found(self, temp_state_dir):
        """Lookup by nonexistent key returns None."""
        result = lookup_hit("a" * 32, state_dir=temp_state_dir)
        assert result is None

    def test_lookup_by_key_uppercase(self, temp_state_dir):
        """Lookup by uppercase hit key works."""
        seq = "ACGTACGTACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = lookup_hit(record.hit_key.upper(), state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == record.hit_key

    def test_lookup_by_sequence_lowercase(self, temp_state_dir):
        """Lookup by lowercase sequence works."""
        seq = "ACGTACGTACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = lookup_hit(seq.lower(), state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == record.hit_key

    def test_lookup_not_found_key(self, temp_state_dir):
        """Lookup by nonexistent key returns None."""
        result = lookup_hit("a" * 32, state_dir=temp_state_dir)
        assert result is None

    def test_lookup_not_found_sequence(self, temp_state_dir):
        """Lookup by nonexistent sequence returns None."""
        result = lookup_hit("ACGTACGTACGT", state_dir=temp_state_dir)
        assert result is None

    def test_lookup_strips_whitespace(self, temp_state_dir):
        """Lookup strips leading/trailing whitespace."""
        seq = "ACGTACGTACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = lookup_hit(f"  {seq}  \n", state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == record.hit_key

    def test_lookup_key_with_whitespace(self, temp_state_dir):
        """Lookup key with whitespace is stripped."""
        seq = "ACGTACGTACGT"
        record = make_hit_record(seq, "set_001", "sample_a", "2024-01-01T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = lookup_hit(f"  {record.hit_key}  ", state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == record.hit_key


class TestIsValidHitKey:
    """Tests for is_valid_hit_key()."""

    def test_valid_key(self):
        """32 lowercase hex characters is valid."""
        assert is_valid_hit_key("a" * 32) is True
        assert is_valid_hit_key("0123456789abcdef" * 2) is True
        assert is_valid_hit_key("deadbeef" * 4) is True

    def test_wrong_length_short(self):
        """Keys shorter than 32 characters are invalid."""
        assert is_valid_hit_key("a" * 31) is False
        assert is_valid_hit_key("a" * 16) is False
        assert is_valid_hit_key("a") is False
        assert is_valid_hit_key("") is False

    def test_wrong_length_long(self):
        """Keys longer than 32 characters are invalid."""
        assert is_valid_hit_key("a" * 33) is False
        assert is_valid_hit_key("a" * 64) is False

    def test_uppercase_invalid(self):
        """Uppercase hex characters are invalid."""
        assert is_valid_hit_key("A" * 32) is False
        assert is_valid_hit_key("ABCDEF" + "a" * 26) is False
        assert is_valid_hit_key("a" * 31 + "A") is False

    def test_non_hex_invalid(self):
        """Non-hex characters are invalid."""
        assert is_valid_hit_key("g" * 32) is False
        assert is_valid_hit_key("z" * 32) is False
        assert is_valid_hit_key("-" * 32) is False
        assert is_valid_hit_key(" " * 32) is False
        assert is_valid_hit_key("a" * 31 + "g") is False


class TestGetSampleParquetPath:
    """Tests for get_sample_parquet_path()."""

    def test_path_structure(self, temp_state_dir):
        """Returns correct path structure."""
        path = get_sample_parquet_path("sample_a", "set_001", temp_state_dir)

        # hits/month=NULL/{sample_set_id}/{sample_id}/data.parquet
        assert path.name == "data.parquet"
        assert path.parent.name == "sample_a"
        assert path.parent.parent.name == "set_001"
        assert path.parent.parent.parent.name == "month=NULL"
        assert "hits" in str(path)

    def test_does_not_create_directories(self, temp_state_dir):
        """Does not create the file or directories."""
        path = get_sample_parquet_path("sample_a", "set_001", temp_state_dir)

        assert not path.exists()
        assert not path.parent.exists()

    def test_different_samples_different_paths(self, temp_state_dir):
        """Different sample_ids produce different paths."""
        path_a = get_sample_parquet_path("sample_a", "set_001", temp_state_dir)
        path_b = get_sample_parquet_path("sample_b", "set_001", temp_state_dir)

        assert path_a != path_b
        # Same filename (data.parquet), different parent directories (sample_a vs sample_b)
        assert path_a.name == path_b.name
        assert path_a.parent != path_b.parent
        assert path_a.parent.parent == path_b.parent.parent  # Same sample set directory

    def test_different_sets_different_paths(self, temp_state_dir):
        """Different sample_set_ids produce different paths."""
        path_1 = get_sample_parquet_path("sample_a", "set_001", temp_state_dir)
        path_2 = get_sample_parquet_path("sample_a", "set_002", temp_state_dir)

        assert path_1 != path_2
        assert path_1.name == path_2.name  # Same filename (data.parquet)
        assert path_1.parent.name == path_2.parent.name  # Same sample_id directory name
        assert (
            path_1.parent.parent != path_2.parent.parent
        )  # Different sample_set directories


class TestSampleHitsExist:
    """Tests for sample_hits_exist()."""

    def test_exists_after_write(self, temp_state_dir):
        """Returns True after writing parquet file."""
        record = make_hit_record(
            "ACGTACGT",
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        assert sample_hits_exist("sample_a", "set_001", temp_state_dir) is True

    def test_not_exists_before_write(self, temp_state_dir):
        """Returns False when no parquet file exists."""
        assert sample_hits_exist("sample_a", "set_001", temp_state_dir) is False

    def test_not_exists_wrong_sample(self, temp_state_dir):
        """Returns False for different sample_id."""
        record = make_hit_record(
            "ACGTACGT",
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        assert sample_hits_exist("sample_b", "set_001", temp_state_dir) is False

    def test_not_exists_wrong_set(self, temp_state_dir):
        """Returns False for different sample_set_id."""
        record = make_hit_record(
            "ACGTACGT",
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        assert sample_hits_exist("sample_a", "set_002", temp_state_dir) is False

    def test_false_after_delete(self, temp_state_dir):
        """Returns False after deleting the parquet file."""
        record = make_hit_record(
            "ACGTACGT",
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
        )
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        assert sample_hits_exist("sample_a", "set_001", temp_state_dir) is True

        delete_sample_hits("sample_a", "set_001", temp_state_dir)

        assert sample_hits_exist("sample_a", "set_001", temp_state_dir) is False


class TestCountSampleSetObservations:
    """Tests for count_sample_set_observations()."""

    def test_empty_sample_set(self, temp_state_dir):
        """Returns 0 for nonexistent sample set."""
        assert count_sample_set_observations("nonexistent", temp_state_dir) == 0

    def test_counts_observations_single_file(self, temp_state_dir):
        """Counts all observations in a single parquet file."""
        records = [
            make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("CCCC", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
            make_hit_record("GGGG", "set_001", "sample_a", "2024-01-01T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        assert count_sample_set_observations("set_001", temp_state_dir) == 3

    def test_counts_across_multiple_samples(self, temp_state_dir):
        """Counts observations across multiple sample files in same set."""
        write_hits_parquet(
            [make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [
                make_hit_record("CCCC", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
                make_hit_record("GGGG", "set_001", "sample_b", "2024-01-01T00:00:00Z"),
            ],
            "sample_b",
            "set_001",
            temp_state_dir,
        )

        assert count_sample_set_observations("set_001", temp_state_dir) == 3

    def test_only_counts_specified_set(self, temp_state_dir):
        """Only counts observations in the specified sample set."""
        write_hits_parquet(
            [make_hit_record("AAAA", "set_001", "sample_a", "2024-01-01T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [
                make_hit_record("CCCC", "set_002", "sample_b", "2024-01-01T00:00:00Z"),
                make_hit_record("GGGG", "set_002", "sample_b", "2024-01-01T00:00:00Z"),
            ],
            "sample_b",
            "set_002",
            temp_state_dir,
        )

        assert count_sample_set_observations("set_001", temp_state_dir) == 1
        assert count_sample_set_observations("set_002", temp_state_dir) == 2


class TestCompactHits:
    """Tests for compact_hits()."""

    def test_compact_empty_returns_empty_result(self, temp_state_dir):
        """Compacting with no uncompacted data returns empty result."""
        result = compact_hits(temp_state_dir)

        assert result.dry_run is False
        assert result.months == []
        assert result.total_observations == 0
        assert result.total_source_files == 0

    def test_dry_run_does_not_modify(self, temp_state_dir):
        """Dry run reports what would be compacted without modifying."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = compact_hits(temp_state_dir, dry_run=True)

        assert result.dry_run is True
        assert len(result.months) == 1
        assert result.months[0].month == "2024-01"
        assert result.total_observations == 1

        # Verify uncompacted file still exists
        uncompacted_path = (
            temp_state_dir
            / "hits"
            / "month=NULL"
            / "set_001"
            / "sample_a"
            / "data.parquet"
        )
        assert uncompacted_path.exists()

        # Verify compacted file was NOT created
        compacted_path = temp_state_dir / "hits" / "month=2024-01" / "data.parquet"
        assert not compacted_path.exists()

    def test_basic_compaction(self, temp_state_dir):
        """Basic compaction moves data from month=NULL to month=YYYY-MM."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = compact_hits(temp_state_dir)

        assert result.dry_run is False
        assert len(result.months) == 1
        assert result.months[0].month == "2024-01"
        assert result.total_observations == 1
        assert result.errors == []

        # Verify compacted file exists
        compacted_path = temp_state_dir / "hits" / "month=2024-01" / "data.parquet"
        assert compacted_path.exists()

        # Verify uncompacted file was deleted (default behavior)
        uncompacted_path = (
            temp_state_dir
            / "hits"
            / "month=NULL"
            / "set_001"
            / "sample_a"
            / "data.parquet"
        )
        assert not uncompacted_path.exists()

    def test_keep_source_preserves_uncompacted(self, temp_state_dir):
        """keep_source=True preserves uncompacted files after compaction."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = compact_hits(temp_state_dir, keep_source=True)

        assert result.errors == []

        # Both files should exist
        compacted_path = temp_state_dir / "hits" / "month=2024-01" / "data.parquet"
        uncompacted_path = (
            temp_state_dir
            / "hits"
            / "month=NULL"
            / "set_001"
            / "sample_a"
            / "data.parquet"
        )
        assert compacted_path.exists()
        assert uncompacted_path.exists()

    def test_compaction_groups_by_month(self, temp_state_dir):
        """Data from different months goes to different partitions."""
        jan_record = make_hit_record(
            "AAAA",
            "set_001",
            "sample_a",
            "2024-01-15T00:00:00Z",
        )
        feb_record = make_hit_record(
            "CCCC",
            "set_001",
            "sample_b",
            "2024-02-15T00:00:00Z",
        )
        write_hits_parquet([jan_record], "sample_a", "set_001", temp_state_dir)
        write_hits_parquet([feb_record], "sample_b", "set_001", temp_state_dir)

        result = compact_hits(temp_state_dir)

        assert len(result.months) == 2
        months = {m.month for m in result.months}
        assert months == {"2024-01", "2024-02"}

        # Both month partitions should exist
        assert (temp_state_dir / "hits" / "month=2024-01" / "data.parquet").exists()
        assert (temp_state_dir / "hits" / "month=2024-02" / "data.parquet").exists()

    def test_month_filter_only_compacts_specified_month(self, temp_state_dir):
        """month parameter filters compaction to a single month."""
        jan_record = make_hit_record(
            "AAAA",
            "set_001",
            "sample_a",
            "2024-01-15T00:00:00Z",
        )
        feb_record = make_hit_record(
            "CCCC",
            "set_001",
            "sample_b",
            "2024-02-15T00:00:00Z",
        )
        write_hits_parquet([jan_record], "sample_a", "set_001", temp_state_dir)
        write_hits_parquet([feb_record], "sample_b", "set_001", temp_state_dir)

        result = compact_hits(temp_state_dir, month="2024-01")

        assert len(result.months) == 1
        assert result.months[0].month == "2024-01"

        # Only January should be compacted
        assert (temp_state_dir / "hits" / "month=2024-01" / "data.parquet").exists()
        assert not (temp_state_dir / "hits" / "month=2024-02" / "data.parquet").exists()

        # February uncompacted file should still exist
        feb_uncompacted = (
            temp_state_dir
            / "hits"
            / "month=NULL"
            / "set_001"
            / "sample_b"
            / "data.parquet"
        )
        assert feb_uncompacted.exists()

    def test_incremental_compaction_merges_with_existing(self, temp_state_dir):
        """Compacting into existing month file merges the data."""
        # First compaction
        record1 = make_hit_record("AAAA", "set_001", "sample_a", "2024-01-10T00:00:00Z")
        write_hits_parquet([record1], "sample_a", "set_001", temp_state_dir)
        compact_hits(temp_state_dir)

        # Second batch of data for same month
        record2 = make_hit_record("CCCC", "set_002", "sample_b", "2024-01-20T00:00:00Z")
        write_hits_parquet([record2], "sample_b", "set_002", temp_state_dir)
        result = compact_hits(temp_state_dir)

        assert result.errors == []

        # Query the compacted file - should have both records
        import duckdb

        compacted_path = temp_state_dir / "hits" / "month=2024-01" / "data.parquet"
        con = duckdb.connect()
        result = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{compacted_path}')",
        ).fetchone()
        assert result is not None
        assert result[0] == 2

    def test_compacted_data_is_sorted_by_hit_key(self, temp_state_dir):
        """Compacted data is sorted by hit_key, then run_date."""
        # Create records with different hit_keys (AAAA < CCCC < GGGG lexicographically)
        records = [
            make_hit_record("GGGG", "set_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("AAAA", "set_001", "sample_a", "2024-01-15T00:00:00Z"),
            make_hit_record("CCCC", "set_001", "sample_a", "2024-01-15T00:00:00Z"),
        ]
        write_hits_parquet(records, "sample_a", "set_001", temp_state_dir)

        compact_hits(temp_state_dir)

        # Read compacted file and verify sort order
        import duckdb

        compacted_path = temp_state_dir / "hits" / "month=2024-01" / "data.parquet"
        con = duckdb.connect()
        hit_keys = con.execute(
            f"SELECT hit_key FROM read_parquet('{compacted_path}')",
        ).fetchall()
        hit_keys = [row[0] for row in hit_keys]

        # Should be sorted (hit_key for AAAA < CCCC < GGGG)
        assert hit_keys == sorted(hit_keys)

    def test_compacted_data_queryable_with_hive_partitioning(self, temp_state_dir):
        """Compacted data is queryable via the standard query_hits() interface."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)
        compact_hits(temp_state_dir)

        # Query using the standard interface
        con = query_hits(temp_state_dir)
        result = con.execute(
            "SELECT hit_key, month, effective_month FROM hits",
        ).fetchone()

        assert result is not None
        assert result[1] == "2024-01"  # month extracted from Hive partition
        assert result[2] == "2024-01"  # effective_month matches

    def test_incremental_compaction_handles_schema_evolution(self, temp_state_dir):
        """Compaction merges old-schema compacted files with new-schema uncompacted files."""
        import duckdb

        # Create an "old" compacted file without classification columns
        # (simulating data compacted before the schema change)
        month_dir = temp_state_dir / "hits" / "month=2024-01"
        month_dir.mkdir(parents=True)
        old_compacted = month_dir / "data.parquet"

        con = duckdb.connect()
        con.execute(f"""
            COPY (
                SELECT
                    'abc123def456abc123def456abc123de' as hit_key,
                    100 as sequence_length,
                    'AAAA'::BLOB as sequence_compressed,
                    0.5 as gc_content,
                    'set-old' as sample_set_id,
                    'sample-old' as sample_id,
                    '2024-01-15'::DATE as run_date,
                    'contig-old' as contig_id
            )
            TO '{old_compacted}' (FORMAT parquet)
        """)

        # Write new uncompacted data WITH classification columns
        record = HitRecord(
            hit_key="def456abc123def456abc123def456ab",
            sequence_length=200,
            sequence_compressed=b"TTTT",
            gc_content=0.6,
            sample_set_id="set-new",
            sample_id="sample-new",
            run_date="2024-01-20",
            contig_id="contig-new",
            adjusted_taxid=12345,
            adjusted_taxid_name="Test virus",
            adjusted_taxid_rank="species",
        )
        write_hits_parquet([record], "sample-new", "set-new", temp_state_dir)

        # Compact - this should merge old (no classification cols) with new (has them)
        result = compact_hits(state_dir=temp_state_dir, month="2024-01")

        assert result.total_observations == 1  # Only the new uncompacted file
        assert len(result.months) == 1
        assert result.errors == []

        # Verify the compacted file has both records with correct schema
        compacted = con.execute(f"""
            SELECT hit_key, adjusted_taxid, adjusted_taxid_name, adjusted_taxid_rank
            FROM read_parquet('{old_compacted}')
            ORDER BY hit_key
        """).fetchall()

        assert len(compacted) == 2
        # Old record should have NULLs for classification columns
        assert compacted[0] == ("abc123def456abc123def456abc123de", None, None, None)
        # New record should have the classification data
        assert compacted[1] == (
            "def456abc123def456abc123def456ab",
            12345,
            "Test virus",
            "species",
        )


class TestGetStatsForSampleSet:
    """Tests for get_stats_for_sample_set()."""

    def test_empty_database_returns_none(self, temp_state_dir):
        """Empty database returns None."""
        result = get_stats_for_sample_set("nonexistent", temp_state_dir)
        assert result is None

    def test_nonexistent_sample_set_returns_none(self, temp_state_dir):
        """Sample set with no data returns None."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = get_stats_for_sample_set("set_002", temp_state_dir)
        assert result is None

    def test_single_sample_single_hit(self, temp_state_dir):
        """Single sample with single hit returns correct stats."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = get_stats_for_sample_set("set_001", temp_state_dir)

        assert result is not None
        assert result.sample_set_id == "set_001"
        assert result.total_observations == 1
        assert result.unique_hits == 1
        assert result.unique_samples == 1
        assert result.date_first == "2024-01-15T00:00:00Z"
        assert result.date_last == "2024-01-15T00:00:00Z"

    def test_multiple_samples_multiple_hits(self, temp_state_dir):
        """Multiple samples with multiple hits returns correct stats."""
        # Use sequences that produce distinct hit keys (not reverse complements)
        write_hits_parquet(
            [
                make_hit_record(
                    "AAACCC",
                    "set_001",
                    "sample_a",
                    "2024-01-10T00:00:00Z",
                ),
                make_hit_record(
                    "AAAGGG",
                    "set_001",
                    "sample_a",
                    "2024-01-10T00:00:00Z",
                ),
            ],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [
                make_hit_record(
                    "AAACCC",
                    "set_001",
                    "sample_b",
                    "2024-01-15T00:00:00Z",
                ),
                make_hit_record(
                    "AAATTT",
                    "set_001",
                    "sample_b",
                    "2024-01-15T00:00:00Z",
                ),
            ],
            "sample_b",
            "set_001",
            temp_state_dir,
        )

        result = get_stats_for_sample_set("set_001", temp_state_dir)

        assert result is not None
        assert result.sample_set_id == "set_001"
        assert result.total_observations == 4
        assert result.unique_hits == 3  # AAACCC, AAAGGG, AAATTT (AAACCC appears twice)
        assert result.unique_samples == 2  # sample_a, sample_b
        assert result.date_first == "2024-01-10T00:00:00Z"
        assert result.date_last == "2024-01-15T00:00:00Z"

    def test_filters_to_specified_sample_set(self, temp_state_dir):
        """Only counts data from the specified sample set."""
        write_hits_parquet(
            [make_hit_record("AAACCC", "set_001", "sample_a", "2024-01-10T00:00:00Z")],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        write_hits_parquet(
            [
                make_hit_record(
                    "AAAGGG",
                    "set_002",
                    "sample_b",
                    "2024-01-15T00:00:00Z",
                ),
                make_hit_record(
                    "AAATTT",
                    "set_002",
                    "sample_b",
                    "2024-01-15T00:00:00Z",
                ),
            ],
            "sample_b",
            "set_002",
            temp_state_dir,
        )

        result_001 = get_stats_for_sample_set("set_001", temp_state_dir)
        result_002 = get_stats_for_sample_set("set_002", temp_state_dir)

        assert result_001 is not None
        assert result_001.total_observations == 1
        assert result_001.unique_hits == 1

        assert result_002 is not None
        assert result_002.total_observations == 2
        assert result_002.unique_hits == 2

    def test_same_hit_multiple_observations_same_sample(self, temp_state_dir):
        """Same hit observed multiple times in one sample counts correctly."""
        write_hits_parquet(
            [
                make_hit_record(
                    "AAACCC",
                    "set_001",
                    "sample_a",
                    "2024-01-10T00:00:00Z",
                    "contig_1",
                ),
                make_hit_record(
                    "AAACCC",
                    "set_001",
                    "sample_a",
                    "2024-01-10T00:00:00Z",
                    "contig_2",
                ),
            ],
            "sample_a",
            "set_001",
            temp_state_dir,
        )

        result = get_stats_for_sample_set("set_001", temp_state_dir)

        assert result is not None
        assert result.total_observations == 2  # Two observations
        assert result.unique_hits == 1  # But only one unique hit
        assert result.unique_samples == 1  # And one sample

    def test_works_with_compacted_data(self, temp_state_dir):
        """Works correctly after compaction."""
        write_hits_parquet(
            [
                make_hit_record(
                    "AAACCC",
                    "set_001",
                    "sample_a",
                    "2024-01-10T00:00:00Z",
                ),
                make_hit_record(
                    "AAAGGG",
                    "set_001",
                    "sample_a",
                    "2024-01-15T00:00:00Z",
                ),
            ],
            "sample_a",
            "set_001",
            temp_state_dir,
        )
        compact_hits(temp_state_dir)

        result = get_stats_for_sample_set("set_001", temp_state_dir)

        assert result is not None
        assert result.sample_set_id == "set_001"
        assert result.total_observations == 2
        assert result.unique_hits == 2
        assert result.unique_samples == 1
        assert result.date_first == "2024-01-10T00:00:00Z"
        assert result.date_last == "2024-01-15T00:00:00Z"

    def test_returns_sample_set_stats_type(self, temp_state_dir):
        """Returns a SampleSetStats instance."""
        record = make_hit_record("ACGT", "set_001", "sample_a", "2024-01-15T00:00:00Z")
        write_hits_parquet([record], "sample_a", "set_001", temp_state_dir)

        result = get_stats_for_sample_set("set_001", temp_state_dir)

        assert isinstance(result, SampleSetStats)

    def test_raises_on_empty_sample_set_id(self, temp_state_dir):
        """Raises AssertionError for empty sample_set_id."""
        with pytest.raises(AssertionError, match="cannot be empty"):
            get_stats_for_sample_set("", temp_state_dir)
