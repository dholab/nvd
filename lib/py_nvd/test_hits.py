"""Tests for py_nvd.hits module.

These tests verify:
- Sequence operations (reverse complement, canonical form)
- Hit key computation (deterministic, strand-agnostic)
- Sequence compression/decompression (round-trip, N handling)
- GC content calculation
- Hit registration and observation recording (state management)
- Hit queries
"""

import tempfile
from pathlib import Path

import pytest

from py_nvd.hits import (
    calculate_gc_content,
    canonical_sequence,
    compress_sequence,
    compute_hit_key,
    count_hit_observations,
    count_hits,
    decompress_sequence,
    delete_observations_for_sample_set,
    delete_orphaned_hits,
    get_discovery_timeline,
    get_hit,
    get_hit_sequence,
    get_hit_stats,
    get_recurring_hits,
    list_hit_observations,
    list_hits,
    list_hits_with_observations,
    lookup_hit,
    prune_hits_for_sample_set,
    record_hit_observation,
    register_hit,
    reverse_complement,
)


# =============================================================================
# Sequence Operations Tests
# =============================================================================


class TestReverseComplement:
    def test_palindrome(self):
        # ACGT is its own reverse complement
        assert reverse_complement("ACGT") == "ACGT"

    def test_asymmetric(self):
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("CCCC") == "GGGG"
        assert reverse_complement("AAACCC") == "GGGTTT"

    def test_palindrome_longer(self):
        # GCGC is a palindrome
        assert reverse_complement("GCGC") == "GCGC"

    def test_with_n(self):
        # N complements to N
        assert reverse_complement("ACNGT") == "ACNGT"
        assert reverse_complement("NNNN") == "NNNN"

    def test_lowercase(self):
        # Should handle lowercase
        assert reverse_complement("acgt") == "acgt"
        assert reverse_complement("aaaa") == "tttt"

    def test_mixed_case(self):
        assert reverse_complement("AcGt") == "aCgT"


class TestCanonicalSequence:
    def test_already_canonical(self):
        # "AAAA" < "TTTT" lexicographically, so AAAA is canonical
        assert canonical_sequence("AAAA") == "AAAA"

    def test_needs_revcomp(self):
        # "TTTT" > "AAAA", so canonical is AAAA
        assert canonical_sequence("TTTT") == "AAAA"

    def test_palindrome(self):
        # ACGT is its own reverse complement, so it's canonical
        assert canonical_sequence("ACGT") == "ACGT"

    def test_lowercase_normalized(self):
        # Should uppercase before comparison
        assert canonical_sequence("aaaa") == "AAAA"
        assert canonical_sequence("tttt") == "AAAA"

    def test_complex_sequence(self):
        # AAACCCGGG vs reverse complement CCCGGGTTT
        # "AAACCCGGG" < "CCCGGGTTT", so original is canonical
        seq = "AAACCCGGG"
        assert canonical_sequence(seq) == "AAACCCGGG"

    def test_complex_sequence_reverse(self):
        # CCCGGGTTT vs reverse complement AAACCCGGG
        # "AAACCCGGG" < "CCCGGGTTT", so revcomp is canonical
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
        assert len(key) == 32  # 128 bits = 32 hex chars

    def test_hex_format(self):
        key = compute_hit_key("ACGT")
        # Should be valid hex
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
        # R, Y, S, W are ambiguous IUPAC codes - should become N
        seq = "ACRYSWGT"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == "ACNNNNGT"

    def test_all_iupac_ambiguous(self):
        # All ambiguous codes: R, Y, S, W, K, M, B, D, H, V
        seq = "RYSWKMBDHV"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == "NNNNNNNNNN"

    def test_compression_ratio(self):
        seq = "ACGT" * 1000  # 4000 bp
        compressed = compress_sequence(seq)
        # Should be roughly 25% or better (2-bit encoding + zlib)
        assert len(compressed) < len(seq) * 0.30

    def test_lowercase_normalized(self):
        seq = "acgtacgt"
        compressed = compress_sequence(seq)
        decompressed = decompress_sequence(compressed, len(seq))
        assert decompressed == "ACGTACGT"

    def test_odd_length(self):
        # Test sequences that don't divide evenly by 4
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
        # N's don't count as GC
        assert calculate_gc_content("GCNN") == 0.5

    def test_lowercase(self):
        assert calculate_gc_content("gcgc") == 1.0

    def test_typical_viral(self):
        # ~40% GC is typical for many viruses
        seq = "AATTAATTGC"  # 2 GC out of 10
        assert calculate_gc_content(seq) == 0.2


# =============================================================================
# State Management Tests
# =============================================================================


class TestRegisterHit:
    """Tests for register_hit()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_register_new_hit(self, temp_state_dir):
        """Registering a new hit succeeds and returns is_new=True."""
        seq = "ACGTACGTACGT"
        hit, is_new = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        assert is_new is True
        assert hit.hit_key == compute_hit_key(seq)
        assert hit.sequence_length == len(seq)
        assert hit.gc_content == 0.5
        assert hit.first_seen_date == "2024-01-01T00:00:00Z"

    def test_register_duplicate_hit(self, temp_state_dir):
        """Registering the same sequence twice returns is_new=False."""
        seq = "ACGTACGTACGT"

        hit1, is_new1 = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)
        hit2, is_new2 = register_hit(seq, "2024-01-02T00:00:00Z", temp_state_dir)

        assert is_new1 is True
        assert is_new2 is False
        assert hit1.hit_key == hit2.hit_key
        # First seen date should be preserved
        assert hit2.first_seen_date == "2024-01-01T00:00:00Z"

    def test_register_strand_agnostic(self, temp_state_dir):
        """Forward and reverse complement sequences produce same hit."""
        seq = "AAACCCGGG"
        revcomp = reverse_complement(seq)

        hit1, is_new1 = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)
        hit2, is_new2 = register_hit(revcomp, "2024-01-02T00:00:00Z", temp_state_dir)

        assert is_new1 is True
        assert is_new2 is False
        assert hit1.hit_key == hit2.hit_key

    def test_sequence_round_trip(self, temp_state_dir):
        """Registered sequence can be recovered via decompression."""
        seq = "ACGTNNNACGT"
        hit, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        recovered = get_hit_sequence(hit)
        assert recovered == seq


class TestRecordHitObservation:
    """Tests for record_hit_observation()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_record_new_observation(self, temp_state_dir):
        """Recording a new observation succeeds."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        obs = record_hit_observation(
            hit_key=hit.hit_key,
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
            contig_id="NODE_1",
            state_dir=temp_state_dir,
        )

        assert obs is not None
        assert obs.hit_key == hit.hit_key
        assert obs.sample_set_id == "set_001"
        assert obs.sample_id == "sample_a"
        assert obs.contig_id == "NODE_1"

    def test_record_duplicate_observation(self, temp_state_dir):
        """Recording the same observation twice returns None."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        obs1 = record_hit_observation(
            hit_key=hit.hit_key,
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        obs2 = record_hit_observation(
            hit_key=hit.hit_key,
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        assert obs1 is not None
        assert obs2 is None  # Duplicate

    def test_same_hit_different_samples(self, temp_state_dir):
        """Same hit can be observed in different samples."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        obs1 = record_hit_observation(
            hit_key=hit.hit_key,
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        obs2 = record_hit_observation(
            hit_key=hit.hit_key,
            sample_set_id="set_001",
            sample_id="sample_b",
            run_date="2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        assert obs1 is not None
        assert obs2 is not None
        assert obs1.sample_id != obs2.sample_id


class TestGetHit:
    """Tests for get_hit()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_get_existing_hit(self, temp_state_dir):
        """Getting an existing hit returns it."""
        seq = "ACGTACGTACGT"
        registered, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        retrieved = get_hit(registered.hit_key, temp_state_dir)

        assert retrieved is not None
        assert retrieved.hit_key == registered.hit_key
        assert retrieved.sequence_length == registered.sequence_length

    def test_get_nonexistent_hit(self, temp_state_dir):
        """Getting a nonexistent hit returns None."""
        result = get_hit("nonexistent_key", temp_state_dir)
        assert result is None


class TestListHits:
    """Tests for list_hits()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_list_empty(self, temp_state_dir):
        """Listing hits from empty database returns empty list."""
        hits = list_hits(state_dir=temp_state_dir)
        assert hits == []

    def test_list_multiple_hits(self, temp_state_dir):
        """Listing hits returns all registered hits."""
        # Use sequences that are NOT reverse complements of each other
        register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        register_hit("AAAGGG", "2024-01-02T00:00:00Z", temp_state_dir)
        register_hit("AAATTT", "2024-01-03T00:00:00Z", temp_state_dir)

        hits = list_hits(state_dir=temp_state_dir)

        assert len(hits) == 3

    def test_list_with_limit(self, temp_state_dir):
        """Listing hits respects limit parameter."""
        # Use sequences that are NOT reverse complements of each other
        register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        register_hit("AAAGGG", "2024-01-02T00:00:00Z", temp_state_dir)
        register_hit("AAATTT", "2024-01-03T00:00:00Z", temp_state_dir)

        hits = list_hits(limit=2, state_dir=temp_state_dir)

        assert len(hits) == 2

    def test_list_ordered_by_date_desc(self, temp_state_dir):
        """Hits are ordered by first_seen_date descending."""
        # Use sequences that are NOT reverse complements of each other
        register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        register_hit("AAAGGG", "2024-01-03T00:00:00Z", temp_state_dir)
        register_hit("AAATTT", "2024-01-02T00:00:00Z", temp_state_dir)

        hits = list_hits(state_dir=temp_state_dir)

        # Most recent first
        assert hits[0].first_seen_date == "2024-01-03T00:00:00Z"
        assert hits[1].first_seen_date == "2024-01-02T00:00:00Z"
        assert hits[2].first_seen_date == "2024-01-01T00:00:00Z"


class TestListHitObservations:
    """Tests for list_hit_observations()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_list_empty(self, temp_state_dir):
        """Listing observations from empty database returns empty list."""
        obs = list_hit_observations(state_dir=temp_state_dir)
        assert obs == []

    def test_list_all_observations(self, temp_state_dir):
        """Listing without filters returns all observations."""
        hit1, _ = register_hit("AAAA", "2024-01-01T00:00:00Z", temp_state_dir)
        hit2, _ = register_hit("CCCC", "2024-01-01T00:00:00Z", temp_state_dir)

        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit2.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        obs = list_hit_observations(state_dir=temp_state_dir)

        assert len(obs) == 3

    def test_filter_by_hit_key(self, temp_state_dir):
        """Filtering by hit_key returns only matching observations."""
        hit1, _ = register_hit("AAAA", "2024-01-01T00:00:00Z", temp_state_dir)
        hit2, _ = register_hit("CCCC", "2024-01-01T00:00:00Z", temp_state_dir)

        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit2.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        obs = list_hit_observations(hit_key=hit1.hit_key, state_dir=temp_state_dir)

        assert len(obs) == 2
        assert all(o.hit_key == hit1.hit_key for o in obs)

    def test_filter_by_sample_id(self, temp_state_dir):
        """Filtering by sample_id returns only matching observations."""
        hit1, _ = register_hit("AAAA", "2024-01-01T00:00:00Z", temp_state_dir)

        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        obs = list_hit_observations(sample_id="sample_a", state_dir=temp_state_dir)

        assert len(obs) == 1
        assert obs[0].sample_id == "sample_a"

    def test_filter_by_sample_set_id(self, temp_state_dir):
        """Filtering by sample_set_id returns only matching observations."""
        hit1, _ = register_hit("AAAA", "2024-01-01T00:00:00Z", temp_state_dir)

        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit1.hit_key,
            "set_002",
            "sample_a",
            "2024-01-02T00:00:00Z",
            state_dir=temp_state_dir,
        )

        obs = list_hit_observations(sample_set_id="set_001", state_dir=temp_state_dir)

        assert len(obs) == 1
        assert obs[0].sample_set_id == "set_001"


class TestCountFunctions:
    """Tests for count_hits() and count_hit_observations()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_count_hits_empty(self, temp_state_dir):
        """Counting hits in empty database returns 0."""
        assert count_hits(temp_state_dir) == 0

    def test_count_hits(self, temp_state_dir):
        """Counting hits returns correct count."""
        # Use sequences that are NOT reverse complements of each other
        register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        register_hit("AAAGGG", "2024-01-01T00:00:00Z", temp_state_dir)
        register_hit("AAATTT", "2024-01-01T00:00:00Z", temp_state_dir)

        assert count_hits(temp_state_dir) == 3

    def test_count_observations_empty(self, temp_state_dir):
        """Counting observations in empty database returns 0."""
        assert count_hit_observations(temp_state_dir) == 0

    def test_count_observations(self, temp_state_dir):
        """Counting observations returns correct count."""
        hit, _ = register_hit("AAAA", "2024-01-01T00:00:00Z", temp_state_dir)

        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        assert count_hit_observations(temp_state_dir) == 2


class TestGetHitSequence:
    """Tests for get_hit_sequence()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_get_sequence_simple(self, temp_state_dir):
        """Getting sequence from hit returns original sequence."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        recovered = get_hit_sequence(hit)

        assert recovered == seq

    def test_get_sequence_with_n(self, temp_state_dir):
        """Getting sequence preserves N positions."""
        seq = "ACGTNNNACGT"
        hit, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        recovered = get_hit_sequence(hit)

        assert recovered == seq

    def test_get_sequence_long(self, temp_state_dir):
        """Getting long sequence works correctly."""
        seq = "ACGT" * 1000 + "NNN" + "TGCA" * 500
        hit, _ = register_hit(seq, "2024-01-01T00:00:00Z", temp_state_dir)

        recovered = get_hit_sequence(hit)

        assert recovered == seq


class TestListHitsWithObservations:
    """Tests for list_hits_with_observations()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_empty_database(self, temp_state_dir):
        """Empty database returns empty list."""
        result = list_hits_with_observations(state_dir=temp_state_dir)
        assert result == []

    def test_returns_joined_data(self, temp_state_dir):
        """Returns hits joined with their observations."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        result = list_hits_with_observations(state_dir=temp_state_dir)

        assert len(result) == 1
        returned_hit, returned_obs = result[0]
        assert returned_hit.hit_key == hit.hit_key
        assert returned_obs.sample_id == "sample_a"

    def test_one_hit_multiple_observations(self, temp_state_dir):
        """Hit observed in multiple samples appears multiple times."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        result = list_hits_with_observations(state_dir=temp_state_dir)

        assert len(result) == 2
        sample_ids = {obs.sample_id for _, obs in result}
        assert sample_ids == {"sample_a", "sample_b"}

    def test_respects_limit(self, temp_state_dir):
        """Limit parameter restricts number of rows."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_c",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        result = list_hits_with_observations(limit=2, state_dir=temp_state_dir)

        assert len(result) == 2


class TestDeleteObservationsForSampleSet:
    """Tests for delete_observations_for_sample_set()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_deletes_observations(self, temp_state_dir):
        """Deletes all observations for the sample set."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        deleted = delete_observations_for_sample_set("set_001", temp_state_dir)

        assert deleted == 2
        assert count_hit_observations(temp_state_dir) == 0

    def test_only_deletes_matching_sample_set(self, temp_state_dir):
        """Only deletes observations for the specified sample set."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit.hit_key,
            "set_002",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        deleted = delete_observations_for_sample_set("set_001", temp_state_dir)

        assert deleted == 1
        assert count_hit_observations(temp_state_dir) == 1

    def test_returns_zero_for_nonexistent(self, temp_state_dir):
        """Returns 0 when no observations match."""
        deleted = delete_observations_for_sample_set("nonexistent", temp_state_dir)
        assert deleted == 0


class TestDeleteOrphanedHits:
    """Tests for delete_orphaned_hits()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_deletes_orphaned_hits(self, temp_state_dir):
        """Deletes hits with no observations."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        # Delete the observation, leaving the hit orphaned
        delete_observations_for_sample_set("set_001", temp_state_dir)
        assert count_hits(temp_state_dir) == 1

        deleted = delete_orphaned_hits(temp_state_dir)

        assert deleted == 1
        assert count_hits(temp_state_dir) == 0

    def test_preserves_hits_with_observations(self, temp_state_dir):
        """Does not delete hits that still have observations."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        deleted = delete_orphaned_hits(temp_state_dir)

        assert deleted == 0
        assert count_hits(temp_state_dir) == 1


class TestPruneHitsForSampleSet:
    """Tests for prune_hits_for_sample_set()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_deletes_observations_and_orphaned_hits(self, temp_state_dir):
        """Deletes observations and cleans up orphaned hits atomically."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        obs_deleted, hits_deleted = prune_hits_for_sample_set("set_001", temp_state_dir)

        assert obs_deleted == 1
        assert hits_deleted == 1
        assert count_hit_observations(temp_state_dir) == 0
        assert count_hits(temp_state_dir) == 0

    def test_preserves_shared_hits(self, temp_state_dir):
        """Hits observed in other sample sets are preserved."""
        hit, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        record_hit_observation(
            hit.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit.hit_key,
            "set_002",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        obs_deleted, hits_deleted = prune_hits_for_sample_set("set_001", temp_state_dir)

        assert obs_deleted == 1
        assert hits_deleted == 0  # Hit still has observation in set_002
        assert count_hit_observations(temp_state_dir) == 1
        assert count_hits(temp_state_dir) == 1

    def test_handles_multiple_hits(self, temp_state_dir):
        """Correctly handles multiple hits, some orphaned some not."""
        hit1, _ = register_hit("AAACCC", "2024-01-01T00:00:00Z", temp_state_dir)
        hit2, _ = register_hit("AAAGGG", "2024-01-01T00:00:00Z", temp_state_dir)

        # hit1 observed in both sets, hit2 only in set_001
        record_hit_observation(
            hit1.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit1.hit_key,
            "set_002",
            "sample_b",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )
        record_hit_observation(
            hit2.hit_key,
            "set_001",
            "sample_a",
            "2024-01-01T00:00:00Z",
            state_dir=temp_state_dir,
        )

        obs_deleted, hits_deleted = prune_hits_for_sample_set("set_001", temp_state_dir)

        assert obs_deleted == 2  # Both observations in set_001
        assert hits_deleted == 1  # Only hit2 becomes orphaned
        assert count_hit_observations(temp_state_dir) == 1
        assert count_hits(temp_state_dir) == 1


# =============================================================================
# Statistics / Aggregation Tests
# =============================================================================


class TestGetHitStats:
    """Tests for get_hit_stats()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

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

    def test_single_hit_no_observations(self, temp_state_dir):
        """Single hit with no observations."""
        register_hit("ACGTACGT", "2024-01-15", temp_state_dir)

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 1
        assert stats.total_observations == 0
        assert stats.unique_samples == 0
        assert stats.unique_runs == 0
        assert stats.length_min == 8
        assert stats.length_max == 8
        assert stats.length_median == 8.0
        assert stats.gc_min == 0.5
        assert stats.gc_max == 0.5
        assert stats.gc_median == 0.5
        assert stats.date_first == "2024-01-15"
        assert stats.date_last == "2024-01-15"

    def test_single_hit_with_observations(self, temp_state_dir):
        """Single hit with observations in multiple samples."""
        hit, _ = register_hit("ACGTACGT", "2024-01-15", temp_state_dir)
        record_hit_observation(
            hit.hit_key, "run_001", "sample_a", "2024-01-15", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_001", "sample_b", "2024-01-15", state_dir=temp_state_dir
        )

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 1
        assert stats.total_observations == 2
        assert stats.unique_samples == 2
        assert stats.unique_runs == 1

    def test_multiple_hits_length_distribution(self, temp_state_dir):
        """Multiple hits with varying lengths."""
        # Lengths: 8, 12, 16 -> median = 12
        register_hit("ACGTACGT", "2024-01-01", temp_state_dir)  # 8 bp
        register_hit("ACGTACGTACGT", "2024-01-02", temp_state_dir)  # 12 bp
        register_hit("ACGTACGTACGTACGT", "2024-01-03", temp_state_dir)  # 16 bp

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 3
        assert stats.length_min == 8
        assert stats.length_max == 16
        assert stats.length_median == 12.0

    def test_multiple_hits_gc_distribution(self, temp_state_dir):
        """Multiple hits with varying GC content."""
        # GC contents: 0.0, 0.5, 1.0 -> median = 0.5
        register_hit("AAAAAAAA", "2024-01-01", temp_state_dir)  # 0% GC
        register_hit("ACGTACGT", "2024-01-02", temp_state_dir)  # 50% GC
        register_hit("GCGCGCGC", "2024-01-03", temp_state_dir)  # 100% GC

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 3
        assert stats.gc_min == 0.0
        assert stats.gc_max == 1.0
        assert stats.gc_median == 0.5

    def test_date_range(self, temp_state_dir):
        """Date range spans first and last seen dates."""
        register_hit("AAACCC", "2024-01-15", temp_state_dir)
        register_hit("AAAGGG", "2024-06-20", temp_state_dir)
        register_hit("AAATTT", "2024-03-10", temp_state_dir)

        stats = get_hit_stats(temp_state_dir)

        assert stats.date_first == "2024-01-15"
        assert stats.date_last == "2024-06-20"

    def test_multiple_runs_and_samples(self, temp_state_dir):
        """Correctly counts unique runs and samples."""
        hit1, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)
        hit2, _ = register_hit("AAAGGG", "2024-01-01", temp_state_dir)

        # hit1 in run_001: sample_a, sample_b
        record_hit_observation(
            hit1.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit1.hit_key, "run_001", "sample_b", "2024-01-01", state_dir=temp_state_dir
        )
        # hit1 in run_002: sample_a (same sample, different run)
        record_hit_observation(
            hit1.hit_key, "run_002", "sample_a", "2024-01-02", state_dir=temp_state_dir
        )
        # hit2 in run_001: sample_c
        record_hit_observation(
            hit2.hit_key, "run_001", "sample_c", "2024-01-01", state_dir=temp_state_dir
        )

        stats = get_hit_stats(temp_state_dir)

        assert stats.total_hits == 2
        assert stats.total_observations == 4
        assert stats.unique_samples == 3  # sample_a, sample_b, sample_c
        assert stats.unique_runs == 2  # run_001, run_002

    def test_median_even_count(self, temp_state_dir):
        """Median calculation for even number of values."""
        # Lengths: 8, 12, 16, 20 -> median = (12 + 16) / 2 = 14
        register_hit("ACGTACGT", "2024-01-01", temp_state_dir)  # 8 bp
        register_hit("ACGTACGTACGT", "2024-01-02", temp_state_dir)  # 12 bp
        register_hit("ACGTACGTACGTACGT", "2024-01-03", temp_state_dir)  # 16 bp
        register_hit("ACGTACGTACGTACGTACGT", "2024-01-04", temp_state_dir)  # 20 bp

        stats = get_hit_stats(temp_state_dir)

        assert stats.length_median == 14.0


class TestGetRecurringHits:
    """Tests for get_recurring_hits()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_empty_database(self, temp_state_dir):
        """Empty database returns empty list."""
        result = get_recurring_hits(state_dir=temp_state_dir)
        assert result == []

    def test_no_recurring_hits(self, temp_state_dir):
        """Hits with only one sample each don't appear."""
        hit1, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)
        hit2, _ = register_hit("AAAGGG", "2024-01-01", temp_state_dir)

        record_hit_observation(
            hit1.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit2.hit_key, "run_001", "sample_b", "2024-01-01", state_dir=temp_state_dir
        )

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)
        assert result == []

    def test_finds_recurring_hit(self, temp_state_dir):
        """Hit appearing in multiple samples is returned."""
        hit, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)

        record_hit_observation(
            hit.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_001", "sample_b", "2024-01-02", state_dir=temp_state_dir
        )

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].hit_key == hit.hit_key
        assert result[0].sample_count == 2
        assert result[0].run_count == 1
        assert result[0].observation_count == 2

    def test_min_samples_filter(self, temp_state_dir):
        """min_samples parameter filters correctly."""
        hit, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)

        # Observed in 3 samples
        record_hit_observation(
            hit.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_001", "sample_b", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_001", "sample_c", "2024-01-01", state_dir=temp_state_dir
        )

        # Should appear with min_samples=2 or 3, not 4
        assert len(get_recurring_hits(min_samples=2, state_dir=temp_state_dir)) == 1
        assert len(get_recurring_hits(min_samples=3, state_dir=temp_state_dir)) == 1
        assert len(get_recurring_hits(min_samples=4, state_dir=temp_state_dir)) == 0

    def test_min_runs_filter(self, temp_state_dir):
        """min_runs parameter filters correctly."""
        hit, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)

        # Observed in 2 samples across 2 runs
        record_hit_observation(
            hit.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_002", "sample_b", "2024-01-02", state_dir=temp_state_dir
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
        hit1, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)
        hit2, _ = register_hit("AAAGGG", "2024-01-01", temp_state_dir)

        # Both hits in 2 samples
        for hit in [hit1, hit2]:
            record_hit_observation(
                hit.hit_key,
                "run_001",
                "sample_a",
                "2024-01-01",
                state_dir=temp_state_dir,
            )
            record_hit_observation(
                hit.hit_key,
                "run_001",
                "sample_b",
                "2024-01-01",
                state_dir=temp_state_dir,
            )

        assert len(get_recurring_hits(min_samples=2, state_dir=temp_state_dir)) == 2
        assert (
            len(get_recurring_hits(min_samples=2, limit=1, state_dir=temp_state_dir))
            == 1
        )

    def test_ordered_by_sample_count(self, temp_state_dir):
        """Results ordered by sample_count descending."""
        hit1, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)
        hit2, _ = register_hit("AAAGGG", "2024-01-01", temp_state_dir)

        # hit1 in 2 samples, hit2 in 3 samples
        record_hit_observation(
            hit1.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit1.hit_key, "run_001", "sample_b", "2024-01-01", state_dir=temp_state_dir
        )

        record_hit_observation(
            hit2.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit2.hit_key, "run_001", "sample_b", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit2.hit_key, "run_001", "sample_c", "2024-01-01", state_dir=temp_state_dir
        )

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 2
        assert result[0].hit_key == hit2.hit_key  # 3 samples first
        assert result[0].sample_count == 3
        assert result[1].hit_key == hit1.hit_key  # 2 samples second
        assert result[1].sample_count == 2

    def test_last_seen_date(self, temp_state_dir):
        """last_seen_date is the most recent observation."""
        hit, _ = register_hit("AAACCC", "2024-01-01", temp_state_dir)

        record_hit_observation(
            hit.hit_key, "run_001", "sample_a", "2024-01-15", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_001", "sample_b", "2024-06-20", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_002", "sample_c", "2024-03-10", state_dir=temp_state_dir
        )

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].first_seen_date == "2024-01-01"
        assert result[0].last_seen_date == "2024-06-20"

    def test_returns_correct_metadata(self, temp_state_dir):
        """Returned RecurringHit has correct hit metadata."""
        hit, _ = register_hit("GCGCGCGC", "2024-01-01", temp_state_dir)  # 100% GC, 8 bp

        record_hit_observation(
            hit.hit_key, "run_001", "sample_a", "2024-01-01", state_dir=temp_state_dir
        )
        record_hit_observation(
            hit.hit_key, "run_001", "sample_b", "2024-01-01", state_dir=temp_state_dir
        )

        result = get_recurring_hits(min_samples=2, state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].sequence_length == 8
        assert result[0].gc_content == 1.0


class TestGetDiscoveryTimeline:
    """Tests for get_discovery_timeline()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_empty_database(self, temp_state_dir):
        """Empty database returns empty list."""
        result = get_discovery_timeline(state_dir=temp_state_dir)
        assert result == []

    def test_single_hit(self, temp_state_dir):
        """Single hit returns one bucket."""
        register_hit("AAACCC", "2024-01-15", temp_state_dir)

        result = get_discovery_timeline(granularity="month", state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].period == "2024-01"
        assert result[0].new_hits == 1

    def test_multiple_hits_same_month(self, temp_state_dir):
        """Multiple hits in same month are grouped."""
        register_hit("AAACCC", "2024-01-10", temp_state_dir)
        register_hit("AAAGGG", "2024-01-20", temp_state_dir)
        register_hit("AAATTT", "2024-01-25", temp_state_dir)

        result = get_discovery_timeline(granularity="month", state_dir=temp_state_dir)

        assert len(result) == 1
        assert result[0].period == "2024-01"
        assert result[0].new_hits == 3

    def test_multiple_months(self, temp_state_dir):
        """Hits across months create separate buckets."""
        register_hit("AAACCC", "2024-01-15", temp_state_dir)
        register_hit("AAAGGG", "2024-01-20", temp_state_dir)
        register_hit("AAATTT", "2024-03-10", temp_state_dir)
        register_hit("CCCGGG", "2024-06-01", temp_state_dir)

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
        register_hit("AAAGGG", "2024-06-01", temp_state_dir)
        register_hit("AAACCC", "2024-01-15", temp_state_dir)
        register_hit("AAATTT", "2024-03-10", temp_state_dir)

        result = get_discovery_timeline(granularity="month", state_dir=temp_state_dir)

        periods = [b.period for b in result]
        assert periods == ["2024-01", "2024-03", "2024-06"]

    def test_granularity_day(self, temp_state_dir):
        """Day granularity groups by exact date."""
        register_hit("AAACCC", "2024-01-15", temp_state_dir)
        register_hit("AAAGGG", "2024-01-15", temp_state_dir)
        register_hit("AAATTT", "2024-01-16", temp_state_dir)

        result = get_discovery_timeline(granularity="day", state_dir=temp_state_dir)

        assert len(result) == 2
        assert result[0].period == "2024-01-15"
        assert result[0].new_hits == 2
        assert result[1].period == "2024-01-16"
        assert result[1].new_hits == 1

    def test_granularity_year(self, temp_state_dir):
        """Year granularity groups by year."""
        register_hit("AAACCC", "2023-06-15", temp_state_dir)
        register_hit("AAAGGG", "2024-01-15", temp_state_dir)
        register_hit("AAATTT", "2024-12-31", temp_state_dir)

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
        register_hit("AAACCC", "2024-01-01", temp_state_dir)
        register_hit("AAAGGG", "2024-01-07", temp_state_dir)  # Same week
        register_hit("AAATTT", "2024-01-08", temp_state_dir)  # Next week

        result = get_discovery_timeline(granularity="week", state_dir=temp_state_dir)

        assert len(result) == 2
        # First bucket should have 2 hits, second should have 1
        assert result[0].new_hits == 2
        assert result[1].new_hits == 1


class TestLookupHit:
    """Tests for lookup_hit()."""

    @pytest.fixture
    def temp_state_dir(self):
        """Create a temporary directory for state database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_lookup_by_key(self, temp_state_dir):
        """Lookup by hit key returns the hit."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", temp_state_dir)

        result = lookup_hit(hit.hit_key, state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == hit.hit_key

    def test_lookup_by_key_uppercase(self, temp_state_dir):
        """Lookup by uppercase hit key works."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", temp_state_dir)

        result = lookup_hit(hit.hit_key.upper(), state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == hit.hit_key

    def test_lookup_by_sequence(self, temp_state_dir):
        """Lookup by sequence returns the hit."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", temp_state_dir)

        result = lookup_hit(seq, state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == hit.hit_key

    def test_lookup_by_sequence_lowercase(self, temp_state_dir):
        """Lookup by lowercase sequence works."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", temp_state_dir)

        result = lookup_hit(seq.lower(), state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == hit.hit_key

    def test_lookup_by_reverse_complement(self, temp_state_dir):
        """Lookup by reverse complement finds the same hit."""
        seq = "AAACCCGGG"
        hit, _ = register_hit(seq, "2024-01-01", temp_state_dir)

        revcomp = reverse_complement(seq)
        result = lookup_hit(revcomp, state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == hit.hit_key

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
        hit, _ = register_hit(seq, "2024-01-01", temp_state_dir)

        result = lookup_hit(f"  {seq}  \n", state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == hit.hit_key

    def test_lookup_key_with_whitespace(self, temp_state_dir):
        """Lookup key with whitespace is stripped."""
        seq = "ACGTACGTACGT"
        hit, _ = register_hit(seq, "2024-01-01", temp_state_dir)

        result = lookup_hit(f"  {hit.hit_key}  ", state_dir=temp_state_dir)

        assert result is not None
        assert result.hit_key == hit.hit_key
