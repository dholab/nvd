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
    get_hit,
    get_hit_sequence,
    list_hit_observations,
    list_hits,
    list_hits_with_observations,
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
