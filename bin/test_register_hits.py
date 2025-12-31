"""
Tests for register_hits.py I/O parsing and orchestration.

The core sequence operations and state management are tested in
lib/py_nvd/test_hits.py. This file tests the thin wrapper's I/O functions.
"""

import tempfile
from pathlib import Path

import pytest

from register_hits import (
    HitRegistrationContext,
    parse_blast_contig_ids,
    parse_fasta,
)


class TestHitRegistrationContext:
    """Tests for HitRegistrationContext validation."""

    def test_valid_context(self):
        """Valid context is created successfully."""
        ctx = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )
        assert ctx.sample_set_id == "set_001"
        assert ctx.sample_id == "sample_a"

    def test_empty_sample_set_id_raises(self):
        """Empty sample_set_id raises ValueError."""
        with pytest.raises(ValueError, match="sample_set_id cannot be empty"):
            HitRegistrationContext(
                state_dir=Path("/tmp"),
                sample_set_id="",
                sample_id="sample_a",
                run_date="2024-01-01T00:00:00Z",
            )

    def test_empty_sample_id_raises(self):
        """Empty sample_id raises ValueError."""
        with pytest.raises(ValueError, match="sample_id cannot be empty"):
            HitRegistrationContext(
                state_dir=Path("/tmp"),
                sample_set_id="set_001",
                sample_id="",
                run_date="2024-01-01T00:00:00Z",
            )

    def test_empty_run_date_raises(self):
        """Empty run_date raises ValueError."""
        with pytest.raises(ValueError, match="run_date cannot be empty"):
            HitRegistrationContext(
                state_dir=Path("/tmp"),
                sample_set_id="set_001",
                sample_id="sample_a",
                run_date="",
            )


class TestParseFasta:
    """Tests for parse_fasta()."""

    def test_parse_simple_fasta(self):
        """Parses a simple FASTA file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">contig1\n")
            f.write("ACGTACGT\n")
            f.write(">contig2\n")
            f.write("GGGGCCCC\n")
            fasta_path = Path(f.name)

        try:
            contigs = parse_fasta(fasta_path)
            assert len(contigs) == 2
            assert contigs["contig1"] == "ACGTACGT"
            assert contigs["contig2"] == "GGGGCCCC"
        finally:
            fasta_path.unlink()

    def test_parse_multiline_sequence(self):
        """Parses FASTA with multi-line sequences."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">contig1\n")
            f.write("ACGT\n")
            f.write("ACGT\n")
            fasta_path = Path(f.name)

        try:
            contigs = parse_fasta(fasta_path)
            assert contigs["contig1"] == "ACGTACGT"
        finally:
            fasta_path.unlink()

    def test_parse_empty_fasta(self):
        """Parses empty FASTA file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            fasta_path = Path(f.name)

        try:
            contigs = parse_fasta(fasta_path)
            assert contigs == {}
        finally:
            fasta_path.unlink()


class TestParseBlastContigIds:
    """Tests for parse_blast_contig_ids()."""

    def test_parse_blast_results(self):
        """Parses BLAST results TSV and extracts qseqid."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample\ttask\tqseqid\tsseqid\tpident\n")
            f.write("s1\tblast\tNODE_1\tNC_001\t99.5\n")
            f.write("s1\tblast\tNODE_2\tNC_002\t98.0\n")
            f.write("s1\tblast\tNODE_1\tNC_003\t97.0\n")  # duplicate qseqid
            tsv_path = Path(f.name)

        try:
            contig_ids = parse_blast_contig_ids(tsv_path)
            assert contig_ids == {"NODE_1", "NODE_2"}
        finally:
            tsv_path.unlink()

    def test_missing_qseqid_column_raises(self):
        """Raises ValueError if qseqid column is missing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample\ttask\tsseqid\n")
            f.write("s1\tblast\tNC_001\n")
            tsv_path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="Expected column 'qseqid' not found"):
                parse_blast_contig_ids(tsv_path)
        finally:
            tsv_path.unlink()

    def test_empty_results(self):
        """Parses TSV with header only (no data rows)."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample\ttask\tqseqid\tsseqid\n")
            tsv_path = Path(f.name)

        try:
            contig_ids = parse_blast_contig_ids(tsv_path)
            assert contig_ids == set()
        finally:
            tsv_path.unlink()
