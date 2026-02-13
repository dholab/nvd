"""
Tests for register_hits.py I/O parsing and orchestration.

The core sequence operations and state management are tested in
lib/py_nvd/test_hits.py. This file tests the thin wrapper's I/O functions.
"""

import tempfile
from pathlib import Path

import pytest
from register_hits import (
    ContigClassification,
    HitRegistrationContext,
    build_hit_records,
    parse_blast_classifications,
    parse_fasta,
    write_hits_to_path,
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


class TestParseBlastClassifications:
    """Tests for parse_blast_classifications()."""

    def test_parse_blast_results_with_classifications(self):
        """Parses BLAST results TSV and extracts classifications per contig."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(
                "sample\ttask\tqseqid\tsseqid\tpident\tadjusted_taxid\tadjusted_taxid_name\tadjusted_taxid_rank\n",
            )
            f.write("s1\tblast\tNODE_1\tNC_001\t99.5\t12345\tHuman virus\tspecies\n")
            f.write("s1\tblast\tNODE_2\tNC_002\t98.0\t67890\tOther virus\tgenus\n")
            f.write(
                "s1\tblast\tNODE_1\tNC_003\t97.0\t12345\tHuman virus\tspecies\n",
            )  # duplicate qseqid
            tsv_path = Path(f.name)

        try:
            classifications = parse_blast_classifications(tsv_path)
            assert set(classifications.keys()) == {"NODE_1", "NODE_2"}
            assert classifications["NODE_1"].adjusted_taxid == 12345
            assert classifications["NODE_1"].adjusted_taxid_name == "Human virus"
            assert classifications["NODE_1"].adjusted_taxid_rank == "species"
            assert classifications["NODE_2"].adjusted_taxid == 67890
        finally:
            tsv_path.unlink()

    def test_missing_qseqid_column_raises(self):
        """Raises ValueError if qseqid column is missing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("sample\ttask\tsseqid\tadjusted_taxid\n")
            f.write("s1\tblast\tNC_001\t12345\n")
            tsv_path = Path(f.name)

        try:
            with pytest.raises(ValueError, match="Expected column 'qseqid' not found"):
                parse_blast_classifications(tsv_path)
        finally:
            tsv_path.unlink()

    def test_empty_results(self):
        """Parses TSV with header only (no data rows)."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(
                "sample\ttask\tqseqid\tsseqid\tadjusted_taxid\tadjusted_taxid_name\tadjusted_taxid_rank\n",
            )
            tsv_path = Path(f.name)

        try:
            classifications = parse_blast_classifications(tsv_path)
            assert classifications == {}
        finally:
            tsv_path.unlink()

    def test_handles_null_classification_values(self):
        """Handles rows with null/empty classification values."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write(
                "sample\ttask\tqseqid\tsseqid\tpident\tadjusted_taxid\tadjusted_taxid_name\tadjusted_taxid_rank\n",
            )
            f.write("s1\tblast\tNODE_1\tNC_001\t99.5\t\t\t\n")
            tsv_path = Path(f.name)

        try:
            classifications = parse_blast_classifications(tsv_path)
            assert "NODE_1" in classifications
            assert classifications["NODE_1"].adjusted_taxid is None
            assert classifications["NODE_1"].adjusted_taxid_name is None
        finally:
            tsv_path.unlink()


def _make_classification(
    contig_id: str,
    taxid: int | None = 12345,
    name: str | None = "Test virus",
    rank: str | None = "species",
) -> ContigClassification:
    """Helper to create ContigClassification for tests."""
    return ContigClassification(
        contig_id=contig_id,
        adjusted_taxid=taxid,
        adjusted_taxid_name=name,
        adjusted_taxid_rank=rank,
    )


class TestBuildHitRecords:
    """Tests for build_hit_records()."""

    def test_builds_records_for_matching_contigs(self):
        """Builds HitRecords for contigs that have BLAST hits."""
        contigs = {"NODE_1": "ACGTACGT", "NODE_2": "GGGGCCCC", "NODE_3": "AAAATTTT"}
        classifications = {
            "NODE_1": _make_classification("NODE_1"),
            "NODE_3": _make_classification("NODE_3"),
        }
        context = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )

        records = build_hit_records(contigs, classifications, context)

        assert len(records) == 2
        sequence_ids = {r.sequence_id for r in records}
        assert sequence_ids == {"NODE_1", "NODE_3"}

    def test_skips_missing_contigs(self):
        """Skips contig IDs that aren't in the FASTA."""
        contigs = {"NODE_1": "ACGTACGT"}
        classifications = {
            "NODE_1": _make_classification("NODE_1"),
            "NODE_MISSING": _make_classification("NODE_MISSING"),
        }
        context = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )

        records = build_hit_records(contigs, classifications, context)

        assert len(records) == 1
        assert records[0].sequence_id == "NODE_1"

    def test_empty_hits_returns_empty_list(self):
        """Returns empty list when no contigs have hits."""
        contigs = {"NODE_1": "ACGTACGT"}
        classifications: dict[str, ContigClassification] = {}
        context = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )

        records = build_hit_records(contigs, classifications, context)

        assert records == []

    def test_populates_metadata_correctly(self):
        """HitRecords have correct metadata from context."""
        contigs = {"NODE_1": "ACGTACGT"}
        classifications = {"NODE_1": _make_classification("NODE_1")}
        context = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )

        records = build_hit_records(contigs, classifications, context)

        assert len(records) == 1
        record = records[0]
        assert record.sample_set_id == "set_001"
        assert record.sample_id == "sample_a"
        assert record.run_date == "2024-01-01T00:00:00Z"
        assert record.sequence_length == 8
        assert len(record.hit_key) == 32

    def test_populates_classification_fields(self):
        """HitRecords include classification data from ContigClassification."""
        contigs = {"NODE_1": "ACGTACGT"}
        classifications = {
            "NODE_1": ContigClassification(
                contig_id="NODE_1",
                adjusted_taxid=99999,
                adjusted_taxid_name="Special virus",
                adjusted_taxid_rank="genus",
            ),
        }
        context = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )

        records = build_hit_records(contigs, classifications, context)

        assert len(records) == 1
        record = records[0]
        assert record.adjusted_taxid == 99999
        assert record.adjusted_taxid_name == "Special virus"
        assert record.adjusted_taxid_rank == "genus"


class TestWriteHitsToPath:
    """Tests for write_hits_to_path()."""

    def test_writes_parquet_to_specified_path(self):
        """Writes parquet file to the specified output path."""
        import polars as pl

        contigs = {"NODE_1": "ACGTACGT"}
        classifications = {"NODE_1": _make_classification("NODE_1")}
        context = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )
        records = build_hit_records(contigs, classifications, context)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.parquet"
            result_path = write_hits_to_path(records, output_path)

            assert result_path == output_path
            assert output_path.exists()

            # Verify contents
            df = pl.read_parquet(output_path)
            assert len(df) == 1
            assert df["sample_id"][0] == "sample_a"

    def test_creates_parent_directories(self):
        """Creates parent directories if they don't exist."""
        contigs = {"NODE_1": "ACGTACGT"}
        classifications = {"NODE_1": _make_classification("NODE_1")}
        context = HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )
        records = build_hit_records(contigs, classifications, context)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "nested" / "dir" / "output.parquet"
            result_path = write_hits_to_path(records, output_path)

            assert result_path == output_path
            assert output_path.exists()

    def test_writes_empty_parquet_with_schema(self):
        """Writes valid parquet file even with empty records list."""
        import polars as pl

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "empty.parquet"
            result_path = write_hits_to_path([], output_path)

            assert result_path == output_path
            assert output_path.exists()

            # Verify it's a valid parquet with correct schema
            df = pl.read_parquet(output_path)
            assert len(df) == 0
            assert "hit_key" in df.columns
            assert "sample_id" in df.columns
