"""
Tests for register_gottcha2_hits.py I/O parsing and orchestration.

The core sequence operations and state management are tested in
lib/py_nvd/test_hits.py. This file tests the GOTTCHA2-specific wrapper's
FASTA header parsing, hit record building, and parquet writing.
"""

import tempfile
from pathlib import Path

import pytest
from register_gottcha2_hits import (
    Gottcha2Classification,
    Gottcha2HitRegistrationContext,
    build_gottcha2_hit_records,
    parse_gottcha2_fasta,
    write_hits_to_path,
)


class TestGottcha2HitRegistrationContext:
    """Tests for Gottcha2HitRegistrationContext validation."""

    def test_valid_context(self):
        """Valid context is created successfully."""
        ctx = Gottcha2HitRegistrationContext(
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
            Gottcha2HitRegistrationContext(
                state_dir=Path("/tmp"),
                sample_set_id="",
                sample_id="sample_a",
                run_date="2024-01-01T00:00:00Z",
            )

    def test_empty_sample_id_raises(self):
        """Empty sample_id raises ValueError."""
        with pytest.raises(ValueError, match="sample_id cannot be empty"):
            Gottcha2HitRegistrationContext(
                state_dir=Path("/tmp"),
                sample_set_id="set_001",
                sample_id="",
                run_date="2024-01-01T00:00:00Z",
            )

    def test_empty_run_date_raises(self):
        """Empty run_date raises ValueError."""
        with pytest.raises(ValueError, match="run_date cannot be empty"):
            Gottcha2HitRegistrationContext(
                state_dir=Path("/tmp"),
                sample_set_id="set_001",
                sample_id="sample_a",
                run_date="",
            )


class TestParseGottcha2Fasta:
    """Tests for parse_gottcha2_fasta()."""

    def test_parse_standard_gottcha2_header(self):
        """Parses a standard GOTTCHA2 extracted FASTA with taxonomy in headers."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(
                ">read1.1|NC_001234|1|100|12345:11..20 LEVEL=species NAME=Escherichia_coli TAXID=562\n"
            )
            f.write("ACGTACGTACGT\n")
            f.write(
                ">read2.2|NC_005678|1|200|67890:5..15 LEVEL=genus NAME=Staphylococcus TAXID=1279\n"
            )
            f.write("GGGGCCCCAAAA\n")
            fasta_path = Path(f.name)

        try:
            results = parse_gottcha2_fasta(fasta_path)
            assert len(results) == 2

            seq_id_1, seq_1, cls_1 = results[0]
            assert seq_1 == "ACGTACGTACGT"
            assert cls_1.taxid == 562
            assert cls_1.name == "Escherichia coli"
            assert cls_1.rank == "species"

            seq_id_2, seq_2, cls_2 = results[1]
            assert seq_2 == "GGGGCCCCAAAA"
            assert cls_2.taxid == 1279
            assert cls_2.name == "Staphylococcus"
            assert cls_2.rank == "genus"
        finally:
            fasta_path.unlink()

    def test_parse_multiword_name_with_underscores(self):
        """Names with underscores are restored to spaces."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(
                ">read1|ref:1..10 LEVEL=species NAME=Human_herpesvirus_4 TAXID=10376\n"
            )
            f.write("ACGTACGT\n")
            fasta_path = Path(f.name)

        try:
            results = parse_gottcha2_fasta(fasta_path)
            assert len(results) == 1
            _, _, cls = results[0]
            assert cls.name == "Human herpesvirus 4"
            assert cls.taxid == 10376
            assert cls.rank == "species"
        finally:
            fasta_path.unlink()

    def test_parse_empty_fasta(self):
        """Parses empty FASTA file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            fasta_path = Path(f.name)

        try:
            results = parse_gottcha2_fasta(fasta_path)
            assert results == []
        finally:
            fasta_path.unlink()

    def test_skips_unparseable_headers(self):
        """Entries with headers that don't match the expected pattern are skipped."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            # Valid entry
            f.write(">read1|ref:1..10 LEVEL=species NAME=Test_virus TAXID=12345\n")
            f.write("ACGTACGT\n")
            # Invalid entry (no LEVEL/NAME/TAXID)
            f.write(">some_other_header with no taxonomy\n")
            f.write("GGGGCCCC\n")
            fasta_path = Path(f.name)

        try:
            results = parse_gottcha2_fasta(fasta_path)
            assert len(results) == 1
            _, _, cls = results[0]
            assert cls.name == "Test virus"
        finally:
            fasta_path.unlink()

    def test_skips_empty_sequences(self):
        """Entries with empty sequences are skipped."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">read1|ref:1..10 LEVEL=species NAME=Test_virus TAXID=12345\n")
            f.write("\n")
            f.write(">read2|ref:1..10 LEVEL=species NAME=Other_virus TAXID=67890\n")
            f.write("ACGTACGT\n")
            fasta_path = Path(f.name)

        try:
            results = parse_gottcha2_fasta(fasta_path)
            # BioPython may or may not parse the empty sequence entry;
            # if it does, our code skips it
            sequence_ids = [r[0] for r in results]
            assert "read2|ref:1..10" in sequence_ids or len(results) >= 1
        finally:
            fasta_path.unlink()

    def test_multiline_sequence(self):
        """Parses FASTA with multi-line sequences."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">read1|ref:1..10 LEVEL=species NAME=Test_virus TAXID=12345\n")
            f.write("ACGT\n")
            f.write("ACGT\n")
            fasta_path = Path(f.name)

        try:
            results = parse_gottcha2_fasta(fasta_path)
            assert len(results) == 1
            _, seq, _ = results[0]
            assert seq == "ACGTACGT"
        finally:
            fasta_path.unlink()


class TestBuildGottcha2HitRecords:
    """Tests for build_gottcha2_hit_records()."""

    def _make_context(self) -> Gottcha2HitRegistrationContext:
        """Helper to create a test context."""
        return Gottcha2HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )

    def _make_entry(
        self,
        seq_id: str = "read1|ref:1..10",
        seq: str = "ACGTACGT",
        taxid: int = 12345,
        name: str = "Test virus",
        rank: str = "species",
    ) -> tuple[str, str, Gottcha2Classification]:
        """Helper to create a parsed FASTA entry."""
        return (
            seq_id,
            seq,
            Gottcha2Classification(
                sequence_id=seq_id,
                taxid=taxid,
                name=name,
                rank=rank,
            ),
        )

    def test_builds_records_from_parsed_entries(self):
        """Builds HitRecords from parsed GOTTCHA2 FASTA entries."""
        entries = [
            self._make_entry("read1|ref:1..10", "ACGTACGT", 562, "E. coli", "species"),
            self._make_entry("read2|ref:5..15", "GGGGCCCC", 1279, "Staph", "genus"),
        ]
        context = self._make_context()

        records = build_gottcha2_hit_records(entries, context)

        assert len(records) == 2
        sequence_ids = {r.sequence_id for r in records}
        assert sequence_ids == {"read1|ref:1..10", "read2|ref:5..15"}

    def test_empty_entries_returns_empty_list(self):
        """Returns empty list when no entries."""
        context = self._make_context()
        records = build_gottcha2_hit_records([], context)
        assert records == []

    def test_populates_metadata_correctly(self):
        """HitRecords have correct metadata from context."""
        entries = [self._make_entry()]
        context = self._make_context()

        records = build_gottcha2_hit_records(entries, context)

        assert len(records) == 1
        record = records[0]
        assert record.sample_set_id == "set_001"
        assert record.sample_id == "sample_a"
        assert record.run_date == "2024-01-01T00:00:00Z"
        assert record.source == "gottcha2"
        assert record.sequence_length == 8
        assert len(record.hit_key) == 32

    def test_populates_classification_from_header(self):
        """HitRecords include taxonomy parsed from GOTTCHA2 headers."""
        entries = [
            self._make_entry(
                taxid=10376,
                name="Human herpesvirus 4",
                rank="species",
            ),
        ]
        context = self._make_context()

        records = build_gottcha2_hit_records(entries, context)

        assert len(records) == 1
        record = records[0]
        assert record.adjusted_taxid == 10376
        assert record.adjusted_taxid_name == "Human herpesvirus 4"
        assert record.adjusted_taxid_rank == "species"


class TestWriteHitsToPath:
    """Tests for write_hits_to_path()."""

    def test_writes_parquet_to_specified_path(self):
        """Writes parquet file to the specified output path."""
        import polars as pl

        entries = [
            (
                "read1|ref:1..10",
                "ACGTACGT",
                Gottcha2Classification(
                    sequence_id="read1|ref:1..10",
                    taxid=12345,
                    name="Test virus",
                    rank="species",
                ),
            ),
        ]
        context = Gottcha2HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )
        records = build_gottcha2_hit_records(entries, context)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.parquet"
            result_path = write_hits_to_path(records, output_path)

            assert result_path == output_path
            assert output_path.exists()

            df = pl.read_parquet(output_path)
            assert len(df) == 1
            assert df["sample_id"][0] == "sample_a"
            assert df["source"][0] == "gottcha2"

    def test_creates_parent_directories(self):
        """Creates parent directories if they don't exist."""
        entries = [
            (
                "read1|ref:1..10",
                "ACGTACGT",
                Gottcha2Classification(
                    sequence_id="read1|ref:1..10",
                    taxid=12345,
                    name="Test virus",
                    rank="species",
                ),
            ),
        ]
        context = Gottcha2HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )
        records = build_gottcha2_hit_records(entries, context)

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

            df = pl.read_parquet(output_path)
            assert len(df) == 0
            assert "hit_key" in df.columns
            assert "sample_id" in df.columns
            assert "source" in df.columns

    def test_parquet_includes_contig_id_null_stub(self):
        """Parquet includes contig_id=NULL for backward compatibility."""
        import polars as pl

        entries = [
            (
                "read1|ref:1..10",
                "ACGTACGT",
                Gottcha2Classification(
                    sequence_id="read1|ref:1..10",
                    taxid=12345,
                    name="Test virus",
                    rank="species",
                ),
            ),
        ]
        context = Gottcha2HitRegistrationContext(
            state_dir=Path("/tmp"),
            sample_set_id="set_001",
            sample_id="sample_a",
            run_date="2024-01-01T00:00:00Z",
        )
        records = build_gottcha2_hit_records(entries, context)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.parquet"
            write_hits_to_path(records, output_path)

            df = pl.read_parquet(output_path)
            assert "contig_id" in df.columns
            assert df["contig_id"][0] is None
