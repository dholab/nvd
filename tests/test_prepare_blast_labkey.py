"""Behavioral tests for prepare_blast_labkey.py, especially hit-less samples.

A sample with zero BLAST hits (water/negative controls) must still emit a
header-bearing CSV. LABKEY_CONCAT_ALL_SAMPLE_BLAST_RESULTS must read prepared
CSVs under one all-String schema so header-only files can mix with populated ones.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "prepare_blast_labkey.py"

# Columns as emitted by finalize_blast_results.py (ADD_READ_COUNTS_TO_BLAST).
BLAST_COLUMNS = [
    "qseqid",
    "task",
    "sample",
    "rank",
    "pident",
    "mapped_reads",
    "total_reads",
    "blast_db_version",
    "virus_index_version",
    "nextflow_run_id",
]

EXPECTED_HEADER = [
    "experiment",
    "qseqid",
    "blast_task",
    "sample_id",
    "tax_rank",
    "pident",
    "mapped_reads",
    "total_reads",
    "blast_db_version",
    "virus_index_version",
    "nextflow_run_id",
]

DATA_ROW = [
    "contig_1",
    "megablast",
    "Sample_S34",
    "species",
    "99.8",
    "1234",
    "500000",
    "nt_2026_07",
    "not_used",
    "run_abc",
]


def write_tsv(path: Path, rows: list[list[str]]) -> Path:
    """Write a BLAST final TSV with a header and the given data rows."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(BLAST_COLUMNS)
        writer.writerows(rows)
    return path


def run_prepare(tmp_path: Path, tsv: Path, meta: str) -> Path:
    """Run prepare_blast_labkey.py and return the output CSV path."""
    output = tmp_path / f"{meta}_blast_labkey.csv"
    result = subprocess.run(  # noqa: S603
        [
            sys.executable,
            str(SCRIPT),
            "--blast-csv",
            str(tsv),
            "--output",
            str(output),
            "--meta",
            meta,
            "--experiment-id",
            "32495",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, f"script failed: {result.stderr}"
    return output


def test_hitless_sample_writes_header_not_empty_file(tmp_path: Path) -> None:
    """A header-only input must yield a header-only CSV, never a 0-byte file."""
    tsv = write_tsv(tmp_path / "Water_S11_blast.final.tsv", [])
    output = run_prepare(tmp_path, tsv, "Water_S11")

    assert output.stat().st_size > 0, "hit-less sample produced a 0-byte CSV"
    rows = list(csv.reader(output.open()))
    assert rows == [EXPECTED_HEADER]


def test_populated_sample_unchanged(tmp_path: Path) -> None:
    """Samples with hits keep their existing header and rows."""
    tsv = write_tsv(tmp_path / "Sample_S34_blast.final.tsv", [DATA_ROW])
    output = run_prepare(tmp_path, tsv, "Sample_S34")

    rows = list(csv.reader(output.open()))
    assert rows[0] == EXPECTED_HEADER
    assert rows[1] == ["32495", *DATA_ROW]


def test_hitless_and_populated_headers_match(tmp_path: Path) -> None:
    """Both branches must emit identical headers for the concat schema."""
    empty_out = run_prepare(
        tmp_path,
        write_tsv(tmp_path / "Water_S11_blast.final.tsv", []),
        "Water_S11",
    )
    full_out = run_prepare(
        tmp_path,
        write_tsv(tmp_path / "Sample_S34_blast.final.tsv", [DATA_ROW]),
        "Sample_S34",
    )

    empty_header = next(csv.reader(empty_out.open()))
    full_header = next(csv.reader(full_out.open()))
    assert empty_header == full_header


def test_header_only_rows_are_skipped(tmp_path: Path) -> None:
    """A duplicated header smuggled in as data is still dropped."""
    tsv = write_tsv(tmp_path / "Sample_S34_blast.final.tsv", [BLAST_COLUMNS])
    output = run_prepare(tmp_path, tsv, "Sample_S34")

    rows = list(csv.reader(output.open()))
    assert rows == [EXPECTED_HEADER]


def test_headerless_input_fails_loudly(tmp_path: Path) -> None:
    """A truly empty input has no recoverable schema and must not exit 0."""
    tsv = tmp_path / "Broken_S99_blast.final.tsv"
    tsv.write_text("")
    output = tmp_path / "Broken_S99_blast_labkey.csv"

    result = subprocess.run(  # noqa: S603
        [
            sys.executable,
            str(SCRIPT),
            "--blast-csv",
            str(tsv),
            "--output",
            str(output),
            "--meta",
            "Broken_S99",
            "--experiment-id",
            "32495",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert not output.exists() or output.stat().st_size == 0


def test_concat_of_hitless_and_populated_succeeds(tmp_path: Path) -> None:
    """Reproduces LABKEY_CONCAT_ALL_SAMPLE_BLAST_RESULTS over a mixed batch."""
    pl = pytest.importorskip("polars")

    blast_results = tmp_path / "blast_results"
    blast_results.mkdir()
    for meta, rows in (("Water_S11", []), ("Sample_S34", [DATA_ROW])):
        output = run_prepare(tmp_path, write_tsv(tmp_path / f"{meta}.tsv", rows), meta)
        output.rename(blast_results / output.name)

    files = sorted(blast_results.glob("*.csv"))
    concatenated = tmp_path / "32495_blast_concatenated.csv"
    pl.concat([pl.scan_csv(f, infer_schema_length=0) for f in files]).sink_csv(
        concatenated,
    )

    rows = list(csv.reader(concatenated.open()))
    assert rows[0] == EXPECTED_HEADER
    assert rows[1] == ["32495", *DATA_ROW]
    assert len(rows) == 2  # noqa: PLR2004
