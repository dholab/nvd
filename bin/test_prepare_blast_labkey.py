"""Tests for the final-BLAST to LabKey projection boundary."""

import csv
import sys
from pathlib import Path

import pytest
from prepare_blast_labkey import main


def run_projection(
    monkeypatch: pytest.MonkeyPatch,
    blast_tsv: Path,
    output_csv: Path,
    *,
    experiment_id: str,
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "prepare_blast_labkey.py",
            "--blast-csv",
            str(blast_tsv),
            "--output",
            str(output_csv),
            "--meta",
            "sample-1",
            "--experiment-id",
            experiment_id,
        ],
    )
    main()


def test_query_class_column_survives_into_labkey_csv(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    blast_tsv = tmp_path / "sample.blast.final.tsv"
    blast_tsv.write_text(
        "sample\tqseqid\tquery_class\ttask\trank\n"
        "sample-1\tq1\tsingle_read\tmegablast\tspecies:x\n",
        encoding="utf-8",
    )
    output_csv = tmp_path / "sample.blast.labkey.csv"

    run_projection(
        monkeypatch,
        blast_tsv,
        output_csv,
        experiment_id="7",
    )

    with output_csv.open(newline="", encoding="utf-8") as handle:
        [row] = csv.DictReader(handle)
    assert row["experiment"] == "7"
    assert row["query_class"] == "single_read"
    assert row["sample_id"] == "sample-1"
    assert row["blast_task"] == "megablast"


def test_labkey_projection_omits_reference_placement_columns(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    blast_tsv = tmp_path / "sample.blast.final.tsv"
    output_csv = tmp_path / "sample.blast.labkey.csv"
    placement = {
        "saccver": "NC_000001.1",
        "qstart": "1",
        "qend": "100",
        "slen": "1000",
        "sstart": "900",
        "send": "801",
        "sstrand": "minus",
    }
    row = {
        "task": "megablast",
        "sample": "sample-1",
        "qseqid": "query-1",
        "sseqid": "ref|NC_000001.1|",
        **placement,
    }
    with blast_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    run_projection(
        monkeypatch,
        blast_tsv,
        output_csv,
        experiment_id="42",
    )

    with output_csv.open(newline="", encoding="utf-8") as handle:
        [projected] = csv.DictReader(handle)
    assert projected == {
        "experiment": "42",
        "blast_task": "megablast",
        "sample_id": "sample-1",
        "qseqid": "query-1",
        "sseqid": "ref|NC_000001.1|",
    }
    assert placement.keys().isdisjoint(projected)
