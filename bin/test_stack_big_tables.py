"""Tests for stacking per-sample Big Table TSVs."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

import pytest
from stack_big_tables import main

if TYPE_CHECKING:
    from pathlib import Path


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stack_big_tables_preserves_schema_and_sorts_inputs(tmp_path: Path) -> None:
    input_dir = tmp_path / "sample_big_tables"
    input_dir.mkdir()
    output = tmp_path / "query_big_table.tsv"
    columns = ["sample_id", "qseqid", "support_tier"]
    write_tsv(
        input_dir / "sample-b.query_big_table.tsv",
        [{"sample_id": "sample-b", "qseqid": "q2", "support_tier": "moderate"}],
        columns,
    )
    write_tsv(
        input_dir / "sample-a.query_big_table.tsv",
        [{"sample_id": "sample-a", "qseqid": "q1", "support_tier": "strong"}],
        columns,
    )

    main(["--input-dir", str(input_dir), "--output", str(output)])

    rows = read_tsv(output)
    assert list(rows[0]) == columns
    assert [row["sample_id"] for row in rows] == ["sample-a", "sample-b"]


def test_stack_big_tables_rejects_column_order_mismatch(tmp_path: Path) -> None:
    input_dir = tmp_path / "sample_big_tables"
    input_dir.mkdir()
    output = tmp_path / "taxon_big_table.tsv"
    write_tsv(
        input_dir / "sample-a.taxon_big_table.tsv",
        [{"sample_id": "sample-a", "taxid": "111", "support_tier": "strong"}],
        ["sample_id", "taxid", "support_tier"],
    )
    write_tsv(
        input_dir / "sample-b.taxon_big_table.tsv",
        [{"sample_id": "sample-b", "support_tier": "moderate", "taxid": "222"}],
        ["sample_id", "support_tier", "taxid"],
    )

    with pytest.raises(ValueError, match="schema does not match"):
        main(["--input-dir", str(input_dir), "--output", str(output)])


def test_stack_big_tables_writes_empty_file_for_empty_input_dir(tmp_path: Path) -> None:
    input_dir = tmp_path / "sample_big_tables"
    input_dir.mkdir()
    output = tmp_path / "query_big_table.tsv"

    main(["--input-dir", str(input_dir), "--output", str(output)])

    assert output.read_text(encoding="utf-8") == ""
