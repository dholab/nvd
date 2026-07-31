"""Tests for deterministic, header-aware TSV concatenation."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

import pytest
from concatenate_tsvs import main

if TYPE_CHECKING:
    from pathlib import Path


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_concatenate_tsvs_preserves_header_and_sorts_inputs(tmp_path: Path) -> None:
    first_stage = tmp_path / "batch000001"
    second_stage = tmp_path / "batch000002"
    first_stage.mkdir()
    second_stage.mkdir()
    output = tmp_path / "query_big_table.tsv"
    columns = ["sample_id", "qseqid", "support_tier"]
    sample_b = first_stage / "sample-b.query_big_table.tsv"
    sample_a = second_stage / "sample-a.query_big_table.tsv"
    write_tsv(
        sample_b,
        [{"sample_id": "sample-b", "qseqid": "q2", "support_tier": "moderate"}],
        columns,
    )
    write_tsv(
        sample_a,
        [{"sample_id": "sample-a", "qseqid": "q1", "support_tier": "strong"}],
        columns,
    )

    main(["--input", str(sample_b), str(sample_a), "--output", str(output)])

    assert output.read_text(encoding="utf-8") == (
        "sample_id\tqseqid\tsupport_tier\n"
        "sample-a\tq1\tstrong\n"
        "sample-b\tq2\tmoderate\n"
    )


@pytest.mark.parametrize(
    "second_contents",
    [
        "qseqid\tsample\nq2\ts2\n",
        "sample\ttaxid\ns2\t2\n",
    ],
)
def test_concatenate_tsvs_rejects_nonmatching_headers(
    tmp_path: Path,
    second_contents: str,
) -> None:
    first = tmp_path / "a.tsv"
    second = tmp_path / "b.tsv"
    output = tmp_path / "combined.tsv"
    first.write_text("sample\tqseqid\ns1\tq1\n", encoding="utf-8")
    second.write_text(second_contents, encoding="utf-8")

    with pytest.raises(ValueError, match="matching input TSV schemas") as error:
        main(["--input", str(first), str(second), "--output", str(output)])

    assert str(second) in str(error.value)


def test_concatenate_tsvs_writes_empty_file_for_empty_input_dir(tmp_path: Path) -> None:
    input_dir = tmp_path / "sample_big_tables"
    input_dir.mkdir()
    output = tmp_path / "query_big_table.tsv"

    main(["--input-dir", str(input_dir), "--output", str(output)])

    assert output.read_text(encoding="utf-8") == ""


def test_concatenate_tsvs_preserves_one_header_for_header_only_inputs(
    tmp_path: Path,
) -> None:
    first = tmp_path / "a.tsv"
    second = tmp_path / "b.tsv"
    output = tmp_path / "combined.tsv"
    write_tsv(first, [], ["sample_id", "qseqid"])
    write_tsv(second, [], ["sample_id", "qseqid"])

    main(["--input", str(second), str(first), "--output", str(output)])

    assert output.read_text(encoding="utf-8") == "sample_id\tqseqid\n"


@pytest.mark.parametrize("empty_name", ["a.tsv", "b.tsv"])
def test_concatenate_tsvs_rejects_empty_input_file(
    tmp_path: Path,
    empty_name: str,
) -> None:
    valid_name = "b.tsv" if empty_name == "a.tsv" else "a.tsv"
    valid = tmp_path / valid_name
    empty = tmp_path / empty_name
    output = tmp_path / "combined.tsv"
    valid.write_text("sample\tqseqid\ns1\tq1\n", encoding="utf-8")
    empty.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="matching input TSV schemas") as error:
        main(["--input", str(valid), str(empty), "--output", str(output)])

    assert str(empty) in str(error.value)


def test_concatenate_tsvs_rejects_missing_input_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing.tsv"
    output = tmp_path / "combined.tsv"

    with pytest.raises(FileNotFoundError) as error:
        main(["--input", str(missing), "--output", str(output)])

    assert str(missing) in str(error.value)


def test_concatenate_tsvs_rejects_missing_input_directory(tmp_path: Path) -> None:
    missing = tmp_path / "missing"
    output = tmp_path / "combined.tsv"

    with pytest.raises(
        FileNotFoundError,
        match="input directory does not exist",
    ) as error:
        main(["--input-dir", str(missing), "--output", str(output)])

    assert str(missing) in str(error.value)


def test_concatenate_tsvs_rejects_invalid_numeric_values(
    tmp_path: Path,
) -> None:
    invalid = tmp_path / "invalid.tsv"
    output = tmp_path / "combined.tsv"
    invalid.write_text("sample\tbitscore\ns1\tnot-a-number\n", encoding="utf-8")

    with pytest.raises(ValueError, match="could not concatenate input TSVs") as error:
        main(["--input", str(invalid), "--output", str(output)])

    assert str(invalid) in str(error.value)
    assert "not-a-number" in str(error.value)
    assert "bitscore" in str(error.value)
