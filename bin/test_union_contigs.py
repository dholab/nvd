"""Tests for long-read assembler contig union behavior."""

from __future__ import annotations

import csv
import json
import sqlite3
from typing import TYPE_CHECKING

import pytest
from union_contigs import (
    FinalizedUnionPaths,
    finalize_union,
    prepare_union,
    read_fasta,
    reverse_complement,
)

if TYPE_CHECKING:
    from pathlib import Path

EXPECTED_REPRESENTATIVE_COUNT = 3


def write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")
    return path


def write_manifest(path: Path, rows: list[tuple[str, Path]]) -> Path:
    lines = ["assembler\tsource_fasta"]
    lines.extend(f"{assembler}\t{fasta}" for assembler, fasta in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def rows(sqlite_path: Path, query: str) -> list[tuple[object, ...]]:
    with sqlite3.connect(sqlite_path) as connection:
        return list(connection.execute(query))


def prepare_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    first = write_fasta(
        tmp_path / "first.fasta",
        [
            ("alpha original", "AAAACCCCGGGG"),
            ("beta", "TTTTCCCC"),
        ],
    )
    second = write_fasta(
        tmp_path / "second.fasta",
        [
            ("alpha duplicate", "AAAACCCCGGGG"),
            ("beta reverse duplicate", reverse_complement("TTTTCCCC")),
            ("gamma", "GGGGAAAATTTT"),
        ],
    )
    manifest = write_manifest(
        tmp_path / "inputs.tsv",
        [("zeta", first), ("alpha", second)],
    )
    reps = tmp_path / "representatives.fasta"
    sqlite_path = tmp_path / "union.sqlite"
    meta = tmp_path / "union_meta.json"

    prepare_union(
        sample_id="Sample 1",
        inputs_manifest=manifest,
        representatives_fasta=reps,
        sqlite_path=sqlite_path,
        meta_json=meta,
    )
    return reps, sqlite_path, meta


def test_prepare_writes_representatives_and_sidecar_metadata(tmp_path: Path) -> None:
    reps, sqlite_path, meta = prepare_fixture(tmp_path)

    records = read_fasta(reps)
    metadata = json.loads(meta.read_text(encoding="utf-8"))

    assert len(records) == EXPECTED_REPRESENTATIVE_COUNT
    assert all(record.header.startswith("Sample_1__union_") for record in records)
    assert metadata == {
        "max_contig_length": 12,
        "representative_count": 3,
        "sample_id": "Sample 1",
    }
    assert rows(sqlite_path, "select count(*) from contig") == [(5,)]
    assert rows(sqlite_path, "select count(*) from duplicate") == [(2,)]
    assert rows(
        sqlite_path,
        "select name from sqlite_master where type = 'table' order by name",
    ) == [("containment",), ("contig",), ("duplicate",)]
    assert "sequence" not in {
        row[1] for row in rows(sqlite_path, "pragma table_info(contig)")
    }


def test_prepare_records_reverse_complement_duplicate_orientation(
    tmp_path: Path,
) -> None:
    _reps, sqlite_path, _meta = prepare_fixture(tmp_path)

    orientations = rows(
        sqlite_path,
        "select orientation, count(*) from duplicate group by orientation order by orientation",
    )

    assert orientations == [("reverse_complement", 1), ("same", 1)]


def test_finalize_exports_duplicate_provenance(tmp_path: Path) -> None:
    reps, sqlite_path, _meta = prepare_fixture(tmp_path)
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text("", encoding="utf-8")
    provenance = tmp_path / "provenance.tsv"
    summary = tmp_path / "summary.json"

    finalize_union(
        sample_id="Sample 1",
        representatives_fasta=reps,
        candidates_tsv=candidates,
        sqlite_path=sqlite_path,
        outputs=FinalizedUnionPaths(
            fasta=tmp_path / "final.fasta",
            provenance_tsv=provenance,
            summary_json=summary,
        ),
    )

    with provenance.open(newline="", encoding="utf-8") as handle:
        provenance_rows = list(csv.DictReader(handle, delimiter="\t"))
    duplicates = [row for row in provenance_rows if row["disposition"] == "duplicate"]
    assert {row["reason"] for row in duplicates} == {
        "canonical_sha256_exact_duplicate",
    }
    assert {row["orientation"] for row in duplicates} == {
        "same",
        "reverse_complement",
    }
    assert json.loads(summary.read_text(encoding="utf-8")) == {
        "emitted_contig_count": 3,
        "emitted_bases": 32,
        "exact_containment_count": 0,
        "exact_duplicate_count": 2,
        "representative_count": 3,
        "sample_id": "Sample 1",
        "schema_version": "nvd.long-read-union-summary/v1",
        "source_bases": 52,
        "source_contig_count": 5,
        "by_assembler": [
            {"producer": "alpha", "source_count": 3, "source_bases": 32},
            {"producer": "zeta", "source_count": 2, "source_bases": 20},
        ],
    }


def test_finalize_exports_empty_union_summary(tmp_path: Path) -> None:
    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.write_text("", encoding="utf-8")
    manifest = write_manifest(tmp_path / "inputs.tsv", [("assembler", empty_fasta)])
    reps = tmp_path / "representatives.fasta"
    sqlite_path = tmp_path / "union.sqlite"
    prepare_union(
        sample_id="empty-sample",
        inputs_manifest=manifest,
        representatives_fasta=reps,
        sqlite_path=sqlite_path,
        meta_json=tmp_path / "prepare.json",
    )
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text("", encoding="utf-8")
    provenance = tmp_path / "provenance.tsv"
    summary = tmp_path / "summary.json"

    finalize_union(
        sample_id="empty-sample",
        representatives_fasta=reps,
        candidates_tsv=candidates,
        sqlite_path=sqlite_path,
        outputs=FinalizedUnionPaths(
            fasta=tmp_path / "final.fasta",
            provenance_tsv=provenance,
            summary_json=summary,
        ),
    )

    with provenance.open(newline="", encoding="utf-8") as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []
    assert json.loads(summary.read_text(encoding="utf-8")) == {
        "emitted_contig_count": 0,
        "emitted_bases": 0,
        "exact_containment_count": 0,
        "exact_duplicate_count": 0,
        "representative_count": 0,
        "sample_id": "empty-sample",
        "schema_version": "nvd.long-read-union-summary/v1",
        "source_bases": 0,
        "source_contig_count": 0,
        "by_assembler": [],
    }


def test_finalize_drops_only_python_verified_exact_containment(tmp_path: Path) -> None:
    container_sequence = "ATGCCCAAATTTGG"
    forward_contained = "GCCCAAA"
    reverse_contained = reverse_complement("CCCAAAT")
    near_miss = "GCCCAAT"
    container = write_fasta(
        tmp_path / "container.fasta",
        [("container", container_sequence)],
    )
    fragments = write_fasta(
        tmp_path / "fragments.fasta",
        [
            ("contained", forward_contained),
            ("rc contained", reverse_contained),
            ("near miss", near_miss),
        ],
    )
    manifest = write_manifest(
        tmp_path / "inputs.tsv",
        [("a", container), ("b", fragments)],
    )
    reps = tmp_path / "representatives.fasta"
    sqlite_path = tmp_path / "union.sqlite"
    meta = tmp_path / "union_meta.json"
    prepare_union(
        sample_id="S1",
        inputs_manifest=manifest,
        representatives_fasta=reps,
        sqlite_path=sqlite_path,
        meta_json=meta,
    )
    rep_records = read_fasta(reps)
    by_sequence = {record.sequence: record.header for record in rep_records}
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text(
        "\n".join(
            [
                candidate_row(
                    by_sequence[forward_contained],
                    by_sequence[container_sequence],
                ),
                candidate_row(
                    by_sequence[reverse_contained],
                    by_sequence[container_sequence],
                ),
                candidate_row(by_sequence[near_miss], by_sequence[container_sequence]),
            ],
        )
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "final.fasta"
    provenance = tmp_path / "provenance.tsv"
    summary = tmp_path / "summary.json"

    finalize_union(
        sample_id="S1",
        representatives_fasta=reps,
        candidates_tsv=candidates,
        sqlite_path=sqlite_path,
        outputs=FinalizedUnionPaths(
            fasta=output,
            provenance_tsv=provenance,
            summary_json=summary,
        ),
    )

    emitted_sequences = {record.sequence for record in read_fasta(output)}
    assert emitted_sequences == {container_sequence, near_miss}
    assert rows(
        sqlite_path,
        "select orientation, start_1based, end_1based from containment order by orientation",
    ) == [
        ("reverse_complement", 4, 10),
        ("same", 3, 9),
    ]
    with provenance.open(newline="", encoding="utf-8") as handle:
        provenance_by_header = {
            row["original_header"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert set(provenance_by_header) == {
        "container",
        "contained",
        "rc contained",
        "near miss",
    }
    assert provenance_by_header["container"]["disposition"] == "emitted"
    assert provenance_by_header["near miss"]["disposition"] == "emitted"
    assert provenance_by_header["contained"] | {
        "contig_id": "ignored",
        "sha256": "ignored",
        "representative_id": "ignored",
        "container_id": "ignored",
    } == {
        "sample_id": "S1",
        "assembler": "b",
        "original_header": "contained",
        "contig_id": "ignored",
        "length": "7",
        "sha256": "ignored",
        "representative_id": "ignored",
        "disposition": "contained",
        "reason": "python_verified_exact_containment",
        "container_id": "ignored",
        "orientation": "same",
        "start_1based": "3",
        "end_1based": "9",
    }
    assert provenance_by_header["rc contained"]["orientation"] == ("reverse_complement")
    assert json.loads(summary.read_text(encoding="utf-8")) == {
        "emitted_contig_count": 2,
        "emitted_bases": 21,
        "exact_containment_count": 2,
        "exact_duplicate_count": 0,
        "representative_count": 4,
        "sample_id": "S1",
        "schema_version": "nvd.long-read-union-summary/v1",
        "source_bases": 35,
        "source_contig_count": 4,
        "by_assembler": [
            {"producer": "a", "source_count": 1, "source_bases": 14},
            {"producer": "b", "source_count": 3, "source_bases": 21},
        ],
    }


def test_finalize_warns_on_equal_length_duplicate_candidate(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    first = write_fasta(tmp_path / "first.fasta", [("one", "AAAACCCC")])
    manifest = write_manifest(tmp_path / "inputs.tsv", [("a", first)])
    reps = tmp_path / "representatives.fasta"
    sqlite_path = tmp_path / "union.sqlite"
    meta = tmp_path / "union_meta.json"
    prepare_union(
        sample_id="S1",
        inputs_manifest=manifest,
        representatives_fasta=reps,
        sqlite_path=sqlite_path,
        meta_json=meta,
    )
    rep_ids = ["S1__union_left", "S1__union_right"]
    write_fasta(reps, [(rep_ids[0], "AAAACCCC"), (rep_ids[1], "AAAACCCC")])
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text(
        candidate_row(rep_ids[0], rep_ids[1]) + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "final.fasta"
    finalize_union(
        sample_id="S1",
        representatives_fasta=reps,
        candidates_tsv=candidates,
        sqlite_path=sqlite_path,
        outputs=FinalizedUnionPaths(
            fasta=output,
            provenance_tsv=tmp_path / "provenance.tsv",
            summary_json=tmp_path / "summary.json",
        ),
    )

    captured = capsys.readouterr()
    assert "Equal-length exact duplicate candidate" in captured.err
    assert rows(sqlite_path, "select count(*) from containment") == [(0,)]
    assert {record.header for record in read_fasta(output)} == set(rep_ids)


def test_source_contig_ids_do_not_depend_on_source_fasta_path(tmp_path: Path) -> None:
    sequence = "AAAACCCCGGGG"
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"
    first_dir.mkdir()
    second_dir.mkdir()
    first = write_fasta(first_dir / "contigs.fasta", [("same header", sequence)])
    second = write_fasta(second_dir / "renamed.fasta", [("same header", sequence)])

    first_manifest = write_manifest(tmp_path / "first.tsv", [("asm", first)])
    second_manifest = write_manifest(tmp_path / "second.tsv", [("asm", second)])
    first_db = tmp_path / "first.sqlite"
    second_db = tmp_path / "second.sqlite"

    prepare_union(
        sample_id="S1",
        inputs_manifest=first_manifest,
        representatives_fasta=tmp_path / "first.reps.fasta",
        sqlite_path=first_db,
        meta_json=tmp_path / "first.json",
    )
    prepare_union(
        sample_id="S1",
        inputs_manifest=second_manifest,
        representatives_fasta=tmp_path / "second.reps.fasta",
        sqlite_path=second_db,
        meta_json=tmp_path / "second.json",
    )

    assert rows(first_db, "select contig_id from contig") == rows(
        second_db,
        "select contig_id from contig",
    )


def test_source_contig_ids_do_not_depend_on_manifest_order(tmp_path: Path) -> None:
    first = write_fasta(tmp_path / "first.fasta", [("one", "AAAACCCC")])
    second = write_fasta(tmp_path / "second.fasta", [("two", "GGGGTTTT")])
    forward_manifest = write_manifest(
        tmp_path / "forward.tsv",
        [("b", second), ("a", first)],
    )
    reversed_manifest = write_manifest(
        tmp_path / "reversed.tsv",
        [("a", first), ("b", second)],
    )
    forward_db = tmp_path / "forward.sqlite"
    reversed_db = tmp_path / "reversed.sqlite"

    prepare_union(
        sample_id="S1",
        inputs_manifest=forward_manifest,
        representatives_fasta=tmp_path / "forward.fasta",
        sqlite_path=forward_db,
        meta_json=tmp_path / "forward.json",
    )
    prepare_union(
        sample_id="S1",
        inputs_manifest=reversed_manifest,
        representatives_fasta=tmp_path / "reversed.fasta",
        sqlite_path=reversed_db,
        meta_json=tmp_path / "reversed.json",
    )

    query = "select contig_id, representative_id from contig order by contig_id"
    assert rows(forward_db, query) == rows(reversed_db, query)


def test_prepare_and_finalize_handle_all_empty_inputs(tmp_path: Path) -> None:
    empty_one = write_fasta(tmp_path / "empty-one.fasta", [])
    empty_two = write_fasta(tmp_path / "empty-two.fasta", [])
    manifest = write_manifest(
        tmp_path / "inputs.tsv",
        [("a", empty_one), ("b", empty_two)],
    )
    reps = tmp_path / "representatives.fasta"
    sqlite_path = tmp_path / "union.sqlite"
    meta = tmp_path / "union_meta.json"

    prepare_union(
        sample_id="S1",
        inputs_manifest=manifest,
        representatives_fasta=reps,
        sqlite_path=sqlite_path,
        meta_json=meta,
    )
    candidates = tmp_path / "candidates.tsv"
    candidates.write_text("", encoding="utf-8")
    output = tmp_path / "final.fasta"
    provenance = tmp_path / "provenance.tsv"
    summary = tmp_path / "summary.json"
    finalize_union(
        sample_id="S1",
        representatives_fasta=reps,
        candidates_tsv=candidates,
        sqlite_path=sqlite_path,
        outputs=FinalizedUnionPaths(
            fasta=output,
            provenance_tsv=provenance,
            summary_json=summary,
        ),
    )

    assert reps.read_text(encoding="utf-8") == ""
    assert output.read_text(encoding="utf-8") == ""
    assert provenance.read_text(encoding="utf-8").count("\n") == 1
    assert json.loads(summary.read_text(encoding="utf-8")) == {
        "emitted_contig_count": 0,
        "emitted_bases": 0,
        "exact_containment_count": 0,
        "exact_duplicate_count": 0,
        "representative_count": 0,
        "sample_id": "S1",
        "schema_version": "nvd.long-read-union-summary/v1",
        "source_bases": 0,
        "source_contig_count": 0,
        "by_assembler": [],
    }
    assert json.loads(meta.read_text(encoding="utf-8")) == {
        "max_contig_length": 0,
        "representative_count": 0,
        "sample_id": "S1",
    }
    assert rows(sqlite_path, "select count(*) from contig") == [(0,)]


def test_prepare_ignores_empty_inputs_when_other_assemblers_emit_contigs(
    tmp_path: Path,
) -> None:
    empty = write_fasta(tmp_path / "empty.fasta", [])
    non_empty = write_fasta(tmp_path / "non-empty.fasta", [("one", "AAAACCCC")])
    manifest = write_manifest(tmp_path / "inputs.tsv", [("a", empty), ("b", non_empty)])
    reps = tmp_path / "representatives.fasta"

    prepare_union(
        sample_id="S1",
        inputs_manifest=manifest,
        representatives_fasta=reps,
        sqlite_path=tmp_path / "union.sqlite",
        meta_json=tmp_path / "union_meta.json",
    )

    assert [record.sequence for record in read_fasta(reps)] == ["AAAACCCC"]


def test_read_fasta_rejects_sequence_before_header(tmp_path: Path) -> None:
    fasta = tmp_path / "bad.fasta"
    fasta.write_text("ACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="before first header"):
        read_fasta(fasta)


def test_read_fasta_rejects_unsupported_symbols(tmp_path: Path) -> None:
    fasta = write_fasta(tmp_path / "bad.fasta", [("contig", "ACGTZ")])

    with pytest.raises(ValueError, match="Unsupported FASTA symbols"):
        read_fasta(fasta)


def test_finalize_rejects_duplicate_representative_headers(tmp_path: Path) -> None:
    reps = write_fasta(tmp_path / "reps.fasta", [("dup", "AAAA"), ("dup", "CCCC")])
    db = tmp_path / "union.sqlite"
    manifest = write_manifest(
        tmp_path / "inputs.tsv",
        [("asm", write_fasta(tmp_path / "source.fasta", [("one", "AAAA")]))],
    )
    prepare_union(
        sample_id="S1",
        inputs_manifest=manifest,
        representatives_fasta=tmp_path / "prepared.fasta",
        sqlite_path=db,
        meta_json=tmp_path / "meta.json",
    )

    with pytest.raises(ValueError, match="Duplicate representative FASTA header"):
        finalize_union(
            sample_id="S1",
            representatives_fasta=reps,
            candidates_tsv=tmp_path / "empty.tsv",
            sqlite_path=db,
            outputs=FinalizedUnionPaths(
                fasta=tmp_path / "final.fasta",
                provenance_tsv=tmp_path / "provenance.tsv",
                summary_json=tmp_path / "summary.json",
            ),
        )


def candidate_row(query: str, target: str) -> str:
    return f"{query}\t{target}\t100.000\t8\t0\t0\t1\t8\t1\t8\t8\t16\t1.000\t0.500\t0\t1"
