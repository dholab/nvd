"""Tests for the experiment-wide sequence-flow ledger."""

from __future__ import annotations

import csv
import json
from typing import TYPE_CHECKING

from build_sequence_flow import OUTPUT_FIELDS, build_rows, write_rows

if TYPE_CHECKING:
    from pathlib import Path

QUERY_SEQUENCE_COUNT = 2
SOURCE_CONTIG_COUNT = 7
EMITTED_CONTIG_COUNT = 4
EXPECTED_ROW_COUNT = 5


def write(path: Path, content: str) -> Path:
    path.write_text(content, encoding="utf-8")
    return path


def test_builds_deterministic_rows_from_stage_specific_evidence(tmp_path: Path) -> None:
    query_batches = write(
        tmp_path / "S1.blast_query_batches.tsv",
        "sample_id\tplatform\tquery_class\tquery_source\tn_query_sequences\tquery_fasta_present\tquery_lookup_present\tquery_fasta\tquery_lookup\n"
        "S1\tillumina\tsingle_read\tread_query\t2\ttrue\ttrue\tq.fa\tq.db\n"
        "S1\tillumina\tlong_assembly_contig\tcontig\t0\tfalse\tfalse\t\t\n",
    )
    partition = write(
        tmp_path / "S1.single_read.megablast_query_partition.tsv",
        "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\n"
        "S1\tsingle_read\t2\t1\t1\n",
    )
    union_summary = tmp_path / "S2.long_read_union.summary.json"
    union_summary.write_text(
        json.dumps(
            {
                "sample_id": "S2",
                "source_contig_count": 7,
                "emitted_contig_count": 4,
            },
        ),
        encoding="utf-8",
    )

    rows = build_rows([union_summary, partition, query_batches])

    assert [row.stage for row in rows] == [
        "blast_query_preparation",
        "blast_query_preparation",
        "megablast_query_partition",
        "megablast_query_partition",
        "long_read_contig_union",
    ]
    assert rows[1].sequences_out == QUERY_SEQUENCE_COUNT
    assert rows[2].output_class == "blastn_candidates"
    assert rows[3].output_class == "megablast_accounted"
    assert rows[4].sequences_in == SOURCE_CONTIG_COUNT
    assert rows[4].sequences_out == EMITTED_CONTIG_COUNT

    output = tmp_path / "sequence_flow.tsv"
    write_rows(output, rows)
    with output.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == OUTPUT_FIELDS
        assert len(list(reader)) == EXPECTED_ROW_COUNT
