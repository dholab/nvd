"""Tests for typed BLAST task-breakdown report inputs."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from py_nvd.multiqc_blast import (
    InvalidMegablastQueryPartition,
    ParsedMegablastQueryPartition,
    blast_task_rows,
    collect_megablast_query_partitions,
)
from py_nvd.multiqc_packages import MegablastQueryPartitionReceipt, ReportPackage

if TYPE_CHECKING:
    from pathlib import Path


def test_megablast_partition_package_becomes_task_breakdown_row(
    tmp_path: Path,
) -> None:
    summary = tmp_path / "summary.payload"
    summary.write_text(
        "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\tproducer_note\n"
        "sample_A\tshort_assembly_contig\t4\t2\t2\tignored\n",
        encoding="utf-8",
    )
    receipt = MegablastQueryPartitionReceipt(
        sample_id="sample_A",
        query_class="short_assembly_contig",
    )
    package = ReportPackage(root=tmp_path, receipt=receipt)

    [partition] = collect_megablast_query_partitions((package,))
    assert isinstance(partition, ParsedMegablastQueryPartition)

    assert blast_task_rows((partition,), sample_order={"sample_A": 0})[0].model_dump(
        exclude={"query_class_code"},
    ) == {
        "sample_id": "sample_A",
        "query_class": "Short assembly contig",
        "queries_searched": 4,
        "matched_by_megablast": 2,
        "forwarded_to_blastn": 2,
        "forwarded_to_blastn_percentage": 50.0,
        "problem": None,
        "validation_details": None,
    }


def test_zero_denominator_and_malformed_payload_remain_explicit(
    tmp_path: Path,
) -> None:
    valid_root = tmp_path / "valid"
    valid_root.mkdir()
    (valid_root / "summary.payload").write_text(
        "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\n"
        "sample_A\tsingle_read\t0\t0\t0\n",
        encoding="utf-8",
    )
    valid_receipt = MegablastQueryPartitionReceipt(
        sample_id="sample_A",
        query_class="single_read",
    )

    invalid_root = tmp_path / "invalid"
    invalid_root.mkdir()
    (invalid_root / "summary.payload").write_text(
        "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\n"
        "sample_B\tlong_assembly_contig\t3\t1\t2\n"
        "sample_B\tlong_assembly_contig\t3\t1\t2\n",
        encoding="utf-8",
    )
    invalid_receipt = MegablastQueryPartitionReceipt(
        sample_id="sample_B",
        query_class="long_assembly_contig",
    )

    partitions = collect_megablast_query_partitions(
        (
            ReportPackage(root=valid_root, receipt=valid_receipt),
            ReportPackage(root=invalid_root, receipt=invalid_receipt),
        ),
    )
    assert isinstance(partitions[1], InvalidMegablastQueryPartition)

    rows = blast_task_rows(
        partitions,
        sample_order={"sample_A": 0, "sample_B": 1},
    )
    assert rows[0].forwarded_to_blastn_percentage is None
    assert rows[1].query_class_code == "long_assembly_contig"
    assert rows[1].problem == "BLAST task breakdown could not be read"
    assert "exactly one row" in (rows[1].validation_details or "")


@pytest.mark.parametrize(
    "summary_text",
    [
        pytest.param(
            "sample_id\tquery_class\tqueries_in\tmegablast_accounted\n"
            "sample_A\tshort_assembly_contig\t3\t1\n",
            id="missing-required-column",
        ),
        pytest.param(
            "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\n"
            "sample_B\tshort_assembly_contig\t3\t1\t2\n",
            id="sample-mismatch",
        ),
        pytest.param(
            "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\n"
            "sample_A\tlong_assembly_contig\t3\t1\t2\n",
            id="query-class-mismatch",
        ),
        pytest.param(
            "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\n"
            "sample_A\tunknown_class\t3\t1\t2\n",
            id="unknown-query-class",
        ),
        pytest.param(
            "sample_id\tquery_class\tqueries_in\tmegablast_accounted\tblastn_candidates\n"
            "sample_A\tshort_assembly_contig\t3\t-1\t2\n",
            id="negative-count",
        ),
    ],
)
def test_malformed_partition_payload_is_localized_to_its_receipt_coordinate(
    tmp_path: Path,
    summary_text: str,
) -> None:
    (tmp_path / "summary.payload").write_text(summary_text, encoding="utf-8")
    receipt = MegablastQueryPartitionReceipt(
        sample_id="sample_A",
        query_class="short_assembly_contig",
    )

    [partition] = collect_megablast_query_partitions(
        (ReportPackage(root=tmp_path, receipt=receipt),),
    )
    assert isinstance(partition, InvalidMegablastQueryPartition)

    [row] = blast_task_rows((partition,), sample_order={"sample_A": 0})
    assert row.query_class_code == "short_assembly_contig"
    assert row.problem == "BLAST task breakdown could not be read"
    assert row.validation_details
