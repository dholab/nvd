"""Tests for typed prepared-query report inputs."""

from __future__ import annotations

import csv
from typing import TYPE_CHECKING

import pytest

from py_nvd.multiqc_packages import (
    PreparedQueryBatchesReceipt,
    ReportPackage,
    write_package,
)
from py_nvd.multiqc_query_preparation import (
    InvalidPreparedQuerySummary,
    ParsedPreparedQuerySummary,
    collect_prepared_query_summaries,
    prepared_query_rows,
)

if TYPE_CHECKING:
    from pathlib import Path


def summary_rows() -> list[dict[str, object]]:
    return [
        {
            "sample_id": "sample_A",
            "platform": "illumina",
            "query_class": "short_assembly_contig",
            "query_source": "contig",
            "n_query_sequences": 2,
            "query_fasta_present": True,
            "query_lookup_present": True,
            "producer_note": "ignored additive metadata",
        },
        {
            "sample_id": "sample_A",
            "platform": "illumina",
            "query_class": "long_assembly_contig",
            "query_source": "contig",
            "n_query_sequences": 0,
            "query_fasta_present": False,
            "query_lookup_present": False,
            "producer_note": "ignored additive metadata",
        },
        {
            "sample_id": "sample_A",
            "platform": "illumina",
            "query_class": "overlap_merged_pair",
            "query_source": "read_query",
            "n_query_sequences": 0,
            "query_fasta_present": False,
            "query_lookup_present": False,
            "producer_note": "ignored additive metadata",
        },
        {
            "sample_id": "sample_A",
            "platform": "illumina",
            "query_class": "single_read",
            "query_source": "read_query",
            "n_query_sequences": 0,
            "query_fasta_present": False,
            "query_lookup_present": False,
            "producer_note": "ignored additive metadata",
        },
    ]


def write_summary(path: Path, rows: list[dict[str, object]]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def packaged_summary(tmp_path: Path, rows: list[dict[str, object]]) -> ReportPackage:
    receipt = PreparedQueryBatchesReceipt(sample_id="sample_A")
    source = write_summary(tmp_path / "summary.tsv", rows)
    package = write_package(receipt, (source,), output_root=tmp_path)
    return ReportPackage(root=package, receipt=receipt)


def test_projection_ignores_additive_columns_and_uses_report_plan_availability(
    tmp_path: Path,
) -> None:
    package = packaged_summary(tmp_path, summary_rows())
    [summary] = collect_prepared_query_summaries((package,))

    assert isinstance(summary, ParsedPreparedQuerySummary)
    assert [batch.query_class for batch in summary.batches] == [
        "short_assembly_contig",
        "long_assembly_contig",
        "overlap_merged_pair",
        "single_read",
    ]

    standard_rows = prepared_query_rows(
        (summary,),
        sample_order={"sample_A": 0},
        sample_platforms={"sample_A": "illumina"},
        assembly_enabled=True,
        unassembled_read_querying_enabled=False,
    )
    assert [row.availability for row in standard_rows] == [
        "enabled",
        "disabled",
        "disabled",
    ]
    assert [row.query_sequences for row in standard_rows] == [2, None, None]
    assert [row.sequence_type for row in standard_rows] == [
        "Assembly contig",
        "Unassembled read",
        "Unassembled read",
    ]

    experimental_rows = prepared_query_rows(
        (summary,),
        sample_order={"sample_A": 0},
        sample_platforms={"sample_A": "illumina"},
        assembly_enabled=True,
        unassembled_read_querying_enabled=True,
    )
    assert [row.query_sequences for row in experimental_rows] == [2]


def test_enabled_empty_batch_requires_success_artifacts(tmp_path: Path) -> None:
    rows = summary_rows()
    rows[-1].update(query_fasta_present=True, query_lookup_present=True)
    package = packaged_summary(tmp_path, rows)
    [summary] = collect_prepared_query_summaries((package,))

    report_rows = prepared_query_rows(
        (summary,),
        sample_order={"sample_A": 0},
        sample_platforms={"sample_A": "illumina"},
        assembly_enabled=True,
        unassembled_read_querying_enabled=True,
    )

    assert [row.query_class_code for row in report_rows] == [
        "short_assembly_contig",
        "single_read",
    ]
    assert report_rows[-1].availability == "enabled"
    assert report_rows[-1].query_sequences == 0


def test_artifactless_zeroes_do_not_claim_enabled_observations(tmp_path: Path) -> None:
    rows = summary_rows()
    rows[0].update(
        n_query_sequences=0,
        query_fasta_present=False,
        query_lookup_present=False,
    )
    package = packaged_summary(tmp_path, rows)
    [summary] = collect_prepared_query_summaries((package,))

    standard_rows = prepared_query_rows(
        (summary,),
        sample_order={"sample_A": 0},
        sample_platforms={"sample_A": "illumina"},
        assembly_enabled=True,
        unassembled_read_querying_enabled=False,
    )
    experimental_rows = prepared_query_rows(
        (summary,),
        sample_order={"sample_A": 0},
        sample_platforms={"sample_A": "illumina"},
        assembly_enabled=True,
        unassembled_read_querying_enabled=True,
    )

    assert [row.query_class_code for row in standard_rows] == [
        "overlap_merged_pair",
        "single_read",
    ]
    assert all(row.availability == "disabled" for row in standard_rows)
    assert experimental_rows == ()


@pytest.mark.parametrize(
    "mutate",
    [
        lambda rows: rows.pop(),
        lambda rows: rows.append(rows[0].copy()),
        lambda rows: rows[0].update(query_class="unknown"),
        lambda rows: rows[0].update(query_source="read_query"),
        lambda rows: rows[0].update(query_lookup_present=False),
    ],
    ids=(
        "missing-class",
        "duplicate-class",
        "unknown-class",
        "wrong-source",
        "inconsistent-presence",
    ),
)
def test_semantic_summary_defects_become_one_invalid_sample_row(
    tmp_path: Path,
    mutate: object,
) -> None:
    rows = summary_rows()
    mutate(rows)  # type: ignore[operator]
    package = packaged_summary(tmp_path, rows)
    [summary] = collect_prepared_query_summaries((package,))

    assert isinstance(summary, InvalidPreparedQuerySummary)
    report_rows = prepared_query_rows(
        (summary,),
        sample_order={"sample_A": 0},
        sample_platforms={"sample_A": "illumina"},
        assembly_enabled=True,
        unassembled_read_querying_enabled=True,
    )
    assert len(report_rows) == 1
    assert report_rows[0].availability == "invalid"
    assert report_rows[0].query_sequences is None
    assert report_rows[0].problem == "Prepared query summary could not be read"
    assert report_rows[0].validation_details
