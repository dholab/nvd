"""Tests for typed report-package construction."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from py_nvd.multiqc_packages import (
    DeaconReceipt,
    Domain,
    EligibilityReceipt,
    LongReadEligibilityReceipt,
    ProfileReceipt,
    UnionReceipt,
    parse_receipt,
    write_package,
)

if TYPE_CHECKING:
    from pathlib import Path


@pytest.mark.parametrize(
    ("receipt", "payload_count"),
    [
        (
            DeaconReceipt(
                domain=Domain.DEPLETION,
                sample_id="sample_A",
                query_class="reads",
            ),
            1,
        ),
        (
            ProfileReceipt(
                domain=Domain.FASTX,
                sample_id="sample_A",
                stage="post_depletion",
                quality_histogram="quality_histogram.payload",
            ),
            3,
        ),
        (
            EligibilityReceipt(
                sample_id="sample_A",
                producer="spades",
                decision="run",
                sequence_count=120,
                threshold=100,
                reason="short_read_minimum_sequence_count",
            ),
            0,
        ),
        (LongReadEligibilityReceipt(sample_id="sample_A"), 1),
        (UnionReceipt(sample_id="sample_A"), 1),
    ],
)
def test_receipts_round_trip_through_materialized_packages(
    tmp_path: Path,
    receipt: DeaconReceipt
    | ProfileReceipt
    | EligibilityReceipt
    | LongReadEligibilityReceipt
    | UnionReceipt,
    payload_count: int,
) -> None:
    payloads = tuple(tmp_path / f"source_{index}" for index in range(payload_count))
    for payload in payloads:
        payload.write_text("payload", encoding="utf-8")

    package = write_package(receipt, payloads, output_root=tmp_path)

    assert package.name == receipt.package_name
    assert parse_receipt(package / "receipt.json") == receipt
    assert {
        path.name for path in package.iterdir() if path.name != "receipt.json"
    } == set(receipt.payload_names)
