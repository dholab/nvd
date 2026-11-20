"""Tests for the full-size-vs-reduced e2e parity check.

The reduced e2e samplesheet carries one SRA run per distinct code path instead
of the original five. These tests pin the comparison that proves the reduction
did not change what the pipeline detects: a recorded summary of the full run is
compared against whatever the current run produced.
"""

from __future__ import annotations

from tests.e2e_parity import find_parity_violations, summarize_detections


def test_summarize_detections_groups_by_sample() -> None:
    rows = [
        {"sample": "sample_a", "adjusted_taxid": "10258", "saccver": "NC_005336.1"},
        {"sample": "sample_a", "adjusted_taxid": "10258", "saccver": "NC_005336.1"},
        {"sample": "sample_b", "adjusted_taxid": "10244", "saccver": "NC_063383.1"},
    ]

    assert summarize_detections(rows) == {
        "sample_a": {"taxids": ["10258"], "references": ["NC_005336.1"]},
        "sample_b": {"taxids": ["10244"], "references": ["NC_063383.1"]},
    }


def test_summarize_detections_ignores_blank_values() -> None:
    """Rows for a sample with no hit carry empty taxid and accession fields."""
    rows = [
        {"sample": "sample_a", "adjusted_taxid": "", "saccver": ""},
        {"sample": "sample_a", "adjusted_taxid": "10258", "saccver": "NC_005336.1"},
    ]

    assert summarize_detections(rows) == {
        "sample_a": {"taxids": ["10258"], "references": ["NC_005336.1"]},
    }


def test_parity_holds_for_subset() -> None:
    """Detecting fewer things than the full run is fine; detecting more is not."""
    full = {"big_sample": {"taxids": ["10244", "10258"], "references": ["NC_005336.1"]}}
    small = {"small_sample": {"taxids": ["10258"], "references": ["NC_005336.1"]}}

    errors, warnings = find_parity_violations(
        small,
        full,
        {"small_sample": "big_sample"},
    )

    assert errors == []
    assert warnings == []


def test_parity_reports_extra_taxid() -> None:
    """A detection the full run never made means the reduced run changed behaviour."""
    full = {"big_sample": {"taxids": ["10258"], "references": ["NC_005336.1"]}}
    small = {"small_sample": {"taxids": ["10258", "99999"], "references": ["NC_005336.1"]}}

    errors, _ = find_parity_violations(small, full, {"small_sample": "big_sample"})

    assert len(errors) == 1
    assert "small_sample" in errors[0]
    assert "99999" in errors[0]


def test_parity_extra_reference_is_warning_only() -> None:
    """Which genome inside a taxon wins is allowed to drift between BLAST DBs.

    The Mpox reference swap (NC_003310.1 -> NC_063383.1) is the real case: the
    taxon is unchanged, so this must not fail the suite.
    """
    full = {"big_sample": {"taxids": ["10244"], "references": ["NC_003310.1"]}}
    small = {"small_sample": {"taxids": ["10244"], "references": ["NC_063383.1"]}}

    errors, warnings = find_parity_violations(
        small,
        full,
        {"small_sample": "big_sample"},
    )

    assert errors == []
    assert len(warnings) == 1
    assert "NC_063383.1" in warnings[0]


def test_parity_reports_unmapped_sample() -> None:
    """A sample id with no recorded counterpart is a typo, not a pass."""
    full = {"big_sample": {"taxids": ["10258"], "references": []}}
    small = {"mystery_sample": {"taxids": ["10258"], "references": []}}

    errors, _ = find_parity_violations(small, full, {})

    assert len(errors) == 1
    assert "mystery_sample" in errors[0]


def test_parity_reports_a_mapped_sample_missing_from_the_run() -> None:
    """If a mapped sample produced no rows at all, the run is not comparable."""
    full = {"big_sample": {"taxids": ["10258"], "references": []}}

    errors, _ = find_parity_violations({}, full, {"small_sample": "big_sample"})

    assert len(errors) == 1
    assert "small_sample" in errors[0]


def test_parity_reports_a_mapping_to_an_unrecorded_sample() -> None:
    """Guard against a map that points at a sample the recording does not have."""
    small = {"small_sample": {"taxids": ["10258"], "references": []}}

    errors, _ = find_parity_violations(small, {}, {"small_sample": "absent_sample"})

    assert len(errors) == 1
    assert "absent_sample" in errors[0]
