"""Tests for preserving taxonomy absence at the LabKey upload boundary."""

import polars as pl
from labkey_upload_blast_results import (
    apply_reference_cutoff,
    combo_already_uploaded,
    dataframe_to_records,
    validate_dataframe,
)


class _FakeQuery:
    def __init__(self, rows):
        self._rows = rows
        self.calls = []

    def select_rows(self, **kwargs):
        self.calls.append(kwargs)
        return {"rows": self._rows}


def test_combo_present_is_detected() -> None:
    q = _FakeQuery(rows=[{"Key": 1}])
    assert combo_already_uploaded(q, "lists", "hits", 7, "s1", "single_read") is True


def test_combo_absent_is_detected() -> None:
    q = _FakeQuery(rows=[])
    assert combo_already_uploaded(q, "lists", "hits", 7, "s1", "single_read") is False
    # filters must include all three identity columns
    (call,) = q.calls
    filtered = {f.column_name if hasattr(f, "column_name") else str(f) for f in call["filter_array"]}
    assert any("experiment" in str(f) for f in filtered)
    assert any("sample_id" in str(f) for f in filtered)
    assert any("query_class" in str(f) for f in filtered)


def test_unavailable_staxid_remains_absent_in_labkey_record() -> None:
    logs: list[str] = []
    frame = pl.DataFrame(
        {
            "experiment": [1],
            "sample_id": ["sample-1"],
            "qseqid": ["query-1"],
            "staxids": [None],
            "adjusted_taxid": [9606],
        },
        schema_overrides={"staxids": pl.Int64},
    )

    [record] = dataframe_to_records(validate_dataframe(frame, "blast.csv", logs))

    assert record["staxids"] is None
    assert record["staxids"] != 0


def test_reference_cutoff_bounds_tied_references_per_contig() -> None:
    """A contig whose references tie at the top bitscore keeps only top_k of them.

    core_nt redundancy makes many references share one bitscore; the production list
    should keep an ordinal top_k per contig rather than every tied reference.
    """
    frame = pl.DataFrame(
        {
            "sample_id": ["s1"] * 8,
            "qseqid": ["q1"] * 8,
            "sseqid": [f"ref-{index:02d}" for index in range(8)],
            "bitscore": [100.0] * 8,
        },
    )

    trimmed = apply_reference_cutoff(frame, top_k=5)

    kept = set(trimmed["sseqid"].to_list())
    assert len(kept) == 5
    assert kept == {f"ref-{index:02d}" for index in range(5)}  # deterministic: lowest sseqid


def test_reference_cutoff_is_per_group_and_keeps_all_taxid_rows() -> None:
    """The cutoff applies within each (sample_id, qseqid); a multi-taxid reference is one hit."""
    frame = pl.DataFrame(
        {
            "sample_id": ["s1", "s1", "s1", "s1", "s1", "s1"],
            "qseqid": ["q1", "q1", "q1", "q1", "q2", "q2"],
            "sseqid": ["ref-a", "ref-a", "ref-b", "ref-c", "ref-x", "ref-y"],
            "bitscore": [100.0, 100.0, 99.0, 98.0, 50.0, 49.0],
            "staxids": [10, 11, 12, 13, 20, 21],
        },
    )

    trimmed = apply_reference_cutoff(frame, top_k=2)

    q1 = trimmed.filter(pl.col("qseqid") == "q1")
    q2 = trimmed.filter(pl.col("qseqid") == "q2")
    # q1: top 2 references = ref-a (two taxid rows) + ref-b (one) -> 3 rows, 2 references.
    assert set(q1["sseqid"].to_list()) == {"ref-a", "ref-b"}
    assert q1.height == 3
    # q2: only two references, both kept, unaffected by q1's cutoff.
    assert set(q2["sseqid"].to_list()) == {"ref-x", "ref-y"}
