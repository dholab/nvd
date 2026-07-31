"""Tests for preserving taxonomy absence at the LabKey upload boundary."""

import sys

import labkey.api_wrapper
import polars as pl
import pytest
from labkey_upload_blast_results import (
    apply_reference_cutoff,
    combo_already_uploaded,
    dataframe_to_records,
    insert_records,
    main,
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


def test_insert_records_propagates_api_failure() -> None:
    """A failed insert_rows call must not be swallowed by insert_records."""

    class _FailingQuery:
        def insert_rows(self, **kwargs):
            raise RuntimeError("insert boom")

    with pytest.raises(RuntimeError, match="insert boom"):
        insert_records(_FailingQuery(), "lists", "hits", [{"experiment": 1}])


class _FakeQueryAPI:
    """Stand-in for the LabKey query API: configurable presence check + insert."""

    def __init__(self, rows, insert_error=None):
        self._rows = rows
        self._insert_error = insert_error
        self.inserted_rows = None

    def select_rows(self, **kwargs):
        return {"rows": self._rows}

    def insert_rows(self, **kwargs):
        if self._insert_error is not None:
            raise self._insert_error
        self.inserted_rows = kwargs.get("rows")
        return {"rows": [{"Key": 1}]}


def _fake_api_wrapper_class(fake_query: "_FakeQueryAPI"):
    """Build a fake APIWrapper class whose .query is the given fake, never the network."""

    class _FakeAPIWrapper:
        def __init__(self, *args, **kwargs):
            self.query = fake_query

    return _FakeAPIWrapper


def _cli_argv(csv_path, table_name: str = "hits") -> list[str]:
    return [
        "labkey_upload_blast_results.py",
        "--experiment-id",
        "7",
        "--sample-id",
        "s1",
        "--query-class",
        "single_read",
        "--csv",
        str(csv_path),
        "--labkey-server",
        "https://example.org",
        "--labkey-project-name",
        "proj",
        "--labkey-api-key",
        "key",
        "--labkey-schema",
        "lists",
        "--table-name",
        table_name,
    ]


def test_main_exits_nonzero_when_insert_fails(tmp_path, monkeypatch, capsys) -> None:
    """A failed atomic insert must hard-fail the run, not report success."""
    csv_path = tmp_path / "batch.csv"
    csv_path.write_text(
        "experiment,sample_id,qseqid,sseqid,bitscore,staxids\n7,s1,q1,ref-1,100.0,9606\n",
    )

    fake_query = _FakeQueryAPI(rows=[], insert_error=RuntimeError("insert boom"))
    monkeypatch.setattr(
        labkey.api_wrapper, "APIWrapper", _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_argv(csv_path))

    with pytest.raises(SystemExit) as exc_info:
        main()

    assert exc_info.value.code == 1
    log_text = (tmp_path / "blast_labkey_upload.log").read_text()
    assert "insert boom" in log_text
    assert "BLAST UPLOAD COMPLETE" not in log_text
    assert "insert boom" in capsys.readouterr().err


def test_main_skip_present_combo_exits_zero(tmp_path, monkeypatch) -> None:
    """A combo already present in the destination list is a no-op success, not a failure."""
    csv_path = tmp_path / "batch.csv"
    csv_path.write_text("experiment,sample_id,qseqid\n7,s1,q1\n")

    fake_query = _FakeQueryAPI(rows=[{"Key": 1}])
    monkeypatch.setattr(
        labkey.api_wrapper, "APIWrapper", _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_argv(csv_path))

    main()  # must return normally (implicit exit 0), never raise SystemExit

    log_text = (tmp_path / "blast_labkey_upload.log").read_text()
    assert "SKIP: combo already uploaded" in log_text
    assert fake_query.inserted_rows is None


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
