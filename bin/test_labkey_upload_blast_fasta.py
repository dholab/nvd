"""Tests for deduping the contig FASTA LabKey upload by (experiment, sample_id)."""

import sys

import labkey.api_wrapper
import pytest
from labkey_upload_blast_fasta import insert_records, main, sample_already_uploaded


class _FakeQuery:
    def __init__(self, rows):
        self._rows = rows
        self.calls = []

    def select_rows(self, **kwargs):
        self.calls.append(kwargs)
        return {"rows": self._rows}


def test_sample_present_is_detected() -> None:
    assert sample_already_uploaded(_FakeQuery([{"Key": 1}]), "lists", "fasta", 7, "s1") is True


def test_sample_absent_is_detected() -> None:
    q = _FakeQuery([])
    assert sample_already_uploaded(q, "lists", "fasta", 7, "s1") is False
    (call,) = q.calls
    assert any("experiment" in str(f) for f in call["filter_array"])
    assert any("sample_id" in str(f) for f in call["filter_array"])


def test_insert_records_propagates_api_failure() -> None:
    """A failed insert_rows call must not be swallowed by insert_records."""

    class _FailingQuery:
        def insert_rows(self, **kwargs):
            raise RuntimeError("insert boom")

    with pytest.raises(RuntimeError, match="insert boom"):
        insert_records(_FailingQuery(), "lists", "fasta", [{"experiment": 1}])


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


def _cli_argv(table_name: str = "fasta") -> list[str]:
    return [
        "labkey_upload_blast_fasta.py",
        "--experiment-id",
        "7",
        "--sample-id",
        "s1",
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
    csv_path = tmp_path / "s1_fasta.csv"
    csv_path.write_text(
        "experiment,sample_id,contig_id,contig_sequence,notes,nextflow_run_id\n"
        "7,s1,contig_1,ACGT,,run1\n",
    )

    fake_query = _FakeQueryAPI(rows=[], insert_error=RuntimeError("insert boom"))
    monkeypatch.setattr(
        labkey.api_wrapper, "APIWrapper", _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_argv())

    with pytest.raises(SystemExit) as exc_info:
        main()

    assert exc_info.value.code == 1
    log_text = (tmp_path / "fasta_labkey_upload.log").read_text()
    assert "insert boom" in log_text
    assert "FASTA UPLOAD COMPLETE" not in log_text
    assert "insert boom" in capsys.readouterr().err


def test_main_skip_present_sample_exits_zero(tmp_path, monkeypatch) -> None:
    """A sample already present in the destination list is a no-op success, not a failure."""
    csv_path = tmp_path / "s1_fasta.csv"
    csv_path.write_text(
        "experiment,sample_id,contig_id,contig_sequence,notes,nextflow_run_id\n"
        "7,s1,contig_1,ACGT,,run1\n",
    )

    fake_query = _FakeQueryAPI(rows=[{"Key": 1}])
    monkeypatch.setattr(
        labkey.api_wrapper, "APIWrapper", _fake_api_wrapper_class(fake_query),
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", _cli_argv())

    main()  # must return normally (implicit exit 0), never raise SystemExit

    log_text = (tmp_path / "fasta_labkey_upload.log").read_text()
    assert "SKIP: sample already uploaded" in log_text
    assert fake_query.inserted_rows is None
