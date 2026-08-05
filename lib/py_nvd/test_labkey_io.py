"""Tests for the shared LabKey list I/O helpers."""

from __future__ import annotations

import pytest

from py_nvd.labkey_io import insert_records, rows_present


class _FakeQuery:
    def __init__(
        self,
        rows: list[dict[str, object]] | None = None,
        *,
        raise_on_insert: bool = False,
    ) -> None:
        self._rows = rows or []
        self.raise_on_insert = raise_on_insert
        self.select_calls: list[dict[str, object]] = []
        self.insert_calls: list[dict[str, object]] = []

    def select_rows(self, **kwargs: object) -> dict[str, object]:
        self.select_calls.append(kwargs)
        return {"rows": self._rows}

    def insert_rows(self, **kwargs: object) -> dict[str, object]:
        self.insert_calls.append(kwargs)
        if self.raise_on_insert:
            msg = "insert failed"
            raise RuntimeError(msg)
        return {"rows": kwargs["rows"]}


def test_rows_present_true_when_any_row_returned() -> None:
    query = _FakeQuery(rows=[{"Key": 1}])
    assert (
        rows_present(query, "lists", "t", {"experiment": 7, "sample_id": "s1"}) is True
    )


def test_rows_present_false_when_empty() -> None:
    query = _FakeQuery(rows=[])
    assert (
        rows_present(query, "lists", "t", {"experiment": 7, "sample_id": "s1"}) is False
    )


def test_rows_present_builds_one_eq_filter_per_column() -> None:
    query = _FakeQuery(rows=[])
    rows_present(
        query,
        "lists",
        "t",
        {"experiment": 7, "sample_id": "s1", "query_class": "single_read"},
    )
    (call,) = query.select_calls
    filter_array = call["filter_array"]
    assert isinstance(filter_array, list)
    assert len(filter_array) == 3
    rendered = " ".join(str(f) for f in filter_array)
    assert "experiment" in rendered
    assert "sample_id" in rendered
    assert "query_class" in rendered


def test_insert_records_inserts_all_rows_in_one_call() -> None:
    query = _FakeQuery()
    records = [{"a": 1}, {"a": 2}]
    insert_records(query, "lists", "t", records)
    (call,) = query.insert_calls
    assert call["rows"] == records
    assert call["schema_name"] == "lists"
    assert call["query_name"] == "t"


def test_insert_records_propagates_api_failure() -> None:
    query = _FakeQuery(raise_on_insert=True)
    with pytest.raises(RuntimeError):
        insert_records(query, "lists", "t", [{"a": 1}])
