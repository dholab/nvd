"""Shared LabKey list I/O helpers for the BLAST hits and FASTA upload scripts.

Both upload scripts dedup against their destination LabKey list (the list is its
own completion ledger) and insert each unit atomically. This module centralizes
the two pieces of client logic they would otherwise duplicate:

- ``rows_present`` — a presence check with the ``QueryFilter`` -> filter-string
  fallback so it still works against older labkey-api-python versions.
- ``insert_records`` — a single-call insert that deliberately lets failures
  propagate so callers can hard-fail rather than log-and-continue.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from collections.abc import Mapping


class LabKeyQueryApi(Protocol):
    """The subset of the LabKey APIWrapper ``query`` interface used here."""

    def select_rows(self, **kwargs: object) -> dict[str, object]: ...

    def insert_rows(self, **kwargs: object) -> object: ...


def rows_present(
    query_api: LabKeyQueryApi,
    schema: str,
    table: str,
    filters: Mapping[str, object],
) -> bool:
    """Return whether the list holds at least one row matching every filter.

    ``filters`` maps column name to the value it must equal; all are ANDed. The
    destination list is its own completion ledger, so any matching row means the
    unit was already uploaded and the caller should skip re-inserting it.

    Uses ``QueryFilter`` objects when available, falling back to filter-string
    syntax (``col~eq=value``) on labkey-api-python versions that lack them.
    """
    try:
        # Imported lazily so the filter-string fallback can handle older
        # labkey-api-python versions that do not ship QueryFilter.
        from labkey.query import QueryFilter  # noqa: PLC0415

        result = query_api.select_rows(
            schema_name=schema,
            query_name=table,
            filter_array=[
                QueryFilter(column, value, "eq") for column, value in filters.items()
            ],
            max_rows=1,
        )
    except (ImportError, AttributeError):
        result = query_api.select_rows(
            schema_name=schema,
            query_name=table,
            filter_array=[f"{column}~eq={value}" for column, value in filters.items()],
        )
    return bool(result and result.get("rows"))


def insert_records(
    query_api: LabKeyQueryApi,
    schema: str,
    table: str,
    records: list[dict[str, object]],
) -> None:
    """Insert every record in one atomic call.

    Deliberately does not catch anything: a failed insert must propagate so the
    caller can hard-fail the run rather than log an error and report success.
    """
    query_api.insert_rows(
        schema_name=schema,
        query_name=table,
        rows=records,
    )
