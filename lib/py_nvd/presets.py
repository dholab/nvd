"""Preset storage for reusable NVD run parameter bundles."""

from __future__ import annotations

import json
import sqlite3
from datetime import UTC, datetime
from pathlib import Path
from typing import Protocol

from py_nvd.models import Preset
from py_nvd.paths import get_preset_db_path


def utc_now_iso() -> str:
    """Return the current UTC time as an ISO8601 string with a Z suffix."""
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


class PresetStore(Protocol):
    """Storage boundary for named parameter presets."""

    def get(self, name: str) -> Preset | None: ...

    def list(self) -> list[Preset]: ...

    def upsert(
        self,
        name: str,
        params: dict[str, object],
        description: str | None = None,
    ) -> Preset: ...

    def delete(self, name: str) -> bool: ...


class SQLitePresetStore:
    """SQLite-backed preset store using only a narrow presets schema."""

    def __init__(self, db_path: Path | str | None = None) -> None:
        self.db_path = Path(db_path) if db_path is not None else get_preset_db_path()

    def _connect(self) -> sqlite3.Connection:
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS presets (
                name TEXT PRIMARY KEY,
                params TEXT NOT NULL,
                description TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
            """,
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_presets_created_at ON presets(created_at)",
        )
        return conn

    def get(self, name: str) -> Preset | None:
        with self._connect() as conn:
            row = conn.execute(
                "SELECT * FROM presets WHERE name = ?",
                (name,),
            ).fetchone()
        return Preset(**dict(row)) if row is not None else None

    def list(self) -> list[Preset]:
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT * FROM presets ORDER BY created_at DESC",
            ).fetchall()
        return [Preset(**dict(row)) for row in rows]

    def upsert(
        self,
        name: str,
        params: dict[str, object],
        description: str | None = None,
    ) -> Preset:
        params_json = json.dumps(params, sort_keys=True)
        now = utc_now_iso()

        with self._connect() as conn:
            existing = conn.execute(
                "SELECT created_at FROM presets WHERE name = ?",
                (name,),
            ).fetchone()
            if existing is None:
                conn.execute(
                    """
                    INSERT INTO presets (name, params, description, created_at, updated_at)
                    VALUES (?, ?, ?, ?, ?)
                    """,
                    (name, params_json, description, now, now),
                )
            else:
                conn.execute(
                    """
                    UPDATE presets
                    SET params = ?, description = ?, updated_at = ?
                    WHERE name = ?
                    """,
                    (params_json, description, now, name),
                )

        preset = self.get(name)
        assert preset is not None, "preset upsert should make preset readable"
        return preset

    def delete(self, name: str) -> bool:
        with self._connect() as conn:
            cursor = conn.execute("DELETE FROM presets WHERE name = ?", (name,))
            return cursor.rowcount > 0


def get_preset_store() -> PresetStore:
    """Return the configured preset store.

    The store is SQLite-backed for now, but callers depend only on the small
    PresetStore interface so a Turso-backed implementation can replace it later.
    """
    return SQLitePresetStore()
