"""Tests for ensure_taxonomy.py preflight policy behavior."""

import sys
from pathlib import Path

import pytest


def test_read_only_mode_fails_fast_when_taxonomy_is_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Explicit read-only preflight failures should not be downgraded to warnings."""
    from ensure_taxonomy import main

    taxonomy_dir = tmp_path / "taxdump"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ensure_taxonomy.py",
            "--taxonomy-dir",
            str(taxonomy_dir),
            "--taxonomy-mode",
            "read_only",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        main()

    assert exc_info.value.code == 1
    assert not (taxonomy_dir / "taxonomy.sqlite").exists()


def test_legacy_offline_env_fails_fast_when_taxonomy_is_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Legacy offline preflight should still require prepared taxonomy."""
    from ensure_taxonomy import main

    taxonomy_dir = tmp_path / "taxdump"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ensure_taxonomy.py",
            "--taxonomy-dir",
            str(taxonomy_dir),
        ],
    )
    monkeypatch.setenv("NVD_TAXONOMY_OFFLINE", "1")

    with pytest.raises(SystemExit) as exc_info:
        main()

    assert exc_info.value.code == 1
    assert not (taxonomy_dir / "taxonomy.sqlite").exists()
