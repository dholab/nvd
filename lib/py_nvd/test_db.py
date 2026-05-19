"""Path and resource-warning tests retained during the v3 state removal."""

from __future__ import annotations

from typing import TYPE_CHECKING

from py_nvd.db import StateUnavailableError, format_state_warning
from py_nvd.paths import (
    get_config_dir,
    get_config_path,
    get_preset_db_path,
    get_setup_conf_path,
    get_state_dir,
    get_taxdump_dir,
)

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def test_state_unavailable_error_keeps_db_path_alias(tmp_path: Path) -> None:
    db_path = tmp_path / "legacy-state.sqlite"
    error = StateUnavailableError(db_path, "operation", "reason")
    assert error.db_path == db_path
    assert error.resource_path == db_path


def test_format_state_warning_includes_context(tmp_path: Path) -> None:
    warning = format_state_warning(
        operation="Legacy operation",
        context="testing",
        error=PermissionError("denied"),
        db_path=tmp_path / "legacy-state.sqlite",
        consequences=["legacy coordination unavailable"],
    )
    assert "Legacy operation" in warning
    assert "legacy coordination unavailable" in warning


def test_config_dir_env_sets_config_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("NVD_CONFIG_DIR", str(tmp_path))
    monkeypatch.delenv("NVD_CONFIG", raising=False)

    assert get_config_dir() == tmp_path
    assert get_config_path() == tmp_path / "user.config"
    assert get_setup_conf_path() == tmp_path / "setup.conf"
    assert get_preset_db_path() == tmp_path / "presets.sqlite"
    assert get_taxdump_dir() == tmp_path / "taxdump"


def test_explicit_config_path_wins(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("NVD_CONFIG_DIR", str(tmp_path / "root"))
    explicit = tmp_path / "custom.config"
    assert get_config_path(explicit) == explicit


def test_nvd_config_file_env_wins(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_file = tmp_path / "env.config"
    monkeypatch.setenv("NVD_CONFIG", str(config_file))
    assert get_config_path() == config_file


def test_legacy_state_dir_alias_still_resolves_local_data(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    legacy_state = tmp_path / "legacy"
    monkeypatch.setenv("NVD_STATE_DIR", str(legacy_state))
    assert get_state_dir() == legacy_state
    assert get_taxdump_dir() == legacy_state / "taxdump"


def test_taxonomy_env_wins_over_config_dir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    taxonomy_dir = tmp_path / "shared-taxdump"
    monkeypatch.setenv("NVD_CONFIG_DIR", str(tmp_path / "config"))
    monkeypatch.setenv("NVD_TAXONOMY_DB", str(taxonomy_dir))
    assert get_taxdump_dir() == taxonomy_dir
