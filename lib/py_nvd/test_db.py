"""Path resolution tests retained after the v3 state removal."""

from __future__ import annotations

from typing import TYPE_CHECKING

from py_nvd.paths import (
    get_config_dir,
    get_config_path,
    get_preset_db_path,
    get_setup_conf_path,
    get_taxdump_dir,
)

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def test_config_dir_env_sets_config_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("NVD_CONFIG_DIR", str(tmp_path))
    monkeypatch.delenv("NVD_CONFIG", raising=False)
    monkeypatch.delenv("NVD_PRESET_STORE", raising=False)

    assert get_config_dir() == tmp_path
    assert get_config_path() == tmp_path / "user.config"
    assert get_setup_conf_path() == tmp_path / "setup.conf"
    assert get_preset_db_path() == tmp_path / "presets.sqlite"
    assert get_taxdump_dir() == tmp_path / "taxdump"


def test_preset_store_env_wins_over_config_dir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_root = tmp_path / "config"
    shared_store = tmp_path / "shared" / "presets.sqlite"
    monkeypatch.setenv("NVD_CONFIG_DIR", str(config_root))
    monkeypatch.setenv("NVD_PRESET_STORE", str(shared_store))

    assert get_preset_db_path() == shared_store
    assert get_config_path() == config_root / "user.config"
    assert get_setup_conf_path() == config_root / "setup.conf"
    assert get_taxdump_dir() == config_root / "taxdump"


def test_explicit_preset_store_path_wins(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("NVD_PRESET_STORE", str(tmp_path / "env.sqlite"))
    explicit = tmp_path / "explicit.sqlite"

    assert get_preset_db_path(explicit) == explicit


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


def test_state_dir_env_does_not_affect_taxonomy_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_root = tmp_path / "config"
    old_state = tmp_path / "old-state"
    monkeypatch.setenv("NVD_CONFIG_DIR", str(config_root))
    monkeypatch.setenv("NVD_STATE_DIR", str(old_state))

    assert get_taxdump_dir() == config_root / "taxdump"


def test_taxonomy_env_wins_over_config_dir(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    taxonomy_dir = tmp_path / "shared-taxdump"
    monkeypatch.setenv("NVD_CONFIG_DIR", str(tmp_path / "config"))
    monkeypatch.setenv("NVD_TAXONOMY_DB", str(taxonomy_dir))
    assert get_taxdump_dir() == taxonomy_dir
