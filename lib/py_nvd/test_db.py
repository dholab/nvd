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

# Every environment variable that path resolution consults. A developer machine
# with a real NVD install exports several of these, and they would otherwise
# override the values a test monkeypatches in.
PATH_RESOLUTION_ENV_VARS = (
    "NVD_CONFIG",
    "NVD_CONFIG_DIR",
    "NVD_PRESET_STORE",
    "NVD_STATE_DIR",
    "NVD_TAXONOMY_DB",
)


@pytest.fixture(autouse=True)
def isolate_path_resolution_from_the_developer_machine(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Make these tests depend only on what each one sets for itself.

    Two things on a real machine leak into path resolution. The first is the
    environment: a developer who actually runs NVD exports NVD_TAXONOMY_DB and
    friends, which outrank a monkeypatched NVD_CONFIG_DIR. The second is the
    legacy home directory: ``get_config_path`` falls back to
    ``~/.nvd/user.config`` whenever that file exists and the configured one does
    not, so anyone upgraded from an older NVD gets the legacy path instead of
    the one under test. Pointing the legacy home at an empty temporary
    directory switches that fallback off without disabling the code path, which
    ``test_legacy_config_is_used_when_no_configured_file_exists`` still covers.
    """
    for variable_name in PATH_RESOLUTION_ENV_VARS:
        monkeypatch.delenv(variable_name, raising=False)
    monkeypatch.setattr("py_nvd.paths.LEGACY_NVD_HOME", tmp_path / "absent-legacy-home")


def test_config_dir_env_sets_config_root(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("NVD_CONFIG_DIR", str(tmp_path))

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


def test_legacy_config_is_used_when_no_configured_file_exists(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An upgraded install keeps reading ~/.nvd/user.config until it moves.

    This is the fallback that the autouse fixture switches off for every other
    test in this file, so it is pinned here explicitly.
    """
    legacy_home = tmp_path / "legacy-home"
    legacy_home.mkdir()
    legacy_config = legacy_home / "user.config"
    legacy_config.write_text("params {}\n", encoding="utf-8")
    monkeypatch.setattr("py_nvd.paths.LEGACY_NVD_HOME", legacy_home)

    configured_root = tmp_path / "config"
    monkeypatch.setenv("NVD_CONFIG_DIR", str(configured_root))

    assert get_config_path() == legacy_config

    # Once the configured file exists it wins, so the fallback is not sticky.
    # The directory is already there: get_config_dir() creates it as a side
    # effect of the call above, which is why exist_ok is needed here.
    configured_root.mkdir(exist_ok=True)
    (configured_root / "user.config").write_text("params {}\n", encoding="utf-8")
    assert get_config_path() == configured_root / "user.config"
