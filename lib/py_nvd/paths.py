"""Path and environment-variable resolution for NVD local resources."""

from __future__ import annotations

import os
from pathlib import Path

# The NVD config directory: user-facing configuration and small control-plane
# data live here. The repository installer may still use ~/.nvd for versioned
# checkouts; that install location is deliberately separate from this runtime
# config root.
DEFAULT_CONFIG_DIR = Path.home() / ".config" / "nvd"

# Legacy install location retained only as a read fallback for user.config.
LEGACY_NVD_HOME = Path.home() / ".nvd"
DEFAULT_NVD_HOME = LEGACY_NVD_HOME

# Environment variables for overriding default paths.
ENV_VAR_CONFIG = "NVD_CONFIG"
ENV_VAR_CONFIG_DIR = "NVD_CONFIG_DIR"
ENV_VAR_PRESET_STORE = "NVD_PRESET_STORE"
ENV_VAR_TAXONOMY = "NVD_TAXONOMY_DB"

# Default paths.
DEFAULT_CONFIG_PATH = DEFAULT_CONFIG_DIR / "user.config"


def get_config_dir(explicit_path: Path | str | None = None) -> Path:
    """Resolve the NVD configuration directory.

    Priority:
        1. Explicit path argument
        2. NVD_CONFIG_DIR environment variable
        3. Default: ~/.config/nvd

    The directory is created because this is the canonical home for setup.conf,
    user.config, shell completions, and the temporary SQLite preset store.
    """
    if explicit_path is not None:
        config_dir = Path(explicit_path)
    elif ENV_VAR_CONFIG_DIR in os.environ:
        config_dir = Path(os.environ[ENV_VAR_CONFIG_DIR])
    else:
        config_dir = DEFAULT_CONFIG_DIR

    config_dir.mkdir(parents=True, exist_ok=True)
    return config_dir


def get_config_path(explicit_path: Path | str | None = None) -> Path:
    """
    Resolve the NVD config file path using hierarchical fallback.

    Priority:
        1. Explicit path argument (from CLI)
        2. NVD_CONFIG environment variable
        3. {NVD_CONFIG_DIR}/user.config, defaulting to ~/.config/nvd/user.config
        4. Legacy fallback: ~/.nvd/user.config if it already exists

    Unlike state/taxonomy directories, this does not create the file if it does
    not exist. Config files should be created explicitly.
    """
    if explicit_path is not None:
        return Path(explicit_path)
    if ENV_VAR_CONFIG in os.environ:
        return Path(os.environ[ENV_VAR_CONFIG])
    config_path = get_config_dir() / "user.config"
    legacy_path = LEGACY_NVD_HOME / "user.config"
    if not config_path.exists() and legacy_path.exists():
        return legacy_path
    return config_path


def get_taxdump_dir(
    taxonomy_dir: Path | str | None = None,
) -> Path:
    """
    Get the taxdump directory containing .dmp files and taxonomy.sqlite.

    Priority:
        1. Explicit taxonomy_dir argument (from CLI or Nextflow param)
        2. NVD_TAXONOMY_DB environment variable (for shared cluster installs)
        3. {NVD_CONFIG_DIR}/taxdump, defaulting to ~/.config/nvd/taxdump
    """
    if taxonomy_dir is not None:
        return Path(taxonomy_dir)
    if ENV_VAR_TAXONOMY in os.environ:
        return Path(os.environ[ENV_VAR_TAXONOMY])
    return get_config_dir() / "taxdump"


def get_setup_conf_path() -> Path:
    """Return the canonical setup.conf path."""
    return get_config_dir() / "setup.conf"


def get_preset_db_path(explicit_path: Path | str | None = None) -> Path:
    """Return the SQLite preset store path.

    Priority:
        1. Explicit path argument
        2. NVD_PRESET_STORE environment variable (for shared mutable stores)
        3. {NVD_CONFIG_DIR}/presets.sqlite, defaulting to ~/.config/nvd/presets.sqlite
    """
    if explicit_path is not None:
        return Path(explicit_path)
    if ENV_VAR_PRESET_STORE in os.environ:
        return Path(os.environ[ENV_VAR_PRESET_STORE])
    return get_config_dir() / "presets.sqlite"


def get_taxonomy_db_path() -> Path:
    """Get the taxonomy SQLite database path inside the taxdump directory."""
    return get_taxdump_dir() / "taxonomy.sqlite"
