"""Path and environment-variable resolution for NVD local resources."""

from __future__ import annotations

import os
from pathlib import Path

# The NVD config directory: user-facing configuration lives here.
DEFAULT_NVD_HOME = Path.home() / ".nvd"

# The NVD cache directory: mutable local data lives here. The name is retained
# as state for compatibility with existing parameters and environment variables,
# even as v3 removes the old state database product.
DEFAULT_NVD_CACHE = Path.home() / ".cache" / "nvd"

# Environment variables for overriding default paths.
ENV_VAR_CONFIG = "NVD_CONFIG"
ENV_VAR = "NVD_STATE_DIR"  # Keep short name for backward compatibility.
ENV_VAR_TAXONOMY = "NVD_TAXONOMY_DB"

# Default paths.
DEFAULT_CONFIG_PATH = DEFAULT_NVD_HOME / "user.config"
DEFAULT_STATE_DIR = DEFAULT_NVD_CACHE


def get_config_path(explicit_path: Path | str | None = None) -> Path:
    """
    Resolve the NVD config file path using hierarchical fallback.

    Priority:
        1. Explicit path argument (from CLI)
        2. NVD_CONFIG environment variable
        3. Default: ~/.nvd/user.config

    Unlike state/taxonomy directories, this does not create the file if it does
    not exist. Config files should be created explicitly.
    """
    if explicit_path is not None:
        return Path(explicit_path)
    if ENV_VAR_CONFIG in os.environ:
        return Path(os.environ[ENV_VAR_CONFIG])
    return DEFAULT_CONFIG_PATH


def get_state_dir(explicit_path: Path | str | None = None) -> Path:
    """
    Resolve the NVD state/local-data directory using hierarchical fallback.

    Priority:
        1. Explicit path argument (from CLI or Nextflow param)
        2. NVD_STATE_DIR environment variable
        3. Default: ~/.cache/nvd/

    The directory is created if it does not exist.
    """
    if explicit_path is not None:
        state_dir = Path(explicit_path)
    elif ENV_VAR in os.environ:
        state_dir = Path(os.environ[ENV_VAR])
    else:
        state_dir = DEFAULT_STATE_DIR

    state_dir.mkdir(parents=True, exist_ok=True)
    return state_dir


def get_taxdump_dir(
    state_dir: Path | str | None = None,
    taxonomy_dir: Path | str | None = None,
) -> Path:
    """
    Get the taxdump directory containing .dmp files and taxonomy.sqlite.

    Priority:
        1. Explicit taxonomy_dir argument (from CLI or Nextflow param)
        2. NVD_TAXONOMY_DB environment variable (for shared cluster installs)
        3. {state_dir}/taxdump (default)
    """
    if taxonomy_dir is not None:
        return Path(taxonomy_dir)
    if ENV_VAR_TAXONOMY in os.environ:
        return Path(os.environ[ENV_VAR_TAXONOMY])
    return get_state_dir(state_dir) / "taxdump"


def get_taxonomy_db_path(state_dir: Path | str | None = None) -> Path:
    """Get the taxonomy SQLite database path inside the taxdump directory."""
    return get_taxdump_dir(state_dir) / "taxonomy.sqlite"
