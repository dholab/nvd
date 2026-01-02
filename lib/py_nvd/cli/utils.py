"""
Shared utilities for the NVD CLI.

This module provides common functions, constants, and the Rich console
instance used across all CLI commands.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import typer
from rich.console import Console

from py_nvd.db import DEFAULT_CONFIG_PATH, get_config_path, get_state_db_path
from py_nvd.fingerprint import is_dev_mode, verify_pipeline

# Re-export for backward compatibility (prefer get_config_path() for new code)
DEFAULT_CONFIG = DEFAULT_CONFIG_PATH

# Resume file for `nvd resume` command (stored in launch directory)
RESUME_FILE = Path(".nfresume")

# Maximum directory depth to search when looking for pipeline root
MAX_PIPELINE_ROOT_DEPTH = 10


def _has_pipeline_files(candidate: Path) -> bool:
    """Check if a directory contains the core pipeline files."""
    return (candidate / "main.nf").exists() and (candidate / "nextflow.config").exists()


def _is_verified_nvd_pipeline(candidate: Path) -> bool:
    """
    Verify a candidate directory is the authentic NVD pipeline.

    In dev mode (editable install or NVD_DEV_MODE set), only checks for
    file presence. In production, verifies blake3 hashes match the
    fingerprint file.
    """
    if not _has_pipeline_files(candidate):
        return False

    # In dev mode, skip hash verification (files change during development)
    # In production, verify hashes match
    return verify_pipeline(candidate, strict=not is_dev_mode())


def _find_pipeline_root() -> Path:
    """
    Find the NVD pipeline root directory.

    Strategy (in order of precedence):
    1. NVD_PIPELINE_ROOT environment variable (explicit override)
    2. Current working directory (if main.nf exists there)
    3. Walk up from this file's location (development/editable install)

    Each candidate is verified against the pipeline fingerprint to ensure
    it's actually the NVD pipeline and not some other Nextflow project.
    In dev mode, verification is relaxed to allow for local edits.

    The cwd check (strategy 2) enables running the CLI from within a container
    when the repo is bind-mounted to the working directory.

    Returns:
        Path to the pipeline root directory.

    Raises:
        RuntimeError: If the pipeline root cannot be found.
    """
    # Strategy 1: Explicit override via environment variable
    if env_root := os.environ.get("NVD_PIPELINE_ROOT"):
        candidate = Path(env_root).resolve()
        if _is_verified_nvd_pipeline(candidate):
            return candidate
        # If explicitly set but invalid, warn and continue to other strategies
        # (could also raise here, but graceful fallback seems friendlier)

    # Strategy 2: Current working directory (container-friendly)
    cwd = Path.cwd()
    if _is_verified_nvd_pipeline(cwd):
        return cwd

    # Strategy 3: Walk up from installed package location (original behavior)
    current = Path(__file__).resolve().parent
    for _ in range(MAX_PIPELINE_ROOT_DEPTH):
        if _is_verified_nvd_pipeline(current):
            return current
        if current.parent == current:  # Hit filesystem root
            break
        current = current.parent

    msg = (
        "Could not find NVD pipeline root (main.nf not found or fingerprint mismatch). "
        "Either run from the repository root, set NVD_PIPELINE_ROOT, "
        "or ensure the CLI is installed from the cloned repository."
    )
    raise RuntimeError(msg)


PIPELINE_ROOT = _find_pipeline_root()


def get_editor() -> str:
    """
    Get the user's preferred text editor.

    Resolution order (following git's convention):
    1. $VISUAL
    2. $EDITOR
    3. vi (POSIX-mandated fallback)

    Returns:
        The editor command to use.
    """
    return os.environ.get("VISUAL") or os.environ.get("EDITOR") or "vi"


VALID_TOOLS = [
    "stat_blast",
    "nvd",
    "stat",
    "blast",
    "stast",
    "gottcha",
    "all",
    "clumpify",
]
MAX_PREVIEW_ITEMS = 5  # Max items to show before "... and N more"

# Help panel names for organizing --help output
PANEL_CORE = "Core Options"
PANEL_PREPROCESSING = "Read Preprocessing"
PANEL_DATABASES = "Database Paths"
PANEL_ANALYSIS = "Analysis Parameters"
PANEL_SRA = "SRA Submission"
PANEL_LABKEY = "LabKey Integration"

console = Console()


def error(message: str, exit_code: int = 1) -> None:
    """Print error message and exit."""
    console.print(f"[red]✗ Error:[/red] {message}")
    sys.exit(exit_code)


def success(message: str) -> None:
    """Print success message."""
    console.print(f"[green]✓[/green] {message}")


def info(message: str) -> None:
    """Print info message."""
    console.print(f"[cyan]i[/cyan]  {message}")


def warning(message: str) -> None:
    """Print warning message."""
    console.print(f"[yellow]⚠[/yellow]  {message}")


def ensure_db_exists(json_output: bool = False) -> Path:
    """
    Ensure the state database exists, exiting with a helpful error if not.

    This is for CLI commands that require an existing database (as opposed
    to commands like `state init` that create it).

    Args:
        json_output: If True, output error as JSON instead of Rich formatted text.

    Returns:
        Path to the database file.

    Raises:
        typer.Exit: If the database does not exist.
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        if json_output:
            console.print_json(
                json.dumps({"error": "Database not found", "path": str(db_path)}),
            )
            raise typer.Exit(1)
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)  # unreachable, but satisfies type checker
    return db_path


def find_config_file(custom_path: Path | None = None) -> Path | None:
    """
    Locate configuration file with fallback logic.

    Priority:
        1. Explicit custom_path argument
        2. NVD_CONFIG environment variable
        3. Default: ~/.nvd/user.config

    Returns None if no config file exists at the resolved path.
    """
    resolved = get_config_path(custom_path)

    if custom_path is not None and not resolved.exists():
        # Explicit path was given but doesn't exist - that's an error
        error(f"Config file not found: {custom_path}")

    if resolved.exists():
        return resolved

    return None


def parse_nextflow_config(config_path: Path) -> dict[str, str]:
    """
    Parse Nextflow config file for parameter values.
    Simple regex-based extraction of params { } block.
    """
    params: dict[str, str] = {}

    try:
        with open(config_path) as f:
            content = f.read()

        # Find params block
        params_match = re.search(r"params\s*\{([^}]+)\}", content, re.DOTALL)
        if not params_match:
            return params

        params_block = params_match.group(1)

        # Extract key = value pairs (simple string values)
        for raw_line in params_block.split("\n"):
            line = raw_line.strip()
            # Match: key = "value" or key = 'value' or key = value
            match = re.match(r'(\w+)\s*=\s*["\']?([^"\'\n]+)["\']?', line)
            if match:
                key, value = match.groups()
                params[key.strip()] = value.strip().strip("\"'")

    except OSError as e:
        warning(f"Failed to parse config file: {e}")
        return params
    else:
        return params


def check_command_exists(command: str) -> bool:
    """Check if a command exists in PATH."""
    return shutil.which(command) is not None


def docker_is_running() -> bool:
    """Check if Docker daemon is running."""
    try:
        result = subprocess.run(
            ["docker", "info"],  # noqa: S607
            capture_output=True,
            check=True,
            timeout=2,
        )
    except (
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
        FileNotFoundError,
    ):
        return False
    else:
        return result.returncode == 0


def auto_detect_profile() -> str:
    """Auto-detect available execution environment."""
    if check_command_exists("docker") and docker_is_running():
        return "docker"
    if check_command_exists("apptainer"):
        return "apptainer"
    if check_command_exists("singularity"):
        return "apptainer"  # Singularity is legacy name for Apptainer
    return "local"
