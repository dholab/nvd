"""
Shared utilities for the NVD CLI.

This module provides common functions, constants, and the Rich console
instance used across all CLI commands.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from rich.console import Console

from py_nvd.db import DEFAULT_CONFIG_PATH, get_config_path

# ============================================================================
# CONSTANTS
# ============================================================================

# Re-export for backward compatibility (prefer get_config_path() for new code)
DEFAULT_CONFIG = DEFAULT_CONFIG_PATH

# Resume file for `nvd resume` command (stored in launch directory)
RESUME_FILE = Path(".nfresume")


def _find_pipeline_root() -> Path:
    """
    Find the NVD pipeline root directory.

    Strategy:
    1. Walk up from this file's location looking for main.nf
    2. Verify nextflow.config also exists (sanity check)

    This is more robust than counting parent directories, which breaks
    if the package structure changes.

    Returns:
        Path to the pipeline root directory.

    Raises:
        RuntimeError: If the pipeline root cannot be found.
    """
    current = Path(__file__).resolve().parent
    max_depth = 10  # Safety limit to avoid infinite loops

    for _ in range(max_depth):
        if (current / "main.nf").exists() and (current / "nextflow.config").exists():
            return current
        if current.parent == current:  # Hit filesystem root
            break
        current = current.parent

    # This should never happen if installed correctly
    msg = (
        "Could not find NVD pipeline root (main.nf not found). "
        "Ensure the CLI is installed from the cloned repository."
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
VALID_PROFILES = ["docker", "apptainer", "local", "chtc_hpc"]
MAX_PREVIEW_ITEMS = 5  # Max items to show before "... and N more"

# Help panel names for organizing --help output
PANEL_CORE = "Core Options"
PANEL_PREPROCESSING = "Read Preprocessing"
PANEL_DATABASES = "Database Paths"
PANEL_ANALYSIS = "Analysis Parameters"
PANEL_SRA = "SRA Submission"
PANEL_LABKEY = "LabKey Integration"

# ============================================================================
# CONSOLE AND OUTPUT UTILITIES
# ============================================================================

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


# ============================================================================
# CONFIG MANAGEMENT
# ============================================================================


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


# ============================================================================
# ENVIRONMENT DETECTION
# ============================================================================


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
