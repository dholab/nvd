"""
Shared utilities for the NVD CLI.

This module provides common functions, constants, and the Rich console
instance used across all CLI commands.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

from rich.console import Console

# ============================================================================
# CONSTANTS
# ============================================================================

DEFAULT_CONFIG = Path.home() / ".nvd2" / "config" / "user.config"
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
    """Locate configuration file with fallback logic."""
    if custom_path:
        if custom_path.exists():
            return custom_path
        error(f"Config file not found: {custom_path}")

    if DEFAULT_CONFIG.exists():
        return DEFAULT_CONFIG

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
