"""
Shared utilities for the NVD CLI.

This module provides common functions, constants, and the Rich console
instance used across all CLI commands.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import NoReturn

import typer
from rich.console import Console

from py_nvd.db import get_state_db_path
from py_nvd.paths import DEFAULT_CONFIG_PATH, get_config_path

# Re-export for backward compatibility (prefer get_config_path() for new code)
DEFAULT_CONFIG = DEFAULT_CONFIG_PATH

# Resume file for `nvd resume` command (stored in launch directory)
RESUME_FILE = Path(".nfresume")

# Files and directories that distinguish NVD from an arbitrary Nextflow project.
NVD_PIPELINE_SENTINELS = (
    "main.nf",
    "nextflow.config",
    "pyproject.toml",
    "lib/py_nvd",
    "modules",
    "subworkflows",
)

# Time constants for duration formatting
SECONDS_PER_MINUTE = 60
SECONDS_PER_HOUR = 3600


def _has_nvd_sentinels(candidate: Path) -> bool:
    """Check whether a directory looks like an NVD source checkout."""
    return all((candidate / sentinel).exists() for sentinel in NVD_PIPELINE_SENTINELS)


def _find_parent_pipeline_root(start: Path) -> Path | None:
    """Walk upward from a path until an NVD pipeline root is found."""
    current = start.resolve()
    if current.is_file():
        current = current.parent

    while True:
        if _has_nvd_sentinels(current):
            return current
        if current.parent == current:
            break
        current = current.parent

    return None


def _pipeline_root_error() -> RuntimeError:
    """Build the standard missing pipeline root error."""
    msg = (
        "Could not find NVD pipeline root. "
        "Either run from the repository root, set NVD_PIPELINE_ROOT, "
        "or ensure the CLI is installed from the cloned repository."
    )
    return RuntimeError(msg)


def _find_pipeline_root() -> Path:
    """
    Find the NVD pipeline root directory.

    Strategy (in order of precedence):
    1. NVD_PIPELINE_ROOT environment variable (explicit override)
    2. Walk up from this file's location (development/editable install)
    3. Current working directory (container bind-mount fallback)

    Candidates must contain NVD-specific sentinel paths so the CLI does not
    accidentally run an unrelated Nextflow project that merely has main.nf and
    nextflow.config.

    Returns:
        Path to the pipeline root directory.

    Raises:
        RuntimeError: If the pipeline root cannot be found.
    """
    # Strategy 1: Explicit override via environment variable
    if env_root := os.environ.get("NVD_PIPELINE_ROOT"):
        candidate = Path(env_root).resolve()
        if _has_nvd_sentinels(candidate):
            return candidate

        msg = (
            "NVD_PIPELINE_ROOT is set, but does not point to an NVD pipeline root: "
            f"{candidate}"
        )
        raise RuntimeError(msg)

    # Strategy 2: Walk up from installed package location.
    package_root = _find_parent_pipeline_root(Path(__file__))
    if package_root is not None:
        return package_root

    # Strategy 3: Current working directory fallback for container bind mounts.
    cwd_root = _find_parent_pipeline_root(Path.cwd())
    if cwd_root is not None:
        return cwd_root

    raise _pipeline_root_error()


def _find_pipeline_root_optional() -> Path | None:
    """Find the pipeline root when available without making CLI imports fail."""
    try:
        return _find_pipeline_root()
    except RuntimeError as exc:
        if os.environ.get("NVD_PIPELINE_ROOT"):
            raise
        if str(exc) != str(_pipeline_root_error()):
            raise
        return None


def get_pipeline_root() -> Path:
    """Return the pipeline root, failing only when a command actually needs it."""
    if PIPELINE_ROOT is not None:
        return PIPELINE_ROOT
    return _find_pipeline_root()


PIPELINE_ROOT = _find_pipeline_root_optional()


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


MAX_PREVIEW_ITEMS = 5  # Max items to show before "... and N more"

# Help panel names for organizing --help output
PANEL_CORE = "Core Options"
PANEL_PREPROCESSING = "Read Preprocessing"
PANEL_DATABASES = "Database Paths"
PANEL_ANALYSIS = "Analysis Parameters"
PANEL_LABKEY = "LabKey Integration"
PANEL_NOTIFICATIONS = "Notifications"

console = Console()


def error(message: str, exit_code: int = 1) -> NoReturn:
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


def format_duration(seconds: float) -> str:
    """
    Format seconds as human-readable duration.

    Truncates fractional seconds (no rounding). For durations >= 1 hour,
    seconds are omitted for readability.

    Examples:
        format_duration(45) -> "45s"
        format_duration(125) -> "2m 5s"
        format_duration(3725) -> "1h 2m"
        format_duration(7200) -> "2h 0m"

    Args:
        seconds: Duration in seconds (must be non-negative and finite).

    Returns:
        Human-readable duration string.

    Raises:
        AssertionError: If seconds is negative, NaN, or infinite.
    """
    assert not math.isnan(seconds), "duration cannot be NaN"
    assert not math.isinf(seconds), "duration cannot be infinite"
    assert seconds >= 0, f"duration must be non-negative, got {seconds}"

    if seconds < SECONDS_PER_MINUTE:
        return f"{int(seconds)}s"
    if seconds < SECONDS_PER_HOUR:
        minutes = int(seconds // SECONDS_PER_MINUTE)
        secs = int(seconds % SECONDS_PER_MINUTE)
        return f"{minutes}m {secs}s"
    hours = int(seconds // SECONDS_PER_HOUR)
    minutes = int((seconds % SECONDS_PER_HOUR) // SECONDS_PER_MINUTE)
    return f"{hours}h {minutes}m"


def ensure_db_exists(json_output: bool = False) -> Path:  # noqa: FBT001, FBT002
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


def get_default_profile() -> str | None:
    """
    Read default profile from ~/.nvd/setup.conf if configured.

    Returns:
        The configured default profile, or None if not set.
    """
    setup_conf = Path.home() / ".nvd" / "setup.conf"
    if not setup_conf.exists():
        return None

    for raw_line in setup_conf.read_text().splitlines():
        line = raw_line.strip()
        if line.startswith("NVD_DEFAULT_PROFILE="):
            value = line.split("=", 1)[1].strip()
            if value:
                return value
    return None


def auto_detect_profile() -> str:
    """Auto-detect available execution environment."""
    if check_command_exists("docker") and docker_is_running():
        return "docker"
    if check_command_exists("apptainer"):
        return "apptainer"
    if check_command_exists("singularity"):
        return "apptainer"  # Singularity is legacy name for Apptainer
    return "local"


_MIN_CMD_PARTS_FOR_CONTINUATION = 3


def format_command_for_display(cmd: list[str]) -> str:
    """
    Format a command list for readable display with line continuations.

    Each argument pair (--flag value) gets its own line, indented and
    with shell continuation characters for copy-paste compatibility.
    """
    if len(cmd) < _MIN_CMD_PARTS_FOR_CONTINUATION:
        return " ".join(cmd)

    lines = []
    # First line: nextflow run <pipeline_root>
    lines.append(f"{cmd[0]} {cmd[1]} {cmd[2]} \\")

    # Process remaining args in pairs where possible
    i = 3
    while i < len(cmd):
        arg = cmd[i]

        # Check if this is a flag that takes a value (not a standalone flag like -resume)
        if i + 1 < len(cmd) and not cmd[i + 1].startswith("-"):
            # Flag with value: --param value
            value = cmd[i + 1]
            if i + 2 < len(cmd):
                lines.append(f"    {arg} {value} \\")
            else:
                lines.append(f"    {arg} {value}")
            i += 2
        else:
            # Standalone flag like -resume
            if i + 1 < len(cmd):
                lines.append(f"    {arg} \\")
            else:
                lines.append(f"    {arg}")
            i += 1

    return "\n".join(lines)
