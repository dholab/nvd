"""
Version command for the NVD CLI.

Commands:
    nvd version  - Show version information and dependency status
"""

from __future__ import annotations

import subprocess

from rich.panel import Panel

from py_nvd import __version__
from py_nvd.cli.utils import (
    check_command_exists,
    console,
    docker_is_running,
)


def version() -> None:  # noqa: PLR0912
    """Show version information."""
    console.print(
        Panel.fit(
            f"[bold cyan]NVD2 CLI[/bold cyan]\n"
            f"Version: {__version__}\n\n"
            f"A user-friendly wrapper for the NVD2 Nextflow pipeline.\n"
            f"Pipeline: dhoconno/nvd",
            title="Version Info",
            border_style="cyan",
        ),
    )

    # Show dependency versions if available
    deps = []

    if check_command_exists("nextflow"):
        try:
            result = subprocess.run(
                ["nextflow", "-version"],  # noqa: S607
                check=False,
                capture_output=True,
                text=True,
                timeout=3,
            )
            nf_version = result.stdout.strip().split("\n")[0]
            deps.append(f"✓ {nf_version}")
        except (OSError, subprocess.TimeoutExpired):
            deps.append("✓ Nextflow (version unknown)")
    else:
        deps.append("✗ Nextflow not found")

    if check_command_exists("java"):
        try:
            result = subprocess.run(
                ["java", "-version"],  # noqa: S607
                check=False,
                capture_output=True,
                text=True,
                timeout=3,
            )
            java_output = result.stderr.strip().split("\n")[0]
            deps.append(f"✓ {java_output}")
        except (OSError, subprocess.TimeoutExpired):
            deps.append("✓ Java (version unknown)")
    else:
        deps.append("✗ Java not found")

    # Docker/Apptainer
    if check_command_exists("docker"):
        if docker_is_running():
            deps.append("✓ Docker (running)")
        else:
            deps.append("⚠ Docker (daemon not running)")
    elif check_command_exists("apptainer"):
        deps.append("✓ Apptainer")
    elif check_command_exists("singularity"):
        deps.append("✓ Singularity")
    else:
        deps.append("⚠ No container runtime found")

    console.print("\n[bold]Dependencies:[/bold]")
    for dep in deps:
        console.print(f"  {dep}")

    console.print()
