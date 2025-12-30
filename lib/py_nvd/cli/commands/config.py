# ruff: noqa: B008
"""
Configuration management commands for the NVD CLI.

Commands:
    nvd config show  - Display current configuration
    nvd config path  - Show configuration file location
"""

from __future__ import annotations

import sys
from pathlib import Path

import typer
from rich.table import Table

from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    console,
    find_config_file,
    info,
    parse_nextflow_config,
    success,
    warning,
)

config_app = typer.Typer(
    name="config",
    help="Configuration management commands",
    no_args_is_help=True,
)


@config_app.command("show")
@config_app.command("s", hidden=True)  # Alias
def config_show(
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Config file to display (default: ~/.nvd2/config/user.config)",
    ),
) -> None:
    """Display current configuration."""
    # Find config
    config_file = find_config_file(config)
    if not config_file:
        warning(f"No config file found at {DEFAULT_CONFIG}")
        info("Run install.sh to create a configuration file")
        sys.exit(1)

    console.print(f"\n[bold]Configuration File:[/bold] {config_file}\n")

    # Parse config
    params = parse_nextflow_config(config_file)

    if not params:
        warning("No parameters found in config file")
        sys.exit(1)

    # Display as table
    table = Table(
        title="Configuration Parameters",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Parameter", style="cyan", width=30)
    table.add_column("Value", style="white")

    for key, value in sorted(params.items()):
        table.add_row(key, value)

    console.print(table)
    console.print()


@config_app.command("path")
@config_app.command("p", hidden=True)  # Alias
def config_path() -> None:
    """Show configuration file location."""
    console.print(f"\n[bold]Default config location:[/bold] {DEFAULT_CONFIG}")

    if DEFAULT_CONFIG.exists():
        success(f"Config file exists: {DEFAULT_CONFIG}")
        file_size = DEFAULT_CONFIG.stat().st_size
        console.print(f"  Size: {file_size} bytes")
    else:
        warning("Config file does not exist")
        info("Run install.sh to create a configuration file")

    console.print()
