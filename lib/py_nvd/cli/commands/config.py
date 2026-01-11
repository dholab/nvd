# ruff: noqa: B008
"""
Configuration management commands for the NVD CLI.

Commands:
    nvd config show  - Display current configuration
    nvd config path  - Show configuration file location
    nvd config edit  - Open configuration file in editor
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import typer
from rich.table import Table

from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    console,
    error,
    find_config_file,
    get_editor,
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
        help="Config file to display (default: ~/.nvd/user.config or NVD_CONFIG)",
    ),
    params: bool = typer.Option(
        False,
        "--params",
        "-p",
        help="Parse and display parameters as a table instead of raw file",
    ),
) -> None:
    """
    Display current configuration.

    By default, prints the raw config file contents. Use --params to parse
    and display as a formatted table of parameters.
    """
    config_file = find_config_file(config)
    if not config_file:
        warning(f"No config file found at {DEFAULT_CONFIG}")
        info("Run 'nvd setup' to create a configuration file")
        sys.exit(1)

    console.print(f"\n[bold]Configuration File:[/bold] {config_file}\n")

    if params:
        # Parse and display as table
        parsed = parse_nextflow_config(config_file)

        if not parsed:
            warning("No parameters found in config file")
            sys.exit(1)

        table = Table(
            title="Configuration Parameters",
            show_header=True,
            header_style="bold cyan",
        )
        table.add_column("Parameter", style="cyan", width=30)
        table.add_column("Value", style="white")

        for key, value in sorted(parsed.items()):
            table.add_row(key, value)

        console.print(table)
    else:
        # Print raw file contents
        console.print(config_file.read_text())

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


@config_app.command("edit")
@config_app.command("e", hidden=True)  # Alias
def config_edit(
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Config file to edit (default: ~/.nvd/user.config or NVD_CONFIG)",
    ),
) -> None:
    """
    Open configuration file in editor.

    Uses $VISUAL if set, otherwise $EDITOR, falling back to vi.

    Examples:

        # Edit default config
        nvd config edit

        # Edit specific config file
        nvd config edit --config /path/to/config
    """
    config_file = find_config_file(config)
    if not config_file:
        error(f"No config file found at {DEFAULT_CONFIG}")
        info("Run 'nvd setup' to create a configuration file")
        raise typer.Exit(1)

    editor = get_editor()
    info(f"Opening {config_file} in {editor}...")

    try:
        subprocess.run([editor, str(config_file)], check=True)
    except FileNotFoundError:
        error(f"Editor '{editor}' not found")
        info("Set $VISUAL or $EDITOR to your preferred editor")
        raise typer.Exit(1)
    except subprocess.CalledProcessError as e:
        error(f"Editor exited with error: {e.returncode}")
        raise typer.Exit(e.returncode)
