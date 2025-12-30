"""
NVD2 Command Line Interface

A user-friendly wrapper around the Nextflow-based NVD2 pipeline for
metagenomic virus detection and taxonomic classification.

Usage:
    nvd run --samplesheet samples.csv --experiment-id exp001
    nvd version

This module defines the Typer app and registers all commands.
Command implementations are in py_nvd.cli.commands.* modules.
"""

from __future__ import annotations

import sys
from datetime import datetime

import typer

from py_nvd.cli.commands.config import config_app
from py_nvd.cli.commands.params import params_app
from py_nvd.cli.commands.preset import preset_app
from py_nvd.cli.commands.resume import resume
from py_nvd.cli.commands.run import run
from py_nvd.cli.commands.samplesheet import samplesheet_app
from py_nvd.cli.commands.secrets import secrets_app
from py_nvd.cli.commands.state import state_app
from py_nvd.cli.commands.validate import validate_app
from py_nvd.cli.commands.version import version
from py_nvd.cli.utils import console

# ============================================================================
# TYPER APP
# ============================================================================

app = typer.Typer(
    name="nvd",
    help="NVD2 - Metagenomic virus detection pipeline",
    add_completion=True,
    no_args_is_help=True,
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
)

# ============================================================================
# COMMAND REGISTRATION
# ============================================================================

# Run command
app.command("run")(run)
app.command("r", hidden=True)(run)  # Alias

# Resume command
app.command("resume")(resume)

# Version command
app.command("version")(version)
app.command("v", hidden=True)(version)  # Alias

# Validation commands
app.add_typer(validate_app, name="validate")
app.add_typer(validate_app, name="val", hidden=True)  # Alias

# Params commands
app.add_typer(params_app, name="params")

# Preset commands
app.add_typer(preset_app, name="preset")

# Config commands
app.add_typer(config_app, name="config")
app.add_typer(config_app, name="cfg", hidden=True)  # Alias

# State commands
app.add_typer(state_app, name="state")

# Samplesheet commands
app.add_typer(samplesheet_app, name="samplesheet")
app.add_typer(samplesheet_app, name="sheet", hidden=True)  # Alias
app.add_typer(samplesheet_app, name="ss", hidden=True)  # Short alias

# Secrets commands
app.add_typer(secrets_app, name="secrets")

# Wrapped command (seasonal easter egg - Dec 15 to Jan 15)
_today = datetime.now()
_is_wrapped_season = (_today.month == 12 and _today.day >= 15) or (
    _today.month == 1 and _today.day <= 15
)

if _is_wrapped_season:
    from py_nvd.cli.commands.wrapped import wrapped

    app.command("wrapped")(wrapped)


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================


def main() -> None:
    """Main entry point for the CLI."""
    try:
        app()
    except KeyboardInterrupt:
        console.print("\n[yellow]Interrupted by user[/yellow]")
        sys.exit(130)
    except (typer.Exit, typer.Abort):
        # Normal exit from Typer
        raise
    except OSError as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
