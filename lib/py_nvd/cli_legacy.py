"""
NVD2 Command Line Interface

A user-friendly wrapper around the Nextflow-based NVD2 pipeline for
metagenomic virus detection and taxonomic classification.

Usage:
    nvd run --samplesheet samples.csv --experiment-id exp001
    nvd version

This module serves as the orchestration layer for the CLI. All command
implementations have been extracted to py_nvd.cli.commands.* modules.
"""

from __future__ import annotations

import sys

import typer

from py_nvd.cli.utils import console

# ============================================================================
# TYPER APP & COMMANDS
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
# RUN COMMAND (extracted to cli/commands/run.py)
# ============================================================================

from py_nvd.cli.commands.run import run

app.command("run")(run)
app.command("r", hidden=True)(run)  # Alias


# ============================================================================
# VERSION COMMAND (extracted to cli/commands/version.py)
# ============================================================================

from py_nvd.cli.commands.version import version

app.command("version")(version)
app.command("v", hidden=True)(version)  # Alias


# ============================================================================
# VALIDATION COMMANDS (extracted to cli/commands/validate.py)
# ============================================================================

from py_nvd.cli.commands.validate import validate_app

app.add_typer(validate_app, name="validate")
app.add_typer(validate_app, name="val", hidden=True)  # Alias


# ============================================================================
# PARAMS COMMANDS (extracted to cli/commands/params.py)
# ============================================================================

from py_nvd.cli.commands.params import params_app

app.add_typer(params_app, name="params")


# ============================================================================
# PRESET COMMANDS (extracted to cli/commands/preset.py)
# ============================================================================

from py_nvd.cli.commands.preset import preset_app

app.add_typer(preset_app, name="preset")


# ============================================================================
# CONFIG COMMANDS (extracted to cli/commands/config.py)
# ============================================================================

from py_nvd.cli.commands.config import config_app

app.add_typer(config_app, name="config")
app.add_typer(config_app, name="cfg", hidden=True)  # Alias


# ============================================================================
# STATE COMMANDS (cli/commands/state.py)
# ============================================================================

from py_nvd.cli.commands.state import state_app

app.add_typer(state_app, name="state")


# ============================================================================
# WRAPPED COMMAND (seasonal easter egg - Dec 15 to Jan 15)
# ============================================================================

from datetime import datetime

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
