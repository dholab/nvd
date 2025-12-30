"""
Typer application instance for the NVD CLI.

This module defines the main Typer app and shared configuration.
Commands are registered via the commands subpackage.
"""

import typer

# The main Typer application instance
app = typer.Typer(
    name="nvd",
    help="NVD: Novel Virus Detection pipeline CLI.",
    rich_markup_mode="rich",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
