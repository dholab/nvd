"""
Command modules for the NVD CLI.

Each submodule defines one or more Typer commands that are registered
with the main app in py_nvd/cli/__init__.py.
"""

# Commands are imported here as they are extracted from cli_legacy.py.
# Importing them triggers registration with the main app.
from py_nvd.cli.commands import (
    config,
    params,
    preset,
    run,
    state,
    validate,
    version,
    wrapped,
)

__all__: list[str] = [
    "config",
    "params",
    "preset",
    "run",
    "state",
    "validate",
    "version",
    "wrapped",
]
