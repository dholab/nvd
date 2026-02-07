"""
NVD CLI - A Typer-based command-line interface for the NVD pipeline.

This package provides a user-friendly wrapper around the Nextflow-based
NVD pipeline for metagenomic virus detection and taxonomic classification.

Usage:
    nvd run --samplesheet samples.csv
    nvd run --samplesheet samples.csv --experiment-id exp001 --labkey
    nvd version
    nvd --help
"""

from py_nvd.cli.app import app, main
from py_nvd.cli.commands.config import config_app
from py_nvd.cli.commands.params import params_app
from py_nvd.cli.commands.preset import preset_app
from py_nvd.cli.commands.validate import validate_app
from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    MAX_PREVIEW_ITEMS,
    PANEL_ANALYSIS,
    PANEL_CORE,
    PANEL_DATABASES,
    PANEL_LABKEY,
    PANEL_PREPROCESSING,
    PANEL_SRA,
    VALID_TOOLS,
    console,
    error,
    info,
    success,
    warning,
)

__all__ = [
    "DEFAULT_CONFIG",
    "MAX_PREVIEW_ITEMS",
    "PANEL_ANALYSIS",
    "PANEL_CORE",
    "PANEL_DATABASES",
    "PANEL_LABKEY",
    "PANEL_PREPROCESSING",
    "PANEL_SRA",
    "VALID_TOOLS",
    "app",
    "config_app",
    "console",
    "error",
    "info",
    "main",
    "params_app",
    "preset_app",
    "success",
    "validate_app",
    "warning",
]
