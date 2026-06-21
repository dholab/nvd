"""CLI presentation helpers for samplesheet validation results."""

from __future__ import annotations

from py_nvd.cli.utils import console, error, success, warning
from py_nvd.samplesheet_validation import REQUIRED_COLUMNS, SamplesheetValidationResult


def render_samplesheet_validation_result(result: SamplesheetValidationResult) -> None:
    """Render a samplesheet validation result and exit non-zero on errors."""
    if result.header_valid:
        columns = result.header_columns or list(REQUIRED_COLUMNS)
        success(f"Header valid: {', '.join(columns)}")

    if result.errors:
        console.print("\n[red]Validation errors:[/red]")
        for err in result.errors:
            console.print(f"  • {err}")
        console.print()
        error(f"Found {len(result.errors)} validation error(s)")

    for warn in result.warnings:
        warning(warn)

    success(f"Found {len(result.samples)} valid samples")
    console.print()
    success("Samplesheet is valid")
