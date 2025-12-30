# ruff: noqa: B008
"""
Parameter file management commands for the NVD CLI.

Commands:
    nvd params init     - Generate a template params file with schema reference
    nvd params validate - Validate a params file against the NVD schema
"""

from __future__ import annotations

from pathlib import Path

import typer

from py_nvd.cli.utils import (
    console,
    error,
    info,
    success,
    warning,
)

params_app = typer.Typer(
    name="params",
    help="Parameter file management (schema validation, templates)",
    no_args_is_help=True,
)


@params_app.command("init")
def params_init(
    output: Path = typer.Argument(
        ...,
        help="Output file path (.yaml or .json)",
    ),
    format: str | None = typer.Option(
        None,
        "--format",
        "-f",
        help="Output format: yaml or json (auto-detected from extension if not specified)",
    ),
) -> None:
    """
    Generate a template params file with schema reference.

    The generated file includes common parameters with defaults and
    comments explaining each option. Edit the file in your IDE to get
    autocomplete and validation.

    Examples:

        # Generate YAML template (default)
        nvd params init my-params.yaml

        # Generate JSON template
        nvd params init my-params.json

        # Specify format explicitly
        nvd params init my-params --format yaml
    """
    from py_nvd import params

    # Determine format from extension or option
    if format is None:
        if output.suffix in (".json",):
            format = "json"
        else:
            format = "yaml"  # Default to YAML

    # Ensure correct extension
    if format == "yaml" and output.suffix not in (".yaml", ".yml"):
        output = output.with_suffix(".yaml")
    elif format == "json" and output.suffix != ".json":
        output = output.with_suffix(".json")

    # Check if file already exists
    if output.exists():
        warning(f"File already exists: {output}")
        if not typer.confirm("Overwrite?"):
            raise typer.Abort()

    try:
        params.generate_template(output, format=format)
        success(f"Created {output}")
        info("Edit this file in your IDE for autocomplete and validation")
        info(f"Use with: nextflow run dhoconno/nvd -params-file {output}")
    except Exception as e:
        error(f"Failed to generate template: {e}")


@params_app.command("validate")
@params_app.command("val", hidden=True)  # Alias
def params_validate_cmd(
    path: Path = typer.Argument(
        ...,
        help="Params file to validate (.yaml or .json)",
        exists=True,
    ),
) -> None:
    """
    Validate a params file against the NVD schema.

    Checks that all parameter names are valid and values match expected types.
    Returns exit code 0 if valid, 1 if errors found.

    Examples:

        nvd params validate my-params.yaml
        nvd params validate config.json
    """
    from py_nvd import params

    console.print(f"\n[bold]Validating:[/bold] {path}\n")

    try:
        errors = params.validate_params_file(path)
    except FileNotFoundError:
        error(f"File not found: {path}")
    except Exception as e:
        error(f"Failed to parse file: {e}")

    if errors:
        console.print(f"[red]✗ {path} has {len(errors)} error(s):[/red]\n")
        for err in errors:
            console.print(f"  • {err}")
        console.print()
        raise typer.Exit(1)
    else:
        success(f"{path} is valid")
        console.print()
