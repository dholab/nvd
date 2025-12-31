# ruff: noqa: B008
"""
Parameter file management commands for the NVD CLI.

Commands:
    nvd params init  - Generate a template params file with schema reference
    nvd params check - Validate and display params with source tracking
    nvd params merge - Merge multiple params files into one
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.panel import Panel
from rich.table import Table

from py_nvd.cli.utils import (
    console,
    error,
    info,
    success,
    warning,
)
from py_nvd.models import (
    PARAM_CATEGORIES,
    ParamSource,
    TracedParams,
    get_field_category,
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


@params_app.command("check")
@params_app.command("validate", hidden=True)  # Backwards compat alias
def params_check(
    path: Path = typer.Argument(
        ...,
        help="Params file to check (.yaml or .json)",
        exists=True,
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        "-p",
        help="Merge with a preset to see combined result",
    ),
    show_defaults: bool = typer.Option(
        False,
        "--show-defaults",
        "-d",
        help="Show parameters using default values",
    ),
) -> None:
    """
    Validate a params file and show merged result with sources.

    Validates parameter types and values using the NvdParams model.
    When --preset is specified, shows the merged result with annotations
    indicating where each value came from.

    Precedence (highest to lowest):
        1. Params file values
        2. Preset values
        3. Pipeline defaults

    Examples:

        # Check a params file
        nvd params check my-params.yaml

        # Check merged with a preset
        nvd params check my-params.yaml --preset production

        # Show all params including defaults
        nvd params check my-params.yaml --show-defaults
    """
    from pydantic import ValidationError

    from py_nvd import state
    from py_nvd.models import ParamSource, trace_merge
    from py_nvd.params import load_params_file

    # Load params file
    try:
        file_params = load_params_file(path)
    except FileNotFoundError:
        error(f"File not found: {path}")
        raise typer.Exit(1) from None
    except Exception as e:
        error(f"Failed to parse file: {e}")
        raise typer.Exit(1) from None

    # Load preset if specified
    preset_params: dict = {}
    preset_name = None
    if preset:
        preset_obj = state.get_preset(preset)
        if not preset_obj:
            error(f"Preset not found: {preset}")
            available = state.list_presets()
            if available:
                console.print("\n[cyan]Available presets:[/cyan]")
                for p in available:
                    console.print(f"  • {p.name}")
            raise typer.Exit(1)
        preset_params = preset_obj.params
        preset_name = preset_obj.name

    # Merge and validate
    try:
        traced = trace_merge(
            ("preset", preset_params),
            ("file", file_params),
        )
    except ValidationError as e:
        console.print(f"\n[red]✗ Validation failed:[/red]\n")
        for err in e.errors():
            field = ".".join(str(loc) for loc in err["loc"])
            msg = err["msg"]
            console.print(f"  • [bold]{field}[/bold]: {msg}")
        console.print()
        raise typer.Exit(1) from None

    # Build display
    title = f"Params Check: {path.name}"
    if preset_name:
        title += f" + preset:{preset_name}"

    # Group sources by category for display
    categories = _categorize_sources(traced.sources, show_defaults)

    # Build the table
    table = Table(
        show_header=False,
        box=None,
        padding=(0, 2),
        collapse_padding=True,
    )
    table.add_column("Field", style="cyan", no_wrap=True)
    table.add_column("Value", style="white")
    table.add_column("Source", style="dim")

    non_default_count = 0
    for category, sources in categories.items():
        if not sources:
            continue

        # Category header
        table.add_row(f"[bold]{category}[/bold]", "", "")

        for src in sources:
            if src.source != "default":
                non_default_count += 1

            # Format value (truncate long values)
            value_str = _format_value(src.value)

            # Format source with override info
            source_str = _format_source(src)

            table.add_row(f"  {src.field}", value_str, source_str)

        table.add_row("", "", "")  # Spacer between categories

    # Summary line
    summary = f"[green]✓ Valid[/green] ({non_default_count} params set"
    if traced.defaults_used > 0:
        summary += f", {traced.defaults_used} defaults"
    summary += ")"

    # Wrap in panel
    panel = Panel(
        table,
        title=f"[bold]{title}[/bold]",
        subtitle=summary,
        border_style="blue",
        padding=(1, 2),
    )

    console.print()
    console.print(panel)
    console.print()


def _categorize_sources(
    sources: list[ParamSource],
    show_defaults: bool,
) -> dict[str, list[ParamSource]]:
    """
    Group sources by category for display.

    Categories are derived from NvdParams field metadata (json_schema_extra),
    ensuring a single source of truth for field categorization.
    """
    # Initialize categories in display order
    categories: dict[str, list[ParamSource]] = {cat: [] for cat in PARAM_CATEGORIES}

    for src in sources:
        # Skip defaults unless requested
        if src.source == "default" and not show_defaults:
            continue

        # Skip None values
        if src.value is None:
            continue

        category = get_field_category(src.field)
        categories[category].append(src)

    return categories


def _format_value(value) -> str:
    """Format a value for display, truncating if needed."""
    if isinstance(value, list):
        if len(value) > 3:
            return f"[{', '.join(str(v) for v in value[:3])}, ...]"
        return f"[{', '.join(str(v) for v in value)}]"

    str_value = str(value)
    if len(str_value) > 50:
        return str_value[:47] + "..."
    return str_value


def _format_source(src: ParamSource) -> str:
    """Format source with override info."""
    if src.source == "default":
        return "[dim]default[/dim]"

    result = f"← {src.source}"
    if src.overridden_value is not None:
        # Truncate overridden value if needed
        ov = str(src.overridden_value)
        if len(ov) > 20:
            ov = ov[:17] + "..."
        result += f" [dim](was {src.overridden_source}: {ov})[/dim]"

    return result


@params_app.command("merge")
def params_merge(
    files: list[Path] = typer.Argument(
        ...,
        help="Params files to merge (later files take precedence)",
        exists=True,
    ),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output file path (.yaml or .json)",
    ),
    format: str | None = typer.Option(
        None,
        "--format",
        "-f",
        help="Output format: yaml or json (auto-detected from extension if not specified)",
    ),
    include_defaults: bool = typer.Option(
        False,
        "--include-defaults",
        "-d",
        help="Include all default values in output (default: only merged values)",
    ),
) -> None:
    """
    Merge multiple params files into one.

    Files are merged in order, with later files taking precedence over earlier
    ones. The merged result is validated against the NvdParams model before
    being written.

    By default, only explicitly set values are included in the output. Use
    --include-defaults to include all parameters with their default values.

    Examples:

        # Merge two files (settings.yaml values override base.yaml)
        nvd params merge base.yaml settings.yaml -o merged.yaml

        # Merge multiple files with explicit format
        nvd params merge a.yaml b.json c.yaml -o result.json

        # Include all defaults in output
        nvd params merge base.yaml overrides.yaml -o full.yaml --include-defaults
    """
    import json as json_module

    import yaml
    from pydantic import ValidationError

    from py_nvd.models import trace_merge
    from py_nvd.params import get_schema_url, load_params_file

    if len(files) < 2:
        error("At least two files are required for merging")
        raise typer.Exit(1)

    # Load all files with labels for trace_merge
    labeled_sources: list[tuple[str, dict]] = []
    for path in files:
        try:
            params_dict = load_params_file(path)
            labeled_sources.append((path.name, params_dict))
        except FileNotFoundError:
            error(f"File not found: {path}")
            raise typer.Exit(1) from None
        except Exception as e:
            error(f"Failed to parse {path}: {e}")
            raise typer.Exit(1) from None

    # Merge with source tracking
    try:
        traced = trace_merge(*labeled_sources)
    except ValidationError as e:
        console.print("\n[red]✗ Validation failed after merge:[/red]\n")
        for err in e.errors():
            field = ".".join(str(loc) for loc in err["loc"])
            msg = err["msg"]
            console.print(f"  • [bold]{field}[/bold]: {msg}")
        console.print()
        raise typer.Exit(1) from None

    # Determine output format
    if format is None:
        if output.suffix == ".json":
            format = "json"
        else:
            format = "yaml"

    # Ensure correct extension
    if format == "yaml" and output.suffix not in (".yaml", ".yml"):
        output = output.with_suffix(".yaml")
    elif format == "json" and output.suffix != ".json":
        output = output.with_suffix(".json")

    # Check if output exists
    if output.exists():
        warning(f"File already exists: {output}")
        if not typer.confirm("Overwrite?"):
            raise typer.Abort()

    # Build output data - only include explicitly set values unless --include-defaults
    explicitly_set = set()
    for _, params_dict in labeled_sources:
        explicitly_set.update(params_dict.keys())

    if include_defaults:
        output_data = traced.params.model_dump(exclude_none=True)
    else:
        output_data = {
            k: v
            for k, v in traced.params.model_dump(exclude_none=True).items()
            if k in explicitly_set
        }

    # Convert Path objects to strings for serialization
    for key, value in output_data.items():
        if isinstance(value, Path):
            output_data[key] = str(value)

    # Write output
    schema_url = get_schema_url()

    if format == "yaml":
        lines = [
            f"# yaml-language-server: $schema={schema_url}",
            "#",
            f"# Merged from: {', '.join(f.name for f in files)}",
            "#",
        ]
        with open(output, "w") as f:
            f.write("\n".join(lines) + "\n\n")
            yaml.dump(output_data, f, default_flow_style=False, sort_keys=True)
    else:
        output_data_with_schema = {"$schema": schema_url, **output_data}
        with open(output, "w") as f:
            json_module.dump(output_data_with_schema, f, indent=2)
            f.write("\n")

    # Display merge summary with Rich Panel
    _display_merge_summary(files, traced, output, include_defaults)


def _display_merge_summary(
    files: list[Path],
    traced: TracedParams,
    output: Path,
    show_defaults: bool,
) -> None:
    """Display a Rich Panel summarizing the merge result."""
    # Build title from file names
    file_names = " + ".join(f.name for f in files)
    title = f"Params Merge: {file_names}"

    # Group sources by category
    categories = _categorize_sources(traced.sources, show_defaults)

    # Build the table
    table = Table(
        show_header=False,
        box=None,
        padding=(0, 2),
        collapse_padding=True,
    )
    table.add_column("Field", style="cyan", no_wrap=True)
    table.add_column("Value", style="white")
    table.add_column("Source", style="dim")

    param_count = 0
    for category, sources in categories.items():
        if not sources:
            continue

        # Category header
        table.add_row(f"[bold]{category}[/bold]", "", "")

        for src in sources:
            if src.source != "default":
                param_count += 1

            # Format value (truncate long values)
            value_str = _format_value(src.value)

            # Format source with override info
            source_str = _format_source(src)

            table.add_row(f"  {src.field}", value_str, source_str)

        table.add_row("", "", "")  # Spacer between categories

    # Summary line
    summary = f"[green]✓ Written to {output.name}[/green] ({param_count} params)"

    # Wrap in panel
    panel = Panel(
        table,
        title=f"[bold]{title}[/bold]",
        subtitle=summary,
        border_style="blue",
        padding=(1, 2),
    )

    console.print()
    console.print(panel)
    console.print()
