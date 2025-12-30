# ruff: noqa: B008
"""
Preset management commands for the NVD CLI.

Commands:
    nvd preset list     - List all registered presets
    nvd preset show     - Show details of a specific preset
    nvd preset register - Register a new preset or update existing
    nvd preset export   - Export a preset to YAML or JSON
    nvd preset import   - Import a preset from a params file
    nvd preset delete   - Delete a preset
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.table import Table

from py_nvd.cli.utils import (
    console,
    error,
    info,
    success,
)

preset_app = typer.Typer(
    name="preset",
    help="Manage parameter presets (named collections of settings)",
    no_args_is_help=True,
)


@preset_app.command("list")
@preset_app.command("ls", hidden=True)  # Alias
def preset_list() -> None:
    """
    List all registered presets.

    Shows preset names, descriptions, parameter counts, and creation dates.
    """
    from py_nvd import state

    presets = state.list_presets()

    if not presets:
        info("No presets registered yet.")
        info("Create one with: nvd preset register <name> --from-file params.yaml")
        return

    console.print(f"\n[bold]Registered Presets ({len(presets)})[/bold]\n")

    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("Name", style="cyan")
    table.add_column("Description")
    table.add_column("Params", justify="right")
    table.add_column("Created", style="dim")

    for preset in presets:
        desc = preset.description or "[dim]—[/dim]"
        # Format date nicely (just date part)
        created = (
            preset.created_at.split("T")[0]
            if "T" in preset.created_at
            else preset.created_at.split(" ")[0]
        )
        table.add_row(preset.name, desc, str(len(preset.params)), created)

    console.print(table)
    console.print()


@preset_app.command("show")
def preset_show(
    name: str = typer.Argument(..., help="Name of the preset to show"),
) -> None:
    """
    Show details of a specific preset.

    Displays all parameters and their values.
    """
    from py_nvd import state

    preset = state.get_preset(name)

    if not preset:
        error(f"Preset not found: {name}")

    console.print(f"\n[bold cyan]Preset: {preset.name}[/bold cyan]\n")

    if preset.description:
        console.print(f"[dim]Description:[/dim] {preset.description}")

    console.print(f"[dim]Created:[/dim]     {preset.created_at}")
    console.print(f"[dim]Updated:[/dim]     {preset.updated_at}")
    console.print()

    console.print(f"[bold]Parameters ({len(preset.params)}):[/bold]")
    for key, value in sorted(preset.params.items()):
        console.print(f"  {key}: [green]{value}[/green]")

    console.print()


@preset_app.command("register")
@preset_app.command("reg", hidden=True)  # Alias
def preset_register(  # noqa: PLR0913
    name: str = typer.Argument(..., help="Name for the preset"),
    description: str | None = typer.Option(
        None,
        "--description",
        "-d",
        help="Description of the preset's purpose",
    ),
    from_file: Path | None = typer.Option(
        None,
        "--from-file",
        "-f",
        help="Import parameters from a YAML or JSON file",
        exists=True,
    ),
    # Common inline parameter options
    cutoff_percent: float | None = typer.Option(
        None, "--cutoff-percent", help="Minimum abundance threshold"
    ),
    tax_stringency: float | None = typer.Option(
        None, "--tax-stringency", help="Taxonomy stringency"
    ),
    entropy: float | None = typer.Option(None, "--entropy", help="Entropy threshold"),
    preprocess: bool | None = typer.Option(
        None, "--preprocess/--no-preprocess", help="Enable preprocessing"
    ),
    dedup: bool | None = typer.Option(
        None, "--dedup/--no-dedup", help="Deduplicate reads"
    ),
    trim_adapters: bool | None = typer.Option(
        None, "--trim-adapters/--no-trim-adapters", help="Trim adapters"
    ),
    scrub_host_reads: bool | None = typer.Option(
        None, "--scrub-host-reads/--no-scrub-host-reads", help="Scrub host reads"
    ),
    filter_reads: bool | None = typer.Option(
        None, "--filter-reads/--no-filter-reads", help="Filter reads"
    ),
    min_gottcha_reads: int | None = typer.Option(
        None, "--min-gottcha-reads", help="Minimum GOTTCHA2 reads"
    ),
    max_blast_targets: int | None = typer.Option(
        None, "--max-blast-targets", help="Maximum BLAST targets"
    ),
    blast_retention_count: int | None = typer.Option(
        None, "--blast-retention-count", help="BLAST hits to retain"
    ),
) -> None:
    """
    Register a new preset or update an existing one.

    Parameters can be imported from a YAML/JSON file or specified inline.
    Inline parameters override file parameters if both are provided.
    All parameters are validated against the NVD schema.

    Examples:

        # Register from a params file
        nvd preset register production --from-file prod.yaml -d "Production settings"

        # Register with inline parameters
        nvd preset register quick-test --cutoff-percent 0.01 --preprocess

        # Update existing preset
        nvd preset register production --from-file updated.yaml
    """
    from py_nvd import params, state

    preset_params: dict[str, str | int | float | bool | None] = {}

    # Load from file if specified
    if from_file:
        preset_params = params.load_params_file(from_file)
        info(f"Loaded {len(preset_params)} parameters from {from_file}")

    # Override/add inline params
    inline_params = {
        "cutoff_percent": cutoff_percent,
        "tax_stringency": tax_stringency,
        "entropy": entropy,
        "preprocess": preprocess,
        "dedup": dedup,
        "trim_adapters": trim_adapters,
        "scrub_host_reads": scrub_host_reads,
        "filter_reads": filter_reads,
        "min_gottcha_reads": min_gottcha_reads,
        "max_blast_targets": max_blast_targets,
        "blast_retention_count": blast_retention_count,
    }
    for key, value in inline_params.items():
        if value is not None:
            preset_params[key] = value

    if not preset_params:
        error("No parameters specified. Use --from-file or inline options.")

    # Validate against schema
    errors = params.validate_params(preset_params)
    if errors:
        console.print("[red]Invalid parameters:[/red]")
        for err in errors:
            console.print(f"  • {err}")
        raise typer.Exit(1)

    # Check if updating existing
    existing = state.get_preset(name)

    preset = state.register_preset(name, preset_params, description)

    if existing:
        success(f"Updated preset '{preset.name}' ({len(preset.params)} parameters)")
    else:
        success(f"Registered preset '{preset.name}' ({len(preset.params)} parameters)")


@preset_app.command("export")
def preset_export(
    name: str = typer.Argument(..., help="Name of the preset to export"),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output file path",
    ),
    format: str | None = typer.Option(
        None,
        "--format",
        "-f",
        help="Output format: yaml or json (auto-detected from extension)",
    ),
) -> None:
    """
    Export a preset to a YAML or JSON file.

    The exported file includes a schema reference for IDE support and
    can be used directly with Nextflow's -params-file option.

    Examples:

        nvd preset export production -o production.yaml
        nvd preset export production -o production.json --format json
    """
    import json as json_module
    from datetime import datetime

    import yaml

    from py_nvd import params, state

    preset = state.get_preset(name)
    if not preset:
        error(f"Preset not found: {name}")

    # Determine format
    if format is None:
        if output.suffix == ".json":
            format = "json"
        else:
            format = "yaml"

    schema_url = params.get_schema_url()

    if format == "yaml":
        lines = [
            f"# yaml-language-server: $schema={schema_url}",
            "#",
            f"# NVD Pipeline Preset: {preset.name}",
            f"# Exported: {datetime.now().isoformat()}",
            "#",
            f"# Usage: nextflow run dhoconno/nvd -params-file {output.name}",
        ]
        if preset.description:
            lines.append(f"# {preset.description}")
        lines.append("")

        with open(output, "w") as f:
            f.write("\n".join(lines) + "\n")
            yaml.dump(dict(preset.params), f, default_flow_style=False, sort_keys=True)
    else:
        export_data = {
            "$schema": schema_url,
            **preset.params,
        }
        with open(output, "w") as f:
            json_module.dump(export_data, f, indent=2)
            f.write("\n")

    success(f"Exported preset '{name}' to {output}")


@preset_app.command("import")
def preset_import(
    path: Path = typer.Argument(
        ...,
        help="Path to params file to import",
        exists=True,
    ),
    name: str | None = typer.Option(
        None,
        "--name",
        "-n",
        help="Name for the preset (default: filename without extension)",
    ),
    description: str | None = typer.Option(
        None,
        "--description",
        "-d",
        help="Description of the preset",
    ),
) -> None:
    """
    Import a preset from a YAML or JSON params file.

    If no name is specified, the filename (without extension) is used.

    Examples:

        nvd preset import shared-settings.yaml
        nvd preset import settings.yaml --name production -d "From shared config"
    """
    from py_nvd import params, state

    # Determine name from filename if not specified
    if name is None:
        name = path.stem

    # Load and validate
    preset_params = params.load_params_file(path)

    if not preset_params:
        error(f"No parameters found in {path}")

    errors = params.validate_params(preset_params)
    if errors:
        console.print(f"[red]Invalid parameters in {path}:[/red]")
        for err in errors:
            console.print(f"  • {err}")
        raise typer.Exit(1)

    # Check if updating existing
    existing = state.get_preset(name)

    preset = state.register_preset(name, preset_params, description)

    if existing:
        success(
            f"Updated preset '{preset.name}' from {path} ({len(preset.params)} parameters)"
        )
    else:
        success(
            f"Imported preset '{preset.name}' from {path} ({len(preset.params)} parameters)"
        )


@preset_app.command("delete")
@preset_app.command("rm", hidden=True)  # Alias
def preset_delete(
    name: str = typer.Argument(..., help="Name of the preset to delete"),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Delete without confirmation",
    ),
) -> None:
    """
    Delete a preset.

    Examples:

        nvd preset delete old-preset
        nvd preset rm old-preset --force
    """
    from py_nvd import state

    # Check if exists
    preset = state.get_preset(name)
    if not preset:
        error(f"Preset not found: {name}")

    # Confirm unless forced
    if not force:
        console.print(f"Delete preset '{name}' ({len(preset.params)} parameters)?")
        if not typer.confirm("Continue?"):
            raise typer.Abort()

    deleted = state.delete_preset(name)
    if deleted:
        success(f"Deleted preset '{name}'")
    else:
        error(f"Failed to delete preset '{name}'")
