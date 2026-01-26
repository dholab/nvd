# ruff: noqa: B008
"""
Preset management commands for the NVD CLI.

Commands:
    nvd preset list     - List all registered presets
    nvd preset show     - Show details of a specific preset
    nvd preset register - Register a new preset or update existing
    nvd preset merge    - Merge multiple presets into a new preset
    nvd preset diff     - Show differences between two presets
    nvd preset export   - Export a preset to YAML or JSON
    nvd preset import   - Import a preset from a params file
    nvd preset delete   - Delete a preset
"""

from __future__ import annotations

from pathlib import Path

import typer
from pydantic import ValidationError
from rich.panel import Panel
from rich.table import Table

from py_nvd import params, state
from py_nvd.cli.utils import (
    console,
    error,
    info,
    success,
)
from py_nvd.models import (
    PARAM_CATEGORIES,
    NvdParams,
    TracedParams,
    get_field_category,
    trace_merge,
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

    # Validate with Pydantic model
    try:
        NvdParams.model_validate(preset_params)
    except ValidationError as e:
        console.print("[red]Invalid parameters:[/red]")
        for err in e.errors():
            field = ".".join(str(loc) for loc in err["loc"])
            console.print(f"  • [bold]{field}[/bold]: {err['msg']}")
        raise typer.Exit(1) from None

    # Check if updating existing
    existing = state.get_preset(name)

    preset = state.register_preset(name, preset_params, description)

    if existing:
        success(f"Updated preset '{preset.name}' ({len(preset.params)} parameters)")
    else:
        success(f"Registered preset '{preset.name}' ({len(preset.params)} parameters)")


@preset_app.command("merge")
def preset_merge(
    presets: list[str] = typer.Argument(
        ...,
        help="Presets to merge (later presets take precedence)",
    ),
    name: str = typer.Option(
        ...,
        "--name",
        "-n",
        help="Name for the new merged preset",
    ),
    description: str | None = typer.Option(
        None,
        "--description",
        "-d",
        help="Description for the new preset",
    ),
) -> None:
    """
    Merge multiple presets into a new preset.

    Presets are merged in order, with later presets taking precedence over
    earlier ones. The merged result is validated before being registered.

    Examples:

        # Merge two presets (production values override base)
        nvd preset merge base production --name prod-merged

        # Merge with description
        nvd preset merge base tweaks --name custom -d "Base with custom tweaks"

        # Merge three presets
        nvd preset merge defaults team-settings my-overrides --name final
    """
    if len(presets) < 2:
        error("At least two presets are required for merging")
        raise typer.Exit(1)

    # Load all presets with labels for trace_merge
    labeled_sources: list[tuple[str, dict]] = []
    for preset_name in presets:
        preset_obj = state.get_preset(preset_name)
        if not preset_obj:
            error(f"Preset not found: {preset_name}")
            available = state.list_presets()
            if available:
                console.print("\n[cyan]Available presets:[/cyan]")
                for p in available:
                    console.print(f"  • {p.name}")
            raise typer.Exit(1)
        labeled_sources.append((preset_name, dict(preset_obj.params)))

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

    # Check if name already exists
    existing = state.get_preset(name)
    if existing:
        console.print(f"[yellow]Preset '{name}' already exists.[/yellow]")
        if not typer.confirm("Overwrite?"):
            raise typer.Abort()

    # Build merged params dict (only explicitly set values, not defaults)
    explicitly_set = set()
    for _, params_dict in labeled_sources:
        explicitly_set.update(params_dict.keys())

    merged_params = {
        k: v
        for k, v in traced.params.model_dump(exclude_none=True).items()
        if k in explicitly_set
    }

    # Register the new preset
    new_preset = state.register_preset(name, merged_params, description)

    # Display merge summary
    _display_preset_merge_summary(presets, traced, new_preset.name, explicitly_set)


def _display_preset_merge_summary(
    preset_names: list[str],
    traced: TracedParams,
    new_name: str,
    explicitly_set: set[str],
) -> None:
    """Display a Rich Panel summarizing the preset merge result."""
    # Build title from preset names
    names_str = " + ".join(preset_names)
    title = f"Preset Merge: {names_str}"

    # Group sources by category (only explicitly set values)
    categories: dict[str, list] = {cat: [] for cat in PARAM_CATEGORIES}

    for src in traced.sources:
        if src.source == "default":
            continue
        if src.value is None:
            continue
        if src.field not in explicitly_set:
            continue

        category = get_field_category(src.field)
        categories[category].append(src)

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
            param_count += 1

            # Format value (truncate long values)
            value_str = str(src.value)
            if len(value_str) > 40:
                value_str = value_str[:37] + "..."

            # Format source with override info
            source_str = f"← {src.source}"
            if src.overridden_value is not None:
                ov = str(src.overridden_value)
                if len(ov) > 20:
                    ov = ov[:17] + "..."
                source_str += f" [dim](was {src.overridden_source}: {ov})[/dim]"

            table.add_row(f"  {src.field}", value_str, source_str)

        table.add_row("", "", "")  # Spacer between categories

    # Summary line
    summary = f"[green]✓ Registered as '{new_name}'[/green] ({param_count} params)"

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


@preset_app.command("export")
def preset_export(
    name: str = typer.Argument(..., help="Name of the preset to export"),
    output: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Output file path",
    ),
    output_format: str | None = typer.Option(
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

    import yaml

    from py_nvd.db import utc_now_iso

    preset = state.get_preset(name)
    if not preset:
        error(f"Preset not found: {name}")

    # Determine format
    if output_format is None:
        if output.suffix == ".json":
            output_format = "json"
        else:
            output_format = "yaml"

    schema_url = params.get_schema_url()

    if output_format == "yaml":
        lines = [
            f"# yaml-language-server: $schema={schema_url}",
            "#",
            f"# NVD Pipeline Preset: {preset.name}",
            f"# Exported: {utc_now_iso()}",
            "#",
            f"# Usage: nextflow run dhoconno/nvd -params-file {output.name}",
        ]
        if preset.description:
            lines.append(f"# {preset.description}")
        lines.append("")

        with open(output, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
            yaml.dump(dict(preset.params), f, default_flow_style=False, sort_keys=True)
    else:
        export_data = {
            "$schema": schema_url,
            **preset.params,
        }
        with open(output, "w", encoding="utf-8") as f:
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
    # Determine name from filename if not specified
    if name is None:
        name = path.stem

    # Load and validate
    preset_params = params.load_params_file(path)

    if not preset_params:
        error(f"No parameters found in {path}")

    try:
        NvdParams.model_validate(preset_params)
    except ValidationError as e:
        console.print(f"[red]Invalid parameters in {path}:[/red]")
        for err in e.errors():
            field = ".".join(str(loc) for loc in err["loc"])
            console.print(f"  • [bold]{field}[/bold]: {err['msg']}")
        raise typer.Exit(1) from None

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


@preset_app.command("diff")
def preset_diff(
    preset1: str = typer.Argument(..., help="First preset (base)"),
    preset2: str = typer.Argument(..., help="Second preset (to compare)"),
) -> None:
    """
    Show differences between two presets.

    Compares the two presets and shows which parameters differ, which are
    only in one preset, and which are the same. Output is grouped by category.

    Examples:

        nvd preset diff base production
        nvd preset diff defaults my-settings
    """
    # Load both presets
    p1 = state.get_preset(preset1)
    if not p1:
        error(f"Preset not found: {preset1}")
        raise typer.Exit(1)

    p2 = state.get_preset(preset2)
    if not p2:
        error(f"Preset not found: {preset2}")
        raise typer.Exit(1)

    params1 = dict(p1.params)
    params2 = dict(p2.params)

    # Find differences
    all_keys = set(params1.keys()) | set(params2.keys())

    # Categorize differences
    changed: list[tuple[str, str, object, object]] = []  # (category, field, val1, val2)
    only_in_first: list[tuple[str, str, object]] = []  # (category, field, val)
    only_in_second: list[tuple[str, str, object]] = []  # (category, field, val)

    for key in sorted(all_keys):
        category = get_field_category(key)
        in_first = key in params1
        in_second = key in params2

        if in_first and in_second:
            if params1[key] != params2[key]:
                changed.append((category, key, params1[key], params2[key]))
        elif in_first:
            only_in_first.append((category, key, params1[key]))
        else:
            only_in_second.append((category, key, params2[key]))

    # Check if identical
    if not changed and not only_in_first and not only_in_second:
        info(f"Presets '{preset1}' and '{preset2}' are identical")
        return

    # Build display
    title = f"Diff: {preset1} → {preset2}"

    table = Table(
        show_header=False,
        box=None,
        padding=(0, 2),
        collapse_padding=True,
    )
    table.add_column("Field", style="cyan", no_wrap=True)
    table.add_column("Change", style="white")

    # Group by category for display
    def add_section(
        items,  # list of tuples: (category, field, ...)
        format_fn,
    ):
        """Add items grouped by category."""
        by_category: dict[str, list] = {cat: [] for cat in PARAM_CATEGORIES}
        for item in items:
            by_category[item[0]].append(item)

        for category in PARAM_CATEGORIES:
            cat_items = by_category[category]
            if not cat_items:
                continue
            table.add_row(f"[bold]{category}[/bold]", "")
            for item in cat_items:
                field = item[1]
                change_str = format_fn(item)
                table.add_row(f"  {field}", change_str)
            table.add_row("", "")

    # Changed values
    if changed:
        table.add_row("[bold yellow]Changed[/bold yellow]", "")
        table.add_row("", "")
        add_section(
            changed,
            lambda x: f"[red]{_format_diff_value(x[2])}[/red] → [green]{_format_diff_value(x[3])}[/green]",
        )

    # Only in first
    if only_in_first:
        table.add_row(f"[bold red]Only in {preset1}[/bold red]", "")
        table.add_row("", "")
        add_section(
            only_in_first,
            lambda x: f"[red]{_format_diff_value(x[2])}[/red]",
        )

    # Only in second
    if only_in_second:
        table.add_row(f"[bold green]Only in {preset2}[/bold green]", "")
        table.add_row("", "")
        add_section(
            only_in_second,
            lambda x: f"[green]{_format_diff_value(x[2])}[/green]",
        )

    # Summary
    parts = []
    if changed:
        parts.append(f"{len(changed)} changed")
    if only_in_first:
        parts.append(f"{len(only_in_first)} only in {preset1}")
    if only_in_second:
        parts.append(f"{len(only_in_second)} only in {preset2}")
    summary = ", ".join(parts)

    panel = Panel(
        table,
        title=f"[bold]{title}[/bold]",
        subtitle=f"[dim]{summary}[/dim]",
        border_style="blue",
        padding=(1, 2),
    )

    console.print()
    console.print(panel)
    console.print()


def _format_diff_value(value: object) -> str:
    """Format a value for diff display, truncating if needed."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, list):
        if len(value) > 3:
            return f"[{', '.join(str(v) for v in value[:3])}, ...]"
        return f"[{', '.join(str(v) for v in value)}]"
    s = str(value)
    if len(s) > 30:
        return s[:27] + "..."
    return s
