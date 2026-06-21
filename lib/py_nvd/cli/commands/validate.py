"""
Validation commands for the NVD CLI.

Commands:
    nvd validate deps       - Validate dependencies (Java, Nextflow, Docker/Apptainer)
    nvd validate samplesheet - Validate samplesheet format and content
    nvd validate databases  - Validate reference paths from config file
    nvd validate all        - Run all validation checks
    nvd validate params     - Validate a params file (alias for 'nvd params check')
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import typer

from py_nvd.cli.commands.params import params_check
from py_nvd.cli.samplesheet_validation import render_samplesheet_validation_result
from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    check_command_exists,
    console,
    docker_is_running,
    error,
    find_config_file,
    info,
    parse_nextflow_config,
    success,
    warning,
)
from py_nvd.samplesheet_validation import validate_samplesheet

validate_app = typer.Typer(
    name="validate",
    help="Validation commands to check system setup",
    no_args_is_help=True,
)


@validate_app.command("deps")
@validate_app.command("d", hidden=True)  # Alias
def validate_deps() -> None:
    """Validate dependencies (Java, Nextflow, Docker/Apptainer)."""
    console.print("\n[bold]Checking Dependencies[/bold]\n")

    all_good = True

    # Check Java
    if check_command_exists("java"):
        try:
            result = subprocess.run(
                ["java", "-version"],  # noqa: S607
                capture_output=True,
                text=True,
                timeout=3,
                check=False,
            )
            # Parse Java version
            version_line = result.stderr.strip().split("\n")[0]
            # Extract version number (handles both old 1.8.0 and new 17.0.0 formats)
            version_match = re.search(r'version "(\d+)\.?(\d+)?', version_line)
            if version_match:
                major = int(version_match.group(1))
                if major == 1 and version_match.group(2):
                    major = int(version_match.group(2))  # Old format: 1.8 -> 8

                if major >= 11:  # noqa: PLR2004
                    success(f"Java {version_line}")
                else:
                    warning(f"Java {version_line} (need Java 11+)")
                    all_good = False
            else:
                success("Java installed (version unknown)")
        except (OSError, subprocess.TimeoutExpired):
            warning("Java found but couldn't determine version")
    else:
        error("Java not found (required, need Java 11+)", exit_code=0)
        all_good = False

    # Check Nextflow
    if check_command_exists("nextflow"):
        try:
            result = subprocess.run(
                ["nextflow", "-version"],  # noqa: S607
                capture_output=True,
                text=True,
                timeout=3,
                check=False,
            )
            version_line = result.stdout.strip().split("\n")[0]
            success(f"Nextflow {version_line}")
        except (OSError, subprocess.TimeoutExpired):
            success("Nextflow installed")
    else:
        error("Nextflow not found (required)", exit_code=0)
        all_good = False

    # Check execution environment
    has_runtime = False
    if check_command_exists("docker"):
        if docker_is_running():
            success("Docker (running)")
            has_runtime = True
        else:
            warning("Docker installed but daemon not running")
    if check_command_exists("apptainer"):
        success("Apptainer")
        has_runtime = True
    elif check_command_exists("singularity"):
        success("Singularity")
        has_runtime = True

    if not has_runtime:
        warning("No container runtime found (Docker/Apptainer)")
        info("Pipeline can run in 'local' mode but containerization is recommended")

    console.print()
    if all_good:
        success("All critical dependencies are installed")
    else:
        console.print("[yellow]⚠ Some dependencies are missing or outdated[/yellow]")
        sys.exit(1)


@validate_app.command("samplesheet")
@validate_app.command("s", hidden=True)  # Alias
def validate_samplesheet_cmd(
    samplesheet: Path = typer.Argument(..., help="Path to samplesheet CSV"),
) -> None:
    """Validate samplesheet format and content."""

    console.print(f"\n[bold]Validating Samplesheet[/bold]: {samplesheet}\n")
    render_samplesheet_validation_result(validate_samplesheet(samplesheet))


@validate_app.command("databases")
@validate_app.command("db", hidden=True)  # Alias
def validate_databases(
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Config file to validate (default: NVD_CONFIG_DIR/user.config or NVD_CONFIG)",
    ),
) -> None:
    """Validate reference paths from config file."""
    console.print("\n[bold]Validating Reference Paths[/bold]\n")

    # Find config
    config_file = find_config_file(config)
    if config_file is None:
        warning(f"No config file found at {DEFAULT_CONFIG}")
        info("Run install.sh to set up reference configuration")
        sys.exit(1)

    assert config_file is not None  # Narrowing hint: guarded by sys.exit above
    info(f"Using config: {config_file}")

    # Parse config
    params = parse_nextflow_config(config_file)

    if not params:
        warning("No parameters found in config file")
        sys.exit(1)

    # Reference parameters to check
    reference_params = {
        "blast_db": "BLAST database directory",
    }

    console.print()
    all_found = True

    for param_name, description in reference_params.items():
        value = params.get(param_name)

        if not value or value == "null":
            warning(f"{description:30s} - not configured")
            all_found = False
            continue

        # Expand path
        db_path = Path(value).expanduser()

        if db_path.exists():
            success(f"{description:30s} - {db_path}")
        else:
            error(
                f"{description:30s} - NOT FOUND: {db_path}",
                exit_code=0,
            )
            all_found = False

    console.print()
    if all_found:
        success("All reference paths are valid")
    else:
        console.print("[yellow]⚠ Some reference paths are missing or invalid[/yellow]")
        info("Run install.sh to download and configure references")
        sys.exit(1)


@validate_app.command("all")
@validate_app.command("a", hidden=True)  # Alias
def validate_all(
    samplesheet: Path | None = typer.Option(
        None,
        "--samplesheet",
        "-s",
        help="Optionally validate a samplesheet",
    ),
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Config file to validate",
    ),
) -> None:
    """Run all validation checks."""
    console.print("\n[bold cyan]Running All Validation Checks[/bold cyan]")

    # Track failures
    failed = []

    # 1. Dependencies
    console.print("\n" + "=" * 70)
    try:
        validate_deps()
    except SystemExit as e:
        if e.code != 0:
            failed.append("Dependencies")

    # 2. References
    console.print("\n" + "=" * 70)
    try:
        validate_databases(config)
    except SystemExit as e:
        if e.code != 0:
            failed.append("References")

    # 3. Samplesheet (if provided)
    if samplesheet:
        console.print("\n" + "=" * 70)
        try:
            validate_samplesheet_cmd(samplesheet)
        except SystemExit as e:
            if e.code != 0:
                failed.append("Samplesheet")

    # Summary
    console.print("\n" + "=" * 70)
    console.print("\n[bold]Validation Summary[/bold]\n")

    if failed:
        console.print(f"[red]✗ Failed checks:[/red] {', '.join(failed)}")
        console.print()
        sys.exit(1)
    else:
        success("All validation checks passed")
        console.print()


@validate_app.command("params")
def validate_params(
    path: Path = typer.Argument(..., help="Path to params file to validate"),
    preset: str | None = typer.Option(
        None,
        "--preset",
        "-p",
        help="Merge with a preset before validating",
    ),
    show_defaults: bool = typer.Option(
        False,
        "--show-defaults",
        "-d",
        help="Show all parameters including defaults",
    ),
) -> None:
    """
    Validate a params file.

    This is an alias for 'nvd params check'. Validates the params file
    and displays the merged result with source tracking.

    Examples:

        nvd validate params my-params.yaml
        nvd validate params my-params.yaml --preset production
    """
    params_check(path=path, preset=preset, show_defaults=show_defaults)
