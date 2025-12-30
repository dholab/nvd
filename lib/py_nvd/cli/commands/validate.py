# ruff: noqa: B008
"""
Validation commands for the NVD CLI.

Commands:
    nvd validate deps       - Validate dependencies (Java, Nextflow, Docker/Apptainer)
    nvd validate samplesheet - Validate samplesheet format and content
    nvd validate databases  - Validate database paths from config file
    nvd validate all        - Run all validation checks
    nvd validate params     - Validate CLI coverage of Nextflow parameters
"""

from __future__ import annotations

import ast
import csv
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import typer

from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    MAX_PREVIEW_ITEMS,
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

# ============================================================================
# PARAMETER COVERAGE VALIDATION HELPERS
# ============================================================================


@dataclass
class ParamInfo:
    """Represents a Nextflow parameter."""

    name: str
    default_value: str | None
    category: str  # "configurable" or "internal"
    line_number: int


@dataclass
class CliOptionInfo:
    """Represents a CLI option."""

    name: str  # Python parameter name (e.g., "gottcha2_db")
    flags: list[str]  # CLI flags (e.g., ["--gottcha2-db"])
    command: str  # Command where it appears (e.g., "run")


# Known internal parameters that should NOT be exposed in CLI
INTERNAL_PARAMS = {
    "date",
    "resources",
    "refman_registry",
    "monoimage",
    "gottcha2_db_version",
    "blast_db_version",
    "nvd_files",
    "human_virus_families",  # Complex list type
}


def categorize_param(name: str, value: str, preceding_comment: str) -> str:  # noqa: PLR0911
    """
    Determine if a parameter should be CLI-exposed.

    Args:
        name: Parameter name
        value: Parameter value
        preceding_comment: Comments before this parameter

    Returns:
        "internal" or "configurable"
    """
    # Explicit internal list
    if name in INTERNAL_PARAMS:
        return "internal"

    # Version tracking
    if name.endswith("_version"):
        return "internal"

    # Complex types (lists)
    if value.strip().startswith("["):
        return "internal"

    # Java expressions
    if "new java" in value:
        return "internal"

    # Project/launch directory paths
    if "${projectDir}" in value or "${launchDir}" in value:
        return "internal"

    # Check for "SHOULD NOT BE CHANGED" in comments
    if "SHOULD NOT BE CHANGED" in preceding_comment:
        return "internal"

    return "configurable"


def extract_nextflow_params(config_path: Path) -> dict[str, ParamInfo]:
    """
    Parse nextflow.config and extract all params from params { } block.

    Returns:
        dict mapping parameter name to ParamInfo
    """
    params = {}

    try:
        with open(config_path) as f:
            lines = f.readlines()

        in_params_block = False
        brace_depth = 0
        preceding_comment = ""

        for line_num, line in enumerate(lines, start=1):
            stripped = line.strip()

            # Track comments
            if stripped.startswith("//"):
                preceding_comment += stripped + "\n"
                continue

            # Detect params block start
            if not in_params_block and re.match(r"params\s*\{", stripped):
                in_params_block = True
                brace_depth = 1
                continue

            if not in_params_block:
                continue

            # Track brace depth
            brace_depth += stripped.count("{") - stripped.count("}")

            # Exit params block
            if brace_depth == 0:
                break

            # Parse parameter line: name = value
            match = re.match(r"(\w+)\s*=\s*(.+)", stripped)
            if match:
                param_name = match.group(1)
                param_value = match.group(2).strip()

                # Remove trailing comments
                param_value = re.sub(r"\s*//.*$", "", param_value)

                category = categorize_param(param_name, param_value, preceding_comment)

                params[param_name] = ParamInfo(
                    name=param_name,
                    default_value=param_value,
                    category=category,
                    line_number=line_num,
                )

                # Reset comment accumulator
                preceding_comment = ""

    except OSError as e:
        warning(f"Failed to parse Nextflow config: {e}")
        return {}

    return params


def extract_cli_options() -> dict[str, CliOptionInfo]:  # noqa: C901, PLR0912
    """
    Extract all CLI options from the run() command.

    Uses AST parsing to analyze the run() function definition.

    Returns:
        dict mapping parameter name to CliOptionInfo
    """
    options = {}

    try:
        # Read run.py where run() is defined
        source_path = Path(__file__).parent / "run.py"
        with open(source_path) as f:
            source = f.read()

        # Parse AST
        tree = ast.parse(source)

        # Find the run() function
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            if node.name != "run":
                continue

            # Analyze function arguments
            for arg in node.args.args:
                param_name = arg.arg

                # Look for default value with typer.Option
                defaults_offset = len(node.args.args) - len(node.args.defaults)
                arg_index = node.args.args.index(arg)
                default_index = arg_index - defaults_offset

                if default_index < 0 or default_index >= len(node.args.defaults):
                    continue

                default = node.args.defaults[default_index]

                # Check if it's a typer.Option() call
                if not isinstance(default, ast.Call):
                    continue

                # Extract flags from Option call
                flags = []
                for kw in default.keywords:
                    if isinstance(kw.value, ast.Constant) and isinstance(
                        kw.value.value,
                        str,
                    ):
                        value = kw.value.value
                        if value.startswith(("--", "-")):
                            flags.append(value)

                # Also check positional args for flags
                for call_arg in default.args:
                    if isinstance(call_arg, ast.Constant) and isinstance(
                        call_arg.value,
                        str,
                    ):
                        value = call_arg.value
                        if value.startswith(("--", "-")):
                            flags.append(value)

                if flags or param_name:
                    options[param_name] = CliOptionInfo(
                        name=param_name,
                        flags=flags,
                        command="run",
                    )

            break  # Found run(), no need to continue

    except OSError as e:
        warning(f"Failed to parse CLI source: {e}")
        return {}

    return options


def compare_coverage(
    nf_params: dict[str, ParamInfo],
    cli_options: dict[str, CliOptionInfo],
) -> tuple[list[ParamInfo], list[str], list[ParamInfo]]:
    """
    Compare Nextflow params to CLI options.

    Returns:
        (missing, extra, internal)
        - missing: Configurable NF params not in CLI
        - extra: CLI options not in NF params (informational)
        - internal: Internal params correctly excluded
    """
    # Get configurable params
    configurable = {
        name: info
        for name, info in nf_params.items()
        if info.category == "configurable"
    }

    # Get internal params
    internal = [info for info in nf_params.values() if info.category == "internal"]

    # CLI option names (normalize underscores)
    cli_param_names = set(cli_options.keys())

    # Find missing
    missing = [
        info for name, info in configurable.items() if name not in cli_param_names
    ]

    # Find extra (CLI options not in Nextflow)
    extra = [name for name in cli_param_names if name not in nf_params]

    return missing, extra, internal


# ============================================================================
# VALIDATION COMMANDS
# ============================================================================

validate_app = typer.Typer(
    name="validate",
    help="Validation commands to check system setup",
    no_args_is_help=True,
)


@validate_app.command("deps")
@validate_app.command("d", hidden=True)  # Alias
def validate_deps() -> None:  # noqa: PLR0915, PLR0912, C901
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
def validate_samplesheet_cmd(  # noqa: PLR0912, C901
    samplesheet: Path = typer.Argument(..., help="Path to samplesheet CSV"),
) -> None:
    """Validate samplesheet format and content."""

    console.print(f"\n[bold]Validating Samplesheet[/bold]: {samplesheet}\n")

    if not samplesheet.exists():
        error(f"Samplesheet not found: {samplesheet}")

    # Read and validate
    required_columns = {"sample_id", "srr", "platform", "fastq1", "fastq2"}
    valid_platforms = {"illumina", "ont", ""}

    try:
        with samplesheet.open() as f:
            reader = csv.DictReader(f)

            if not reader.fieldnames:
                error("Samplesheet is empty or has no header")

            # Check for required columns
            fieldnames = set(reader.fieldnames) if reader.fieldnames else set()
            missing_cols = required_columns - fieldnames
            if missing_cols:
                error(f"Missing required columns: {', '.join(missing_cols)}")

            success(f"Header valid: {', '.join(reader.fieldnames or [])}")

            # Validate rows
            samples = []
            errors_found = []
            for i, row in enumerate(reader, start=2):  # Start at 2 (header is line 1)
                sample_id = row.get("sample_id", "").strip()

                # Skip comment lines
                if sample_id.startswith("#"):
                    continue

                # Check sample_id not empty
                if not sample_id:
                    errors_found.append(f"Line {i}: sample_id is empty")
                    continue

                # Check platform is valid
                platform = row.get("platform", "").strip().lower()
                if platform and platform not in valid_platforms:
                    errors_found.append(
                        f"Line {i} ({sample_id}): invalid platform '{platform}'",
                    )

                # Check that either SRR or FASTQ files are provided
                srr = row.get("srr", "").strip()
                fastq1 = row.get("fastq1", "").strip()

                if not srr and not fastq1:
                    msg = (
                        f"Line {i} ({sample_id}): must provide either 'srr' or 'fastq1'"
                    )
                    errors_found.append(msg)

                samples.append(sample_id)

            if errors_found:
                console.print("\n[red]Validation errors:[/red]")
                for err in errors_found:
                    console.print(f"  • {err}")
                console.print()
                error(f"Found {len(errors_found)} validation error(s)")

            # Check for duplicate sample IDs
            duplicates = [s for s in set(samples) if samples.count(s) > 1]
            if duplicates:
                warning(f"Duplicate sample_id found: {', '.join(duplicates)}")

            success(f"Found {len(samples)} valid samples")

            console.print()
            success("Samplesheet is valid")

    except csv.Error as e:
        error(f"CSV parsing error: {e}")
    except OSError as e:
        error(f"Failed to read samplesheet: {e}")


@validate_app.command("databases")
@validate_app.command("db", hidden=True)  # Alias
def validate_databases(
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Config file to validate (default: ~/.nvd2/config/user.config)",
    ),
) -> None:
    """Validate database paths from config file."""
    console.print("\n[bold]Validating Database Paths[/bold]\n")

    # Find config
    config_file = find_config_file(config)
    if not config_file:
        warning(f"No config file found at {DEFAULT_CONFIG}")
        info("Run install.sh to set up database configuration")
        sys.exit(1)

    info(f"Using config: {config_file}")

    # Parse config
    params = parse_nextflow_config(config_file)

    if not params:
        warning("No parameters found in config file")
        sys.exit(1)

    # Database parameters to check
    db_params = {
        "stat_index": "STAT index file",
        "stat_dbss": "STAT dbss file",
        "stat_annotation": "STAT annotation file",
        "human_virus_taxlist": "Human virus taxlist",
        "blast_db": "BLAST database directory",
        "gottcha2_db": "GOTTCHA2 database file",
    }

    console.print()
    all_found = True

    for param_name, description in db_params.items():
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
        success("All database paths are valid")
    else:
        console.print("[yellow]⚠ Some database paths are missing or invalid[/yellow]")
        info("Run install.sh to download and configure databases")
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

    # 2. Databases
    console.print("\n" + "=" * 70)
    try:
        validate_databases(config)
    except SystemExit as e:
        if e.code != 0:
            failed.append("Databases")

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
def validate_params_coverage(  # noqa: C901, PLR0912, PLR0915
    strict: bool = typer.Option(
        False,
        "--strict",
        help="Fail on any discrepancies (for CI)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Show all covered parameters",
    ),
) -> None:
    """
    Validate that CLI options cover all configurable Nextflow parameters.

    Checks that all user-configurable parameters in nextflow.config have
    corresponding CLI options in the 'nvd run' command.

    Exit codes:
    - 0: Full coverage (or warnings only without --strict)
    - 1: Missing coverage found
    """
    console.print("\n[bold cyan]Validating Parameter Coverage[/bold cyan]\n")

    # Locate nextflow.config
    config_path = Path.cwd() / "nextflow.config"
    if not config_path.exists():
        error(f"nextflow.config not found at: {config_path}")

    cli_path = Path(__file__).parent / "run.py"

    info(f"Nextflow config: {config_path}")
    info(f"CLI source:      {cli_path}")
    console.print()

    # Extract parameters
    console.print("[dim]Parsing Nextflow config...[/dim]")
    nf_params = extract_nextflow_params(config_path)

    console.print("[dim]Analyzing CLI options...[/dim]")
    cli_options = extract_cli_options()

    if not nf_params:
        error("Failed to parse Nextflow parameters")
    if not cli_options:
        error("Failed to parse CLI options")

    console.print()

    # Compare
    missing, extra, internal = compare_coverage(nf_params, cli_options)

    # Display results
    console.print("━" * 70)
    console.print("[bold]COVERAGE ANALYSIS[/bold]")
    console.print("━" * 70)
    console.print()

    # Covered parameters
    configurable = [
        param_info
        for param_info in nf_params.values()
        if param_info.category == "configurable"
    ]
    covered = [p for p in configurable if p.name in cli_options]

    console.print(f"[green]✓ Covered Parameters ({len(covered)}):[/green]")
    if verbose:
        for param in sorted(covered, key=lambda p: p.name):
            cli_opt = cli_options.get(param.name)
            flags = ", ".join(cli_opt.flags) if cli_opt and cli_opt.flags else "N/A"
            console.print(f"  • {param.name:30s} → {flags}")
    else:
        # Show first few
        for param in sorted(covered, key=lambda p: p.name)[:MAX_PREVIEW_ITEMS]:
            cli_opt = cli_options.get(param.name)
            flags = ", ".join(cli_opt.flags) if cli_opt and cli_opt.flags else "N/A"
            console.print(f"  • {param.name:30s} → {flags}")
        if len(covered) > MAX_PREVIEW_ITEMS:
            remaining = len(covered) - MAX_PREVIEW_ITEMS
            console.print(f"  [dim]... and {remaining} more[/dim]")
        console.print("  [dim](use --verbose to see all)[/dim]")

    console.print()

    # Missing parameters
    if missing:
        console.print(f"[yellow]⚠ Missing from CLI ({len(missing)}):[/yellow]")
        for param in sorted(missing, key=lambda p: p.name):
            console.print(f"  • [bold]{param.name}[/bold] (line {param.line_number})")
            console.print(f"    Default: {param.default_value}")
        console.print()

    # Internal parameters (informational)
    if internal:
        console.print(
            f"[cyan]i Internal Parameters (Correctly Excluded) ({len(internal)}):[/cyan]",
        )
        if verbose:
            for param in sorted(internal, key=lambda p: p.name):
                console.print(f"  • {param.name}")
        else:
            for param in sorted(internal, key=lambda p: p.name)[:MAX_PREVIEW_ITEMS]:
                console.print(f"  • {param.name}")
            if len(internal) > MAX_PREVIEW_ITEMS:
                remaining = len(internal) - MAX_PREVIEW_ITEMS
                console.print(f"  [dim]... and {remaining} more[/dim]")
        console.print()

    # Extra CLI options (informational, only in verbose)
    if extra and verbose:
        console.print(
            f"[dim]i CLI-only Options (Not in nextflow.config) ({len(extra)}):[/dim]",
        )
        for opt in sorted(extra):
            console.print(f"  • [dim]{opt}[/dim]")
        console.print()

    # Summary
    console.print("━" * 70)
    console.print("[bold]SUMMARY[/bold]")
    console.print("━" * 70)
    console.print()

    total_configurable = len(configurable)
    total_covered = len(covered)
    coverage_pct = (
        (total_covered / total_configurable * 100) if total_configurable > 0 else 0
    )

    console.print(
        f"Coverage: {total_covered}/{total_configurable} ({coverage_pct:.1f}%)",
    )
    console.print()

    if missing:
        console.print(
            f"[yellow]⚠ Status: Incomplete coverage ({len(missing)} parameters missing)[/yellow]",
        )
        console.print()
        if not strict:
            warning(
                "Missing parameters found (use --strict to fail on missing coverage)",
            )
        else:
            error(
                f"Validation failed: {len(missing)} configurable parameters are missing CLI options",
            )
    else:
        console.print("[green]✓ Status: Full coverage[/green]")
        console.print()
        success("All configurable parameters are covered by CLI options")
