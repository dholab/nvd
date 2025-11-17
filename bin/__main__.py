#!/usr/bin/env python3
# ruff: noqa: B008, FBT001, FBT002, FBT003
"""
NVD2 Command Line Interface

A user-friendly wrapper around the Nextflow-based NVD2 pipeline for
metagenomic virus detection and taxonomic classification.

Usage:
    nvd run --samplesheet samples.csv --experiment-id exp001
    nvd version

Note: B008 (function calls in defaults) and FBT (boolean traps) are disabled
because typer.Option() in function defaults is the idiomatic pattern for Typer,
and boolean CLI flags are standard and expected.
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

# ============================================================================
# CONSTANTS
# ============================================================================

VERSION = "0.1.0"
DEFAULT_CONFIG = Path.home() / ".nvd2" / "config" / "user.config"
VALID_TOOLS = [
    "stat_blast",
    "nvd",
    "stat",
    "blast",
    "stast",
    "gottcha",
    "all",
    "clumpify",
]
VALID_PROFILES = ["docker", "apptainer", "local", "chtc_hpc"]

# ============================================================================
# UTILITIES
# ============================================================================

console = Console()


def error(message: str, exit_code: int = 1) -> None:
    """Print error message and exit."""
    console.print(f"[red]✗ Error:[/red] {message}")
    sys.exit(exit_code)


def success(message: str) -> None:
    """Print success message."""
    console.print(f"[green]✓[/green] {message}")


def info(message: str) -> None:
    """Print info message."""
    console.print(f"[cyan]i[/cyan]  {message}")


def warning(message: str) -> None:
    """Print warning message."""
    console.print(f"[yellow]⚠[/yellow]  {message}")


# ============================================================================
# CONFIG MANAGEMENT
# ============================================================================


def find_config_file(custom_path: Path | None = None) -> Path | None:
    """Locate configuration file with fallback logic."""
    if custom_path:
        if custom_path.exists():
            return custom_path
        error(f"Config file not found: {custom_path}")

    if DEFAULT_CONFIG.exists():
        return DEFAULT_CONFIG

    return None


def parse_nextflow_config(config_path: Path) -> dict[str, str]:
    """
    Parse Nextflow config file for parameter values.
    Simple regex-based extraction of params { } block.
    """
    params = {}

    try:
        with open(config_path) as f:
            content = f.read()

        # Find params block
        params_match = re.search(r"params\s*\{([^}]+)\}", content, re.DOTALL)
        if not params_match:
            return params

        params_block = params_match.group(1)

        # Extract key = value pairs (simple string values)
        for raw_line in params_block.split("\n"):
            line = raw_line.strip()
            # Match: key = "value" or key = 'value' or key = value
            match = re.match(r'(\w+)\s*=\s*["\']?([^"\'\n]+)["\']?', line)
            if match:
                key, value = match.groups()
                params[key.strip()] = value.strip().strip("\"'")

    except OSError as e:
        warning(f"Failed to parse config file: {e}")
        return params
    else:
        return params


# ============================================================================
# ENVIRONMENT DETECTION
# ============================================================================


def check_command_exists(command: str) -> bool:
    """Check if a command exists in PATH."""
    return shutil.which(command) is not None


def docker_is_running() -> bool:
    """Check if Docker daemon is running."""
    try:
        result = subprocess.run(
            ["docker", "info"],  # noqa: S607
            capture_output=True,
            check=True,
            timeout=2,
        )
    except (
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
        FileNotFoundError,
    ):
        return False
    else:
        return result.returncode == 0


def auto_detect_profile() -> str:
    """Auto-detect available execution environment."""
    if check_command_exists("docker") and docker_is_running():
        return "docker"
    if check_command_exists("apptainer"):
        return "apptainer"
    if check_command_exists("singularity"):
        return "apptainer"  # Singularity is legacy name for Apptainer
    return "local"


# ============================================================================
# NEXTFLOW COMMAND CONSTRUCTION
# ============================================================================


def build_nextflow_command(  # noqa: PLR0913
    samplesheet: Path,
    experiment_id: str,
    tools: str,
    profile: str,
    config: Path | None,
    results: Path | None = None,
    resume: bool = False,
    cleanup: bool = False,
    work_dir: Path | None = None,
    **extra_params: str | Path | float | bool,
) -> list[str]:
    """Construct the nextflow run command."""

    cmd = ["nextflow", "run", "dhoconno/nvd"]

    # Profile
    cmd.extend(["-profile", profile])

    # Config file
    if config:
        cmd.extend(["-c", str(config)])

    # Work directory
    if work_dir:
        cmd.extend(["-work-dir", str(work_dir)])

    # Resume flag
    if resume:
        cmd.append("-resume")

    # Required pipeline parameters
    cmd.extend(
        [
            "--samplesheet",
            str(samplesheet),
            "--experiment_id",
            experiment_id,
            "--tools",
            tools,
        ],
    )

    # Optional pipeline parameters
    if results:
        cmd.extend(["--results", str(results)])

    if cleanup:
        cmd.append("--cleanup")

    # Add any extra parameters passed as database overrides
    for key, value in extra_params.items():
        if value is not None:
            # Convert Python arg names (gottcha2_db) to Nextflow params (--gottcha2-db)
            param_name = key.replace("_", "-")
            cmd.extend([f"--{param_name}", str(value)])

    return cmd


# ============================================================================
# TYPER APP & COMMANDS
# ============================================================================

app = typer.Typer(
    name="nvd",
    help="NVD2 - Metagenomic virus detection pipeline",
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    context_settings={"help_option_names": ["-h", "--help"]},
)


@app.command("run")
@app.command("r", hidden=True)  # Alias
def run(  # noqa: PLR0913, PLR0912, PLR0915, C901
    # Complexity/boolean/B008 warnings are acceptable for CLI with many options
    samplesheet: Path = typer.Option(
        ...,
        "--samplesheet",
        "-s",
        help="Path to samplesheet CSV",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    experiment_id: str = typer.Option(
        ...,
        "--experiment-id",
        "-e",
        help="Unique experiment identifier",
    ),
    tools: str = typer.Option(
        "all",
        "--tools",
        "-t",
        help="Tools to run: stat_blast, gottcha, or all",
    ),
    results: Path | None = typer.Option(
        None,
        "--results",
        "-r",
        help="Results directory (default: ./results)",
    ),
    profile: str | None = typer.Option(
        None,
        "--profile",
        "-p",
        help="Execution profile: docker, apptainer, local (auto-detect if not set)",
    ),
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Custom config file (default: ~/.nvd2/config/user.config)",
        exists=True,
    ),
    resume: bool = typer.Option(
        False,
        "--resume",
        help="Resume from checkpoint",
    ),
    cleanup: bool = typer.Option(
        False,
        "--cleanup",
        help="Clean work directory after success",
    ),
    work_dir: Path | None = typer.Option(
        None,
        "--work-dir",
        "-w",
        help="Nextflow work directory",
    ),
    # Database overrides
    gottcha2_db: Path | None = typer.Option(
        None,
        "--gottcha2-db",
        help="Override GOTTCHA2 database path",
    ),
    blast_db: Path | None = typer.Option(
        None,
        "--blast-db",
        help="Override BLAST database directory",
    ),
    blast_db_prefix: str | None = typer.Option(
        None,
        "--blast-db-prefix",
        help="Override BLAST database prefix",
    ),
    stat_index: Path | None = typer.Option(
        None,
        "--stat-index",
        help="Override STAT index file",
    ),
    stat_dbss: Path | None = typer.Option(
        None,
        "--stat-dbss",
        help="Override STAT dbss file",
    ),
    stat_annotation: Path | None = typer.Option(
        None,
        "--stat-annotation",
        help="Override STAT annotation file",
    ),
    human_virus_taxlist: Path | None = typer.Option(
        None,
        "--human-virus-taxlist",
        help="Override human virus taxlist file",
    ),
    # Analysis parameters
    cutoff_percent: float | None = typer.Option(
        None,
        "--cutoff-percent",
        help="Cutoff percentage (default: 0.001)",
    ),
    tax_stringency: float | None = typer.Option(
        None,
        "--tax-stringency",
        help="Taxonomy stringency (default: 0.7)",
    ),
    entropy: float | None = typer.Option(
        None,
        "--entropy",
        help="Entropy threshold (default: 0.9)",
    ),
    min_gottcha_reads: int | None = typer.Option(
        None,
        "--min-gottcha-reads",
        help="Minimum reads for GOTTCHA2 (default: 250)",
    ),
    max_blast_targets: int | None = typer.Option(
        None,
        "--max-blast-targets",
        help="Maximum BLAST targets (default: 100)",
    ),
    # LabKey
    labkey: bool = typer.Option(
        False,
        "--labkey",
        help="Enable LabKey integration",
    ),
    # Dry run
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show command without executing",
    ),
) -> None:
    """
    Run the NVD2 pipeline.

    This command wraps 'nextflow run dhoconno/nvd' with a simpler interface.
    Database paths and settings are loaded from ~/.nvd2/config/user.config
    unless overridden with command-line options.

    Examples:

        # Simple run with auto-detected profile
        nvd run -s samples.csv -e exp001

        # Specify profile and results directory
        nvd run -s samples.csv -e exp002 -p docker -r ./my_results

        # Run only GOTTCHA2 workflow
        nvd run -s samples.csv -e exp003 -t gottcha

        # Resume a failed run
        nvd run -s samples.csv -e exp003 --resume
    """

    # Validate tools option
    if tools not in VALID_TOOLS:
        error(
            f"Invalid tools option: {tools}. Must be one of: {', '.join(VALID_TOOLS)}",
        )

    # Auto-detect profile if not specified
    if not profile:
        profile = auto_detect_profile()
        info(f"Auto-detected execution profile: {profile}")
    elif profile not in VALID_PROFILES:
        valid_profiles = ", ".join(VALID_PROFILES)
        error(f"Invalid profile: {profile}. Must be one of: {valid_profiles}")

    # Find config file
    config_file = find_config_file(config)
    if config_file:
        info(f"Using config: {config_file}")
    else:
        warning("No config file found. Using command-line parameters only.")
        warning(f"Consider creating a config file at: {DEFAULT_CONFIG}")

    # Check Nextflow is installed
    if not check_command_exists("nextflow"):
        error(
            "Nextflow not found. Please install Nextflow:\n"
            "  curl -s https://get.nextflow.io | bash\n"
            "  Or visit: https://www.nextflow.io/docs/latest/getstarted.html",
        )

    # Gather database overrides and extra params
    extra_params = {}
    if gottcha2_db:
        extra_params["gottcha2_db"] = gottcha2_db
    if blast_db:
        extra_params["blast_db"] = blast_db
    if blast_db_prefix:
        extra_params["blast_db_prefix"] = blast_db_prefix
    if stat_index:
        extra_params["stat_index"] = stat_index
    if stat_dbss:
        extra_params["stat_dbss"] = stat_dbss
    if stat_annotation:
        extra_params["stat_annotation"] = stat_annotation
    if human_virus_taxlist:
        extra_params["human_virus_taxlist"] = human_virus_taxlist
    if cutoff_percent is not None:
        extra_params["cutoff_percent"] = cutoff_percent
    if tax_stringency is not None:
        extra_params["tax_stringency"] = tax_stringency
    if entropy is not None:
        extra_params["entropy"] = entropy
    if min_gottcha_reads is not None:
        extra_params["min_gottcha_reads"] = min_gottcha_reads
    if max_blast_targets is not None:
        extra_params["max_blast_targets"] = max_blast_targets
    if labkey:
        extra_params["labkey"] = "true"

    # Build command
    cmd = build_nextflow_command(
        samplesheet=samplesheet,
        experiment_id=experiment_id,
        tools=tools,
        profile=profile,
        config=config_file,
        results=results,
        resume=resume,
        cleanup=cleanup,
        work_dir=work_dir,
        **extra_params,
    )

    # Show command
    console.print("\n[bold]Executing command:[/bold]")
    console.print(f"[dim]{' '.join(cmd)}[/dim]\n")

    if dry_run:
        success("Dry-run mode: command shown above but not executed")
        return

    # Execute nextflow
    try:
        result = subprocess.run(cmd, check=False)  # noqa: S603
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        console.print("\n[yellow]Pipeline interrupted by user[/yellow]")
        sys.exit(130)
    except OSError as e:
        error(f"Failed to execute Nextflow: {e}")


@app.command("version")
@app.command("v", hidden=True)  # Alias
def version() -> None:  # noqa: PLR0912
    """Show version information."""

    console.print(
        Panel.fit(
            f"[bold cyan]NVD2 CLI[/bold cyan]\n"
            f"Version: {VERSION}\n\n"
            f"A user-friendly wrapper for the NVD2 Nextflow pipeline.\n"
            f"Pipeline: dhoconno/nvd",
            title="Version Info",
            border_style="cyan",
        ),
    )

    # Show dependency versions if available
    deps = []

    if check_command_exists("nextflow"):
        try:
            result = subprocess.run(
                ["nextflow", "-version"],  # noqa: S607
                check=False,
                capture_output=True,
                text=True,
                timeout=3,
            )
            nf_version = result.stdout.strip().split("\n")[0]
            deps.append(f"✓ {nf_version}")
        except (OSError, subprocess.TimeoutExpired):
            deps.append("✓ Nextflow (version unknown)")
    else:
        deps.append("✗ Nextflow not found")

    if check_command_exists("java"):
        try:
            result = subprocess.run(
                ["java", "-version"],  # noqa: S607
                check=False,
                capture_output=True,
                text=True,
                timeout=3,
            )
            java_output = result.stderr.strip().split("\n")[0]
            deps.append(f"✓ {java_output}")
        except (OSError, subprocess.TimeoutExpired):
            deps.append("✓ Java (version unknown)")
    else:
        deps.append("✗ Java not found")

    # Docker/Apptainer
    if check_command_exists("docker"):
        if docker_is_running():
            deps.append("✓ Docker (running)")
        else:
            deps.append("⚠ Docker (daemon not running)")
    elif check_command_exists("apptainer"):
        deps.append("✓ Apptainer")
    elif check_command_exists("singularity"):
        deps.append("✓ Singularity")
    else:
        deps.append("⚠ No container runtime found")

    console.print("\n[bold]Dependencies:[/bold]")
    for dep in deps:
        console.print(f"  {dep}")

    console.print()


# ============================================================================
# VALIDATION COMMANDS
# ============================================================================

validate_app = typer.Typer(
    name="validate",
    help="Validation commands to check system setup",
    no_args_is_help=True,
)
app.add_typer(validate_app, name="validate")
app.add_typer(validate_app, name="val", hidden=True)  # Alias


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
                    msg = f"Line {i} ({sample_id}): must provide either 'srr' or 'fastq1'"
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


# ============================================================================
# CONFIG COMMANDS
# ============================================================================

config_app = typer.Typer(
    name="config",
    help="Configuration management commands",
    no_args_is_help=True,
)
app.add_typer(config_app, name="config")
app.add_typer(config_app, name="cfg", hidden=True)  # Alias


@config_app.command("show")
@config_app.command("s", hidden=True)  # Alias
def config_show(
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Config file to display (default: ~/.nvd2/config/user.config)",
    ),
) -> None:
    """Display current configuration."""
    # Find config
    config_file = find_config_file(config)
    if not config_file:
        warning(f"No config file found at {DEFAULT_CONFIG}")
        info("Run install.sh to create a configuration file")
        sys.exit(1)

    console.print(f"\n[bold]Configuration File:[/bold] {config_file}\n")

    # Parse config
    params = parse_nextflow_config(config_file)

    if not params:
        warning("No parameters found in config file")
        sys.exit(1)

    # Display as table
    table = Table(
        title="Configuration Parameters",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Parameter", style="cyan", width=30)
    table.add_column("Value", style="white")

    for key, value in sorted(params.items()):
        table.add_row(key, value)

    console.print(table)
    console.print()


@config_app.command("path")
@config_app.command("p", hidden=True)  # Alias
def config_path() -> None:
    """Show configuration file location."""
    console.print(f"\n[bold]Default config location:[/bold] {DEFAULT_CONFIG}")

    if DEFAULT_CONFIG.exists():
        success(f"Config file exists: {DEFAULT_CONFIG}")
        file_size = DEFAULT_CONFIG.stat().st_size
        console.print(f"  Size: {file_size} bytes")
    else:
        warning("Config file does not exist")
        info("Run install.sh to create a configuration file")

    console.print()


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================


def main() -> None:
    """Main entry point for the CLI."""
    try:
        app()
    except KeyboardInterrupt:
        console.print("\n[yellow]Interrupted by user[/yellow]")
        sys.exit(130)
    except (typer.Exit, typer.Abort):
        # Normal exit from Typer
        raise
    except OSError as e:
        console.print(f"[red]Error:[/red] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
