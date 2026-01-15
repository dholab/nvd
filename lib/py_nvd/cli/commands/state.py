# ruff: noqa: B008
"""
State management commands for the NVD CLI.

Commands for inspecting and managing the NVD pipeline state database,
which tracks runs, samples, uploads, databases, and presets.

Commands:
    nvd state path  - Show resolved state directory path
    nvd state init  - Initialize state database
    nvd state info  - Show database overview
    nvd state runs  - List pipeline runs
    nvd state run   - Show details for a specific run
    nvd state move  - Move state directory to a new location
"""

from __future__ import annotations

import json
import os
from datetime import date, datetime
from pathlib import Path
from typing import Annotated

import typer
from rich.table import Table

from py_nvd import __version__, state
from py_nvd.cli.utils import (
    console,
    ensure_db_exists,
    error,
    format_duration,
    info,
    success,
    warning,
)
from py_nvd.db import (
    DEFAULT_STATE_DIR,
    ENV_VAR,
    EXPECTED_VERSION,
    connect,
    get_state_db_path,
    get_state_dir,
    get_taxdump_dir,
)
from py_nvd.hits import (
    count_hit_observations,
    count_hits,
    count_sample_set_observations,
    delete_sample_set_hits,
)
from py_nvd.models import ProcessedSampleStatus, Status, UploadType

state_app = typer.Typer(
    name="state",
    help="Inspect and manage pipeline state (runs, samples, uploads, databases)",
    no_args_is_help=True,
)

# Database subcommand group
database_app = typer.Typer(
    name="database",
    help="Manage reference database registrations",
    no_args_is_help=True,
)
state_app.add_typer(database_app, name="database")


def _format_size(size_bytes: int) -> str:
    """Format bytes as human-readable size."""
    size = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def _output_json(data: dict | list) -> None:
    """Print data as formatted JSON."""
    console.print_json(json.dumps(data, default=str))


def _validate_path_for_sql(path: Path) -> Path:
    """
    Validate that a path can be safely used in SQL ATTACH statement.

    SQLite's ATTACH DATABASE requires the path as a string literal,
    which cannot be parameterized. This validates the path is safe.

    Args:
        path: The path to validate

    Returns:
        The resolved absolute path

    Raises:
        ValueError: If path contains characters unsafe for SQL string literals
    """
    resolved = path.resolve()
    path_str = str(resolved)
    if "'" in path_str:
        raise ValueError(
            f"Path contains single quote which cannot be safely used in SQL: {path_str}\n"
            "Please rename the file to remove the quote.",
        )
    return resolved


@state_app.command("path")
def state_path(
    json_output: bool = typer.Option(
        False,
        "--json",
        help="Output as JSON",
    ),
) -> None:
    """
    Show resolved state directory path and precedence rule used.

    Displays which method was used to determine the state directory:
    explicit path, environment variable, or default location.
    """
    # Determine which precedence rule applies
    env_value = os.environ.get(ENV_VAR)

    if env_value:
        source = "environment"
        source_detail = f"{ENV_VAR}={env_value}"
        resolved_path = Path(env_value)
    else:
        source = "default"
        source_detail = f"~/.nvd/ (no {ENV_VAR} set)"
        resolved_path = DEFAULT_STATE_DIR

    db_path = get_state_db_path()
    db_exists = db_path.exists()
    db_size = db_path.stat().st_size if db_exists else 0

    taxdump_dir = get_taxdump_dir()
    taxdump_exists = taxdump_dir.exists()

    if json_output:
        _output_json(
            {
                "state_dir": str(resolved_path),
                "source": source,
                "source_detail": source_detail,
                "database": {
                    "path": str(db_path),
                    "exists": db_exists,
                    "size_bytes": db_size,
                },
                "taxdump": {
                    "path": str(taxdump_dir),
                    "exists": taxdump_exists,
                },
            },
        )
        return

    console.print(f"\n[bold]State Directory[/bold]: {resolved_path}")
    console.print(f"[dim]Resolved via:[/dim] {source_detail}")
    console.print()

    if db_exists:
        success(f"Database exists: {db_path}")
        console.print(f"  Size: {_format_size(db_size)}")
    else:
        warning(f"Database not found: {db_path}")
        info("Run 'nvd state init' to create the database")

    if taxdump_exists:
        success(f"Taxdump directory exists: {taxdump_dir}")
    else:
        info(f"Taxdump directory not found: {taxdump_dir}")

    console.print()


@state_app.command("init")
def state_init(
    path: Path | None = typer.Option(
        None,
        "--path",
        "-p",
        help="Explicit state directory path (overrides NVD_STATE_DIR)",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Reinitialize even if database exists",
    ),
) -> None:
    """
    Initialize the state database.

    Creates the state directory and initializes an empty database.
    If no --path is specified, uses NVD_STATE_DIR or the default location.

    Examples:

        # Initialize at default location
        nvd state init

        # Initialize at custom location
        nvd state init --path /scratch/nvd-state

        # Reinitialize (recreate) existing database
        nvd state init --force
    """
    # Resolve state directory
    state_dir = get_state_dir(path)
    db_path = get_state_db_path(path)

    if db_path.exists() and not force:
        warning(f"Database already exists: {db_path}")
        info("Use --force to reinitialize")
        raise typer.Exit(1)

    if db_path.exists() and force:
        info(f"Removing existing database: {db_path}")
        db_path.unlink()

    # Initialize by connecting (this creates schema)
    with connect(path):
        pass  # Schema is created on connect

    success(f"Initialized state database: {db_path}")
    console.print()

    # Show how to persist this location
    if path:
        console.print("[bold]To use this location by default:[/bold]")
        console.print(f"  export {ENV_VAR}={state_dir}")
        console.print()
        console.print(
            "[dim]Add this to your shell profile (.bashrc, .zshrc, etc.)[/dim]",
        )
        console.print()


@state_app.command("info")
def state_info(
    json_output: bool = typer.Option(
        False,
        "--json",
        help="Output as JSON",
    ),
    skip_hits: bool = typer.Option(
        False,
        "--skip-hits",
        help="Skip querying hits data (faster when hits database is large)",
    ),
) -> None:
    """
    Show overview of state database contents.

    Displays counts of runs, samples, uploads, databases, and presets,
    along with database size and schema version.

    Examples:

        # Show full overview
        nvd state info

        # Skip hits query for faster output
        nvd state info --skip-hits

        # JSON output
        nvd state info --json
    """
    db_path = ensure_db_exists(json_output)

    db_size = db_path.stat().st_size

    with connect() as conn:
        # Get schema version
        version = conn.execute("PRAGMA user_version").fetchone()[0]

        # Get counts
        runs_count = conn.execute("SELECT COUNT(*) FROM runs").fetchone()[0]
        samples_count = conn.execute(
            "SELECT COUNT(*) FROM processed_samples",
        ).fetchone()[0]
        uploads_count = conn.execute("SELECT COUNT(*) FROM uploads").fetchone()[0]
        databases_count = conn.execute("SELECT COUNT(*) FROM databases").fetchone()[0]
        presets_count = conn.execute("SELECT COUNT(*) FROM presets").fetchone()[0]
        taxonomy_count = conn.execute(
            "SELECT COUNT(*) FROM taxonomy_versions",
        ).fetchone()[0]

        # Get lock counts (handle case where table doesn't exist yet)
        try:
            locks_count = conn.execute("SELECT COUNT(*) FROM sample_locks").fetchone()[
                0
            ]
            expired_locks_count = conn.execute(
                "SELECT COUNT(*) FROM sample_locks WHERE expires_at < datetime('now')",
            ).fetchone()[0]
            active_locks_count = locks_count - expired_locks_count
        except Exception:
            locks_count = 0
            expired_locks_count = 0
            active_locks_count = 0

        # Get recent activity
        latest_run = conn.execute(
            "SELECT run_id, started_at FROM runs ORDER BY started_at DESC LIMIT 1",
        ).fetchone()

    # Get hits counts from parquet files (not SQLite)
    if skip_hits:
        hits_count = None
        hit_observations_count = None
    else:
        hits_count = count_hits()
        hit_observations_count = count_hit_observations()

    # Check taxdump
    taxdump_dir = get_taxdump_dir()
    taxdump_exists = taxdump_dir.exists()
    taxdump_files = list(taxdump_dir.glob("*.dmp")) if taxdump_exists else []

    if json_output:
        _output_json(
            {
                "database": {
                    "path": str(db_path),
                    "size_bytes": db_size,
                    "schema_version": version,
                    "expected_version": EXPECTED_VERSION,
                },
                "counts": {
                    "runs": runs_count,
                    "samples": samples_count,
                    "uploads": uploads_count,
                    "databases": databases_count,
                    "presets": presets_count,
                    "taxonomy_versions": taxonomy_count,
                    "hits": hits_count,
                    "hit_observations": hit_observations_count,
                    "locks_active": active_locks_count,
                    "locks_expired": expired_locks_count,
                },
                "latest_run": {
                    "run_id": latest_run["run_id"] if latest_run else None,
                    "started_at": latest_run["started_at"] if latest_run else None,
                },
                "taxdump": {
                    "path": str(taxdump_dir),
                    "exists": taxdump_exists,
                    "file_count": len(taxdump_files),
                },
            },
        )
        return

    console.print(f"\n[bold]State Database:[/bold] {db_path}")
    console.print(f"  Size: {_format_size(db_size)}")
    console.print(f"  Schema version: {version}", end="")
    if version == EXPECTED_VERSION:
        console.print(" [green](current)[/green]")
    else:
        console.print(f" [yellow](expected: {EXPECTED_VERSION})[/yellow]")
    console.print()

    # Counts table
    table = Table(title="Contents", show_header=True, header_style="bold cyan")
    table.add_column("Entity", style="cyan")
    table.add_column("Count", justify="right")

    table.add_row("Runs", str(runs_count))
    table.add_row("Processed Samples", str(samples_count))
    table.add_row("Uploads", str(uploads_count))
    table.add_row("Databases", str(databases_count))
    table.add_row("Presets", str(presets_count))
    table.add_row("Taxonomy Versions", str(taxonomy_count))
    hits_display = "[dim]skipped[/dim]" if hits_count is None else str(hits_count)
    obs_display = (
        "[dim]skipped[/dim]"
        if hit_observations_count is None
        else str(hit_observations_count)
    )
    table.add_row("Hits", hits_display)
    table.add_row("Hit Observations", obs_display)
    if active_locks_count > 0 or expired_locks_count > 0:
        lock_str = str(active_locks_count)
        if expired_locks_count > 0:
            lock_str += f" [dim](+{expired_locks_count} expired)[/dim]"
        table.add_row("Sample Locks", lock_str)

    console.print(table)
    console.print()

    if latest_run:
        console.print(f"[bold]Latest Run:[/bold] {latest_run['run_id']}")
        console.print(f"  Started: {latest_run['started_at']}")
        console.print()

    if taxdump_exists:
        console.print(f"[bold]Taxonomy Cache:[/bold] {taxdump_dir}")
        console.print(f"  Files: {len(taxdump_files)} .dmp files")
    else:
        console.print(f"[dim]Taxonomy cache not found: {taxdump_dir}[/dim]")

    console.print()


def _parse_since(value: str) -> date | datetime:
    """
    Parse a date or datetime string in ISO format.

    Accepts:
        - Date: "2024-12-01"
        - Datetime: "2024-12-01T14:30:00" or "2024-12-01 14:30:00"

    Raises:
        typer.BadParameter: If the string cannot be parsed
    """
    try:
        # Try datetime first (more specific)
        return datetime.fromisoformat(value)
    except ValueError:
        pass

    try:
        return date.fromisoformat(value)
    except ValueError:
        raise typer.BadParameter(
            f"Invalid date format: {value!r}. Use ISO format: YYYY-MM-DD or YYYY-MM-DDTHH:MM:SS",
        )


# Status values for CLI validation
_STATUS_VALUES = ["running", "completed", "failed"]


@state_app.command("runs")
def state_runs(
    status: Annotated[
        str | None,
        typer.Option(
            "--status",
            "-s",
            help=f"Filter by status ({', '.join(_STATUS_VALUES)})",
        ),
    ] = None,
    since: Annotated[
        str | None,
        typer.Option(
            "--since",
            help="Filter to runs started on or after this date (ISO format: YYYY-MM-DD)",
        ),
    ] = None,
    limit: Annotated[
        int | None,
        typer.Option(
            "--limit",
            "-n",
            help="Maximum number of runs to show",
        ),
    ] = None,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    List pipeline runs.

    Shows all recorded pipeline runs with their status and timing information.
    Results are ordered by start time, most recent first.

    Examples:

        # List all runs
        nvd state runs

        # List only completed runs
        nvd state runs --status completed

        # List runs from the last week
        nvd state runs --since 2024-12-22

        # List the 5 most recent failed runs
        nvd state runs --status failed --limit 5
    """
    db_path = ensure_db_exists(json_output)

    # Validate and convert status
    status_filter: Status | None = None
    if status is not None:
        if status not in _STATUS_VALUES:
            error(
                f"Invalid status: {status!r}. Must be one of: {', '.join(_STATUS_VALUES)}",
            )
            raise typer.Exit(1)
        status_filter = status  # type: ignore[assignment]

    # Parse since date
    since_filter: date | datetime | None = None
    if since is not None:
        since_filter = _parse_since(since)

    runs = state.list_runs(status=status_filter, since=since_filter, limit=limit)

    if json_output:
        _output_json(
            [
                {
                    "run_id": r.run_id,
                    "sample_set_id": r.sample_set_id,
                    "status": r.status,
                    "started_at": r.started_at,
                    "completed_at": r.completed_at,
                    "experiment_id": r.experiment_id,
                }
                for r in runs
            ],
        )
        return

    if not runs:
        info("No runs found")
        if status or since:
            console.print("[dim]Try adjusting your filters[/dim]")
        return

    table = Table(title="Pipeline Runs", show_header=True, header_style="bold cyan")
    table.add_column("Run ID", style="cyan")
    table.add_column("Status")
    table.add_column("Started")
    table.add_column("Completed")
    table.add_column("Experiment ID", justify="right")

    for run in runs:
        # Color-code status
        if run.status == "completed":
            status_str = "[green]completed[/green]"
        elif run.status == "failed":
            status_str = "[red]failed[/red]"
        else:
            status_str = "[yellow]running[/yellow]"

        table.add_row(
            run.run_id,
            status_str,
            run.started_at,
            run.completed_at or "-",
            str(run.experiment_id) if run.experiment_id else "-",
        )

    console.print()
    console.print(table)
    console.print(f"\n[dim]Showing {len(runs)} run(s)[/dim]")
    console.print()


@state_app.command("run")
def state_run(
    run_id: Annotated[
        str,
        typer.Argument(help="The run ID to show details for"),
    ],
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Show details for a specific run.

    Displays run metadata, associated samples, and taxonomy version information.

    Examples:

        # Show run details
        nvd state run friendly_turing

        # Output as JSON for scripting
        nvd state run friendly_turing --json
    """
    db_path = ensure_db_exists(json_output)

    run = state.get_run(run_id)
    if run is None:
        if json_output:
            _output_json({"error": "Run not found", "run_id": run_id})
            raise typer.Exit(1)
        error(f"Run not found: {run_id}")

    samples = state.get_samples_for_run(run_id)
    taxonomy = state.get_taxonomy_version(run_id)

    if json_output:
        _output_json(
            {
                "run": {
                    "run_id": run.run_id,
                    "sample_set_id": run.sample_set_id,
                    "status": run.status,
                    "started_at": run.started_at,
                    "completed_at": run.completed_at,
                    "duration_seconds": run.duration_seconds,
                    "experiment_id": run.experiment_id,
                },
                "samples": [
                    {
                        "sample_id": s.sample_id,
                        "status": s.status,
                        "processed_at": s.processed_at,
                        "blast_db_version": s.blast_db_version,
                        "stat_db_version": s.stat_db_version,
                        "taxonomy_hash": s.taxonomy_hash,
                    }
                    for s in samples
                ],
                "taxonomy": {
                    "file_hash": taxonomy.file_hash,
                    "ncbi_rebuild_timestamp": taxonomy.ncbi_rebuild_timestamp,
                    "file_size": taxonomy.file_size,
                    "downloaded_at": taxonomy.downloaded_at,
                }
                if taxonomy
                else None,
            },
        )
        return

    # Human-readable output
    console.print()

    # Run details
    if run.status == "completed":
        status_str = "[green]completed[/green]"
    elif run.status == "failed":
        status_str = "[red]failed[/red]"
    else:
        status_str = "[yellow]running[/yellow]"

    console.print(f"[bold]Run:[/bold] {run.run_id}")
    console.print(f"  Status: {status_str}")
    console.print(f"  Sample Set ID: {run.sample_set_id}")
    console.print(f"  Started: {run.started_at}")
    if run.completed_at:
        console.print(f"  Completed: {run.completed_at}")
    if run.duration_seconds is not None:
        console.print(f"  Duration: {format_duration(run.duration_seconds)}")
    if run.experiment_id:
        console.print(f"  Experiment ID: {run.experiment_id}")
    console.print()

    # Samples
    if samples:
        table = Table(
            title=f"Samples ({len(samples)})",
            show_header=True,
            header_style="bold cyan",
        )
        table.add_column("Sample ID", style="cyan")
        table.add_column("Status")
        table.add_column("Processed At")
        table.add_column("BLAST DB")
        table.add_column("STAT DB")

        for sample in samples:
            if sample.status == "uploaded":
                sample_status = "[green]uploaded[/green]"
            elif sample.status == "completed":
                sample_status = "[cyan]completed[/cyan]"
            elif sample.status == "failed":
                sample_status = "[red]failed[/red]"
            else:
                sample_status = f"[yellow]{sample.status}[/yellow]"

            table.add_row(
                sample.sample_id,
                sample_status,
                sample.processed_at,
                sample.blast_db_version or "-",
                sample.stat_db_version or "-",
            )

        console.print(table)
        console.print()
    else:
        info("No samples recorded for this run")
        console.print()

    # Taxonomy
    if taxonomy:
        console.print("[bold]Taxonomy Version:[/bold]")
        console.print(f"  File Hash: {taxonomy.file_hash[:16]}...")
        if taxonomy.ncbi_rebuild_timestamp:
            console.print(f"  NCBI Rebuild: {taxonomy.ncbi_rebuild_timestamp}")
        if taxonomy.file_size:
            console.print(f"  File Size: {_format_size(taxonomy.file_size)}")
        console.print(f"  Downloaded: {taxonomy.downloaded_at}")
        console.print()


# Sample status values for CLI validation
_SAMPLE_STATUS_VALUES = ["completed", "uploaded", "failed"]


@state_app.command("samples")
def state_samples(
    run_id: Annotated[
        str | None,
        typer.Option(
            "--run",
            "-r",
            help="Filter to samples from a specific run",
        ),
    ] = None,
    status: Annotated[
        str | None,
        typer.Option(
            "--status",
            "-s",
            help=f"Filter by status ({', '.join(_SAMPLE_STATUS_VALUES)})",
        ),
    ] = None,
    limit: Annotated[
        int | None,
        typer.Option(
            "--limit",
            "-n",
            help="Maximum number of samples to show",
        ),
    ] = None,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    List processed samples.

    Shows all processed samples with their status and run information.
    Results are ordered by processing time, most recent first.

    Examples:

        # List all samples
        nvd state samples

        # List samples from a specific run
        nvd state samples --run friendly_turing

        # List failed samples
        nvd state samples --status failed

        # List the 10 most recent samples
        nvd state samples --limit 10
    """
    db_path = ensure_db_exists(json_output)

    # Validate status
    status_filter: ProcessedSampleStatus | None = None
    if status is not None:
        if status not in _SAMPLE_STATUS_VALUES:
            error(
                f"Invalid status: {status!r}. Must be one of: {', '.join(_SAMPLE_STATUS_VALUES)}",
            )
            raise typer.Exit(1)
        status_filter = status  # type: ignore[assignment]

    samples = state.list_samples(run_id=run_id, status=status_filter, limit=limit)

    if json_output:
        _output_json(
            [
                {
                    "sample_id": s.sample_id,
                    "run_id": s.run_id,
                    "status": s.status,
                    "processed_at": s.processed_at,
                    "blast_db_version": s.blast_db_version,
                    "stat_db_version": s.stat_db_version,
                    "taxonomy_hash": s.taxonomy_hash,
                }
                for s in samples
            ],
        )
        return

    if not samples:
        info("No samples found")
        if run_id or status:
            console.print("[dim]Try adjusting your filters[/dim]")
        return

    table = Table(title="Processed Samples", show_header=True, header_style="bold cyan")
    table.add_column("Sample ID", style="cyan")
    table.add_column("Run ID")
    table.add_column("Status")
    table.add_column("Processed At")
    table.add_column("BLAST DB")

    for sample in samples:
        if sample.status == "uploaded":
            status_str = "[green]uploaded[/green]"
        elif sample.status == "completed":
            status_str = "[cyan]completed[/cyan]"
        elif sample.status == "failed":
            status_str = "[red]failed[/red]"
        else:
            status_str = f"[yellow]{sample.status}[/yellow]"

        table.add_row(
            sample.sample_id,
            sample.run_id,
            status_str,
            sample.processed_at,
            sample.blast_db_version or "-",
        )

    console.print()
    console.print(table)
    console.print(f"\n[dim]Showing {len(samples)} sample(s)[/dim]")
    console.print()


@state_app.command("sample")
def state_sample(
    sample_id: Annotated[
        str,
        typer.Argument(help="The sample ID to show history for"),
    ],
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Show processing and upload history for a sample.

    Displays all times this sample was processed (across different runs)
    and its complete upload history.

    Examples:

        # Show sample history
        nvd state sample sample_a

        # Output as JSON for scripting
        nvd state sample sample_a --json
    """
    db_path = ensure_db_exists(json_output)

    processing_history = state.get_sample_history(sample_id)
    upload_history = state.get_uploads_for_sample(sample_id)

    if not processing_history and not upload_history:
        if json_output:
            _output_json({"error": "Sample not found", "sample_id": sample_id})
            raise typer.Exit(1)
        error(f"Sample not found: {sample_id}")

    if json_output:
        _output_json(
            {
                "sample_id": sample_id,
                "processing_history": [
                    {
                        "run_id": p.run_id,
                        "status": p.status,
                        "processed_at": p.processed_at,
                        "blast_db_version": p.blast_db_version,
                        "stat_db_version": p.stat_db_version,
                        "taxonomy_hash": p.taxonomy_hash,
                    }
                    for p in processing_history
                ],
                "upload_history": [
                    {
                        "upload_type": u.upload_type,
                        "upload_target": u.upload_target,
                        "uploaded_at": u.uploaded_at,
                        "content_hash": u.content_hash,
                        "target_metadata": u.target_metadata,
                    }
                    for u in upload_history
                ],
            },
        )
        return

    # Human-readable output
    console.print()
    console.print(f"[bold]Sample:[/bold] {sample_id}")
    console.print()

    # Processing history
    if processing_history:
        table = Table(
            title=f"Processing History ({len(processing_history)})",
            show_header=True,
            header_style="bold cyan",
        )
        table.add_column("Run ID", style="cyan")
        table.add_column("Status")
        table.add_column("Processed At")
        table.add_column("BLAST DB")
        table.add_column("STAT DB")

        for record in processing_history:
            if record.status == "uploaded":
                status_str = "[green]uploaded[/green]"
            elif record.status == "completed":
                status_str = "[cyan]completed[/cyan]"
            elif record.status == "failed":
                status_str = "[red]failed[/red]"
            else:
                status_str = f"[yellow]{record.status}[/yellow]"

            table.add_row(
                record.run_id,
                status_str,
                record.processed_at,
                record.blast_db_version or "-",
                record.stat_db_version or "-",
            )

        console.print(table)
        console.print()
    else:
        info("No processing records found")
        console.print()

    # Upload history
    if upload_history:
        table = Table(
            title=f"Upload History ({len(upload_history)})",
            show_header=True,
            header_style="bold cyan",
        )
        table.add_column("Type", style="cyan")
        table.add_column("Target")
        table.add_column("Uploaded At")
        table.add_column("Content Hash")

        for upload in upload_history:
            # Truncate hash for display
            short_hash = (
                upload.content_hash[:12] + "..." if upload.content_hash else "-"
            )

            table.add_row(
                upload.upload_type,
                upload.upload_target,
                upload.uploaded_at,
                short_hash,
            )

        console.print(table)
        console.print()
    else:
        info("No upload records found")
        console.print()


# Upload type values for CLI validation
_UPLOAD_TYPE_VALUES = ["blast", "blast_fasta", "gottcha2", "gottcha2_fasta"]


@state_app.command("uploads")
def state_uploads(
    sample_id: Annotated[
        str | None,
        typer.Option(
            "--sample",
            "-s",
            help="Filter to uploads for a specific sample",
        ),
    ] = None,
    upload_type: Annotated[
        str | None,
        typer.Option(
            "--type",
            "-t",
            help=f"Filter by upload type ({', '.join(_UPLOAD_TYPE_VALUES)})",
        ),
    ] = None,
    target: Annotated[
        str | None,
        typer.Option(
            "--target",
            help="Filter by upload target (e.g., labkey, local)",
        ),
    ] = None,
    limit: Annotated[
        int | None,
        typer.Option(
            "--limit",
            "-n",
            help="Maximum number of uploads to show",
        ),
    ] = None,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    List uploads.

    Shows all upload records with their type, target, and timing information.
    Results are ordered by upload time, most recent first.

    Examples:

        # List all uploads
        nvd state uploads

        # List uploads for a specific sample
        nvd state uploads --sample sample_a

        # List only BLAST uploads to LabKey
        nvd state uploads --type blast --target labkey

        # List the 10 most recent uploads
        nvd state uploads --limit 10
    """
    db_path = ensure_db_exists(json_output)

    # Validate upload_type
    type_filter: UploadType | None = None
    if upload_type is not None:
        if upload_type not in _UPLOAD_TYPE_VALUES:
            error(
                f"Invalid upload type: {upload_type!r}. Must be one of: {', '.join(_UPLOAD_TYPE_VALUES)}",
            )
            raise typer.Exit(1)
        type_filter = upload_type  # type: ignore[assignment]

    uploads = state.list_uploads(
        sample_id=sample_id,
        upload_type=type_filter,
        upload_target=target,
        limit=limit,
    )

    if json_output:
        _output_json(
            [
                {
                    "sample_id": u.sample_id,
                    "upload_type": u.upload_type,
                    "upload_target": u.upload_target,
                    "uploaded_at": u.uploaded_at,
                    "content_hash": u.content_hash,
                    "target_metadata": u.target_metadata,
                }
                for u in uploads
            ],
        )
        return

    if not uploads:
        info("No uploads found")
        if sample_id or upload_type or target:
            console.print("[dim]Try adjusting your filters[/dim]")
        return

    table = Table(title="Uploads", show_header=True, header_style="bold cyan")
    table.add_column("Sample ID", style="cyan")
    table.add_column("Type")
    table.add_column("Target")
    table.add_column("Uploaded At")
    table.add_column("Content Hash")

    for upload in uploads:
        short_hash = upload.content_hash[:12] + "..." if upload.content_hash else "-"

        table.add_row(
            upload.sample_id,
            upload.upload_type,
            upload.upload_target,
            upload.uploaded_at,
            short_hash,
        )

    console.print()
    console.print(table)
    console.print(f"\n[dim]Showing {len(uploads)} upload(s)[/dim]")
    console.print()


@database_app.command("list")
def database_list(
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Output as JSON"),
    ] = False,
) -> None:
    """
    List registered reference databases.

    Shows all databases registered with the state system, including their
    type, version, path, and validation status.

    Examples:

        nvd state database list
        nvd state database list --json
    """
    db_path = ensure_db_exists(json_output)

    databases = state.list_databases()

    if json_output:
        _output_json(
            [
                {
                    "db_type": db.db_type,
                    "version": db.version,
                    "path": db.path,
                    "registered_at": db.registered_at,
                    "path_exists": Path(db.path).exists(),
                }
                for db in databases
            ],
        )
        return

    if not databases:
        info("No databases registered")
        console.print("[dim]Use 'nvd state database register' to add databases[/dim]")
        return

    table = Table(
        title="Registered Databases",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Type", style="cyan")
    table.add_column("Version")
    table.add_column("Path")
    table.add_column("Status")
    table.add_column("Registered At")

    for db in databases:
        path_exists = Path(db.path).exists()
        status_str = "[green]valid[/green]" if path_exists else "[red]missing[/red]"

        # Truncate path if too long
        display_path = db.path
        if len(display_path) > 40:
            display_path = "..." + display_path[-37:]

        table.add_row(
            db.db_type,
            db.version,
            display_path,
            status_str,
            db.registered_at,
        )

    console.print()
    console.print(table)
    console.print(f"\n[dim]Showing {len(databases)} database(s)[/dim]")
    console.print()


@database_app.command("show")
def database_show(
    db_type: Annotated[
        str,
        typer.Argument(help="Database type (blast, stat, gottcha2)"),
    ],
    version: Annotated[
        str | None,
        typer.Argument(help="Version (optional, shows latest if omitted)"),
    ] = None,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Output as JSON"),
    ] = False,
) -> None:
    """
    Show details for a specific database.

    Displays database metadata and validates that the path exists.

    Examples:

        nvd state database show blast
        nvd state database show blast core-nt_2024-12
        nvd state database show blast --json
    """
    db_path = ensure_db_exists(json_output)

    # Validate db_type
    valid_types = ["blast", "stat", "gottcha2"]
    if db_type not in valid_types:
        if json_output:
            _output_json(
                {
                    "error": f"Invalid database type: {db_type}",
                    "valid_types": valid_types,
                },
            )
            raise typer.Exit(1)
        error(
            f"Invalid database type: {db_type!r}. Must be one of: {', '.join(valid_types)}",
        )

    db = state.get_database_by_version(db_type, version)  # type: ignore[arg-type]

    if db is None:
        if json_output:
            _output_json(
                {"error": "Database not found", "db_type": db_type, "version": version},
            )
            raise typer.Exit(1)
        if version:
            error(f"Database not found: {db_type} version {version}")
        else:
            error(f"No {db_type} database registered")

    path_exists = Path(db.path).exists()

    if json_output:
        _output_json(
            {
                "db_type": db.db_type,
                "version": db.version,
                "path": db.path,
                "registered_at": db.registered_at,
                "path_exists": path_exists,
            },
        )
        return

    console.print()
    console.print(f"[bold]Database:[/bold] {db.db_type}")
    console.print(f"  Version: {db.version}")
    console.print(f"  Path: {db.path}")
    console.print(f"  Registered: {db.registered_at}")
    console.print()

    if path_exists:
        success("Path exists and is accessible")
    else:
        warning("Path does not exist or is not accessible")

    console.print()


@database_app.command("register")
def database_register(
    db_type: Annotated[
        str,
        typer.Argument(help="Database type (blast, stat, gottcha2)"),
    ],
    path: Annotated[
        Path,
        typer.Argument(help="Path to the database"),
    ],
    version: Annotated[
        str,
        typer.Option("--version", "-v", help="Version label for this database"),
    ],
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Output as JSON"),
    ] = False,
) -> None:
    """
    Register a reference database.

    Associates a database path with a version label for tracking and
    auto-resolution during pipeline runs.

    Examples:

        nvd state database register blast /data/blast/core-nt -v core-nt_2025-01
        nvd state database register gottcha2 /data/gottcha2/RefSeq-r220 -v RefSeq-r220
        nvd state database register stat /data/stat/index.k31 -v 2024-12
    """
    # Validate db_type
    valid_types = ["blast", "stat", "gottcha2"]
    if db_type not in valid_types:
        if json_output:
            _output_json(
                {
                    "error": f"Invalid database type: {db_type}",
                    "valid_types": valid_types,
                },
            )
            raise typer.Exit(1)
        error(
            f"Invalid database type: {db_type!r}. Must be one of: {', '.join(valid_types)}",
        )
        raise typer.Exit(1)

    # Warn if path doesn't exist (but allow registration)
    if not path.exists():
        if not json_output:
            warning(f"Path does not exist: {path}")
            warning(
                "Registering anyway - ensure path is correct before running pipeline",
            )

    # Register the database (path is canonicalized by register_database)
    db = state.register_database(
        db_type=db_type,  # type: ignore[arg-type]
        version=version,
        path=str(path),
    )

    if json_output:
        _output_json(
            {
                "registered": True,
                "db_type": db.db_type,
                "version": db.version,
                "path": db.path,
                "registered_at": db.registered_at,
            },
        )
        return

    success(f"Registered {db.db_type} database: {db.version}")
    console.print(f"  Path: {db.path}")
    console.print()


@database_app.command("unregister")
def database_unregister(
    db_type: Annotated[
        str,
        typer.Argument(help="Database type (blast, stat, gottcha2)"),
    ],
    version: Annotated[
        str,
        typer.Argument(help="Version to unregister"),
    ],
    yes: Annotated[
        bool,
        typer.Option("--yes", "-y", help="Skip confirmation prompt"),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Output as JSON"),
    ] = False,
) -> None:
    """
    Unregister a reference database.

    Removes a database registration from the state. This does not delete
    the actual database files.

    Examples:

        nvd state database unregister blast core-nt_2024-01
        nvd state database unregister blast core-nt_2024-01 --yes
    """
    # Validate db_type
    valid_types = ["blast", "stat", "gottcha2"]
    if db_type not in valid_types:
        if json_output:
            _output_json(
                {
                    "error": f"Invalid database type: {db_type}",
                    "valid_types": valid_types,
                },
            )
            raise typer.Exit(1)
        error(
            f"Invalid database type: {db_type!r}. Must be one of: {', '.join(valid_types)}",
        )
        raise typer.Exit(1)

    # Check if database exists
    db = state.get_database_by_version(db_type, version)  # type: ignore[arg-type]
    if db is None:
        if json_output:
            _output_json(
                {"error": "Database not found", "db_type": db_type, "version": version},
            )
            raise typer.Exit(1)
        error(f"Database not found: {db_type} version {version}")

    # Confirm unless --yes
    if not yes and not json_output:
        console.print("\n[bold]Database to unregister:[/bold]")
        console.print(f"  Type: {db.db_type}")
        console.print(f"  Version: {db.version}")
        console.print(f"  Path: {db.path}")
        console.print()
        warning("This will remove the registration (not the actual files).")

        confirm = typer.confirm("Proceed with unregistration?")
        if not confirm:
            info("Cancelled")
            raise typer.Exit(0)

    # Delete the registration
    with connect() as conn:
        conn.execute(
            "DELETE FROM databases WHERE db_type = ? AND version = ?",
            (db_type, version),
        )
        conn.commit()

    if json_output:
        _output_json({"unregistered": True, "db_type": db_type, "version": version})
        return

    success(f"Unregistered {db_type} database: {version}")
    console.print()


@state_app.command("taxonomy")
def state_taxonomy(
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Show taxonomy cache status.

    Displays information about the local taxonomy cache, including
    file locations and sizes.

    Examples:

        # Show taxonomy status
        nvd state taxonomy

        # Output as JSON
        nvd state taxonomy --json
    """
    taxdump_dir = get_taxdump_dir()
    taxonomy_db = taxdump_dir / "taxonomy.sqlite"

    taxdump_exists = taxdump_dir.exists()
    taxonomy_db_exists = taxonomy_db.exists()

    # Check for .dmp files
    dmp_files: dict[str, int] = {}
    if taxdump_exists:
        for dmp_file in taxdump_dir.glob("*.dmp"):
            dmp_files[dmp_file.name] = dmp_file.stat().st_size

    taxonomy_db_size = taxonomy_db.stat().st_size if taxonomy_db_exists else 0

    if json_output:
        _output_json(
            {
                "taxdump_dir": str(taxdump_dir),
                "taxdump_exists": taxdump_exists,
                "taxonomy_db": {
                    "path": str(taxonomy_db),
                    "exists": taxonomy_db_exists,
                    "size_bytes": taxonomy_db_size,
                },
                "dmp_files": dmp_files,
            },
        )
        return

    console.print()
    console.print(f"[bold]Taxonomy Cache:[/bold] {taxdump_dir}")
    console.print()

    if not taxdump_exists:
        warning("Taxonomy cache directory does not exist")
        info("The taxonomy cache will be created when you run the pipeline")
        console.print()
        return

    # Taxonomy database
    if taxonomy_db_exists:
        success(f"taxonomy.sqlite: {_format_size(taxonomy_db_size)}")
    else:
        info("taxonomy.sqlite: not found")

    # DMP files
    if dmp_files:
        console.print()
        console.print("[bold]Taxonomy dump files:[/bold]")
        for name, size in sorted(dmp_files.items()):
            console.print(f"  {name}: {_format_size(size)}")
    else:
        console.print()
        info("No .dmp files found")

    console.print()


def _parse_duration(value: str) -> int:
    """
    Parse a duration string into days.

    Accepts:
        - "30d" — 30 days
        - "12w" — 12 weeks (84 days)
        - "6m" — 6 months (180 days, approximate)

    Returns:
        Number of days

    Raises:
        typer.BadParameter: If the string cannot be parsed
    """
    import re

    match = re.match(r"^(\d+)([dwm])$", value.lower())
    if not match:
        raise typer.BadParameter(
            f"Invalid duration format: {value!r}. Use format like '30d', '12w', or '6m'",
        )

    amount = int(match.group(1))
    unit = match.group(2)

    if unit == "d":
        return amount
    if unit == "w":
        return amount * 7
    if unit == "m":
        return amount * 30  # Approximate
    raise typer.BadParameter(f"Unknown duration unit: {unit}")


@state_app.command("prune")
def state_prune(
    older_than: Annotated[
        str,
        typer.Option(
            "--older-than",
            "-o",
            help="Remove runs older than this duration (e.g., 90d, 12w, 6m)",
        ),
    ],
    status: Annotated[
        str | None,
        typer.Option(
            "--status",
            "-s",
            help="Only prune runs with this status (completed, failed)",
        ),
    ] = None,
    include_running: Annotated[
        bool,
        typer.Option(
            "--include-running",
            help="Include runs with 'running' status (normally skipped)",
        ),
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            "-n",
            help="Show what would be deleted without deleting",
        ),
    ] = False,
    yes: Annotated[
        bool,
        typer.Option(
            "--yes",
            "-y",
            help="Skip confirmation prompt",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Remove old runs and their associated data.

    Prunes runs older than the specified duration, along with their
    associated samples, uploads, and taxonomy versions.

    Runs with 'running' status are skipped by default (use --include-running
    to include them).

    Examples:

        # Show what would be deleted (dry run)
        nvd state prune --older-than 90d --dry-run

        # Prune runs older than 90 days
        nvd state prune --older-than 90d

        # Prune only failed runs older than 30 days
        nvd state prune --older-than 30d --status failed

        # Prune without confirmation
        nvd state prune --older-than 90d --yes
    """
    from datetime import timedelta

    db_path = ensure_db_exists(json_output)

    # Parse duration
    try:
        days = _parse_duration(older_than)
    except typer.BadParameter as e:
        if json_output:
            _output_json({"error": str(e)})
            raise typer.Exit(1)
        error(str(e))

    # Validate status filter
    valid_statuses = ["completed", "failed"]
    if status is not None and status not in valid_statuses:
        if not include_running or status != "running":
            if json_output:
                _output_json(
                    {
                        "error": f"Invalid status: {status}",
                        "valid_statuses": valid_statuses,
                    },
                )
                raise typer.Exit(1)
            error(
                f"Invalid status: {status!r}. Must be one of: {', '.join(valid_statuses)}",
            )
            raise typer.Exit(1)

    # Calculate cutoff date
    cutoff = datetime.now() - timedelta(days=days)
    cutoff_str = cutoff.strftime("%Y-%m-%d %H:%M:%S")

    with connect() as conn:
        # Build query for runs to prune
        query = "SELECT * FROM runs WHERE started_at < ?"
        params: list[str] = [cutoff_str]

        if not include_running:
            query += " AND status != 'running'"

        if status is not None:
            query += " AND status = ?"
            params.append(status)

        query += " ORDER BY started_at"

        runs_to_prune = conn.execute(query, params).fetchall()

        if not runs_to_prune:
            if json_output:
                _output_json(
                    {
                        "message": "No runs match the criteria",
                        "cutoff_date": cutoff_str,
                        "status_filter": status,
                        "include_running": include_running,
                    },
                )
                return
            info("No runs match the criteria")
            console.print(f"[dim]Cutoff date: {cutoff_str}[/dim]")
            return

        # Gather statistics for each run
        prune_details = []
        total_samples = 0
        total_uploads = 0
        total_taxonomy = 0
        total_observations = 0

        for run in runs_to_prune:
            run_id = run["run_id"]
            sample_set_id = run["sample_set_id"]

            # Count samples for this run
            sample_count = conn.execute(
                "SELECT COUNT(*) FROM processed_samples WHERE run_id = ?",
                (run_id,),
            ).fetchone()[0]

            # Count uploads for samples in this run's sample_set
            upload_count = conn.execute(
                "SELECT COUNT(*) FROM uploads WHERE sample_set_id = ?",
                (sample_set_id,),
            ).fetchone()[0]

            # Count hit observations for this sample_set (from parquet files)
            observation_count = count_sample_set_observations(sample_set_id)

            # Check for taxonomy version
            has_taxonomy = (
                conn.execute(
                    "SELECT 1 FROM taxonomy_versions WHERE run_id = ?",
                    (run_id,),
                ).fetchone()
                is not None
            )

            prune_details.append(
                {
                    "run_id": run_id,
                    "started_at": run["started_at"],
                    "status": run["status"],
                    "sample_count": sample_count,
                    "upload_count": upload_count,
                    "observation_count": observation_count,
                    "has_taxonomy": has_taxonomy,
                },
            )

            total_samples += sample_count
            total_uploads += upload_count
            total_observations += observation_count
            if has_taxonomy:
                total_taxonomy += 1

        # JSON output
        if json_output:
            _output_json(
                {
                    "dry_run": dry_run,
                    "cutoff_date": cutoff_str,
                    "status_filter": status,
                    "include_running": include_running,
                    "runs": prune_details,
                    "totals": {
                        "runs": len(prune_details),
                        "samples": total_samples,
                        "uploads": total_uploads,
                        "taxonomy_versions": total_taxonomy,
                        "hit_observations": total_observations,
                    },
                },
            )
            if dry_run:
                return
            # For JSON mode with --yes, proceed to delete
            if not yes:
                raise typer.Exit(0)

        # Human-readable output
        if not json_output:
            console.print()
            console.print(
                f"[bold]Runs older than {older_than}[/bold] (before {cutoff.strftime('%Y-%m-%d')}):",
            )
            console.print()

            table = Table(show_header=True, header_style="bold cyan")
            table.add_column("Run ID", style="cyan")
            table.add_column("Started")
            table.add_column("Status")
            table.add_column("Samples", justify="right")
            table.add_column("Uploads", justify="right")

            for detail in prune_details:
                if detail["status"] == "completed":
                    status_str = "[green]completed[/green]"
                elif detail["status"] == "failed":
                    status_str = "[red]failed[/red]"
                else:
                    status_str = "[yellow]running[/yellow]"

                table.add_row(
                    detail["run_id"],
                    detail["started_at"],
                    status_str,
                    str(detail["sample_count"]),
                    str(detail["upload_count"]),
                )

            console.print(table)
            console.print()

            console.print("[bold]Summary:[/bold]")
            console.print(f"  {len(prune_details)} run(s)")
            console.print(f"  {total_samples} processed sample(s)")
            console.print(f"  {total_uploads} upload(s)")
            console.print(f"  {total_taxonomy} taxonomy version(s)")
            console.print(f"  {total_observations} hit observation(s)")
            console.print()

            if dry_run:
                info("Dry run — no changes made")
                return

            # Confirmation
            if not yes:
                warning("This action cannot be undone.")
                confirm = typer.confirm("Proceed?", default=False)
                if not confirm:
                    info("Aborted")
                    raise typer.Exit(0)

        # Perform the deletion
        sample_set_ids_to_prune = []
        for detail in prune_details:
            run_id = detail["run_id"]

            # Get sample_set_id for this run
            run_row = conn.execute(
                "SELECT sample_set_id FROM runs WHERE run_id = ?",
                (run_id,),
            ).fetchone()
            sample_set_id = run_row["sample_set_id"]
            sample_set_ids_to_prune.append(sample_set_id)

            # Delete uploads for this sample_set
            conn.execute(
                "DELETE FROM uploads WHERE sample_set_id = ?",
                (sample_set_id,),
            )

            # Delete processed samples for this run
            conn.execute(
                "DELETE FROM processed_samples WHERE run_id = ?",
                (run_id,),
            )

            # Delete taxonomy version for this run
            conn.execute(
                "DELETE FROM taxonomy_versions WHERE run_id = ?",
                (run_id,),
            )

            # Delete the run itself
            conn.execute(
                "DELETE FROM runs WHERE run_id = ?",
                (run_id,),
            )

        conn.commit()

    # Delete parquet hit files for pruned sample sets
    total_hit_files_deleted = 0
    total_months_rewritten = 0
    for sample_set_id in sample_set_ids_to_prune:
        result = delete_sample_set_hits(sample_set_id)
        total_hit_files_deleted += result.uncompacted_files_deleted
        total_months_rewritten += result.compacted_months_rewritten

    if not json_output:
        hit_summary = f"{total_hit_files_deleted} hit file(s)"
        if total_months_rewritten > 0:
            hit_summary += f", {total_months_rewritten} month(s) rewritten"
        success(
            f"Pruned {len(prune_details)} run(s), {total_samples} sample(s), "
            f"{total_uploads} upload(s), {total_taxonomy} taxonomy version(s), "
            f"{hit_summary}",
        )


def _get_table_counts(conn) -> dict[str, int]:
    """Get counts for all tables in the state database."""
    return {
        "runs": conn.execute("SELECT COUNT(*) FROM runs").fetchone()[0],
        "processed_samples": conn.execute(
            "SELECT COUNT(*) FROM processed_samples",
        ).fetchone()[0],
        "uploads": conn.execute("SELECT COUNT(*) FROM uploads").fetchone()[0],
        "databases": conn.execute("SELECT COUNT(*) FROM databases").fetchone()[0],
        "presets": conn.execute("SELECT COUNT(*) FROM presets").fetchone()[0],
        "taxonomy_versions": conn.execute(
            "SELECT COUNT(*) FROM taxonomy_versions",
        ).fetchone()[0],
        "sra_cache": conn.execute("SELECT COUNT(*) FROM sra_cache").fetchone()[0],
        "sample_locks": conn.execute(
            "SELECT COUNT(*) FROM sample_locks",
        ).fetchone()[0],
    }


def _build_export_manifest(conn) -> dict:
    """Build the manifest for a state export."""
    import socket

    schema_version = conn.execute("PRAGMA user_version").fetchone()[0]
    counts = _get_table_counts(conn)

    return {
        "schema_version": schema_version,
        "exported_at": datetime.now().isoformat(),
        "nvd_version": __version__,
        "hostname": socket.gethostname(),
        "counts": counts,
    }


@state_app.command("export")
def state_export(
    path: Annotated[
        Path | None,
        typer.Argument(
            help="Output path for the archive (auto-appends .nvd-state if needed)",
        ),
    ] = None,
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            "-f",
            help="Overwrite existing file",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output manifest to stdout instead of creating archive (inspection mode)",
        ),
    ] = False,
) -> None:
    """
    Export state database to a portable archive.

    Creates a .nvd-state archive containing the database and a manifest
    with metadata about the export. Use 'nvd state import' to restore.

    The --json flag outputs the manifest to stdout without creating an
    archive, useful for inspecting what would be exported.

    Examples:

        # Export to an archive
        nvd state export backup.nvd-state

        # Export with auto-appended extension
        nvd state export backup
        # Creates: backup.nvd-state

        # Inspect what would be exported (no file created)
        nvd state export --json

        # Overwrite existing archive
        nvd state export backup.nvd-state --force
    """
    import shutil
    import tarfile
    import tempfile

    # Validate arguments
    if json_output and path is not None:
        error("Cannot specify both --json and a path")
        info("Use --json alone to inspect, or provide a path to create an archive")
        raise typer.Exit(1)

    if not json_output and path is None:
        error("Path is required unless using --json")
        info("Usage: nvd state export <path> or nvd state export --json")
        raise typer.Exit(1)

    db_path = ensure_db_exists(json_output)

    # Build manifest
    with connect() as conn:
        manifest = _build_export_manifest(conn)

    # JSON inspection mode
    if json_output:
        _output_json(manifest)
        return

    # Ensure path has .nvd-state extension
    assert path is not None  # Already validated above
    output_path = path
    if not str(output_path).endswith(".nvd-state"):
        output_path = Path(str(output_path) + ".nvd-state")

    # Check if output exists
    if output_path.exists() and not force:
        error(f"File already exists: {output_path}")
        info("Use --force to overwrite")
        raise typer.Exit(1)

    # Create the archive
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Copy database to temp directory
        tmp_db = tmpdir_path / "state.sqlite"
        shutil.copy2(db_path, tmp_db)

        # Write manifest
        tmp_manifest = tmpdir_path / "manifest.json"
        tmp_manifest.write_text(json.dumps(manifest, indent=2, default=str))

        # Create tar.gz archive
        with tarfile.open(output_path, "w:gz") as tar:
            tar.add(tmp_db, arcname="state.sqlite")
            tar.add(tmp_manifest, arcname="manifest.json")

    archive_size = output_path.stat().st_size

    success(f"Exported state to: {output_path}")
    console.print(f"  Archive size: {_format_size(archive_size)}")
    console.print(f"  Schema version: {manifest['schema_version']}")
    console.print(f"  NVD version: {manifest['nvd_version']}")
    console.print()
    console.print("[bold]Contents:[/bold]")
    console.print(f"  {manifest['counts']['runs']} run(s)")
    console.print(f"  {manifest['counts']['processed_samples']} processed sample(s)")
    console.print(f"  {manifest['counts']['uploads']} upload(s)")
    console.print(f"  {manifest['counts']['databases']} database(s)")
    console.print(f"  {manifest['counts']['presets']} preset(s)")
    console.print(f"  {manifest['counts']['taxonomy_versions']} taxonomy version(s)")
    console.print(f"  {manifest['counts']['sra_cache']} SRA cache entries")
    console.print()


# Table definitions for import conflict detection
# Each tuple: (table_name, primary_key_columns, display_columns)
_IMPORT_TABLES = [
    ("runs", ["run_id"], ["started_at", "status"]),
    ("processed_samples", ["sample_id", "sample_set_id"], ["run_id", "status"]),
    (
        "uploads",
        ["sample_id", "sample_set_id", "upload_type", "upload_target"],
        ["uploaded_at"],
    ),
    ("databases", ["db_type", "version"], ["path", "registered_at"]),
    ("presets", ["name"], ["created_at"]),
    ("taxonomy_versions", ["run_id"], ["file_hash", "recorded_at"]),
    ("sra_cache", ["srr_accession"], ["download_path", "downloaded_at"]),
    ("sample_locks", ["sample_id"], ["run_id", "locked_at", "expires_at"]),
]


def _validate_archive(archive_path: Path) -> tuple[dict, Path]:
    """
    Validate and extract an nvd-state archive.

    Returns:
        Tuple of (manifest dict, path to extracted temp directory)

    Raises:
        typer.Exit: If archive is invalid
    """
    import tarfile
    import tempfile

    if not archive_path.exists():
        error(f"Archive not found: {archive_path}")

    # Create temp directory (caller is responsible for cleanup)
    tmpdir = tempfile.mkdtemp(prefix="nvd-import-")
    tmpdir_path = Path(tmpdir)

    try:
        with tarfile.open(archive_path, "r:gz") as tar:
            # Check for required files
            names = tar.getnames()
            if "state.sqlite" not in names:
                error("Invalid archive: missing state.sqlite")
            if "manifest.json" not in names:
                error("Invalid archive: missing manifest.json")

    except tarfile.TarError as e:
        error(f"Failed to read archive: {e}")

    # Parse manifest
    manifest_path = tmpdir_path / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (json.JSONDecodeError, OSError) as e:
        error(f"Failed to read manifest: {e}")

    # Validate manifest has required fields
    required_fields = ["schema_version", "exported_at", "nvd_version", "counts"]
    for field in required_fields:
        if field not in manifest:
            error(f"Invalid manifest: missing '{field}'")

    if json_output:
        _output_json(
            {
                "error": "Merge aborted due to conflicts",
                "conflicts": conflicts,
                "conflict_count": sum(len(v) for v in conflicts.values()),
            },
        )
        return

    total = sum(len(v) for v in conflicts.values())
    # Don't use error() here as it calls sys.exit()
    console.print(
        f"[red]✗ Error:[/red] Merge aborted due to {total} conflict(s) "
        f"in {len(conflicts)} table(s)",
    )
    console.print()

    for table_name, table_conflicts in conflicts.items():
        console.print(
            f"[bold]{table_name}[/bold] ({len(table_conflicts)} conflict(s)):",
        )

        for conflict in table_conflicts:
            # Format primary key
            pk_str = ", ".join(f"{k}={v}" for k, v in conflict["pk"].items())
            console.print(f"  [cyan]{pk_str}[/cyan]")

            # Show local vs incoming
            local_str = ", ".join(f"{k}={v}" for k, v in conflict["local"].items())
            incoming_str = ", ".join(
                f"{k}={v}" for k, v in conflict["incoming"].items()
            )
            console.print(f"    Local:    {local_str}")
            console.print(f"    Incoming: {incoming_str}")

        console.print()

    console.print("[dim]Resolve conflicts manually, then retry import.[/dim]")


def _do_merge_import(conn, backup_db_path: Path) -> dict[str, int]:
    """
    Perform merge import from backup database.

    Assumes conflicts have already been checked.

    Returns:
        Dict mapping table names to count of records imported.
    """
    # Validate path is safe for SQL (ATTACH DATABASE cannot use parameterized queries)
    validated_path = _validate_path_for_sql(backup_db_path)

    imported: dict[str, int] = {}

    conn.execute(f"ATTACH DATABASE '{validated_path}' AS backup")

    try:
        for table_name, pk_cols, _ in _IMPORT_TABLES:
            # Build WHERE NOT EXISTS to skip records that already exist
            pk_match = " AND ".join(
                f"main.{table_name}.{col} = backup.{table_name}.{col}"
                for col in pk_cols
            )

            # Get column names for this table
            cols_info = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
            col_names = [col["name"] for col in cols_info]
            cols_str = ", ".join(col_names)

            # Insert records that don't exist in main
            query = f"""
                INSERT INTO main.{table_name} ({cols_str})
                SELECT {cols_str} FROM backup.{table_name}
                WHERE NOT EXISTS (
                    SELECT 1 FROM main.{table_name}
                    WHERE {pk_match}
                )
            """

            cursor = conn.execute(query)
            imported[table_name] = cursor.rowcount

        conn.commit()

    finally:
        conn.execute("DETACH DATABASE backup")

    return imported


def _do_replace_import(backup_db_path: Path, current_db_path: Path) -> None:
    """
    Perform replace import by swapping database files.

    Uses atomic rename for safety.
    """
    import shutil

    # Backup current DB just in case (to .bak)
    backup_path = current_db_path.with_suffix(".sqlite.bak")
    if current_db_path.exists():
        shutil.copy2(current_db_path, backup_path)

    try:
        # Copy new DB over current
        shutil.copy2(backup_db_path, current_db_path)

        # Remove backup on success
        if backup_path.exists():
            backup_path.unlink()

    except OSError:
        # Restore from backup on failure
        if backup_path.exists():
            shutil.copy2(backup_path, current_db_path)
            backup_path.unlink()
        raise


@state_app.command("import")
def state_import(
    path: Annotated[
        Path,
        typer.Argument(help="Path to the .nvd-state archive to import"),
    ],
    replace: Annotated[
        bool,
        typer.Option(
            "--replace",
            help="Replace current state entirely (default: merge)",
        ),
    ] = False,
    yes: Annotated[
        bool,
        typer.Option(
            "--yes",
            "-y",
            help="Skip confirmation prompt",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Import state from an archive.

    By default, performs a merge import that adds records from the archive
    that don't exist locally. If any primary key conflicts are found, the
    import is aborted and a conflict report is shown.

    Use --replace for a full restore that discards current state.

    Examples:

        # Merge import (default) — adds non-conflicting records
        nvd state import backup.nvd-state

        # Full restore — replaces current state entirely
        nvd state import backup.nvd-state --replace

        # Skip confirmation prompt
        nvd state import backup.nvd-state --yes
        nvd state import backup.nvd-state --replace --yes

        # JSON output (useful for scripting)
        nvd state import backup.nvd-state --json
    """
    import shutil

    # Validate archive and extract to temp directory
    manifest, tmpdir_path = _validate_archive(path)

    try:
        backup_db_path = tmpdir_path / "state.sqlite"

        # Check schema version compatibility
        incoming_version = manifest.get("schema_version")
        if incoming_version != EXPECTED_VERSION:
            if json_output:
                _output_json(
                    {
                        "error": "Schema version mismatch",
                        "incoming_version": incoming_version,
                        "current_version": EXPECTED_VERSION,
                    },
                )
            else:
                # Don't use error() here as it calls sys.exit()
                console.print("[red]✗ Error:[/red] Schema version mismatch")
                console.print(f"  Archive version: {incoming_version}")
                console.print(f"  Current version: {EXPECTED_VERSION}")
                console.print()
                if incoming_version is not None and incoming_version < EXPECTED_VERSION:
                    info(
                        "The archive is from an older NVD version. "
                        "Export from a newer version to upgrade.",
                    )
                else:
                    info(
                        "The archive is from a newer NVD version. "
                        "Upgrade NVD to import this archive.",
                    )
            raise typer.Exit(1)

        db_path = get_state_db_path()

        # Ensure current database exists (for merge mode)
        if not replace and not db_path.exists():
            if json_output:
                _output_json(
                    {
                        "error": "Database not found",
                        "path": str(db_path),
                        "hint": "Use --replace to create from archive, or run 'nvd state init' first",
                    },
                )
            else:
                error(f"Database not found: {db_path}")
                info(
                    "Use --replace to create from archive, or run 'nvd state init' first",
                )
            raise typer.Exit(1)

        if replace:
            # Replace mode: full restore
            incoming_counts = manifest.get("counts", {})

            # Get current counts if DB exists
            current_counts: dict[str, int] = {}
            if db_path.exists():
                with connect() as conn:
                    current_counts = _get_table_counts(conn)

            if json_output:
                if not yes:
                    _output_json(
                        {
                            "mode": "replace",
                            "current_counts": current_counts,
                            "incoming_counts": incoming_counts,
                            "requires_confirmation": True,
                        },
                    )
                    raise typer.Exit(0)
            else:
                console.print()
                console.print(
                    "[bold]Replace mode[/bold] — current state will be discarded",
                )
                console.print()

                # Show comparison table
                table = Table(show_header=True, header_style="bold cyan")
                table.add_column("Table")
                table.add_column("Current", justify="right")
                table.add_column("Incoming", justify="right")

                for key in [
                    "runs",
                    "processed_samples",
                    "uploads",
                    "databases",
                    "presets",
                    "taxonomy_versions",
                    "sra_cache",
                ]:
                    table.add_row(
                        key,
                        str(current_counts.get(key, 0)),
                        str(incoming_counts.get(key, 0)),
                    )

                console.print(table)
                console.print()

                if not yes:
                    warning("This will discard all current state data.")
                    confirm = typer.confirm("Proceed?", default=False)
                    if not confirm:
                        info("Aborted")
                        raise typer.Exit(0)

            # Perform replace
            _do_replace_import(backup_db_path, db_path)

            if json_output:
                _output_json(
                    {
                        "mode": "replace",
                        "success": True,
                        "imported_counts": incoming_counts,
                    },
                )
            else:
                success("State replaced from archive")
                console.print(f"  Source: {path}")
                console.print(
                    f"  Exported at: {manifest.get('exported_at', 'unknown')}",
                )
                console.print(f"  From host: {manifest.get('hostname', 'unknown')}")
                console.print()

        else:
            # Merge mode: check for conflicts first
            with connect() as conn:
                conflicts = _find_conflicts(conn, backup_db_path)

                if conflicts:
                    _print_conflicts(conflicts, json_output)
                    raise typer.Exit(1)

                # No conflicts — show what will be imported
                incoming_counts = manifest.get("counts", {})
                current_counts = _get_table_counts(conn)

                if not json_output:
                    console.print()
                    console.print(
                        "[bold]Merge mode[/bold] — adding non-conflicting records",
                    )
                    console.print()

                    # Show what will be added
                    table = Table(show_header=True, header_style="bold cyan")
                    table.add_column("Table")
                    table.add_column("Current", justify="right")
                    table.add_column("In Archive", justify="right")

                    for key in [
                        "runs",
                        "processed_samples",
                        "uploads",
                        "databases",
                        "presets",
                        "taxonomy_versions",
                        "sra_cache",
                    ]:
                        table.add_row(
                            key,
                            str(current_counts.get(key, 0)),
                            str(incoming_counts.get(key, 0)),
                        )

                    console.print(table)
                    console.print()

                    if not yes:
                        confirm = typer.confirm("Proceed with merge?", default=True)
                        if not confirm:
                            info("Aborted")
                            raise typer.Exit(0)

                # Perform merge
                imported = _do_merge_import(conn, backup_db_path)

            if json_output:
                _output_json(
                    {
                        "mode": "merge",
                        "success": True,
                        "imported": imported,
                        "total_imported": sum(imported.values()),
                    },
                )
            else:
                total = sum(imported.values())
                if total > 0:
                    success(f"Merged {total} record(s) from archive")
                    for table_name, count in imported.items():
                        if count > 0:
                            console.print(f"  {table_name}: {count}")
                else:
                    info("No new records to import (all records already exist)")
                console.print()

    finally:
        # Clean up temp directory
        shutil.rmtree(tmpdir_path, ignore_errors=True)


@state_app.command("reset")
def state_reset(
    yes_i_really_mean_it: Annotated[
        bool,
        typer.Option(
            "--yes-i-really-mean-it",
            help="Confirm deletion without interactive prompt",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Clear all state data.

    Deletes all records from the state database while preserving the schema.
    This is a destructive operation that cannot be undone.

    Consider using 'nvd state export' to create a backup before resetting.

    Examples:

        # Show what would be deleted (requires confirmation)
        nvd state reset

        # Reset without interactive prompt
        nvd state reset --yes-i-really-mean-it

        # JSON output for scripting
        nvd state reset --json --yes-i-really-mean-it
    """
    db_path = ensure_db_exists(json_output)

    with connect() as conn:
        counts = _get_table_counts(conn)

    total = sum(counts.values())

    if total == 0:
        if json_output:
            _output_json({"message": "Database is already empty", "counts": counts})
            return
        info("Database is already empty")
        return

    # Show what will be deleted
    if json_output:
        if not yes_i_really_mean_it:
            _output_json(
                {
                    "counts": counts,
                    "total": total,
                    "requires_confirmation": True,
                    "hint": "Use --yes-i-really-mean-it to confirm",
                },
            )
            raise typer.Exit(0)
    else:
        console.print()
        console.print("[bold red]Warning:[/bold red] This will delete all state data!")
        console.print()

        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("Table")
        table.add_column("Records", justify="right")

        for table_name, count in counts.items():
            if count > 0:
                table.add_row(table_name, str(count))

        console.print(table)
        console.print()
        console.print(f"[bold]Total: {total} record(s) will be deleted[/bold]")
        console.print()

        if not yes_i_really_mean_it:
            warning("This action cannot be undone.")
            info("Consider running 'nvd state export' first to create a backup.")
            console.print()
            confirm = typer.confirm(
                "Are you sure you want to delete all state data?",
                default=False,
            )
            if not confirm:
                info("Aborted")
                raise typer.Exit(0)

    # Perform deletion
    with connect() as conn:
        deleted = {}
        for table_name, _, _ in _IMPORT_TABLES:
            cursor = conn.execute(f"DELETE FROM {table_name}")  # noqa: S608
            deleted[table_name] = cursor.rowcount

        conn.commit()

    if json_output:
        _output_json(
            {
                "success": True,
                "deleted": deleted,
                "total_deleted": sum(deleted.values()),
            },
        )
    else:
        total_deleted = sum(deleted.values())
        success(f"Deleted {total_deleted} record(s)")
        for table_name, count in deleted.items():
            if count > 0:
                console.print(f"  {table_name}: {count}")
        console.print()


# =============================================================================
# Sample Lock Commands
# =============================================================================


def _format_time_remaining(expires_at: str) -> str:
    """Format time remaining until lock expires."""
    from datetime import timezone

    try:
        expires = datetime.fromisoformat(expires_at.replace("Z", "+00:00"))
        if expires.tzinfo is None:
            expires = expires.replace(tzinfo=timezone.utc)
        now = datetime.now(timezone.utc)
        remaining = expires - now

        if remaining.total_seconds() <= 0:
            return "[red]expired[/red]"

        hours = int(remaining.total_seconds() // 3600)
        minutes = int((remaining.total_seconds() % 3600) // 60)

        if hours > 24:
            days = hours // 24
            hours = hours % 24
            return f"{days}d {hours}h"
        elif hours > 0:
            return f"{hours}h {minutes}m"
        else:
            return f"{minutes}m"
    except (ValueError, TypeError):
        return "unknown"


@state_app.command("locks")
def state_locks(
    run_id: Annotated[
        str | None,
        typer.Option(
            "--run-id",
            "-r",
            help="Filter by run ID",
        ),
    ] = None,
    include_expired: Annotated[
        bool,
        typer.Option(
            "--expired",
            "-e",
            help="Include expired locks",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    List current sample locks.

    Shows samples that are currently locked for processing, preventing
    duplicate work across concurrent runs.

    Examples:

        # List all active locks
        nvd state locks

        # List locks for a specific run
        nvd state locks --run-id happy_volta

        # Include expired locks
        nvd state locks --expired

        # JSON output for scripting
        nvd state locks --json
    """
    ensure_db_exists(json_output)

    if run_id:
        locks = state.get_locks_for_run(run_id, include_expired=include_expired)
    else:
        locks = state.list_locks(include_expired=include_expired)

    if json_output:
        _output_json(
            [
                {
                    "sample_id": lock.sample_id,
                    "run_id": lock.run_id,
                    "hostname": lock.hostname,
                    "username": lock.username,
                    "locked_at": lock.locked_at,
                    "expires_at": lock.expires_at,
                    "is_expired": lock.is_expired,
                }
                for lock in locks
            ]
        )
        return

    if not locks:
        if run_id:
            info(f"No locks found for run '{run_id}'")
        else:
            info("No active locks")
        return

    console.print()
    table = Table(show_header=True, header_style="bold cyan")
    table.add_column("Sample ID")
    table.add_column("Run ID")
    table.add_column("Host")
    table.add_column("User")
    table.add_column("Locked At")
    table.add_column("Expires In")

    for lock in locks:
        time_remaining = _format_time_remaining(lock.expires_at)
        table.add_row(
            lock.sample_id,
            lock.run_id,
            lock.hostname,
            lock.username,
            lock.locked_at[:19],  # Trim to datetime without timezone
            time_remaining,
        )

    console.print(table)
    console.print()
    console.print(f"[bold]Total: {len(locks)} lock(s)[/bold]")
    console.print()


@state_app.command("unlock")
def state_unlock(
    sample_id: Annotated[
        str | None,
        typer.Option(
            "--sample-id",
            "-s",
            help="Unlock specific sample",
        ),
    ] = None,
    run_id: Annotated[
        str | None,
        typer.Option(
            "--run-id",
            "-r",
            help="Unlock all samples for a run",
        ),
    ] = None,
    force: Annotated[
        bool,
        typer.Option(
            "--force",
            "-f",
            help="Skip confirmation prompt",
        ),
    ] = False,
    steal: Annotated[
        bool,
        typer.Option(
            "--steal",
            help="Steal locks even if held by another run (requires --force)",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Manually release sample locks.

    Use this to unlock samples that are blocking a run due to a crashed
    or abandoned previous run.

    Examples:

        # Unlock a specific sample (with confirmation)
        nvd state unlock --sample-id sample_A

        # Unlock all samples for a run
        nvd state unlock --run-id happy_volta

        # Force unlock without confirmation
        nvd state unlock --run-id happy_volta --force

        # Steal a lock held by another run
        nvd state unlock --sample-id sample_A --steal --force
    """
    ensure_db_exists(json_output)

    if not sample_id and not run_id:
        if json_output:
            _output_json({"error": "Must specify --sample-id or --run-id"})
            raise typer.Exit(1)
        error("Must specify --sample-id or --run-id")
        raise typer.Exit(1)

    if steal and not force:
        if json_output:
            _output_json({"error": "--steal requires --force"})
            raise typer.Exit(1)
        error("--steal requires --force")
        raise typer.Exit(1)

    # Collect locks to release
    locks_to_release: list[state.SampleLock] = []

    if sample_id:
        lock = state.get_lock(sample_id, include_expired=True)
        if lock:
            locks_to_release.append(lock)
        else:
            if json_output:
                _output_json({"message": f"No lock found for sample '{sample_id}'"})
                return
            info(f"No lock found for sample '{sample_id}'")
            return

    if run_id:
        run_locks = state.get_locks_for_run(run_id, include_expired=True)
        locks_to_release.extend(run_locks)

    if not locks_to_release:
        if json_output:
            _output_json({"message": "No locks to release"})
            return
        info("No locks to release")
        return

    # Deduplicate by sample_id
    seen = set()
    unique_locks = []
    for lock in locks_to_release:
        if lock.sample_id not in seen:
            seen.add(lock.sample_id)
            unique_locks.append(lock)
    locks_to_release = unique_locks

    # Show what will be unlocked
    if not json_output and not force:
        console.print()
        console.print(
            f"[bold yellow]Warning:[/bold yellow] About to release {len(locks_to_release)} lock(s):"
        )
        console.print()

        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("Sample ID")
        table.add_column("Run ID")
        table.add_column("Host")
        table.add_column("User")

        for lock in locks_to_release:
            table.add_row(
                lock.sample_id,
                lock.run_id,
                lock.hostname,
                lock.username,
            )

        console.print(table)
        console.print()

        confirm = typer.confirm("Release these locks?", default=False)
        if not confirm:
            info("Aborted")
            raise typer.Exit(0)

    # Release locks
    released_count = 0
    with connect() as conn:
        for lock in locks_to_release:
            if steal:
                # Delete regardless of run_id
                cursor = conn.execute(
                    "DELETE FROM sample_locks WHERE sample_id = ?",
                    (lock.sample_id,),
                )
            else:
                # Only delete if we "own" it (same run_id as specified, or any if only sample_id given)
                if run_id:
                    cursor = conn.execute(
                        "DELETE FROM sample_locks WHERE sample_id = ? AND run_id = ?",
                        (lock.sample_id, run_id),
                    )
                else:
                    cursor = conn.execute(
                        "DELETE FROM sample_locks WHERE sample_id = ?",
                        (lock.sample_id,),
                    )
            released_count += cursor.rowcount
        conn.commit()

    if json_output:
        _output_json(
            {
                "success": True,
                "released_count": released_count,
                "samples": [lock.sample_id for lock in locks_to_release],
            }
        )
    else:
        success(f"Released {released_count} lock(s)")


@state_app.command("cleanup")
def state_cleanup(
    stale_after: Annotated[
        int,
        typer.Option(
            "--stale-after",
            help="Hours after which locks/runs are considered stale",
        ),
    ] = state.DEFAULT_LOCK_TTL_HOURS,
    mark_failed: Annotated[
        bool,
        typer.Option(
            "--mark-failed/--no-mark-failed",
            help="Mark stale 'running' runs as 'failed'",
        ),
    ] = True,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            "-n",
            help="Show what would be cleaned without doing it",
        ),
    ] = False,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Clean up expired locks and stale runs.

    Removes expired sample locks and optionally marks runs that have been
    'running' for longer than the TTL as 'failed'.

    Examples:

        # Clean up with default TTL (72 hours)
        nvd state cleanup

        # Use custom stale threshold
        nvd state cleanup --stale-after 24

        # Dry run - show what would be cleaned
        nvd state cleanup --dry-run

        # Only clean locks, don't mark runs as failed
        nvd state cleanup --no-mark-failed
    """
    ensure_db_exists(json_output)

    results = {
        "expired_locks_removed": 0,
        "stale_runs_marked_failed": 0,
    }

    if dry_run:
        # Count what would be cleaned
        locks = state.list_locks(include_expired=True)
        expired_locks = [lock for lock in locks if lock.is_expired]
        results["expired_locks_removed"] = len(expired_locks)

        if mark_failed:
            # Count stale running runs
            from datetime import timedelta, timezone

            cutoff = datetime.now(timezone.utc) - timedelta(hours=stale_after)
            runs = state.list_runs(status="running")
            stale_runs = []
            for run in runs:
                started = datetime.fromisoformat(run.started_at.replace("Z", "+00:00"))
                if started.tzinfo is None:
                    started = started.replace(tzinfo=timezone.utc)
                if started < cutoff:
                    stale_runs.append(run)
            results["stale_runs_marked_failed"] = len(stale_runs)

        if json_output:
            _output_json({"dry_run": True, **results})
        else:
            console.print()
            info("Dry run - no changes made")
            console.print()
            console.print(
                f"  Expired locks that would be removed: {results['expired_locks_removed']}"
            )
            if mark_failed:
                console.print(
                    f"  Stale runs that would be marked failed: {results['stale_runs_marked_failed']}"
                )
            console.print()
        return

    # Actually clean up
    results["expired_locks_removed"] = state.cleanup_expired_locks()

    if mark_failed:
        results["stale_runs_marked_failed"] = state.mark_stale_runs_failed(
            ttl_hours=stale_after
        )

    if json_output:
        _output_json({"success": True, **results})
    else:
        total = sum(results.values())
        if total == 0:
            info("Nothing to clean up")
        else:
            success("Cleanup complete")
            if results["expired_locks_removed"] > 0:
                console.print(
                    f"  Expired locks removed: {results['expired_locks_removed']}"
                )
            if results["stale_runs_marked_failed"] > 0:
                console.print(
                    f"  Stale runs marked failed: {results['stale_runs_marked_failed']}"
                )
        console.print()


def _get_directory_stats(directory: Path) -> tuple[int, int, int]:
    """
    Get statistics for a directory.

    Returns:
        Tuple of (file_count, dir_count, total_size_bytes)
    """
    file_count = 0
    dir_count = 0
    total_size = 0

    for item in directory.rglob("*"):
        if item.is_file():
            file_count += 1
            total_size += item.stat().st_size
        elif item.is_dir():
            dir_count += 1

    return file_count, dir_count, total_size


def _get_file_metadata(directory: Path) -> dict[str, int]:
    """
    Get metadata (relative path -> size) for all files in a directory.

    Returns:
        Dict mapping relative paths to file sizes.
    """
    metadata = {}
    for item in directory.rglob("*"):
        if item.is_file():
            rel_path = str(item.relative_to(directory))
            metadata[rel_path] = item.stat().st_size
    return metadata


def _verify_copy(source: Path, dest: Path) -> tuple[bool, str]:
    """
    Verify that a copy operation succeeded by comparing metadata.

    Compares file count, relative paths, and sizes between source and dest.

    Returns:
        Tuple of (success, message)
    """
    source_meta = _get_file_metadata(source)
    dest_meta = _get_file_metadata(dest)

    if len(source_meta) != len(dest_meta):
        return False, (
            f"File count mismatch: source has {len(source_meta)}, "
            f"dest has {len(dest_meta)}"
        )

    missing_in_dest = set(source_meta.keys()) - set(dest_meta.keys())
    if missing_in_dest:
        return (
            False,
            f"Missing files in destination: {', '.join(sorted(missing_in_dest)[:5])}...",
        )

    extra_in_dest = set(dest_meta.keys()) - set(source_meta.keys())
    if extra_in_dest:
        return (
            False,
            f"Extra files in destination: {', '.join(sorted(extra_in_dest)[:5])}...",
        )

    size_mismatches = []
    for path, source_size in source_meta.items():
        dest_size = dest_meta[path]
        if source_size != dest_size:
            size_mismatches.append(f"{path}: {source_size} vs {dest_size}")

    if size_mismatches:
        return False, f"Size mismatches: {', '.join(size_mismatches[:5])}..."

    return True, f"Verified {len(source_meta)} files"


def _update_setup_conf(new_state_dir: Path) -> bool:
    """
    Update ~/.nvd/setup.conf to point to the new state directory.

    Returns:
        True if setup.conf was updated, False if it doesn't exist.
    """
    from py_nvd.cli.commands.setup import NVD_HOME

    setup_conf = NVD_HOME / "setup.conf"
    if not setup_conf.exists():
        return False

    lines = setup_conf.read_text().splitlines()
    new_lines = []
    updated = False

    for line in lines:
        if line.startswith("NVD_STATE_DIR="):
            new_lines.append(f"NVD_STATE_DIR={new_state_dir}")
            updated = True
        else:
            new_lines.append(line)

    if updated:
        setup_conf.write_text("\n".join(new_lines) + "\n")

    return updated


def _count_uncompacted_hits(state_dir: Path) -> int:
    """Count uncompacted hit files in the state directory."""
    hits_dir = state_dir / "hits" / "month=NULL"
    if not hits_dir.exists():
        return 0
    return len(list(hits_dir.rglob("*.parquet")))


@state_app.command("move")
def state_move(
    destination: Annotated[
        Path,
        typer.Argument(help="New location for the state directory"),
    ],
    force: bool = typer.Option(
        False,
        "--force",
        "-f",
        help="Overwrite destination if it exists",
    ),
    keep_source: bool = typer.Option(
        False,
        "--keep-source",
        help="Copy instead of move (keep original)",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        "-n",
        help="Show what would happen without doing it",
    ),
) -> None:
    """
    Move the state directory to a new location.

    Copies the entire state directory to the destination, verifies the copy
    succeeded, updates ~/.nvd/setup.conf to point to the new location, and
    (unless --keep-source) removes the original.

    Examples:

        # Move to a new location
        nvd state move /scratch/nvd-state

        # Preview what would happen
        nvd state move /scratch/nvd-state --dry-run

        # Copy instead of move (keep original)
        nvd state move /backup/nvd-state --keep-source

        # Overwrite existing destination
        nvd state move /scratch/nvd-state --force
    """
    import shutil

    source = get_state_dir()
    dest = destination.expanduser().resolve()

    # Validation
    if not source.exists():
        error(f"Source state directory does not exist: {source}")
        raise typer.Exit(1)

    if source == dest:
        error("Source and destination are the same")
        raise typer.Exit(1)

    # Check for nested paths
    try:
        dest.relative_to(source)
        error("Destination cannot be inside source directory")
        raise typer.Exit(1)
    except ValueError:
        pass  # Not nested, good

    try:
        source.relative_to(dest)
        error("Source cannot be inside destination directory")
        raise typer.Exit(1)
    except ValueError:
        pass  # Not nested, good

    if dest.exists() and not force:
        error(f"Destination already exists: {dest}")
        info("Use --force to overwrite")
        raise typer.Exit(1)

    # Check destination parent exists and is writable
    dest_parent = dest.parent
    if not dest_parent.exists():
        error(f"Destination parent directory does not exist: {dest_parent}")
        raise typer.Exit(1)

    if not os.access(dest_parent, os.W_OK):
        error(f"Destination parent directory is not writable: {dest_parent}")
        raise typer.Exit(1)

    # Gather stats
    file_count, dir_count, total_size = _get_directory_stats(source)
    uncompacted_hits = _count_uncompacted_hits(source)

    # Check for setup.conf
    from py_nvd.cli.commands.setup import NVD_HOME

    setup_conf_exists = (NVD_HOME / "setup.conf").exists()

    if dry_run:
        console.print()
        info("Dry run - no changes will be made")
        console.print()
        console.print(f"  Source: {source}")
        console.print(f"  Destination: {dest}")
        console.print(f"  Operation: {'copy' if keep_source else 'move'}")
        console.print()
        console.print(f"  Files: {file_count:,}")
        console.print(f"  Directories: {dir_count:,}")
        console.print(f"  Total size: {_format_size(total_size)}")
        console.print()

        if uncompacted_hits > 100:
            warning(
                f"Found {uncompacted_hits:,} uncompacted hit files. "
                "Consider running 'nvd hits compact' before moving to reduce transfer time."
            )
            console.print()

        if setup_conf_exists:
            console.print("  [green]✓[/green] ~/.nvd/setup.conf will be updated")
        else:
            console.print(
                "  [yellow]![/yellow] No ~/.nvd/setup.conf found (config won't be updated)"
            )

        if dest.exists():
            console.print(
                "  [yellow]![/yellow] Destination exists and will be overwritten (--force)"
            )

        console.print()
        return

    # Actually perform the move
    console.print()
    operation = "Copying" if keep_source else "Moving"
    console.print(f"{operation} state directory...")
    console.print(f"  From: {source}")
    console.print(f"  To: {dest}")
    console.print(f"  Size: {_format_size(total_size)} ({file_count:,} files)")
    console.print()

    if uncompacted_hits > 100:
        warning(
            f"Found {uncompacted_hits:,} uncompacted hit files. This may take a while."
        )
        console.print()

    # Remove existing destination if force
    if dest.exists() and force:
        info(f"Removing existing destination: {dest}")
        shutil.rmtree(dest)

    # Copy the directory
    try:
        shutil.copytree(source, dest)
    except Exception as e:
        error(f"Copy failed: {e}")
        raise typer.Exit(1)

    # Verify the copy
    info("Verifying copy...")

    verified, verify_msg = _verify_copy(source, dest)
    if not verified:
        error(f"Copy verification failed: {verify_msg}")
        warning("Destination may be incomplete. Source was not modified.")
        raise typer.Exit(1)

    success(verify_msg)

    # Update setup.conf
    config_updated = _update_setup_conf(dest)
    if config_updated:
        success("Updated ~/.nvd/setup.conf")
    else:
        warning("No ~/.nvd/setup.conf found - config not updated")
        console.print(f"  Set NVD_STATE_DIR={dest} in your environment")

    # Remove source if not keeping
    if not keep_source:
        info("Removing source directory...")
        try:
            shutil.rmtree(source)
            success(f"Removed: {source}")
        except Exception as e:
            warning(f"Could not remove source: {e}")
            console.print("  Move completed but source remains")

    console.print()
    success("State directory moved successfully!")
    console.print()
    console.print("[bold]Next steps:[/bold]")
    console.print("  1. Restart your shell or run: source ~/.bashrc (or ~/.zshrc)")
    console.print("  2. Verify with: nvd state path")
    console.print()
