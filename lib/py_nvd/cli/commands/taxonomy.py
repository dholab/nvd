"""Taxonomy cache commands for the NVD CLI."""

from __future__ import annotations

import json
from pathlib import Path  # noqa: TC003  # Typer needs Path available for CLI parsing.
from typing import Annotated

import typer

from py_nvd import taxonomy as taxonomy_runtime
from py_nvd.cli.utils import console, info, success, warning
from py_nvd.paths import get_taxdump_dir

_SIZE_UNIT_THRESHOLD = 1024

taxonomy_app = typer.Typer(
    name="taxonomy",
    help="Inspect and prepare the NCBI taxonomy cache",
    no_args_is_help=True,
)


def _format_size(size_bytes: int) -> str:
    """Format bytes as a human-readable size."""
    size = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if size < _SIZE_UNIT_THRESHOLD:
            return f"{size:.1f} {unit}"
        size /= _SIZE_UNIT_THRESHOLD
    return f"{size:.1f} TB"


def _taxonomy_status(
    *,
    taxonomy_dir: Path | None,
) -> dict[str, object]:
    """Return taxonomy cache status without creating or downloading files."""
    taxdump_dir = get_taxdump_dir(taxonomy_dir=taxonomy_dir)
    taxonomy_db = taxdump_dir / "taxonomy.sqlite"

    dmp_files: dict[str, int] = {}
    if taxdump_dir.exists():
        for dmp_file in taxdump_dir.glob("*.dmp"):
            dmp_files[dmp_file.name] = dmp_file.stat().st_size

    taxonomy_db_size = taxonomy_db.stat().st_size if taxonomy_db.exists() else 0
    return {
        "taxdump_dir": str(taxdump_dir),
        "taxdump_exists": taxdump_dir.exists(),
        "taxonomy_db": {
            "path": str(taxonomy_db),
            "exists": taxonomy_db.exists(),
            "size_bytes": taxonomy_db_size,
        },
        "dmp_files": dmp_files,
    }


def _print_taxonomy_status(status: dict[str, object]) -> None:
    """Print taxonomy cache status for humans."""
    taxdump_dir = status["taxdump_dir"]
    taxonomy_db = status["taxonomy_db"]
    assert isinstance(taxonomy_db, dict)
    dmp_files = status["dmp_files"]
    assert isinstance(dmp_files, dict)

    console.print()
    console.print(f"[bold]Taxonomy Cache:[/bold] {taxdump_dir}")
    console.print()

    if not status["taxdump_exists"]:
        warning("Taxonomy cache directory does not exist")
        info("Run 'nvd taxonomy ensure' to create it")
        console.print()
        return

    if taxonomy_db["exists"]:
        size_bytes = taxonomy_db["size_bytes"]
        assert isinstance(size_bytes, int)
        success(f"taxonomy.sqlite: {_format_size(size_bytes)}")
    else:
        info("taxonomy.sqlite: not found")

    if dmp_files:
        console.print()
        console.print("[bold]Taxonomy dump files:[/bold]")
        for name, size in sorted(dmp_files.items()):
            assert isinstance(name, str)
            assert isinstance(size, int)
            console.print(f"  {name}: {_format_size(size)}")
    else:
        console.print()
        info("No .dmp files found")

    console.print()


@taxonomy_app.command("status")
def taxonomy_status(
    taxonomy_dir: Annotated[
        Path | None,
        typer.Option(
            "--taxonomy-dir",
            help="Explicit taxonomy directory; overrides config-directory fallback",
        ),
    ] = None,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Output as JSON"),
    ] = False,
) -> None:
    """Show taxonomy cache status."""
    status = _taxonomy_status(taxonomy_dir=taxonomy_dir)
    if json_output:
        console.print_json(json.dumps(status))
        return
    _print_taxonomy_status(status)


@taxonomy_app.command("ensure")
def taxonomy_ensure(
    taxonomy_dir: Annotated[
        Path | None,
        typer.Option(
            "--taxonomy-dir",
            help="Explicit taxonomy directory; overrides config-directory fallback",
        ),
    ] = None,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Output as JSON"),
    ] = False,
    taxonomy_refresh: Annotated[
        str,
        typer.Option(
            "--taxonomy-refresh",
            help="Refresh policy: missing, stale, or force",
        ),
    ] = taxonomy_runtime.TaxonomyRefresh.MISSING.value,
    taxonomy_max_age_days: Annotated[
        int,
        typer.Option(
            "--taxonomy-max-age-days",
            help="Freshness threshold for stale refreshes and warnings",
        ),
    ] = taxonomy_runtime.DEFAULT_TAXONOMY_MAX_AGE_DAYS,
) -> None:
    """Download or rebuild the taxonomy cache if needed."""
    taxdump_dir = taxonomy_runtime.ensure_taxonomy_available(
        taxonomy_dir=taxonomy_dir,
        refresh=taxonomy_refresh,
        max_age_days=taxonomy_max_age_days,
    )
    status = _taxonomy_status(taxonomy_dir=taxdump_dir)
    if json_output:
        console.print_json(json.dumps(status))
        return
    success(f"Taxonomy cache ready: {taxdump_dir}")
