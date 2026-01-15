"""
CLI commands for managing and querying metagenomic hits.

Commands:
    nvd hits export    - Export hits to various formats
    nvd hits stats     - Show summary statistics
    nvd hits lookup    - Look up a hit by key or sequence
    nvd hits trace     - Show full details for a hit including sequence
    nvd hits recur     - Find recurring hits across samples
    nvd hits timeline  - Show discovery timeline as histogram
    nvd hits compact   - Compact uncompacted hits into monthly partitions
"""

from __future__ import annotations

import sys
from dataclasses import asdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Annotated

import polars as pl
import typer
from rich.table import Table

from py_nvd.cli.utils import console, error, info
from py_nvd.db import get_state_db_path
from py_nvd.hits import (
    compact_hits,
    get_discovery_timeline,
    get_hit,
    get_hit_sequence,
    get_hit_stats,
    get_recurring_hits,
    is_valid_hit_key,
    list_hit_observations,
    list_hits,
    list_hits_with_observations,
    lookup_hit,
)
from py_nvd.models import TimelineBucket
from py_nvd.state import get_run_ids_for_sample_sets

hits_app = typer.Typer(
    name="hits",
    help="Export and query metagenomic hits",
    no_args_is_help=True,
)


@hits_app.command("stats")
def hits_stats(
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
    top: Annotated[
        int,
        typer.Option(
            "--top",
            "-n",
            help="Number of top recurring hits to show",
        ),
    ] = 5,
) -> None:
    """
    Show summary statistics for hits.

    Displays counts, length/GC distributions, date range, and top recurring
    hits (sequences appearing in multiple samples).

    Examples:

        # Show stats with default formatting
        nvd hits stats

        # Output as JSON (for scripting)
        nvd hits stats --json

        # Show top 10 recurring hits
        nvd hits stats --top 10
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    # Get data from library
    stats = get_hit_stats()
    top_recurring = get_recurring_hits(min_samples=2, limit=top) if top > 0 else []

    # Handle empty database
    if stats.total_hits == 0:
        if json_output:
            console.print_json(data=_stats_to_dict(stats, []))
        else:
            info("No hits in database")
        return

    # Output
    if json_output:
        console.print_json(data=_stats_to_dict(stats, top_recurring))
    else:
        _print_stats_rich(stats, top_recurring)


def _stats_to_dict(stats, top_recurring: list) -> dict:
    """Convert stats and recurring hits to JSON-serializable dict."""
    result = asdict(stats)
    result["top_recurring"] = [asdict(hit) for hit in top_recurring]
    return result


def _print_stats_rich(stats, top_recurring: list) -> None:
    """Print stats with Rich formatting."""
    console.print("\n[bold]Hit Statistics[/bold]")
    console.print("─" * 40)

    # Basic counts
    console.print(f"  Total unique hits:     {stats.total_hits:,}")
    console.print(f"  Total observations:    {stats.total_observations:,}")
    console.print(f"  Samples represented:   {stats.unique_samples:,}")
    console.print(f"  Runs represented:      {stats.unique_runs:,}")

    # Date range
    if stats.date_first and stats.date_last:
        console.print(
            f"  Date range:            {stats.date_first} → {stats.date_last}"
        )

    # Length distribution
    if stats.length_min is not None:
        console.print("\n[bold]Length distribution:[/bold]")
        console.print(
            f"  Min: {stats.length_min:,} bp  │  "
            f"Median: {stats.length_median:,.0f} bp  │  "
            f"Max: {stats.length_max:,} bp"
        )

    # GC distribution
    if stats.gc_min is not None:
        console.print("\n[bold]GC content:[/bold]")
        console.print(
            f"  Min: {stats.gc_min * 100:.1f}%  │  "
            f"Median: {stats.gc_median * 100:.1f}%  │  "
            f"Max: {stats.gc_max * 100:.1f}%"
        )

    # Top recurring hits
    console.print(f"\n[bold]Top {len(top_recurring)} recurring hits:[/bold]")
    if not top_recurring:
        console.print(
            "  [dim]No recurring hits found (all hits appear in only 1 sample)[/dim]"
        )
    else:
        table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
        table.add_column("Hit Key", style="cyan")
        table.add_column("Length", justify="right")
        table.add_column("Samples", justify="right")
        table.add_column("Observations", justify="right")

        for hit in top_recurring:
            table.add_row(
                hit.hit_key,
                f"{hit.sequence_length:,} bp",
                str(hit.sample_count),
                str(hit.observation_count),
            )

        console.print(table)

    console.print()


@hits_app.command("lookup")
def hits_lookup(
    query: Annotated[
        str,
        typer.Argument(
            help="Hit key (32 hex chars) or DNA sequence. Use '-' for stdin.",
        ),
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
    Look up a hit by key or sequence.

    Searches for a hit matching the given key or sequence. If found, shows
    the hit metadata and all observations (samples/runs where it appeared).

    For sequences, computes the canonical hash (strand-agnostic) and looks
    up by that key. This means you can paste a sequence in either orientation.

    Examples:

        # Look up by hit key
        nvd hits lookup a1b2c3d4e5f67890...

        # Look up by sequence
        nvd hits lookup ACGTACGTACGT...

        # Look up from stdin (e.g., from a file)
        cat sequence.txt | nvd hits lookup -

        # Output as JSON
        nvd hits lookup ACGTACGT... --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    # Handle stdin
    if query == "-":
        query = sys.stdin.read().strip()
        if not query:
            error("No input provided on stdin")
            raise typer.Exit(1)

    # Look up the hit
    hit = lookup_hit(query)

    if hit is None:
        if json_output:
            console.print_json(data={"found": False, "query": query[:50]})
        else:
            info(
                f"Hit not found for query: {query[:50]}{'...' if len(query) > 50 else ''}"
            )
        return

    # Get observations
    observations = list_hit_observations(hit_key=hit.hit_key)

    # Output
    if json_output:
        console.print_json(data=_lookup_to_dict(hit, observations))
    else:
        _print_lookup_rich(hit, observations)


def _lookup_to_dict(hit, observations: list) -> dict:
    """Convert lookup result to JSON-serializable dict."""
    result = {
        "found": True,
        "hit": asdict(hit),
        "observations": [asdict(obs) for obs in observations],
    }
    # Remove compressed sequence from JSON (not useful, binary data)
    del result["hit"]["sequence_compressed"]
    return result


def _print_lookup_rich(hit, observations: list) -> None:
    """Print lookup result with Rich formatting."""
    console.print(f"\n[bold green]Hit found:[/bold green] [cyan]{hit.hit_key}[/cyan]")
    console.print()
    console.print(
        f"  Length: {hit.sequence_length:,} bp  │  "
        f"GC: {hit.gc_content * 100:.1f}%  │  "
        f"First seen: {hit.first_seen_date}"
    )

    # Observations
    console.print(f"\n[bold]Observations ({len(observations)} total):[/bold]")
    if observations:
        # Look up run_ids for all sample_set_ids
        sample_set_ids = list({obs.sample_set_id for obs in observations})
        run_id_map = get_run_ids_for_sample_sets(sample_set_ids)

        table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
        table.add_column("Date", style="dim")
        table.add_column("Sample")
        table.add_column("Sample Set ID")
        table.add_column("Run")
        table.add_column("Contig", style="dim")

        for obs in observations:
            run_id = run_id_map.get(obs.sample_set_id, "-")
            table.add_row(
                obs.run_date,
                obs.sample_id,
                obs.sample_set_id,
                run_id,
                obs.contig_id or "-",
            )

        console.print(table)
    else:
        console.print("  [dim]No observations recorded[/dim]")

    console.print()


@hits_app.command("trace")
def hits_trace(
    hit_key: Annotated[
        str,
        typer.Argument(
            help="Hit key (32 hex chars)",
        ),
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
    Show full details for a hit including sequence.

    Deep-dive on a known hit key. Shows metadata, all observations,
    and the full sequence (wrapped at 80 characters).

    Examples:

        # Trace a hit
        nvd hits trace a1b2c3d4e5f67890...

        # Output as JSON (includes full sequence)
        nvd hits trace a1b2c3d4e5f67890... --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    # Validate key format (32 hex chars)
    hit_key = hit_key.strip().lower()
    if not is_valid_hit_key(hit_key):
        error(f"Invalid hit key: expected 32 hex characters, got {len(hit_key)}")
        raise typer.Exit(1)

    # Look up the hit
    hit = get_hit(hit_key)

    if hit is None:
        if json_output:
            console.print_json(data={"found": False, "hit_key": hit_key})
        else:
            info(f"Hit not found: {hit_key}")
        return

    # Get observations and sequence
    observations = list_hit_observations(hit_key=hit.hit_key)
    sequence = get_hit_sequence(hit)

    # Output
    if json_output:
        console.print_json(data=_trace_to_dict(hit, observations, sequence))
    else:
        _print_trace_rich(hit, observations, sequence)


def _trace_to_dict(hit, observations: list, sequence: str) -> dict:
    """Convert trace result to JSON-serializable dict."""
    result = {
        "found": True,
        "hit": asdict(hit),
        "observations": [asdict(obs) for obs in observations],
        "sequence": sequence,
    }
    del result["hit"]["sequence_compressed"]
    return result


def _print_trace_rich(hit, observations: list, sequence: str) -> None:
    """Print trace result with Rich formatting."""
    # Reuse lookup display for metadata + observations
    _print_lookup_rich(hit, observations)

    # Add sequence section
    console.print("[bold]Sequence:[/bold]")
    for i in range(0, len(sequence), 80):
        console.print(f"  {sequence[i : i + 80]}")
    console.print()


@hits_app.command("timeline")
def hits_timeline(
    granularity: Annotated[
        str,
        typer.Option(
            "--granularity",
            "-g",
            help="Time period grouping: day, week, month, year",
        ),
    ] = "month",
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Show discovery timeline as ASCII histogram.

    Displays when new hits were first discovered, grouped by time period.
    Useful for understanding discovery rate trends over time.

    Examples:

        # Show monthly timeline (default)
        nvd hits timeline

        # Show weekly timeline
        nvd hits timeline --granularity week

        # Output as JSON
        nvd hits timeline --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    valid_granularities = ("day", "week", "month", "year")
    if granularity not in valid_granularities:
        error(
            f"Invalid granularity: {granularity!r}. "
            f"Must be one of: {', '.join(valid_granularities)}"
        )
        raise typer.Exit(1)

    buckets = get_discovery_timeline(granularity=granularity)

    if json_output:
        console.print_json(data=[asdict(b) for b in buckets])
    else:
        filled = _fill_timeline_gaps(buckets, granularity)
        _print_timeline_rich(filled, granularity)


def _fill_timeline_gaps(buckets: list, granularity: str) -> list[TimelineBucket]:
    """Fill gaps in timeline with zero-count buckets."""
    if not buckets:
        return []

    existing = {b.period: b.new_hits for b in buckets}
    all_periods = _generate_period_range(
        buckets[0].period, buckets[-1].period, granularity
    )

    return [TimelineBucket(period=p, new_hits=existing.get(p, 0)) for p in all_periods]


def _generate_period_range(start: str, end: str, granularity: str) -> list[str]:
    """Generate all periods between start and end (inclusive)."""
    if granularity == "year":
        start_year, end_year = int(start), int(end)
        return [str(y) for y in range(start_year, end_year + 1)]

    elif granularity == "month":
        # Parse "YYYY-MM"
        start_dt = datetime.strptime(start, "%Y-%m")
        end_dt = datetime.strptime(end, "%Y-%m")
        periods = []
        current = start_dt
        while current <= end_dt:
            periods.append(current.strftime("%Y-%m"))
            # Advance to next month
            if current.month == 12:
                current = current.replace(year=current.year + 1, month=1)
            else:
                current = current.replace(month=current.month + 1)
        return periods

    elif granularity == "week":
        # Parse "YYYY-Www" (ISO week)
        start_dt = datetime.strptime(start + "-1", "%Y-W%W-%w")
        end_dt = datetime.strptime(end + "-1", "%Y-W%W-%w")
        periods = []
        current = start_dt
        while current <= end_dt:
            periods.append(current.strftime("%Y-W%W"))
            current += timedelta(weeks=1)
        return periods

    elif granularity == "day":
        # Parse "YYYY-MM-DD"
        start_dt = datetime.strptime(start, "%Y-%m-%d")
        end_dt = datetime.strptime(end, "%Y-%m-%d")
        periods = []
        current = start_dt
        while current <= end_dt:
            periods.append(current.strftime("%Y-%m-%d"))
            current += timedelta(days=1)
        return periods

    return []


def _print_timeline_rich(buckets: list, granularity: str) -> None:
    """Print timeline as ASCII histogram."""
    if not buckets:
        info("No hits in database")
        return

    max_count = max(b.new_hits for b in buckets)
    max_bar_width = 30

    # Determine period width for alignment
    period_width = max(len(b.period) for b in buckets)

    console.print(f"\n[bold]Discovery Timeline[/bold] (by {granularity})")
    console.print("─" * 50)

    for bucket in buckets:
        bar_len = (
            int((bucket.new_hits / max_count) * max_bar_width) if max_count > 0 else 0
        )
        bar = "█" * bar_len
        padded_period = bucket.period.ljust(period_width)
        console.print(f"  {padded_period}  [cyan]{bar}[/cyan]  {bucket.new_hits}")

    total = sum(b.new_hits for b in buckets)
    console.print(f"\n[dim]Total: {total:,} unique hits discovered[/dim]\n")


@hits_app.command("recur")
def hits_recur(
    min_samples: Annotated[
        int,
        typer.Option(
            "--min-samples",
            "-s",
            help="Minimum number of distinct samples",
        ),
    ] = 2,
    min_runs: Annotated[
        int | None,
        typer.Option(
            "--min-runs",
            "-r",
            help="Minimum number of distinct runs",
        ),
    ] = None,
    limit: Annotated[
        int | None,
        typer.Option(
            "--limit",
            "-n",
            help="Maximum number of results",
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
    Find recurring hits (sequences appearing in multiple samples).

    Recurring hits may indicate persistent infections, contamination,
    or lab reagent artifacts. Use this to identify sequences that
    warrant further investigation.

    Examples:

        # Show all hits appearing in 2+ samples
        nvd hits recur

        # Show hits appearing in 5+ samples
        nvd hits recur --min-samples 5

        # Show hits appearing across multiple runs
        nvd hits recur --min-runs 2

        # Limit results
        nvd hits recur --limit 20

        # Output as JSON
        nvd hits recur --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    # Get recurring hits
    hits = get_recurring_hits(
        min_samples=min_samples,
        min_runs=min_runs,
        limit=limit,
    )

    # Output
    if json_output:
        console.print_json(data=[asdict(h) for h in hits])
    else:
        _print_recur_rich(hits, min_samples, min_runs)


def _print_recur_rich(hits: list, min_samples: int, min_runs: int | None) -> None:
    """Print recurring hits with Rich formatting."""
    # Build filter description
    filters = [f"{min_samples}+ samples"]
    if min_runs:
        filters.append(f"{min_runs}+ runs")
    filter_desc = ", ".join(filters)

    console.print(f"\n[bold]Recurring Hits[/bold] ({filter_desc})")
    console.print("─" * 50)

    if not hits:
        console.print("[dim]No recurring hits found matching criteria[/dim]\n")
        return

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Hit Key", style="cyan")
    table.add_column("Length", justify="right")
    table.add_column("Samples", justify="right")
    table.add_column("Runs", justify="right")
    table.add_column("Observations", justify="right")
    table.add_column("First Seen", style="dim")
    table.add_column("Last Seen", style="dim")

    for hit in hits:
        table.add_row(
            hit.hit_key,
            f"{hit.sequence_length:,} bp",
            str(hit.sample_count),
            str(hit.run_count),
            str(hit.observation_count),
            hit.first_seen_date,
            hit.last_seen_date,
        )

    console.print(table)
    console.print(f"\n[dim]{len(hits)} recurring hit(s) found[/dim]\n")


@hits_app.command("export")
def hits_export(
    format: Annotated[
        str,
        typer.Option(
            "--format",
            "-f",
            help="Output format: tsv, csv, parquet, fasta",
        ),
    ] = "tsv",
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output file path (default: stdout for tsv/csv/fasta)",
        ),
    ] = None,
    include_sequence: Annotated[
        bool,
        typer.Option(
            "--include-sequence",
            help="Include full sequence in tabular output (large!)",
        ),
    ] = False,
) -> None:
    """
    Export hits from the state database.

    Exports all registered hits with their observation metadata.
    For tabular formats (tsv, csv, parquet), joins hits with observations
    so each row represents one observation of a hit.

    Examples:

        # Export to TSV (stdout)
        nvd hits export

        # Export to CSV file
        nvd hits export --format csv --output hits.csv

        # Export to Parquet (requires --output)
        nvd hits export --format parquet --output hits.parquet

        # Export sequences as FASTA
        nvd hits export --format fasta --output hits.fasta

        # Include full sequences in tabular output
        nvd hits export --include-sequence --output hits_with_seqs.tsv
    """
    valid_formats = ["tsv", "csv", "parquet", "fasta"]
    if format not in valid_formats:
        error(f"Invalid format: {format!r}. Must be one of: {', '.join(valid_formats)}")
        raise typer.Exit(1)

    if format == "parquet" and output is None:
        error("Parquet format requires --output file path")
        raise typer.Exit(1)

    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    if format == "fasta":
        _export_fasta(output)
    else:
        _export_tabular(format, output, include_sequence)


def _export_tabular(
    format: str,
    output: Path | None,
    include_sequence: bool,
) -> None:
    """Export hits as tabular data using Polars."""
    rows = list_hits_with_observations()

    if not rows:
        info("No hits found in database")
        return

    # Look up run_ids for all sample_set_ids
    sample_set_ids = list({obs.sample_set_id for _, obs in rows})
    run_id_map = get_run_ids_for_sample_sets(sample_set_ids)

    # Build data for DataFrame
    data: dict[str, list] = {
        "hit_key": [],
        "sequence_length": [],
        "gc_content": [],
        "first_seen_date": [],
        "sample_set_id": [],
        "run_id": [],
        "sample_id": [],
        "contig_id": [],
        "run_date": [],
    }

    if include_sequence:
        data["sequence"] = []

    for hit, obs in rows:
        data["hit_key"].append(hit.hit_key)
        data["sequence_length"].append(hit.sequence_length)
        data["gc_content"].append(hit.gc_content)
        data["first_seen_date"].append(hit.first_seen_date)
        data["sample_set_id"].append(obs.sample_set_id)
        data["run_id"].append(run_id_map.get(obs.sample_set_id))
        data["sample_id"].append(obs.sample_id)
        data["contig_id"].append(obs.contig_id)
        data["run_date"].append(obs.run_date)

        if include_sequence:
            data["sequence"].append(get_hit_sequence(hit))

    df = pl.DataFrame(data)

    # Output
    if format == "tsv":
        content = df.write_csv(separator="\t")
        if output:
            output.write_text(content)
            info(f"Wrote {len(df)} rows to {output}")
        else:
            console.print(content, end="")
    elif format == "csv":
        content = df.write_csv()
        if output:
            output.write_text(content)
            info(f"Wrote {len(df)} rows to {output}")
        else:
            console.print(content, end="")
    elif format == "parquet":
        assert output is not None  # Validated in hits_export
        df.write_parquet(output)
        info(f"Wrote {len(df)} rows to {output}")


def _export_fasta(output: Path | None) -> None:
    """Export hits as FASTA with minimal headers."""
    hits = list_hits()

    if not hits:
        info("No hits found in database")
        return

    # Build FASTA content
    lines = []
    for hit in hits:
        sequence = get_hit_sequence(hit)
        lines.append(f">{hit.hit_key}")
        # Wrap sequence at 80 characters
        for i in range(0, len(sequence), 80):
            lines.append(sequence[i : i + 80])

    content = "\n".join(lines) + "\n"

    if output:
        output.write_text(content)
        info(f"Wrote {len(hits)} sequences to {output}")
    else:
        console.print(content, end="")


@hits_app.command("compact")
def hits_compact(
    month: Annotated[
        str | None,
        typer.Option(
            "--month",
            "-m",
            help="Only compact data for this month (YYYY-MM format)",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            "-n",
            help="Show what would be compacted without doing it",
        ),
    ] = False,
    keep_source: Annotated[
        bool,
        typer.Option(
            "--keep-source",
            help="Keep uncompacted files after compaction",
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
    Compact uncompacted hits into monthly partitions.

    Consolidates small per-sample parquet files from hits/month=NULL/ into
    larger sorted files in hits/month=YYYY-MM/ partitions. This improves
    query performance through:

    - Fewer files to open (better for DuckDB)
    - Sorted data (better compression, faster aggregations)
    - Time-based partitions (enables partition pruning)

    Examples:

        # Show what would be compacted (dry run)
        nvd hits compact --dry-run

        # Compact all uncompacted data
        nvd hits compact

        # Compact only January 2024 data
        nvd hits compact --month 2024-01

        # Compact but keep source files (for testing)
        nvd hits compact --keep-source
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    # Validate month format
    if month is not None:
        if len(month) != 7 or month[4] != "-":
            error(f"Invalid month format: {month!r}. Expected YYYY-MM (e.g., 2024-01)")
            raise typer.Exit(1)
        try:
            year, mon = int(month[:4]), int(month[5:])
            if not (1 <= mon <= 12):
                raise ValueError("Month must be 1-12")
        except ValueError as e:
            error(f"Invalid month format: {month!r}. {e}")
            raise typer.Exit(1)

    # Run compaction
    result = compact_hits(month=month, dry_run=dry_run, keep_source=keep_source)

    # Output
    if json_output:
        console.print_json(data=_compact_to_dict(result))
    else:
        _print_compact_rich(result, dry_run, keep_source)


def _compact_to_dict(result) -> dict:
    """Convert compaction result to JSON-serializable dict."""
    return {
        "dry_run": result.dry_run,
        "months": [
            {
                "month": m.month,
                "observation_count": m.observation_count,
                "unique_hits": m.unique_hits,
                "sample_sets": m.sample_sets,
                "source_files": m.source_files,
            }
            for m in result.months
        ],
        "total_observations": result.total_observations,
        "total_source_files": result.total_source_files,
        "errors": result.errors,
    }


def _print_compact_rich(result, dry_run: bool, keep_source: bool) -> None:
    """Print compaction result with Rich formatting."""
    if not result.months:
        info("No uncompacted data to compact")
        return

    action = "Would compact" if dry_run else "Compacted"
    console.print(f"\n[bold]{action} {len(result.months)} month(s)[/bold]")
    console.print("─" * 50)

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Month", style="cyan")
    table.add_column("Observations", justify="right")
    table.add_column("Unique Hits", justify="right")
    table.add_column("Sample Sets", justify="right")
    table.add_column("Source Files", justify="right")

    for m in sorted(result.months, key=lambda x: x.month, reverse=True):
        table.add_row(
            m.month,
            f"{m.observation_count:,}",
            f"{m.unique_hits:,}",
            str(m.sample_sets),
            str(m.source_files),
        )

    console.print(table)

    console.print(f"\n  Total observations: {result.total_observations:,}")
    console.print(f"  Total source files: {result.total_source_files:,}")

    if dry_run:
        console.print("\n[yellow]Dry run - no changes made[/yellow]")
    elif keep_source:
        console.print("\n[dim]Source files preserved (--keep-source)[/dim]")
    else:
        console.print("\n[green]Source files deleted[/green]")

    if result.errors:
        console.print("\n[red bold]Errors:[/red bold]")
        for err in result.errors:
            console.print(f"  [red]• {err}[/red]")

    console.print()
