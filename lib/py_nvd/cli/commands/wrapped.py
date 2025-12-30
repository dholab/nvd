# ruff: noqa: B008
"""
NVD Wrapped - Your year in metagenomics!

A fun easter egg command that summarizes your pipeline activity
for the year, Spotify Wrapped style.
"""

from __future__ import annotations

from datetime import datetime
from typing import Annotated

import typer
from rich.console import Group
from rich.panel import Panel
from rich.text import Text

from py_nvd.cli.utils import console, info
from py_nvd.db import get_state_db_path


def _format_duration(seconds: float) -> str:
    """Format seconds as human-readable duration."""
    if seconds < 60:
        return f"{int(seconds)}s"
    elif seconds < 3600:
        minutes = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{minutes}m {secs}s"
    else:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        return f"{hours}h {minutes}m"


def _make_bar(value: int, max_value: int, width: int = 20) -> str:
    """Create a simple text progress bar."""
    if max_value == 0:
        return " " * width
    filled = int((value / max_value) * width)
    return "\u2588" * filled + "\u2591" * (width - filled)


def wrapped(
    year: Annotated[
        int | None,
        typer.Option(
            "--year",
            "-y",
            help="Year to summarize (defaults to current year)",
        ),
    ] = None,
) -> None:
    """
    Your year in metagenomics - NVD Wrapped!

    See your pipeline statistics for the year, Spotify Wrapped style.

    Examples:

        # Show this year's wrapped
        nvd wrapped

        # Show 2024's wrapped
        nvd wrapped --year 2024
    """
    from py_nvd.db import connect

    db_path = get_state_db_path()
    if not db_path.exists():
        info("No pipeline data yet - run some samples and come back!")
        console.print("[dim]Your wrapped awaits...[/dim]")
        return

    target_year = year or datetime.now().year
    year_start = f"{target_year}-01-01"
    year_end = f"{target_year}-12-31 23:59:59"

    with connect() as conn:
        # Get runs for the year
        runs = conn.execute(
            """SELECT * FROM runs 
               WHERE started_at >= ? AND started_at <= ?
               ORDER BY started_at""",
            (year_start, year_end),
        ).fetchall()

        if not runs:
            _show_empty_wrapped(target_year)
            return

        # Calculate run stats
        total_runs = len(runs)
        completed_runs = sum(1 for r in runs if r["status"] == "completed")
        failed_runs = sum(1 for r in runs if r["status"] == "failed")
        success_rate = (completed_runs / total_runs * 100) if total_runs > 0 else 0

        # Find longest run
        longest_run = None
        longest_duration = 0
        for run in runs:
            if run["completed_at"] and run["started_at"]:
                try:
                    start = datetime.fromisoformat(run["started_at"])
                    end = datetime.fromisoformat(run["completed_at"])
                    duration = (end - start).total_seconds()
                    if duration > longest_duration:
                        longest_duration = duration
                        longest_run = run
                except ValueError:
                    pass

        # Find busiest month
        month_counts: dict[int, int] = {}
        for run in runs:
            try:
                month = datetime.fromisoformat(run["started_at"]).month
                month_counts[month] = month_counts.get(month, 0) + 1
            except ValueError:
                pass

        busiest_month = None
        busiest_month_count = 0
        if month_counts:
            busiest_month = max(month_counts, key=lambda m: month_counts[m])
            busiest_month_count = month_counts[busiest_month]

        month_names = [
            "",
            "January",
            "February",
            "March",
            "April",
            "May",
            "June",
            "July",
            "August",
            "September",
            "October",
            "November",
            "December",
        ]

        # Get sample stats
        samples = conn.execute(
            """SELECT ps.* FROM processed_samples ps
               JOIN runs r ON ps.run_id = r.run_id
               WHERE r.started_at >= ? AND r.started_at <= ?""",
            (year_start, year_end),
        ).fetchall()

        total_samples = len(samples)
        unique_samples = len({s["sample_id"] for s in samples})

        # Find top sample
        sample_counts: dict[str, int] = {}
        for s in samples:
            sample_counts[s["sample_id"]] = sample_counts.get(s["sample_id"], 0) + 1

        top_sample = None
        top_sample_count = 0
        if sample_counts:
            top_sample = max(sample_counts, key=lambda s: sample_counts[s])
            top_sample_count = sample_counts[top_sample]

        # Get upload stats
        uploads = conn.execute(
            """SELECT u.* FROM uploads u
               JOIN runs r ON u.sample_set_id = (
                   SELECT sample_set_id FROM processed_samples 
                   WHERE sample_id = u.sample_id AND sample_set_id = u.sample_set_id
                   LIMIT 1
               )
               WHERE u.uploaded_at >= ? AND u.uploaded_at <= ?""",
            (year_start, year_end),
        ).fetchall()

        # Simpler upload query - just count by date
        uploads = conn.execute(
            """SELECT * FROM uploads
               WHERE uploaded_at >= ? AND uploaded_at <= ?""",
            (year_start, year_end),
        ).fetchall()

        total_uploads = len(uploads)
        labkey_uploads = sum(1 for u in uploads if u["upload_target"] == "labkey")

        # Get database versions used
        db_versions = conn.execute(
            """SELECT DISTINCT db_type, version FROM databases
               WHERE registered_at >= ? AND registered_at <= ?""",
            (year_start, year_end),
        ).fetchall()

    # Build the display
    _show_wrapped(
        year=target_year,
        total_runs=total_runs,
        completed_runs=completed_runs,
        failed_runs=failed_runs,
        success_rate=success_rate,
        total_samples=total_samples,
        unique_samples=unique_samples,
        top_sample=top_sample,
        top_sample_count=top_sample_count,
        busiest_month=month_names[busiest_month] if busiest_month else None,
        busiest_month_count=busiest_month_count,
        longest_run=longest_run,
        longest_duration=longest_duration,
        total_uploads=total_uploads,
        labkey_uploads=labkey_uploads,
        db_versions=len(db_versions),
    )


def _show_empty_wrapped(year: int) -> None:
    """Show wrapped for a year with no data."""
    console.print()

    header = Text()
    header.append("\n")
    header.append("      ", style="bold")
    header.append("\U0001f9ec NVD WRAPPED ", style="bold magenta")
    header.append(str(year), style="bold cyan")
    header.append(" \U0001f9ec", style="bold")
    header.append("\n")

    console.print(
        Panel(
            header,
            border_style="magenta",
            padding=(0, 2),
        )
    )

    console.print()
    console.print("  [dim]No pipeline runs found for this year.[/dim]")
    console.print()
    console.print(
        f"  [italic]{year + 1} is your year! Get sequencing! \U0001f680[/italic]"
    )
    console.print()


def _show_wrapped(
    year: int,
    total_runs: int,
    completed_runs: int,
    failed_runs: int,
    success_rate: float,
    total_samples: int,
    unique_samples: int,
    top_sample: str | None,
    top_sample_count: int,
    busiest_month: str | None,
    busiest_month_count: int,
    longest_run: dict | None,
    longest_duration: float,
    total_uploads: int,
    labkey_uploads: int,
    db_versions: int,
) -> None:
    """Display the wrapped summary."""
    console.print()

    # Header
    header = Text()
    header.append("\n")
    header.append("          ", style="bold")
    header.append("\U0001f9ec NVD WRAPPED ", style="bold magenta")
    header.append(str(year), style="bold cyan")
    header.append(" \U0001f9ec", style="bold")
    header.append("\n")

    console.print(
        Panel(
            header,
            border_style="magenta",
            padding=(0, 2),
        )
    )

    console.print()
    console.print("  [bold italic]You had a big year in metagenomics![/bold italic]")
    console.print()

    # Main stats
    max_stat = max(total_runs, total_samples, 1)

    stats_text = Text()
    stats_text.append(f"    {total_runs:,} runs", style="bold cyan")
    stats_text.append("       ")
    stats_text.append(f"{total_samples:,} samples processed", style="bold green")
    stats_text.append("\n")
    stats_text.append(f"    {_make_bar(total_runs, max_stat, 12)}", style="cyan")
    stats_text.append("       ")
    stats_text.append(f"{_make_bar(total_samples, max_stat, 20)}", style="green")

    console.print(
        Panel(
            stats_text,
            border_style="blue",
            padding=(1, 2),
        )
    )

    console.print()

    # Superlatives
    if top_sample and top_sample_count > 1:
        console.print(
            f"  \U0001f3c6 [bold]Top Sample:[/bold] {top_sample} [dim](appeared in {top_sample_count} runs)[/dim]"
        )
        console.print()

    if busiest_month:
        console.print(
            f"  \U0001f4c5 [bold]Busiest Month:[/bold] {busiest_month} [dim]({busiest_month_count} runs)[/dim]"
        )
        console.print()

    if longest_run and longest_duration > 0:
        console.print(
            f"  \U000023f1\U0000fe0f  [bold]Longest Run:[/bold] {longest_run['run_id']} [dim]({_format_duration(longest_duration)})[/dim]"
        )
        console.print()

    if total_uploads > 0:
        if labkey_uploads > 0:
            console.print(
                f"  \U0001f4e4 [bold]Uploads:[/bold] {total_uploads:,} results [dim]({labkey_uploads:,} to LabKey)[/dim]"
            )
        else:
            console.print(
                f"  \U0001f4e4 [bold]Uploads:[/bold] {total_uploads:,} results"
            )
        console.print()

    # Success rate with appropriate emoji
    if success_rate >= 95:
        rate_emoji = "\U0001f389"  # party popper
        rate_style = "bold green"
    elif success_rate >= 80:
        rate_emoji = "\U00002705"  # check mark
        rate_style = "green"
    elif success_rate >= 50:
        rate_emoji = "\U0001f914"  # thinking face
        rate_style = "yellow"
    else:
        rate_emoji = "\U0001f4aa"  # flexed bicep (keep trying!)
        rate_style = "yellow"

    console.print(
        f"  {rate_emoji} [bold]Success Rate:[/bold] [{rate_style}]{success_rate:.0f}%[/{rate_style}] of runs completed"
    )
    console.print()

    if unique_samples != total_samples:
        console.print(
            f"  [dim]({unique_samples:,} unique samples across all runs)[/dim]"
        )
        console.print()

    # Footer
    next_year = year + 1
    footer = Text()
    footer.append("\n")
    footer.append(f"        See you in {next_year}! Keep sequencing. ", style="italic")
    footer.append("\U0001f52c", style="")
    footer.append("\n")

    console.print(
        Panel(
            footer,
            border_style="magenta",
            padding=(0, 2),
        )
    )

    console.print()
