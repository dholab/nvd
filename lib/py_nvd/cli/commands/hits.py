"""
CLI commands for managing and querying metagenomic hits.

Commands:
    nvd hits report    - Show comprehensive report for a run
    nvd hits top       - Show top taxa findings for a run
    nvd hits samples   - Show per-sample breakdown for a run
    nvd hits categories - Show findings grouped by virus family
    nvd hits quality   - Show contig quality metrics for a run
    nvd hits novel     - Show novel taxa (first-time sightings) for a run
    nvd hits movers    - Show top movers (biggest prevalence changes) for a run
    nvd hits rare      - Show rare taxa in a run
    nvd hits history   - Show history of a specific taxon
    nvd hits compare   - Compare run metrics to historical norms
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
    get_contig_quality,
    get_discovery_timeline,
    get_hit,
    get_hit_sequence,
    get_hit_stats,
    get_negative_samples,
    get_novel_taxa,
    get_rare_taxa,
    get_recurring_hits,
    get_run_comparison,
    get_run_report,
    get_sample_summaries,
    get_taxa_by_category,
    get_taxon_history,
    get_top_movers,
    get_top_taxa,
    is_valid_hit_key,
    list_hit_observations,
    list_hits,
    list_hits_with_observations,
    lookup_hit,
)
from py_nvd.models import (
    CategorySummary,
    ContigQuality,
    HitStats,
    NovelTaxon,
    RareTaxon,
    RecurringHit,
    RunComparison,
    RunReport,
    SampleSummary,
    TaxonChange,
    TaxonHistory,
    TaxonSummary,
    TimelineBucket,
)
from py_nvd.state import get_run, get_run_ids_for_sample_sets


def resolve_run_identifier(identifier: str) -> str:
    """
    Resolve a run identifier to a sample_set_id.

    Accepts either:
    - A run_id (e.g., "evil_heisenberg") - looked up in the runs table
    - A sample_set_id (e.g., "4eb47dfdadfa0d33") - returned as-is

    This allows users to use the human-friendly Nextflow run names from
    `nvd state runs` when querying hits.

    Args:
        identifier: Either a run_id or sample_set_id

    Returns:
        The sample_set_id for use in hits queries
    """
    # Try to look up as run_id first
    run = get_run(identifier)
    if run is not None:
        return run.sample_set_id

    # Not found as run_id, assume it's already a sample_set_id
    return identifier


def _format_classification(
    name: str | None,
    rank: str | None,
    taxid: int | None = None,
    *,
    include_taxid: bool = False,
) -> str:
    """Format taxonomic classification for display.

    Args:
        name: Taxonomic name (e.g., "Influenza A virus").
        rank: Taxonomic rank (e.g., "species").
        taxid: Taxonomic ID (e.g., 12345).
        include_taxid: If True, append "[taxid: N]" with dim styling.

    Returns:
        Formatted string like "Name (rank)" or "Name (rank) [dim][taxid: N][/dim]",
        or "-" if name is None.
    """
    if not name:
        return "-"
    result = name
    if rank:
        result += f" ({rank})"
    if include_taxid and taxid:
        result += f" [dim][taxid: {taxid}][/dim]"
    return result


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


def _print_stats_rich(stats: HitStats, top_recurring: list[RecurringHit]) -> None:
    """Print stats with Rich formatting."""
    console.print("\n[bold]Hit Statistics[/bold]")
    console.print("─" * 40)


@hits_app.command("report")
def hits_report(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
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
    Show comprehensive report for a pipeline run.

    Displays summary counts, top findings, quality metrics, and negative
    samples for a specific run. This is the primary reporting command for
    understanding run results.

    Examples:

        # Show report for a run (using run ID from nvd state runs)
        nvd hits report --run evil_heisenberg

        # Or using sample set ID directly
        nvd hits report --run 4eb47dfdadfa0d33

        # Output as JSON (for scripting)
        nvd hits report --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    report = get_run_report(sample_set_id)

    if report is None:
        if json_output:
            console.print_json(data={"found": False, "sample_set_id": run})
        else:
            info(f"No hits found for run: {run}")
        return

    if json_output:
        console.print_json(data=_report_to_dict(report))
    else:
        _print_report_rich(report)


def _report_to_dict(report: RunReport) -> dict:
    """Convert RunReport to JSON-serializable dict."""
    return {
        "sample_set_id": report.sample_set_id,
        "run_date": report.run_date,
        "samples_analyzed": report.samples_analyzed,
        "samples_with_viral_hits": report.samples_with_viral_hits,
        "samples_negative": report.samples_negative,
        "unique_hits": report.unique_hits,
        "viral_taxa_found": report.viral_taxa_found,
        "median_contig_length": report.median_contig_length,
        "contigs_over_500bp": report.contigs_over_500bp,
        "top_findings": [
            {
                "taxid": t.taxid,
                "name": t.name,
                "rank": t.rank,
                "sample_count": t.sample_count,
                "hit_count": t.hit_count,
                "avg_contig_length": t.avg_contig_length,
            }
            for t in report.top_findings
        ],
        "sample_summaries": [
            {
                "sample_id": s.sample_id,
                "total_hits": s.total_hits,
                "viral_taxa_count": s.viral_taxa_count,
                "taxa_detected": list(s.taxa_detected),
            }
            for s in report.sample_summaries
        ],
    }


def _print_report_rich(report: RunReport) -> None:
    """Print run report with Rich formatting."""
    # Header
    console.print(f"\n[bold]Run Report:[/bold] [cyan]{report.sample_set_id}[/cyan]")
    if report.run_date:
        # Extract just the date part
        date_str = (
            report.run_date[:10] if len(report.run_date) >= 10 else report.run_date
        )
        console.print(f"[dim]Date: {date_str}[/dim]")
    console.print()

    # Summary section
    console.print("[bold]Summary[/bold]")
    console.print("─" * 40)
    console.print(f"  Samples analyzed:      {report.samples_analyzed:,}")
    if report.samples_analyzed > 0:
        pct = report.samples_with_viral_hits / report.samples_analyzed * 100
        console.print(
            f"  Samples with hits:     {report.samples_with_viral_hits:,} ({pct:.0f}%)"
        )
    else:
        console.print(f"  Samples with hits:     {report.samples_with_viral_hits:,}")
    console.print(f"  Samples negative:      {report.samples_negative:,}")
    console.print(f"  Unique contigs:        {report.unique_hits:,}")
    console.print(f"  Taxa identified:       {report.viral_taxa_found:,}")

    # Top findings
    console.print(f"\n[bold]Top Findings[/bold] (by sample prevalence)")
    console.print("─" * 40)
    if report.top_findings:
        table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
        table.add_column("Taxon")
        table.add_column("Samples", justify="right")
        table.add_column("Contigs", justify="right")
        table.add_column("Avg Length", justify="right")

        for t in report.top_findings:
            name = t.name or "-"
            if t.rank:
                name += f" [dim]({t.rank})[/dim]"
            table.add_row(
                name,
                str(t.sample_count),
                str(t.hit_count),
                f"{t.avg_contig_length:.0f} bp",
            )

        console.print(table)
    else:
        console.print("  [dim]No viral taxa identified[/dim]")

    # Quality metrics
    console.print("\n[bold]Contig Quality[/bold]")
    console.print("─" * 40)
    console.print(f"  Median length:         {report.median_contig_length:,} bp")
    if report.unique_hits > 0:
        pct_500 = report.contigs_over_500bp / report.unique_hits * 100
        console.print(
            f"  Contigs ≥500 bp:       {report.contigs_over_500bp:,} ({pct_500:.1f}%)"
        )
    else:
        console.print(f"  Contigs ≥500 bp:       {report.contigs_over_500bp:,}")

    # Negative samples
    if report.samples_negative > 0:
        negative = get_negative_samples(report.sample_set_id)
        console.print("\n[bold]Samples with No Viral Hits[/bold]")
        console.print("─" * 40)
        # Display as comma-separated list, wrapping if needed
        if len(negative) <= 10:
            console.print(f"  {', '.join(negative)}")
        else:
            # Show first 10 with count
            console.print(f"  {', '.join(negative[:10])}")
            console.print(f"  [dim]...and {len(negative) - 10} more[/dim]")

    console.print()


@hits_app.command("top")
def hits_top(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
        ),
    ],
    limit: Annotated[
        int,
        typer.Option(
            "--limit",
            "-n",
            help="Number of top taxa to show",
        ),
    ] = 10,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Show top taxa findings for a run.

    Displays taxa ordered by sample prevalence (number of samples where
    each taxon was detected). Use --limit to control how many to show.

    Examples:

        # Show top 10 taxa (default)
        nvd hits top --run evil_heisenberg

        # Show top 20 taxa
        nvd hits top --run evil_heisenberg --limit 20

        # Output as JSON
        nvd hits top --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    taxa = get_top_taxa(sample_set_id=sample_set_id, limit=limit)

    if json_output:
        console.print_json(data=[_taxon_to_dict(t) for t in taxa])
    else:
        _print_top_taxa_rich(taxa, sample_set_id, limit)


def _taxon_to_dict(t: TaxonSummary) -> dict:
    """Convert TaxonSummary to JSON-serializable dict."""
    return {
        "taxid": t.taxid,
        "name": t.name,
        "rank": t.rank,
        "sample_count": t.sample_count,
        "hit_count": t.hit_count,
        "avg_contig_length": t.avg_contig_length,
    }


def _print_top_taxa_rich(taxa: list[TaxonSummary], run: str, limit: int) -> None:
    """Print top taxa with Rich formatting."""
    console.print(f"\n[bold]Top {limit} Taxa:[/bold] [cyan]{run}[/cyan]")
    console.print("─" * 50)

    if not taxa:
        console.print("  [dim]No viral taxa identified[/dim]\n")
        return

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Taxon")
    table.add_column("Samples", justify="right")
    table.add_column("Contigs", justify="right")
    table.add_column("Avg Length", justify="right")

    for t in taxa:
        name = t.name or "-"
        if t.rank:
            name += f" [dim]({t.rank})[/dim]"
        table.add_row(
            name,
            str(t.sample_count),
            str(t.hit_count),
            f"{t.avg_contig_length:.0f} bp",
        )

    console.print(table)
    console.print()


@hits_app.command("samples")
def hits_samples(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
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
    Show per-sample breakdown for a run.

    Displays each sample with its hit count and detected taxa. Samples
    are ordered by viral taxa count (most interesting first).

    Examples:

        # Show per-sample breakdown
        nvd hits samples --run evil_heisenberg

        # Output as JSON
        nvd hits samples --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    summaries = get_sample_summaries(sample_set_id)

    if json_output:
        console.print_json(data=[_sample_summary_to_dict(s) for s in summaries])
    else:
        _print_samples_rich(summaries, sample_set_id)


def _sample_summary_to_dict(s: SampleSummary) -> dict:
    """Convert SampleSummary to JSON-serializable dict."""
    return {
        "sample_id": s.sample_id,
        "total_hits": s.total_hits,
        "viral_taxa_count": s.viral_taxa_count,
        "taxa_detected": list(s.taxa_detected),
    }


def _print_samples_rich(summaries: list[SampleSummary], run: str) -> None:
    """Print sample summaries with Rich formatting."""
    console.print(f"\n[bold]Per-Sample Breakdown:[/bold] [cyan]{run}[/cyan]")
    console.print("─" * 60)

    if not summaries:
        console.print("  [dim]No samples found[/dim]\n")
        return

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Sample")
    table.add_column("Hits", justify="right")
    table.add_column("Taxa", justify="right")
    table.add_column("Detected Taxa")

    for s in summaries:
        # Truncate taxa list if too long
        taxa_str = ", ".join(s.taxa_detected[:5])
        if len(s.taxa_detected) > 5:
            taxa_str += f" [dim](+{len(s.taxa_detected) - 5} more)[/dim]"
        elif not s.taxa_detected:
            taxa_str = "[dim]-[/dim]"

        table.add_row(
            s.sample_id,
            str(s.total_hits),
            str(s.viral_taxa_count),
            taxa_str,
        )

    console.print(table)
    console.print(f"\n[dim]{len(summaries)} sample(s)[/dim]\n")


@hits_app.command("categories")
def hits_categories(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
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
    Show findings grouped by virus family.

    Displays a high-level triage view with taxa grouped by virus family
    (e.g., Orthoherpesviridae, Picornaviridae). Useful for quickly
    understanding the types of viruses detected.

    Examples:

        # Show category breakdown
        nvd hits categories --run evil_heisenberg

        # Output as JSON
        nvd hits categories --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    categories = get_taxa_by_category(sample_set_id)

    if json_output:
        console.print_json(data=[_category_to_dict(c) for c in categories])
    else:
        _print_categories_rich(categories, sample_set_id)


def _category_to_dict(c: CategorySummary) -> dict:
    """Convert CategorySummary to JSON-serializable dict."""
    return {
        "category": c.category,
        "sample_count": c.sample_count,
        "hit_count": c.hit_count,
        "taxa": list(c.taxa),
    }


def _print_categories_rich(categories: list[CategorySummary], run: str) -> None:
    """Print category summaries with Rich formatting."""
    console.print(f"\n[bold]Findings by Virus Family:[/bold] [cyan]{run}[/cyan]")
    console.print("─" * 60)

    if not categories:
        console.print("  [dim]No viral taxa identified[/dim]\n")
        return

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Family")
    table.add_column("Samples", justify="right")
    table.add_column("Contigs", justify="right")
    table.add_column("Taxa")

    for c in categories:
        # Truncate taxa list if too long
        taxa_str = ", ".join(c.taxa[:3])
        if len(c.taxa) > 3:
            taxa_str += f" [dim](+{len(c.taxa) - 3} more)[/dim]"

        table.add_row(
            c.category,
            str(c.sample_count),
            str(c.hit_count),
            taxa_str,
        )

    console.print(table)
    console.print(f"\n[dim]{len(categories)} family/families[/dim]\n")


@hits_app.command("quality")
def hits_quality(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
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
    Show contig quality metrics for a run.

    Displays assembly quality statistics including length distribution
    and GC content. Useful for QC and troubleshooting.

    Examples:

        # Show quality metrics
        nvd hits quality --run evil_heisenberg

        # Output as JSON
        nvd hits quality --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    quality = get_contig_quality(sample_set_id)

    if quality is None:
        if json_output:
            console.print_json(data={"found": False, "sample_set_id": sample_set_id})
        else:
            info(f"No hits found for run: {run}")
        return

    if json_output:
        console.print_json(data=_quality_to_dict(quality))
    else:
        _print_quality_rich(quality, sample_set_id)


def _quality_to_dict(q: ContigQuality) -> dict:
    """Convert ContigQuality to JSON-serializable dict."""
    return {
        "total_contigs": q.total_contigs,
        "total_bases": q.total_bases,
        "mean_length": q.mean_length,
        "median_length": q.median_length,
        "max_length": q.max_length,
        "contigs_500bp_plus": q.contigs_500bp_plus,
        "contigs_1kb_plus": q.contigs_1kb_plus,
        "mean_gc": q.mean_gc,
    }


def _print_quality_rich(quality: ContigQuality, run: str) -> None:
    """Print quality metrics with Rich formatting."""
    console.print(f"\n[bold]Contig Quality:[/bold] [cyan]{run}[/cyan]")
    console.print("─" * 40)

    console.print(f"  Total contigs:         {quality.total_contigs:,}")
    console.print(f"  Total bases:           {quality.total_bases:,}")
    console.print()
    console.print("[bold]Length Distribution[/bold]")
    console.print(f"  Mean:                  {quality.mean_length:,.0f} bp")
    console.print(f"  Median:                {quality.median_length:,} bp")
    console.print(f"  Max:                   {quality.max_length:,} bp")
    console.print()

    # Percentage calculations
    if quality.total_contigs > 0:
        pct_500 = quality.contigs_500bp_plus / quality.total_contigs * 100
        pct_1k = quality.contigs_1kb_plus / quality.total_contigs * 100
        console.print(
            f"  Contigs ≥500 bp:       {quality.contigs_500bp_plus:,} ({pct_500:.1f}%)"
        )
        console.print(
            f"  Contigs ≥1 kb:         {quality.contigs_1kb_plus:,} ({pct_1k:.1f}%)"
        )
    else:
        console.print(f"  Contigs ≥500 bp:       {quality.contigs_500bp_plus:,}")
        console.print(f"  Contigs ≥1 kb:         {quality.contigs_1kb_plus:,}")

    console.print()
    console.print("[bold]GC Content[/bold]")
    console.print(f"  Mean GC:               {quality.mean_gc * 100:.1f}%")
    console.print()


@hits_app.command("novel")
def hits_novel(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
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
    Show novel taxa (first-time sightings) for a run.

    Identifies taxa in this run that have never appeared in any previous run.
    Useful for detecting emerging pathogens or new contamination sources.

    Note: "Novel" means not seen in the NVD state database, not necessarily
    novel to science. Historical data availability affects what's considered novel.

    Examples:

        # Show novel taxa for a run
        nvd hits novel --run evil_heisenberg

        # Output as JSON
        nvd hits novel --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    taxa = get_novel_taxa(sample_set_id)

    if json_output:
        console.print_json(data=[_novel_taxon_to_dict(t) for t in taxa])
    else:
        _print_novel_taxa_rich(taxa, sample_set_id)


def _novel_taxon_to_dict(t: NovelTaxon) -> dict:
    """Convert NovelTaxon to JSON-serializable dict."""
    return {
        "taxid": t.taxid,
        "name": t.name,
        "rank": t.rank,
        "sample_count": t.sample_count,
        "hit_count": t.hit_count,
    }


def _print_novel_taxa_rich(taxa: list[NovelTaxon], run: str) -> None:
    """Print novel taxa with Rich formatting."""
    console.print(
        f"\n[bold]Novel Taxa (First-Time Sightings):[/bold] [cyan]{run}[/cyan]"
    )
    console.print("─" * 50)

    if not taxa:
        console.print("  [dim]No novel taxa identified[/dim]")
        console.print("  [dim](All taxa in this run have been seen before)[/dim]\n")
        return

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Taxon")
    table.add_column("Samples", justify="right")
    table.add_column("Contigs", justify="right")

    for t in taxa:
        name = t.name or "-"
        if t.rank:
            name += f" [dim]({t.rank})[/dim]"
        table.add_row(
            name,
            str(t.sample_count),
            str(t.hit_count),
        )

    console.print(table)
    console.print(f"\n[dim]{len(taxa)} novel taxon/taxa[/dim]\n")


@hits_app.command("movers")
def hits_movers(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
        ),
    ],
    limit: Annotated[
        int,
        typer.Option(
            "--limit",
            "-n",
            help="Number of top movers to show",
        ),
    ] = 10,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Show top movers (biggest prevalence changes) for a run.

    Compares each taxon's prevalence in this run against its historical
    average across all prior runs. Returns taxa sorted by absolute change.

    Positive change means higher prevalence than usual; negative means lower.

    Examples:

        # Show top 10 movers (default)
        nvd hits movers --run evil_heisenberg

        # Show top 20 movers
        nvd hits movers --run evil_heisenberg --limit 20

        # Output as JSON
        nvd hits movers --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    movers = get_top_movers(sample_set_id, limit=limit)

    if json_output:
        console.print_json(data=[_taxon_change_to_dict(t) for t in movers])
    else:
        _print_movers_rich(movers, sample_set_id, limit)


def _taxon_change_to_dict(t: TaxonChange) -> dict:
    """Convert TaxonChange to JSON-serializable dict."""
    return {
        "taxid": t.taxid,
        "name": t.name,
        "rank": t.rank,
        "current_prevalence_pct": t.current_prevalence_pct,
        "baseline_prevalence_pct": t.baseline_prevalence_pct,
        "change_pct": t.change_pct,
        "current_sample_count": t.current_sample_count,
        "baseline_run_count": t.baseline_run_count,
    }


def _print_movers_rich(movers: list[TaxonChange], run: str, limit: int) -> None:
    """Print top movers with Rich formatting."""
    console.print(f"\n[bold]Top {limit} Movers:[/bold] [cyan]{run}[/cyan]")
    console.print("─" * 60)

    if not movers:
        console.print("  [dim]No movers to show[/dim]")
        console.print(
            "  [dim](This may be the first run, or no taxa in this run)[/dim]\n"
        )
        return

    # Show baseline info
    if movers:
        console.print(
            f"  [dim]Baseline: {movers[0].baseline_run_count} prior run(s)[/dim]\n"
        )

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Taxon")
    table.add_column("Current", justify="right")
    table.add_column("Baseline", justify="right")
    table.add_column("Change", justify="right")
    table.add_column("Samples", justify="right")

    for t in movers:
        name = t.name or "-"
        if t.rank:
            name += f" [dim]({t.rank})[/dim]"

        # Color the change based on direction
        if t.change_pct > 0:
            change_str = f"[green]+{t.change_pct:.1f}%[/green]"
        elif t.change_pct < 0:
            change_str = f"[red]{t.change_pct:.1f}%[/red]"
        else:
            change_str = "0.0%"

        table.add_row(
            name,
            f"{t.current_prevalence_pct:.1f}%",
            f"{t.baseline_prevalence_pct:.1f}%",
            change_str,
            str(t.current_sample_count),
        )

    console.print(table)
    console.print()


@hits_app.command("rare")
def hits_rare(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
        ),
    ],
    threshold: Annotated[
        float,
        typer.Option(
            "--threshold",
            "-t",
            help="Maximum historical run percentage to be considered rare (default: 5%)",
        ),
    ] = 5.0,
    json_output: Annotated[
        bool,
        typer.Option(
            "--json",
            help="Output as JSON",
        ),
    ] = False,
) -> None:
    """
    Show rare taxa in a run.

    Identifies taxa present in this run that have historically appeared
    in few runs (less than threshold percentage). Useful for flagging
    unusual findings that warrant closer inspection.

    Examples:

        # Show taxa appearing in <5% of runs (default)
        nvd hits rare --run evil_heisenberg

        # Show taxa appearing in <10% of runs
        nvd hits rare --run evil_heisenberg --threshold 10

        # Output as JSON
        nvd hits rare --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    if threshold <= 0 or threshold > 100:
        error("Threshold must be between 0 and 100")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    taxa = get_rare_taxa(sample_set_id, threshold_pct=threshold)

    if json_output:
        console.print_json(data=[_rare_taxon_to_dict(t) for t in taxa])
    else:
        _print_rare_taxa_rich(taxa, sample_set_id, threshold)


def _rare_taxon_to_dict(t: RareTaxon) -> dict:
    """Convert RareTaxon to JSON-serializable dict."""
    return {
        "taxid": t.taxid,
        "name": t.name,
        "rank": t.rank,
        "historical_run_count": t.historical_run_count,
        "historical_run_pct": t.historical_run_pct,
        "current_sample_count": t.current_sample_count,
    }


def _print_rare_taxa_rich(taxa: list[RareTaxon], run: str, threshold: float) -> None:
    """Print rare taxa with Rich formatting."""
    console.print(
        f"\n[bold]Rare Taxa (<{threshold:.0f}% of runs):[/bold] [cyan]{run}[/cyan]"
    )
    console.print("─" * 50)

    if not taxa:
        console.print("  [dim]No rare taxa identified[/dim]")
        console.print(f"  [dim](All taxa appear in ≥{threshold:.0f}% of runs)[/dim]\n")
        return

    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Taxon")
    table.add_column("Historical", justify="right")
    table.add_column("Samples", justify="right")

    for t in taxa:
        name = t.name or "-"
        if t.rank:
            name += f" [dim]({t.rank})[/dim]"
        table.add_row(
            name,
            f"{t.historical_run_count} runs ({t.historical_run_pct:.1f}%)",
            str(t.current_sample_count),
        )

    console.print(table)
    console.print(f"\n[dim]{len(taxa)} rare taxon/taxa[/dim]\n")


@hits_app.command("history")
def hits_history(
    taxon: Annotated[
        str,
        typer.Option(
            "--taxon",
            "-t",
            help="Taxon name to look up",
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
    Show history of a specific taxon across all runs.

    Provides a detailed track record showing when and how often a taxon
    has been detected. Useful for understanding whether a finding is
    typical or unusual for this dataset.

    Examples:

        # Show history for a taxon
        nvd hits history --taxon "Human gammaherpesvirus 4"

        # Output as JSON
        nvd hits history --taxon "Rhinovirus A" --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    history = get_taxon_history(taxon)

    if history is None:
        if json_output:
            console.print_json(data={"found": False, "taxon": taxon})
        else:
            info(f"Taxon not found: {taxon}")
        return

    if json_output:
        console.print_json(data=_taxon_history_to_dict(history))
    else:
        _print_taxon_history_rich(history)


def _taxon_history_to_dict(h: TaxonHistory) -> dict:
    """Convert TaxonHistory to JSON-serializable dict."""
    return {
        "found": True,
        "taxid": h.taxid,
        "name": h.name,
        "rank": h.rank,
        "total_runs_seen": h.total_runs_seen,
        "total_runs_in_system": h.total_runs_in_system,
        "overall_prevalence_pct": h.overall_prevalence_pct,
        "first_seen_date": h.first_seen_date,
        "last_seen_date": h.last_seen_date,
        "run_history": [
            {
                "sample_set_id": r.sample_set_id,
                "run_date": r.run_date,
                "sample_count": r.sample_count,
                "prevalence_pct": r.prevalence_pct,
            }
            for r in h.run_history
        ],
    }


def _print_taxon_history_rich(history: TaxonHistory) -> None:
    """Print taxon history with Rich formatting."""
    console.print(f"\n[bold]Taxon History:[/bold] [cyan]{history.name}[/cyan]")
    if history.rank:
        console.print(f"[dim]Rank: {history.rank}[/dim]")
    if history.taxid:
        console.print(f"[dim]Taxid: {history.taxid}[/dim]")
    console.print("─" * 50)

    # Summary
    console.print(
        f"  Runs with taxon:       {history.total_runs_seen} / {history.total_runs_in_system}"
    )
    console.print(
        f"  Overall prevalence:    {history.overall_prevalence_pct:.1f}% of runs"
    )
    console.print(f"  First seen:            {history.first_seen_date[:10]}")
    console.print(f"  Last seen:             {history.last_seen_date[:10]}")

    # Run history table
    console.print(f"\n[bold]Run History ({len(history.run_history)} runs):[/bold]")

    if history.run_history:
        table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
        table.add_column("Date", style="dim")
        table.add_column("Run")
        table.add_column("Samples", justify="right")
        table.add_column("Prevalence", justify="right")

        for r in history.run_history:
            date_str = r.run_date[:10] if len(r.run_date) >= 10 else r.run_date
            table.add_row(
                date_str,
                r.sample_set_id,
                str(r.sample_count),
                f"{r.prevalence_pct:.1f}%",
            )

        console.print(table)

    console.print()


@hits_app.command("compare")
def hits_compare(
    run: Annotated[
        str,
        typer.Option(
            "--run",
            "-r",
            help="Run ID (e.g., 'evil_heisenberg') or sample set ID",
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
    Compare run metrics to historical norms.

    Provides context for interpreting a run by comparing its key metrics
    (sample count, hit count, taxa count, contig length) against the
    historical average across all prior runs.

    Examples:

        # Compare a run to historical norms
        nvd hits compare --run evil_heisenberg

        # Output as JSON
        nvd hits compare --run evil_heisenberg --json
    """
    db_path = get_state_db_path()
    if not db_path.exists():
        error(f"Database not found: {db_path}\nRun 'nvd state init' to create it.")
        raise typer.Exit(1)

    sample_set_id = resolve_run_identifier(run)
    comparison = get_run_comparison(sample_set_id)

    if comparison is None:
        if json_output:
            console.print_json(data={"found": False, "sample_set_id": sample_set_id})
        else:
            info(f"Run not found: {run}")
        return

    if json_output:
        console.print_json(data=_run_comparison_to_dict(comparison))
    else:
        _print_run_comparison_rich(comparison)


def _run_comparison_to_dict(c: RunComparison) -> dict:
    """Convert RunComparison to JSON-serializable dict."""
    return {
        "found": True,
        "sample_set_id": c.sample_set_id,
        "run_date": c.run_date,
        "samples_analyzed": c.samples_analyzed,
        "unique_hits": c.unique_hits,
        "taxa_found": c.taxa_found,
        "median_contig_length": c.median_contig_length,
        "historical_run_count": c.historical_run_count,
        "avg_samples": c.avg_samples,
        "avg_hits": c.avg_hits,
        "avg_taxa": c.avg_taxa,
        "avg_median_length": c.avg_median_length,
    }


def _print_run_comparison_rich(comparison: RunComparison) -> None:
    """Print run comparison with Rich formatting."""
    console.print(
        f"\n[bold]Run Comparison:[/bold] [cyan]{comparison.sample_set_id}[/cyan]"
    )
    if comparison.run_date:
        date_str = (
            comparison.run_date[:10]
            if len(comparison.run_date) >= 10
            else comparison.run_date
        )
        console.print(f"[dim]Date: {date_str}[/dim]")
    console.print("─" * 50)

    if comparison.historical_run_count == 0:
        console.print("  [dim]No prior runs for comparison[/dim]\n")
        console.print("[bold]Current Run Metrics[/bold]")
        console.print(f"  Samples analyzed:      {comparison.samples_analyzed:,}")
        console.print(f"  Unique hits:           {comparison.unique_hits:,}")
        console.print(f"  Taxa found:            {comparison.taxa_found:,}")
        console.print(
            f"  Median contig length:  {comparison.median_contig_length:,} bp"
        )
        console.print()
        return

    console.print(
        f"  [dim]Baseline: {comparison.historical_run_count} prior run(s)[/dim]\n"
    )

    # Build comparison table
    table = Table(show_header=True, header_style="bold", box=None, padding=(0, 2))
    table.add_column("Metric")
    table.add_column("This Run", justify="right")
    table.add_column("Avg", justify="right")
    table.add_column("Diff", justify="right")

    def _format_diff(current: float, avg: float) -> str:
        """Format difference with color."""
        if avg == 0:
            return "-"
        diff_pct = ((current - avg) / avg) * 100
        if diff_pct > 10:
            return f"[green]+{diff_pct:.0f}%[/green]"
        elif diff_pct < -10:
            return f"[red]{diff_pct:.0f}%[/red]"
        else:
            return f"{diff_pct:+.0f}%"

    table.add_row(
        "Samples",
        f"{comparison.samples_analyzed:,}",
        f"{comparison.avg_samples:.1f}",
        _format_diff(comparison.samples_analyzed, comparison.avg_samples),
    )
    table.add_row(
        "Unique hits",
        f"{comparison.unique_hits:,}",
        f"{comparison.avg_hits:.1f}",
        _format_diff(comparison.unique_hits, comparison.avg_hits),
    )
    table.add_row(
        "Taxa found",
        f"{comparison.taxa_found:,}",
        f"{comparison.avg_taxa:.1f}",
        _format_diff(comparison.taxa_found, comparison.avg_taxa),
    )
    table.add_row(
        "Median length",
        f"{comparison.median_contig_length:,} bp",
        f"{comparison.avg_median_length:.0f} bp",
        _format_diff(comparison.median_contig_length, comparison.avg_median_length),
    )

    console.print(table)
    console.print()


def _print_stats_rich(stats: HitStats, top_recurring: list[RecurringHit]) -> None:
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
        table.add_column("Top Classification")

        for hit in top_recurring:
            table.add_row(
                hit.hit_key,
                f"{hit.sequence_length:,} bp",
                str(hit.sample_count),
                str(hit.observation_count),
                _format_classification(hit.top_taxid_name, hit.top_taxid_rank),
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

    # Classification
    classification = _format_classification(
        hit.top_taxid_name, hit.top_taxid_rank, hit.top_taxid, include_taxid=True
    )
    if hit.top_taxid_name:
        console.print(f"  Classification: {classification}")
    else:
        console.print("  Classification: [dim]None[/dim]")

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


def _print_recur_rich(
    hits: list[RecurringHit], min_samples: int, min_runs: int | None
) -> None:
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
    table.add_column("Top Classification")
    table.add_column("First Seen", style="dim")
    table.add_column("Last Seen", style="dim")

    for hit in hits:
        table.add_row(
            hit.hit_key,
            f"{hit.sequence_length:,} bp",
            str(hit.sample_count),
            str(hit.run_count),
            str(hit.observation_count),
            _format_classification(hit.top_taxid_name, hit.top_taxid_rank),
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
        "top_taxid": [],
        "top_taxid_name": [],
        "top_taxid_rank": [],
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
        data["top_taxid"].append(hit.top_taxid)
        data["top_taxid_name"].append(hit.top_taxid_name)
        data["top_taxid_rank"].append(hit.top_taxid_rank)
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
    try:
        result = compact_hits(month=month, dry_run=dry_run, keep_source=keep_source)
    except RuntimeError as e:
        error(str(e))
        raise typer.Exit(1)

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
