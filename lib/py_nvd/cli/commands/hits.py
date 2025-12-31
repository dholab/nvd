"""
CLI commands for exporting hits from the state database.

Commands:
    nvd hits export  - Export hits to various formats
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from py_nvd.cli.utils import console, error, info
from py_nvd.db import get_state_db_path

hits_app = typer.Typer(
    name="hits",
    help="Export and query metagenomic hits",
    no_args_is_help=True,
)


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
    import polars as pl

    from py_nvd.hits import get_hit_sequence, list_hits_with_observations

    rows = list_hits_with_observations()

    if not rows:
        info("No hits found in database")
        return

    # Build data for DataFrame
    data: dict[str, list] = {
        "hit_key": [],
        "sequence_length": [],
        "gc_content": [],
        "first_seen_date": [],
        "sample_set_id": [],
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
    from py_nvd.hits import get_hit_sequence, list_hits

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
