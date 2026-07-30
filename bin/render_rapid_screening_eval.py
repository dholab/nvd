#!/usr/bin/env python3
"""Render the rapid-screening eval report."""

from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path
from typing import TYPE_CHECKING

import duckdb

if TYPE_CHECKING:
    from collections.abc import Sequence


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def query_rows(database: Path, sql: str) -> list[dict[str, object]]:
    con = duckdb.connect(str(database), read_only=True)
    try:
        result = con.execute(sql)
        columns = [description[0] for description in result.description]
        return [dict(zip(columns, row, strict=True)) for row in result.fetchall()]
    finally:
        con.close()


def esc(value: object) -> str:
    if value is None:
        return ""
    return html.escape(str(value))


def fmt_float(value: object) -> str:
    if value in (None, ""):
        return ""
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return esc(value)


def table(rows: Sequence[dict[str, object]], columns: Sequence[str]) -> str:
    if not rows:
        return '<p class="muted">No rows.</p>'
    header = "".join(f"<th>{esc(column)}</th>" for column in columns)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{esc(row.get(column, ''))}</td>" for column in columns)
        body_rows.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def numeric_table(rows: Sequence[dict[str, object]], columns: Sequence[str]) -> str:
    if not rows:
        return '<p class="muted">No rows.</p>'
    header = "".join(f"<th>{esc(column)}</th>" for column in columns)
    body_rows = []
    for row in rows:
        cells = []
        for column in columns:
            value = row.get(column, "")
            text = fmt_float(value) if isinstance(value, float | int) else esc(value)
            cells.append(f"<td>{text}</td>")
        body_rows.append(f"<tr>{''.join(cells)}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def render_html(
    *,
    database: Path,
    followup_tsv: Path,
    without_followup_tsv: Path,
) -> str:
    samples = query_rows(
        database,
        """
        select sample_id,
               total_reads,
               has_sourmash_tax_summary,
               has_nvd_crumbs_taxa,
               sourmash_found_fraction,
               sourmash_unclassified_fraction
        from samples
        order by sample_id
        """,
    )
    bin_totals = query_rows(
        database,
        """
        select followup_bin,
               count(*) as n_signals,
               sum(primary_screening_mass) as total_primary_screening_mass
        from screening_detections
        group by followup_bin
        order by followup_bin
        """,
    )
    rank_totals = query_rows(
        database,
        """
        select rank as comparison_rank,
               followup_bin,
               count(*) as n_signals,
               sum(primary_screening_mass) as total_primary_screening_mass
        from screening_detections
        group by rank, followup_bin
        order by rank, followup_bin
        """,
    )
    followup_rows = read_tsv(followup_tsv)
    without_followup_rows = read_tsv(without_followup_tsv)[:100]
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Rapid screening eval</title>
  <style>
    body {{ margin: 0; font: 15px/1.5 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17202a; background: #f8fafc; }}
    header {{ padding: 2rem; background: #ffffff; border-bottom: 1px solid #d9e2ec; }}
    main {{ max-width: 1180px; margin: 0 auto; padding: 1.5rem; }}
    h1 {{ margin: 0 0 .5rem; font-size: clamp(1.8rem, 4vw, 3rem); letter-spacing: -0.04em; }}
    h2 {{ margin-top: 2rem; }}
    .subtitle, .muted {{ color: #5f6f83; }}
    .card {{ background: #fff; border: 1px solid #d9e2ec; border-radius: 14px; padding: 1rem; margin: 1rem 0; box-shadow: 0 1px 2px rgba(15, 23, 42, .04); }}
    table {{ border-collapse: collapse; width: 100%; background: #fff; margin: .75rem 0; }}
    th, td {{ border: 1px solid #d9e2ec; padding: .45rem .55rem; vertical-align: top; }}
    th {{ text-align: left; background: #f1f5f9; }}
    code {{ background: #eef2f7; border-radius: 4px; padding: .05rem .25rem; }}
  </style>
</head>
<body>
  <header>
    <h1>Rapid screening eval</h1>
    <p class="subtitle">The pipeline performs rapid k-mer-based taxonomic decomposition before assembly to surface potentially important signals while the more expensive downstream analysis is still running. This experimental eval records how those preliminary sourmash signals relate to later assembly-, alignment-, and mapback-based evidence, supporting continuous calibration of the screening stage independently of the particular screening engine.</p>
    <p class="muted">This is not a reciprocal precision/recall report, and later pipeline evidence is not treated as perfect ground truth.</p>
  </header>
  <main>
    <section class="card">
      <h2>Run/sample status</h2>
      {numeric_table(samples, ["sample_id", "total_reads", "has_sourmash_tax_summary", "has_nvd_crumbs_taxa", "sourmash_found_fraction", "sourmash_unclassified_fraction"])}
    </section>

    <section class="card">
      <h2>Overall follow-up bins</h2>
      {numeric_table(bin_totals, ["followup_bin", "n_signals", "total_primary_screening_mass"])}
    </section>

    <section class="card">
      <h2>Follow-up by rank</h2>
      {numeric_table(rank_totals, ["comparison_rank", "followup_bin", "n_signals", "total_primary_screening_mass"])}
    </section>

    <section class="card">
      <h2>Follow-up by sample and rank</h2>
      <p class="muted">Feeds <code>screening_signal_followup_by_sample_rank.tsv</code>.</p>
      {table(followup_rows, ["sample_id", "screening_engine", "comparison_rank", "followup_bin", "n_signals", "total_primary_screening_mass"])}
    </section>

    <section class="card">
      <h2>Signals without same-rank follow-up evidence</h2>
      <p class="muted">Top rows from <code>screening_signals_without_same_rank_followup.tsv</code>, sorted by raw primary screening mass. Absence of follow-up evidence does not prove that a screening signal is false; it means later NVD assembly, BLAST/LCA, and CRUMBS evidence did not report the taxon at the same rank.</p>
      {table(without_followup_rows, ["sample_id", "screening_engine", "comparison_rank", "screening_taxon_name", "primary_screening_mass", "sourmash_fraction", "sourmash_f_weighted_at_rank", "lineage_source", "screening_lineage"])}
    </section>
  </main>
</body>
</html>
"""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render the rapid-screening eval HTML")
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--followup-tsv", type=Path, required=True)
    parser.add_argument("--without-followup-tsv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    args.output.write_text(
        render_html(
            database=args.database,
            followup_tsv=args.followup_tsv,
            without_followup_tsv=args.without_followup_tsv,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
