#!/usr/bin/env python3
"""
Enrich merged BLAST results with pipeline metadata columns.

Appends per-contig mapped read counts (joined by qseqid), sample-level
total reads, BLAST database version, and workflow run identifier to each
row. The output TSV is the canonical final BLAST result — complete
regardless of whether LabKey is enabled.
"""

import argparse
import csv
import sqlite3
from pathlib import Path

MIN_CONTIG_COUNT_COLUMNS = 2


def load_contig_counts(contig_counts_path: Path) -> dict[str, str]:
    """Load per-contig mapped read counts from a two-column TSV."""
    contig_counts = {}
    if contig_counts_path.exists() and contig_counts_path.stat().st_size > 0:
        with open(contig_counts_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= MIN_CONTIG_COUNT_COLUMNS:
                    contig_counts[parts[0]] = parts[1]
    return contig_counts


def load_contig_metadata(contig_lookup_path: Path | None) -> dict[str, dict[str, str]]:
    """Load contig metadata keyed by stable BLAST query ID."""
    if contig_lookup_path is None:
        return {}
    with sqlite3.connect(contig_lookup_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(
            """
            select qseqid, evidence_class, producer, contig_id
            from contigs
            """,
        ).fetchall()
    return {row["qseqid"]: dict(row) for row in rows}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Enrich merged BLAST results with pipeline metadata.",
    )
    parser.add_argument(
        "--blast-tsv",
        type=Path,
        required=True,
        help="Input merged BLAST TSV (with LCA annotations)",
    )
    parser.add_argument(
        "--contig-counts",
        type=Path,
        required=True,
        help="Two-column TSV of contig_id and mapped read count",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output enriched BLAST TSV",
    )
    parser.add_argument(
        "--contig-lookup",
        type=Path,
        default=None,
        help="Optional per-sample SQLite lookup keyed by qseqid",
    )
    parser.add_argument(
        "--total-reads",
        required=True,
        help="Sample-level total input read count",
    )
    parser.add_argument(
        "--blast-db-version",
        required=True,
        help="BLAST database version string",
    )
    parser.add_argument(
        "--virus-index-version",
        required=True,
        help="STAT k-mer database version used for virus enrichment index",
    )
    parser.add_argument(
        "--run-id",
        required=True,
        help="Workflow run identifier",
    )
    args = parser.parse_args(argv)

    contig_counts = load_contig_counts(args.contig_counts)
    contigs_by_qseqid = load_contig_metadata(args.contig_lookup)

    metadata_columns = [
        "evidence_class",
        "producer",
        "contig_id",
        "mapped_reads",
        "total_reads",
        "blast_db_version",
        "virus_index_version",
        "nextflow_run_id",
    ]

    with open(args.blast_tsv) as infile, open(args.output, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        writer = csv.DictWriter(
            outfile,
            fieldnames=reader.fieldnames + metadata_columns,
            delimiter="\t",
        )
        writer.writeheader()
        for row in reader:
            qseqid = row.get("qseqid", "")
            contig = contigs_by_qseqid.get(qseqid, {})
            row["evidence_class"] = contig.get("evidence_class", "")
            row["producer"] = contig.get("producer", "")
            row["contig_id"] = contig.get("contig_id", "")
            row["mapped_reads"] = contig_counts.get(qseqid, "0")
            row["total_reads"] = args.total_reads
            row["blast_db_version"] = args.blast_db_version
            row["virus_index_version"] = args.virus_index_version
            row["nextflow_run_id"] = args.run_id
            writer.writerow(row)


if __name__ == "__main__":
    main()
