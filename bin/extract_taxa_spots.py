#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "snakemake",
# ]
# ///
"""
Taxonomic Read Filter for Snakemake
This script filters reads based on taxonomic classifications, specifically tailored for
identifying reads potentially originating from human viruses. It includes (whitelist)
reads that match specified taxonomic groups, with adjustable stringency.
This script is designed to be used with Snakemake's script directive.
"""

import argparse
import logging
import sys
from collections import Counter
from contextlib import redirect_stderr
from typing import Literal

from hits_to_report import TaxonomyDatabase, ensure_taxonomy_file, get_lineage

__version__ = "1.7"
logger = logging.getLogger("taxonomic_read_filter")

# snakemake setup, if applicable
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"


def resolve_taxa(taxa_list, tax_db):
    """Resolve the input taxa list to a set of valid taxonomic IDs."""
    resolved_taxa = set()
    for taxon in taxa_list:
        cur = tax_db.conn.cursor()
        cur.execute("SELECT tax_id FROM taxons WHERE scientific_name = ?", (taxon,))
        result = cur.fetchone()
        if result:
            resolved_taxa.add(result[0])
        else:
            logger.warning(f"Taxon '{taxon}' not found in database.")
    return resolved_taxa


def is_in_target_lineage(tax_id, target_taxa, lineage_cache, tax_db, include_children):
    """Check if a given tax_id is in the target lineage based on the filtering options."""
    lineage = get_lineage(tax_id, lineage_cache, tax_db)

    if include_children:
        return any(taxon in target_taxa for taxon in lineage)
    return tax_id in target_taxa


def process_input(input_file, output_file, target_taxa, tax_db, include_children, stringency):
    """Process the input file, apply filtering, and write results to the output file."""
    lineage_cache = {}
    total_spots = 0
    filtered_spots = 0

    try:
        with open(input_file) as infile:
            first_line = infile.readline().strip()
            if not first_line:
                logger.warning("Input file is empty. No spots to process.")
                return

            # Reset file pointer to the beginning of the file
            infile.seek(0)

            with open(output_file, "w") as outfile:
                for line in infile:
                    total_spots += 1
                    parts = line.strip().split("\t")
                    spot = parts[0]

                    tax_counts = Counter()
                    for hit in parts[1:]:
                        tax_id, count = hit.split("x") if "x" in hit else (hit, "1")
                        tax_id = int(tax_id)
                        count = int(count)
                        tax_counts[tax_id] += count

                    total_hits = sum(tax_counts.values())
                    target_hits = sum(
                        count
                        for tax_id, count in tax_counts.items()
                        if is_in_target_lineage(
                            tax_id,
                            target_taxa,
                            lineage_cache,
                            tax_db,
                            include_children,
                        )
                    )

                    if total_hits == 0:
                        continue

                    target_ratio = target_hits / total_hits

                    if target_ratio >= stringency:
                        filtered_spots += 1
                        outfile.write(f"{spot}\n")

    except FileNotFoundError:
        logger.exception(f"Input file {input_file} not found. Exiting.")
        return

    logger.info(f"Total spots processed: {total_spots}")
    logger.info(f"Spots after filtering: {filtered_spots}")
    logger.info(
        f"Filtering rate: {filtered_spots/total_spots:.2%}"
        if total_spots > 0
        else "No spots processed.",
    )


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Taxa filtering pipeline")
    parser.add_argument("--gettax_sqlite_path", required=True, help="Path to SQLite taxonomy DB")
    parser.add_argument("--hits_file", required=True, help="Path to hits file")
    parser.add_argument("--output_file", required=True, help="Where to write filtered results")
    parser.add_argument("--taxa", required=True, nargs="+", help="List of taxa to include")
    parser.add_argument("--stringency", type=float, default=0.9, help="Stringency threshold")
    parser.add_argument(
        "--include_children",
        action="store_true",
        help="Include child taxa in filter",
    )

    return parser.parse_args()


def execute(
    gettax_sqlite_path: str,
    hits_file: str,
    output_file: str,
    taxa: list[str],
    stringency: float,
    include_children: bool,
) -> None:
    try:
        logger.info(f"Using gettax_sqlite_path: {gettax_sqlite_path}")

        sqlite_cache = ensure_taxonomy_file(gettax_sqlite_path)

        with TaxonomyDatabase(sqlite_cache) as tax_db:
            resolved_taxa = resolve_taxa(taxa, tax_db)
            if not resolved_taxa:
                logger.error("No valid taxa specified. Exiting.")
                msg = "No valid taxa specified"
                raise ValueError(msg)

            logger.info(f"Filtering based on the following taxa: {resolved_taxa}")
            logger.info(f"Using stringency threshold: {stringency}")

            process_input(
                hits_file,
                output_file,
                resolved_taxa,
                tax_db,
                include_children,
                stringency,
            )
    except Exception:
        logger.exception("An error occurred. Detailed traceback:")
        raise


def main() -> None:
    if MODE == "snakemake":
        with open(snakemake.log[0], "w") as log_file, redirect_stderr(log_file):
            logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s - %(levelname)s - %(message)s",
                stream=sys.stderr,
            )
            execute(
                gettax_sqlite_path=snakemake.config["global"]["gettax_sqlite_path"],
                hits_file=snakemake.input.hits,
                output_file=snakemake.output.filtered,
                taxa=snakemake.params.taxa,
                stringency=snakemake.params.stringency,
                include_children=snakemake.params.include_children,
            )
        return

    args = parse_command_line_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        stream=sys.stderr,
    )

    execute(
        gettax_sqlite_path=args.gettax_sqlite_path,
        hits_file=args.hits_file,
        output_file=args.output_file,
        taxa=args.taxa,
        stringency=args.stringency,
        include_children=args.include_children,
    )


if __name__ == "__main__":
    main()
