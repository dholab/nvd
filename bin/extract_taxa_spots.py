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
from pathlib import Path
from typing import Literal

from py_nvd import taxonomy

__version__ = "1.8"
logger = logging.getLogger("taxonomic_read_filter")

# snakemake setup, if applicable
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"


def resolve_taxa(taxa_list: list[str], tax) -> set[int]:
    """
    Resolve the input taxa list to a set of valid taxonomic IDs.

    Looks up scientific names in the taxonomy database and returns
    corresponding taxids.

    Args:
        taxa_list: List of scientific names to resolve
        tax: TaxonomyDB instance

    Returns:
        Set of resolved taxids
    """
    resolved_taxa = set()
    for taxon_name in taxa_list:
        # Query by scientific name
        # We need to access the connection directly for name lookup
        row = tax._conn.execute(
            "SELECT tax_id FROM taxons WHERE scientific_name = ?", (taxon_name,)
        ).fetchone()
        if row:
            resolved_taxa.add(row[0])
        else:
            logger.warning(f"Taxon '{taxon_name}' not found in database.")
    return resolved_taxa


def is_in_target_lineage(
    tax_id: int, target_taxa: set[int], tax, include_children: bool
) -> bool:
    """
    Check if a given tax_id is in the target lineage based on the filtering options.

    Args:
        tax_id: The taxid to check
        target_taxa: Set of target taxids
        tax: TaxonomyDB instance
        include_children: If True, check if any ancestor is in target_taxa

    Returns:
        True if the taxid matches the filter criteria
    """
    if include_children:
        # Check if any taxid in the lineage is in target_taxa
        lineage_ids = tax.get_lineage_ids(tax_id)
        return any(tid in target_taxa for tid in lineage_ids)
    return tax_id in target_taxa


def process_input(
    input_file: str,
    output_file: str,
    target_taxa: set[int],
    tax,
    include_children: bool,
    stringency: float,
) -> None:
    """Process the input file, apply filtering, and write results to the output file."""
    total_spots = 0
    filtered_spots = 0

    try:
        with open(input_file) as infile:
            first_line = infile.readline().strip()
            if not first_line:
                logger.warning("Input file is empty. No spots to process.")
                Path(output_file).touch()
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
                            tax,
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
        f"Filtering rate: {filtered_spots / total_spots:.2%}"
        if total_spots > 0
        else "No spots processed.",
    )


def parse_command_line_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Taxa filtering pipeline")
    # Keep --gettax_sqlite_path for backward compatibility with Nextflow,
    # but it's no longer used (taxonomy.open() handles DB location)
    parser.add_argument(
        "--gettax_sqlite_path",
        required=False,
        default=None,
        help="(Deprecated) Path to SQLite taxonomy DB - now handled automatically",
    )
    parser.add_argument("--hits_file", required=True, help="Path to hits file")
    parser.add_argument(
        "--output_file", required=True, help="Where to write filtered results"
    )
    parser.add_argument(
        "--taxa", required=True, nargs="+", help="List of taxa to include"
    )
    parser.add_argument(
        "--stringency", type=float, default=0.9, help="Stringency threshold"
    )
    parser.add_argument(
        "--include_children",
        action="store_true",
        help="Include child taxa in filter",
    )

    return parser.parse_args()


def execute(
    hits_file: str,
    output_file: str,
    taxa: list[str],
    stringency: float,
    include_children: bool,
) -> None:
    """Execute the taxa filtering pipeline."""
    try:
        with taxonomy.open() as tax:
            resolved_taxa = resolve_taxa(taxa, tax)
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
                tax,
                include_children,
                stringency,
            )
    except Exception:
        logger.exception("An error occurred. Detailed traceback:")
        raise


def main() -> None:
    """Main entry point."""
    if MODE == "snakemake":
        with open(snakemake.log[0], "w") as log_file, redirect_stderr(log_file):
            logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s - %(levelname)s - %(message)s",
                stream=sys.stderr,
            )
            # gettax_sqlite_path from snakemake is ignored - taxonomy.open() handles it
            execute(
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

    # args.gettax_sqlite_path is ignored - taxonomy.open() handles DB location
    execute(
        hits_file=args.hits_file,
        output_file=args.output_file,
        taxa=args.taxa,
        stringency=args.stringency,
        include_children=args.include_children,
    )


if __name__ == "__main__":
    main()
