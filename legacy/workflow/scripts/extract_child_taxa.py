#!/usr/bin/env python3
"""
Taxa Lookup Script
This script takes a list of taxonomic names or IDs, resolves them into taxonomic IDs,
and recursively retrieves all descendant taxa (children) for each specified taxon.
It outputs the results to a file with one taxon per line.
"""
import logging
import argparse
from hits_to_report import TaxonomyDatabase, ensure_taxonomy_file

__version__ = '1.8'
logger = logging.getLogger('taxa_lookup')

def resolve_taxa(taxa_list, tax_db):
    """Resolve the input taxa list to a set of valid taxonomic IDs."""
    resolved_taxa = set()
    for taxon in taxa_list:
        cur = tax_db.conn.cursor()
        cur.execute("SELECT tax_id FROM taxons WHERE scientific_name = ? OR tax_id = ?", (taxon, taxon))
        result = cur.fetchone()
        if result:
            resolved_taxa.add(result[0])
        else:
            logger.warning(f"Taxon '{taxon}' not found in database.")
    return resolved_taxa

def get_all_descendants(tax_id, tax_db):
    """Recursively fetch all child taxonomic IDs for the given tax_id."""
    descendants = set()
    cur = tax_db.conn.cursor()
    cur.execute("SELECT tax_id FROM taxons WHERE parent_tax_id = ?", (tax_id,))
    children = cur.fetchall()
    for child in children:
        child_tax_id = child[0]
        descendants.add(child_tax_id)
        # Recursively get all descendants of this child
        descendants.update(get_all_descendants(child_tax_id, tax_db))
    return descendants

def process_taxa(taxa_list, tax_db, include_children):
    """Process taxa, optionally including all child taxa."""
    resolved_taxa = resolve_taxa(taxa_list, tax_db)
    if include_children:
        all_taxa = set(resolved_taxa)  # Start with the resolved taxa
        for taxon in resolved_taxa:
            all_taxa.update(get_all_descendants(taxon, tax_db))
        return all_taxa
    return resolved_taxa

def main():
    parser = argparse.ArgumentParser(description="Taxa Lookup")
    
    # Required arguments
    parser.add_argument("--output", required=True, help="Output file for taxa list")
    parser.add_argument("--taxa", nargs="+", required=True, help="List of taxonomic names or IDs to resolve")
    parser.add_argument("--include-children", action="store_true", help="Include all child taxa for each specified taxon")
    parser.add_argument("--gettax-sqlite-path", required=True, help="Path to gettax SQLite database file")
    parser.add_argument("--log", help="Log file (default: stdout)", default=None)

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filename=args.log if args.log else None
    )

    try:
        logger.info(f"Using gettax_sqlite_path: {args.gettax_sqlite_path}")
        sqlite_cache = ensure_taxonomy_file(args.gettax_sqlite_path)

        with TaxonomyDatabase(sqlite_cache) as tax_db:
            all_taxa = process_taxa(args.taxa, tax_db, args.include_children)

            if not all_taxa:
                logger.error("No valid taxa specified. Exiting.")
                raise ValueError("No valid taxa specified")
            
            logger.info(f"Resolved the following taxa: {all_taxa}")

            # Write the results to the output file
            with open(args.output, 'w') as outfile:
                for taxon in sorted(all_taxa):
                    outfile.write(f"{taxon}\n")

    except Exception as e:
        logger.error(f"An error occurred: {e}")
        logger.exception("Detailed traceback:")
        raise

if __name__ == "__main__":
    main()
