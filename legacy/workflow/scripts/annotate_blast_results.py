#!/usr/bin/env python3

import csv
import logging
import os
import sqlite3
import sys
from snakemake.logging import logger as snakemake_logger

__version__ = '1.7'

# Set up logging
logger = logging.getLogger(__name__)

class TaxonomyDatabase:
    def __init__(self, sqlite_cache):
        self.sqlite_cache = sqlite_cache
        self.conn = None

    def connect(self):
        self.conn = sqlite3.connect(self.sqlite_cache)
        return self

    def close(self):
        if self.conn:
            self.conn.close()

    def __enter__(self):
        return self.connect()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def get_lineage(self, tax_id):
        cur = self.conn.cursor()
        lineage = []
        while tax_id != 1:
            cur.execute('SELECT parent_tax_id, rank, scientific_name FROM taxons WHERE tax_id = ?', [int(tax_id)])
            row = cur.fetchone()
            if row is None:
                break
            parent_tax_id, rank, scientific_name = row
            if rank:
                lineage.append(f"{rank}:{scientific_name}")
            tax_id = parent_tax_id
        cur.close()
        return "; ".join(reversed(lineage))

    def check_database_structure(self):
        cur = self.conn.cursor()
        cur.execute("PRAGMA table_info(taxons)")
        columns = [column[1] for column in cur.fetchall()]
        required_columns = ['tax_id', 'parent_tax_id', 'rank', 'scientific_name']
        missing_columns = set(required_columns) - set(columns)
        if missing_columns:
            raise ValueError(f"Missing required columns in taxons table: {', '.join(missing_columns)}")

def main(snakemake):
    sqlite_cache = snakemake.params.gettax_sqlite_path
    input_file = snakemake.input.blast_results
    output_file = snakemake.output[0]
    sample_name = snakemake.params.sample
    task = snakemake.params.task  # Get task from Snakemake params

    try:
        if not os.path.exists(sqlite_cache):
            logger.error(f"SQLite cache file not found: {sqlite_cache}")
            sys.exit(1)

        if not os.path.exists(input_file):
            logger.error(f"Input file not found: {input_file}")
            sys.exit(1)

        if os.path.getsize(input_file) == 0:
            logger.warning("Input file is empty. Creating an empty output file.")
            open(output_file, 'w').close()
            return

        with TaxonomyDatabase(sqlite_cache) as tax_db:
            tax_db.check_database_structure()

            with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
                reader = csv.reader(infile, delimiter='\t')
                writer = csv.writer(outfile, delimiter='\t')

                # Write header
                header = ["task", "sample", "qseqid", "qlen", "sseqid", "stitle", "length", "pident", "evalue", "bitscore", "sscinames", "staxids", "rank"]
                writer.writerow(header)

                for row in reader:
                    qseqid, qlen, sseqid, stitle, length, pident, evalue, bitscore, sscinames, staxids = row
                    
                    # Get lineages for all taxids
                    taxids = staxids.split(';')
                    lineages = [tax_db.get_lineage(taxid) for taxid in taxids]
                    full_lineage = "; ".join(lineages)

                    # Create new row with additional information
                    new_row = [task, sample_name] + row + [full_lineage]
                    writer.writerow(new_row)

        logger.info(f"Annotated results written to {output_file}")

    except sqlite3.OperationalError as e:
        logger.error(f"SQLite error: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.error(f"Database structure error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        logger.exception("Detailed traceback:")
        sys.exit(1)

if __name__ == '__main__':
    # Use Snakemake's logger
    logger = snakemake_logger
    main(snakemake)