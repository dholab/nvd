#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "snakemake",
# ]
# ///

import argparse
from collections import defaultdict
from contextlib import redirect_stderr
from pathlib import Path
from typing import Literal

# snakemake setup, if applicable
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"


def parse_taxids(hits_file: str | Path, output_file: str | Path) -> None:
    # Dictionary to store the count of each tax_id
    tax_hits = defaultdict(int)

    # Open and read the hits file
    with open(hits_file) as f:
        for line in f:
            _spot_id, spot_hits = line.strip().split("\t", 1)
            for tax_id in spot_hits.split():
                tax_id_int = int(tax_id.split("x")[0])
                tax_hits[tax_id_int] += 1

    # Write the unique list of taxids to the output file
    with open(output_file, "w") as tax_list:
        for tax_id in tax_hits:
            tax_list.write(f"{tax_id}\n")
        tax_list.flush()


def main() -> None:
    # Check if the script is being run by Snakemake or from the command line
    if "snakemake" in globals() and MODE == "snakemake":
        # Snakemake provides the `snakemake` object
        hits_file = str(snakemake.input.first_pass_stat_file)
        output_file = str(snakemake.output.tax_list)
        log_file = str(snakemake.log[0])

        # Redirect stderr to the log file
        with open(log_file, "w") as log, redirect_stderr(log):
            parse_taxids(hits_file, output_file)

        # early return to reduce nesting below
        return

    # Set up argument parsing for command-line usage
    parser = argparse.ArgumentParser(description="Parse taxids from hits file")
    parser.add_argument("hits_file", help="Path to the hits file")
    parser.add_argument("output_file", help="Path to the output file")
    args = parser.parse_args()

    hits_file = args.hits_file
    output_file = args.output_file

    # Run the main function
    parse_taxids(hits_file, output_file)


if __name__ == "__main__":
    main()
