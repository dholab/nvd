#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pandas",
#     "snakemake",
# ]
# ///

import os
import sys
from typing import Literal

import pandas as pd

# snakemake setup
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"


def usage() -> None:
    print("Usage: python3 filter_non_virus_blast_hits.py <INPUT TXT> <OUTPUT NAME>")  # noqa: T201


def check_args(args: list[str]) -> None:
    if len(args) == 0 or (len(args) == 2 and args[1] == "-h"):  # noqa: PLR2004
        usage()
        sys.exit(0)
    elif len(args) != 3:  # noqa: PLR2004
        print("Error! Incorrect argument provided.")  # noqa: T201
        usage()
        sys.exit(1)


def contains_non_phage_viruses(group: str) -> bool:
    virus_hits = group[group["rank"].str.contains("superkingdom:Viruses", na=False)]
    non_phage_viruses = virus_hits[~virus_hits["stitle"].str.contains("phage", case=False)]
    return len(non_phage_viruses) > 0


def main() -> None:
    check_args(sys.argv)

    input_file = snakemake.input[0] if MODE == "snakemake" else sys.argv[1]
    output_file = snakemake.output[1] if MODE == "snakemake" else sys.argv[2]

    # Check if the input file is empty
    # TODO(@Nick): This is necessary for snakemake, but could actually cause problems for Nextflow
    if os.path.getsize(input_file) == 0:
        # Create an empty output file if the input is empty
        open(output_file, "w").close()
        print("Input file is empty. Created an empty output file.", file=sys.stderr)  # noqa: T201
        return

    # Read the input file
    all_hits_df = pd.read_csv(input_file, sep="\t")

    # Filter the dataframe
    filtered_df = all_hits_df.groupby("qseqid").filter(contains_non_phage_viruses)

    # Write the filtered data to the output file
    filtered_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
