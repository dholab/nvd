#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "pandas",
# ]
# ///

import argparse
import os

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Keep BLAST query groups with non-phage viral hits",
    )
    parser.add_argument("input_file", help="Annotated BLAST TSV")
    parser.add_argument("output_file", help="Filtered output TSV")
    return parser.parse_args()


def contains_non_phage_viruses(group: pd.DataFrame) -> bool:
    virus_hits = group[group["rank"].str.contains("superkingdom:Viruses", na=False)]
    non_phage_viruses = virus_hits[
        ~virus_hits["stitle"].str.contains("phage", case=False)
    ]
    return len(non_phage_viruses) > 0


def filter_non_virus_blast_nodes(input_file: str, output_file: str) -> None:

    # BLAST processes may emit zero-byte files when no hits are found.
    if os.path.getsize(input_file) == 0:
        # Create an empty output file if the input is empty
        open(output_file, "w").close()
        return

    # Read the input file
    all_hits_df = pd.read_csv(input_file, sep="\t")

    # Filter the dataframe
    filtered_df = all_hits_df.groupby("qseqid").filter(contains_non_phage_viruses)

    # Write the filtered data to the output file
    filtered_df.to_csv(output_file, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    filter_non_virus_blast_nodes(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
