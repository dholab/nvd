#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "biopython",
#     "snakemake",
# ]
# ///


"""
Script to extract contigs that are not in either megablast or blastn results.
"""

import argparse
import subprocess
import sys
from typing import Literal

from Bio import SeqIO

# snakemake setup
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"


def get_classified_nodes(blast_results):
    """
    Extract classified node IDs from BLAST results file.
    """
    try:
        with open(blast_results) as f:
            # Extract the first column (node ID) from each line
            return set(line.split("\t")[0] for line in f)
    except Exception as e:
        print(f"Error processing BLAST results file {blast_results}: {e}", file=sys.stderr)
        return set()


def get_unclassified_contigs(megablast_nodes, blastn_nodes, fasta_contigs, output_file):
    """
    Identify and extract unclassified contigs.
    """
    assembly_contigs = [record.id for record in SeqIO.parse(fasta_contigs, "fasta")]
    classified_nodes = megablast_nodes | blastn_nodes
    unclassified_nodes = set(assembly_contigs) - classified_nodes

    if unclassified_nodes:
        unclassified_node_names = ",".join(unclassified_nodes)
        cmd = [
            "filterbyname.sh",
            "include=t",
            f"in={fasta_contigs}",
            f"names={unclassified_node_names}",
            f"out={output_file}",
            "ow=t",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running filterbyname.sh: {e}", file=sys.stderr)
            print(f"Command output: {e.output}", file=sys.stderr)
            sys.exit(1)
    else:
        # Create an empty output file if no unclassified contigs
        open(output_file, "w").close()


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify contigs by blast result presence")
    parser.add_argument("--contigs", required=True, help="Input contigs FASTA file")
    parser.add_argument("--megablast_results", required=True, help="Input megablast results file")
    parser.add_argument("--blastn_results", required=True, help="Input blastn results file")
    parser.add_argument(
        "--unclassified",
        required=True,
        help="Output file for unclassified contigs",
    )

    return parser.parse_args()


def main():
    """
    Main function to process inputs and extract unclassified contigs.
    """
    if MODE == "snakemake":
        input_contigs = snakemake.input.contigs
        input_megablast = snakemake.input.megablast_results
        input_blastn = snakemake.input.blastn_results
        output_file = snakemake.output.unclassified
    else:
        args = parse_command_line_args()
        input_contigs = args.contigs
        input_megablast = args.megablast_results
        input_blastn = args.blastn_results
        output_file = args.unclassified

    # Create output file in case input files are empty
    open(output_file, "w").close()

    megablast_nodes = get_classified_nodes(input_megablast)
    blastn_nodes = get_classified_nodes(input_blastn)

    if megablast_nodes or blastn_nodes:
        get_unclassified_contigs(megablast_nodes, blastn_nodes, input_contigs, output_file)
    else:
        print("Both BLAST results files are empty. All contigs are unclassified.", file=sys.stderr)
        # Copy all contigs to output file as they are all unclassified
        cmd = ["cp", input_contigs, output_file]
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
