#!/usr/bin/env python3
"""
Script to extract contigs that are not in either megablast or blastn results.
"""

import os
import sys
import subprocess

try:
    import pandas as pd
    from Bio import SeqIO
except ImportError as e:
    print(f"Error importing required libraries: {e}", file=sys.stderr)
    print("Please ensure pandas and biopython are installed in your conda environment.", file=sys.stderr)
    sys.exit(1)

def get_classified_nodes(blast_results):
    """
    Extract classified node IDs from BLAST results file.
    """
    try:
        with open(blast_results, 'r') as f:
            # Extract the first column (node ID) from each line
            return set(line.split('\t')[0] for line in f)
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
        unclassified_node_names = ','.join(unclassified_nodes)
        cmd = [
            'filterbyname.sh', 'include=t', 
            f'in={fasta_contigs}', 
            f'names={unclassified_node_names}', 
            f'out={output_file}', 'ow=t'
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running filterbyname.sh: {e}", file=sys.stderr)
            print(f"Command output: {e.output}", file=sys.stderr)
            sys.exit(1)
    else:
        # Create an empty output file if no unclassified contigs
        open(output_file, 'w').close()

def main(input_contigs, input_megablast, input_blastn, output_file):
    """
    Main function to process inputs and extract unclassified contigs.
    """
    # Create output file in case input files are empty
    open(output_file, 'w').close()

    megablast_nodes = get_classified_nodes(input_megablast)
    blastn_nodes = get_classified_nodes(input_blastn)

    if megablast_nodes or blastn_nodes:
        get_unclassified_contigs(megablast_nodes, blastn_nodes, input_contigs, output_file)
    else:
        print("Both BLAST results files are empty. All contigs are unclassified.", file=sys.stderr)
        # Copy all contigs to output file as they are all unclassified
        cmd = ['cp', input_contigs, output_file]
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main(snakemake.input.contigs, 
         snakemake.input.megablast_results, 
         snakemake.input.blastn_results, 
         snakemake.output.unclassified)