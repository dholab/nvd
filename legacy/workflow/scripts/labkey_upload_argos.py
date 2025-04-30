#!/usr/bin/env python3

"""
Script to aggregate results and upload to LabKey.
"""

import sys
import os
import pandas as pd
import gzip
from Bio import SeqIO
from labkey.api_wrapper import APIWrapper
from labkey.exceptions import ServerContextError, RequestError
import more_itertools

def main():
    # Get inputs and outputs from snakemake
    gc_files = snakemake.input.gc_files
    mapped_reads_files = snakemake.input.mapped_reads_files
    bed_files = snakemake.input.bed_files
    total_reads_files = snakemake.input.total_reads_files
    filtered_fasta_files = snakemake.input.filtered_fasta_files
    output_csv = snakemake.output.aggregated_csv
    experiment = snakemake.params.experiment
    argos_version = snakemake.params.argos_version
    snakemake_run_id = snakemake.params.snakemake_run_id
    labkey_server = snakemake.params.labkey_server
    project_name = snakemake.params.project_name
    api_key = snakemake.params.api_key

    # Initialize an empty list to collect data
    aggregated_data = []

    # List of sample names
    samples = [s['name'] for s in snakemake.config["samples"]]

    # Map sample names to input files
    sample_gc = dict(zip(samples, gc_files))
    sample_mapped_reads = dict(zip(samples, mapped_reads_files))
    sample_bed = dict(zip(samples, bed_files))
    sample_total_reads = dict(zip(samples, total_reads_files))
    sample_fasta = dict(zip(samples, filtered_fasta_files))

    for sample in samples:
        print(f"Processing sample: {sample}")

        # Get total_reads from total_reads file
        total_reads_file = sample_total_reads[sample]
        try:
            with open(total_reads_file, 'r') as f:
                total_reads_str = f.read().strip()
                if total_reads_str:
                    total_reads = int(float(total_reads_str))
                else:
                    total_reads = 0
                    print(f"Warning: Total reads file {total_reads_file} is empty for sample {sample}. Setting to 0.")
        except (FileNotFoundError, ValueError):
            total_reads = 0
            print(f"Warning: Unable to read total_reads for sample {sample}. Setting to 0.")

        # Read genome coverage
        gc_file = sample_gc[sample]
        try:
            gc_df = pd.read_csv(
                gc_file, sep=' ', header=None, names=['argos_id', 'genome_coverage']
            )
        except pd.errors.EmptyDataError:
            gc_df = pd.DataFrame(columns=['argos_id', 'genome_coverage'])
            print(f"Warning: Genome coverage file {gc_file} is empty for sample {sample}.")

        # Read mapped reads
        mapped_reads_file = sample_mapped_reads[sample]
        try:
            mapped_reads_df = pd.read_csv(
                mapped_reads_file, sep=' ', header=None, names=['argos_id', 'mapped_reads']
            )
        except pd.errors.EmptyDataError:
            mapped_reads_df = pd.DataFrame(columns=['argos_id', 'mapped_reads'])
            print(f"Warning: Mapped reads file {mapped_reads_file} is empty for sample {sample}.")

        # Read BED file to get contig lengths
        bed_file = sample_bed[sample]
        try:
            bed_df = pd.read_csv(
                bed_file, sep='\t', header=None, names=['argos_id', 'start', 'end']
            )
            bed_df['argos_contig_length'] = bed_df['end'] - bed_df['start']
        except pd.errors.EmptyDataError:
            bed_df = pd.DataFrame(columns=['argos_id', 'start', 'end', 'argos_contig_length'])
            print(f"Warning: BED file {bed_file} is empty for sample {sample}.")

        # Read filtered_fasta to get argos_contig mapping
        fasta_file = sample_fasta[sample]
        contig_to_desc = {}
        try:
            with gzip.open(fasta_file, 'rt') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    header = record.description  # full header
                    # Split the header into argos_id and argos_contig
                    parts = header.split(None, 1)  # Split on first whitespace
                    if len(parts) == 2:
                        argos_id, argos_contig_desc = parts
                    elif len(parts) == 1:
                        argos_id = parts[0]
                        argos_contig_desc = parts[0]
                    else:
                        # Handle cases with no header
                        continue
                    contig_to_desc[argos_id] = argos_contig_desc
        except FileNotFoundError:
            print(f"Warning: FASTA file {fasta_file} not found for sample {sample}.")
        except gzip.BadGzipFile:
            print(f"Warning: FASTA file {fasta_file} is not a valid gzip file for sample {sample}.")

        # Assign 'argos_contig' based on 'argos_id'
        if not mapped_reads_df.empty:
            # Assign argos_contig using the mapping
            mapped_reads_df['argos_contig'] = mapped_reads_df['argos_id'].map(contig_to_desc)
        else:
            mapped_reads_df['argos_contig'] = pd.Series(dtype='object')

        # Merge dataframes on 'argos_id'
        df = pd.merge(gc_df, mapped_reads_df, on='argos_id', how='outer')
        df = pd.merge(df, bed_df[['argos_id', 'argos_contig_length']], on='argos_id', how='outer')

        # Fill missing values with zeros
        df['genome_coverage'] = df['genome_coverage'].fillna(0).astype(int)
        df['mapped_reads'] = df['mapped_reads'].fillna(0).astype(int)
        df['argos_contig_length'] = df['argos_contig_length'].fillna(0).astype(int)

        # Add additional columns
        df['sample_id'] = sample
        df['experiment'] = experiment
        df['argos_version'] = argos_version
        df['snakemake_run_id'] = snakemake_run_id
        df['total_reads'] = total_reads

        # Reorder columns
        df = df[[
            'experiment',
            'sample_id',
            'argos_contig',
            'argos_id',
            'genome_coverage',
            'mapped_reads',
            'total_reads',
            'argos_contig_length',
            'argos_version',
            'snakemake_run_id'
        ]]

        # Append to aggregated_data
        aggregated_data.append(df)

    if aggregated_data:
        # Concatenate all dataframes
        aggregated_df = pd.concat(aggregated_data, ignore_index=True)
    else:
        # Create an empty DataFrame with the required columns if no data is present
        aggregated_df = pd.DataFrame(columns=[
            'experiment',
            'sample_id',
            'argos_contig',
            'argos_id',
            'genome_coverage',
            'mapped_reads',
            'total_reads',
            'argos_contig_length',
            'argos_version',
            'snakemake_run_id'
        ])

    # Drop rows where 'argos_id' or 'argos_contig' is missing
    aggregated_df = aggregated_df.dropna(subset=['argos_id', 'argos_contig'])

    # Ensure 'argos_id' and 'argos_contig' are strings
    aggregated_df['argos_id'] = aggregated_df['argos_id'].astype(str)
    aggregated_df['argos_contig'] = aggregated_df['argos_contig'].astype(str)

    # Remove duplicate rows
    aggregated_df = aggregated_df.drop_duplicates()

    # Sort by sample_id, then by largest number of mapped reads
    aggregated_df = aggregated_df.sort_values(['sample_id', 'mapped_reads'], ascending=[True, False])

    # Write to CSV
    aggregated_df.to_csv(output_csv, index=False)
    print(f"Aggregated results saved to {output_csv}", file=sys.stderr)

    # Now, upload to LabKey
    # Only proceed if API key is provided
    if api_key:
        # Prepare data for LabKey
        labkey_rows = aggregated_df.to_dict(orient='records')
        # Use APIWrapper to insert rows into 'metagenomics_argos' list
        try:
            api = APIWrapper(labkey_server, project_name, api_key=api_key, use_ssl=True)

            # Chunk the data if necessary
            labkey_row_chunks = list(more_itertools.chunked(labkey_rows, 1000))
            for counter, chunk in enumerate(labkey_row_chunks):
                try:
                    api.query.insert_rows(schema_name='lists', query_name='metagenomics_argos', rows=chunk)
                    number_processed = (counter + 1) * len(chunk)
                    print(f"{number_processed} aggregated results rows added to LabKey")
                except (ServerContextError, RequestError) as e:
                    print(f"Error inserting chunk {counter+1}: {str(e)}", file=sys.stderr)
                    print(f"Server response: {e.response.text if hasattr(e, 'response') else 'No response text'}", file=sys.stderr)
                    print(f"Continuing with next chunk...", file=sys.stderr)
        except Exception as e:
            print(f"Error uploading to LabKey: {str(e)}", file=sys.stderr)
            raise

if __name__ == "__main__":
    main()