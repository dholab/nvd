#!/usr/bin/env python3
"""
Script to upload classified contig information to LabKey and save to CSV regardless.
"""

import sys
import pandas as pd
import more_itertools
from labkey.api_wrapper import APIWrapper
from labkey.exceptions import ServerContextError, RequestError
import os
import gzip

def count_total_reads(fastq_file):
    """
    Count the number of reads in a gzipped FASTQ file. Each read is represented by 4 lines in the FASTQ format.
    """
    read_count = 0
    with gzip.open(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:  # Each read starts on every 4th line
                read_count += 1
    return read_count

def insert_blast_results(experiment, blast_results, mapped_reads_file, stat_db_version, sample_name, blast_db_name, snakemake_run_id, labkey_server, project_name, api_key, out_dir, token_file, fastq_file):
    """
    Insert records into LabKey and save to CSV regardless.
    """
    try:
        # Check if the blast_results file is empty
        if os.stat(blast_results).st_size == 0:
            print(f"Warning: The file {blast_results} is empty.", file=sys.stderr)
            return
        
        # Read the BLAST results CSV file
        df = pd.read_csv(blast_results, sep='\t')
        
        if df.empty:
            print(f"Warning: No data found in the file {blast_results}.", file=sys.stderr)
            return

        # Read the mapped reads file into a dictionary
        mapped_reads_dict = {}
        if os.path.exists(mapped_reads_file) and os.stat(mapped_reads_file).st_size > 0:
            with open(mapped_reads_file, 'r') as f:
                for line in f:
                    ref_name, mapped_count = line.strip().split()
                    mapped_reads_dict[ref_name] = int(mapped_count)
        else:
            print(f"Warning: The mapped reads file {mapped_reads_file} is missing or empty.", file=sys.stderr)
        
        # Count total reads in the input FASTQ file
        total_reads = count_total_reads(fastq_file)

        blast_rows = []
        
        for r in df.itertuples(index=False):
            try:
                normalized_qseqid = r.qseqid.strip()  # Normalize qseqid
                mapped_reads_value = mapped_reads_dict.get(normalized_qseqid, 0)
                if mapped_reads_value == 0:
                    print(f"Warning: No mapped reads found for qseqid '{normalized_qseqid}'", file=sys.stderr)
                
                blast_row = {
                    'experiment': str(experiment),
                    'blast_task': str(r.task),
                    'sample_id': str(r.sample),
                    'qseqid': normalized_qseqid,  
                    'qlen': int(r.qlen),
                    'sseqid': str(r.sseqid),
                    'stitle': str(r.stitle),
                    'tax_rank': str(r.rank),
                    'length': int(r.length),
                    'pident': float(r.pident),
                    'evalue': float(r.evalue),
                    'bitscore': float(r.bitscore),
                    'sscinames': str(r.sscinames),
                    'blast_db_version': str(blast_db_name),
                    'snakemake_run_id': str(snakemake_run_id),
                    'mapped_reads': mapped_reads_value,
                    'total_reads': total_reads,
                    'stat_db_version': str(stat_db_version)
                }
                blast_rows.append(blast_row)
            except AttributeError as e:
                print(f"AttributeError for row: {r}", file=sys.stderr)
                print(f"Error details: {str(e)}", file=sys.stderr)
            except ValueError as e:
                print(f"ValueError for row: {r}", file=sys.stderr)
                print(f"Error details: {str(e)}", file=sys.stderr)

        # Always save results to CSV
        if not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        output_csv_file = os.path.join(out_dir, f"{sample_name}.csv")
        pd.DataFrame(blast_rows).to_csv(output_csv_file, index=False)
        print(f"Results saved to CSV: {output_csv_file}", file=sys.stderr)
        
        # Insert into LabKey if API key is present
        if api_key:
            blast_row_chunks = list(more_itertools.chunked(blast_rows, 1000))
            api = APIWrapper(labkey_server, project_name, api_key=api_key, use_ssl=True)
            for counter, chunk in enumerate(blast_row_chunks):
                try:
                    api.query.insert_rows(schema_name='lists', query_name='metagenomic_hits', rows=chunk)
                    number_processed = (counter + 1) * len(chunk)
                    print(f"{number_processed} BLAST output rows added for sample {sample_name}", file=sys.stderr)
                except (ServerContextError, RequestError) as e:
                    print(f"Error inserting chunk {counter+1}: {str(e)}", file=sys.stderr)
                    print(f"Server response: {e.response.text if hasattr(e, 'response') else 'No response text'}", file=sys.stderr)
                    print(f"Continuing with next chunk...", file=sys.stderr)

        # Create the token file to signal completion
        with open(token_file, 'w') as f:
            f.write('Completed')
        print(f"Token file created at: {token_file}", file=sys.stderr)

    except Exception as e:
        print(f"Error in insert_blast_results: {str(e)}", file=sys.stderr)
        raise

if __name__ == "__main__":
    insert_blast_results(
        experiment=snakemake.params.experiment,
        blast_results=snakemake.input.blast_results,
        mapped_reads_file=snakemake.input.mapped_reads,
        stat_db_version=snakemake.params.stat_db_version,
        sample_name=snakemake.wildcards.sample,
        blast_db_name=snakemake.params.blast_db_name,
        snakemake_run_id=snakemake.params.snakemake_run_id,
        labkey_server=snakemake.params.labkey_server,
        project_name=snakemake.params.project_name,
        api_key=snakemake.params.api_key,
        out_dir=snakemake.params.out_dir,
        token_file=snakemake.output.token,
        fastq_file=snakemake.input.fastq_file  # Pass the FASTQ file to count total reads
    )
