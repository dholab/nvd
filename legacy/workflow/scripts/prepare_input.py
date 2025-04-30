#!/usr/bin/env python3
"""
prepare_input.py

This script is used as part of a Snakemake workflow to prepare input files for further processing.
It handles various input types including paired-end reads, SRA accessions, interleaved FASTQ,
and Oxford Nanopore Technologies (ONT) FASTQ files.

The script performs different operations based on the input type:
- For paired-end reads: It runs a repair operation using repair.sh
- For SRA accessions: It downloads and processes the data using fasterq-dump
- For interleaved and ONT FASTQ: It copies the input file to the output location

Usage:
    This script is intended to be called by Snakemake and not run directly.
    It expects certain parameters to be provided by Snakemake.

Dependencies:
    - snakemake
    - repair.sh (for paired-end repair)
    - fasterq-dump (for SRA accession processing)
    - pigz (for compressing output)
"""

import os
from snakemake.shell import shell

# Extract parameters passed from Snakemake
input_type = snakemake.params.input_type
sra_accession = snakemake.params.sra_accession
temp_dir = snakemake.params.temp_dir
output_file = snakemake.output[0]
log_file = snakemake.log[0]
threads = snakemake.threads
input_files = snakemake.input.files

# Debug information
print(f"Debug: input_type = {input_type}")
print(f"Debug: input_files = {input_files}")
print(f"Debug: SRA Accession = {sra_accession}")

if input_type == "paired_end":
    if len(input_files) != 2:
        raise ValueError(f"Paired-end input requires exactly two input files. Got: {input_files}")
    command = (f"repair.sh in={input_files[0]} in2={input_files[1]} "
               f"out={output_file} threads={threads} 2> {log_file}")
    print(f"Running paired-end repair: {command}")
    shell(command)

elif input_type == "sra" and sra_accession:
    command = (
        f"fasterq-dump.3.1.1 {sra_accession} -O {temp_dir} "
        f"--split-3 --verbose --threads {threads} && "
        f"cat {temp_dir}/*.fastq | pigz -c > {output_file} 2> {log_file} && "
        f"rm -rf {temp_dir} 2>> {log_file}"
    )
    print(f"Running SRA accession processing: {command}")
    shell(command)

elif input_type in ["interleaved", "ont"]:
    if len(input_files) != 1:
        raise ValueError(f"{input_type} input requires exactly one input file. Got: {input_files}")
    if not os.path.exists(input_files[0]):
        raise FileNotFoundError(f"Input file does not exist: {input_files[0]}")
    command = f"cp {input_files[0]} {output_file} 2> {log_file}"
    print(f"Copying {input_type}: {command}")
    shell(command)

else:
    raise ValueError(f"Unknown input type or missing sra_accession: {input_type}")

# Ensure the log file exists (create it if it doesn't)
if not os.path.exists(log_file):
    open(log_file, 'a').close()

print(f"Input preparation completed. Output saved to: {output_file}")