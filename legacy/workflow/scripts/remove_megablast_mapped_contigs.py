#!/usr/bin/env python3

import subprocess
import logging

def setup_logger(log_file):
    logger = logging.getLogger('remove_megablast_mapped_contigs')
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

def extract_unique_ids(input_file, output_file, logger):
    unique_ids = set()
    logger.info(f"Reading megablast results from {input_file}")
    with open(input_file, 'r') as f:
        for line in f:
            unique_ids.add(line.rstrip('\n').split('\t')[0] + '\n')
    
    logger.info(f"Writing classified contigs to {output_file}")
    with open(output_file, 'w') as f:
        for id in unique_ids:
            f.write(id)
    logger.info(f"Wrote {len(unique_ids)} unique contig IDs")

def run_filterbyname(input_fasta, output_fasta, names_file, logger):
    command = [
        "filterbyname.sh",
        f"in={input_fasta}",
        f"out={output_fasta}",
        f"names={names_file}"
    ]
    logger.info(f"Running command: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        logger.info(result.stdout)
        if result.stderr:
            logger.warning(result.stderr)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running filterbyname.sh: {e}")
        logger.error(e.stderr)
        raise

def main(snakemake):
    logger = setup_logger(snakemake.log[0])

    logger.info("Starting remove_megablast_mapped_contigs script")

    try:
        extract_unique_ids(snakemake.input.megablast_results, 
                           snakemake.output.classified_contigs, 
                           logger)

        run_filterbyname(snakemake.input.contigs_fasta, 
                         snakemake.output.pruned_contigs, 
                         snakemake.output.classified_contigs, 
                         logger)

        logger.info("Script completed successfully")
    except Exception as e:
        logger.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main(snakemake)