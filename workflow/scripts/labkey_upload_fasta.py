import os
import more_itertools
from labkey.api_wrapper import APIWrapper

def reformat_fasta(input_fasta, output_fasta):
    """
    Use reformat.sh to remove newlines in the FASTA sequences.
    """
    command = f"reformat.sh in={input_fasta} out={output_fasta} fastawrap=0 ow=t"
    os.system(command)
    return output_fasta

def parse_fasta(fasta_file):
    """
    Parse FASTA file without using Biopython.
    """
    with open(fasta_file, 'r') as file:
        header = ''
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    yield header, ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, ''.join(sequence)

def insert_fasta_records(experiment, fasta_file, sample_name, api_key, snakemake_run_id):
    """
    Insert records into LabKey if API key is provided.
    """
    if not api_key:
        print(f"Skipping LabKey upload because API key is missing.", file=sys.stderr)
        return  # Exit function if API key is missing

    fasta_records = []
    for header, sequence in parse_fasta(fasta_file):
        fasta_records.append({
            'experiment': str(experiment),
            'sample_id': str(sample_name),
            'contig_id': str(header),
            'contig_sequence': str(sequence),
            'notes': '',
            'snakemake_run_id': str(snakemake_run_id)
        })

    # LabKey has a hard time handling large numbers of inserts at a time
    # Break apart fasta_records into smaller chunks of 1000 rows each
    fasta_record_chunks = list(more_itertools.chunked(fasta_records, 1000))

    api = APIWrapper('dholk.primate.wisc.edu', 'dho/projects/FuturePath/InfinitePath', api_key=api_key)

    for chunk in fasta_record_chunks:
        api.query.insert_rows(
            schema_name='lists',
            query_name='fasta_hits',
            rows=chunk
        )

# Main execution
if __name__ == "__main__":
    # Reformat the FASTA file to remove newlines
    formatted_fasta = reformat_fasta(snakemake.input.fasta, snakemake.params.singleline_fasta_path)

    # Run LabKey insert function if API key is present
    insert_fasta_records(
        experiment=snakemake.params.experiment,
        fasta_file=formatted_fasta,
        sample_name=snakemake.wildcards.sample,
        api_key=snakemake.params.api_key,
        snakemake_run_id=snakemake.params.snakemake_run_id
    )

    # Ensure the token file is created even if the LabKey upload is skipped
    with open(snakemake.output.token, 'w') as token_file:
        token_file.write('Completed')
    print(f"Token file created at: {snakemake.output.token}", file=sys.stderr)
