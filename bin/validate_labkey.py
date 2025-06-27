#!/usr/bin/env python3

import argparse
from labkey.api_wrapper import APIWrapper

# Define the expected schema for the BLAST hit lists
blast_fields = [
    'Experiment', 'Blast Task', 'Sample Id', 'Qseqid', 'Qlen', 'Sseqid',
    'Stitle', 'Tax Rank', 'Length', 'Pident', 'Evalue', 'Bitscore',
    'Blast Db Version', 'Snakemake Run Id', 'Sscinames', 'Mapped Reads',
    'Stat Db Version', 'Total Reads', 'Exclude'
]

# Define the expected schema for the FASTA output lists
blast_fasta_fields = [
    'Experiment', 'Sample Id', 'Contig Id', 'Contig Sequence',
    'Notes', 'Snakemake Run Id'
]

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--server')              # LabKey server URL
    parser.add_argument('--container')           # LabKey container path (e.g., project folder)
    parser.add_argument('--list')                # Name of the LabKey list to validate
    parser.add_argument('--api_key')             # API key with insert/delete permissions
    parser.add_argument('--type', choices=['blast', 'fasta'], required=True)  # List type selector
    args = parser.parse_args()

    # Initialize LabKey API wrapper
    labkey = APIWrapper(args.server, args.container, api_key=args.api_key)

    # Use a unique experiment ID that's unlikely to conflict with real data
    unique_experiment_id = 123456789

    # Construct a dummy row with properly typed fields matching the expected schema
    if args.type == 'blast':
        dummy = {
            'Experiment': unique_experiment_id,
            'Blast Task': '__DUMMY__TASK__',
            'Sample Id': '__DUMMY__SAMPLE__',
            'Qseqid': '__DUMMY__QSEQID__',
            'Qlen': 0,
            'Sseqid': '__DUMMY__SSEQID__',
            'Stitle': '__DUMMY__STITLE__',
            'Tax Rank': '__DUMMY__TAX__',
            'Length': 0,
            'Pident': 0.0,
            'Evalue': 0.0,
            'Bitscore': 0,
            'Blast Db Version': 'unknown',
            'Snakemake Run Id': '__DUMMY__SMK__',
            'Sscinames': 0.0,
            'Mapped Reads': 0,
            'Stat Db Version': 'unknown',
            'Total Reads': 0.0,
            'Exclude': False
        }
    else:
        dummy = {
            'Experiment': unique_experiment_id,
            'Sample Id': '__DUMMY__SAMPLE__',
            'Contig Id': '__DUMMY__CONTIG__',
            'Contig Sequence': 'ATGCATGC',
            'Notes': 0.0,
            'Snakemake Run Id': '__DUMMY__SMK__'
        }

    try:
        print(f"[INFO] Trying insert to list '{args.list}'...")

        # Attempt to insert the dummy row into the specified LabKey list
        insert_result = labkey.query.insert_rows("lists", args.list, [dummy])
        inserted_row = insert_result['rows'][0]

        # Get the Key field if it exists for logging
        row_key = inserted_row.get('Key', 'N/A')
        print(f"[OK] Inserted dummy row with Key={row_key}, Experiment={unique_experiment_id}")

        # Delete the dummy row using the complete inserted row data
        # LabKey requires the full row data (including auto-generated fields like Key)
        delete_result = labkey.query.delete_rows("lists", args.list, [inserted_row])
        deleted_key = delete_result['rows'][0].get('Key', 'N/A')
        print(f"[OK] Deleted dummy row successfully with Key={deleted_key}")

    except Exception as e:
        print(f"[FAIL] Insert/Delete failed: {e}")
        exit(1)

if __name__ == '__main__':
    main()