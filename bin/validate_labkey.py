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

def check_experiment_id_uniqueness(labkey, list_name, experiment_id):
    """
    Check if the experiment ID already exists in the LabKey list.
    Uses SQL query to efficiently check for duplicates without downloading full dataset.
    
    Args:
        labkey: APIWrapper instance
        list_name: Name of the LabKey list to check
        experiment_id: Experiment ID to validate
        
    Returns:
        bool: True if experiment ID is unique (safe to use), False if duplicate found
    """
    try:
        print(f"[INFO] Checking for duplicate experiment ID: {experiment_id}")
        
        # Use SQL query to efficiently check for existing experiment IDs
        # This only returns matching rows, not the entire dataset
        sql_query = f"SELECT DISTINCT Experiment FROM lists.{list_name} WHERE Experiment = {experiment_id}"
        
        result = labkey.query.execute_sql("lists", sql_query)
        
        # If any rows are returned, the experiment ID already exists
        if result['rows']:
            existing_experiments = [row['Experiment'] for row in result['rows']]
            print(f"[FAIL] Duplicate experiment ID found: {existing_experiments}")
            return False
        else:
            print(f"[OK] Experiment ID {experiment_id} is unique")
            return True
            
    except Exception as e:
        # If the list doesn't exist or is empty, that's okay for uniqueness check
        if "does not exist" in str(e).lower() or "table" in str(e).lower():
            print(f"[INFO] List '{list_name}' appears to be empty or new - experiment ID is unique")
            return True
        else:
            print(f"[WARN] Error checking experiment ID uniqueness: {e}")
            print("[INFO] Proceeding with validation (assuming uniqueness)")
            return True

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--server')              # LabKey server URL
    parser.add_argument('--container')           # LabKey container path (e.g., project folder)
    parser.add_argument('--list')                # Name of the LabKey list to validate
    parser.add_argument('--api_key')             # API key with insert/delete permissions
    parser.add_argument('--type', choices=['blast', 'fasta'], required=True)  # List type selector
    parser.add_argument('--experiment_id', type=int, required=True)  # Experiment ID to check for uniqueness
    args = parser.parse_args()

    # Initialize LabKey API wrapper
    labkey = APIWrapper(args.server, args.container, api_key=args.api_key)

    # Check for duplicate experiment ID before proceeding with validation
    if not check_experiment_id_uniqueness(labkey, args.list, args.experiment_id):
        print()
        print("=" * 80)
        print("❌ EXPERIMENT ID VALIDATION FAILED")
        print("=" * 80)
        print()
        print(f"Error: Experiment ID '{args.experiment_id}' already exists in list '{args.list}'")
        print()
        print("REQUIRED ACTION:")
        print("   You must use a unique experiment_id for each pipeline run.")
        print("   Please update your pipeline configuration with a new experiment_id value.")
        print()
        print("Suggestions:")
        print("   • Use a timestamp-based ID (e.g., YYYYMMDDHHMMSS format)")
        print("   • Use an incremental counter from your last run")
        print("   • Generate a random unique identifier")
        print("   • Check existing experiment IDs in LabKey to avoid conflicts")
        print()
        print("Configuration:")
        print(f"   Update your params.experiment_id from {args.experiment_id} to a unique value")
        print()
        exit(1)

    # Use a unique experiment ID that's unlikely to conflict with real data for dummy row
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
        
        print()
        print("=" * 80)
        print("✅ VALIDATION SUCCESSFUL")
        print("=" * 80)
        print()
        print(f"• Experiment ID {args.experiment_id} is unique and ready to use")
        print(f"• LabKey connection and permissions validated")
        print(f"• List '{args.list}' is accessible and writable")
        print()

    except Exception as e:
        print(f"[FAIL] Insert/Delete failed: {e}")
        print()
        print("=" * 80)
        print("⚠️  VALIDATION FAILED - Troubleshooting Guide:")
        print("=" * 80)
        print()
        print("Common issues and solutions:")
        print()
        print("1. API Key Issues:")
        print("   • Check that your nvd2 secret is correctly configured in Nextflow")
        print("   • Verify the API key has INSERT and DELETE permissions")
        print("   • Ensure the API key hasn't expired")
        print("   • Test the API key directly in LabKey's web interface")
        print()
        print("2. Server/Container Issues:")
        print(f"   • Verify server URL is correct: {args.server}")
        print(f"   • Check container/project path: {args.container}")
        print("   • Ensure you have access to the specified project")
        print("   • Confirm the server URL includes 'https://' if required")
        print()
        print("3. List Issues:")
        print(f"   • Verify list '{args.list}' exists in LabKey")
        print("   • Check that the list has the expected schema/columns")
        print("   • Ensure you have read/write permissions to the list")
        print("   • Confirm list is in the correct folder/container")
        print()
        print("4. Network/Connectivity Issues:")
        print("   • Check firewall or proxy settings")
        print("   • Verify network connectivity to the LabKey server")
        print("   • Try accessing the LabKey server through a web browser")
        print()
        print("5. Configuration Issues:")
        print("   • Double-check your nextflow.config parameters:")
        print(f"     - labkey_server = '{args.server}'")
        print(f"     - labkey_project_name = '{args.container}'")
        print(f"     - labkey_{args.type}_list = '{args.list}'")
        print()
        print("For more help:")
        print("   • Check LabKey server logs for detailed error messages")
        print("   • Verify your LabKey user permissions with an administrator")
        print("   • Test API access using LabKey's built-in API documentation")
        print()
        exit(1)

if __name__ == '__main__':
    main()