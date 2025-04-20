#### filepath: scripts/labkey_upload_gottcha.py
#!/usr/bin/env python3
"""
Script to upload GOTTCHA2 results and validation summary to LabKey.
"""

import sys
import pandas as pd
import more_itertools
from labkey.api_wrapper import APIWrapper
from labkey.exceptions import ServerContextError, RequestError
import os
import argparse

# --- Constants ---
# Define target LabKey list names
SPECIES_LIST_NAME = "GOTTCHAResults"
VALIDATION_SUMMARY_LIST_NAME = "GOTTCHAValidationSummary"
CHUNK_SIZE = 1000

# --- Helper Functions ---

def print_error(*args, **kwargs):
    """Prints messages to stderr."""
    print(*args, file=sys.stderr, **kwargs)

def upload_dataframe(api, schema_name, query_name, df, description="data"):
    """Uploads a pandas DataFrame to LabKey in chunks."""
    if df.empty:
        print_error(f"Skipping upload: DataFrame for {description} is empty.")
        return 0 # Return 0 rows uploaded

    rows_to_upload = df.to_dict('records')
    row_chunks = list(more_itertools.chunked(rows_to_upload, CHUNK_SIZE))
    total_uploaded = 0

    print_error(f"Attempting to upload {len(rows_to_upload)} rows for {description} to {schema_name}.{query_name}...")

    for i, chunk in enumerate(row_chunks):
        try:
            result = api.query.insert_rows(schema_name=schema_name, query_name=query_name, rows=chunk)
            uploaded_in_chunk = result.get('rowsAffected', len(chunk)) # Use rowsAffected if available
            total_uploaded += uploaded_in_chunk
            print_error(f"  Chunk {i+1}/{len(row_chunks)}: Uploaded {uploaded_in_chunk} rows for {description}.")
        except (ServerContextError, RequestError) as e:
            print_error(f"Error inserting chunk {i+1} for {description}: {e}")
            print_error(f"  Server response: {e.response.text if hasattr(e, 'response') else 'No response text'}")
            print_error(f"  Problematic chunk data (first row): {chunk[0] if chunk else 'N/A'}")
            # Decide whether to continue or stop on error
            # print_error("  Continuing with next chunk...")
            print_error("  Stopping upload due to error.")
            return -1 # Indicate failure
        except Exception as e:
            print_error(f"An unexpected error occurred during chunk {i+1} upload for {description}: {e}")
            print_error(f"  Problematic chunk data (first row): {chunk[0] if chunk else 'N/A'}")
            print_error("  Stopping upload due to error.")
            return -1 # Indicate failure

    print_error(f"Finished uploading {total_uploaded} rows for {description}.")
    return total_uploaded


# --- Main Upload Functions ---

def upload_species_results(api, species_tsv, sample_id, experiment, db_version):
    """Reads GOTTCHA2 species TSV, adds metadata, and uploads."""
    print_error(f"\n--- Processing Species Results ({SPECIES_LIST_NAME}) ---")
    if not os.path.exists(species_tsv) or os.path.getsize(species_tsv) == 0:
        print_error(f"Warning: Species TSV file '{species_tsv}' is missing or empty. Skipping upload.")
        return 0

    try:
        df = pd.read_csv(species_tsv, sep='\t')
        if df.empty:
             print_error(f"Warning: Species TSV file '{species_tsv}' contained no data. Skipping upload.")
             return 0

        # Add standard columns
        df['sample_id'] = sample_id
        df['experiment'] = experiment
        df['gottcha2_db_version'] = db_version

        # Rename columns to match expected LabKey list fields (adjust as needed)
        # Example: df.rename(columns={'TAXA': 'TaxonName', 'ROLLUP_COUNT': 'ReadCount'}, inplace=True)
        # Make sure column names in the DataFrame match the LabKey list exactly (case-sensitive)

        # Perform the upload
        return upload_dataframe(api, 'lists', SPECIES_LIST_NAME, df, description="GOTTCHA2 species results")

    except pd.errors.EmptyDataError:
        print_error(f"Warning: Species TSV file '{species_tsv}' was empty or invalid. Skipping upload.")
        return 0
    except Exception as e:
        print_error(f"Error processing species TSV file '{species_tsv}': {e}")
        return -1 # Indicate failure

def upload_validation_summary(api, summary_tsv, sample_id, experiment):
    """Reads validation summary TSV and uploads."""
    print_error(f"\n--- Processing Validation Summary ({VALIDATION_SUMMARY_LIST_NAME}) ---")
    if not os.path.exists(summary_tsv) or os.path.getsize(summary_tsv) == 0:
        print_error(f"Warning: Validation summary file '{summary_tsv}' is missing or empty. Skipping upload.")
        return 0

    try:
        df = pd.read_csv(summary_tsv, sep='\t')
        if df.empty:
             print_error(f"Warning: Validation summary file '{summary_tsv}' contained no data. Skipping upload.")
             return 0

        # Add standard columns if they weren't added during aggregation
        if 'sample_id' not in df.columns:
             df['sample_id'] = sample_id
        if 'experiment' not in df.columns:
             df['experiment'] = experiment

        # Rename columns if necessary to match LabKey list fields
        # Example: df.rename(columns={'validation_taxid': 'TargetTaxID'}, inplace=True)

        # Perform the upload
        return upload_dataframe(api, 'lists', VALIDATION_SUMMARY_LIST_NAME, df, description="GOTTCHA2 validation summary")

    except pd.errors.EmptyDataError:
        print_error(f"Warning: Validation summary file '{summary_tsv}' was empty or invalid. Skipping upload.")
        return 0
    except Exception as e:
        print_error(f"Error processing validation summary file '{summary_tsv}': {e}")
        return -1 # Indicate failure


# --- Main Execution ---

def main(args):
    """Connects to LabKey and orchestrates uploads."""
    print_error("Starting GOTTCHA2 LabKey upload script...")
    print_error(f"  Sample: {args.sample_id}")
    print_error(f"  Experiment: {args.experiment}")
    print_error(f"  Species TSV: {args.species_tsv}")
    print_error(f"  Validation Summary: {args.validation_summary}")
    print_error(f"  DB Version: {args.db_version}")

    success = True

    if not args.api_key:
        print_error("Warning: No LabKey API key provided. Skipping LabKey uploads.")
        # Still create token file to allow workflow to proceed? Or fail?
        # Let's assume we proceed but indicate skipped upload.
        success = False # Mark as not fully successful regarding LabKey
    else:
        try:
            api = APIWrapper(args.labkey_server, args.project_name, api_key=args.api_key, use_ssl=True)
            print_error("LabKey connection configured.")

            # Upload Species Results
            species_upload_status = upload_species_results(api, args.species_tsv, args.sample_id, args.experiment, args.db_version)
            if species_upload_status == -1:
                success = False # Upload failed

            # Upload Validation Summary
            validation_upload_status = upload_validation_summary(api, args.validation_summary, args.sample_id, args.experiment)
            if validation_upload_status == -1:
                success = False # Upload failed

        except Exception as e:
            print_error(f"Failed to initialize LabKey API connection or during upload process: {e}")
            success = False

    # Create the token file to signal completion, regardless of LabKey success?
    # Or only if success is True? Let's create it regardless for now,
    # but log the status.
    try:
        with open(args.token_file, 'w') as f:
            f.write('Completed')
        print_error(f"\nToken file created: {args.token_file}")
        if not success:
             print_error("Note: Token file created, but LabKey upload may have been skipped or encountered errors.")
    except Exception as e:
        print_error(f"Error creating token file '{args.token_file}': {e}")
        sys.exit(1) # Exit with error if token cannot be created

    if not success:
        print_error("Script finished with warnings or errors during LabKey upload.")
        # Optionally exit with a non-zero status code if LabKey upload is critical
        # sys.exit(1)
    else:
        print_error("Script finished successfully.")


if __name__ == "__main__":
    # Check if running under Snakemake
    if 'snakemake' in globals():
        # Use Snakemake inputs/outputs/params directly
        parser = argparse.ArgumentParser(description="Dummy parser for Snakemake execution.")
        # Add arguments just to store Snakemake object attributes
        parser.add_argument('--species_tsv', default=snakemake.input.tsv_species)
        parser.add_argument('--validation_summary', default=snakemake.input.validation_summary)
        parser.add_argument('--token_file', default=snakemake.output.token)
        parser.add_argument('--sample_id', default=snakemake.wildcards.sample)
        parser.add_argument('--experiment', default=snakemake.params.experiment)
        parser.add_argument('--labkey_server', default=snakemake.params.labkey_server)
        parser.add_argument('--project_name', default=snakemake.params.project_name)
        parser.add_argument('--api_key', default=snakemake.params.api_key)
        parser.add_argument('--db_version', default=snakemake.params.gottcha2_db_version)
        args = parser.parse_args([]) # Pass empty list as Snakemake provides globals
    else:
        # Use argparse for standalone execution (optional, for testing)
        parser = argparse.ArgumentParser(description="Upload GOTTCHA2 results to LabKey.")
        parser.add_argument('--species_tsv', required=True, help="Path to GOTTCHA2 species TSV file.")
        parser.add_argument('--validation_summary', required=True, help="Path to validation summary TSV file.")
        parser.add_argument('--token_file', required=True, help="Path to output token file.")
        parser.add_argument('--sample_id', required=True, help="Sample identifier.")
        parser.add_argument('--experiment', required=True, help="Experiment identifier.")
        parser.add_argument('--labkey_server', required=True, help="LabKey server URL (e.g., https://your.labkey.org).")
        parser.add_argument('--project_name', required=True, help="LabKey project path.")
        parser.add_argument('--api_key', help="LabKey API key (omit to skip upload).")
        parser.add_argument('--db_version', required=True, help="GOTTCHA2 DB version string.")
        args = parser.parse_args()

    main(args)