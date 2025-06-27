#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Upload GOTTCHA2 full TSV data to LabKey")
    parser.add_argument('--input-tsv', required=True, help='Path to GOTTCHA2 full TSV file')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--experiment', required=True, help='Experiment name')
    parser.add_argument('--labkey-server', required=True, help='LabKey server address (e.g., labkey.org)')
    parser.add_argument('--labkey-schema', required=True)
    parser.add_argument('--labkey-project-name', required=True, help='LabKey project name')
    parser.add_argument('--labkey-api-key', required=True, help='LabKey API key')
    parser.add_argument('--db-version', required=True, help='GOTTCHA2 DB version')
    parser.add_argument('--table-name', required=True, help="LabKey table to insert data into")
    return parser.parse_args()

def main():
    args = parse_args()

    log_entries = [
        f"GOTTCHA2 Upload Log - {datetime.now()}",
        f"Sample: {args.sample}",
        f"Experiment: {args.experiment}",
        f"Server: {args.labkey_server}",
        f"Project: {args.labkey_project_name}",
        f"Schema: {args.labkey_schema}",
        f"Target Table: {args.table_name}",
        "=" * 80
    ]

    # Check if upload should proceed
    upload_enabled = bool(args.labkey_server and args.labkey_project_name and args.labkey_api_key)
    api = None

    if upload_enabled:
        try:
            from labkey.api_wrapper import APIWrapper
            import more_itertools
            api = APIWrapper(
                args.labkey_server,
                args.labkey_project_name,
                api_key=args.labkey_api_key
)
            log_entries.append("LabKey API wrapper found - proceeding with upload.")
        except ImportError as e:
            log_entries.append(f"ERROR: LabKey API wrapper or dependencies not available - {str(e)}")
            log_entries.append("Falling back to simulation mode.")
            upload_enabled = False
    else:
        log_entries.append("Missing LabKey credentials - running in simulation mode.")

    # Read the input TSV
    if not os.path.exists(args.input_tsv):
        log_entries.append(f"ERROR: Input file not found - {args.input_tsv}")
        print('\n'.join(log_entries))
        sys.exit(1)

    df = pd.read_csv(args.input_tsv, sep='\t')

    if df.empty:
        log_entries.append("No data found in TSV file. Exiting.")
        print('\n'.join(log_entries))
        sys.exit(0)

    # Add metadata
    df['SAMPLE'] = args.sample
    df['EXPERIMENT'] = args.experiment
    df['GOTTCHA2_DB_VERSION'] = args.db_version

    # Perform upload or simulation
    total_uploaded = 0
    total_records = len(df)
    log_entries.append(f"Total records to process: {total_records}")

    if upload_enabled:
        batch_size = 1000
        record_chunks = list(more_itertools.chunked(df.to_dict("records"), batch_size))
        for i, chunk in enumerate(record_chunks, 1):
            try:
                api.query.insert_rows(
                    schema_name=args.labkey_schema,
                    query_name=args.table_name,
                    rows=chunk
                )
                log_entries.append(f"Batch {i}/{len(record_chunks)}: SUCCESS ({len(chunk)} records)")
                total_uploaded += len(chunk)
            except Exception as e:
                log_entries.append(f"Batch {i}/{len(record_chunks)}: ERROR - {str(e)}")
        log_entries.append(f"Total successfully uploaded: {total_uploaded} / {total_records}")
    else:
        batch_count = (total_records + 999) // 1000
        log_entries.append(f"SIMULATION MODE: Would upload {total_records} records in {batch_count} batches.")

    log_entries.append("=" * 80)
    log_entries.append("GOTTCHA2 UPLOAD COMPLETE" if upload_enabled else "GOTTCHA2 SIMULATION COMPLETE")

    # Save log
    with open("fasta_labkey_upload.log", "w") as f:
        f.write('\n'.join(log_entries))

    # Also print to stdout
    print('\n'.join(log_entries))

if __name__ == '__main__':
    main()
