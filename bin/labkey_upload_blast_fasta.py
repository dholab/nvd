#!/usr/bin/env python3

import os
import csv
import sys
import argparse
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(description="Upload FASTA CSVs to LabKey.")
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--labkey-server", required=True)
    parser.add_argument("--labkey-project-name", required=True)
    parser.add_argument("--labkey-api-key", required=True)
    parser.add_argument("--labkey-schema", required=True)
    parser.add_argument("--table-name", default="fasta_hits_test_nvd2")
    args = parser.parse_args()

    log_entries = [
        f"LabKey FASTA Upload Log - {datetime.now()}",
        f"Experiment ID: {args.experiment_id}",
        f"Run ID: {args.run_id}",
        f"Server: {args.labkey_server}",
        f"Project: {args.labkey_project_name}",
        f"Target Table: {args.table_name}",
        "=" * 80
    ]

    upload_enabled = bool(args.labkey_server and args.labkey_project_name and args.labkey_api_key)
    if upload_enabled:
        try:
            from labkey.api_wrapper import APIWrapper
            import more_itertools
            log_entries.append("LabKey API wrapper found - attempting real upload")
            api = APIWrapper(args.labkey_server, args.labkey_project_name, api_key=args.labkey_api_key)
        except ImportError as e:
            log_entries.append(f"ERROR: LabKey API wrapper not available - {str(e)}")
            log_entries.append("Falling back to simulation mode")
            upload_enabled = False
    else:
        log_entries.append("No LabKey credentials provided - running in simulation mode")

    csv_files = [f for f in os.listdir('.') if f.endswith('.csv') and 'fasta' in f.lower()]
    log_entries.append(f"Found {len(csv_files)} FASTA CSV files: {csv_files}")

    total_records_processed = 0
    total_records_uploaded = 0

    for csv_file in sorted(csv_files):
        log_entries.append(f"\nProcessing FASTA file: {csv_file}")
        record_count = 0

        if os.path.getsize(csv_file) > 0:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                records = list(reader)
                record_count = len(records)
                total_records_processed += record_count

                if record_count > 0:
                    log_entries.append(f"  Records: {record_count}")
                    log_entries.append(f"  Sample fields: {', '.join(list(records[0].keys())[:5])}")

                    if upload_enabled:
                        batch_size = 1000
                        record_chunks = list(more_itertools.chunked(records, batch_size))
                        success_count = 0
                        for i, chunk in enumerate(record_chunks, 1):
                            try:
                                api.query.insert_rows(schema_name=args.labkey_schema, query_name=args.table_name, rows=chunk)
                                log_entries.append(f"    Batch {i}/{len(record_chunks)}: SUCCESS ({len(chunk)} records)")
                                success_count += len(chunk)
                            except Exception as e:
                                log_entries.append(f"    Batch {i}/{len(record_chunks)}: ERROR - {str(e)}")
                        total_records_uploaded += success_count
                    else:
                        num_batches = (record_count + 999) // 1000
                        log_entries.append(f"  Would upload in {num_batches} batch(es) (SIMULATION)")
                        total_records_uploaded += record_count
                else:
                    log_entries.append("  No records found in file")
        else:
            log_entries.append("  Empty file - no records to upload")

    log_entries += [
        "\n" + "=" * 80,
        "FASTA DATA UPLOAD SUMMARY",
        f"Files processed: {len(csv_files)}",
        f"Total records processed: {total_records_processed}",
        f"Total records uploaded: {total_records_uploaded}" if upload_enabled else f"Total records that would be uploaded: {total_records_uploaded}",
        "FASTA UPLOAD COMPLETE" if upload_enabled else "FASTA SIMULATION COMPLETE - No actual upload performed"
    ]

    with open('fasta_labkey_upload.log', 'w') as f:
        f.write('\n'.join(log_entries))

if __name__ == "__main__":
    main()
