#!/usr/bin/env python3

import argparse
import csv
from labkey.api_wrapper import APIWrapper

def upload_fasta_to_labkey(server, container, list_name, api_key, input_tsv):
    print(f"[INFO] Connecting to LabKey at {server} / {container}")
    labkey = APIWrapper(server, container, api_key=api_key)

    try:
        # Read TSV input file
        print(f"[INFO] Reading TSV file: {input_tsv}")
        with open(input_tsv, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            rows = list(reader)

        if not rows:
            print("[WARN] No rows found in the input file. Aborting.")
            return

        print(f"[INFO] Parsed {len(rows)} rows for upload.")
        response = labkey.query.insert_rows("lists", list_name, rows)

        print("[✅] Upload complete.")
        for row in response.get("rows", []):
            key = row.get("Key", "N/A")
            sample = row.get("Sample Id", "N/A")
            print(f"   - Inserted: Sample Id={sample}, Key={key}")

    except Exception as e:
        print(f"[❌ ERROR] Upload failed: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Upload GOTTCHA2 FASTA TSV to LabKey list.")
    parser.add_argument('--server', required=True, help="LabKey server URL")
    parser.add_argument('--container', required=True, help="LabKey container/project path")
    parser.add_argument('--list', required=True, help="LabKey list name to upload to")
    parser.add_argument('--api_key', required=True, help="LabKey user API key with insert permissions")
    parser.add_argument('--input_tsv', required=True, help="Path to TSV file to upload")

    args = parser.parse_args()

    upload_fasta_to_labkey(
        server=args.server,
        container=args.container,
        list_name=args.list,
        api_key=args.api_key,
        input_tsv=args.input_tsv
    )
