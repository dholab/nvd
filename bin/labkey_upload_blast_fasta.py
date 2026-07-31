#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "labkey>=4.0.1",
# ]
# ///

import argparse
import csv
import os
import sys
from datetime import datetime
from typing import Any


def sample_already_uploaded(query_api, schema, table, experiment, sample_id) -> bool:
    """True if the FASTA list already holds this sample's contigs (its ledger).

    FASTA is contig-only: one unit per sample, so presence is checked by
    (experiment, sample_id) alone (no query_class). Reuses the guard script's
    filter fallbacks (QueryFilter objects, then filter-string syntax) so this
    still works against older labkey-api-python versions.
    """
    try:
        from labkey.query import QueryFilter

        result = query_api.select_rows(
            schema_name=schema,
            query_name=table,
            filter_array=[
                QueryFilter("experiment", experiment, "eq"),
                QueryFilter("sample_id", sample_id, "eq"),
            ],
            max_rows=1,
        )
    except (ImportError, AttributeError):
        result = query_api.select_rows(
            schema_name=schema,
            query_name=table,
            filter_array=[f"experiment~eq={experiment}", f"sample_id~eq={sample_id}"],
        )
    return bool(result and result.get("rows"))


def insert_records(
    query_api,
    schema: str,
    table: str,
    records: list[dict[str, Any]],
) -> None:
    """Insert every record for a sample's contigs in one atomic call.

    Deliberately does not catch anything: a failed insert must propagate so the
    caller can hard-fail the run rather than log an error and report success.
    """
    query_api.insert_rows(
        schema_name=schema,
        query_name=table,
        rows=records,
    )


def _write_log(log_entries: list[str]) -> None:
    """Write the accumulated log entries to the upload log file and stdout."""
    text = "\n".join(log_entries)
    with open("fasta_labkey_upload.log", "w") as f:
        f.write(text)
    print(text)


def main():
    parser = argparse.ArgumentParser(description="Upload FASTA CSVs to LabKey.")
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--labkey-server", required=True)
    parser.add_argument("--labkey-project-name", required=True)
    parser.add_argument("--labkey-api-key", required=True)
    parser.add_argument("--labkey-schema", required=True)
    parser.add_argument("--table-name", default="fasta_hits_test_nvd2")
    args = parser.parse_args()

    log_entries = [
        f"LabKey FASTA Upload Log - {datetime.now()}",
        f"Experiment ID: {args.experiment_id}",
        f"Sample: {args.sample_id}",
        f"Server: {args.labkey_server}",
        f"Project: {args.labkey_project_name}",
        f"Target Table: {args.table_name}",
        "=" * 80,
    ]

    upload_enabled = bool(
        args.labkey_server and args.labkey_project_name and args.labkey_api_key,
    )
    if upload_enabled:
        try:
            from labkey.api_wrapper import APIWrapper

            log_entries.append("LabKey API wrapper found - attempting real upload")
            api = APIWrapper(
                args.labkey_server,
                args.labkey_project_name,
                api_key=args.labkey_api_key,
            )
        except ImportError as e:
            log_entries.append(f"ERROR: LabKey API wrapper not available - {e!s}")
            log_entries.append("Falling back to simulation mode")
            upload_enabled = False
    else:
        log_entries.append(
            "No LabKey credentials provided - running in simulation mode",
        )

    # The destination list is its own completion ledger. If this sample already
    # has rows there, a prior run already uploaded its contigs: skip
    # re-inserting rather than risk duplicating the FASTA list.
    if upload_enabled and sample_already_uploaded(
        api.query,
        args.labkey_schema,
        args.table_name,
        int(args.experiment_id),
        args.sample_id,
    ):
        log_entries.append(
            f"SKIP: sample already uploaded (exp={args.experiment_id}, "
            f"sample={args.sample_id}); no insert.",
        )
        _write_log(log_entries)
        return

    csv_files = [
        f for f in os.listdir(".") if f.endswith(".csv") and "fasta" in f.lower()
    ]
    log_entries.append(f"Found {len(csv_files)} FASTA CSV files: {csv_files}")

    total_records_processed = 0
    total_records_uploaded = 0

    for csv_file in sorted(csv_files):
        log_entries.append(f"\nProcessing FASTA file: {csv_file}")
        record_count = 0

        if os.path.getsize(csv_file) > 0:
            with open(csv_file) as f:
                reader = csv.DictReader(f)
                records = list(reader)

                record_count = len(records)
                total_records_processed += record_count

                if record_count > 0:
                    log_entries.append(f"  Records: {record_count}")
                    log_entries.append(
                        f"  Sample fields: {', '.join(list(records[0].keys())[:5])}",
                    )

                    if upload_enabled:
                        # Insert the whole file atomically: chunking only adds
                        # partial-failure risk (some rows committed, some not)
                        # without a throughput benefit at this batch size.
                        try:
                            insert_records(
                                api.query, args.labkey_schema, args.table_name, records,
                            )
                        except Exception as e:
                            log_entries.append(f"  Upload: ERROR - {e!s}")
                            log_entries.append(
                                "\nFASTA UPLOAD FAILED - insert error, see above (no rows committed)",
                            )
                            _write_log(log_entries)
                            print(
                                f"ERROR: LabKey insert failed for "
                                f"sample={args.sample_id}: {e!s}",
                                file=sys.stderr,
                            )
                            sys.exit(1)
                        log_entries.append(f"  Upload: SUCCESS ({len(records)} records)")
                        total_records_uploaded += record_count

                    else:
                        num_batches = (record_count + 999) // 1000
                        log_entries.append(
                            f"  Would upload in {num_batches} batch(es) (SIMULATION)",
                        )
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
        f"Total records uploaded: {total_records_uploaded}"
        if upload_enabled
        else f"Total records that would be uploaded: {total_records_uploaded}",
        "FASTA UPLOAD COMPLETE"
        if upload_enabled
        else "FASTA SIMULATION COMPLETE - No actual upload performed",
    ]

    _write_log(log_entries)


if __name__ == "__main__":
    main()
