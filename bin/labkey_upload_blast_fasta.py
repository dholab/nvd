#!/usr/bin/env python3

import argparse
import csv
import os
from datetime import datetime

import more_itertools


def main():
    parser = argparse.ArgumentParser(description="Upload FASTA CSVs to LabKey.")
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument(
        "--sample-set-id", required=False, help="Sample set ID for upload tracking"
    )
    parser.add_argument(
        "--state-dir",
        required=False,
        help="State directory for upload tracking database",
    )
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
        "=" * 80,
    ]

    upload_enabled = bool(
        args.labkey_server and args.labkey_project_name and args.labkey_api_key
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
            "No LabKey credentials provided - running in simulation mode"
        )

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

                # Filter out samples that were already uploaded (cross-run deduplication)
                # Justification for local import: script must work standalone without py_nvd
                if args.sample_set_id and records and "sample_id" in records[0]:
                    try:
                        import py_nvd.state as nvd_state

                        unique_samples = set(r.get("sample_id") for r in records)
                        samples_to_skip = set()
                        for sample_id in unique_samples:
                            if sample_id and nvd_state.was_sample_ever_uploaded(
                                str(sample_id),
                                upload_type="blast_fasta",
                                upload_target="labkey",
                                state_dir=args.state_dir,
                            ):
                                samples_to_skip.add(sample_id)
                                log_entries.append(
                                    f"  SKIP: Sample '{sample_id}' was already uploaded in a previous run",
                                )
                        if samples_to_skip:
                            records = [
                                r
                                for r in records
                                if r.get("sample_id") not in samples_to_skip
                            ]
                            log_entries.append(
                                f"  Filtered out {len(samples_to_skip)} already-uploaded samples",
                            )
                    except ImportError:
                        log_entries.append(
                            "  WARNING: py_nvd not available, skipping cross-run check",
                        )

                record_count = len(records)
                total_records_processed += record_count

                if record_count > 0:
                    log_entries.append(f"  Records: {record_count}")
                    log_entries.append(
                        f"  Sample fields: {', '.join(list(records[0].keys())[:5])}"
                    )

                    if upload_enabled:
                        batch_size = 1000
                        record_chunks = list(
                            more_itertools.chunked(records, batch_size)
                        )
                        success_count = 0
                        for i, chunk in enumerate(record_chunks, 1):
                            try:
                                api.query.insert_rows(
                                    schema_name=args.labkey_schema,
                                    query_name=args.table_name,
                                    rows=chunk,
                                )
                                log_entries.append(
                                    f"    Batch {i}/{len(record_chunks)}: SUCCESS ({len(chunk)} records)",
                                )
                                success_count += len(chunk)
                            except Exception as e:
                                log_entries.append(
                                    f"    Batch {i}/{len(record_chunks)}: ERROR - {e!s}",
                                )
                        total_records_uploaded += success_count

                        # Record successful uploads in state database (per-sample)
                        # Justification for local import: script must work standalone without py_nvd
                        if (
                            args.sample_set_id
                            and success_count > 0
                            and records
                            and "sample_id" in records[0]
                        ):
                            try:
                                import py_nvd.state as nvd_state

                                uploaded_samples = set(
                                    r.get("sample_id")
                                    for r in records
                                    if r.get("sample_id")
                                )
                                for sample_id in uploaded_samples:
                                    sample_records = [
                                        r
                                        for r in records
                                        if r.get("sample_id") == sample_id
                                    ]
                                    content_hash = nvd_state.hash_upload_content(
                                        sample_records
                                    )
                                    nvd_state.record_upload(
                                        sample_id=str(sample_id),
                                        sample_set_id=args.sample_set_id,
                                        upload_type="blast_fasta",
                                        upload_target="labkey",
                                        content_hash=content_hash,
                                        target_metadata={
                                            "experiment_id": args.experiment_id,
                                            "table_name": args.table_name,
                                            "records_uploaded": len(sample_records),
                                            "source_file": csv_file,
                                        },
                                        state_dir=args.state_dir,
                                    )
                                    # Mark sample as uploaded (terminal state for LabKey runs)
                                    nvd_state.mark_sample_uploaded(
                                        sample_id=str(sample_id),
                                        sample_set_id=args.sample_set_id,
                                        state_dir=args.state_dir,
                                    )
                                    # Release sample lock after successful upload
                                    if args.run_id:
                                        nvd_state.release_sample_lock(
                                            sample_id=str(sample_id),
                                            run_id=args.run_id,
                                            state_dir=args.state_dir,
                                        )
                                log_entries.append(
                                    f"  Recorded uploads for {len(uploaded_samples)} samples in state database",
                                )
                            except ImportError:
                                log_entries.append(
                                    "  WARNING: py_nvd.state not available, uploads not recorded",
                                )
                            except Exception as e:
                                log_entries.append(
                                    f"  WARNING: Failed to record uploads: {e!s}"
                                )
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

    with open("fasta_labkey_upload.log", "w") as f:
        f.write("\n".join(log_entries))


if __name__ == "__main__":
    main()
