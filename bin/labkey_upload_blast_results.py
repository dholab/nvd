#!/usr/bin/env python3

import os
import sys
import argparse
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

try:
    import polars as pl
except ImportError:
    print("ERROR: Polars is not installed. Please install it with: pip install polars")
    sys.exit(1)


def validate_dataframe(df: pl.DataFrame, csv_file: str, log_entries: List[str]) -> pl.DataFrame:
    """
    Validate and clean the DataFrame.
    
    Args:
        df: Polars DataFrame to validate
        csv_file: Name of the source file for logging
        log_entries: List to append log messages to
    
    Returns:
        Cleaned DataFrame
    """
    # Log initial shape
    log_entries.append(f"  Initial shape: {df.shape}")
    
    # Remove any completely null rows
    initial_rows = len(df)
    df = df.filter(~pl.all_horizontal(pl.all().is_null()))
    rows_removed = initial_rows - len(df)
    if rows_removed > 0:
        log_entries.append(f"  Removed {rows_removed} completely empty rows")
    
    # Log column info
    log_entries.append(f"  Columns: {', '.join(df.columns[:10])}{'...' if len(df.columns) > 10 else ''}")
    
    # Check for critical columns and log missing ones
    critical_columns = ['experiment', 'sample_id', 'qseqid']
    missing_columns = [col for col in critical_columns if col not in df.columns]
    if missing_columns:
        log_entries.append(f"  WARNING: Missing critical columns: {missing_columns}")
    
    # Convert numeric columns to proper types if they exist
    numeric_columns = {
        'length': pl.Int64,
        'pident': pl.Float64,
        'evalue': pl.Float64,
        'bitscore': pl.Float64,
        'staxids': pl.Int64,
        'adjusted_taxid': pl.Int64,
        'mapped_reads': pl.Int64,
        'total_reads': pl.Int64,
        'experiment': pl.Int64,
    }
    
    for col, dtype in numeric_columns.items():
        if col in df.columns:
            try:
                # Try to cast to the appropriate type
                df = df.with_columns(pl.col(col).cast(dtype, strict=False))
                log_entries.append(f"  Converted '{col}' to {dtype}")
            except Exception as e:
                log_entries.append(f"  WARNING: Could not convert '{col}' to {dtype}: {str(e)}")
    
    # Fill null values with appropriate defaults
    for col in df.columns:
        if df[col].dtype in [pl.Float32, pl.Float64]:
            df = df.with_columns(pl.col(col).fill_null(0.0))
        elif df[col].dtype in [pl.Int8, pl.Int16, pl.Int32, pl.Int64]:
            df = df.with_columns(pl.col(col).fill_null(0))
        elif df[col].dtype == pl.Utf8:
            df = df.with_columns(pl.col(col).fill_null(""))
    
    return df


def get_sample_stats(df: pl.DataFrame) -> Dict[str, Any]:
    """
    Get statistical summary of the DataFrame.
    
    Args:
        df: Polars DataFrame
    
    Returns:
        Dictionary with statistics
    """
    stats = {
        "total_rows": len(df),
        "unique_samples": df["sample_id"].n_unique() if "sample_id" in df.columns else 0,
        "unique_queries": df["qseqid"].n_unique() if "qseqid" in df.columns else 0,
    }
    
    # Get numeric column stats
    if "pident" in df.columns:
        stats["avg_pident"] = df["pident"].mean()
        stats["max_pident"] = df["pident"].max()
    
    if "bitscore" in df.columns:
        stats["avg_bitscore"] = df["bitscore"].mean()
        stats["max_bitscore"] = df["bitscore"].max()
    
    if "evalue" in df.columns:
        stats["min_evalue"] = df["evalue"].min()
    
    return stats


def dataframe_to_records(df: pl.DataFrame) -> List[Dict[str, Any]]:
    """
    Convert Polars DataFrame to list of dictionaries for LabKey upload.
    
    Args:
        df: Polars DataFrame
    
    Returns:
        List of dictionaries
    """
    # Convert to dictionaries, handling null values
    records = df.to_dicts()
    
    # Clean up records - convert None to empty string for LabKey
    cleaned_records = []
    for record in records:
        cleaned_record = {}
        for key, value in record.items():
            if value is None:
                cleaned_record[key] = ""
            elif isinstance(value, float) and value != value:  # Check for NaN
                cleaned_record[key] = ""
            else:
                cleaned_record[key] = value
        cleaned_records.append(cleaned_record)
    
    return cleaned_records


def process_csv_batch(df: pl.DataFrame, batch_size: int = 1000) -> List[pl.DataFrame]:
    """
    Split DataFrame into batches for upload.
    
    Args:
        df: Polars DataFrame
        batch_size: Number of records per batch
    
    Returns:
        List of DataFrame batches
    """
    n_batches = (len(df) + batch_size - 1) // batch_size
    batches = []
    
    for i in range(n_batches):
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, len(df))
        batch = df.slice(start_idx, end_idx - start_idx)
        batches.append(batch)
    
    return batches


def main():
    parser = argparse.ArgumentParser(description="Upload BLAST CSVs to LabKey using Polars.")
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--labkey-server", required=True)
    parser.add_argument("--labkey-project-name", required=True)
    parser.add_argument("--labkey-api-key", required=True)
    parser.add_argument("--labkey-schema", required=True)
    parser.add_argument("--table-name", default="metagenomic_hits_test_nvd2")
    parser.add_argument("--batch-size", type=int, default=1000, 
                        help="Number of records per upload batch")
    parser.add_argument("--validate-only", action="store_true",
                        help="Only validate data without uploading")
    args = parser.parse_args()

    log_entries = [
        f"LabKey BLAST Upload Log (Polars Version) - {datetime.now()}",
        f"Experiment ID: {args.experiment_id}",
        f"Run ID: {args.run_id}",
        f"Server: {args.labkey_server}",
        f"Project: {args.labkey_project_name}",
        f"Target Table: {args.table_name}",
        f"Batch Size: {args.batch_size}",
        "=" * 80,
    ]

    upload_enabled = bool(
        args.labkey_server and args.labkey_project_name and args.labkey_api_key
    ) and not args.validate_only
    
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
            log_entries.append(f"ERROR: LabKey API wrapper not available - {str(e)}")
            log_entries.append("Falling back to simulation mode")
            upload_enabled = False
    elif args.validate_only:
        log_entries.append("Running in validation-only mode")
    else:
        log_entries.append(
            "No LabKey credentials provided - running in simulation mode"
        )

    # Find all BLAST CSV files
    csv_files = [
        f for f in os.listdir(".") if f.endswith(".csv") and "blast" in f.lower()
    ]
    log_entries.append(f"Found {len(csv_files)} BLAST CSV files: {csv_files}")

    total_records_processed = 0
    total_records_uploaded = 0
    all_stats = {}

    for csv_file in sorted(csv_files):
        log_entries.append(f"\nProcessing BLAST file: {csv_file}")
        
        file_path = Path(csv_file)
        if file_path.stat().st_size == 0:
            log_entries.append("  Empty file - skipping")
            continue
        
        try:
            # Read CSV with Polars
            df = pl.read_csv(
                csv_file,
                separator="\t" if csv_file.endswith(".tsv") else ",",
                has_header=True,
                ignore_errors=True,  # Continue on parsing errors
                null_values=["", "NA", "N/A", "null", "NULL", "None"],
                infer_schema_length=10000,  # Infer types from more rows
            )
            
            # Validate and clean the DataFrame
            df = validate_dataframe(df, csv_file, log_entries)
            
            record_count = len(df)
            total_records_processed += record_count
            
            if record_count > 0:
                # Get statistics
                stats = get_sample_stats(df)
                all_stats[csv_file] = stats
                
                log_entries.append(f"  Records: {record_count}")
                log_entries.append(f"  Statistics:")
                for key, value in stats.items():
                    if isinstance(value, float):
                        log_entries.append(f"    {key}: {value:.4f}")
                    else:
                        log_entries.append(f"    {key}: {value}")
                
                # Show sample of first record
                if len(df) > 0:
                    first_record = df.head(1).to_dicts()[0]
                    sample_fields = list(first_record.keys())[:5]
                    log_entries.append(f"  Sample fields: {', '.join(sample_fields)}")
                    log_entries.append(f"  First record sample:")
                    for field in sample_fields:
                        value = first_record.get(field, "")
                        if isinstance(value, float):
                            log_entries.append(f"    {field}: {value:.4f}")
                        else:
                            log_entries.append(f"    {field}: {value}")
                
                if upload_enabled:
                    # Process in batches using Polars
                    batches = process_csv_batch(df, args.batch_size)
                    success_count = 0
                    
                    for i, batch_df in enumerate(batches, 1):
                        try:
                            # Convert batch to records
                            batch_records = dataframe_to_records(batch_df)
                            
                            # Upload to LabKey
                            api.query.insert_rows(
                                schema_name=args.labkey_schema,
                                query_name=args.table_name,
                                rows=batch_records,
                            )
                            log_entries.append(
                                f"    Batch {i}/{len(batches)}: SUCCESS ({len(batch_records)} records)"
                            )
                            success_count += len(batch_records)
                        except Exception as e:
                            log_entries.append(
                                f"    Batch {i}/{len(batches)}: ERROR - {str(e)}"
                            )
                            # Log first record of failed batch for debugging
                            if len(batch_df) > 0:
                                log_entries.append(f"      Failed batch first record: {batch_df.head(1).to_dicts()[0]}")
                    
                    total_records_uploaded += success_count
                else:
                    # Simulation mode - just count batches
                    num_batches = (record_count + args.batch_size - 1) // args.batch_size
                    log_entries.append(
                        f"  Would upload in {num_batches} batch(es) (SIMULATION)"
                    )
                    total_records_uploaded += record_count
            else:
                log_entries.append("  No valid records found in file after cleaning")
                
        except pl.exceptions.ComputeError as e:
            log_entries.append(f"  ERROR: Polars compute error - {str(e)}")
        except Exception as e:
            log_entries.append(f"  ERROR: Failed to process file - {str(e)}")
            import traceback
            log_entries.append(f"  Traceback: {traceback.format_exc()}")

    # Summary statistics using Polars
    log_entries += [
        "\n" + "=" * 80,
        "BLAST DATA UPLOAD SUMMARY",
        f"Files processed: {len(csv_files)}",
        f"Total records processed: {total_records_processed}",
        f"Total records uploaded: {total_records_uploaded}"
        if upload_enabled
        else f"Total records that would be uploaded: {total_records_uploaded}",
    ]
    
    # Add aggregate statistics if available
    if all_stats:
        log_entries.append("\nAGGREGATE STATISTICS:")
        total_unique_samples = len(set(
            sample 
            for stats in all_stats.values() 
            for sample in [stats.get("unique_samples", 0)]
        ))
        log_entries.append(f"  Total unique samples across all files: {total_unique_samples}")
        
        # Calculate overall averages
        if any("avg_pident" in stats for stats in all_stats.values()):
            avg_pidents = [stats.get("avg_pident", 0) for stats in all_stats.values() if "avg_pident" in stats]
            if avg_pidents:
                overall_avg_pident = sum(avg_pidents) / len(avg_pidents)
                log_entries.append(f"  Overall average percent identity: {overall_avg_pident:.2f}%")
    
    log_entries.append(
        "\nBLAST UPLOAD COMPLETE"
        if upload_enabled
        else "\nBLAST SIMULATION COMPLETE - No actual upload performed"
    )

    # Write log file
    with open("blast_labkey_upload.log", "w") as f:
        f.write("\n".join(log_entries))
    
    # Also print to console
    print("\n".join(log_entries))


if __name__ == "__main__":
    main()