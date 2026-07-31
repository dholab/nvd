#!/usr/bin/env python3

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import polars as pl
except ImportError:
    print("ERROR: Polars is not installed. Please install it with: pip install polars")
    sys.exit(1)


# Taxid columns held strictly to a clean integer or a passable null so LabKey's
# integer columns stay strict. A row whose taxid is neither (e.g. a stray
# "11676;11706" or a label) is skipped before upload and logged, never coerced
# or sent as something the schema would reject. Absent taxids reach LabKey as a
# true SQL NULL (never a fabricated 0 or "").
TAXID_COLUMNS = {"staxids", "adjusted_taxid"}

# Identifier columns a row must have to be uploadable; a row missing any of them
# is skipped and logged rather than sent for LabKey to reject.
REQUIRED_COLUMNS = {"experiment", "sample_id", "qseqid"}


def drop_invalid_rows(
    df: pl.DataFrame,
    taxid_columns: set[str],
    required_columns: set[str],
    log_entries: list[str],
) -> pl.DataFrame:
    """Drop rows LabKey's strict schema would reject, logging each one.

    A row is skipped when a taxid column holds a value that is neither a clean
    integer nor null, or when a required identifier column is null/empty. Kept
    rows have their taxid columns cast to integers; a genuinely absent taxid
    stays null (a passable null LabKey accepts). Skipped rows are recorded in
    ``log_entries`` and are not uploaded, so a strict integer column never has
    to reject them at insert time.
    """
    taxids = [c for c in sorted(taxid_columns) if c in df.columns]
    required = [c for c in sorted(required_columns) if c in df.columns]

    # Per taxid column: normalized text, "present" (empty -> null), parsed int.
    # A value is bad if it is present but does not parse as an integer.
    text = {c: pl.col(c).cast(pl.Utf8, strict=False).str.strip_chars() for c in taxids}
    present = {
        c: pl.when(text[c].str.len_chars() == 0).then(None).otherwise(text[c])
        for c in taxids
    }
    as_int = {
        c: present[c].str.to_integer(strict=False).cast(pl.Int64, strict=False)
        for c in taxids
    }

    annotations = {f"__present_{c}": present[c] for c in taxids}
    annotations.update(
        {f"__bad_{c}": present[c].is_not_null() & as_int[c].is_null() for c in taxids},
    )
    annotations.update(
        {
            f"__missing_{c}": (
                pl.col(c).is_null()
                | (
                    pl.col(c).cast(pl.Utf8, strict=False).str.strip_chars().str.len_chars()
                    == 0
                )
            )
            for c in required
        },
    )

    if not annotations:
        return df

    annotated = df.with_columns(**annotations)
    flag_cols = [f"__bad_{c}" for c in taxids] + [f"__missing_{c}" for c in required]
    row_invalid = pl.any_horizontal([pl.col(flag) for flag in flag_cols])

    invalid = annotated.filter(row_invalid)
    for row in invalid.iter_rows(named=True):
        reasons = [
            f"{c}={row['__present_' + c]!r} (not an integer)"
            for c in taxids
            if row[f"__bad_{c}"]
        ]
        reasons += [f"{c} missing" for c in required if row[f"__missing_{c}"]]
        ident = ", ".join(
            f"{key}={row.get(key)!r}" for key in ("sample_id", "qseqid") if key in row
        )
        log_entries.append(f"  SKIPPED ROW ({ident}): {'; '.join(reasons)}")

    if invalid.height > 0:
        log_entries.append(
            f"  SKIPPED {invalid.height} row(s) that would violate the LabKey "
            f"schema; not uploaded (see entries above).",
        )

    return (
        annotated.filter(~row_invalid)
        .with_columns(**{c: as_int[c] for c in taxids})
        .drop(list(annotations.keys()))
    )


def validate_dataframe(
    df: pl.DataFrame,
    csv_file: str,
    log_entries: list[str],
) -> pl.DataFrame:
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
    log_entries.append(
        f"  Columns: {', '.join(df.columns[:10])}{'...' if len(df.columns) > 10 else ''}",
    )

    # Check for critical columns and log missing ones
    critical_columns = ["experiment", "sample_id", "qseqid"]
    missing_columns = [col for col in critical_columns if col not in df.columns]
    if missing_columns:
        log_entries.append(f"  WARNING: Missing critical columns: {missing_columns}")

    # Convert numeric columns to proper types if they exist
    numeric_columns = {
        "length": pl.Int64,
        "pident": pl.Float64,
        "evalue": pl.Float64,
        "bitscore": pl.Float64,
        "mapped_reads": pl.Int64,
        "total_reads": pl.Int64,
        "experiment": pl.Int64,
    }

    for col, dtype in numeric_columns.items():
        if col in df.columns:
            try:
                # Try to cast to the appropriate type
                df = df.with_columns(pl.col(col).cast(dtype, strict=False))
                log_entries.append(f"  Converted '{col}' to {dtype}")
            except Exception as e:
                log_entries.append(
                    f"  WARNING: Could not convert '{col}' to {dtype}: {e!s}",
                )

    # Skip (and log) rows LabKey's strict schema would reject: a taxid that is
    # neither a clean integer nor a passable null, or a missing identifier. Kept
    # rows carry integer taxids or a true null.
    df = drop_invalid_rows(df, TAXID_COLUMNS, REQUIRED_COLUMNS, log_entries)

    # Fill null values with appropriate defaults
    for col in df.columns:
        if df[col].dtype in [pl.Float32, pl.Float64]:
            df = df.with_columns(pl.col(col).fill_null(0.0))
        elif (
            df[col].dtype in [pl.Int8, pl.Int16, pl.Int32, pl.Int64]
            and col not in TAXID_COLUMNS
        ):
            df = df.with_columns(pl.col(col).fill_null(0))
        elif df[col].dtype == pl.Utf8 and col not in TAXID_COLUMNS:
            df = df.with_columns(pl.col(col).fill_null(""))

    return df


def apply_reference_cutoff(df: pl.DataFrame, top_k: int) -> pl.DataFrame:
    """Keep only the top_k highest-bitscore reference sequences per (sample_id, qseqid).

    core_nt contains many redundant reference sequences, so a single contig can match
    dozens of references at an identical bitscore. Inserting every tied reference bloats
    the LabKey production list. This bounds the list to the top_k references per
    (sample_id, qseqid) by descending bitscore, breaking ties deterministically on
    sseqid; every taxid row of a retained reference is kept.

    Only the production-list insert is trimmed here. The raw enriched BLAST result is
    uploaded to WebDAV separately (upstream of this step) and is left untouched.

    Args:
        df: Cleaned BLAST DataFrame (expects sample_id, qseqid, sseqid, bitscore).
        top_k: Maximum reference sequences to keep per (sample_id, qseqid).

    Returns:
        Rows whose reference is within the top_k for its group. Returned unchanged when
        top_k <= 0 or the required columns are absent.
    """
    required_columns = {"sample_id", "qseqid", "sseqid", "bitscore"}
    if top_k <= 0 or not required_columns.issubset(df.columns):
        return df

    retained_references = (
        df.group_by(["sample_id", "qseqid", "sseqid"])
        .agg(pl.col("bitscore").max().alias("_bitscore"))
        # Best bitscore first, then a deterministic sseqid tiebreak so a truncated tie
        # is reproducible run-to-run.
        .sort(
            ["sample_id", "qseqid", "_bitscore", "sseqid"],
            descending=[False, False, True, False],
        )
        # 0-based ordinal position within each (sample_id, qseqid) group.
        .with_columns(
            pl.int_range(pl.len()).over(["sample_id", "qseqid"]).alias("_position"),
        )
        .filter(pl.col("_position") < top_k)
        .select(["sample_id", "qseqid", "sseqid"])
    )

    return df.join(
        retained_references,
        on=["sample_id", "qseqid", "sseqid"],
        how="semi",
    )


def get_sample_stats(df: pl.DataFrame) -> dict[str, Any]:
    """
    Get statistical summary of the DataFrame.

    Args:
        df: Polars DataFrame

    Returns:
        Dictionary with statistics
    """
    stats = {
        "total_rows": len(df),
        "unique_samples": df["sample_id"].n_unique()
        if "sample_id" in df.columns
        else 0,
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


def dataframe_to_records(df: pl.DataFrame) -> list[dict[str, Any]]:
    """
    Convert Polars DataFrame to list of dictionaries for LabKey upload.

    Args:
        df: Polars DataFrame

    Returns:
        List of dictionaries
    """
    # Convert to dictionaries, handling null values
    records = df.to_dicts()

    # Clean up records for LabKey. Taxid columns preserve absence as a true null
    # (LabKey stores SQL NULL); every other column keeps the legacy empty-string
    # form for a missing value.
    cleaned_records = []
    for record in records:
        cleaned_record = {}
        for key, value in record.items():
            is_missing = value is None or (isinstance(value, float) and value != value)
            if is_missing:
                cleaned_record[key] = None if key in TAXID_COLUMNS else ""
            else:
                cleaned_record[key] = value
        cleaned_records.append(cleaned_record)

    return cleaned_records


def combo_already_uploaded(
    query_api,
    schema: str,
    table: str,
    experiment,
    sample_id: str,
    query_class: str,
) -> bool:
    """True if the hits list already holds any row for this combo.

    The destination list is its own completion ledger: presence of any row for
    (experiment, sample_id, query_class) means "already uploaded", so the caller
    skips the insert rather than risk a duplicate. Reuses the guard script's
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
                QueryFilter("query_class", query_class, "eq"),
            ],
            max_rows=1,
        )
    except (ImportError, AttributeError):
        result = query_api.select_rows(
            schema_name=schema,
            query_name=table,
            filter_array=[
                f"experiment~eq={experiment}",
                f"sample_id~eq={sample_id}",
                f"query_class~eq={query_class}",
            ],
        )
    return bool(result and result.get("rows"))


def insert_records(
    query_api,
    schema: str,
    table: str,
    records: list[dict[str, Any]],
) -> None:
    """Insert every record for a (sample_id, query_class) combo in one atomic call.

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
    with open("blast_labkey_upload.log", "w") as f:
        f.write(text)
    print(text)


def main():
    parser = argparse.ArgumentParser(
        description="Upload a (sample, query_class) BLAST CSV batch to LabKey using Polars.",
    )
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--query-class", required=True)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--labkey-server", required=True)
    parser.add_argument("--labkey-project-name", required=True)
    parser.add_argument("--labkey-api-key", required=True)
    parser.add_argument("--labkey-schema", required=True)
    parser.add_argument("--table-name", default="metagenomic_hits_test_nvd2")
    parser.add_argument(
        "--blast-retention-count",
        type=int,
        default=5,
        help="Max reference sequences to keep per (sample, contig) in the "
        "production list (default: 5). Trims only the insert, not the WebDAV upload.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate data without uploading",
    )
    args = parser.parse_args()

    log_entries = [
        f"LabKey BLAST Upload Log (Polars Version) - {datetime.now()}",
        f"Experiment ID: {args.experiment_id}",
        f"Sample: {args.sample_id}",
        f"Query class: {args.query_class}",
        f"Server: {args.labkey_server}",
        f"Project: {args.labkey_project_name}",
        f"Target Table: {args.table_name}",
        "=" * 80,
    ]

    upload_enabled = (
        bool(args.labkey_server and args.labkey_project_name and args.labkey_api_key)
        and not args.validate_only
    )

    api = None
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
    elif args.validate_only:
        log_entries.append("Running in validation-only mode")
    else:
        log_entries.append(
            "No LabKey credentials provided - running in simulation mode",
        )

    # The destination list is its own completion ledger. If this combo already
    # has rows there, a prior run already uploaded it: skip re-inserting rather
    # than risk duplicating the hits list.
    if upload_enabled and combo_already_uploaded(
        api.query,
        args.labkey_schema,
        args.table_name,
        int(args.experiment_id),
        args.sample_id,
        args.query_class,
    ):
        log_entries.append(
            f"SKIP: combo already uploaded (exp={args.experiment_id}, "
            f"sample={args.sample_id}, query_class={args.query_class}); no insert.",
        )
        _write_log(log_entries)
        return

    log_entries.append(f"\nProcessing BLAST file: {args.csv}")

    csv_path = Path(args.csv)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        log_entries.append("  Empty or missing file - skipping")
        log_entries.append(
            "\nBLAST UPLOAD COMPLETE"
            if upload_enabled
            else "\nBLAST SIMULATION COMPLETE - No actual upload performed",
        )
        _write_log(log_entries)
        return

    try:
        # Read CSV with Polars
        df = pl.read_csv(
            args.csv,
            separator="\t" if args.csv.endswith(".tsv") else ",",
            has_header=True,
            ignore_errors=True,  # Continue on parsing errors
            null_values=["", "NA", "N/A", "null", "NULL", "None"],
            infer_schema_length=10000,  # Infer types from more rows
        )

        # Validate and clean the DataFrame
        df = validate_dataframe(df, args.csv, log_entries)

        # Trim the production list to the top references per contig. core_nt
        # redundancy can otherwise insert dozens of tied references per contig;
        # this bounds the inserted rows without touching the WebDAV upload.
        rows_before_cutoff = len(df)
        df = apply_reference_cutoff(df, args.blast_retention_count)
        if len(df) != rows_before_cutoff:
            log_entries.append(
                f"  Production-list cutoff (<= {args.blast_retention_count} "
                f"references per contig): {rows_before_cutoff} -> {len(df)} rows",
            )

        record_count = len(df)

        if record_count > 0:
            # Get statistics
            stats = get_sample_stats(df)

            log_entries.append(f"  Records: {record_count}")
            log_entries.append("  Statistics:")
            for key, value in stats.items():
                if isinstance(value, float):
                    log_entries.append(f"    {key}: {value:.4f}")
                else:
                    log_entries.append(f"    {key}: {value}")

            if upload_enabled:
                # Insert the whole combo atomically: one (sample_id, query_class)
                # batch is small enough that chunked inserts only add partial-
                # failure risk without a throughput benefit.
                records = dataframe_to_records(df)
                try:
                    insert_records(
                        api.query, args.labkey_schema, args.table_name, records,
                    )
                except Exception as e:
                    log_entries.append(f"  Upload: ERROR - {e!s}")
                    log_entries.append(
                        f"    Failed batch first record: {df.head(1).to_dicts()[0]}",
                    )
                    log_entries.append(
                        "\nBLAST UPLOAD FAILED - insert error, see above (no rows committed)",
                    )
                    _write_log(log_entries)
                    print(
                        f"ERROR: LabKey insert failed for sample={args.sample_id}, "
                        f"query_class={args.query_class}: {e!s}",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                log_entries.append(f"  Upload: SUCCESS ({len(records)} records)")
            else:
                log_entries.append(
                    f"  Would upload {record_count} record(s) (SIMULATION)",
                )
        else:
            log_entries.append("  No valid records found in file after cleaning")

    except pl.exceptions.ComputeError as e:
        log_entries.append(f"  ERROR: Polars compute error - {e!s}")
    except Exception as e:
        log_entries.append(f"  ERROR: Failed to process file - {e!s}")
        import traceback

        log_entries.append(f"  Traceback: {traceback.format_exc()}")

    log_entries.append(
        "\nBLAST UPLOAD COMPLETE"
        if upload_enabled
        else "\nBLAST SIMULATION COMPLETE - No actual upload performed",
    )

    _write_log(log_entries)


if __name__ == "__main__":
    main()
