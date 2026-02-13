#!/usr/bin/env python3
# ruff: noqa: TC003
"""
labkey_upload_gottcha_fasta.py — keep one taxonomic level per read: STRAIN if present, else SPECIES

Why this shape?
---------------
GOTTCHA2 duplicates each hit across multiple taxonomy levels using the *same read header*.
We just want all the reads, but only at the lowest relevant level:
    - If a read has any STRAIN entries -> keep *only* those STRAIN entries for that read
    - Else, if it has SPECIES -> keep those SPECIES entries
    - Ignore other levels

This script:
  1) Reads a plain-text FASTA (file path only; no stdin, no compression).
  2) Keeps only LEVEL ∈ {species, strain}, selecting per-read as above.
  3) Writes a LabKey TSV with columns:
       "Experiment", "Sample Id", "Header", "Sequence", "Notes", "Snakemake Run Id", "Upload Type"
     (Header is the FASTA header without the leading '>')
  4) Uploads TSV rows to a LabKey list in **batches** (unless --no_upload).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections.abc import Generator, Iterable, Sequence
from typing import TextIO

# Lazy import so --no_upload works without LabKey client installed
try:
    from labkey.api_wrapper import APIWrapper
except Exception:
    APIWrapper = None  # ty:ignore[invalid-assignment]


# ----------------------------- FASTA parsing -----------------------------


def iterate_fasta(handle: TextIO) -> Generator[tuple[str, list[str]], None, None]:
    """
    Yield (header_line_including_>, [sequence_lines]) for each FASTA record.
    """
    current_header: str | None = None
    current_seq: list[str] = []
    for raw in handle:
        line = raw.rstrip("\r\n")
        if line.startswith(">"):
            if current_header is not None:
                yield current_header, current_seq
            current_header = line
            current_seq = []
        else:
            if current_header is None:
                continue  # ignore stray sequence lines before first header
            current_seq.append(line)
    if current_header is not None:
        yield current_header, current_seq


LEVEL_RE = re.compile(r"LEVEL=(\w+)", re.IGNORECASE)


def extract_level(header: str) -> str | None:
    """
    Return 'species' or 'strain' if present as LEVEL=..., else None.
    """
    m = LEVEL_RE.search(header)
    return m.group(1).lower() if m else None


def read_key_from_header(header: str) -> str:
    """
    GOTTCHA2 duplicates the SAME read across levels with the same leading token.
    Read key is everything after '>' up to the first '|'.
    """
    h = header.lstrip(">")
    return h.split("|", 1)[0]


# ------------------------ Selection (per-read) ---------------------------


def select_records_for_upload(fasta_path: str) -> list[tuple[str, str]]:
    """
    Parse the FASTA and return a list of (header_without_gt, concatenated_sequence)
    keeping ONLY species/strain entries and preferring strain per read.

    Algorithm:
      - Map: read_key -> {'strain': [(hdr, seq)], 'species': [(hdr, seq)]}
      - For each read_key:
            if any 'strain' -> take ALL 'strain' entries for that read
            elif any 'species' -> take ALL 'species' entries
            else -> take nothing
    """
    per_read: dict[str, dict[str, list[tuple[str, str]]]] = {}

    # 1) Build the map
    with open(fasta_path, encoding="utf-8") as fh:
        for header, seq_lines in iterate_fasta(fh):
            level = extract_level(header)
            if level not in {"species", "strain"}:
                continue  # ignore other taxonomic levels entirely

            read_key = read_key_from_header(header)
            seq = "".join(seq_lines)
            if not seq:
                continue

            bucket = per_read.setdefault(read_key, {"strain": [], "species": []})
            bucket[level].append((header.lstrip(">"), seq))

    # 2) For each read, keep only the desired level
    selected: list[tuple[str, str]] = []
    for bucket in per_read.values():
        if bucket["strain"]:
            selected.extend(bucket["strain"])
        elif bucket["species"]:
            selected.extend(bucket["species"])
        # else: neither level present -> nothing for this read

    return selected


# --------------------------- LabKey batching -----------------------------


def upload_in_batches(
    lk: APIWrapper,
    list_name: str,
    rows_iterable: Iterable[dict],
    batch_size: int,
) -> tuple[int, int]:
    """
    Upload rows to LabKey in batches. Returns (total_inserted, num_batches).
    """
    batch: list[dict] = []
    total = 0
    batches = 0

    def flush():
        nonlocal batch, total, batches
        if not batch:
            return
        resp = lk.query.insert_rows("lists", list_name, batch)
        inserted = len(resp.get("rows", [])) if isinstance(resp, dict) else 0
        total += inserted
        batches += 1
        for r in resp.get("rows", []):
            key = r.get("Key", "N/A")
            samp = r.get("Sample Id", "N/A")
            print(f"   - Inserted: Sample Id={samp}, Key={key}")
        batch = []

    for row in rows_iterable:
        batch.append(row)
        if len(batch) >= batch_size:
            flush()
    flush()
    return total, batches


# ------------------------------- CLI ------------------------------------

TSV_COLUMNS = [
    "Experiment",
    "Sample Id",
    "Header",
    "Sequence",
    "Notes",
    "Snakemake Run Id",
    "Upload Type",
]


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Upload GOTTCHA2 FASTA hits at lowest relevant level per read (strain if present, else species).",
    )
    p.add_argument(
        "--fasta",
        required=True,
        help="Path to input FASTA (plain text; no stdin).",
    )
    p.add_argument("--sample_id", required=True, help="Value for 'Sample Id'.")
    p.add_argument("--experiment_id", required=True, help="Value for 'Experiment'.")
    p.add_argument(
        "--sample_set_id",
        required=False,
        help="Sample set ID for upload tracking",
    )
    p.add_argument(
        "--state_dir",
        required=False,
        help="State directory for upload tracking database",
    )
    p.add_argument(
        "--output_tsv",
        required=True,
        help="Path to write TSV (file path only).",
    )

    # LabKey settings (omit with --no_upload)
    p.add_argument("--server", help="LabKey server URL")
    p.add_argument("--container", help="LabKey container/project path")
    p.add_argument("--list", dest="list_name", help="LabKey list name")
    p.add_argument("--api_key", help="LabKey API key")
    p.add_argument("--batch_size", type=int, default=1000, help="Upload batch size")
    p.add_argument(
        "--no_upload",
        action="store_true",
        help="Only write TSV; skip upload",
    )

    # Fixed column defaults
    p.add_argument("--notes", default="gottcha2", help="Value for 'Notes'")
    p.add_argument(
        "--run_id",
        default="GOTTCHA2_UPLOAD",
        help="Value for 'Snakemake Run Id'",
    )
    p.add_argument("--upload_type", default="fasta", help="Value for 'Upload Type'")

    args = p.parse_args(argv)

    # Validate file paths
    if args.fasta.strip() in {"-", ""}:
        print(
            "[❌ ERROR] --fasta must be a real file path (no stdin).",
            file=sys.stderr,
        )
        return 1
    if not os.path.isfile(args.fasta):
        print(f"[❌ ERROR] FASTA not found: {args.fasta}", file=sys.stderr)
        return 1
    if args.output_tsv.strip() in {"-", ""}:
        print(
            "[❌ ERROR] --output_tsv must be a file path (no stdout).",
            file=sys.stderr,
        )
        return 1

    # Check for cross-run duplicate uploads if sample_set_id is provided
    # Justification for local import: script must work standalone without py_nvd
    if args.sample_set_id and not args.no_upload:
        try:
            import py_nvd.state as nvd_state

            if nvd_state.was_sample_ever_uploaded(
                args.sample_id,
                upload_type="gottcha2_fasta",
                upload_target="labkey",
                state_dir=args.state_dir,
            ):
                print(
                    f"[INFO] SKIP: Sample '{args.sample_id}' was already uploaded to LabKey in a previous run.",
                    file=sys.stderr,
                )
                return 0
        except ImportError:
            print(
                "[WARN] py_nvd not available, skipping cross-run check",
                file=sys.stderr,
            )

    try:
        # Select records exactly once
        selected_records = select_records_for_upload(args.fasta)
        print(f"[INFO] Selected records: {len(selected_records)}", file=sys.stderr)

        # Write TSV
        with open(args.output_tsv, "w", encoding="utf-8", newline="") as out_fh:
            writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
            writer.writerow(TSV_COLUMNS)

            # Optional: set up LabKey client
            lk = None
            if not args.no_upload:
                if APIWrapper is None:
                    raise RuntimeError(
                        "labkey.api_wrapper not available. Install it or pass --no_upload.",
                    )
                missing = [
                    flag
                    for flag, val in {
                        "--server": args.server,
                        "--container": args.container,
                        "--list": args.list_name,
                        "--api_key": args.api_key,
                    }.items()
                    if not val
                ]
                if missing:
                    raise ValueError(f"Missing LabKey parameters: {', '.join(missing)}")
                lk = APIWrapper(args.server, args.container, api_key=args.api_key)  # type: ignore

            # Build all rows for TSV and upload
            all_rows: list[dict] = []
            for header_no_gt, seq in selected_records:
                row = {
                    "Experiment": args.experiment_id,
                    "Sample Id": args.sample_id,
                    "Header": header_no_gt,
                    "Sequence": seq,
                    "Notes": args.notes,
                    "Snakemake Run Id": args.run_id,
                    "Upload Type": args.upload_type,
                }
                writer.writerow([row[c] for c in TSV_COLUMNS])
                all_rows.append(row)

            total_rows = len(all_rows)

            # Compute content hash for idempotency tracking
            # Justification for local import: script must work standalone without py_nvd
            content_hash = None
            if args.sample_set_id:
                try:
                    import py_nvd.state as nvd_state

                    content_hash = nvd_state.hash_upload_content(all_rows)
                    print(
                        f"[INFO] Content hash: {content_hash[:16]}...",
                        file=sys.stderr,
                    )
                except ImportError:
                    pass

            # Upload to LabKey
            total_uploaded = 0
            if lk is not None:
                uploaded, batches = upload_in_batches(
                    lk,
                    args.list_name,
                    all_rows,
                    args.batch_size,
                )
                total_uploaded = uploaded
                print(
                    f"[INFO] Uploaded {uploaded} rows in {batches} batch(es).",
                    file=sys.stderr,
                )

                # Record successful upload in state database
                # Justification for local import: script must work standalone without py_nvd
                if args.sample_set_id and content_hash and total_uploaded > 0:
                    try:
                        import py_nvd.state as nvd_state

                        nvd_state.record_upload(
                            sample_id=args.sample_id,
                            sample_set_id=args.sample_set_id,
                            upload_type="gottcha2_fasta",
                            upload_target="labkey",
                            content_hash=content_hash,
                            target_metadata={
                                "experiment_id": args.experiment_id,
                                "list_name": args.list_name,
                                "records_uploaded": total_uploaded,
                            },
                            state_dir=args.state_dir,
                        )
                        # Mark sample as uploaded (terminal state for LabKey runs)
                        nvd_state.mark_sample_uploaded(
                            sample_id=args.sample_id,
                            sample_set_id=args.sample_set_id,
                            state_dir=args.state_dir,
                        )
                        print(
                            f"[INFO] Recorded upload in state database for sample '{args.sample_id}'",
                            file=sys.stderr,
                        )
                    except ImportError:
                        print(
                            "[WARN] py_nvd.state not available, upload not recorded",
                            file=sys.stderr,
                        )
                    except Exception as e:
                        print(f"[WARN] Failed to record upload: {e}", file=sys.stderr)

        print(f"[INFO] TSV rows written: {total_rows}", file=sys.stderr)
        if args.no_upload:
            print("[INFO] Upload skipped (--no_upload).", file=sys.stderr)
        else:
            print("[✅] Upload complete.", file=sys.stderr)
        return 0

    except Exception as e:
        print(f"[❌ ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
