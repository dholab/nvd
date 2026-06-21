#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars>=1.27.1"]
# ///
"""Prepare WVDB inputs for a combined sourmash reference database."""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

WVDB_ID_PREFIX = "WVDB|"
LINEAGE_COLUMNS = (
    "ident",
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


def ensure_parent(path: Path) -> None:
    """Create the parent directory for an output path when needed."""
    path.parent.mkdir(parents=True, exist_ok=True)


def cleaned_string(name: str) -> pl.Expr:
    """Return a stripped string column expression with blanks converted to null."""
    value = pl.col(name).cast(pl.Utf8).str.strip_chars()
    return pl.when(value.str.len_chars() > 0).then(value).otherwise(None)


def first_nonempty(available_columns: set[str], *columns: str) -> pl.Expr:
    """Return the first non-empty value among columns present in a frame."""
    values = [cleaned_string(name) for name in columns if name in available_columns]
    if not values:
        return pl.lit("")
    return pl.coalesce(values).fill_null("")


def normalize_fasta_ids(source: Path, destination: Path) -> None:
    """Write a FASTA with WVDB-prefixed, whitespace-free sequence identifiers."""
    ensure_parent(destination)
    seen: set[str] = set()

    with (
        source.open(encoding="utf-8") as src,
        destination.open("w", encoding="utf-8") as dst,
    ):
        for line_number, line in enumerate(src, start=1):
            if not line.startswith(">"):
                dst.write(line)
                continue

            raw_id = line[1:].strip().split(maxsplit=1)[0]
            if not raw_id:
                msg = f"Empty FASTA identifier at line {line_number} in {source}"
                raise ValueError(msg)

            normalized_id = f"{WVDB_ID_PREFIX}{raw_id}"
            if normalized_id in seen:
                msg = f"Duplicate FASTA identifier after normalization: {normalized_id}"
                raise ValueError(msg)
            seen.add(normalized_id)

            dst.write(f">{normalized_id}\n")


def scan_csv(path: Path, *, separator: str = ",", skip_rows: int = 0) -> pl.LazyFrame:
    """Read a CSV-like file lazily with all columns as strings."""
    return pl.scan_csv(
        path,
        separator=separator,
        skip_rows=skip_rows,
        infer_schema_length=0,
        null_values=[""],
    )


def require_lazy_columns(
    frame: pl.LazyFrame,
    *,
    path: Path,
    columns: set[str],
) -> pl.LazyFrame:
    """Fail if a lazy CSV frame lacks required columns."""
    missing = columns.difference(frame.collect_schema().names())
    if missing:
        msg = f"{path} is missing required column(s): {', '.join(sorted(missing))}"
        raise ValueError(msg)
    return frame


def scan_annotations(annotations_tsv: Path) -> pl.LazyFrame:
    """Read the WVDB annotation TSV lazily as strings."""
    return scan_csv(
        annotations_tsv,
        separator="\t",
    ).pipe(require_lazy_columns, path=annotations_tsv, columns={"contig"})


def scan_sourmash_manifest(manifest_csv: Path) -> pl.LazyFrame:
    """Read a sourmash manifest CSV, skipping its metadata preamble."""
    return scan_csv(manifest_csv, skip_rows=1).pipe(
        require_lazy_columns,
        path=manifest_csv,
        columns={"name"},
    )


def require_unique_identifiers(lineages: pl.DataFrame) -> pl.DataFrame:
    """Fail if lineage identifiers are duplicated."""
    duplicates = (
        lineages.group_by("ident")
        .len()
        .filter(pl.col("len") > 1)
        .select("ident")
        .get_column("ident")
        .to_list()
    )
    if duplicates:
        examples = ", ".join(duplicates[:10])
        msg = f"Duplicate WVDB lineage identifiers: {examples}"
        raise ValueError(msg)
    return lineages


def select_wvdb_lineages(annotations: pl.LazyFrame) -> pl.LazyFrame:
    """Project WVDB annotation columns into sourmash lineage columns."""
    available_columns = set(annotations.collect_schema().names())
    return annotations.select(
        (
            pl.lit(WVDB_ID_PREFIX) + pl.col("contig").cast(pl.Utf8).str.strip_chars()
        ).alias("ident"),
        pl.lit("Viruses").alias("superkingdom"),
        first_nonempty(available_columns, "Phylum_RdRp", "Phylum_gNd").alias("phylum"),
        first_nonempty(available_columns, "Class_RdRp", "Class_gNd").alias("class"),
        first_nonempty(
            available_columns,
            "Order_consensus",
            "Order_RdRp",
            "Order_gNd",
        ).alias("order"),
        first_nonempty(
            available_columns,
            "Family_consensus",
            "Family_RdRp",
            "Family_gNd",
            "Family_ICTV",
        ).alias("family"),
        first_nonempty(available_columns, "Genus_RdRp").alias("genus"),
        first_nonempty(available_columns, "Species_RdRp").alias("species"),
    )


def build_wvdb_lineages(annotations: pl.LazyFrame) -> pl.DataFrame:
    """Build WVDB sourmash lineages as a validated DataFrame."""
    return (
        annotations.pipe(select_wvdb_lineages)
        .collect()
        .pipe(require_unique_identifiers)
    )


def write_wvdb_lineages(annotations_tsv: Path, lineages_csv: Path) -> None:
    """Write sourmash-compatible WVDB lineages from the WVDB annotation TSV."""
    ensure_parent(lineages_csv)
    build_wvdb_lineages(scan_annotations(annotations_tsv)).write_csv(lineages_csv)


def prepare_inputs(
    *,
    fasta: Path,
    annotations_tsv: Path,
    normalized_fasta: Path,
    lineages_csv: Path,
) -> None:
    """Write normalized FASTA and sourmash lineage CSV outputs."""
    normalize_fasta_ids(fasta, normalized_fasta)
    write_wvdb_lineages(annotations_tsv, lineages_csv)


def find_missing_wvdb_lineages(
    manifest: pl.LazyFrame,
    lineages: pl.LazyFrame,
) -> pl.DataFrame:
    """Return WVDB manifest names that do not have matching lineages."""
    return (
        manifest.filter(pl.col("name").str.starts_with(WVDB_ID_PREFIX))
        .select("name")
        .join(lineages.select("ident"), left_on="name", right_on="ident", how="anti")
        .collect()
    )


def fail_on_missing_wvdb_lineages(missing: pl.DataFrame) -> pl.DataFrame:
    """Fail when the missing-lineage frame is not empty."""

    if missing.height:
        examples = ", ".join(missing.get_column("name").head(10).to_list())
        msg = f"{missing.height} WVDB signatures are missing lineages. Examples: {examples}"
        raise ValueError(msg)
    return missing


def check_manifest_coverage(manifest_csv: Path, lineages_csv: Path) -> None:
    """Verify that every WVDB signature in a sourmash manifest has a lineage."""
    manifest = scan_sourmash_manifest(manifest_csv)
    lineages = scan_csv(lineages_csv).pipe(
        require_lazy_columns,
        path=lineages_csv,
        columns={"ident"},
    )

    find_missing_wvdb_lineages(manifest, lineages).pipe(fail_on_missing_wvdb_lineages)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Prepare WVDB inputs for a combined NCBI Virus + WVDB sourmash reference.",
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    prepare = subcommands.add_parser(
        "prepare-inputs",
        help="Normalize WVDB FASTA IDs and write WVDB sourmash lineages.",
    )
    prepare.add_argument("--fasta", type=Path, required=True)
    prepare.add_argument("--annotations-tsv", type=Path, required=True)
    prepare.add_argument("--normalized-fasta", type=Path, required=True)
    prepare.add_argument("--lineages-csv", type=Path, required=True)

    check = subcommands.add_parser(
        "check-manifest-coverage",
        help="Fail if any WVDB signatures in a sourmash manifest lack lineages.",
    )
    check.add_argument("--manifest-csv", type=Path, required=True)
    check.add_argument("--lineages-csv", type=Path, required=True)

    return parser


def main() -> None:
    """Run the command-line interface."""
    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.command == "prepare-inputs":
            prepare_inputs(
                fasta=args.fasta,
                annotations_tsv=args.annotations_tsv,
                normalized_fasta=args.normalized_fasta,
                lineages_csv=args.lineages_csv,
            )
        elif args.command == "check-manifest-coverage":
            check_manifest_coverage(args.manifest_csv, args.lineages_csv)
    except ValueError as exc:
        parser.exit(1, f"error: {exc}\n")


if __name__ == "__main__":
    main()
