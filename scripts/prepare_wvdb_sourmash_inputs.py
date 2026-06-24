#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars>=1.27.1"]
# ///
"""Prepare WVDB inputs for a combined sourmash reference database."""

from __future__ import annotations

import argparse
import csv
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import polars as pl

WVDB_ID_PREFIX = "WVDB|"
TAXONOMY_RANKS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
)
LINEAGE_COLUMNS = ("ident", *TAXONOMY_RANKS[:-1])
WVDB_TAXONOMY_RANKS = LINEAGE_COLUMNS[1:]
CURATED_LINEAGE_COLUMNS = ("ident", *TAXONOMY_RANKS)
LINEAGE_IDENTIFIER_COLUMNS = ("ident", "identifiers")
SIGNATURE_IDENTITY_COLUMNS = (
    "md5",
    "ksize",
    "moltype",
    "num",
    "scaled",
    "with_abundance",
)
TAXONOMY_CONFLICT_POLICIES = (
    "error",
    "most-specified",
    "ncbi-wins",
    "wvdb-wins",
)
UNKNOWN_TAXA = frozenset({"", "unknown", "na", "n/a", "none"})
DIAGNOSTIC_COLUMNS = (
    "conflict_type",
    "rank",
    "normalized_name",
    "policy",
    "decision",
    "source",
    "ident",
    "signature_identity",
    "specificity",
    "lineage",
)


@dataclass(frozen=True)
class CurationPaths:
    """Input and output paths for combined-reference taxonomy curation."""

    reference_manifest_csv: Path
    wvdb_manifest_csv: Path
    reference_lineages_csv: Path
    wvdb_lineages_csv: Path
    lineages_csv: Path
    diagnostics_tsv: Path


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


def normalize_lineage_identifier_column(
    frame: pl.LazyFrame,
    *,
    path: Path,
) -> pl.LazyFrame:
    """Expose raw and prepared sourmash lineage identifiers as ident."""
    available_columns = set(frame.collect_schema().names())
    for column in LINEAGE_IDENTIFIER_COLUMNS:
        if column not in available_columns:
            continue
        if column == "ident":
            return frame
        return frame.rename({column: "ident"})

    expected = ", ".join(LINEAGE_IDENTIFIER_COLUMNS)
    msg = f"{path} is missing required lineage identifier column (expected one of: {expected})"
    raise ValueError(msg)


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


def fill_internal_lineage_gaps(lineages: pl.LazyFrame) -> pl.LazyFrame:
    """Fill missing intermediate ranks when a lower rank is present."""
    expressions = []
    for index, rank in enumerate(WVDB_TAXONOMY_RANKS):
        lower_ranks = WVDB_TAXONOMY_RANKS[index + 1 :]
        if not lower_ranks:
            continue

        has_lower_rank = pl.any_horizontal(
            cleaned_string(lower_rank).is_not_null() for lower_rank in lower_ranks
        )
        missing_rank = cleaned_string(rank).is_null()
        expressions.append(
            pl.when(missing_rank & has_lower_rank)
            .then(pl.lit(f"unclassified {rank}"))
            .otherwise(pl.col(rank))
            .alias(rank),
        )

    return lineages.with_columns(expressions)


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
    ).pipe(fill_internal_lineage_gaps)


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
        normalize_lineage_identifier_column,
        path=lineages_csv,
    )

    find_missing_wvdb_lineages(manifest, lineages).pipe(fail_on_missing_wvdb_lineages)


def clean_cell(value: object) -> str:
    """Return a stripped string for CSV values that may be missing."""
    return "" if value is None else str(value).strip()


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    """Read a CSV file into string-valued dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return [
            {key: clean_cell(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def read_manifest_dicts(path: Path, *, source: str) -> list[dict[str, str]]:
    """Read sourmash manifest rows and attach their curation source."""
    with path.open(newline="", encoding="utf-8") as handle:
        first_line = handle.readline()
        if not first_line.startswith("# SOURMASH-MANIFEST-VERSION:"):
            handle.seek(0)
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or ())
        required = {"name", *SIGNATURE_IDENTITY_COLUMNS}
        missing = required.difference(columns)
        if missing:
            msg = f"{path} is missing required column(s): {', '.join(sorted(missing))}"
            raise ValueError(msg)
        records = []
        for row in reader:
            record = {key: clean_cell(value) for key, value in row.items()}
            record["ident"] = record["name"].split(maxsplit=1)[0]
            record["source"] = source
            records.append(record)
        return records


def read_lineage_dicts(path: Path, *, source: str) -> list[dict[str, str]]:
    """Read lineage rows into one curation schema with source metadata."""
    rows = []
    for raw in read_csv_dicts(path):
        ident = clean_cell(raw.get("ident") or raw.get("identifiers"))
        if not ident:
            msg = f"{path} contains a lineage row without ident/identifiers"
            raise ValueError(msg)
        row = {rank: clean_cell(raw.get(rank)) for rank in TAXONOMY_RANKS}
        row["ident"] = ident
        row["source"] = source
        rows.append(row)
    return rows


def lineage_by_identifier(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Index unique lineage records by identifier."""
    by_ident: dict[str, dict[str, str]] = {}
    for row in rows:
        ident = row["ident"]
        if ident in by_ident:
            msg = f"Duplicate lineage identifier: {ident}"
            raise ValueError(msg)
        by_ident[ident] = row
    return by_ident


def signature_identity(record: dict[str, str]) -> tuple[str, ...]:
    """Return the complete manifest-level identity of a sourmash sketch."""
    return tuple(record[column] for column in SIGNATURE_IDENTITY_COLUMNS)


def signature_identity_text(record: dict[str, str]) -> str:
    """Format a sketch identity for diagnostics."""
    return ";".join(
        f"{column}={record[column]}" for column in SIGNATURE_IDENTITY_COLUMNS
    )


def normalized_taxon_name(value: object) -> str:
    """Normalize WVDB/NCBI spelling differences for taxonomy comparison."""
    return " ".join(clean_cell(value).replace("_", " ").casefold().split())


def is_informative_taxon(value: object) -> bool:
    """Return whether a rank value provides more than unknown classification."""
    normalized = normalized_taxon_name(value)
    return normalized not in UNKNOWN_TAXA and not normalized.startswith("unclassified")


def lineage_prefix(row: dict[str, str], rank_index: int) -> tuple[str, ...]:
    """Return a normalized lineage prefix through the selected rank."""
    return tuple(
        normalized_taxon_name(row[rank]) for rank in TAXONOMY_RANKS[: rank_index + 1]
    )


def specificity(row: dict[str, str], rank_index: int) -> tuple[int, int]:
    """Score informative ranks through the taxon whose placement is contested."""
    informative = [
        index
        for index, rank in enumerate(TAXONOMY_RANKS[: rank_index + 1])
        if is_informative_taxon(row[rank])
    ]
    return len(informative), max(informative, default=-1)


def lineage_text(row: dict[str, str]) -> str:
    """Format a lineage for actionable diagnostics."""
    return " > ".join(row[rank] or "<blank>" for rank in TAXONOMY_RANKS)


def write_diagnostics(path: Path, rows: list[dict[str, str]]) -> None:
    """Write all automatic and unresolved curation decisions."""
    ensure_parent(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=DIAGNOSTIC_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def fail_on_duplicate_signature_identities(
    manifest_rows: list[dict[str, str]],
    lineages: dict[str, dict[str, str]],
    diagnostics: list[dict[str, str]],
) -> None:
    """Reject duplicate sketch identities instead of relying on input ordering."""
    grouped: defaultdict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for record in manifest_rows:
        grouped[signature_identity(record)].append(record)

    duplicate_groups = [records for records in grouped.values() if len(records) > 1]
    for records in duplicate_groups:
        for record in records:
            lineage = lineages.get(record["ident"], {})
            diagnostics.append(
                {
                    "conflict_type": "signature_identity",
                    "rank": "",
                    "normalized_name": "",
                    "policy": "error",
                    "decision": "ambiguous",
                    "source": record["source"],
                    "ident": record["ident"],
                    "signature_identity": signature_identity_text(record),
                    "specificity": "",
                    "lineage": lineage_text(lineage) if lineage else "<missing>",
                },
            )

    if duplicate_groups:
        records = duplicate_groups[0]
        msg = (
            "Duplicate sourmash signature identity: "
            f"{signature_identity_text(records[0])}; identifiers="
            f"{', '.join(record['ident'] for record in records)}; lineages="
            f"{' || '.join(lineage_text(lineages[record['ident']]) for record in records if record['ident'] in lineages)}"
        )
        raise ValueError(msg)


def preferred_source(policy: str) -> str | None:
    """Return the source selected by an explicit source-precedence policy."""
    if policy == "ncbi-wins":
        return "ncbi"
    if policy == "wvdb-wins":
        return "wvdb"
    return None


def select_canonical_lineage(
    rows: list[dict[str, str]],
    *,
    rank_index: int,
    policy: str,
) -> dict[str, str] | None:
    """Select one canonical lineage or return None for an unresolved conflict."""
    if policy == "error":
        return None

    candidates = rows
    source = preferred_source(policy)
    if source and any(row["source"] == source for row in rows):
        candidates = [row for row in rows if row["source"] == source]

    best_score = max(specificity(row, rank_index) for row in candidates)
    best = [row for row in candidates if specificity(row, rank_index) == best_score]
    prefixes = {lineage_prefix(row, rank_index) for row in best}
    if len(prefixes) != 1:
        return None
    return min(best, key=lambda row: row["ident"])


def taxonomy_conflicts_at_rank(
    rows: list[dict[str, str]],
    rank_index: int,
) -> list[tuple[str, list[dict[str, str]]]]:
    """Return taxon-name groups assigned to incompatible parents at one rank."""
    rank = TAXONOMY_RANKS[rank_index]
    grouped: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        name = normalized_taxon_name(row[rank])
        if name:
            grouped[name].append(row)
    return [
        (name, candidates)
        for name, candidates in grouped.items()
        if len({lineage_prefix(row, rank_index) for row in candidates}) > 1
    ]


def append_taxonomy_diagnostics(
    diagnostics: list[dict[str, str]],
    candidates: list[dict[str, str]],
    *,
    conflict: tuple[int, str, str],
    winner: dict[str, str] | None,
) -> None:
    """Record every candidate participating in one taxonomy-name conflict."""
    rank_index, name, policy = conflict
    rank = TAXONOMY_RANKS[rank_index]
    winner_prefix = lineage_prefix(winner, rank_index) if winner else None
    for row in candidates:
        decision = (
            "ambiguous"
            if winner is None
            else "winner"
            if lineage_prefix(row, rank_index) == winner_prefix
            else "rewritten"
        )
        diagnostics.append(
            {
                "conflict_type": "taxonomy_name",
                "rank": rank,
                "normalized_name": name,
                "policy": policy,
                "decision": decision,
                "source": row["source"],
                "ident": row["ident"],
                "signature_identity": "",
                "specificity": ":".join(
                    str(value) for value in specificity(row, rank_index)
                ),
                "lineage": lineage_text(row),
            },
        )


def curate_taxonomy_names(
    rows: list[dict[str, str]],
    *,
    policy: str,
    diagnostics: list[dict[str, str]],
) -> tuple[int, list[str]]:
    """Canonicalize conflicting taxon placements from broad to specific ranks."""
    resolved = 0
    unresolved: list[str] = []
    for rank_index, rank in enumerate(TAXONOMY_RANKS[1:], start=1):
        for name, candidates in taxonomy_conflicts_at_rank(rows, rank_index):
            winner = select_canonical_lineage(
                candidates,
                rank_index=rank_index,
                policy=policy,
            )
            if winner is None:
                unresolved.append(f"{rank}={name}")

            append_taxonomy_diagnostics(
                diagnostics,
                candidates,
                conflict=(rank_index, name, policy),
                winner=winner,
            )

            if winner is None:
                continue

            resolved += 1
            for row in candidates:
                for ancestor_rank in TAXONOMY_RANKS[: rank_index + 1]:
                    row[ancestor_rank] = winner[ancestor_rank]

    return resolved, unresolved


def find_cross_rank_taxon_names(
    rows: list[dict[str, str]],
    *,
    policy: str,
    diagnostics: list[dict[str, str]],
) -> list[str]:
    """Record names reused at multiple ranks, which Taxburst cannot represent."""
    occurrences: defaultdict[str, list[tuple[str, dict[str, str]]]] = defaultdict(list)
    for row in rows:
        for rank in TAXONOMY_RANKS:
            name = normalized_taxon_name(row[rank])
            if name:
                occurrences[name].append((rank, row))

    conflicts = []
    for name, candidates in occurrences.items():
        ranks = {rank for rank, _row in candidates}
        if len(ranks) <= 1:
            continue
        conflicts.append(f"{name} ({', '.join(sorted(ranks))})")
        for rank, row in candidates:
            diagnostics.append(
                {
                    "conflict_type": "taxonomy_rank",
                    "rank": rank,
                    "normalized_name": name,
                    "policy": policy,
                    "decision": "ambiguous",
                    "source": row["source"],
                    "ident": row["ident"],
                    "signature_identity": "",
                    "specificity": "",
                    "lineage": lineage_text(row),
                },
            )
    return conflicts


def write_curated_lineages(path: Path, rows: list[dict[str, str]]) -> None:
    """Write only the lineage columns consumed by sourmash taxonomy commands."""
    ensure_parent(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CURATED_LINEAGE_COLUMNS,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {column: row[column] for column in CURATED_LINEAGE_COLUMNS} for row in rows
        )


def curate_lineages(
    paths: CurationPaths,
    *,
    policy: str,
) -> None:
    """Validate sketch identities and canonicalize conflicting taxonomy names."""
    manifest_rows = [
        *read_manifest_dicts(paths.reference_manifest_csv, source="ncbi"),
        *read_manifest_dicts(paths.wvdb_manifest_csv, source="wvdb"),
    ]
    lineage_rows = [
        *read_lineage_dicts(paths.reference_lineages_csv, source="ncbi"),
        *read_lineage_dicts(paths.wvdb_lineages_csv, source="wvdb"),
    ]
    lineages = lineage_by_identifier(lineage_rows)
    missing = [
        record["ident"] for record in manifest_rows if record["ident"] not in lineages
    ]
    if missing:
        msg = f"{len(missing)} signatures are missing lineages. Examples: {', '.join(missing[:10])}"
        raise ValueError(msg)

    diagnostics: list[dict[str, str]] = []
    try:
        fail_on_duplicate_signature_identities(manifest_rows, lineages, diagnostics)
    except ValueError:
        write_diagnostics(paths.diagnostics_tsv, diagnostics)
        raise

    retained_rows = [lineages[record["ident"]] for record in manifest_rows]
    resolved, unresolved = curate_taxonomy_names(
        retained_rows,
        policy=policy,
        diagnostics=diagnostics,
    )
    cross_rank_conflicts = find_cross_rank_taxon_names(
        retained_rows,
        policy=policy,
        diagnostics=diagnostics,
    )
    write_diagnostics(paths.diagnostics_tsv, diagnostics)
    if resolved:
        warnings.warn(
            f"resolved {resolved} taxonomy conflict(s) using policy {policy}; "
            f"evidence={paths.diagnostics_tsv}",
            RuntimeWarning,
            stacklevel=2,
        )
    if unresolved or cross_rank_conflicts:
        examples = [*unresolved, *cross_rank_conflicts]
        msg = (
            f"{len(examples)} ambiguous taxonomy conflict(s) under policy {policy}. "
            f"Examples: {', '.join(examples[:10])}. Evidence: {paths.diagnostics_tsv}"
        )
        raise ValueError(msg)

    write_curated_lineages(paths.lineages_csv, retained_rows)


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

    curate = subcommands.add_parser(
        "curate-lineages",
        help="Validate signature identity and resolve conflicting taxonomy names.",
    )
    curate.add_argument("--reference-manifest-csv", type=Path, required=True)
    curate.add_argument("--wvdb-manifest-csv", type=Path, required=True)
    curate.add_argument("--reference-lineages-csv", type=Path, required=True)
    curate.add_argument("--wvdb-lineages-csv", type=Path, required=True)
    curate.add_argument(
        "--taxonomy-conflict-policy",
        choices=TAXONOMY_CONFLICT_POLICIES,
        required=True,
    )
    curate.add_argument("--lineages-csv", type=Path, required=True)
    curate.add_argument("--diagnostics-tsv", type=Path, required=True)

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
        elif args.command == "curate-lineages":
            curate_lineages(
                CurationPaths(
                    reference_manifest_csv=args.reference_manifest_csv,
                    wvdb_manifest_csv=args.wvdb_manifest_csv,
                    reference_lineages_csv=args.reference_lineages_csv,
                    wvdb_lineages_csv=args.wvdb_lineages_csv,
                    lineages_csv=args.lineages_csv,
                    diagnostics_tsv=args.diagnostics_tsv,
                ),
                policy=args.taxonomy_conflict_policy,
            )
    except ValueError as exc:
        parser.exit(1, f"error: {exc}\n")


if __name__ == "__main__":
    main()
