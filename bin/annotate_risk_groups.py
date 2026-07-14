#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["loguru", "polars", "taxopy"]
# ///
"""Annotate Sourmash and BLAST results with directional WHO risk groups.

The nearest classified taxon at or above a result determines its risk group.
Classifications may flow from ancestors to descendants, but never from
descendants to ancestors. Conflicts on one taxid resolve to the highest known
risk group and emit a warning.
"""

from __future__ import annotations

import argparse
import sqlite3
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from loguru import logger
from py_nvd import taxonomy

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from py_nvd.taxonomy import TaxonomyDB

RISK_GROUP_COLUMN = "who_risk_group"
JOIN_KEY_COLUMN = "_risk_group_join_key"
RISK_GROUPS = ("RG1", "RG2", "RG3", "RG4")
RISK_GROUP_PRIORITY = {
    risk_group: priority for priority, risk_group in enumerate(RISK_GROUPS, 1)
}
RiskGroupsByTaxid = dict[int, set[str]]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate taxonomic result tables with WHO risk groups.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    sourmash = commands.add_parser(
        "sourmash-summary",
        help="Annotate a sourmash tax metagenome CSV summary.",
    )
    sourmash.add_argument("--input", type=Path, required=True)
    sourmash.add_argument("--bioboxes", type=Path, required=True)
    sourmash.add_argument("--lookup", type=Path, required=True)
    sourmash.add_argument("--output", type=Path, required=True)
    blast = commands.add_parser(
        "blast",
        help="Annotate BLAST results using the consensus taxid.",
    )
    blast.add_argument("--input", type=Path, required=True)
    blast.add_argument("--lookup", type=Path, required=True)
    blast.add_argument("--taxonomy-dir", type=Path, required=True)
    blast.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def require_columns(
    frame: pl.LazyFrame,
    required: set[str],
    *,
    label: str,
) -> pl.LazyFrame:
    missing = sorted(required - set(frame.collect_schema().names()))
    if missing:
        message = f"{label} missing required columns: {', '.join(missing)}"
        raise ValueError(message)
    return frame


def require_new_risk_group_column(frame: pl.LazyFrame, *, label: str) -> pl.LazyFrame:
    reserved = {RISK_GROUP_COLUMN, JOIN_KEY_COLUMN}
    collisions = sorted(reserved.intersection(frame.collect_schema().names()))
    if collisions:
        message = f"{label} already contains reserved columns: {', '.join(collisions)}"
        raise ValueError(message)
    return frame


def load_risk_groups(path: Path) -> RiskGroupsByTaxid:
    rows = (
        pl.scan_csv(path, separator="\t", infer_schema=False)
        .pipe(
            require_columns,
            {"canonical_taxid", "human_classification"},
            label=str(path),
        )
        .select(
            pl.col("canonical_taxid")
            .str.strip_chars()
            .replace("", None)
            .alias("_taxid_text"),
            pl.col("human_classification")
            .str.strip_chars()
            .str.to_uppercase()
            .replace("", None)
            .alias(RISK_GROUP_COLUMN),
        )
        .with_columns(
            pl.col("_taxid_text").cast(pl.Int64, strict=False).alias("canonical_taxid"),
        )
        .collect()
    )
    malformed_taxids = rows.filter(
        pl.col("_taxid_text").is_not_null() & pl.col("canonical_taxid").is_null(),
    )
    if malformed_taxids.height:
        logger.warning(
            "Ignoring malformed canonical_taxid values in {}: {}",
            path,
            malformed_taxids.get_column("_taxid_text").unique().to_list(),
        )

    by_taxid: RiskGroupsByTaxid = {}
    invalid_groups: set[str] = set()
    nonpositive_taxids: set[int] = set()
    for row in rows.drop_nulls([RISK_GROUP_COLUMN]).iter_rows(named=True):
        group = row[RISK_GROUP_COLUMN]
        if group not in RISK_GROUP_PRIORITY:
            invalid_groups.add(group)
            continue
        tax_id = row["canonical_taxid"]
        if tax_id is None:
            continue
        if tax_id <= 0:
            nonpositive_taxids.add(tax_id)
            continue
        by_taxid.setdefault(row["canonical_taxid"], set()).add(
            group,
        )
    if invalid_groups:
        logger.warning(
            "Ignoring unsupported WHO risk-group values in {}: {}",
            path,
            sorted(invalid_groups),
        )
    if nonpositive_taxids:
        logger.warning(
            "Ignoring non-positive canonical_taxid values in {}: {}",
            path,
            sorted(nonpositive_taxids),
        )
    return by_taxid


def highest_risk_group(groups: set[str], *, tax_id: int) -> str | None:
    valid_groups = groups.intersection(RISK_GROUP_PRIORITY)
    if not valid_groups:
        return None
    if len(valid_groups) == 1:
        return next(iter(valid_groups))
    selected = max(valid_groups, key=RISK_GROUP_PRIORITY.__getitem__)
    logger.warning(
        "Conflicting WHO risk groups for taxid {}: {}; selecting {}",
        tax_id,
        sorted(valid_groups),
        selected,
    )
    return selected


def resolve_taxid_path(
    taxid_path: Sequence[int],
    classifications: Mapping[int, set[str]],
) -> str | None:
    for tax_id in reversed(taxid_path):
        groups = classifications.get(tax_id)
        if groups:
            return highest_risk_group(groups, tax_id=tax_id)
    return None


def load_bioboxes_paths(path: Path) -> dict[str, set[tuple[int, ...]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header_index = next(
        (index for index, line in enumerate(lines) if line.startswith("@@TAXID\t")),
        None,
    )
    if header_index is None:
        message = f"{path} missing BioBoxes @@TAXID header"
        raise ValueError(message)
    rows = pl.read_csv(
        StringIO("\n".join(lines[header_index:])),
        separator="\t",
        infer_schema=False,
    )
    required = {"@@TAXID", "TAXPATH", "TAXPATHSN"}
    missing = sorted(required - set(rows.columns))
    if missing:
        message = f"{path} missing required columns: {', '.join(missing)}"
        raise ValueError(message)

    paths: dict[str, set[tuple[int, ...]]] = {}
    for row in (
        rows.select("@@TAXID", "TAXPATH", "TAXPATHSN")
        .drop_nulls()
        .iter_rows(named=True)
    ):
        lineage = row["TAXPATHSN"].replace("|", ";")
        raw_taxids = row["TAXPATH"].split("|")
        raw_names = row["TAXPATHSN"].split("|")
        try:
            target_taxid = int(row["@@TAXID"])
            taxid_path = tuple(int(value) for value in raw_taxids)
        except ValueError:
            logger.warning(
                "Ignoring malformed BioBoxes taxids for lineage {!r}: taxid={!r}, "
                "taxpath={!r}",
                lineage,
                row["@@TAXID"],
                row["TAXPATH"],
            )
            continue
        if (
            target_taxid <= 0
            or any(tax_id <= 0 for tax_id in taxid_path)
            or any(not value.strip() for value in (*raw_taxids, *raw_names))
            or len(raw_taxids) != len(raw_names)
            or taxid_path[-1] != target_taxid
        ):
            logger.warning(
                "Ignoring inconsistent BioBoxes path for lineage {!r}: "
                "taxid={!r}, taxpath={!r}, taxpathsn={!r}",
                lineage,
                row["@@TAXID"],
                row["TAXPATH"],
                row["TAXPATHSN"],
            )
            continue
        paths.setdefault(lineage, set()).add(taxid_path)
    return paths


def resolve_sourmash_lineage(
    lineage: str,
    paths: Mapping[str, set[tuple[int, ...]]],
    classifications: Mapping[int, set[str]],
) -> str | None:
    taxid_paths = paths.get(lineage)
    if not taxid_paths:
        return None
    if len(taxid_paths) > 1:
        logger.warning(
            "Ambiguous BioBoxes taxpaths for lineage {!r}: {}",
            lineage,
            sorted(taxid_paths),
        )
        return None
    return resolve_taxid_path(next(iter(taxid_paths)), classifications)


def left_join_risk_groups(
    frame: pl.LazyFrame,
    lookup: pl.LazyFrame,
    *,
    join_key: pl.Expr,
) -> pl.LazyFrame:
    return (
        frame.with_columns(join_key.alias(JOIN_KEY_COLUMN))
        .join(lookup, on=JOIN_KEY_COLUMN, how="left", maintain_order="left")
        .drop(JOIN_KEY_COLUMN)
    )


def annotate_sourmash_summary(
    input_path: Path,
    bioboxes_path: Path,
    classifications: Mapping[int, set[str]],
) -> pl.LazyFrame:
    summary = (
        pl.scan_csv(input_path, infer_schema=False)
        .pipe(require_columns, {"lineage"}, label=str(input_path))
        .pipe(require_new_risk_group_column, label=str(input_path))
    )
    input_columns = summary.collect_schema().names()
    output_columns = input_columns.copy()
    output_columns.insert(output_columns.index("lineage") + 1, RISK_GROUP_COLUMN)
    lineages = (
        summary.select("lineage")
        .drop_nulls()
        .unique(maintain_order=True)
        .collect()
        .get_column("lineage")
        .to_list()
    )
    paths = load_bioboxes_paths(bioboxes_path)
    lookup = pl.DataFrame(
        {
            JOIN_KEY_COLUMN: lineages,
            RISK_GROUP_COLUMN: [
                resolve_sourmash_lineage(lineage, paths, classifications)
                for lineage in lineages
            ],
        },
        schema={JOIN_KEY_COLUMN: pl.String, RISK_GROUP_COLUMN: pl.String},
    ).lazy()
    return summary.pipe(
        left_join_risk_groups,
        lookup,
        join_key=pl.col("lineage"),
    ).select(output_columns)


def blast_taxid_path(tax_id: int, tax: TaxonomyDB | None) -> list[int]:
    if tax is None:
        return [tax_id]
    lineage = tax.get_lineage_ids(tax_id)
    return lineage or [tax_id]


def annotate_blast(
    input_path: Path,
    classifications: Mapping[int, set[str]],
    tax: TaxonomyDB | None,
) -> pl.LazyFrame:
    blast = (
        pl.scan_csv(input_path, separator="\t", infer_schema=False)
        .pipe(
            require_columns,
            {"adjusted_taxid", "adjusted_taxid_rank"},
            label=str(input_path),
        )
        .pipe(require_new_risk_group_column, label=str(input_path))
    )
    input_columns = blast.collect_schema().names()
    output_columns = input_columns.copy()
    output_columns.insert(
        output_columns.index("adjusted_taxid_rank") + 1,
        RISK_GROUP_COLUMN,
    )
    parsed_taxid_expr = (
        pl.col("adjusted_taxid").str.strip_chars().cast(pl.Int64, strict=False)
    )
    taxid_expr = pl.when(parsed_taxid_expr > 0).then(parsed_taxid_expr)
    invalid = (
        blast.filter(
            pl.col("adjusted_taxid").str.strip_chars().ne("") & taxid_expr.is_null(),
        )
        .select("adjusted_taxid")
        .unique()
        .collect()
    )
    if invalid.height:
        logger.warning(
            "Leaving WHO risk group blank for malformed adjusted_taxid values: {}",
            invalid.get_column("adjusted_taxid").to_list(),
        )
    taxids = (
        blast.select(taxid_expr.alias(JOIN_KEY_COLUMN))
        .drop_nulls()
        .unique(maintain_order=True)
        .collect()
        .get_column(JOIN_KEY_COLUMN)
        .to_list()
    )
    lookup = pl.DataFrame(
        {
            JOIN_KEY_COLUMN: taxids,
            RISK_GROUP_COLUMN: [
                resolve_taxid_path(blast_taxid_path(tax_id, tax), classifications)
                for tax_id in taxids
            ],
        },
        schema={JOIN_KEY_COLUMN: pl.Int64, RISK_GROUP_COLUMN: pl.String},
    ).lazy()
    return blast.pipe(
        left_join_risk_groups,
        lookup,
        join_key=taxid_expr,
    ).select(output_columns)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    classifications = load_risk_groups(args.lookup)
    if args.command == "sourmash-summary":
        annotate_sourmash_summary(
            args.input,
            args.bioboxes,
            classifications,
        ).sink_csv(args.output, null_value="")
        return

    try:
        with taxonomy.open(taxonomy_dir=args.taxonomy_dir, sync=True) as tax:
            annotated = annotate_blast(args.input, classifications, tax)
    except (
        taxonomy.ResourceUnavailableError,
        taxonomy.TaxonomyBuildError,
        sqlite3.DatabaseError,
    ) as error:
        logger.warning(
            "Taxonomy lineage traversal unavailable; using exact BLAST taxids only: {}",
            error,
        )
        annotated = annotate_blast(args.input, classifications, None)
    annotated.sink_csv(args.output, separator="\t", null_value="")


if __name__ == "__main__":
    main()
