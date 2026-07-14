#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#   "polars>=1.27.1",
# ]
# ///
"""Build the query-level NVD Big Table from retained BLAST evidence."""

# ruff: noqa: COM812

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl
from py_nvd.blast_consensus import ASSIGNMENT_BITSCORE_FRACTION

if TYPE_CHECKING:
    from collections.abc import Sequence


REQUIRED_BLAST_COLUMNS = {
    "sample",
    "qseqid",
    "qlen",
    "sseqid",
    "stitle",
    "length",
    "pident",
    "evalue",
    "bitscore",
    "staxids",
    "rank",
    "adjusted_taxid",
    "adjusted_taxid_name",
    "adjusted_taxid_rank",
    "who_risk_group",
    "adjustment_method",
    "query_class",
    "producer",
    "source_id",
    "support_record_count",
    "mapped_reads",
    "blast_db_version",
    "virus_index_version",
    "nextflow_run_id",
}
OUTPUT_COLUMNS = [
    "sample_id",
    "assigned_taxid_name",
    "assigned_taxid_rank",
    "support_tier",
    "query_class",
    "crumbs_score",
    "qlen",
    "assigned_taxid",
    "who_risk_group",
    "assignment_method",
    "qseqid",
    "best_hit_qcov",
    "best_hit_pident",
    "best_hit_evalue",
    "best_hit_bitscore",
    "retained_reference_count",
    "assignment_reference_count",
    "assignment_taxid_count",
    "support_tier_rule",
    "support_note",
    "support_record_count",
    "mapped_reads",
    "producer",
    "source_id",
    "best_hit_alignment_length",
    "blast_db_version",
    "virus_index_version",
    "nextflow_run_id",
]
SUPPORT_TIER_ORDER = {
    "strong": 1,
    "moderate": 2,
    "review": 3,
    "weak": 4,
    "redacted": 5,
}
QUERY_CLASS_ORDER = {
    "genome_like_assembly_contig": 1,
    "long_assembly_contig": 2,
    "short_assembly_contig": 3,
    "overlap_merged_pair": 4,
    "single_read": 5,
}


@dataclass(frozen=True)
class Thresholds:
    short_qcov: float = 0.4
    long_qcov: float = 0.8


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a one-row-per-query Big Table from final BLAST rows.",
    )
    parser.add_argument("--blast-tsv", type=Path, required=True)
    parser.add_argument("--crumbs-tsv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--short-qcov", type=float, default=Thresholds.short_qcov)
    parser.add_argument("--long-qcov", type=float, default=Thresholds.long_qcov)
    return parser.parse_args(argv)


def scan_tsv(path: Path) -> pl.LazyFrame:
    return pl.scan_csv(path, separator="\t", infer_schema=False)


def require_columns(
    frame: pl.LazyFrame, required: set[str], *, label: str
) -> pl.LazyFrame:
    missing = sorted(required - set(frame.collect_schema().names()))
    if missing:
        message = f"{label} missing required columns: {', '.join(missing)}"
        raise ValueError(message)
    return frame


def scan_crumbs(path: Path | None) -> pl.LazyFrame:
    if path is None:
        return pl.LazyFrame(
            schema={
                "sample_id": pl.String,
                "qseqid": pl.String,
                "crumbs_score": pl.String,
            },
        )
    return (
        scan_tsv(path)
        .pipe(require_columns, {"sample_id", "qseqid", "crumbs_score"}, label=str(path))
        .select("sample_id", "qseqid", "crumbs_score")
        .group_by("sample_id", "qseqid")
        .agg(
            pl.when(pl.col("crumbs_score").n_unique() == 1)
            .then(pl.col("crumbs_score").first())
            .otherwise(None)
            .alias("crumbs_score"),
        )
    )


def with_numeric_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    sample_expr = (
        pl.coalesce(pl.col("sample_id"), pl.col("sample"))
        if "sample_id" in frame.collect_schema().names()
        else pl.col("sample")
    )
    return frame.with_columns(
        sample_expr.alias("sample_id"),
        pl.col("qlen").cast(pl.Float64).alias("_qlen_num"),
        pl.col("length").cast(pl.Float64).alias("_length_num"),
        pl.col("pident").cast(pl.Float64).alias("_pident_num"),
        pl.col("evalue").cast(pl.Float64).alias("_evalue_num"),
        pl.col("bitscore").cast(pl.Float64).alias("_bitscore_num"),
        pl.when(pl.col("staxids").is_null() | (pl.col("staxids") == ""))
        .then(None)
        .otherwise(pl.col("staxids"))
        .alias("_canonical_taxid"),
    )


def with_score_context(frame: pl.LazyFrame) -> pl.LazyFrame:
    return (
        frame.sort(
            [
                "sample_id",
                "qseqid",
                "_bitscore_num",
                "_evalue_num",
                "_pident_num",
                "_length_num",
            ],
            descending=[False, False, True, False, True, True],
        )
        .with_columns(
            pl.max("_bitscore_num").over("sample_id", "qseqid").alias("_best_bitscore"),
        )
        .with_columns(
            (
                pl.col("_bitscore_num")
                >= pl.col("_best_bitscore") * ASSIGNMENT_BITSCORE_FRACTION
            ).alias("_used_for_assignment"),
        )
    )


def group_metrics(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.group_by("sample_id", "qseqid").agg(
        pl.col("sseqid").n_unique().alias("retained_reference_count"),
        pl.col("sseqid")
        .filter(pl.col("_used_for_assignment"))
        .n_unique()
        .alias("assignment_reference_count"),
        pl.col("_canonical_taxid")
        .filter(pl.col("_used_for_assignment"))
        .drop_nulls()
        .n_unique()
        .alias("assignment_taxid_count"),
    )


def with_query_class_family(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.with_columns(
        pl.when(pl.col("query_class") == "long_assembly_contig")
        .then(pl.lit("long_contig"))
        .when(pl.col("query_class") == "short_assembly_contig")
        .then(pl.lit("short_contig"))
        .when(pl.col("query_class") == "overlap_merged_pair")
        .then(pl.lit("merged_pair"))
        .when(pl.col("query_class") == "single_read")
        .then(pl.lit("single_read"))
        .when(pl.col("query_class") == "genome_like_assembly_contig")
        .then(pl.lit("genome_like_contig"))
        .otherwise(pl.lit("query"))
        .alias("_query_class_family"),
    )


def coverage_band(thresholds: Thresholds) -> pl.Expr:
    return (
        pl.when(pl.col("best_hit_qcov") >= thresholds.long_qcov)
        .then(pl.lit("high_qcov"))
        .when(pl.col("best_hit_qcov") < thresholds.short_qcov)
        .then(pl.lit("low_qcov"))
        .otherwise(pl.lit("partial_qcov"))
    )


def support_note() -> pl.Expr:
    alignment = pl.concat_str(
        [
            pl.lit("The best alignment spans "),
            pl.col("length"),
            pl.lit(" of "),
            pl.col("qlen"),
            pl.lit(" query bases ("),
            (pl.col("best_hit_qcov") * 100).round(1).cast(pl.String),
            pl.lit("%). "),
        ],
    )
    retained_references = pl.concat_str(
        [
            pl.col("retained_reference_count").cast(pl.String),
            pl.when(pl.col("retained_reference_count") == 1)
            .then(pl.lit(" retained reference"))
            .otherwise(pl.lit(" retained references")),
        ],
    )
    assignment_references = (
        pl.when(pl.col("assignment_reference_count") == 1)
        .then(pl.lit("1 has a bitscore"))
        .otherwise(
            pl.concat_str(
                [
                    pl.col("assignment_reference_count").cast(pl.String),
                    pl.lit(" have bitscores"),
                ],
            ),
        )
    )
    assignment_taxids = pl.concat_str(
        [
            pl.col("assignment_taxid_count").cast(pl.String),
            pl.when(pl.col("assignment_taxid_count") == 1)
            .then(pl.lit(" taxid"))
            .otherwise(pl.lit(" taxids")),
        ],
    )
    score_context = pl.concat_str(
        [
            alignment,
            pl.lit("Of "),
            retained_references,
            pl.lit(", "),
            assignment_references,
            pl.lit(
                f" at least {ASSIGNMENT_BITSCORE_FRACTION:.0%} of the best bitscore"
            ),
        ],
    )
    assigned_taxon = pl.concat_str(
        [
            pl.col("adjusted_taxid_name"),
            pl.lit(" ("),
            pl.col("adjusted_taxid_rank"),
            pl.lit(")"),
        ],
    )
    metadata_available = (
        pl.col("adjusted_taxid_name").is_not_null()
        & pl.col("adjusted_taxid_rank").is_not_null()
        & pl.col("adjustment_method").is_in(["dominant", "lca"])
    )

    return (
        pl.when(~metadata_available)
        .then(
            pl.concat_str(
                [
                    alignment,
                    pl.lit(
                        "Taxonomic assignment metadata is unavailable for the selected taxid."
                    ),
                ],
            ),
        )
        .when(pl.col("adjustment_method") == "dominant")
        .then(
            pl.concat_str(
                [
                    score_context,
                    pl.lit("; the available taxid from those references resolves to "),
                    assigned_taxon,
                    pl.lit("."),
                ],
            ),
        )
        .otherwise(
            pl.concat_str(
                [
                    score_context,
                    pl.lit(" and map to "),
                    assignment_taxids,
                    pl.lit("; their LCA is "),
                    assigned_taxon,
                    pl.lit("."),
                ],
            ),
        )
    )


def with_support_columns(frame: pl.LazyFrame, thresholds: Thresholds) -> pl.LazyFrame:
    return (
        frame.pipe(with_query_class_family)
        .with_columns(coverage_band(thresholds).alias("_coverage_band"))
        .with_columns(
            pl.when(pl.col("_coverage_band") == "high_qcov")
            .then(pl.lit("strong"))
            .when(pl.col("_coverage_band") == "low_qcov")
            .then(pl.lit("weak"))
            .otherwise(pl.lit("moderate"))
            .alias("support_tier"),
            pl.concat_str(
                [
                    pl.col("_query_class_family"),
                    pl.lit("_"),
                    pl.col("adjustment_method"),
                    pl.lit("_"),
                    pl.col("_coverage_band"),
                ],
            ).alias("support_tier_rule"),
            support_note().alias("support_note"),
        )
    )


def assignment_rows(frame: pl.LazyFrame, thresholds: Thresholds) -> pl.LazyFrame:
    metrics = frame.pipe(group_metrics)
    return (
        frame.unique(
            subset=["sample_id", "qseqid"],
            keep="first",
            maintain_order=True,
        )
        .join(metrics, on=["sample_id", "qseqid"], how="inner")
        .with_columns(
            (pl.col("_length_num") / pl.col("_qlen_num")).alias("best_hit_qcov"),
        )
        .pipe(with_support_columns, thresholds)
    )


def select_output_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.select(
        "sample_id",
        pl.col("adjusted_taxid_name").alias("assigned_taxid_name"),
        pl.col("adjusted_taxid_rank").alias("assigned_taxid_rank"),
        pl.col("support_tier"),
        pl.col("query_class"),
        pl.col("crumbs_score"),
        pl.col("qlen"),
        pl.col("adjusted_taxid").alias("assigned_taxid"),
        pl.col("who_risk_group"),
        pl.col("adjustment_method").alias("assignment_method"),
        "qseqid",
        pl.col("best_hit_qcov"),
        pl.col("pident").alias("best_hit_pident"),
        pl.col("evalue").alias("best_hit_evalue"),
        pl.col("bitscore").alias("best_hit_bitscore"),
        pl.col("retained_reference_count"),
        pl.col("assignment_reference_count"),
        pl.col("assignment_taxid_count"),
        pl.col("support_tier_rule"),
        pl.col("support_note"),
        pl.col("support_record_count"),
        pl.col("mapped_reads"),
        pl.col("producer"),
        pl.col("source_id"),
        pl.col("length").alias("best_hit_alignment_length"),
        pl.col("blast_db_version"),
        pl.col("virus_index_version"),
        pl.col("nextflow_run_id"),
    )


def with_salience_sort_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.with_columns(
        pl.col("support_tier")
        .replace_strict(SUPPORT_TIER_ORDER, default=99)
        .alias("_support_tier_sort"),
        pl.col("query_class")
        .replace_strict(QUERY_CLASS_ORDER, default=99)
        .alias("_query_class_sort"),
        pl.col("crumbs_score")
        .cast(pl.Float64, strict=False)
        .alias("_crumbs_score_num"),
        pl.col("best_hit_evalue")
        .cast(pl.Float64, strict=False)
        .alias("_best_hit_evalue_num"),
    )


def salience_sort(frame: pl.LazyFrame) -> pl.LazyFrame:
    return (
        frame.pipe(with_salience_sort_columns)
        .sort(
            [
                "sample_id",
                "_support_tier_sort",
                "_crumbs_score_num",
                "_query_class_sort",
                "assigned_taxid_name",
                "_best_hit_evalue_num",
                "qseqid",
            ],
            descending=[False, False, True, False, False, False, False],
            nulls_last=True,
        )
        .drop(
            "_support_tier_sort",
            "_query_class_sort",
            "_crumbs_score_num",
            "_best_hit_evalue_num",
        )
    )


def build_query_big_table(
    blast_rows: pl.LazyFrame,
    crumbs: pl.LazyFrame,
    thresholds: Thresholds,
) -> pl.LazyFrame:
    return (
        blast_rows.pipe(require_columns, REQUIRED_BLAST_COLUMNS, label="BLAST TSV")
        .pipe(with_numeric_columns)
        .pipe(with_score_context)
        .pipe(assignment_rows, thresholds)
        .join(crumbs, on=["sample_id", "qseqid"], how="left")
        .pipe(select_output_columns)
        .pipe(salience_sort)
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    thresholds = Thresholds(
        short_qcov=args.short_qcov,
        long_qcov=args.long_qcov,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    build_query_big_table(
        scan_tsv(args.blast_tsv),
        scan_crumbs(args.crumbs_tsv),
        thresholds,
    ).collect().write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
