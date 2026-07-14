#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = [
#   "polars>=1.27.1",
# ]
# ///
"""Build the taxon-level NVD Big Table from query-level Big Tables."""

# ruff: noqa: COM812

from __future__ import annotations

import argparse
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence


REQUIRED_QUERY_COLUMNS = {
    "sample_id",
    "assigned_taxid",
    "assigned_taxid_name",
    "assigned_taxid_rank",
    "who_risk_group",
    "support_tier",
    "query_class",
    "crumbs_score",
    "qlen",
}
OUTPUT_COLUMNS = [
    "sample_id",
    "taxon_name",
    "taxon_rank",
    "support_tier",
    "taxon_crumbs",
    "relative_crumbs_percent",
    "supporting_query_count",
    "taxid",
    "who_risk_group",
    "total_query_span",
    "total_crumbs_score",
    "strong_query_count",
    "moderate_query_count",
    "weak_query_count",
    "review_query_count",
    "redacted_query_count",
    "supporting_genome_like_contig_count",
    "supporting_long_contig_count",
    "supporting_short_contig_count",
    "supporting_merged_pair_count",
    "supporting_single_read_count",
    "support_tier_rule",
    "support_note",
]
SUPPORT_TIER_ORDER = {
    "strong": 1,
    "moderate": 2,
    "review": 3,
    "weak": 4,
    "redacted": 5,
}
MULTI_QUERY_STRONG_MIN = 2


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a one-row-per-taxon Big Table from query Big Tables.",
    )
    parser.add_argument("--query-big-table", type=Path, required=True)
    parser.add_argument("--crumbs-taxa-tsv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
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


def scan_query_big_table(path: Path) -> pl.LazyFrame:
    return scan_tsv(path).pipe(
        require_columns,
        REQUIRED_QUERY_COLUMNS,
        label=str(path),
    )


def scan_crumbs_taxa(path: Path | None) -> pl.LazyFrame:
    schema = {
        "sample_id": pl.String,
        "taxid": pl.String,
        "taxon_crumbs": pl.String,
        "percentage_emitted": pl.String,
        "crumbs_total_crumbs_score": pl.String,
    }
    if path is None:
        return pl.LazyFrame(schema=schema)
    return (
        scan_tsv(path)
        .pipe(
            require_columns,
            {
                "sample_id",
                "taxon_id",
                "total_crumbs_score",
                "taxon_crumbs",
                "percentage_emitted",
            },
            label=str(path),
        )
        .select(
            "sample_id",
            pl.col("taxon_id").alias("taxid"),
            pl.col("taxon_crumbs"),
            pl.col("percentage_emitted"),
            pl.col("total_crumbs_score").alias("crumbs_total_crumbs_score"),
        )
        .group_by("sample_id", "taxid")
        .agg(
            *[
                pl.when(pl.col(column).n_unique() == 1)
                .then(pl.col(column).first())
                .otherwise(None)
                .alias(column)
                for column in (
                    "taxon_crumbs",
                    "percentage_emitted",
                    "crumbs_total_crumbs_score",
                )
            ],
        )
    )


def with_numeric_query_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.with_columns(
        pl.col("qlen").cast(pl.Int64, strict=False).alias("_qlen_num"),
        pl.col("crumbs_score")
        .cast(pl.Float64, strict=False)
        .alias("_crumbs_score_num"),
    )


def query_count(condition: pl.Expr, alias: str) -> pl.Expr:
    return condition.sum().alias(alias)


def aggregate_queries(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.group_by(
        "sample_id",
        pl.col("assigned_taxid").alias("taxid"),
        pl.col("assigned_taxid_name").alias("taxon_name"),
        pl.col("assigned_taxid_rank").alias("taxon_rank"),
    ).agg(
        pl.col("who_risk_group").first(),
        pl.len().alias("supporting_query_count"),
        pl.col("_qlen_num").sum().alias("total_query_span"),
        pl.col("_crumbs_score_num").sum().alias("query_total_crumbs_score"),
        query_count(pl.col("support_tier") == "strong", "strong_query_count"),
        query_count(pl.col("support_tier") == "moderate", "moderate_query_count"),
        query_count(pl.col("support_tier") == "weak", "weak_query_count"),
        query_count(pl.col("support_tier") == "review", "review_query_count"),
        query_count(pl.col("support_tier") == "redacted", "redacted_query_count"),
        query_count(
            pl.col("query_class") == "genome_like_assembly_contig",
            "supporting_genome_like_contig_count",
        ),
        query_count(
            pl.col("query_class") == "long_assembly_contig",
            "supporting_long_contig_count",
        ),
        query_count(
            pl.col("query_class") == "short_assembly_contig",
            "supporting_short_contig_count",
        ),
        query_count(
            pl.col("query_class") == "overlap_merged_pair",
            "supporting_merged_pair_count",
        ),
        query_count(
            pl.col("query_class") == "single_read", "supporting_single_read_count"
        ),
        query_count(
            (pl.col("support_tier") == "strong")
            & (pl.col("query_class") == "genome_like_assembly_contig"),
            "_strong_genome_like_query_count",
        ),
        query_count(
            (pl.col("support_tier") == "strong")
            & (pl.col("query_class") == "long_assembly_contig"),
            "_strong_long_query_count",
        ),
        query_count(
            (pl.col("support_tier") == "strong")
            & (pl.col("query_class") == "short_assembly_contig"),
            "_strong_short_query_count",
        ),
        query_count(
            (pl.col("support_tier") == "strong")
            & pl.col("query_class").is_in(["overlap_merged_pair", "single_read"]),
            "_strong_read_query_count",
        ),
    )


def with_derived_crumbs(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.with_columns(
        pl.when(pl.col("total_query_span") > 0)
        .then(pl.col("query_total_crumbs_score") / pl.col("total_query_span"))
        .otherwise(None)
        .alias("query_taxon_crumbs"),
    )


def with_taxon_crumbs(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.with_columns(
        pl.coalesce("crumbs_total_crumbs_score", "query_total_crumbs_score").alias(
            "total_crumbs_score",
        ),
        pl.coalesce("taxon_crumbs", "query_taxon_crumbs").alias("taxon_crumbs"),
        pl.coalesce("percentage_emitted", pl.lit("")).alias("relative_crumbs_percent"),
    )


def taxon_support_branch() -> pl.Expr:
    single_query = pl.col("supporting_query_count") == 1
    strong_queries = pl.col("strong_query_count")
    strong_long_or_genome_like = pl.col("_strong_genome_like_query_count") + pl.col(
        "_strong_long_query_count"
    )
    return (
        pl.when(pl.col("redacted_query_count") > 0)
        .then(pl.lit("redacted"))
        .when(strong_queries >= MULTI_QUERY_STRONG_MIN)
        .then(pl.lit("multi_query_strong_support"))
        .when(strong_long_or_genome_like > 0)
        .then(pl.lit("strong_long_or_genome_like_query"))
        .when(single_query & (pl.col("_strong_short_query_count") == 1))
        .then(pl.lit("single_short_contig_strong_query_lacks_corroboration"))
        .when(single_query & (pl.col("_strong_read_query_count") == 1))
        .then(pl.lit("single_read_strong_query_lacks_corroboration"))
        .when(strong_queries > 0)
        .then(pl.lit("mixed_support_with_strong_query"))
        .when(pl.col("moderate_query_count") > 0)
        .then(pl.lit("moderate_query_support"))
        .when(pl.col("review_query_count") > 0)
        .then(pl.lit("review_query_support"))
        .when(pl.col("weak_query_count") == pl.col("supporting_query_count"))
        .then(pl.lit("all_weak_queries"))
        .otherwise(pl.lit("mixed_query_support"))
    )


def taxon_label() -> pl.Expr:
    return (
        pl.when(pl.col("taxon_name").is_not_null() & pl.col("taxon_rank").is_not_null())
        .then(
            pl.concat_str(
                [
                    pl.col("taxon_name"),
                    pl.lit(" ("),
                    pl.col("taxon_rank"),
                    pl.lit(")"),
                ],
            ),
        )
        .otherwise(
            pl.concat_str(
                [pl.lit("taxid "), pl.col("taxid"), pl.lit(" (rank unavailable)")],
            ),
        )
    )


def taxon_support_note() -> pl.Expr:
    taxon = taxon_label()
    additional_queries = pl.col("supporting_query_count") - 1
    long_query_kind = (
        pl.when(pl.col("_strong_genome_like_query_count") > 0)
        .then(pl.lit("genome-like-contig"))
        .otherwise(pl.lit("long-contig"))
    )
    strong_long_note = pl.concat_str(
        [
            pl.lit("1 strong "),
            long_query_kind,
            pl.lit(" query assignment supports "),
            taxon,
            pl.when(additional_queries == 0)
            .then(pl.lit("."))
            .otherwise(
                pl.concat_str(
                    [
                        pl.lit("; "),
                        additional_queries.cast(pl.String),
                        pl.when(additional_queries == 1)
                        .then(pl.lit(" additional query is"))
                        .otherwise(pl.lit(" additional queries are")),
                        pl.lit(" assigned to this taxon."),
                    ],
                ),
            ),
        ],
    )

    return (
        pl.when(pl.col("support_tier_rule") == "redacted")
        .then(
            pl.concat_str(
                [
                    pl.lit("At least one query assignment to "),
                    taxon,
                    pl.lit(" is redacted."),
                ]
            )
        )
        .when(pl.col("support_tier_rule") == "multi_query_strong_support")
        .then(
            pl.concat_str(
                [
                    pl.col("strong_query_count").cast(pl.String),
                    pl.lit(" strong query assignments support "),
                    taxon,
                    pl.lit("."),
                ],
            ),
        )
        .when(pl.col("support_tier_rule") == "strong_long_or_genome_like_query")
        .then(strong_long_note)
        .when(
            pl.col("support_tier_rule")
            == "single_short_contig_strong_query_lacks_corroboration"
        )
        .then(
            pl.concat_str(
                [
                    pl.lit("1 strong short-contig query assignment supports "),
                    taxon,
                    pl.lit("; no additional query is assigned to this taxon."),
                ],
            ),
        )
        .when(
            pl.col("support_tier_rule")
            == "single_read_strong_query_lacks_corroboration"
        )
        .then(
            pl.concat_str(
                [
                    pl.lit("1 strong read-derived query assignment supports "),
                    taxon,
                    pl.lit("; no additional query is assigned to this taxon."),
                ],
            ),
        )
        .when(pl.col("support_tier_rule") == "mixed_support_with_strong_query")
        .then(
            pl.concat_str(
                [
                    pl.col("strong_query_count").cast(pl.String),
                    pl.lit(" strong and "),
                    (
                        pl.col("supporting_query_count") - pl.col("strong_query_count")
                    ).cast(pl.String),
                    pl.when(
                        (
                            pl.col("supporting_query_count")
                            - pl.col("strong_query_count")
                        )
                        == 1
                    )
                    .then(pl.lit(" additional query assignment support "))
                    .otherwise(pl.lit(" additional query assignments support ")),
                    taxon,
                    pl.lit("."),
                ],
            ),
        )
        .when(pl.col("support_tier_rule") == "moderate_query_support")
        .then(
            pl.concat_str(
                [
                    pl.col("supporting_query_count").cast(pl.String),
                    pl.lit(" query assignments support "),
                    taxon,
                    pl.lit(", including "),
                    pl.col("moderate_query_count").cast(pl.String),
                    pl.lit(" moderate assignments."),
                ],
            ),
        )
        .when(pl.col("support_tier_rule") == "review_query_support")
        .then(
            pl.concat_str(
                [
                    pl.col("review_query_count").cast(pl.String),
                    pl.lit(" review query assignments support "),
                    taxon,
                    pl.lit("; none is strong or moderate."),
                ],
            ),
        )
        .when(pl.col("support_tier_rule") == "all_weak_queries")
        .then(
            pl.concat_str(
                [
                    pl.lit("All "),
                    pl.col("supporting_query_count").cast(pl.String),
                    pl.lit(" query assignments to "),
                    taxon,
                    pl.lit(" have low best-hit query coverage."),
                ],
            ),
        )
        .otherwise(
            pl.concat_str(
                [
                    pl.col("supporting_query_count").cast(pl.String),
                    pl.lit(" query assignments support "),
                    taxon,
                    pl.lit("."),
                ],
            ),
        )
    )


def with_support_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    branch = taxon_support_branch()
    return frame.with_columns(branch.alias("support_tier_rule")).with_columns(
        pl.when(pl.col("support_tier_rule") == "redacted")
        .then(pl.lit("redacted"))
        .when(
            pl.col("support_tier_rule").is_in(
                ["multi_query_strong_support", "strong_long_or_genome_like_query"]
            )
        )
        .then(pl.lit("strong"))
        .when(pl.col("support_tier_rule") == "all_weak_queries")
        .then(pl.lit("weak"))
        .when(pl.col("support_tier_rule") == "review_query_support")
        .then(pl.lit("review"))
        .otherwise(pl.lit("moderate"))
        .alias("support_tier"),
        taxon_support_note().alias("support_note"),
    )


def select_output_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.select(OUTPUT_COLUMNS)


def with_salience_sort_columns(frame: pl.LazyFrame) -> pl.LazyFrame:
    return frame.with_columns(
        pl.col("support_tier")
        .replace_strict(SUPPORT_TIER_ORDER, default=99)
        .alias("_support_tier_sort"),
        pl.col("taxon_crumbs")
        .cast(pl.Float64, strict=False)
        .alias("_taxon_crumbs_num"),
        pl.col("total_crumbs_score")
        .cast(pl.Float64, strict=False)
        .alias("_total_crumbs_score_num"),
    )


def salience_sort(frame: pl.LazyFrame) -> pl.LazyFrame:
    return (
        frame.pipe(with_salience_sort_columns)
        .sort(
            [
                "sample_id",
                "_support_tier_sort",
                "_taxon_crumbs_num",
                "_total_crumbs_score_num",
                "taxon_name",
                "taxid",
            ],
            descending=[False, False, True, True, False, False],
            nulls_last=True,
        )
        .drop("_support_tier_sort", "_taxon_crumbs_num", "_total_crumbs_score_num")
    )


def build_taxon_big_table(
    query_big_tables: pl.LazyFrame,
    crumbs_taxa: pl.LazyFrame,
) -> pl.LazyFrame:
    return (
        query_big_tables.pipe(
            require_columns, REQUIRED_QUERY_COLUMNS, label="query Big Table"
        )
        .pipe(with_numeric_query_columns)
        .pipe(aggregate_queries)
        .pipe(with_derived_crumbs)
        .join(crumbs_taxa, on=["sample_id", "taxid"], how="left")
        .pipe(with_taxon_crumbs)
        .pipe(with_support_columns)
        .pipe(select_output_columns)
        .pipe(salience_sort)
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    build_taxon_big_table(
        scan_query_big_table(args.query_big_table),
        scan_crumbs_taxa(args.crumbs_taxa_tsv),
    ).collect().write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
