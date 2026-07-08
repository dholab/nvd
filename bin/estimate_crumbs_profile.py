#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "polars",
# ]
# ///

"""Estimate CRUMBS taxonomic profile artifacts from NVD assignments."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence


ASSIGNMENT_COLUMNS = [
    "qseqid",
    "adjusted_taxid",
    "adjusted_taxid_name",
    "adjusted_taxid_rank",
    "adjustment_method",
]
PROFILE_ASSIGNMENT_SCORE_COLUMNS = ["bitscore", "evalue"]
READ_QUERY_CLASSES = {"overlap_merged_pair", "single_read"}
READ_QUERY_PREFIXES = ("nvdMergeReadQuery_", "nvdReadQuery_")
READ_QUERY_METADATA_COLUMNS = ["query_class", "support_record_count", "qlen"]
COVERAGE_COLUMNS = [
    "sample_id",
    "qseqid",
    "contig_length",
    "covered_bases_1x",
    "breadth_1x",
    "raw_aligned_bases",
    "mean_depth_full",
    "median_depth_full",
    "median_depth_positive",
    "depth_p95",
    "depth_p99",
    "max_depth",
    "crumbs_p95",
    "crumbs_p99",
]
PROFILE_TAXONOMY_COLUMNS = [
    "taxon_id",
    "taxon_name",
    "rank",
    "taxpath",
    "taxpathsn",
    "rankpath",
]
CONTIG_OUTPUT_COLUMNS = [
    "sample_id",
    "qseqid",
    "adjusted_taxid",
    "adjusted_taxid_name",
    "adjusted_taxid_rank",
    "adjustment_method",
    "contig_length",
    "covered_bases_1x",
    "breadth_1x",
    "raw_aligned_bases",
    "mean_depth_full",
    "median_depth_full",
    "median_depth_positive",
    "depth_p95",
    "depth_p99",
    "max_depth",
    "crumbs_p95",
    "crumbs_p99",
    "winsorization_percentile",
    "winsorization_cap",
    "crumbs_score",
]
TAXON_OUTPUT_COLUMNS = [
    "sample_id",
    "taxon_id",
    "taxon_name",
    "rank",
    "taxpath",
    "taxpathsn",
    "rankpath",
    "n_contigs",
    "total_contig_length",
    "total_covered_bases_1x",
    "total_raw_aligned_bases",
    "total_crumbs_score",
    "taxon_crumbs",
    "percentage_emitted",
    "n_zero_crumbs_contigs",
    "median_contig_breadth_1x",
    "min_contig_breadth_1x",
    "max_contig_breadth_1x",
    "fraction_crumbs_from_top_contig",
    "fraction_crumbs_from_low_breadth_contigs",
]
BIOBOXES_COLUMNS = ["@@TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"]
LOW_BREADTH_THRESHOLD = 0.1
FLOAT_TOLERANCE = 1e-9
BLAST_SCHEMA_OVERRIDES = {
    "sample": pl.String,
    "qseqid": pl.String,
    "evalue": pl.Float64,
    "bitscore": pl.Float64,
    "adjusted_taxid": pl.String,
    "adjusted_taxid_name": pl.String,
    "adjusted_taxid_rank": pl.String,
    "adjustment_method": pl.String,
    "query_class": pl.String,
    "support_record_count": pl.Int64,
    "qlen": pl.Int64,
}
COVERAGE_SCHEMA_OVERRIDES = {
    "sample_id": pl.String,
    "qseqid": pl.String,
    "contig_length": pl.Int64,
    "covered_bases_1x": pl.Int64,
    "breadth_1x": pl.Float64,
    "raw_aligned_bases": pl.Float64,
    "mean_depth_full": pl.Float64,
    "median_depth_full": pl.Float64,
    "median_depth_positive": pl.Float64,
    "depth_p95": pl.Float64,
    "depth_p99": pl.Float64,
    "max_depth": pl.Float64,
    "crumbs_p95": pl.Float64,
    "crumbs_p99": pl.Float64,
}
PROFILE_TAXONOMY_SCHEMA_OVERRIDES = {
    "taxon_id": pl.String,
    "taxon_name": pl.String,
    "rank": pl.String,
    "taxpath": pl.String,
    "taxpathsn": pl.String,
    "rankpath": pl.String,
}


class CrumbsProfileError(ValueError):
    """Raised when CRUMBS profile estimation cannot proceed safely."""


@dataclass(frozen=True)
class OutputPaths:
    contigs: Path
    taxa: Path
    bioboxes_profile: Path
    qc_json: Path


def output_paths(sample_id: str, output_dir: Path) -> OutputPaths:
    """Return standard CRUMBS output paths for a sample."""
    return OutputPaths(
        contigs=output_dir / f"{sample_id}.crumbs.contigs.tsv",
        taxa=output_dir / f"{sample_id}.crumbs.taxa.tsv",
        bioboxes_profile=output_dir / f"{sample_id}.crumbs.bioboxes.profile.tsv",
        qc_json=output_dir / f"{sample_id}.crumbs.qc.json",
    )


def require_columns(
    frame: pl.DataFrame,
    required: Sequence[str],
    *,
    label: str,
) -> None:
    """Fail when a required input column is absent."""
    missing = [column for column in required if column not in frame.columns]
    if missing:
        message = f"{label} missing required columns: {', '.join(missing)}"
        raise CrumbsProfileError(message)


def read_inputs(
    *,
    blast_tsv: Path,
    coverage_tsv: Path,
    profile_taxonomy_tsv: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Read and validate estimator inputs."""
    blast = pl.read_csv(
        blast_tsv,
        separator="\t",
        schema_overrides=BLAST_SCHEMA_OVERRIDES,
    )
    coverage = pl.read_csv(
        coverage_tsv,
        separator="\t",
        schema_overrides=COVERAGE_SCHEMA_OVERRIDES,
    )
    profile_taxonomy = pl.read_csv(
        profile_taxonomy_tsv,
        separator="\t",
        schema_overrides=PROFILE_TAXONOMY_SCHEMA_OVERRIDES,
    )

    require_columns(blast, ASSIGNMENT_COLUMNS, label="BLAST TSV")
    require_columns(coverage, COVERAGE_COLUMNS, label="coverage TSV")
    require_columns(
        profile_taxonomy,
        PROFILE_TAXONOMY_COLUMNS,
        label="profile taxonomy TSV",
    )
    return blast, coverage, profile_taxonomy


def collapse_assignments(blast: pl.DataFrame, sample_id: str) -> pl.DataFrame:
    """Collapse repeated BLAST hit rows to one assignment per contig."""
    if "sample" in blast.columns:
        blast = blast.filter(pl.col("sample").cast(pl.String) == sample_id)

    has_score_columns = all(
        column in blast.columns for column in PROFILE_ASSIGNMENT_SCORE_COLUMNS
    )
    selected_columns = [
        *ASSIGNMENT_COLUMNS,
        *(PROFILE_ASSIGNMENT_SCORE_COLUMNS if has_score_columns else []),
        *[column for column in READ_QUERY_METADATA_COLUMNS if column in blast.columns],
    ]
    assignments = blast.select(selected_columns).with_columns(
        *[pl.col(column).cast(pl.String) for column in ASSIGNMENT_COLUMNS],
    )

    if assignments.is_empty():
        message = f"no BLAST assignments found for sample {sample_id}"
        raise CrumbsProfileError(message)

    missing_required_values = assignments.filter(
        pl.any_horizontal(
            [
                pl.col(column).is_null() | (pl.col(column).str.strip_chars() == "")
                for column in ASSIGNMENT_COLUMNS
            ],
        ),
    )
    if missing_required_values.height > 0:
        qseqids = (
            missing_required_values.get_column("qseqid").fill_null("<null>").to_list()
        )
        message = (
            "BLAST assignments contain null or blank required values for contigs: "
            + ", ".join(qseqids)
        )
        raise CrumbsProfileError(message)

    if has_score_columns:
        missing_score_values = assignments.filter(
            pl.any_horizontal(
                [pl.col(column).is_null() for column in PROFILE_ASSIGNMENT_SCORE_COLUMNS],
            ),
        )
        if missing_score_values.height > 0:
            qseqids = missing_score_values.get_column("qseqid").to_list()
            message = (
                "BLAST assignments contain null score values for contigs: "
                + ", ".join(qseqids)
            )
            raise CrumbsProfileError(message)

        assignments = (
            assignments.with_columns(
                pl.max("bitscore").over("qseqid").alias("_best_bitscore"),
            )
            .filter(pl.col("bitscore") == pl.col("_best_bitscore"))
            .with_columns(
                pl.min("evalue").over("qseqid").alias("_best_evalue"),
            )
            .filter(pl.col("evalue") == pl.col("_best_evalue"))
            .select(
                [
                    *ASSIGNMENT_COLUMNS,
                    *[column for column in READ_QUERY_METADATA_COLUMNS if column in assignments.columns],
                ],
            )
        )

    conflicts = assignments.unique(maintain_order=True).group_by("qseqid").len()
    conflicting_qseqids = conflicts.filter(pl.col("len") > 1).get_column("qseqid")
    if not conflicting_qseqids.is_empty():
        message = "conflicting assignments for contigs: " + ", ".join(
            conflicting_qseqids.cast(pl.String).to_list(),
        )
        raise CrumbsProfileError(message)

    return assignments.unique(subset=["qseqid"], keep="first", maintain_order=True)


def prepare_coverage(coverage: pl.DataFrame, sample_id: str) -> pl.DataFrame:
    """Select and type CRUMBS coverage rows for one sample."""
    prepared = coverage.filter(pl.col("sample_id") == sample_id).select(
        COVERAGE_COLUMNS,
    )
    if prepared.is_empty():
        message = f"no coverage rows found for sample {sample_id}"
        raise CrumbsProfileError(message)

    duplicates = prepared.group_by("qseqid").len().filter(pl.col("len") > 1)
    if duplicates.height > 0:
        qseqids = duplicates.get_column("qseqid").to_list()
        message = (
            f"coverage TSV contains duplicate qseqid rows for sample {sample_id}: "
            + ", ".join(qseqids)
        )
        raise CrumbsProfileError(message)

    invalid = prepared.filter(
        (pl.col("qseqid").is_null())
        | (pl.col("qseqid").str.strip_chars() == "")
        | pl.any_horizontal(
            [pl.col(column).is_null() for column in COVERAGE_COLUMNS],
        )
        | (pl.col("contig_length") <= 0)
        | (pl.col("covered_bases_1x") < 0)
        | (pl.col("covered_bases_1x") > pl.col("contig_length"))
        | (pl.col("breadth_1x") < 0)
        | (pl.col("breadth_1x") > 1)
        | (
            (
                pl.col("breadth_1x")
                - (pl.col("covered_bases_1x") / pl.col("contig_length"))
            ).abs()
            > FLOAT_TOLERANCE
        )
        | (pl.col("raw_aligned_bases") < 0)
        | (pl.col("mean_depth_full") < 0)
        | (
            (
                pl.col("mean_depth_full")
                - (pl.col("raw_aligned_bases") / pl.col("contig_length"))
            ).abs()
            > FLOAT_TOLERANCE
        )
        | (pl.col("median_depth_full") < 0)
        | (pl.col("median_depth_positive") < 0)
        | (pl.col("depth_p95") < 0)
        | (pl.col("depth_p99") < 0)
        | (pl.col("max_depth") < 0)
        | (pl.col("depth_p95") > pl.col("depth_p99"))
        | (pl.col("depth_p99") > pl.col("max_depth"))
        | (pl.col("crumbs_p95") < 0)
        | (pl.col("crumbs_p99") < 0)
        | (pl.col("crumbs_p95") > pl.col("raw_aligned_bases"))
        | (pl.col("crumbs_p99") > pl.col("raw_aligned_bases"))
        | (pl.col("crumbs_p95") > (pl.col("depth_p95") * pl.col("contig_length")))
        | (pl.col("crumbs_p99") > (pl.col("depth_p99") * pl.col("contig_length"))),
    )
    if invalid.height > 0:
        qseqids = invalid.get_column("qseqid").fill_null("<null>").to_list()
        message = (
            "coverage TSV contains invalid coverage values for contigs: "
            + ", ".join(qseqids)
        )
        raise CrumbsProfileError(message)
    return prepared


def prepare_profile_taxonomy(profile_taxonomy: pl.DataFrame) -> pl.DataFrame:
    """Type and validate explicit profile taxonomy coordinates."""
    typed = profile_taxonomy.select(PROFILE_TAXONOMY_COLUMNS).with_columns(
        pl.all().cast(pl.String),
    )
    duplicate_taxa = typed.group_by("taxon_id").len().filter(pl.col("len") > 1)
    if duplicate_taxa.height > 0:
        message = "profile taxonomy contains duplicate taxon_id values"
        raise CrumbsProfileError(message)

    blank_rows = typed.filter(
        pl.any_horizontal(
            [
                pl.col(column).is_null() | (pl.col(column).str.strip_chars() == "")
                for column in PROFILE_TAXONOMY_COLUMNS
            ],
        ),
    )
    if blank_rows.height > 0:
        message = "profile taxonomy contains null or blank required values"
        raise CrumbsProfileError(message)

    with_paths = typed.with_columns(
        pl.col("taxpath").str.split("|").alias("taxpath_parts"),
        pl.col("taxpathsn").str.split("|").alias("taxpathsn_parts"),
        pl.col("rankpath").str.split("|").alias("rankpath_parts"),
    )
    length_mismatches = with_paths.filter(
        (pl.col("taxpath_parts").list.len() != pl.col("taxpathsn_parts").list.len())
        | (pl.col("taxpath_parts").list.len() != pl.col("rankpath_parts").list.len()),
    )
    if length_mismatches.height > 0:
        taxon_id = length_mismatches.get_column("taxon_id").item(0)
        message = f"profile taxonomy path length mismatch for taxon_id {taxon_id}"
        raise CrumbsProfileError(message)

    taxpath_terminal_mismatches = with_paths.filter(
        pl.col("taxpath_parts").list.last() != pl.col("taxon_id"),
    )
    if taxpath_terminal_mismatches.height > 0:
        taxon_id = taxpath_terminal_mismatches.get_column("taxon_id").item(0)
        message = f"profile taxonomy taxpath for taxon_id {taxon_id} does not end with taxon_id"
        raise CrumbsProfileError(message)

    taxpathsn_terminal_mismatches = with_paths.filter(
        pl.col("taxpathsn_parts").list.last() != pl.col("taxon_name"),
    )
    if taxpathsn_terminal_mismatches.height > 0:
        taxon_id = taxpathsn_terminal_mismatches.get_column("taxon_id").item(0)
        message = f"profile taxonomy taxpathsn for taxon_id {taxon_id} does not end with taxon_name"
        raise CrumbsProfileError(message)
    rankpath_terminal_mismatches = with_paths.filter(
        pl.col("rankpath_parts").list.last() != pl.col("rank"),
    )
    if rankpath_terminal_mismatches.height > 0:
        taxon_id = rankpath_terminal_mismatches.get_column("taxon_id").item(0)
        message = f"profile taxonomy rankpath for taxon_id {taxon_id} does not end with rank"
        raise CrumbsProfileError(message)
    return typed


def validate_assignment_coverage_alignment(
    assignments: pl.DataFrame,
    coverage: pl.DataFrame,
) -> None:
    """Fail when an assigned contig has no coverage summary row."""
    contig_assignments = contig_query_assignments(assignments)
    missing_coverage = contig_assignments.join(
        coverage.select("qseqid"),
        on="qseqid",
        how="anti",
    )
    if missing_coverage.height > 0:
        qseqids = missing_coverage.get_column("qseqid").to_list()
        message = "coverage TSV is missing assigned contigs: " + ", ".join(qseqids)
        raise CrumbsProfileError(message)


def read_query_assignment_mask(assignments: pl.DataFrame) -> pl.Expr:
    """Return an expression selecting read-derived query assignments."""
    prefix_expr = pl.any_horizontal(
        [pl.col("qseqid").str.starts_with(prefix) for prefix in READ_QUERY_PREFIXES],
    )
    if "query_class" not in assignments.columns:
        return prefix_expr
    return pl.col("query_class").is_in(READ_QUERY_CLASSES) | prefix_expr


def read_query_assignments(assignments: pl.DataFrame) -> pl.DataFrame:
    """Select BLAST assignments for read-derived query records."""
    return assignments.filter(read_query_assignment_mask(assignments))


def contig_query_assignments(assignments: pl.DataFrame) -> pl.DataFrame:
    """Select BLAST assignments that still require contig coverage rows."""
    return assignments.filter(~read_query_assignment_mask(assignments))


def build_read_query_table(
    *,
    sample_id: str,
    assignments: pl.DataFrame,
) -> pl.DataFrame:
    """Convert read-derived query assignments into CRUMBS evidence rows."""
    if assignments.is_empty():
        return pl.DataFrame(schema={column: pl.Null for column in CONTIG_OUTPUT_COLUMNS})
    missing_columns = [
        column
        for column in ("support_record_count", "qlen")
        if column not in assignments.columns
    ]
    if missing_columns:
        message = (
            "read-derived BLAST assignments missing required columns: "
            + ", ".join(missing_columns)
        )
        raise CrumbsProfileError(message)

    typed = assignments.with_columns(
        pl.col("support_record_count").cast(pl.Int64),
        pl.col("qlen").cast(pl.Int64),
    )
    invalid = typed.filter(
        pl.col("support_record_count").is_null()
        | (pl.col("support_record_count") <= 0)
        | pl.col("qlen").is_null()
        | (pl.col("qlen") <= 0),
    )
    if invalid.height > 0:
        qseqids = invalid.get_column("qseqid").to_list()
        message = (
            "read-derived BLAST assignments contain invalid support or qlen values: "
            + ", ".join(qseqids)
        )
        raise CrumbsProfileError(message)

    return typed.with_columns(
        pl.lit(sample_id).alias("sample_id"),
        pl.col("qlen").alias("contig_length"),
        pl.lit(0).alias("covered_bases_1x"),
        pl.lit(0.0).alias("breadth_1x"),
        pl.col("support_record_count").cast(pl.Float64).alias("raw_aligned_bases"),
        (pl.col("support_record_count") / pl.col("qlen")).cast(pl.Float64).alias("mean_depth_full"),
        pl.lit(0.0).alias("median_depth_full"),
        pl.lit(0.0).alias("median_depth_positive"),
        pl.lit(0.0).alias("depth_p95"),
        pl.lit(0.0).alias("depth_p99"),
        pl.lit(0.0).alias("max_depth"),
        pl.col("support_record_count").cast(pl.Float64).alias("crumbs_p95"),
        pl.col("support_record_count").cast(pl.Float64).alias("crumbs_p99"),
        pl.lit("read_query_support").alias("winsorization_percentile"),
        pl.col("support_record_count").cast(pl.Float64).alias("winsorization_cap"),
        pl.col("support_record_count").cast(pl.Float64).alias("crumbs_score"),
    ).select(CONTIG_OUTPUT_COLUMNS)


def build_contig_table(
    *,
    sample_id: str,
    assignments: pl.DataFrame,
    coverage: pl.DataFrame,
    profile_taxonomy: pl.DataFrame,
) -> pl.DataFrame:
    """Join assignments, coverage, and explicit taxonomy into contig evidence."""
    validate_assignment_coverage_alignment(assignments, coverage)
    contig_assignments = contig_query_assignments(assignments)
    contigs = (
        contig_assignments.join(coverage, on="qseqid", how="inner")
        .join(
            profile_taxonomy.select(
                "taxon_id",
                pl.lit(1).alias("has_profile_taxonomy"),
            ),
            left_on="adjusted_taxid",
            right_on="taxon_id",
            how="left",
        )
        .with_columns(
            pl.lit(sample_id).alias("sample_id"),
            pl.lit("p99").alias("winsorization_percentile"),
            pl.col("depth_p99").alias("winsorization_cap"),
            pl.col("crumbs_p99").alias("crumbs_score"),
        )
    )
    read_queries = build_read_query_table(
        sample_id=sample_id,
        assignments=read_query_assignments(assignments),
    )
    if read_queries.height > 0:
        contigs = pl.concat(
            [contigs.select(CONTIG_OUTPUT_COLUMNS), read_queries],
            how="vertical_relaxed",
        )
    contigs = contigs.select(CONTIG_OUTPUT_COLUMNS)
    missing_taxonomy = contigs.select("adjusted_taxid").unique().join(
        profile_taxonomy.select(pl.col("taxon_id").alias("adjusted_taxid")),
        on="adjusted_taxid",
        how="anti",
    )
    if missing_taxonomy.height > 0:
        missing = missing_taxonomy.get_column("adjusted_taxid").unique().to_list()
        message = "profile taxonomy is missing adjusted taxids: " + ", ".join(missing)
        raise CrumbsProfileError(message)

    return contigs


def build_taxon_table(
    *,
    sample_id: str,
    contigs: pl.DataFrame,
    profile_taxonomy: pl.DataFrame,
) -> pl.DataFrame:
    """Aggregate contig evidence into Taxon CRUMBS rows."""
    grouped = (
        contigs.group_by("adjusted_taxid")
        .agg(
            pl.len().alias("n_contigs"),
            pl.col("contig_length").sum().alias("total_contig_length"),
            pl.col("covered_bases_1x").sum().alias("total_covered_bases_1x"),
            pl.col("raw_aligned_bases").sum().alias("total_raw_aligned_bases"),
            pl.col("crumbs_score").sum().alias("total_crumbs_score"),
            (pl.col("crumbs_score") == 0).sum().alias("n_zero_crumbs_contigs"),
            pl.col("breadth_1x").median().alias("median_contig_breadth_1x"),
            pl.col("breadth_1x").min().alias("min_contig_breadth_1x"),
            pl.col("breadth_1x").max().alias("max_contig_breadth_1x"),
            pl.col("crumbs_score").max().alias("top_contig_crumbs"),
            pl.when(pl.col("breadth_1x") < LOW_BREADTH_THRESHOLD)
            .then(pl.col("crumbs_score"))
            .otherwise(0)
            .sum()
            .alias("low_breadth_crumbs"),
        )
        .rename({"adjusted_taxid": "taxon_id"})
        .with_columns(
            (pl.col("total_crumbs_score") / pl.col("total_contig_length")).alias(
                "taxon_crumbs",
            ),
        )
    )
    taxon_crumb_rows = list(
        grouped.select("taxon_id", "taxon_crumbs").iter_rows(named=True),
    )
    total_taxon_crumbs = math.fsum(row["taxon_crumbs"] for row in taxon_crumb_rows)
    if total_taxon_crumbs > 0:
        percentages = pl.DataFrame(
            {
                "taxon_id": [row["taxon_id"] for row in taxon_crumb_rows],
                "percentage_emitted": [
                    row["taxon_crumbs"] / total_taxon_crumbs * 100
                    for row in taxon_crumb_rows
                ],
            },
        )
        grouped = grouped.join(percentages, on="taxon_id", how="left")
    else:
        grouped = grouped.with_columns(pl.lit(0.0).alias("percentage_emitted"))

    grouped = grouped.with_columns(
        pl.when(pl.col("total_crumbs_score") > 0)
        .then(pl.col("top_contig_crumbs") / pl.col("total_crumbs_score"))
        .otherwise(0.0)
        .alias("fraction_crumbs_from_top_contig"),
        pl.when(pl.col("total_crumbs_score") > 0)
        .then(pl.col("low_breadth_crumbs") / pl.col("total_crumbs_score"))
        .otherwise(0.0)
        .alias("fraction_crumbs_from_low_breadth_contigs"),
    )

    return (
        grouped.join(profile_taxonomy, on="taxon_id", how="left")
        .with_columns(pl.lit(sample_id).alias("sample_id"))
        .select(TAXON_OUTPUT_COLUMNS)
        .sort("taxon_id")
    )


def bioboxes_rows(taxa: pl.DataFrame) -> pl.DataFrame:
    """Convert nonzero taxon rows to a BioBoxes profile table."""
    return taxa.filter(pl.col("percentage_emitted") > 0).select(
        pl.col("taxon_id").alias("@@TAXID"),
        pl.col("rank").alias("RANK"),
        pl.col("taxpath").alias("TAXPATH"),
        pl.col("taxpathsn").alias("TAXPATHSN"),
        pl.col("percentage_emitted").alias("PERCENTAGE"),
    )


def write_bioboxes_profile(sample_id: str, taxa: pl.DataFrame, path: Path) -> None:
    """Write a CAMI/BioBoxes-style taxonomic profile."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.write(f"@SampleID:{sample_id}\n")
        handle.write("@Version:0.9.1\n")
        bioboxes_rows(taxa).write_csv(handle, separator="\t")


def write_qc(
    sample_id: str,
    contigs: pl.DataFrame,
    taxa: pl.DataFrame,
    path: Path,
) -> None:
    """Write lightweight CRUMBS profile QC metadata."""
    qc = {
        "sample_id": sample_id,
        "crumbs_score_column": "crumbs_p99",
        "winsorization_rule": "per-contig full-depth nearest-rank p99, zeros included",
        "taxon_crumbs_rule": "sum(contig crumbs_p99) / sum(contig_length)",
        "is_genome_copy_abundance": False,
        "n_contigs": contigs.height,
        "n_taxa": taxa.height,
        "n_bioboxes_taxa": taxa.filter(pl.col("percentage_emitted") > 0).height,
        "low_breadth_threshold": LOW_BREADTH_THRESHOLD,
    }
    path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def estimate_profile(
    *,
    sample_id: str,
    blast_tsv: Path,
    coverage_tsv: Path,
    profile_taxonomy_tsv: Path,
    output_dir: Path,
) -> OutputPaths:
    """Estimate CRUMBS profile artifacts for one sample."""
    blast, coverage, profile_taxonomy = read_inputs(
        blast_tsv=blast_tsv,
        coverage_tsv=coverage_tsv,
        profile_taxonomy_tsv=profile_taxonomy_tsv,
    )
    assignments = collapse_assignments(blast, sample_id)
    typed_coverage = prepare_coverage(coverage, sample_id)
    typed_taxonomy = prepare_profile_taxonomy(profile_taxonomy)
    contigs = build_contig_table(
        sample_id=sample_id,
        assignments=assignments,
        coverage=typed_coverage,
        profile_taxonomy=typed_taxonomy,
    )
    taxa = build_taxon_table(
        sample_id=sample_id,
        contigs=contigs,
        profile_taxonomy=typed_taxonomy,
    )

    paths = output_paths(sample_id, output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    contigs.write_csv(paths.contigs, separator="\t")
    taxa.write_csv(paths.taxa, separator="\t")
    write_bioboxes_profile(sample_id, taxa, paths.bioboxes_profile)
    write_qc(sample_id, contigs, taxa, paths.qc_json)
    return paths


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Estimate CRUMBS taxonomic profiles")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--blast-tsv", type=Path, required=True)
    parser.add_argument("--coverage-tsv", type=Path, required=True)
    parser.add_argument("--profile-taxonomy-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    _ = estimate_profile(
        sample_id=args.sample_id,
        blast_tsv=args.blast_tsv,
        coverage_tsv=args.coverage_tsv,
        profile_taxonomy_tsv=args.profile_taxonomy_tsv,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
