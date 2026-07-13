#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
#     "polars",
#     "taxopy",
# ]
# ///

"""
Annotate BLAST hits with taxonomic consensus assignments.

This module resolves retained BLAST references into either a dominant taxid or
a lowest common ancestor (LCA). References scoring at least 95% as highly as
the best retained reference determine the assignment: one distinct taxid is
dominant, while multiple taxids resolve to their strict NCBI LCA.

Key features:
- Relative bitscore selection of references used for assignment
- Dominant assignment when selected references agree on one taxid
- Strict LCA assignment when selected references represent multiple taxids
- Automatic NCBI taxonomy database management via py_nvd.taxonomy
- Exact name and rank annotation for the selected taxid

Graceful Degradation:
    By default, if the pre-cached taxonomy database is unavailable, this script
    will warn and attempt to download fresh taxonomy data from NCBI. Use --sync
    to require pre-cached taxonomy and fail if unavailable.

Usage:
    python annotate_blast_lca.py -i blast_results.txt -o annotated.txt

    # Require pre-cached taxonomy (fail if unavailable)
    python annotate_blast_lca.py -i blast_results.txt -o annotated.txt --sync
"""

import argparse
from pathlib import Path

import polars as pl
import taxopy.utilities as taxutils
from loguru import logger
from py_nvd import taxonomy
from py_nvd.blast_consensus import ASSIGNMENT_BITSCORE_FRACTION
from taxopy import TaxDb, Taxon
from taxopy.exceptions import TaxidError

INPUT_COLUMNS = [
    "task",
    "sample",
    "qseqid",
    "qlen",
    "sseqid",
    "stitle",
    "length",
    "pident",
    "evalue",
    "bitscore",
    "sscinames",
    "staxids",
    "rank",
]
LCA_COLUMNS = [
    "adjusted_taxid",
    "adjusted_taxid_name",
    "adjusted_taxid_rank",
    "adjustment_method",
]
OUTPUT_COLUMNS = [*INPUT_COLUMNS, *LCA_COLUMNS]


def input_has_data_rows(path: Path) -> bool:
    """Return true when a TSV has at least one row after its header."""
    if path.stat().st_size == 0:
        return False
    with path.open(encoding="utf-8") as handle:
        next(handle, None)
        return next(handle, None) is not None


def write_empty_lca_output(path: Path) -> None:
    """Write a header-only LCA TSV for valid no-hit samples."""
    path.write_text("\t".join(OUTPUT_COLUMNS) + "\n", encoding="utf-8")


def with_score_context(hits: pl.LazyFrame) -> pl.LazyFrame:
    """Mark retained references whose bitscores participate in assignment."""
    return hits.with_columns(
        pl.max("bitscore").over(["task", "sample", "qseqid"]).alias("best_bitscore"),
    ).with_columns(
        (
            pl.col("bitscore") >= pl.col("best_bitscore") * ASSIGNMENT_BITSCORE_FRACTION
        ).alias("used_for_assignment"),
    )


def _canonical_taxid(tid: int | None, tx: TaxDb) -> int | None:
    """Resolve one merged taxid, returning null when taxonomy has no match."""
    if tid is None:
        return None

    canonical_tid = tid
    merged_taxids = tx.oldtaxid2newtaxid or {}
    seen: set[int] = set()
    while canonical_tid in merged_taxids and canonical_tid not in seen:
        seen.add(canonical_tid)
        canonical_tid = int(merged_taxids[canonical_tid])

    try:
        return int(Taxon(canonical_tid, tx).taxid)
    except (KeyError, ValueError, TaxidError):
        return None


def with_canonical_taxids(hits: pl.LazyFrame, tx: TaxDb) -> pl.LazyFrame:
    """Normalize retained-row taxids and collapse resulting duplicate rows."""
    return hits.with_columns(
        pl.col("staxids")
        .map_elements(
            lambda tid: _canonical_taxid(tid, tx),
            return_dtype=pl.Int64,
        )
        .alias("staxids"),
    ).unique(maintain_order=True)


def _find_lca(
    canonical_taxids: list[int] | pl.Series,
    tx: TaxDb,
) -> int:
    """Return the strict LCA of canonical, taxonomy-backed taxids."""
    values = (
        canonical_taxids.to_list()
        if isinstance(canonical_taxids, pl.Series)
        else canonical_taxids
    )
    if len(values) == 1:
        return values[0]

    taxons = [Taxon(tid, tx) for tid in values]
    return taxutils.find_lca(taxons, tx).taxid


def _annotate_assignments(assignments: pl.LazyFrame, tx: TaxDb) -> pl.LazyFrame:
    """Add the exact assigned taxid's scientific name and NCBI rank."""
    id2name = {int(t): name for t, name in tx.taxid2name.items()}
    id2rank = {int(t): rank for t, rank in tx.taxid2rank.items()}

    return assignments.with_columns(
        pl.col("adjusted_taxid")
        .replace_strict(id2name, default=None)
        .alias("adjusted_taxid_name"),
        pl.col("adjusted_taxid")
        .replace_strict(id2rank, default=None)
        .alias("adjusted_taxid_rank"),
    )


def assign_taxonomic_consensus(hits: pl.LazyFrame, tx: TaxDb) -> pl.LazyFrame:
    """
    Assign one taxid from references selected by relative bitscore.

    Args:
        hits: Retained BLAST hits with score-context columns
        tx: TaxDb instance with loaded taxonomy

    Returns:
        LazyFrame joined with assignment results, including:
            - adjusted_taxid: Consensus taxid (dominant or LCA)
            - adjusted_taxid_name: Scientific name
            - adjusted_taxid_rank: Refined taxonomic rank
            - adjustment_method: "dominant" or "lca"

    Side effects: None (pure function, returns lazy computation graph)

    Note: This is the core taxonomic consensus algorithm.
    """
    assignment_taxids = (
        hits.filter(pl.col("used_for_assignment"))
        .group_by(["task", "sample", "qseqid"])
        .agg(
            assignment_taxids=pl.col("staxids").drop_nulls().unique().sort(),
        )
        .filter(pl.col("assignment_taxids").list.len() > 0)
        .with_columns(
            pl.col("assignment_taxids")
            .map_elements(lambda tids: _find_lca(tids, tx), return_dtype=pl.Int64)
            .alias("adjusted_taxid"),
            pl.when(pl.col("assignment_taxids").list.len() == 1)
            .then(pl.lit("dominant"))
            .otherwise(pl.lit("lca"))
            .alias("adjustment_method"),
        )
        .select(["task", "sample", "qseqid", "adjusted_taxid", "adjustment_method"])
    )

    assignments = _annotate_assignments(assignment_taxids, tx)

    return hits.join(
        assignments,
        how="left",
        on=["task", "sample", "qseqid"],
    ).drop(
        "used_for_assignment",
        "best_bitscore",
    )


def require_taxonomic_assignments(assignments: pl.LazyFrame) -> None:
    """Fail when retained assignment references have no canonical taxid."""
    missing = (
        assignments.filter(pl.col("adjusted_taxid").is_null())
        .select("task", "sample", "qseqid")
        .unique()
        .with_columns(
            pl.concat_str(
                [
                    pl.col("task"),
                    pl.lit("/"),
                    pl.col("sample"),
                    pl.lit("/"),
                    pl.col("qseqid"),
                ],
            ).alias("query"),
        )
        .select("query")
        .collect()
    )
    if missing.height == 0:
        return

    queries = ", ".join(missing.get_column("query").to_list())
    message = (
        "Assignment references had no canonical taxid for: "
        f"{queries}. Inspect the pre-LCA combined-search rows and taxonomy snapshot."
    )
    raise RuntimeError(message)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for LCA annotation.

    Returns:
        Namespace with input, output, and taxonomy paths

    Side effects:
        - Reads sys.argv
        - May exit program if arguments invalid or --help requested
    """
    parser = argparse.ArgumentParser(
        description=(
            "Annotate retained BLAST references with a dominant taxid or the "
            "lowest common ancestor of near-best taxids."
        ),
    )

    # Input/Output arguments
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="Path to the input BLAST results text file.",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        required=True,
        help="Path to write the annotated BLAST results.",
    )

    parser.add_argument(
        "--taxonomy-dir",
        type=str,
        default=None,
        help="Explicit taxonomy directory",
    )
    parser.add_argument(
        "--sync",
        action="store_true",
        help="Require pre-cached taxonomy database (fail if unavailable)",
    )

    return parser.parse_args()


def main() -> None:
    """
    Main entry point for BLAST LCA annotation.

    Orchestrates taxonomic consensus annotation:
    1. Parse command-line arguments
    2. Load NCBI taxonomy database via py_nvd.taxonomy
    3. Mark references used for assignment from relative bitscores
    4. Assign consensus taxonomic IDs (dominant or LCA)
    5. Write every retained BLAST row with the query-level assignment

    Side effects:
        - Reads input BLAST file
        - May download NCBI taxonomy database (~50MB) on first run
        - Writes output file
        - Logs progress to stderr via loguru

    Raises:
        SystemExit: On invalid arguments or I/O errors
    """
    args = parse_args()

    logger.info("Loading NCBI taxonomy database...")
    input_path = Path(args.input_file)
    output_path = Path(args.output_file)
    if not input_has_data_rows(input_path):
        logger.warning(
            "Input file has no BLAST data rows. Creating header-only output file.",
        )
        write_empty_lca_output(output_path)
        return

    with taxonomy.open(
        taxonomy_dir=args.taxonomy_dir,
        sync=args.sync,
    ) as tax:
        tx = tax.taxopy_db

        logger.info(
            "Assigning taxonomy from references scoring at least {:.0%} of the best hit",
            ASSIGNMENT_BITSCORE_FRACTION,
        )
        hits = (
            pl.scan_csv(args.input_file, separator="\t")
            .pipe(with_canonical_taxids, tx)
            .pipe(with_score_context)
        )
        assignment_lf = assign_taxonomic_consensus(hits, tx)
        require_taxonomic_assignments(assignment_lf)

        logger.info("Writing annotated results to {}", args.output_file)
        assignment_lf.sink_csv(args.output_file, separator="\t")


if __name__ == "__main__":
    main()
