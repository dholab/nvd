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

This module resolves ambiguous BLAST results by computing either dominant taxid
assignments (when one taxid has >80% support) or Least Common Ancestor (LCA)
for near-tie cases. The LCA approach minimizes over-specificity when multiple
closely-scoring hits disagree at the species level.

Key features:
- Quality filtering by alignment length, percent identity, and e-value
- Near-tie detection using bitscore delta windows
- Dual assignment strategy: dominant hits vs. LCA calculation
- Automatic NCBI taxonomy database management via py_nvd.taxonomy
- Taxonomic rank refinement to avoid uninterpretable levels

Usage:
    python annotate_blast_lca.py -i blast_results.txt -o annotated.txt
"""

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import polars as pl
import polars.selectors as cs
import taxopy.utilities as taxutils
from loguru import logger
from py_nvd import taxonomy
from taxopy import TaxDb, Taxon
from taxopy.exceptions import TaxidError

# Parameters (tune as needed)
DEFAULT_MIN_BP = 70  # or 150 to 300 for contigs
DEFAULT_MIN_ID = 90.0  # for blastn; adjust for aa
DEFAULT_MAX_E = 1e-10
DEFAULT_DELTA_S = 5.0  # ΔS window for "near ties"
DEFAULT_MIN_SUPPORT = 0.8  # dominance threshold


@dataclass
class LcaParams:
    """
    Parameters for LCA-based taxonomic assignment.

    Attributes:
        min_bp: Minimum alignment length in base pairs (filters short/spurious hits)
        min_identity: Minimum percent identity threshold (0-100)
        max_evalue: Maximum e-value for considering hits (lower = more stringent)
        delta_s_window: Bitscore window for near-tie detection (ΔS threshold)
        min_support: Minimum support fraction for dominant assignment (0-1)

    Note: This is a pure data class with no side effects.
    """

    min_bp: int = DEFAULT_MIN_BP
    min_identity: float = DEFAULT_MIN_ID
    max_evalue: float = DEFAULT_MAX_E
    delta_s_window: float = DEFAULT_DELTA_S
    min_support: float = DEFAULT_MIN_SUPPORT


def filter_blast_hits(unfiltered_hits: pl.LazyFrame, params: LcaParams) -> pl.LazyFrame:
    """
    Filter BLAST hits by quality thresholds and identify near-tie cases.

    Applies quality filters (length, identity, e-value), computes the maximum
    bitscore (Smax) per query, and flags hits within delta_s_window of Smax
    as "near ties" for potential LCA calculation.

    Args:
        unfiltered_hits: LazyFrame of raw BLAST results
        params: Quality and near-tie detection parameters

    Returns:
        Filtered LazyFrame with added columns:
            - Smax: Maximum bitscore for each query
            - near_tie: Boolean flag for hits within ΔS window

    Side effects: None (pure function, returns lazy computation graph)

    Expected input columns: length, pident, evalue, bitscore, task, sample, qseqid
    """
    return (
        unfiltered_hits
        # Compute Smax per query
        .with_columns(
            pl.max("bitscore").over(["task", "sample", "qseqid"]).alias("Smax"),
        )
        # Near-tie mask: ΔS <= delta_s_window
        .with_columns(
            (pl.col("bitscore") >= (pl.col("Smax") - params.delta_s_window)).alias(
                "near_tie",
            ),
        )
    )


def _find_lca_list(tids: list[int | None] | None, tx: TaxDb) -> int | None:
    """
    Compute the Least Common Ancestor (LCA) for a list of taxids.

    Finds the most recent common ancestor in the NCBI taxonomy tree.
    Handles edge cases: single taxid returns itself, empty list returns None.

    Args:
        tids: List of NCBI taxonomy IDs
        tx: TaxDb instance with loaded taxonomy

    Returns:
        Taxid of LCA, or None if list is empty

    Side effects: None (pure function)

    Note: Private helper function, intended for use within map_elements.
    """
    if tids is None:
        return None
    if len(tids) == 0:
        return None
    if len(tids) == 1:
        return tids[0]

    # Filter out invalid taxids
    valid_taxons = []
    for tid in tids:
        if tid is None:
            continue
        try:
            valid_taxons.append(Taxon(tid, tx))
        except (KeyError, ValueError, TaxidError):
            # Skip invalid taxids - could be merged/obsolete/corrupted
            logger.warning(f"Skipping invalid taxid: {tid}")
            continue

    # Edge cases after filtering
    if len(valid_taxons) == 0:
        logger.warning(f"No valid taxids found in list: {tids}")
        return None
    if len(valid_taxons) == 1:
        return valid_taxons[0].taxid

    return taxutils.find_lca(valid_taxons, tx).taxid


def _nearest_ranked(tid: int | None, tx: TaxDb) -> tuple[int, str] | tuple[None, None]:
    """
    Find the nearest ancestor with a proper taxonomic rank.

    Climbs the taxonomic lineage to find the first ancestor with a standard
    rank (strain to family), avoiding unhelpful ranks like "clade" or "no rank".
    Stops climbing above family level to prevent over-generalization.

    Args:
        tid: NCBI taxonomy ID to refine, or None
        tx: TaxDb instance with loaded taxonomy

    Returns:
        Tuple of (taxid, rank) for nearest ranked ancestor, or (None, None) if:
            - Input tid is None/invalid
            - Taxid not found in database
            - No lineage information available

    Side effects: None (pure function)

    Acceptable ranks: strain, species, subspecies, subgenus, genus, subfamily, family

    Note: Private helper function for rank refinement.
    """
    # Bail early on obviously invalid input
    if tid is None or not isinstance(tid, int):
        return None, None

    try:
        t = Taxon(tid, tx)
    except (KeyError, ValueError, TaxidError):
        # e.g., taxid not found in nodes.dmp or merged.dmp
        return None, None

    if not getattr(t, "taxid_lineage", None):
        return None, None

    # Walk lineage from leaf to root to find first properly ranked node,
    # but stop climbing above family. If none found, keep the original.
    for anc_tid in reversed(t.taxid_lineage):
        anc = Taxon(anc_tid, tx)
        if anc.rank not in {"no rank", "clade"}:
            if anc.rank not in {
                "strain",
                "species",
                "subspecies",
                "subgenus",
                "genus",
                "subfamily",
                "family",
            }:
                # too high: keep the original
                return t.taxid, t.rank
            return anc.taxid, anc.rank

    # Fallback: no ranked ancestor found, return original taxon
    return t.taxid, t.rank


def _refine_assignments(assignments: pl.LazyFrame, tx: TaxDb) -> pl.LazyFrame:
    """
    Refine taxonomic assignments by adding names, ranks, and fixing rank issues.

    Enriches assignments with taxonomic names and ranks. For taxids with
    unhelpful ranks ("clade", "no rank"), climbs to nearest properly ranked
    ancestor using _nearest_ranked.

    Args:
        assignments: LazyFrame with adjusted_taxid column
        tx: TaxDb instance with loaded taxonomy

    Returns:
        LazyFrame with added columns:
            - adjusted_taxid_name: Scientific name for taxid
            - adjusted_taxid_rank: Taxonomic rank (refined to standard levels)

    Side effects: None (pure function, returns lazy computation graph)

    Note: Private helper function for taxonomic enrichment.
    """
    id2name = {int(t): name for t, name in tx.taxid2name.items()}
    id2rank = {int(t): rank for t, rank in tx.taxid2rank.items()}

    needs_nearest = pl.col("adjusted_taxid_rank").is_null() | pl.col(
        "adjusted_taxid_rank",
    ).is_in(
        ["clade", "no rank"],
    )

    return (
        assignments.with_columns(
            pl.col("adjusted_taxid")
            .replace_strict(id2name, default=None)
            .alias("adjusted_taxid_name"),
            pl.col("adjusted_taxid")
            .replace_strict(id2rank, default=None)
            .alias("adjusted_taxid_rank"),
        )
        .with_columns(
            pl.when(needs_nearest)
            .then(
                pl.col("adjusted_taxid").map_elements(
                    lambda t: _nearest_ranked(t, tx)[0],
                    return_dtype=pl.Int64,
                ),
            )
            .otherwise(None)
            .alias("_nearest_ranked_taxid"),
        )
        .with_columns(
            # overwrite rank where needed
            pl.when(needs_nearest)
            .then(pl.col("_nearest_ranked_taxid").replace_strict(id2rank, default=None))
            .otherwise(pl.col("adjusted_taxid_rank"))
            .alias("adjusted_taxid_rank"),
            # overwrite name where needed
            pl.when(needs_nearest)
            .then(pl.col("_nearest_ranked_taxid").replace_strict(id2name, default=None))
            .otherwise(pl.col("adjusted_taxid_name"))
            .alias("adjusted_taxid_name"),
        )
        .drop(cs.starts_with("_"))
    )


def call_least_common_ancestors(
    filtered_lf: pl.LazyFrame,
    tx: TaxDb,
    params: LcaParams,
) -> pl.LazyFrame:
    """
    Assign consensus taxonomic IDs using dominant hits or LCA.

    Implements a two-path assignment strategy:
    1. Dominant path: If one taxid has ≥min_support fraction of total bitscore
       weight, assign it directly (no LCA needed)
    2. LCA path: For queries without dominant taxid, compute LCA of all
       near-tie hits (within delta_s_window of Smax)

    After assignment, refines ranks and adds taxonomic names.

    Args:
        filtered_lf: Filtered BLAST hits with near_tie and Smax columns
        tx: TaxDb instance with loaded taxonomy
        params: Parameters including min_support threshold

    Returns:
        LazyFrame joined with assignment results, including:
            - adjusted_taxid: Consensus taxid (dominant or LCA)
            - adjusted_taxid_name: Scientific name
            - adjusted_taxid_rank: Refined taxonomic rank
            - adjustment_method: "dominant" or "lca"

    Side effects: None (pure function, returns lazy computation graph)

    Note: This is the core taxonomic consensus algorithm.
    """
    no_lca_needed = (
        filtered_lf.group_by(["task", "sample", "qseqid", "staxids"])
        .agg(total_w=pl.col("bitscore").sum())
        .with_columns(
            pl.sum("total_w").over(["task", "sample", "qseqid"]).alias("Wsum"),
        )
        .with_columns(
            (pl.col("total_w") / pl.col("Wsum")).alias("support"),
        )
        .filter(pl.col("support") >= params.min_support)
        .sort(
            ["task", "sample", "qseqid", "support"],
            descending=[False, False, False, True],
        )
        .group_by(["task", "sample", "qseqid"])
        .agg(
            adjusted_taxid=pl.col("staxids").first(),
            adjustment_method=pl.lit("dominant"),
        )
    )

    lca_assigned = (
        filtered_lf.filter(pl.col("near_tie"))  # only the ΔS window
        .group_by(["task", "sample", "qseqid"])
        .agg(staxids=pl.col("staxids").unique().sort())
        .join(
            no_lca_needed.select(["task", "sample", "qseqid"]),
            on=["task", "sample", "qseqid"],
            how="anti",
        )
        .with_columns(
            pl.col("staxids")
            .map_elements(lambda tids: _find_lca_list(tids, tx), return_dtype=pl.Int64)
            .alias("adjusted_taxid"),
            pl.lit("lca").alias("adjustment_method"),
        )
        .select(["task", "sample", "qseqid", "adjusted_taxid", "adjustment_method"])
    )

    assignments = _refine_assignments(pl.concat([no_lca_needed, lca_assigned]), tx)

    return filtered_lf.join(
        assignments,
        how="left",
        on=["task", "sample", "qseqid"],
    ).drop(
        "near_tie",
        "Smax",
    )


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for LCA annotation.

    Returns:
        Namespace with parsed arguments including I/O paths and LCA parameters

    Side effects:
        - Reads sys.argv
        - May exit program if arguments invalid or --help requested
    """
    parser = argparse.ArgumentParser(
        description="Annotate contigs with conflicting BLAST hits with a 'consensus' taxid, "
        "either by taking the best hit or by computing the Least Common Ancestor for near ties.",
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

    # LCA parameter arguments
    parser.add_argument(
        "--min-bp",
        type=int,
        default=DEFAULT_MIN_BP,
        help=f"Minimum alignment length in base pairs (default: {DEFAULT_MIN_BP})",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=DEFAULT_MIN_ID,
        help=f"Minimum percent identity for BLAST hits (default: {DEFAULT_MIN_ID})",
    )
    parser.add_argument(
        "--max-evalue",
        type=float,
        default=DEFAULT_MAX_E,
        help=f"Maximum e-value threshold (default: {DEFAULT_MAX_E})",
    )
    parser.add_argument(
        "--delta-s",
        type=float,
        default=DEFAULT_DELTA_S,
        help=f"Bitscore window for defining near-ties (default: {DEFAULT_DELTA_S})",
    )
    parser.add_argument(
        "--min-support",
        type=float,
        default=DEFAULT_MIN_SUPPORT,
        help=f"Minimum support fraction for dominant taxid assignment (default: {DEFAULT_MIN_SUPPORT})",
    )

    return parser.parse_args()


def main() -> None:
    """
    Main entry point for BLAST LCA annotation.

    Orchestrates the full pipeline:
    1. Parse command-line arguments
    2. Load NCBI taxonomy database via py_nvd.taxonomy
    3. Filter BLAST hits by quality thresholds
    4. Assign consensus taxonomic IDs (dominant or LCA)
    5. Write annotated results to output file

    Side effects:
        - Reads input BLAST file
        - May download NCBI taxonomy database (~50MB) on first run
        - Writes output file
        - Logs progress to stderr via loguru

    Raises:
        SystemExit: On invalid arguments or I/O errors
    """
    args = parse_args()

    # Construct LcaParams from command-line arguments
    params = LcaParams(
        min_bp=args.min_bp,
        min_identity=args.min_identity,
        max_evalue=args.max_evalue,
        delta_s_window=args.delta_s,
        min_support=args.min_support,
    )

    logger.info("Loading NCBI taxonomy database...")
    # Handle empty input file
    if os.path.getsize(args.input_file) == 0:
        logger.warning("Input file is empty. Creating empty output file.")
        Path(args.output_file).touch()
        return

    with taxonomy.open() as tax:
        tx = tax.taxopy_db

        logger.info("Processing BLAST results with parameters: {}", params)
        unfiltered_hits = pl.scan_csv(args.input_file, separator="\t")
        filtered_lf = filter_blast_hits(unfiltered_hits, params)
        assignment_lf = call_least_common_ancestors(filtered_lf, tx, params)

        logger.info("Writing annotated results to {}", args.output_file)
        assignment_lf.sink_csv(args.output_file, separator="\t")


if __name__ == "__main__":
    main()
