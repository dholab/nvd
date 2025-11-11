#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "loguru",
#     "polars",
#     "taxopy",
# ]
# ///


import argparse
import tarfile
import urllib.request
from pathlib import Path

import polars as pl
import polars.selectors as cs
import taxopy.utilities as taxutils
from loguru import logger
from taxopy import TaxDb, Taxon

# URL for NCBI taxdump files
NCBI_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Parameters (tune as needed)
MIN_BP = 70  # or 150 to 300 for contigs
MIN_ID = 90.0  # for blastn; adjust for aa
MAX_E = 1e-10
DELTA_S = 5.0  # ΔS window for “near ties”
MIN_SUPPORT = 0.8  # dominance threshold


def ensure_taxdump(dest_dir: str | Path = "taxdump") -> Path:
    """
    Download and extract NCBI taxdump files if they don't already exist locally.
    Returns the Path to the directory containing nodes.dmp, names.dmp, merged.dmp.
    """
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)

    needed = {"nodes.dmp", "names.dmp", "merged.dmp"}
    existing = {p.name for p in dest.iterdir() if p.is_file()}

    if needed.issubset(existing):
        logger.info("Taxdump files already present in %s", dest)
        return dest

    logger.info("Downloading NCBI taxdump to %s ...", dest)
    tmp_tar = dest / "taxdump.tar.gz"

    with (
        urllib.request.urlopen(NCBI_TAXDUMP_URL) as resp,  # noqa: S310
        open(tmp_tar, "wb") as f,
    ):
        f.write(resp.read())

    with tarfile.open(tmp_tar, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.name in needed]
        tar.extractall(dest, members=members)

    tmp_tar.unlink(missing_ok=True)
    logger.info("Taxdump downloaded and extracted.")
    return dest


def load_taxdb(dest_dir: str | Path = "taxdump") -> TaxDb:
    """
    Ensures the taxonomy dump is available, then loads and returns a TaxDb instance.
    """
    dest = ensure_taxdump(dest_dir)
    return TaxDb(
        nodes_dmp=dest / "nodes.dmp",
        names_dmp=dest / "names.dmp",
        merged_dmp=dest / "merged.dmp",
    )


def filter_blast_hits(unfiltered_hits: pl.LazyFrame) -> pl.LazyFrame:
    return (
        unfiltered_hits.filter(
            pl.col("length") >= MIN_BP,
            pl.col("pident") >= MIN_ID,
            pl.col("evalue") <= MAX_E,
        )
        # Compute Smax per query
        .with_columns(pl.max("bitscore").over(["task", "sample", "qseqid"]).alias("Smax"))
        # Near-tie mask: ΔS <= DELTA_S
        .with_columns((pl.col("bitscore") >= (pl.col("Smax") - DELTA_S)).alias("near_tie"))
    )


def _find_lca_list(tids: list[int], tx: TaxDb) -> int | None:
    if len(tids) == 1:
        return tids[0]
    if len(tids) == 0:
        return None
    tids = [Taxon(tid, tx) for tid in tids]
    return taxutils.find_lca(tids, tx).taxid


def _nearest_ranked(tid: int | None, tx: TaxDb) -> tuple[int, str] | tuple[None, None]:
    # Bail early on obviously invalid input
    if tid is None or not isinstance(tid, int):
        return None, None

    try:
        t = Taxon(tid, tx)
    except (KeyError, ValueError):
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


def refine_assignments(assignments: pl.LazyFrame, tx: TaxDb) -> pl.LazyFrame:
    id2name = {int(t): name for t, name in tx.taxid2name.items()}
    id2rank = {int(t): rank for t, rank in tx.taxid2rank.items()}

    needs_nearest = pl.col("adjusted_taxid_rank").is_null() | pl.col("adjusted_taxid_rank").is_in(
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


def adjust_taxids(filtered_lf: pl.LazyFrame, tx: TaxDb) -> pl.LazyFrame:
    no_lca_needed = (
        filtered_lf.group_by(["task", "sample", "qseqid", "staxids"])
        .agg(total_w=pl.col("bitscore").sum())
        .with_columns(pl.sum("total_w").over(["task", "sample", "qseqid"]).alias("Wsum"))
        .with_columns(
            (pl.col("total_w") / pl.col("Wsum")).alias("support"),
        )
        .filter(pl.col("support") >= MIN_SUPPORT)
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
            .map_elements(_find_lca_list, return_dtype=pl.Int64)
            .alias("adjusted_taxid"),
            pl.lit("lca").alias("adjustment_method"),
        )
        .select(["task", "sample", "qseqid", "adjusted_taxid", "adjustment_method"])
    )

    assignments = refine_assignments(pl.concat([no_lca_needed, lca_assigned]), tx)

    return filtered_lf.join(assignments, how="left", on=["task", "sample", "qseqid"]).drop(
        "near_tie",
        "Smax",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate contigs with conflicting BLAST hits with a 'consensus' taxid, \
        either by taking the best it or by computing the Least Common Ancestor for near ties.",
    )

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
        help="Path to write the annotated BLAST results (top hits).",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tx = load_taxdb()

    unfiltered_hits = pl.scan_csv(args.input_file, separator="\t")
    filtered_lf = filter_blast_hits(unfiltered_hits)
    assignment_lf = adjust_taxids(filtered_lf, tx)

    assignment_lf.sink_csv(args.out_file, separator="\t")


if __name__ == "__main__":
    main()
