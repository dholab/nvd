"""Tests for pruning MEGABLAST-classified query sequences."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING
from unittest.mock import patch

from remove_megablast_mapped_contigs import main

if TYPE_CHECKING:
    from pathlib import Path


def test_prunes_stable_query_ids_from_annotated_fasta_headers(
    tmp_path: Path,
) -> None:
    megablast_results = tmp_path / "megablast.tsv"
    megablast_results.write_text(
        "nvdContig1_sample-1.2_000001\t4\tref|NC_000001.1|\n",
        encoding="utf-8",
    )
    queries = tmp_path / "queries.fasta"
    queries.write_text(
        ">nvdContig1_sample-1.2_000001 evidence_class=assembly_contig producer=spades contig_id=NODE_1\n"
        "AAAA\n"
        ">nvdContig1_sample-1.2_000002 evidence_class=assembly_contig producer=spades contig_id=NODE_2\n"
        "CCCC\n"
        ">nvdContig1_sample-1.2_0000010 evidence_class=assembly_contig producer=spades contig_id=NODE_10\n"
        "GGGG\n",
        encoding="utf-8",
    )
    classified = tmp_path / "classified.txt"
    pruned = tmp_path / "pruned.fasta"
    with patch.object(
        sys,
        "argv",
        [
            "remove_megablast_mapped_contigs.py",
            "--megablast_results",
            str(megablast_results),
            "--contigs_fasta",
            str(queries),
            "--classified_contigs",
            str(classified),
            "--pruned_contigs",
            str(pruned),
        ],
    ):
        main()

    assert classified.read_text(encoding="utf-8") == ("nvdContig1_sample-1.2_000001\n")
    assert pruned.read_text(encoding="utf-8") == (
        ">nvdContig1_sample-1.2_000002 evidence_class=assembly_contig producer=spades contig_id=NODE_2\n"
        "CCCC\n"
        ">nvdContig1_sample-1.2_0000010 evidence_class=assembly_contig producer=spades contig_id=NODE_10\n"
        "GGGG\n"
    )
