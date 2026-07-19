"""Tests for partitioning MEGABLAST query sequences."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING
from unittest.mock import patch

from partition_megablast_queries import main

if TYPE_CHECKING:
    from pathlib import Path


def test_partitions_stable_query_ids_from_annotated_fasta_headers(
    tmp_path: Path,
) -> None:
    megablast_results = tmp_path / "megablast.tsv"
    megablast_results.write_text(
        "nvdContig1_sample-1.2_000001\t4\tref|NC_000001.1|\n",
        encoding="utf-8",
    )
    queries = tmp_path / "queries.fasta"
    queries.write_text(
        ">nvdContig1_sample-1.2_000001 query_class=assembly_contig producer=spades contig_id=NODE_1\n"
        "AAAA\n"
        ">nvdContig1_sample-1.2_000002 query_class=assembly_contig producer=spades contig_id=NODE_2\n"
        "CCCC\n"
        ">nvdContig1_sample-1.2_0000010 query_class=assembly_contig producer=spades contig_id=NODE_10\n"
        "GGGG\n",
        encoding="utf-8",
    )
    accounted_query_ids = tmp_path / "accounted_query_ids.txt"
    blastn_candidates = tmp_path / "blastn_candidates.fasta"
    with patch.object(
        sys,
        "argv",
        [
            "partition_megablast_queries.py",
            "--megablast-results",
            str(megablast_results),
            "--query-fasta",
            str(queries),
            "--accounted-query-ids",
            str(accounted_query_ids),
            "--blastn-candidate-fasta",
            str(blastn_candidates),
        ],
    ):
        main()

    assert accounted_query_ids.read_text(encoding="utf-8") == (
        "nvdContig1_sample-1.2_000001\n"
    )
    assert blastn_candidates.read_text(encoding="utf-8") == (
        ">nvdContig1_sample-1.2_000002 query_class=assembly_contig producer=spades contig_id=NODE_2\n"
        "CCCC\n"
        ">nvdContig1_sample-1.2_0000010 query_class=assembly_contig producer=spades contig_id=NODE_10\n"
        "GGGG\n"
    )
