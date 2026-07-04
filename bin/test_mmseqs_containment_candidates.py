"""Regression tests for MMseqs2 exact-containment candidate generation."""

from __future__ import annotations

import shutil
import subprocess
from typing import TYPE_CHECKING

import pytest
from union_contigs import reverse_complement

if TYPE_CHECKING:
    from pathlib import Path


def write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")
    return path


def run_mmseqs_containment(
    *,
    fasta: Path,
    output: Path,
    tmp_path: Path,
    max_seq_len: int,
) -> None:
    mmseqs = shutil.which("mmseqs")
    if mmseqs is None:
        pytest.skip("MMseqs2 is not available in the test environment")

    command = [
        mmseqs,
        "easy-search",
        str(fasta),
        str(fasta),
        str(output),
        str(tmp_path / "mmseqs_tmp"),
        "--search-type",
        "3",
        "--threads",
        "1",
        "--strand",
        "2",
        "--min-seq-id",
        "1.0",
        "-c",
        "1.0",
        "--cov-mode",
        "2",
        "--alignment-mode",
        "3",
        "--mask",
        "0",
        "--max-seq-len",
        str(max_seq_len),
        "--max-seqs",
        "1000000",
        "--format-output",
        "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,qlen,tlen,qcov,tcov,qframe,tframe",
    ]
    subprocess.run(command, check=True, text=True, capture_output=True)  # noqa: S603


def candidate_pairs(path: Path) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            query, target, *_rest = line.rstrip("\n").split("\t")
            pairs.add((query, target))
    return pairs


def test_mmseqs_generates_exact_containment_candidates_on_both_strands(
    tmp_path: Path,
) -> None:
    forward_contained = "CGTACGTAGCTAGCTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT"
    reverse_source = "TTGGAACCGGTTCCAAGGTTAACCGGTTCCAAGGTTAACCGGTTCCAAGGTTAAC"
    reverse_contained = reverse_complement(reverse_source)
    one_mismatch = forward_contained[:-1] + "A"
    container = "A" * 40 + forward_contained + "T" * 40 + reverse_source + "G" * 40
    fasta = write_fasta(
        tmp_path / "representatives.fasta",
        [
            ("container", container),
            ("forward_contained", forward_contained),
            ("reverse_contained", reverse_contained),
            ("one_mismatch", one_mismatch),
        ],
    )
    output = tmp_path / "candidates.tsv"

    run_mmseqs_containment(
        fasta=fasta,
        output=output,
        tmp_path=tmp_path,
        max_seq_len=1000000,
    )

    pairs = candidate_pairs(output)
    assert ("forward_contained", "container") in pairs
    assert ("reverse_contained", "container") in pairs
    assert ("one_mismatch", "container") not in pairs


def test_mmseqs_max_seq_len_allows_split_boundary_crossing_candidates(
    tmp_path: Path,
) -> None:
    left = "A" * 9980
    contained = "CGTACGTAGCTAGCTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT"
    right = "T" * 6000
    container = left + contained + right
    fasta = write_fasta(
        tmp_path / "representatives.fasta",
        [
            ("container", container),
            ("crosses_boundary", contained),
        ],
    )
    output = tmp_path / "candidates.tsv"

    run_mmseqs_containment(
        fasta=fasta,
        output=output,
        tmp_path=tmp_path,
        max_seq_len=1000000,
    )

    assert ("crosses_boundary", "container") in candidate_pairs(output)
