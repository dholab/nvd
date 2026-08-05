#!/usr/bin/env python3
"""Normalize stage-specific routing evidence into one experiment ledger."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


OUTPUT_FIELDS = (
    "sample_id",
    "stage",
    "input_class",
    "output_class",
    "decision",
    "reason",
    "sequences_in",
    "sequences_out",
    "bases_in",
    "bases_out",
)
MAPBACK_FIELD_COUNT = 2


@dataclass(frozen=True)
class LedgerRow:
    sample_id: str
    stage: str
    input_class: str
    output_class: str
    decision: str
    reason: str
    sequences_in: object = ""
    sequences_out: object = ""
    bases_in: object = ""
    bases_out: object = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def eligibility_rows(path: Path) -> list[LedgerRow]:
    return [
        LedgerRow(
            sample_id=source["sample_id"],
            stage="long_read_assembly_eligibility",
            input_class="long_reads",
            output_class=assembler,
            decision=source[f"{assembler}_decision"],
            reason=source[f"{assembler}_reason"],
            sequences_in=source["sequence_count"],
            sequences_out=source[f"{assembler}_qualifying_reads"],
            bases_in=source["total_bases"],
            bases_out=source[f"{assembler}_qualifying_bases"],
        )
        for source in read_tsv(path)
        for assembler in ("metamdbg", "myloasm", "metaflye")
    ]


def union_summary_rows(path: Path) -> list[LedgerRow]:
    summary = json.loads(path.read_text(encoding="utf-8"))
    return [
        LedgerRow(
            sample_id=summary["sample_id"],
            stage="long_read_contig_union",
            input_class="assembler_contigs",
            output_class="union_contigs",
            decision="collapse",
            reason="exact_duplicate_and_containment_collapse",
            sequences_in=summary["source_contig_count"],
            sequences_out=summary["emitted_contig_count"],
        ),
    ]


def contig_filter_rows(path: Path) -> list[LedgerRow]:
    rows = []
    for source in read_tsv(path):
        phase = source["phase"]
        rows.append(
            LedgerRow(
                sample_id=source["sample_id"],
                stage=f"contig_{phase}",
                input_class="contigs",
                output_class="retained_contigs",
                decision="filter",
                reason=f"passed_{phase}",
                sequences_in=source["sequences_in"],
                sequences_out=source["sequences_out"],
                bases_in=source["bases_in"],
                bases_out=source["bases_out"],
            ),
        )
    return rows


def mapback_rows(path: Path) -> list[LedgerRow]:
    mapped_segments = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != MAPBACK_FIELD_COUNT:
                message = f"Unexpected mapback-count row in {path}: {line!r}"
                raise ValueError(message)
            mapped_segments += int(fields[1])
    return [
        LedgerRow(
            sample_id=path.name.removesuffix("_mapped_counts.txt"),
            stage="contig_read_mapback",
            input_class="preprocessed_reads",
            output_class="contig_mapped_read_segments",
            decision="map",
            reason="aligned_to_retained_contigs",
            sequences_out=mapped_segments,
        ),
    ]


def query_batch_rows(path: Path) -> list[LedgerRow]:
    rows = []
    for source in read_tsv(path):
        count = int(source["n_query_sequences"])
        rows.append(
            LedgerRow(
                sample_id=source["sample_id"],
                stage="blast_query_preparation",
                input_class=source["query_source"],
                output_class=source["query_class"],
                decision="emit" if count else "skip",
                reason="query_batch_emitted" if count else "no_query_sequences",
                sequences_out=count,
            ),
        )
    return rows


def megablast_partition_rows(path: Path) -> list[LedgerRow]:
    rows = []
    for source in read_tsv(path):
        common = {
            "sample_id": source["sample_id"],
            "stage": "megablast_query_partition",
            "input_class": source["query_class"],
            "sequences_in": source["queries_in"],
        }
        rows.extend(
            (
                LedgerRow(
                    **common,
                    output_class="megablast_accounted",
                    decision="account",
                    reason="megablast_hit_present",
                    sequences_out=source["megablast_accounted"],
                ),
                LedgerRow(
                    **common,
                    output_class="blastn_candidates",
                    decision="dispatch",
                    reason="no_megablast_hit",
                    sequences_out=source["blastn_candidates"],
                ),
            ),
        )
    return rows


def blast_filter_rows(path: Path) -> list[LedgerRow]:
    return [
        LedgerRow(
            sample_id=source["sample_id"],
            stage=source["stage"],
            input_class=source["query_class"],
            output_class="virus_only_hits",
            decision="retain",
            reason="viral_taxonomy_match",
            sequences_in=source["queries_in"],
            sequences_out=source["queries_retained"],
        )
        for source in read_tsv(path)
    ]


def rows_for_evidence(path: Path) -> list[LedgerRow]:
    name = path.name
    adapters = (
        ((".long_read_assembly_eligibility.tsv",), eligibility_rows),
        ((".long_read_union.summary.json",), union_summary_rows),
        ((".contig_target_filtering.tsv",), contig_filter_rows),
        (("_mapped_counts.txt",), mapback_rows),
        ((".blast_query_batches.tsv",), query_batch_rows),
        ((".megablast_query_partition.tsv",), megablast_partition_rows),
        (
            (".megablast_query_filtering.tsv", ".blastn_query_filtering.tsv"),
            blast_filter_rows,
        ),
    )
    for suffixes, adapter in adapters:
        if name.endswith(suffixes):
            return adapter(path)
    message = f"Unsupported sequence-flow evidence file: {path}"
    raise ValueError(message)


def build_rows(paths: Iterable[Path]) -> list[LedgerRow]:
    rows = [row for path in paths for row in rows_for_evidence(path)]
    return sorted(
        rows,
        key=lambda row: (
            row.sample_id,
            row.stage,
            row.input_class,
            row.output_class,
        ),
    )


def write_rows(path: Path, rows: Iterable[LedgerRow]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def main() -> None:
    args = parse_args()
    write_rows(args.output, build_rows(args.evidence))


if __name__ == "__main__":
    main()
