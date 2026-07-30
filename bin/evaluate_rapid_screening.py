#!/usr/bin/env python3
"""Build rapid-screening follow-up evidence and eval artifacts."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import duckdb

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence


RANKS = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")
COMPARISON_RANKS = ("family", "genus", "species")
WVDB_ID_PREFIX = "WVDB|"
FOLLOWUP_BINS = (
    "same_rank_followup_evidence",
    "screening_reference_only",
    "no_same_rank_followup_evidence",
    "not_comparable",
)


class RapidScreeningEvalError(ValueError):
    """Raised when rapid-screening eval inputs are not usable."""


@dataclass(frozen=True)
class ScreeningDetection:
    sample_id: str
    screening_engine: str
    source_path: str
    rank: str
    screening_taxon_name: str
    screening_lineage: str
    primary_screening_mass: float
    sourmash_fraction: float | None
    sourmash_f_weighted_at_rank: float | None
    lineage_source: str
    followup_bin: str


@dataclass(frozen=True)
class NvdDetection:
    sample_id: str
    rank: str
    taxon_name: str
    taxon_id: str
    taxon_crumbs: float | None
    percentage_emitted: float | None


@dataclass(frozen=True)
class SampleSummary:
    sample_id: str
    total_reads: str
    has_sourmash_tax_summary: bool
    has_nvd_crumbs_taxa: bool
    sourmash_found_fraction: float | None
    sourmash_unclassified_fraction: float | None


def parse_float(value: object) -> float | None:
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def split_path(value: str, *, separator: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in value.split(separator) if part.strip())


def lineage_key_from_parts(
    rank: str,
    parts: Sequence[str],
) -> tuple[str, tuple[str, ...]] | None:
    if rank not in RANKS:
        return None
    depth = RANKS.index(rank) + 1
    if len(parts) < depth:
        return None
    path = tuple(parts[:depth])
    if any(not part for part in path):
        return None
    return rank, path


def normalize_lineage_identifier(row: dict[str, str]) -> str:
    return (row.get("ident") or row.get("identifiers") or "").strip()


def read_lineage_source_map(
    lineages_csv: Path,
) -> dict[tuple[str, tuple[str, ...]], frozenset[str]]:
    source_sets: defaultdict[tuple[str, tuple[str, ...]], set[str]] = defaultdict(set)
    with lineages_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            msg = f"{lineages_csv}: missing header"
            raise RapidScreeningEvalError(msg)
        for row in reader:
            ident = normalize_lineage_identifier(row)
            source = "wvdb" if ident.startswith(WVDB_ID_PREFIX) else "non_wvdb"
            parts = tuple((row.get(rank) or "").strip() for rank in RANKS)
            for rank in COMPARISON_RANKS:
                key = lineage_key_from_parts(rank, parts)
                if key is not None:
                    source_sets[key].add(source)
    return {key: frozenset(value) for key, value in source_sets.items()}


def read_sample_counts(path: Path) -> dict[str, str]:
    samples: dict[str, str] = {}
    if not path.exists() or path.stat().st_size == 0:
        return samples
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            sample_id = row[0].strip()
            if sample_id:
                samples[sample_id] = row[1].strip() if len(row) > 1 else ""
    return samples


def sample_id_from_sourmash_tax(path: Path) -> str:
    return path.name.split(".sourmash.tax_metagenome", maxsplit=1)[0]


def sample_id_from_gather(path: Path) -> str:
    return path.name.split(".sourmash.gather.csv", maxsplit=1)[0]


def sample_id_from_crumbs_taxa(path: Path) -> str:
    return path.name.split(".crumbs.taxa.tsv", maxsplit=1)[0]


def lineage_source_label(sources: frozenset[str] | None) -> str:
    if sources is None:
        return "missing_lineage_source"
    if sources == frozenset({"wvdb"}):
        return "wvdb"
    if sources == frozenset({"non_wvdb"}):
        return "non_wvdb"
    return "+".join(sorted(sources))


def classify_followup(
    *,
    source_set: frozenset[str] | None,
    has_followup_evidence: bool,
) -> str:
    if source_set is None:
        return "not_comparable"
    if has_followup_evidence:
        return "same_rank_followup_evidence"
    if source_set == frozenset({"wvdb"}):
        return "screening_reference_only"
    return "no_same_rank_followup_evidence"


def read_nvd_detections(crumbs_taxa_paths: Iterable[Path]) -> list[NvdDetection]:
    detections: dict[tuple[str, str, str], NvdDetection] = {}
    for path in crumbs_taxa_paths:
        sample_id = sample_id_from_crumbs_taxa(path)
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                rankpath = split_path(row.get("rankpath", ""), separator="|")
                taxpathsn = split_path(row.get("taxpathsn", ""), separator="|")
                if len(rankpath) != len(taxpathsn):
                    continue
                for rank, taxon_name in zip(rankpath, taxpathsn, strict=True):
                    if rank not in COMPARISON_RANKS or not taxon_name:
                        continue
                    key = (sample_id, rank, taxon_name)
                    previous = detections.get(key)
                    taxon_crumbs = parse_float(row.get("taxon_crumbs"))
                    percentage_emitted = parse_float(row.get("percentage_emitted"))
                    if previous is None:
                        detections[key] = NvdDetection(
                            sample_id=sample_id,
                            rank=rank,
                            taxon_name=taxon_name,
                            taxon_id=row.get("taxon_id", ""),
                            taxon_crumbs=taxon_crumbs,
                            percentage_emitted=percentage_emitted,
                        )
                    else:
                        detections[key] = NvdDetection(
                            sample_id=sample_id,
                            rank=rank,
                            taxon_name=taxon_name,
                            taxon_id=previous.taxon_id,
                            taxon_crumbs=sum_optional(
                                previous.taxon_crumbs,
                                taxon_crumbs,
                            ),
                            percentage_emitted=sum_optional(
                                previous.percentage_emitted,
                                percentage_emitted,
                            ),
                        )
    return list(detections.values())


def sum_optional(left: float | None, right: float | None) -> float | None:
    if left is None:
        return right
    if right is None:
        return left
    return left + right


def followup_evidence_index(
    detections: Iterable[NvdDetection],
) -> set[tuple[str, str, str]]:
    return {(row.sample_id, row.rank, row.taxon_name) for row in detections}


def primary_mass(
    row: dict[str, str],
) -> tuple[float | None, float | None, float | None]:
    fraction = parse_float(row.get("fraction"))
    f_weighted_at_rank = parse_float(row.get("f_weighted_at_rank"))
    primary = f_weighted_at_rank if f_weighted_at_rank is not None else fraction
    return primary, fraction, f_weighted_at_rank


def read_sourmash_detections(
    tax_summary_paths: Iterable[Path],
    *,
    lineage_sources: dict[tuple[str, tuple[str, ...]], frozenset[str]],
    followup_evidence: set[tuple[str, str, str]],
) -> list[ScreeningDetection]:
    detections: list[ScreeningDetection] = []
    for path in tax_summary_paths:
        sample_id = sample_id_from_sourmash_tax(path)
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rank = (row.get("rank") or "").strip()
                lineage = (row.get("lineage") or "").strip()
                if rank not in COMPARISON_RANKS or lineage == "unclassified":
                    continue
                primary, fraction, f_weighted_at_rank = primary_mass(row)
                if primary is None or primary <= 0:
                    continue
                parts = split_path(lineage, separator=";")
                key = lineage_key_from_parts(rank, parts)
                source_set = lineage_sources.get(key) if key is not None else None
                taxon_name = parts[RANKS.index(rank)] if key is not None else ""
                has_followup_evidence = (
                    sample_id,
                    rank,
                    taxon_name,
                ) in followup_evidence
                followup_bin = classify_followup(
                    source_set=source_set,
                    has_followup_evidence=has_followup_evidence,
                )
                detections.append(
                    ScreeningDetection(
                        sample_id=sample_id,
                        screening_engine="sourmash",
                        source_path=str(path),
                        rank=rank,
                        screening_taxon_name=taxon_name,
                        screening_lineage=lineage,
                        primary_screening_mass=primary,
                        sourmash_fraction=fraction,
                        sourmash_f_weighted_at_rank=f_weighted_at_rank,
                        lineage_source=lineage_source_label(source_set),
                        followup_bin=followup_bin,
                    ),
                )
    return detections


def gather_summaries(
    gather_paths: Iterable[Path],
) -> dict[str, tuple[float | None, float | None]]:
    summaries: dict[str, tuple[float | None, float | None]] = {}
    for path in gather_paths:
        sample_id = sample_id_from_gather(path)
        found_values: list[float] = []
        weighted_found: list[float] = []
        total_weighted: list[float] = []
        if path.exists() and path.stat().st_size > 0:
            with path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    value = parse_float(row.get("f_unique_to_query"))
                    if value is not None:
                        found_values.append(value)
                    weighted_value = parse_float(row.get("sum_weighted_found"))
                    total_value = parse_float(row.get("total_weighted_hashes"))
                    if weighted_value is not None:
                        weighted_found.append(weighted_value)
                    if total_value is not None:
                        total_weighted.append(total_value)
        found_fraction: float | None = None
        if weighted_found and total_weighted and max(total_weighted) > 0:
            found_fraction = max(weighted_found) / max(total_weighted)
        elif found_values:
            found_fraction = min(sum(found_values), 1.0)
        unclassified = (
            None if found_fraction is None else max(0.0, 1.0 - found_fraction)
        )
        summaries[sample_id] = (found_fraction, unclassified)
    return summaries


def build_sample_summaries(
    *,
    sample_counts: dict[str, str],
    tax_summary_paths: Sequence[Path],
    crumbs_taxa_paths: Sequence[Path],
    gather_paths: Sequence[Path],
) -> list[SampleSummary]:
    tax_samples = {sample_id_from_sourmash_tax(path) for path in tax_summary_paths}
    crumbs_samples = {sample_id_from_crumbs_taxa(path) for path in crumbs_taxa_paths}
    gather_by_sample = gather_summaries(gather_paths)
    sample_ids = sorted(
        set(sample_counts) | tax_samples | crumbs_samples | set(gather_by_sample),
    )
    return [
        SampleSummary(
            sample_id=sample_id,
            total_reads=sample_counts.get(sample_id, ""),
            has_sourmash_tax_summary=sample_id in tax_samples,
            has_nvd_crumbs_taxa=sample_id in crumbs_samples,
            sourmash_found_fraction=gather_by_sample.get(sample_id, (None, None))[0],
            sourmash_unclassified_fraction=gather_by_sample.get(
                sample_id,
                (None, None),
            )[1],
        )
        for sample_id in sample_ids
    ]


def followup_summary_rows(
    detections: Iterable[ScreeningDetection],
) -> list[dict[str, object]]:
    buckets: defaultdict[tuple[str, str, str, str], dict[str, object]] = defaultdict(
        lambda: {"n_signals": 0, "total_primary_screening_mass": 0.0},
    )
    for detection in detections:
        key = (
            detection.sample_id,
            detection.screening_engine,
            detection.rank,
            detection.followup_bin,
        )
        buckets[key]["n_signals"] = int(buckets[key]["n_signals"]) + 1
        buckets[key]["total_primary_screening_mass"] = (
            float(buckets[key]["total_primary_screening_mass"])
            + detection.primary_screening_mass
        )
    rows = []
    for (sample_id, screening_engine, rank, followup_bin), values in sorted(
        buckets.items(),
    ):
        rows.append(
            {
                "sample_id": sample_id,
                "screening_engine": screening_engine,
                "comparison_rank": rank,
                "followup_bin": followup_bin,
                "n_signals": values["n_signals"],
                "total_primary_screening_mass": values["total_primary_screening_mass"],
            },
        )
    return rows


def without_followup_rows(
    detections: Iterable[ScreeningDetection],
) -> list[dict[str, object]]:
    rows = [
        {
            "sample_id": detection.sample_id,
            "screening_engine": detection.screening_engine,
            "comparison_rank": detection.rank,
            "screening_taxon_name": detection.screening_taxon_name,
            "screening_lineage": detection.screening_lineage,
            "primary_screening_mass": detection.primary_screening_mass,
            "sourmash_fraction": detection.sourmash_fraction,
            "sourmash_f_weighted_at_rank": detection.sourmash_f_weighted_at_rank,
            "lineage_source": detection.lineage_source,
            "source_path": detection.source_path,
        }
        for detection in detections
        if detection.followup_bin == "no_same_rank_followup_evidence"
    ]
    return sorted(
        rows,
        key=lambda row: float(row["primary_screening_mass"]),
        reverse=True,
    )


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_database(
    path: Path,
    *,
    samples: Sequence[SampleSummary],
    screening_detections: Sequence[ScreeningDetection],
    nvd_detections: Sequence[NvdDetection],
    followup_rows: Sequence[dict[str, object]],
    without_followup: Sequence[dict[str, object]],
) -> None:
    if path.exists():
        path.unlink()
    con = duckdb.connect(str(path))
    try:
        con.execute(
            """
            create table samples (
              sample_id varchar,
              total_reads varchar,
              has_sourmash_tax_summary boolean,
              has_nvd_crumbs_taxa boolean,
              sourmash_found_fraction double,
              sourmash_unclassified_fraction double
            )
            """,
        )
        sample_values = [
            (
                row.sample_id,
                row.total_reads,
                row.has_sourmash_tax_summary,
                row.has_nvd_crumbs_taxa,
                row.sourmash_found_fraction,
                row.sourmash_unclassified_fraction,
            )
            for row in samples
        ]
        if sample_values:
            con.executemany(
                "insert into samples values (?, ?, ?, ?, ?, ?)",
                sample_values,
            )
        con.execute(
            """
            create table screening_detections (
              sample_id varchar,
              screening_engine varchar,
              source_path varchar,
              rank varchar,
              screening_taxon_name varchar,
              screening_lineage varchar,
              primary_screening_mass double,
              sourmash_fraction double,
              sourmash_f_weighted_at_rank double,
              lineage_source varchar,
              followup_bin varchar
            )
            """,
        )
        screening_values = [
            (
                row.sample_id,
                row.screening_engine,
                row.source_path,
                row.rank,
                row.screening_taxon_name,
                row.screening_lineage,
                row.primary_screening_mass,
                row.sourmash_fraction,
                row.sourmash_f_weighted_at_rank,
                row.lineage_source,
                row.followup_bin,
            )
            for row in screening_detections
        ]
        if screening_values:
            con.executemany(
                "insert into screening_detections values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                screening_values,
            )
        con.execute(
            """
            create table nvd_detections (
              sample_id varchar,
              rank varchar,
              taxon_name varchar,
              taxon_id varchar,
              taxon_crumbs double,
              percentage_emitted double
            )
            """,
        )
        nvd_values = [
            (
                row.sample_id,
                row.rank,
                row.taxon_name,
                row.taxon_id,
                row.taxon_crumbs,
                row.percentage_emitted,
            )
            for row in nvd_detections
        ]
        if nvd_values:
            con.executemany(
                "insert into nvd_detections values (?, ?, ?, ?, ?, ?)",
                nvd_values,
            )
        write_dict_table(con, "screening_signal_followup_by_rank", followup_rows)
        write_dict_table(
            con,
            "screening_signals_without_same_rank_followup",
            without_followup,
        )
    finally:
        con.close()


def write_dict_table(
    con: duckdb.DuckDBPyConnection,
    table: str,
    rows: Sequence[dict[str, object]],
) -> None:
    if not rows:
        if table == "screening_signal_followup_by_rank":
            con.execute(
                """
                create table screening_signal_followup_by_rank (
                  sample_id varchar,
                  screening_engine varchar,
                  comparison_rank varchar,
                  followup_bin varchar,
                  n_signals bigint,
                  total_primary_screening_mass double
                )
                """,
            )
        else:
            con.execute(
                """
                create table screening_signals_without_same_rank_followup (
                  sample_id varchar,
                  screening_engine varchar,
                  comparison_rank varchar,
                  screening_taxon_name varchar,
                  screening_lineage varchar,
                  primary_screening_mass double,
                  sourmash_fraction double,
                  sourmash_f_weighted_at_rank double,
                  lineage_source varchar,
                  source_path varchar
                )
                """,
            )
        return
    columns = list(rows[0])
    column_sql = ", ".join(f"{column} varchar" for column in columns)
    con.execute(f"create table {table} ({column_sql})")
    placeholders = ", ".join("?" for _ in columns)
    con.executemany(
        f"insert into {table} values ({placeholders})",
        [tuple(row.get(column) for column in columns) for row in rows],
    )


def build_outputs(
    *,
    sample_read_counts: Path,
    lineages_csv: Path,
    sourmash_tax_summaries: Sequence[Path],
    sourmash_gather_csvs: Sequence[Path],
    crumbs_taxa_tsvs: Sequence[Path],
    output_database: Path,
    output_dir: Path,
) -> None:
    sample_counts = read_sample_counts(sample_read_counts)
    lineage_sources = read_lineage_source_map(lineages_csv)
    nvd_detections = read_nvd_detections(crumbs_taxa_tsvs)
    screening_detections = read_sourmash_detections(
        sourmash_tax_summaries,
        lineage_sources=lineage_sources,
        followup_evidence=followup_evidence_index(nvd_detections),
    )
    samples = build_sample_summaries(
        sample_counts=sample_counts,
        tax_summary_paths=sourmash_tax_summaries,
        crumbs_taxa_paths=crumbs_taxa_tsvs,
        gather_paths=sourmash_gather_csvs,
    )
    followup = followup_summary_rows(screening_detections)
    without_followup = without_followup_rows(screening_detections)
    exports_dir = output_dir / "exports"
    write_tsv(
        exports_dir / "screening_signal_followup_by_sample_rank.tsv",
        followup,
        [
            "sample_id",
            "screening_engine",
            "comparison_rank",
            "followup_bin",
            "n_signals",
            "total_primary_screening_mass",
        ],
    )
    write_tsv(
        exports_dir / "screening_signals_without_same_rank_followup.tsv",
        without_followup,
        [
            "sample_id",
            "screening_engine",
            "comparison_rank",
            "screening_taxon_name",
            "screening_lineage",
            "primary_screening_mass",
            "sourmash_fraction",
            "sourmash_f_weighted_at_rank",
            "lineage_source",
            "source_path",
        ],
    )
    write_database(
        output_database,
        samples=samples,
        screening_detections=screening_detections,
        nvd_detections=nvd_detections,
        followup_rows=followup,
        without_followup=without_followup,
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build rapid-screening eval artifacts from sourmash and NVD evidence",
    )
    parser.add_argument("--sample-read-counts", type=Path, required=True)
    parser.add_argument("--lineages-csv", type=Path, required=True)
    parser.add_argument(
        "--sourmash-tax-summary",
        type=Path,
        action="append",
        default=[],
    )
    parser.add_argument("--sourmash-gather-csv", type=Path, action="append", default=[])
    parser.add_argument("--blast-tsv", type=Path, action="append", default=[])
    parser.add_argument("--crumbs-taxa-tsv", type=Path, action="append", default=[])
    parser.add_argument("--crumbs-contigs-tsv", type=Path, action="append", default=[])
    parser.add_argument("--output-database", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    _ = args.blast_tsv
    _ = args.crumbs_contigs_tsv
    build_outputs(
        sample_read_counts=args.sample_read_counts,
        lineages_csv=args.lineages_csv,
        sourmash_tax_summaries=args.sourmash_tax_summary,
        sourmash_gather_csvs=args.sourmash_gather_csv,
        crumbs_taxa_tsvs=args.crumbs_taxa_tsv,
        output_database=args.output_database,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
