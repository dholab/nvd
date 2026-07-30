#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# ///

"""Export CRUMBS taxon evidence in third-party taxonomic report formats."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence


REQUIRED_COLUMNS = [
    "sample_id",
    "taxon_id",
    "taxon_name",
    "rank",
    "taxpath",
    "taxpathsn",
    "rankpath",
    "percentage_emitted",
]
KRONA_RANKS = [
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
    "genome",
]
KREPORT_RANK_CODES = {
    "superkingdom": "D",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
}
DEFAULT_PROFILE_MASS_SCALE = 1_000_000_000


class CrumbsReportExportError(ValueError):
    """Raised when CRUMBS report exports cannot be written safely."""


@dataclass(frozen=True)
class TaxonProfileRow:
    taxon_id: str
    taxon_name: str
    rank: str
    taxpath: tuple[str, ...]
    taxpathsn: tuple[str, ...]
    rankpath: tuple[str, ...]
    percentage_emitted: float

    @property
    def fraction(self) -> float:
        return self.percentage_emitted / 100


@dataclass
class TreeNode:
    taxon_id: str
    taxon_name: str
    rank: str
    direct_fraction: float = 0.0
    children: dict[str, TreeNode] = field(default_factory=dict)

    @property
    def rank_code(self) -> str:
        return KREPORT_RANK_CODES.get(self.rank, "-")

    def clade_fraction(self) -> float:
        return self.direct_fraction + sum(
            child.clade_fraction() for child in self.children.values()
        )


def split_path(value: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in value.split("|") if part.strip())


def require_columns(fieldnames: Sequence[str] | None) -> None:
    if fieldnames is None:
        message = "CRUMBS taxa TSV is missing a header"
        raise CrumbsReportExportError(message)

    missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing:
        message = "CRUMBS taxa TSV missing required columns: " + ", ".join(missing)
        raise CrumbsReportExportError(message)


def parse_percentage(value: str, *, taxon_id: str) -> float:
    try:
        percentage = float(value)
    except ValueError as exc:
        message = f"taxon_id {taxon_id} has invalid percentage_emitted: {value}"
        raise CrumbsReportExportError(message) from exc

    if percentage < 0:
        message = f"taxon_id {taxon_id} has negative percentage_emitted"
        raise CrumbsReportExportError(message)
    return percentage


def parse_taxon_row(row: dict[str, str]) -> TaxonProfileRow:
    taxon_id = row["taxon_id"].strip()
    taxon_name = row["taxon_name"].strip()
    rank = row["rank"].strip()
    taxpath = split_path(row["taxpath"])
    taxpathsn = split_path(row["taxpathsn"])
    rankpath = split_path(row["rankpath"])

    if not taxon_id or not taxon_name or not rank:
        message = "CRUMBS taxa TSV contains blank taxon metadata"
        raise CrumbsReportExportError(message)
    if not taxpath or len(taxpath) != len(taxpathsn) or len(taxpath) != len(rankpath):
        message = f"taxon_id {taxon_id} has misaligned taxpath, taxpathsn, and rankpath"
        raise CrumbsReportExportError(message)
    if taxpath[-1] != taxon_id:
        message = f"taxon_id {taxon_id} taxpath does not end with taxon_id"
        raise CrumbsReportExportError(message)
    if taxpathsn[-1] != taxon_name:
        message = f"taxon_id {taxon_id} taxpathsn does not end with taxon_name"
        raise CrumbsReportExportError(message)
    if rankpath[-1] != rank:
        message = f"taxon_id {taxon_id} rankpath does not end with rank"
        raise CrumbsReportExportError(message)

    return TaxonProfileRow(
        taxon_id=taxon_id,
        taxon_name=taxon_name,
        rank=rank,
        taxpath=taxpath,
        taxpathsn=taxpathsn,
        rankpath=rankpath,
        percentage_emitted=parse_percentage(
            row["percentage_emitted"],
            taxon_id=taxon_id,
        ),
    )


def read_positive_taxa(path: Path, sample_id: str) -> list[TaxonProfileRow]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames)
        rows = [parse_taxon_row(row) for row in reader if row["sample_id"] == sample_id]

    return [row for row in rows if row.percentage_emitted > 0]


def krona_values(row: TaxonProfileRow) -> dict[str, str]:
    values = dict.fromkeys(KRONA_RANKS, "")
    for rank, name in zip(row.rankpath, row.taxpathsn, strict=True):
        if rank in values:
            values[rank] = name
    if not values["superkingdom"]:
        values["superkingdom"] = "unclassified"
    return values


def taxburst_safe_krona_ranks(rows: Sequence[TaxonProfileRow]) -> list[str]:
    if not rows:
        return KRONA_RANKS

    row_values = [krona_values(row) for row in rows]
    ranks: list[str] = []
    for rank in KRONA_RANKS:
        if all(values[rank] for values in row_values):
            ranks.append(rank)
            continue
        break
    return ranks


def krona_row(row: TaxonProfileRow, ranks: Sequence[str]) -> dict[str, str]:
    values = krona_values(row)
    selected_values = {rank: values[rank] for rank in ranks}
    return {"fraction": f"{row.fraction:.12g}", **selected_values}


def write_krona(rows: Sequence[TaxonProfileRow], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        ranks = taxburst_safe_krona_ranks(rows)
        fieldnames = ["fraction", *ranks]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(krona_row(row, ranks))


def add_to_tree(roots: dict[str, TreeNode], row: TaxonProfileRow) -> None:
    siblings = roots
    node: TreeNode | None = None
    for taxon_id, taxon_name, rank in zip(
        row.taxpath,
        row.taxpathsn,
        row.rankpath,
        strict=True,
    ):
        node = siblings.setdefault(
            taxon_id,
            TreeNode(taxon_id=taxon_id, taxon_name=taxon_name, rank=rank),
        )
        if node.taxon_name != taxon_name or node.rank != rank:
            message = f"taxon_id {taxon_id} has conflicting lineage metadata"
            raise CrumbsReportExportError(message)
        siblings = node.children

    if node is None:
        message = f"taxon_id {row.taxon_id} has an empty lineage"
        raise CrumbsReportExportError(message)
    node.direct_fraction += row.fraction


def build_tree(rows: Sequence[TaxonProfileRow]) -> dict[str, TreeNode]:
    roots: dict[str, TreeNode] = {}
    for row in rows:
        add_to_tree(roots, row)
    return roots


def scaled_mass(fraction: float, profile_mass_scale: int) -> int:
    return int(fraction * profile_mass_scale)


def kreport_lines(
    nodes: dict[str, TreeNode],
    *,
    profile_mass_scale: int,
    depth: int = 0,
) -> list[str]:
    lines: list[str] = []
    sorted_nodes = sorted(
        nodes.values(),
        key=lambda node: (-node.clade_fraction(), node.taxon_name, node.taxon_id),
    )
    for node in sorted_nodes:
        clade_fraction = node.clade_fraction()
        direct_fraction = node.direct_fraction
        indent = "  " * depth
        lines.append(
            "\t".join(
                [
                    f"{clade_fraction * 100:.2f}",
                    str(scaled_mass(clade_fraction, profile_mass_scale)),
                    str(scaled_mass(direct_fraction, profile_mass_scale)),
                    node.rank_code,
                    node.taxon_id,
                    f"{indent}{node.taxon_name}",
                ],
            ),
        )
        lines.extend(
            kreport_lines(
                node.children,
                profile_mass_scale=profile_mass_scale,
                depth=depth + 1,
            ),
        )
    return lines


def write_kreport(
    rows: Sequence[TaxonProfileRow],
    path: Path,
    *,
    profile_mass_scale: int,
) -> None:
    roots = build_tree(rows)
    text = "\n".join(
        kreport_lines(roots, profile_mass_scale=profile_mass_scale),
    )
    if text:
        text += "\n"
    path.write_text(text, encoding="utf-8")


@dataclass(frozen=True)
class OutputPaths:
    krona_tsv: Path
    kreport_txt: Path


def output_paths(sample_id: str, output_dir: Path) -> OutputPaths:
    return OutputPaths(
        krona_tsv=output_dir / f"{sample_id}.crumbs.krona.tsv",
        kreport_txt=output_dir / f"{sample_id}.crumbs.kreport.txt",
    )


def export_reports(
    *,
    sample_id: str,
    taxa_tsv: Path,
    output_dir: Path,
    profile_mass_scale: int,
) -> OutputPaths:
    if profile_mass_scale <= 0:
        message = "profile mass scale must be positive"
        raise CrumbsReportExportError(message)

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = output_paths(sample_id, output_dir)
    rows = read_positive_taxa(taxa_tsv, sample_id)
    write_krona(rows, paths.krona_tsv)
    write_kreport(rows, paths.kreport_txt, profile_mass_scale=profile_mass_scale)
    return paths


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export CRUMBS taxon evidence as Krona and Kraken-style reports",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--taxa-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--profile-mass-scale",
        type=int,
        default=DEFAULT_PROFILE_MASS_SCALE,
        help="Integer scale used for kreport count-like CRUMBS profile mass fields.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    export_reports(
        sample_id=args.sample_id,
        taxa_tsv=args.taxa_tsv,
        output_dir=args.output_dir,
        profile_mass_scale=args.profile_mass_scale,
    )


if __name__ == "__main__":
    main()
