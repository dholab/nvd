#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "loguru",
# ]
# ///

"""Prepare an explicit CRUMBS profile taxonomy table from NCBI taxids."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from py_nvd import taxonomy

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import TextIO

    from py_nvd.models import Taxon


OUTPUT_COLUMNS = ["taxon_id", "taxon_name", "rank", "taxpath", "taxpathsn"]


class ProfileTaxonomyError(ValueError):
    """Raised when profile taxonomy rows cannot be prepared safely."""


@dataclass(frozen=True)
class ProfileTaxonRow:
    taxon_id: str
    taxon_name: str
    rank: str
    taxpath: str
    taxpathsn: str

    def as_dict(self) -> dict[str, str]:
        return {
            "taxon_id": self.taxon_id,
            "taxon_name": self.taxon_name,
            "rank": self.rank,
            "taxpath": self.taxpath,
            "taxpathsn": self.taxpathsn,
        }


def read_taxids(path: Path) -> list[int]:
    """Read observed NCBI taxids from a one-column text file."""
    taxids: list[int] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                taxid = int(stripped)
            except ValueError as exc:
                message = f"line {line_number}: expected an integer taxid"
                raise ProfileTaxonomyError(message) from exc
            if taxid <= 0:
                message = f"line {line_number}: taxid must be positive"
                raise ProfileTaxonomyError(message)
            taxids.append(taxid)
    return taxids


def unique_preserving_order(values: Iterable[int]) -> list[int]:
    """Return unique values in first-seen order."""
    seen: set[int] = set()
    unique: list[int] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        unique.append(value)
    return unique


def profile_taxon_row(tax: taxonomy.TaxonomyDB, observed_taxid: int) -> ProfileTaxonRow:
    """Resolve one observed NCBI taxid into a profile taxonomy row."""
    taxon = tax.get_taxon(observed_taxid)
    if taxon is None:
        message = f"taxid {observed_taxid} was not found in the taxonomy database"
        raise ProfileTaxonomyError(message)

    lineage = tax.get_lineage(taxon.tax_id)
    if not lineage or lineage[-1].tax_id != taxon.tax_id:
        lineage = [*lineage, taxon]

    return build_profile_taxon_row(taxon, lineage)


def build_profile_taxon_row(
    taxon: Taxon,
    lineage: Sequence[Taxon],
) -> ProfileTaxonRow:
    """Build and validate a profile taxonomy row from resolved lineage."""
    taxpath_parts = [str(lineage_taxon.tax_id) for lineage_taxon in lineage]
    taxpathsn_parts = [lineage_taxon.scientific_name for lineage_taxon in lineage]
    taxon_id = str(taxon.tax_id)
    taxon_name = taxon.scientific_name

    if not taxpath_parts or taxpath_parts[-1] != taxon_id:
        message = f"lineage for taxid {taxon_id} does not end with the taxon id"
        raise ProfileTaxonomyError(message)
    if not taxpathsn_parts or taxpathsn_parts[-1] != taxon_name:
        message = f"lineage for taxid {taxon_id} does not end with the taxon name"
        raise ProfileTaxonomyError(message)
    if len(taxpath_parts) != len(taxpathsn_parts):
        message = f"lineage ids and names differ in length for taxid {taxon_id}"
        raise ProfileTaxonomyError(message)

    return ProfileTaxonRow(
        taxon_id=taxon_id,
        taxon_name=taxon_name,
        rank=taxon.rank,
        taxpath="|".join(taxpath_parts),
        taxpathsn="|".join(taxpathsn_parts),
    )


def prepare_profile_taxonomy(
    taxids: Iterable[int],
    *,
    taxonomy_dir: Path,
    taxonomy_mode: str,
    taxonomy_max_age_days: int | None,
) -> list[ProfileTaxonRow]:
    """Resolve observed taxids through NVD's NCBI taxonomy cache."""
    rows_by_taxon_id: dict[str, ProfileTaxonRow] = {}
    with taxonomy.open(
        taxonomy_dir=taxonomy_dir,
        taxonomy_mode=taxonomy_mode,
        max_age_days=taxonomy_max_age_days,
    ) as tax:
        for observed_taxid in unique_preserving_order(taxids):
            row = profile_taxon_row(tax, observed_taxid)
            rows_by_taxon_id.setdefault(row.taxon_id, row)
    return list(rows_by_taxon_id.values())


def write_profile_taxonomy(rows: Iterable[ProfileTaxonRow], output: TextIO) -> None:
    """Write profile taxonomy rows as a TSV."""
    writer = csv.DictWriter(output, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(row.as_dict() for row in rows)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare an explicit CRUMBS profile taxonomy TSV from NCBI taxids",
    )
    parser.add_argument("--taxids", type=Path, required=True)
    parser.add_argument("--taxonomy-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--taxonomy-mode",
        default="read_only",
        choices=[mode.value for mode in taxonomy.TaxonomyMode],
    )
    parser.add_argument("--taxonomy-max-age-days", type=int, default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    try:
        rows = prepare_profile_taxonomy(
            read_taxids(args.taxids),
            taxonomy_dir=args.taxonomy_dir,
            taxonomy_mode=args.taxonomy_mode,
            taxonomy_max_age_days=args.taxonomy_max_age_days,
        )
        with args.output.open("w", newline="", encoding="utf-8") as output:
            write_profile_taxonomy(rows, output)
    except (
        OSError,
        ProfileTaxonomyError,
        taxonomy.TaxonomyBuildError,
        taxonomy.TaxonomyOfflineError,
        taxonomy.TaxonomyUnavailableError,
    ) as exc:
        sys.stderr.write(f"{exc}\n")
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
