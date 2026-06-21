"""Behavior tests for WVDB sourmash reference preparation."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_wvdb_sourmash_inputs.py"


def run_script(*args: str | Path) -> subprocess.CompletedProcess[str]:
    """Run the WVDB preparation script through its command-line interface."""
    return subprocess.run(  # noqa: S603
        [sys.executable, str(SCRIPT), *(str(arg) for arg in args)],
        check=False,
        cwd=ROOT,
        text=True,
        capture_output=True,
    )


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV rows as dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_tsv_rows(path: Path, rows: list[list[str]]) -> None:
    """Write TSV rows to a fixture file."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def write_wvdb_fixtures(tmp_path: Path) -> tuple[Path, Path]:
    """Write tiny WVDB-like FASTA and annotation fixtures."""
    fasta = tmp_path / "WVDB_v1.0.fasta"
    annotations = tmp_path / "WVDB_v1.0_annotations.tsv"

    fasta.write_text(
        """\
>contig_with_species original WVDB description
ACGTACGT
>contig_without_species
TGCATGCA
""",
        encoding="utf-8",
    )
    write_tsv_rows(
        annotations,
        [
            [
                "contig",
                "Phylum_gNd",
                "Class_gNd",
                "Order_gNd",
                "Family_gNd",
                "Phylum_RdRp",
                "Class_RdRp",
                "Order_RdRp",
                "Family_RdRp",
                "Genus_RdRp",
                "Species_RdRp",
                "Order_consensus",
                "Family_consensus",
                "Family_ICTV",
            ],
            [
                "contig_with_species",
                "Fallbackphylum",
                "Fallbackclass",
                "Fallbackorder",
                "Fallbackfamily",
                "RdRpphylum",
                "RdRpclass",
                "RdRporder",
                "RdRpfamily",
                "Knownvirus",
                "Knownvirus species",
                "Consensusorder",
                "Consensusfamily",
                "ICTVfamily",
            ],
            [
                "contig_without_species",
                "Fallbackonlyphylum",
                "Fallbackonlyclass",
                "Fallbackonlyorder",
                "Fallbackonlyfamily",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "ICTVfallbackfamily",
            ],
        ],
    )
    return fasta, annotations


def test_prepare_inputs_normalizes_fasta_and_maps_lineages(tmp_path: Path) -> None:
    """prepare-inputs writes normalized FASTA IDs and sourmash lineages."""
    fasta, annotations = write_wvdb_fixtures(tmp_path)
    normalized_fasta = tmp_path / "WVDB_v2.normalized.fasta"
    lineages_csv = tmp_path / "wvdb-v2.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--normalized-fasta",
        normalized_fasta,
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    assert normalized_fasta.read_text(encoding="utf-8") == (
        ">WVDB|contig_with_species\nACGTACGT\n>WVDB|contig_without_species\nTGCATGCA\n"
    )

    rows = read_csv_rows(lineages_csv)
    assert rows == [
        {
            "ident": "WVDB|contig_with_species",
            "superkingdom": "Viruses",
            "phylum": "RdRpphylum",
            "class": "RdRpclass",
            "order": "Consensusorder",
            "family": "Consensusfamily",
            "genus": "Knownvirus",
            "species": "Knownvirus species",
        },
        {
            "ident": "WVDB|contig_without_species",
            "superkingdom": "Viruses",
            "phylum": "Fallbackonlyphylum",
            "class": "Fallbackonlyclass",
            "order": "Fallbackonlyorder",
            "family": "Fallbackonlyfamily",
            "genus": "",
            "species": "",
        },
    ]


def test_manifest_coverage_ignores_non_wvdb_references(tmp_path: Path) -> None:
    """check-manifest-coverage succeeds when WVDB rows are covered."""
    manifest = tmp_path / "manifest.csv"
    lineages = tmp_path / "lineages.csv"

    manifest.write_text(
        "# SOURMASH-MANIFEST-VERSION: 1.0\n"
        "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n"
        "signatures/a.sig.gz,abc,abc,31,DNA,0,50,10,False,NC_005336.1,ref.fa\n"
        "signatures/b.sig.gz,def,def,31,DNA,0,50,10,False,WVDB|contig_1,ref.fa\n",
        encoding="utf-8",
    )
    lineages.write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species\n"
        "WVDB|contig_1,Viruses,,,,,,\n",
        encoding="utf-8",
    )

    result = run_script(
        "check-manifest-coverage",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 0, result.stderr


def test_manifest_coverage_fails_on_missing_wvdb_lineage(tmp_path: Path) -> None:
    """check-manifest-coverage fails when a WVDB signature lacks a lineage."""
    manifest = tmp_path / "manifest.csv"
    lineages = tmp_path / "lineages.csv"

    manifest.write_text(
        "# SOURMASH-MANIFEST-VERSION: 1.0\n"
        "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n"
        "signatures/b.sig.gz,def,def,31,DNA,0,50,10,False,WVDB|missing_contig,ref.fa\n",
        encoding="utf-8",
    )
    lineages.write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species\n"
        "WVDB|other_contig,Viruses,,,,,,\n",
        encoding="utf-8",
    )

    result = run_script(
        "check-manifest-coverage",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "WVDB signatures are missing lineages" in result.stderr
    assert "WVDB|missing_contig" in result.stderr
