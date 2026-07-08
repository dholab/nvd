"""Behavior tests for WVDB sourmash reference preparation."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_wvdb_sourmash_inputs.py"
MANIFEST_HEADER = (
    "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,"
    "with_abundance,name,filename"
)
LINEAGE_HEADER = "ident,superkingdom,phylum,class,order,family,genus,species,strain"
WVDB_LINEAGE_COLUMNS = (
    "ident",
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


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


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read diagnostic TSV rows as dictionaries."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv_rows_by_ident(path: Path) -> dict[str, dict[str, str]]:
    """Read CSV rows keyed by sourmash lineage identifier."""
    return {row["ident"]: row for row in read_csv_rows(path)}


def wvdb_lineage_columns(row: dict[str, str]) -> dict[str, str]:
    """Return lineage-name columns, excluding BioBoxes taxpath metadata."""
    return {column: row[column] for column in WVDB_LINEAGE_COLUMNS}


def write_reference_lineages(path: Path) -> Path:
    """Write a tiny reference lineage CSV with sourmash taxpaths."""
    rows = [
        [
            "ident",
            "taxid",
            "taxpath",
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "strain",
        ],
        [
            "REF_PICORNA",
            "1234567",
            "10239|2732408|2732506|464095|12345|23456|34567|",
            "Viruses",
            "Pisuviricota",
            "Pisoniviricetes",
            "Picornavirales",
            "Picornaviridae",
            "Referencegenus",
            "Reference species",
            "",
        ],
        [
            "REF_KNOWN",
            "8888888",
            "10239|111111|222222|333333|444444|555555|666666|",
            "Viruses",
            "RdRpphylum",
            "RdRpclass",
            "Consensusorder",
            "Consensusfamily",
            "Knownvirus",
            "Knownvirus species",
            "",
        ],
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows(rows)
    return path


def write_tsv_rows(path: Path, rows: list[list[str]]) -> None:
    """Write TSV rows to a fixture file."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def write_manifest(
    path: Path,
    rows: list[tuple[str, str]],
    *,
    scaled: str = "50",
) -> Path:
    """Write a minimal sourmash manifest from identifier/digest pairs."""
    lines = ["# SOURMASH-MANIFEST-VERSION: 1.0", MANIFEST_HEADER]
    lines.extend(
        f"signatures/{md5}.sig.gz,{md5},{md5[:8]},31,DNA,0,{scaled},10,False,{ident},ref.fa"
        for ident, md5 in rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_lineages(path: Path, rows: list[tuple[str, ...]]) -> Path:
    """Write minimal lineage rows for curation tests."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(LINEAGE_HEADER.split(","))
        writer.writerows(rows)
    return path


def run_curation(
    tmp_path: Path,
    *,
    manifests: tuple[Path, Path],
    lineages: tuple[Path, Path],
) -> tuple[subprocess.CompletedProcess[str], Path, Path]:
    """Run lineage curation and return the process and output paths."""
    output = tmp_path / "combined.lineages.csv"
    diagnostics = tmp_path / "taxonomy-conflicts.tsv"
    result = run_script(
        "curate-lineages",
        "--reference-manifest-csv",
        manifests[0],
        "--wvdb-manifest-csv",
        manifests[1],
        "--reference-lineages-csv",
        lineages[0],
        "--wvdb-lineages-csv",
        lineages[1],
        "--lineages-csv",
        output,
        "--diagnostics-tsv",
        diagnostics,
    )
    return result, output, diagnostics


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
>contig_species_without_genus
AACCGGTT
>contig_order_without_phylum_or_class
CCAATTGG
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
            [
                "contig_species_without_genus",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Unclassified",
                "",
                "",
                "",
                "",
                "",
                "Wenzhou_picorna-like_virus_47",
                "",
                "",
                "",
            ],
            [
                "contig_order_without_phylum_or_class",
                "",
                "",
                "Picornavirales",
                "Unclassified",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
            ],
        ],
    )
    return fasta, annotations


def test_prepare_inputs_normalizes_fasta_and_maps_lineages(tmp_path: Path) -> None:
    """prepare-inputs writes normalized FASTA IDs and sourmash lineages."""
    fasta, annotations = write_wvdb_fixtures(tmp_path)
    normalized_fasta = tmp_path / "WVDB_v1.0.normalized.fasta"
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

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
        ">WVDB|contig_with_species\nACGTACGT\n"
        ">WVDB|contig_without_species\nTGCATGCA\n"
        ">WVDB|contig_species_without_genus\nAACCGGTT\n"
        ">WVDB|contig_order_without_phylum_or_class\nCCAATTGG\n"
    )

    rows = read_csv_rows(lineages_csv)
    assert [wvdb_lineage_columns(row) for row in rows] == [
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
        {
            "ident": "WVDB|contig_species_without_genus",
            "superkingdom": "Viruses",
            "phylum": "Pisuviricota",
            "class": "Pisoniviricetes",
            "order": "Picornavirales",
            "family": "unclassified family [under Pisoniviricetes > Picornavirales]",
            "genus": "unclassified genus [under Picornavirales]",
            "species": "Wenzhou_picorna-like_virus_47",
        },
        {
            "ident": "WVDB|contig_order_without_phylum_or_class",
            "superkingdom": "Viruses",
            "phylum": "unclassified phylum [under Viruses]",
            "class": "unclassified class [under Viruses]",
            "order": "Picornavirales",
            "family": "unclassified family [under unclassified class > Picornavirales]",
            "genus": "",
            "species": "",
        },
    ]
    assert all(row["taxpath"].startswith("10239|") for row in rows)


def test_prepare_inputs_contextualizes_wvdb_placeholder_taxa(tmp_path: Path) -> None:
    """WVDB placeholder ranks should be named by stable lineage context."""
    fasta, annotations = write_wvdb_fixtures(tmp_path)
    normalized_fasta = tmp_path / "WVDB_v1.0.normalized.fasta"
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

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
    rows = read_csv_rows_by_ident(lineages_csv)

    assert rows["WVDB|contig_species_without_genus"]["family"] == (
        "unclassified family [under Pisoniviricetes > Picornavirales]"
    )
    assert rows["WVDB|contig_species_without_genus"]["genus"] == (
        "unclassified genus [under Picornavirales]"
    )
    assert rows["WVDB|contig_order_without_phylum_or_class"]["phylum"] == (
        "unclassified phylum [under Viruses]"
    )
    assert rows["WVDB|contig_order_without_phylum_or_class"]["class"] == (
        "unclassified class [under Viruses]"
    )
    assert rows["WVDB|contig_order_without_phylum_or_class"]["family"] == (
        "unclassified family [under unclassified class > Picornavirales]"
    )


def test_prepare_inputs_selects_one_coherent_ancestry_source(tmp_path: Path) -> None:
    """A deeper gNd path must not inherit an incompatible isolated RdRp phylum."""
    fasta = tmp_path / "wvdb.fasta"
    fasta.write_text(">contig_1\nACGT\n", encoding="utf-8")
    annotations = tmp_path / "annotations.tsv"
    write_tsv_rows(
        annotations,
        [
            [
                "contig",
                "Phylum_RdRp",
                "Class_RdRp",
                "Order_RdRp",
                "Family_RdRp",
                "Phylum_gNd",
                "Class_gNd",
                "Order_gNd",
                "Family_gNd",
                "Order_consensus",
                "Family_consensus",
            ],
            [
                "contig_1",
                "Pisuviricota",
                "",
                "",
                "",
                "Kitrinoviricota",
                "Alsuviricetes",
                "Tymovirales",
                "Betaflexiviridae",
                "Tymovirales",
                "Betaflexiviridae",
            ],
        ],
    )
    lineages = tmp_path / "lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--normalized-fasta",
        tmp_path / "normalized.fasta",
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 0, result.stderr
    row = read_csv_rows(lineages)[0]
    assert row["phylum"] == "Kitrinoviricota"
    assert row["class"] == "Alsuviricetes"
    assert row["order"] == "Tymovirales"
    assert row["family"] == "Betaflexiviridae"


def test_prepare_inputs_reuses_reference_taxpaths_for_exact_prefixes(
    tmp_path: Path,
) -> None:
    """WVDB lineages should reuse exact reference taxids where possible."""
    fasta, annotations = write_wvdb_fixtures(tmp_path)
    reference_lineages = write_reference_lineages(tmp_path / "reference.lineages.csv")
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--reference-lineages-csv",
        reference_lineages,
        "--normalized-fasta",
        tmp_path / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(lineages_csv)
    knownvirus_taxpath = rows["WVDB|contig_with_species"]["taxpath"].split("|")
    assert knownvirus_taxpath == [
        "10239",
        "111111",
        "222222",
        "333333",
        "444444",
        "555555",
        "666666",
    ]


def test_prepare_inputs_backfills_unique_reference_ancestors(
    tmp_path: Path,
) -> None:
    """Known descendants should backfill missing ancestors from reference lineages."""
    fasta, annotations = write_wvdb_fixtures(tmp_path)
    reference_lineages = write_reference_lineages(tmp_path / "reference.lineages.csv")
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--reference-lineages-csv",
        reference_lineages,
        "--normalized-fasta",
        tmp_path / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    row = read_csv_rows_by_ident(lineages_csv)[
        "WVDB|contig_order_without_phylum_or_class"
    ]
    assert row["phylum"] == "Pisuviricota"
    assert row["class"] == "Pisoniviricetes"
    assert row["order"] == "Picornavirales"
    assert row["taxpath"].split("|")[:3] == ["10239", "2732408", "2732506"]


def test_prepare_inputs_backfills_from_shallower_known_rank_for_novel_descendant(
    tmp_path: Path,
) -> None:
    """A novel low rank should not block ancestor backfill from a known order."""
    fasta = tmp_path / "WVDB_v1.0.fasta"
    annotations = tmp_path / "WVDB_v1.0_annotations.tsv"
    fasta.write_text(
        ">novel_picorna_like\nACGTACGT\n",
        encoding="utf-8",
    )
    write_tsv_rows(
        annotations,
        [
            ["contig", "Order_gNd", "Species_RdRp"],
            ["novel_picorna_like", "Picornavirales", "Novel picorna-like virus"],
        ],
    )
    reference_lineages = write_reference_lineages(tmp_path / "reference.lineages.csv")
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--reference-lineages-csv",
        reference_lineages,
        "--normalized-fasta",
        tmp_path / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    row = read_csv_rows_by_ident(lineages_csv)["WVDB|novel_picorna_like"]
    assert row["phylum"] == "Pisuviricota"
    assert row["class"] == "Pisoniviricetes"
    assert row["order"] == "Picornavirales"
    assert row["species"] == "Novel picorna-like virus"
    assert row["taxpath"].split("|")[:4] == [
        "10239",
        "2732408",
        "2732506",
        "464095",
    ]


def test_prepare_inputs_placeholder_names_do_not_depend_on_row_order(
    tmp_path: Path,
) -> None:
    """Contextual placeholder names should be stable across annotation ordering."""
    forward_dir = tmp_path / "forward"
    reverse_dir = tmp_path / "reverse"
    forward_dir.mkdir()
    reverse_dir.mkdir()
    fasta, annotations = write_wvdb_fixtures(forward_dir)

    with annotations.open(newline="", encoding="utf-8") as handle:
        annotation_rows = list(csv.reader(handle, delimiter="\t"))
    reversed_annotations = reverse_dir / annotations.name
    write_tsv_rows(
        reversed_annotations,
        [annotation_rows[0], *reversed(annotation_rows[1:])],
    )

    forward_lineages = forward_dir / "wvdb-v1.0.sourmash.lineages.csv"
    reverse_lineages = reverse_dir / "wvdb-v1.0.sourmash.lineages.csv"

    forward_result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--normalized-fasta",
        forward_dir / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        forward_lineages,
    )
    reverse_result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        reversed_annotations,
        "--normalized-fasta",
        reverse_dir / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        reverse_lineages,
    )

    assert forward_result.returncode == 0, forward_result.stderr
    assert reverse_result.returncode == 0, reverse_result.stderr

    assert read_csv_rows_by_ident(forward_lineages) == read_csv_rows_by_ident(
        reverse_lineages,
    )


def test_prepare_inputs_falls_through_sentinel_source_values(
    tmp_path: Path,
) -> None:
    """Sentinels fall through without combining incompatible source paths."""
    fasta = tmp_path / "WVDB_v1.0.fasta"
    annotations = tmp_path / "WVDB_v1.0_annotations.tsv"
    fasta.write_text(
        ">contig_with_fallbacks\nACGTACGT\n",
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
                "contig_with_fallbacks",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Fallbackfamily",
                "Unclassified",
                "Unknown",
                "Picornavirales",
                "Picornaviridae",
                "",
                "",
                "Unclassified",
                "Unclassified",
                "ICTVfamily",
            ],
        ],
    )
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--normalized-fasta",
        tmp_path / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    row = read_csv_rows_by_ident(lineages_csv)["WVDB|contig_with_fallbacks"]
    assert row["phylum"] == "Pisuviricota"
    assert row["class"] == "Pisoniviricetes"
    assert row["order"] == "Picornavirales"
    assert row["family"] == "Fallbackfamily"


def test_prepare_inputs_contextualizes_sentinel_variants(tmp_path: Path) -> None:
    """Common unknown-value spellings should become contextual placeholders."""
    fasta = tmp_path / "WVDB_v1.0.fasta"
    annotations = tmp_path / "WVDB_v1.0_annotations.tsv"
    fasta.write_text(
        ">contig_with_sentinels\nACGTACGT\n",
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
                "Genus_RdRp",
                "Species_RdRp",
            ],
            [
                "contig_with_sentinels",
                "  Unclassified  ",
                "Unknown",
                "NA",
                "n/a",
                "None",
                "Example species",
            ],
        ],
    )
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--normalized-fasta",
        tmp_path / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    row = read_csv_rows_by_ident(lineages_csv)["WVDB|contig_with_sentinels"]
    assert row["phylum"].startswith("unclassified phylum [under Viruses]")
    assert row["class"].startswith("unclassified class [under Viruses]")
    assert row["order"].startswith("unclassified order [under Viruses]")
    assert row["family"].startswith("unclassified family [under Viruses]")
    assert row["genus"].startswith("unclassified genus [under Viruses]")


def test_prepare_inputs_avoids_placeholder_label_collisions(
    tmp_path: Path,
) -> None:
    """Generated placeholder labels should not duplicate real source labels."""
    fasta = tmp_path / "WVDB_v1.0.fasta"
    annotations = tmp_path / "WVDB_v1.0_annotations.tsv"
    fasta.write_text(
        ">contig_real_label\nACGTACGT\n>contig_placeholder\nTGCATGCA\n",
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
            ],
            [
                "contig_real_label",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "unclassified family [under Picornavirales]",
            ],
            [
                "contig_placeholder",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Unclassified",
            ],
        ],
    )
    lineages_csv = tmp_path / "wvdb-v1.0.sourmash.lineages.csv"

    result = run_script(
        "prepare-inputs",
        "--fasta",
        fasta,
        "--annotations-tsv",
        annotations,
        "--normalized-fasta",
        tmp_path / "WVDB_v1.0.normalized.fasta",
        "--lineages-csv",
        lineages_csv,
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(lineages_csv)
    real_label = rows["WVDB|contig_real_label"]["family"]
    generated_label = rows["WVDB|contig_placeholder"]["family"]
    assert real_label == "unclassified family [under Picornavirales]"
    assert generated_label != real_label
    assert generated_label.startswith(
        "unclassified family [under Pisoniviricetes > Picornavirales]",
    )


def test_combine_lineages_preserves_taxpath_column(tmp_path: Path) -> None:
    """combine-lineages writes a publishable taxonomy without stripping taxpaths."""
    reference_lineages = write_reference_lineages(tmp_path / "reference.lineages.csv")
    wvdb_lineages = tmp_path / "wvdb.lineages.csv"
    combined_lineages = tmp_path / "combined.lineages.csv"
    wvdb_lineages.write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species,taxpath\n"
        "WVDB|contig_1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,,,,10239|2732408|2732506|464095\n",
        encoding="utf-8",
    )

    result = run_script(
        "combine-lineages",
        "--reference-lineages-csv",
        reference_lineages,
        "--wvdb-lineages-csv",
        wvdb_lineages,
        "--lineages-csv",
        combined_lineages,
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows(combined_lineages)
    assert rows[0]["ident"] == "REF_PICORNA"
    assert rows[-1]["ident"] == "WVDB|contig_1"
    assert rows[-1]["strain"] == ""
    assert rows[-1]["taxpath"] == "10239|2732408|2732506|464095"


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


def test_manifest_coverage_accepts_prepared_sourmash_taxonomy(
    tmp_path: Path,
) -> None:
    """check-manifest-coverage accepts sourmash tax prepare output."""
    manifest = tmp_path / "manifest.csv"
    lineages = tmp_path / "lineages.csv"

    manifest.write_text(
        "# SOURMASH-MANIFEST-VERSION: 1.0\n"
        "internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename\n"
        "signatures/b.sig.gz,def,def,31,DNA,0,50,10,False,WVDB|contig_1,ref.fa\n",
        encoding="utf-8",
    )
    lineages.write_text(
        "identifiers,superkingdom,phylum,class,order,family,genus,species,strain\n"
        "WVDB|contig_1,Viruses,,,,,,,\n",
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


def test_validate_lineages_rejects_named_strain_without_positional_taxid(
    tmp_path: Path,
) -> None:
    """Production validation reports the historical BioBoxes crash shape."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("GCA_002816655.1", "a" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\n"
        "GCA_002816655.1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,Coxsackievirus A16,"
        "10239|2732408|2732506|464095|12058|12059|138948\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "short taxpaths: 1" in result.stderr
    assert "populated ranks with missing taxids: 1" in result.stderr
    assert "GCA_002816655.1" in result.stderr


def test_validate_lineages_requires_all_declared_rank_columns(tmp_path: Path) -> None:
    """Production validation rejects a lineage schema missing the strain rank."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        "ident,superkingdom,phylum,class,order,family,genus,species,taxpath\n"
        "NCBI_1,Viruses,,,,,,,10239|||||||\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "missing required lineage columns: strain" in result.stderr


def test_validate_lineages_rejects_blank_identifier(tmp_path: Path) -> None:
    """Production validation reports lineage rows without an identifier."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\n,Viruses,,,,,,,,10239|||||||\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "blank identifiers: 1" in result.stderr


def test_validate_lineages_rejects_taxpaths_longer_than_declared_ranks(
    tmp_path: Path,
) -> None:
    """Production validation requires exactly one taxpath position per rank."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\nNCBI_1,Viruses,,,,,,,,10239||||||||extra\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "long taxpaths: 1" in result.stderr
    assert "NCBI_1" in result.stderr


def test_validate_lineages_rejects_taxid_without_rank_name(tmp_path: Path) -> None:
    """Production validation rejects positional taxids for absent rank names."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\n"
        "NCBI_1,Viruses,Pisuviricota,,,,,,,10239|2732408|999|||||\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "taxids without populated ranks: 1" in result.stderr
    assert "NCBI_1:class" in result.stderr


def test_validate_lineages_rejects_non_numeric_taxid(tmp_path: Path) -> None:
    """Production validation rejects taxids outside the decimal BioBoxes form."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\n"
        "NCBI_1,Viruses,Pisuviricota,,,,,,,10239|not-a-taxid||||||\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "non-numeric taxids: 1" in result.stderr
    assert "NCBI_1:phylum=not-a-taxid" in result.stderr


def test_validate_lineages_rejects_duplicate_identifiers(tmp_path: Path) -> None:
    """Production validation reports every repeated lineage identifier."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    row = "NCBI_1,Viruses,,,,,,,,10239|||||||\n"
    lineages.write_text(f"{LINEAGE_HEADER},taxpath\n{row}{row}", encoding="utf-8")

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "duplicate identifiers: 1" in result.stderr
    assert "NCBI_1" in result.stderr


def test_validate_lineages_rejects_duplicate_manifest_identifiers(
    tmp_path: Path,
) -> None:
    """Production validation rejects repeated names in the final manifest."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32), ("NCBI_1", "b" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\nNCBI_1,Viruses,,,,,,,,10239|||||||\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "duplicate manifest identifiers: 1" in result.stderr
    assert "NCBI_1" in result.stderr


def test_validate_lineages_requires_exact_manifest_identifier_coverage(
    tmp_path: Path,
) -> None:
    """Production lineages contain exactly the identifiers in the manifest."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_expected", "a" * 32), ("WVDB|missing", "b" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\n"
        "NCBI_expected,Viruses,,,,,,,,10239|||||||\n"
        "NCBI_stale,Viruses,,,,,,,,10239|||||||\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "manifest identifiers without lineages: 1" in result.stderr
    assert "WVDB|missing" in result.stderr
    assert "lineage identifiers absent from manifest: 1" in result.stderr
    assert "NCBI_stale" in result.stderr


def test_validate_lineages_rejects_conflicting_prefix_taxid_mappings(
    tmp_path: Path,
) -> None:
    """Production validation requires a bijection between nodes and taxids."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [
            ("NCBI_A", "a" * 32),
            ("NCBI_B", "b" * 32),
            ("NCBI_C", "c" * 32),
        ],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\n"
        "NCBI_A,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,,10239|1|2|3|4|5|6|\n"
        "NCBI_B,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,,10239|1|2|3|4|5|7|\n"
        "NCBI_C,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus B,,10239|1|2|3|4|5|6|\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "lineage prefixes with conflicting taxids: 1" in result.stderr
    assert "taxids with conflicting lineage prefixes: 1" in result.stderr


def test_validate_lineages_rejects_semantic_name_conflicts(tmp_path: Path) -> None:
    """The final publication gate repeats taxonomy-name safety checks."""
    manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("REF_1", "a" * 32), ("REF_2", "b" * 32)],
    )
    lineages = tmp_path / "combined.lineages.csv"
    lineages.write_text(
        f"{LINEAGE_HEADER},taxpath\n"
        "REF_1,Viruses,P1,C1,O1,F1,G1,Shared species,,1|2|3|4|5|6|7|\n"
        "REF_2,Viruses,P1,C1,O1,F2,G2,Shared species,,1|2|3|4|8|9|10|\n",
        encoding="utf-8",
    )

    result = run_script(
        "validate-lineages",
        "--manifest-csv",
        manifest,
        "--lineages-csv",
        lineages,
    )

    assert result.returncode == 1
    assert "same-rank taxonomy conflicts: 1" in result.stderr


def test_write_checksums_records_complete_artifact_pair(tmp_path: Path) -> None:
    """The ready marker records both runtime artifacts with SHA-256 digests."""
    signature = tmp_path / "reference.sig.zip"
    lineages = tmp_path / "reference.lineages.csv"
    signature.write_bytes(b"a")
    lineages.write_bytes(b"b")
    checksums = tmp_path / "reference.sha256"

    result = run_script(
        "write-checksums",
        "--output",
        checksums,
        signature,
        lineages,
    )

    assert result.returncode == 0, result.stderr
    assert checksums.read_text(encoding="utf-8") == (
        "ca978112ca1bbdcafac231b39a23dc4da786eff8147c4e72b9807785afee48bb  "
        "reference.sig.zip\n"
        "3e23e8160039594a33894f6564e1b1348bbd7a0088d42c4acb73eeaed59c009d  "
        "reference.lineages.csv\n"
    )


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


def test_curate_lineages_most_specified_rewrites_conflicting_parents(
    tmp_path: Path,
) -> None:
    """The uniquely most-specified lineage should become canonical with a warning."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1 reference title", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            (
                "NCBI_1",
                "Viruses",
                "unclassified Viruses phylum",
                "unclassified Viruses class",
                "unclassified Viruses order",
                "unclassified Viruses family",
                "unclassified Viruses genus",
                "Shared virus",
                "ncbi-strain",
            ),
        ],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Picornaviridae",
                "",
                "Shared_virus",
                "wvdb-strain",
            ),
        ],
    )

    result, output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    output_rows = read_csv_rows(output)
    rows = {row["ident"]: row for row in output_rows}
    assert [row["ident"] for row in output_rows] == ["NCBI_1", "WVDB|1"]
    assert rows["NCBI_1"]["phylum"] == "Pisuviricota"
    assert rows["NCBI_1"]["species"] == "Shared_virus"
    assert rows["NCBI_1"]["strain"] == "ncbi-strain"
    assert rows["WVDB|1"]["strain"] == "wvdb-strain"
    assert "WARNING" in result.stderr
    assert "resolved 1 taxonomy conflict" in result.stderr
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["normalized_name"] for row in diagnostic_rows} == {"shared virus"}
    assert {row["decision"] for row in diagnostic_rows} == {"winner", "rewritten"}


def test_curate_lineages_rewrites_taxpath_with_taxonomy_prefix(tmp_path: Path) -> None:
    """A rewritten taxonomy prefix must retain aligned BioBoxes taxids."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    header = f"{LINEAGE_HEADER},taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,unclassified phylum,unclassified class,"
        "unclassified order,unclassified family,unclassified genus,"
        "Shared virus,ncbi-strain,10239|1|2|3|4|5|6|7\n",
        encoding="utf-8",
    )
    wvdb_lineages = tmp_path / "wvdb.lineages.csv"
    wvdb_lineages.write_text(
        header + "WVDB|1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,,Shared_virus,wvdb-strain,"
        "10239|11|12|13|14||16|17\n",
        encoding="utf-8",
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    assert rows["NCBI_1"]["taxpath"] == "10239|11|12|13|14||6|7"
    assert rows["WVDB|1"]["taxpath"] == "10239|11|12|13|14||6|17"


def test_curate_lineages_assigns_taxid_to_named_strain_after_prefix_rewrite(
    tmp_path: Path,
) -> None:
    """A named strain remains addressable after its taxonomy prefix is rewritten."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    header = f"{LINEAGE_HEADER},taxid,taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,unclassified phylum,unclassified class,"
        "unclassified order,unclassified family,unclassified genus,"
        "Shared virus,Named strain,33757,10239|1|2|3|4|5|6|\n",
        encoding="utf-8",
    )
    wvdb_lineages = tmp_path / "wvdb.lineages.csv"
    wvdb_lineages.write_text(
        header + "WVDB|1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,,Shared_virus,,,10239|11|12|13|14||16|\n",
        encoding="utf-8",
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    row = read_csv_rows_by_ident(output)["NCBI_1"]
    taxids = row["taxpath"].split("|")
    assert len(taxids) == len(LINEAGE_HEADER.split(",")) - 1
    assert row["strain"] == "Named strain"
    assert taxids[-1] == "33757"


def test_curate_lineages_uses_one_taxid_for_equivalent_lineage_prefixes(
    tmp_path: Path,
) -> None:
    """Equivalent finalized lineage prefixes receive one shared taxid path."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    header = f"{LINEAGE_HEADER},taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,,10239|1|2|3|4|5|6|\n",
        encoding="utf-8",
    )
    wvdb_lineages = tmp_path / "wvdb.lineages.csv"
    wvdb_lineages.write_text(
        header + "WVDB|1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,,10239|11|12|13|14|15|16|\n",
        encoding="utf-8",
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    assert rows["NCBI_1"]["taxpath"] == rows["WVDB|1"]["taxpath"]


def test_curate_lineages_synthesizes_conflicting_ncbi_taxids_for_one_prefix(
    tmp_path: Path,
) -> None:
    """Contradictory native IDs produce one explicit synthetic canonical node."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32), ("NCBI_2", "b" * 32)],
    )
    wvdb_manifest = write_manifest(tmp_path / "wvdb.manifest.csv", [])
    header = f"{LINEAGE_HEADER},taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,,10239|1|2|3|4|5|6|\n"
        "NCBI_2,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,,10239|1|2|3|4|5|7|\n",
        encoding="utf-8",
    )
    wvdb_lineages = write_lineages(tmp_path / "wvdb.lineages.csv", [])

    result, output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    assert rows["NCBI_1"]["taxpath"] == rows["NCBI_2"]["taxpath"]
    assert rows["NCBI_1"]["taxpath"].split("|")[6] not in {"6", "7"}
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert diagnostic_rows[-1]["resolution"] == (
        "synthesized-conflicting-native-taxids"
    )


def test_curate_lineages_rejects_one_taxid_for_conflicting_lineage_prefixes(
    tmp_path: Path,
) -> None:
    """One taxid cannot identify two distinct finalized taxonomy nodes."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    header = f"{LINEAGE_HEADER},taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus A,,10239|1|2|3|4|5|999|\n",
        encoding="utf-8",
    )
    wvdb_lineages = tmp_path / "wvdb.lineages.csv"
    wvdb_lineages.write_text(
        header + "WVDB|1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Enterovirus B,,10239|1|2|3|4|5|999|\n",
        encoding="utf-8",
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 1
    assert "taxid 999 represents conflicting lineage prefixes" in result.stderr
    assert not output.exists()


def test_curate_lineages_hashes_delimiter_bearing_prefixes_unambiguously(
    tmp_path: Path,
) -> None:
    """Taxonomy label delimiters cannot alias distinct synthetic-ID prefixes."""
    reference_manifest = write_manifest(tmp_path / "reference.manifest.csv", [])
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "a" * 32), ("WVDB|2", "b" * 32)],
    )
    reference_lineages = write_lineages(tmp_path / "reference.lineages.csv", [])
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            ("WVDB|1", "Viruses", "A|class=B", "", "", "", "", "", ""),
            ("WVDB|2", "Viruses", "A", "B", "", "", "", "", ""),
        ],
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    first_taxids = rows["WVDB|1"]["taxpath"].split("|")
    second_taxids = rows["WVDB|2"]["taxpath"].split("|")
    assert first_taxids[1] != second_taxids[2]


@settings(
    deadline=None,
    max_examples=30,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    optional_labels=st.lists(
        st.one_of(
            st.none(),
            st.text(
                alphabet=st.characters(
                    blacklist_categories=("Cs",),
                    blacklist_characters="\x00\r\n",
                ),
                min_size=1,
                max_size=20,
            ),
        ),
        min_size=7,
        max_size=7,
    ),
)
def test_curate_lineages_materializes_bioboxes_positions_for_sparse_taxonomy(
    tmp_path: Path,
    optional_labels: list[str | None],
) -> None:
    """Every generated rank has one numeric taxid in the declared position."""
    ranks = LINEAGE_HEADER.split(",")[1:]
    names = [
        f"{rank} {label}" if label is not None else ""
        for rank, label in zip(ranks[1:], optional_labels, strict=True)
    ]
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(tmp_path / "wvdb.manifest.csv", [])
    reference_lineages = tmp_path / "reference.lineages.csv"
    with reference_lineages.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([*LINEAGE_HEADER.split(","), "taxpath"])
        writer.writerow(["NCBI_1", "Viruses", *names, "10239|||||||"])
    wvdb_lineages = write_lineages(tmp_path / "wvdb.lineages.csv", [])

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    row = read_csv_rows_by_ident(output)["NCBI_1"]
    taxids = row["taxpath"].split("|")
    assert len(taxids) == len(ranks)
    for rank, taxid in zip(ranks, taxids, strict=True):
        assert bool(row[rank]) == bool(taxid)
        assert not taxid or taxid.isdecimal()


@settings(
    deadline=None,
    max_examples=20,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    labels=st.lists(
        st.text(
            alphabet=st.characters(
                blacklist_categories=("Cs",),
                blacklist_characters="\x00\r\n",
            ),
            min_size=1,
            max_size=20,
        ),
        min_size=2,
        max_size=2,
    ),
)
def test_curated_taxids_are_stable_across_input_order(
    tmp_path: Path,
    labels: list[str],
) -> None:
    """Row order cannot change shared or synthetic lineage-node taxids."""
    rows = [
        (
            "WVDB|1",
            "Viruses",
            "Pisuviricota",
            "Pisoniviricetes",
            "Picornavirales",
            "Picornaviridae",
            "Enterovirus",
            f"Species 1 {labels[0]}",
            "",
        ),
        (
            "WVDB|2",
            "Viruses",
            "Pisuviricota",
            "Pisoniviricetes",
            "Picornavirales",
            "Picornaviridae",
            "Enterovirus",
            f"Species 2 {labels[1]}",
            "",
        ),
    ]
    outputs: list[dict[str, dict[str, str]]] = []
    for name, ordered_rows in (("forward", rows), ("reverse", list(reversed(rows)))):
        case_dir = tmp_path / name
        case_dir.mkdir(exist_ok=True)
        reference_manifest = write_manifest(case_dir / "reference.manifest.csv", [])
        wvdb_manifest = write_manifest(
            case_dir / "wvdb.manifest.csv",
            [
                (row[0], digest * 32)
                for row, digest in zip(ordered_rows, ("a", "b"), strict=True)
            ],
        )
        reference_lineages = write_lineages(case_dir / "reference.lineages.csv", [])
        wvdb_lineages = write_lineages(
            case_dir / "wvdb.lineages.csv",
            ordered_rows,
        )
        result, output, _diagnostics = run_curation(
            case_dir,
            manifests=(reference_manifest, wvdb_manifest),
            lineages=(reference_lineages, wvdb_lineages),
        )
        assert result.returncode == 0, result.stderr
        outputs.append(read_csv_rows_by_ident(output))

    assert outputs[0] == outputs[1]
    first_taxids = outputs[0]["WVDB|1"]["taxpath"].split("|")
    second_taxids = outputs[0]["WVDB|2"]["taxpath"].split("|")
    assert first_taxids[:6] == second_taxids[:6]


def test_curate_lineages_rejects_non_numeric_taxids(tmp_path: Path) -> None:
    """Published BioBoxes taxids must be decimal strings."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(tmp_path / "wvdb.manifest.csv", [])
    header = f"{LINEAGE_HEADER},taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,Pisuviricota,,,,,,,10239|not-a-taxid||||||\n",
        encoding="utf-8",
    )
    wvdb_lineages = write_lineages(tmp_path / "wvdb.lineages.csv", [])

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 1
    assert "non-numeric taxid: not-a-taxid" in result.stderr
    assert not output.exists()


def test_curate_lineages_rejects_source_taxpath_beyond_declared_ranks(
    tmp_path: Path,
) -> None:
    """Curation cannot silently discard source taxpath positions."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(tmp_path / "wvdb.manifest.csv", [])
    header = f"{LINEAGE_HEADER},taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,,,,,,,,10239||||||||999\n",
        encoding="utf-8",
    )
    wvdb_lineages = write_lineages(tmp_path / "wvdb.lineages.csv", [])

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 1
    assert "source taxpath has 9 positions; expected at most 8" in result.stderr
    assert not output.exists()


def test_curate_lineages_uses_informative_ncbi_corroboration(tmp_path: Path) -> None:
    """NCBI resolves a tied placement only when it supports an informative path."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32), ("WVDB|2", "c" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            (
                "NCBI_1",
                "Viruses",
                "Lenarviricota",
                "Leviviricetes",
                "Norzivirales",
                "Fiersviridae",
                "",
                "",
                "",
            ),
        ],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Lenarviricota",
                "Leviviricetes",
                "Norzivirales",
                "Fiersviridae",
                "",
                "",
                "",
            ),
            (
                "WVDB|2",
                "Viruses",
                "Lenarviricota",
                "Leviviricetes",
                "Timlovirales",
                "Fiersviridae",
                "",
                "",
                "",
            ),
        ],
    )

    result, output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    assert {row["order"] for row in rows.values()} == {"Norzivirales"}
    assert {row["family"] for row in rows.values()} == {"Fiersviridae"}
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["resolution"] for row in diagnostic_rows} == {
        "selected-informative-ncbi",
    }
    assert {row["decision"] for row in diagnostic_rows} == {"winner", "rewritten"}


def test_curate_lineages_contextualizes_equally_specified_conflicts(
    tmp_path: Path,
) -> None:
    """Incompatible tied placements retain distinct contextual taxon names."""
    reference_manifest = write_manifest(tmp_path / "reference.manifest.csv", [])
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "a" * 32), ("WVDB|2", "b" * 32)],
    )
    reference_lineages = write_lineages(tmp_path / "reference.lineages.csv", [])
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Picornaviridae",
                "",
                "Shared_virus",
                "",
            ),
            (
                "WVDB|2",
                "Viruses",
                "Pisuviricota",
                "Pisoniviricetes",
                "Picornavirales",
                "Dicistroviridae",
                "",
                "Shared_virus",
                "",
            ),
        ],
    )

    result, output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    assert rows["WVDB|1"]["species"] == "Shared_virus [under Picornaviridae]"
    assert rows["WVDB|2"]["species"] == "Shared_virus [under Dicistroviridae]"
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["conflict_type"] for row in diagnostic_rows} == {"taxonomy_name"}
    assert {row["resolution"] for row in diagnostic_rows} == {
        "contextualized-ambiguity",
    }
    assert {row["decision"] for row in diagnostic_rows} == {"contextualized"}


def test_curate_lineages_treats_placeholder_spellings_as_equivalent(
    tmp_path: Path,
) -> None:
    """Placeholder presentation differences do not create distinct parents."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            (
                "NCBI_1",
                "Viruses",
                "Kitrinoviricota",
                "Flasuviricetes",
                "Amarillovirales",
                "Flaviviridae",
                "unclassified Flaviviridae genus",
                "Apis flavivirus",
                "",
            ),
        ],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Kitrinoviricota",
                "Flasuviricetes",
                "Amarillovirales",
                "Flaviviridae",
                "unclassified genus [under Flaviviridae]",
                "Apis_flavivirus",
                "",
            ),
        ],
    )

    result, output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    assert rows["NCBI_1"]["species"] == "Apis flavivirus"
    assert rows["WVDB|1"]["species"] == "Apis flavivirus"
    assert rows["WVDB|1"]["genus"] == "unclassified Flaviviridae genus"
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["resolution"] for row in diagnostic_rows} == {
        "canonicalized-placeholder",
    }
    assert {row["decision"] for row in diagnostic_rows} == {"winner", "rewritten"}


def test_curate_lineages_reuses_ncbi_taxid_for_normalized_equivalent_name(
    tmp_path: Path,
) -> None:
    """Underscore spelling differences retain the canonical NCBI node ID."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    header = f"{LINEAGE_HEADER},taxpath\n"
    reference_lineages = tmp_path / "reference.lineages.csv"
    reference_lineages.write_text(
        header + "NCBI_1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Example virus,,10239|1|2|3|4|5|6|\n",
        encoding="utf-8",
    )
    wvdb_lineages = tmp_path / "wvdb.lineages.csv"
    wvdb_lineages.write_text(
        header + "WVDB|1,Viruses,Pisuviricota,Pisoniviricetes,Picornavirales,"
        "Picornaviridae,Enterovirus,Example_virus,,10239|11|12|13|14|15|16|\n",
        encoding="utf-8",
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    rows = read_csv_rows_by_ident(output)
    assert rows["WVDB|1"]["species"] == "Example virus"
    assert rows["WVDB|1"]["taxpath"] == rows["NCBI_1"]["taxpath"]
    assert rows["WVDB|1"]["taxpath"].split("|")[6] == "6"
    combined_manifest = write_manifest(
        tmp_path / "combined.manifest.csv",
        [("NCBI_1", "a" * 32), ("WVDB|1", "b" * 32)],
    )
    validation = run_script(
        "validate-lineages",
        "--manifest-csv",
        combined_manifest,
        "--lineages-csv",
        output,
    )
    assert validation.returncode == 0, validation.stderr


def test_curate_lineages_fails_on_duplicate_signature_identity(
    tmp_path: Path,
) -> None:
    """Duplicate sketch identities remain a hard build error with evidence."""
    digest = "c" * 32
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", digest)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", digest)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [("NCBI_1", "Viruses", "Reference phylum", "", "", "", "", "Virus", "")],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [("WVDB|1", "Viruses", "WVDB phylum", "", "", "", "", "Virus", "")],
    )

    result, _output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 1
    assert digest in result.stderr
    assert "NCBI_1" in result.stderr
    assert "WVDB|1" in result.stderr
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert [row["ident"] for row in diagnostic_rows] == ["NCBI_1", "WVDB|1"]
    assert {row["conflict_type"] for row in diagnostic_rows} == {"signature_identity"}


def test_curate_lineages_rejects_duplicate_manifest_identifiers(tmp_path: Path) -> None:
    """Each retained signature identifier must produce exactly one lineage row."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32), ("NCBI_1", "b" * 32)],
    )
    wvdb_manifest = write_manifest(tmp_path / "wvdb.manifest.csv", [])
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [("NCBI_1", "Viruses", "Pisuviricota", "", "", "", "", "", "")],
    )
    wvdb_lineages = write_lineages(tmp_path / "wvdb.lineages.csv", [])

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 1
    assert "Duplicate manifest identifier: NCBI_1" in result.stderr
    assert not output.exists()


def test_curate_lineages_removes_rows_absent_from_manifests(tmp_path: Path) -> None:
    """Published taxonomy should contain only identifiers retained in the collection."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(tmp_path / "wvdb.manifest.csv", [])
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            ("NCBI_1", "Viruses", "Pisuviricota", "", "", "", "", "Virus 1", ""),
            ("NCBI_extra", "Viruses", "Pisuviricota", "", "", "", "", "Virus 2", ""),
        ],
    )
    wvdb_lineages = write_lineages(tmp_path / "wvdb.lineages.csv", [])

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    assert [row["ident"] for row in read_csv_rows(output)] == ["NCBI_1"]


def test_signature_identity_includes_sketch_parameters(tmp_path: Path) -> None:
    """Matching digests with different sketch parameters are distinct identities."""
    digest = "d" * 32
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", digest)],
        scaled="50",
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", digest)],
        scaled="100",
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [("NCBI_1", "Viruses", "Phylum A", "", "", "", "", "Virus A", "")],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [("WVDB|1", "Viruses", "Phylum B", "", "", "", "", "Virus B", "")],
    )

    result, output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 0, result.stderr
    assert [row["ident"] for row in read_csv_rows(output)] == ["NCBI_1", "WVDB|1"]


def test_curate_lineages_rejects_duplicate_identifiers(tmp_path: Path) -> None:
    """Lineage identifiers must be unique across both taxonomy sources."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("SHARED", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [("SHARED", "Viruses", "Phylum A", "", "", "", "", "Virus A", "")],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            ("SHARED", "Viruses", "Phylum B", "", "", "", "", "Virus B", ""),
            ("WVDB|1", "Viruses", "Phylum B", "", "", "", "", "Virus C", ""),
        ],
    )

    result, _output, _diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 1
    assert "Duplicate lineage identifier: SHARED" in result.stderr


def test_curate_lineages_rejects_cross_rank_name_reuse(tmp_path: Path) -> None:
    """Taxburst requires node names to be globally unique across taxonomy ranks."""
    reference_manifest = write_manifest(
        tmp_path / "reference.manifest.csv",
        [("NCBI_1", "a" * 32)],
    )
    wvdb_manifest = write_manifest(
        tmp_path / "wvdb.manifest.csv",
        [("WVDB|1", "b" * 32)],
    )
    reference_lineages = write_lineages(
        tmp_path / "reference.lineages.csv",
        [
            (
                "NCBI_1",
                "Viruses",
                "Phylum A",
                "Class A",
                "Order A",
                "Family A",
                "Repeated taxon",
                "Species A",
                "",
            ),
        ],
    )
    wvdb_lineages = write_lineages(
        tmp_path / "wvdb.lineages.csv",
        [
            (
                "WVDB|1",
                "Viruses",
                "Phylum B",
                "Class B",
                "Order B",
                "Family B",
                "Genus B",
                "Repeated_taxon",
                "",
            ),
        ],
    )

    result, _output, diagnostics = run_curation(
        tmp_path,
        manifests=(reference_manifest, wvdb_manifest),
        lineages=(reference_lineages, wvdb_lineages),
    )

    assert result.returncode == 1
    assert "cross-rank taxonomy name conflict" in result.stderr
    diagnostic_rows = read_tsv_rows(diagnostics)
    assert {row["conflict_type"] for row in diagnostic_rows} == {"taxonomy_rank"}
