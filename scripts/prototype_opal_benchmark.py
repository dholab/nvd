"""Prototype OPAL-style comparison inputs from e2e-like BLAST final results.

This script intentionally uses a tiny synthetic fixture. It exercises the core
scientific rule needed for NVD BLAST final results: repeated BLAST hit rows are
collapsed to one contig assignment before mapped read support is counted.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import shlex
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

BIOBOXES_VERSION = "0.10.0"
BIOBOXES_RANKS = "superkingdom|phylum|class|order|family|genus|species"
BIOBOXES_RANK_LIST = BIOBOXES_RANKS.split("|")
DEFAULT_OUTPUT_DIR = Path("prototype_outputs/opal_cami")
SAMPLE_ID = "cami_complex_synthetic_fixture"


@dataclass(frozen=True)
class TaxonFixture:
    """Explicit synthetic BioBoxes taxpath metadata for a leaf assignment."""

    taxid: str
    rank: str
    taxpath: str
    taxpathsn: str


@dataclass(frozen=True)
class ContigAssignment:
    """One collapsed contig assignment with read-segment support."""

    qseqid: str
    taxid: str
    mapped_reads: int
    selected_bitscore: float
    hit_rows_collapsed: int


TAXONOMY_FIXTURE = {
    "900111": TaxonFixture(
        taxid="900111",
        rank="species",
        taxpath="10239|900001|900010|900011|900012|900013|900111",
        taxpathsn="Viruses|Syntheticviricota|Syntheticviricetes|Syntheticvirales|Syntheticviridae|Syntheticvirus|Synthetic virus alpha",
    ),
    "900222": TaxonFixture(
        taxid="900222",
        rank="species",
        taxpath="10239|900002|900020|900022|900023|900024|900222",
        taxpathsn="Viruses|Fixtureviricota|Fixtureviricetes|Fixturevirales|Fixtureviridae|Fixturevirus|Synthetic virus beta",
    ),
    "900333": TaxonFixture(
        taxid="900333",
        rank="species",
        taxpath="10239|900003|900030|900033|900034|900035|900333",
        taxpathsn="Viruses|Benchmarkviricota|Benchmarkviricetes|Benchmarkvirales|Benchmarkviridae|Benchmarkvirus|Synthetic virus gamma",
    ),
}

SYNTHETIC_BLAST_FINAL_ROWS = [
    {
        "sample": "cami_complex_synthetic_fixture",
        "qseqid": "contig_alpha_001",
        "sseqid": "alpha_ref_a",
        "bitscore": "440.0",
        "adjusted_taxid": "900111",
        "adjusted_taxid_name": "Synthetic virus alpha",
        "mapped_reads": "120",
    },
    {
        "sample": "cami_complex_synthetic_fixture",
        "qseqid": "contig_alpha_001",
        "sseqid": "alpha_ref_b",
        "bitscore": "438.0",
        "adjusted_taxid": "900111",
        "adjusted_taxid_name": "Synthetic virus alpha",
        "mapped_reads": "120",
    },
    {
        "sample": "cami_complex_synthetic_fixture",
        "qseqid": "contig_alpha_002",
        "sseqid": "alpha_ref_c",
        "bitscore": "390.0",
        "adjusted_taxid": "900111",
        "adjusted_taxid_name": "Synthetic virus alpha",
        "mapped_reads": "75",
    },
    {
        "sample": "cami_complex_synthetic_fixture",
        "qseqid": "contig_beta_001",
        "sseqid": "beta_ref_a",
        "bitscore": "410.0",
        "adjusted_taxid": "900222",
        "adjusted_taxid_name": "Synthetic virus beta",
        "mapped_reads": "95",
    },
    {
        "sample": "cami_complex_synthetic_fixture",
        "qseqid": "contig_beta_001",
        "sseqid": "beta_ref_b",
        "bitscore": "407.0",
        "adjusted_taxid": "900222",
        "adjusted_taxid_name": "Synthetic virus beta",
        "mapped_reads": "95",
    },
    {
        "sample": "cami_complex_synthetic_fixture",
        "qseqid": "contig_gamma_001",
        "sseqid": "gamma_ref_a",
        "bitscore": "365.0",
        "adjusted_taxid": "900333",
        "adjusted_taxid_name": "Synthetic virus gamma",
        "mapped_reads": "30",
    },
]

SYNTHETIC_SOURMASH_SUPPORT = {
    "900111": 188,
    "900222": 105,
    "900333": 27,
}

SYNTHETIC_GOLD_SUPPORT = {
    "900111": 200,
    "900222": 100,
    "900333": 25,
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Write prototype OPAL BioBoxes profiles for NVD benchmarking.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for prototype artifacts.",
    )
    parser.add_argument(
        "--opal-command",
        help=(
            "Command prefix used to run OPAL, for example: "
            "'pixi exec -c conda-forge -c bioconda -s cami-opal opal.py'. "
            "When omitted, the script looks for opal.py or opal on PATH."
        ),
    )
    parser.add_argument(
        "--gold-profile",
        type=Path,
        help="External CAMI BioBoxes gold profile to copy into the prototype output.",
    )
    parser.add_argument(
        "--blast-profile",
        type=Path,
        help="External BLAST-derived CAMI BioBoxes profile to compare with OPAL.",
    )
    parser.add_argument(
        "--sourmash-profile",
        type=Path,
        help="External sourmash CAMI BioBoxes profile to compare with OPAL.",
    )
    return parser.parse_args()


def external_profiles_requested(args: argparse.Namespace) -> bool:
    """Return true when any external profile path is set."""
    return any([args.gold_profile, args.blast_profile, args.sourmash_profile])


def copy_external_profiles(args: argparse.Namespace, output_dir: Path) -> list[str]:
    """Copy a complete external profile triplet into the prototype output layout."""
    required_profiles = {
        "--gold-profile": args.gold_profile,
        "--blast-profile": args.blast_profile,
        "--sourmash-profile": args.sourmash_profile,
    }
    missing = [flag for flag, path in required_profiles.items() if path is None]
    if missing:
        msg = "External profile mode requires all profile paths: " + ", ".join(missing)
        raise ValueError(msg)

    output_dir.mkdir(parents=True, exist_ok=True)
    profile_destinations = {
        args.gold_profile: output_dir / "gold.profile",
        args.blast_profile: output_dir / "blast_assembled_read_support.profile",
        args.sourmash_profile: output_dir / "sourmash_read_support.profile",
    }
    for source, destination in profile_destinations.items():
        if not source.is_file():
            msg = f"External profile does not exist or is not a file: {source}"
            raise FileNotFoundError(msg)
        shutil.copy2(source, destination)

    return [
        "External profile triplet copied into prototype output layout.",
        "The script did not validate that external profiles share SampleID, ranks, or taxonomy.",
    ]


def collapse_blast_contigs(
    rows: list[dict[str, str]],
) -> tuple[list[ContigAssignment], list[str]]:
    """Collapse repeated BLAST hit rows to one mapped-read contribution per contig."""
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["qseqid"]].append(row)

    assignments = []
    warnings = []
    for qseqid, hit_rows in sorted(grouped.items()):
        taxids = [row["adjusted_taxid"] for row in hit_rows]
        taxid_counts = Counter(taxids)
        selected = max(hit_rows, key=lambda row: float(row["bitscore"]))
        if len(taxid_counts) > 1:
            warnings.append(
                f"{qseqid} has conflicting adjusted_taxid values; selected highest bitscore row",
            )
        assignments.append(
            ContigAssignment(
                qseqid=qseqid,
                taxid=selected["adjusted_taxid"],
                mapped_reads=int(selected["mapped_reads"]),
                selected_bitscore=float(selected["bitscore"]),
                hit_rows_collapsed=len(hit_rows),
            ),
        )
    return assignments, warnings


def support_by_taxid(assignments: list[ContigAssignment]) -> dict[str, int]:
    """Sum assembled read-segment support by collapsed contig taxid."""
    support: dict[str, int] = defaultdict(int)
    for assignment in assignments:
        support[assignment.taxid] += assignment.mapped_reads
    return dict(sorted(support.items()))


def write_bioboxes_profile(
    path: Path,
    *,
    sample_id: str,
    support: dict[str, int],
) -> None:
    """Write a rank-complete CAMI BioBoxes-compatible taxonomic profile."""
    total_support = sum(support.values())
    if total_support <= 0:
        msg = "Cannot write BioBoxes profile with zero total support"
        raise ValueError(msg)

    support_by_ranked_path: dict[tuple[str, str, str, str], int] = defaultdict(int)
    for leaf_taxid, reads in support.items():
        taxon = TAXONOMY_FIXTURE[leaf_taxid]
        taxpath_parts = taxon.taxpath.split("|")
        taxpathsn_parts = taxon.taxpathsn.split("|")
        if len(taxpath_parts) != len(BIOBOXES_RANK_LIST) or len(taxpathsn_parts) != len(
            BIOBOXES_RANK_LIST,
        ):
            msg = f"Taxpath for {leaf_taxid} is not aligned to declared BioBoxes ranks"
            raise ValueError(msg)
        for rank_index, rank in enumerate(BIOBOXES_RANK_LIST):
            ranked_taxid = taxpath_parts[rank_index]
            ranked_taxpath = "|".join(taxpath_parts[: rank_index + 1])
            ranked_taxpathsn = "|".join(taxpathsn_parts[: rank_index + 1])
            support_by_ranked_path[
                (rank, ranked_taxid, ranked_taxpath, ranked_taxpathsn)
            ] += reads

    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(f"@SampleID:{sample_id}\n")
        handle.write(f"@Version:{BIOBOXES_VERSION}\n")
        handle.write(f"@Ranks:{BIOBOXES_RANKS}\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["@@TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"])
        for rank in BIOBOXES_RANK_LIST:
            rank_rows = [
                (key, reads)
                for key, reads in support_by_ranked_path.items()
                if key[0] == rank
            ]
            for (_rank, taxid, taxpath, taxpathsn), reads in sorted(rank_rows):
                percentage = reads / total_support * 100
                writer.writerow(
                    [
                        taxid,
                        rank,
                        taxpath,
                        taxpathsn,
                        f"{percentage:.6f}",
                    ],
                )


def opal_command_prefix(opal_command: str | None) -> list[str] | None:
    """Resolve an explicit or PATH-discovered OPAL command prefix."""
    if opal_command:
        return shlex.split(opal_command)
    opal = shutil.which("opal.py") or shutil.which("opal")
    if opal is None:
        return None
    return [opal]


def run_opal_if_available(
    output_dir: Path,
    *,
    opal_command: str | None,
) -> tuple[bool, str, str]:
    """Run OPAL when installed; otherwise write an honest prototype HTML report."""
    opal_out = output_dir / "opal_out"
    if opal_out.exists():
        shutil.rmtree(opal_out)
    command_prefix = opal_command_prefix(opal_command)
    command = [
        *(command_prefix or ["opal.py"]),
        "-g",
        str(output_dir / "gold.profile"),
        "-o",
        str(opal_out),
        "-l",
        "blast_assembled_read_support,sourmash_read_support",
        str(output_dir / "blast_assembled_read_support.profile"),
        str(output_dir / "sourmash_read_support.profile"),
    ]
    command_text = shlex.join(command)
    if command_prefix is None:
        write_prototype_html(opal_out / "index.html", opal_command=command_text)
        return (
            False,
            "OPAL executable was not available; wrote prototype HTML report",
            command_text,
        )

    opal_out.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(command, text=True, capture_output=True, check=False)  # noqa: S603
    results_html = opal_out / "results.html"
    index_html = opal_out / "index.html"
    html_report = results_html.exists() or index_html.exists()
    if result.returncode != 0 or not html_report:
        write_prototype_html(
            index_html,
            opal_command=command_text,
            stderr=result.stderr,
        )
        return False, "OPAL command failed; wrote prototype HTML report", command_text
    if results_html.exists() and not index_html.exists():
        shutil.copy2(results_html, index_html)
    return True, "OPAL completed and produced HTML/results artifacts", command_text


def write_prototype_html(
    path: Path,
    *,
    opal_command: str | None,
    stderr: str = "",
) -> None:
    """Write a clearly labeled fallback comparison report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    escaped_command = html.escape(opal_command or "not attempted")
    escaped_stderr = html.escape(stderr)
    path.write_text(
        "\n".join(
            [
                "<!doctype html>",
                '<html lang="en">',
                "<head>",
                '<meta charset="utf-8">',
                "<title>NVD OPAL CAMI Prototype Report</title>",
                "<style>body{font-family:system-ui,sans-serif;max-width:900px;margin:3rem auto;line-height:1.5}code,pre{background:#f4f4f4;padding:.2rem .35rem}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:.4rem .6rem}</style>",
                "</head>",
                "<body>",
                "<h1>NVD OPAL CAMI Prototype Report</h1>",
                "<p><strong>Status:</strong> prototype fallback, not an OPAL-generated benchmark report.</p>",
                "<p>The profiles are CAMI BioBoxes-style fixtures for checking assembled-read-support semantics before wiring the requested future complex CAMI dataset.</p>",
                "<h2>Inputs</h2>",
                "<table><tr><th>Profile</th><th>Description</th></tr>",
                "<tr><td>gold.profile</td><td>Explicit synthetic gold fixture over the same sample ID and taxpaths.</td></tr>",
                "<tr><td>blast_assembled_read_support.profile</td><td>Collapsed e2e-like BLAST final rows; each contig contributes mapped_reads once to adjusted_taxid.</td></tr>",
                "<tr><td>sourmash_read_support.profile</td><td>Comparable synthetic sourmash-style read support over the same explicit fixture taxpaths.</td></tr>",
                "</table>",
                "<h2>OPAL Execution</h2>",
                "<p>If OPAL was unavailable or failed locally, this is the exact command that would have been run.</p>",
                f"<p><code>{escaped_command}</code></p>",
                f"<pre>{escaped_stderr}</pre>",
                "</body>",
                "</html>",
            ],
        )
        + "\n",
        encoding="utf-8",
    )


def write_qc(
    path: Path,
    *,
    assignments: list[ContigAssignment],
    blast_support: dict[str, int],
    opal_result: tuple[bool, str, str],
    warnings: list[str],
) -> None:
    """Write audit metadata for parent-agent review."""
    opal_ran, opal_status, opal_command = opal_result
    total_hit_rows = len(SYNTHETIC_BLAST_FINAL_ROWS)
    qc = {
        "created_at_utc": datetime.now(tz=UTC).isoformat(),
        "dataset": {
            "requested_future_target": "one complex CAMI dataset for real OPAL benchmarking",
            "actual_source": "synthetic e2e-like fixture embedded in scripts/prototype_opal_benchmark.py",
            "prototype_status": "No large CAMI dataset was downloaded or consumed in this iteration.",
            "taxonomy_status": "Synthetic explicit taxids and taxpaths are used; they are not NCBI-validated taxids.",
            "taxpath_status": "Synthetic TAXPATH and TAXPATHSN values are rank-aligned fixtures for the declared CAMI ranks.",
            "sample_id": SAMPLE_ID,
        },
        "blast_semantics": {
            "input_hit_rows": total_hit_rows,
            "collapsed_contigs": len(assignments),
            "repeated_hit_rows_collapsed": total_hit_rows - len(assignments),
            "mass_definition": "V0 assembled read-segment support: each contig contributes mapped_reads once to its consensus adjusted_taxid.",
            "support_by_taxid": blast_support,
            "contig_assignments": [assignment.__dict__ for assignment in assignments],
        },
        "sourmash_semantics": {
            "actual_source": "synthetic comparable fixture profile",
            "support_by_taxid": SYNTHETIC_SOURMASH_SUPPORT,
        },
        "gold_semantics": {
            "actual_source": "explicit synthetic gold fixture profile",
            "support_by_taxid": SYNTHETIC_GOLD_SUPPORT,
        },
        "bioboxes": {
            "version": BIOBOXES_VERSION,
            "ranks": BIOBOXES_RANK_LIST,
            "taxonomy_fixture": {
                key: value.__dict__ for key, value in TAXONOMY_FIXTURE.items()
            },
        },
        "opal": {
            "real_opal_run": opal_ran,
            "status": opal_status,
            "command": opal_command,
            "command_shape": "opal.py -g GOLD_PROFILE -o OUTPUT_DIR -l LABEL1,LABEL2 profile1 profile2",
            "html_report": "opal_out/results.html",
            "index_html": "opal_out/index.html",
        },
        "warnings": warnings,
    }
    path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    """Generate prototype profiles, report, and QC metadata."""
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    assignments: list[ContigAssignment] = []
    warnings: list[str] = []
    blast_support: dict[str, int] = {}
    if external_profiles_requested(args):
        warnings = copy_external_profiles(args, output_dir)
    else:
        assignments, warnings = collapse_blast_contigs(SYNTHETIC_BLAST_FINAL_ROWS)
        blast_support = support_by_taxid(assignments)
        write_bioboxes_profile(
            output_dir / "gold.profile",
            sample_id=SAMPLE_ID,
            support=SYNTHETIC_GOLD_SUPPORT,
        )
        write_bioboxes_profile(
            output_dir / "blast_assembled_read_support.profile",
            sample_id=SAMPLE_ID,
            support=blast_support,
        )
        write_bioboxes_profile(
            output_dir / "sourmash_read_support.profile",
            sample_id=SAMPLE_ID,
            support=SYNTHETIC_SOURMASH_SUPPORT,
        )
    opal_result = run_opal_if_available(output_dir, opal_command=args.opal_command)
    write_qc(
        output_dir / "comparison.qc.json",
        assignments=assignments,
        blast_support=blast_support,
        opal_result=opal_result,
        warnings=warnings,
    )


if __name__ == "__main__":
    main()
