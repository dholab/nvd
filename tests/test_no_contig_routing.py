"""Focused Nextflow contracts for explicit no-contig routing."""

from __future__ import annotations

import gzip
import json
import os
import shutil
import sqlite3
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROCESS_CONTIGS = ROOT / "subworkflows" / "process_contigs"
PREPARE_BLAST_QUERIES = ROOT / "subworkflows" / "prepare_blast_queries"
SHORT_READ_ASSEMBLY = ROOT / "subworkflows" / "short_read_denovo_assembly"
LONG_READ_ENSEMBLE = ROOT / "subworkflows" / "long_read_denovo_ensembly"
FASTX_MODULE = ROOT / "modules" / "fastx"


def write_executable(path: Path, source: str) -> None:
    path.write_text(source, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def workflow_environment(bin_dir: Path) -> dict[str, str]:
    environment = os.environ.copy()
    environment["NXF_ANSI_LOG"] = "false"
    environment["PATH"] = (
        f"{bin_dir}{os.pathsep}{ROOT / 'bin'}{os.pathsep}{environment['PATH']}"
    )
    return environment


def run_nextflow(
    workflow: Path,
    *,
    bin_dir: Path,
    parameters: list[str] | None = None,
) -> subprocess.CompletedProcess[str]:
    nextflow = shutil.which("nextflow")
    assert nextflow is not None
    return subprocess.run(  # noqa: S603
        [
            nextflow,
            "-C",
            "/dev/null",
            "run",
            str(workflow),
            *(parameters or []),
        ],
        cwd=workflow.parent,
        env=workflow_environment(bin_dir),
        text=True,
        capture_output=True,
        check=False,
    )


def write_process_contig_fakes(bin_dir: Path) -> None:
    write_executable(
        bin_dir / "collect_contigs.py",
        """#!/usr/bin/env python3
import shutil
import sys
from pathlib import Path

def arg(name):
    return sys.argv[sys.argv.index(name) + 1]

source = Path(arg("--input-fasta"))
output = Path(arg("--output-fasta"))
if source.stat().st_size:
    shutil.copyfile(source, output)
else:
    output.touch()
Path(arg("--query-lookup")).touch()
""",
    )
    write_executable(
        bin_dir / "bbmask.sh",
        """#!/usr/bin/env python3
import shutil
import sys

values = dict(argument.split("=", 1) for argument in sys.argv[1:] if "=" in argument)
shutil.copyfile(values["in"], values["out"])
""",
    )
    write_executable(
        bin_dir / "reformat.sh",
        """#!/usr/bin/env python3
import shutil
import sys
from pathlib import Path

values = dict(argument.split("=", 1) for argument in sys.argv[1:] if "=" in argument)
output = Path(values["out"])
if "length_empty" in output.name:
    output.touch()
else:
    shutil.copyfile(values["in"], output)
""",
    )


def test_report_only_assembly_profile_failure_does_not_block_sentinel(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    attempts = tmp_path / "profile_attempts.log"
    write_executable(
        bin_dir / "profile_fastx.py",
        f"#!/bin/sh\nprintf 'attempt\\n' >> {json.dumps(str(attempts))}\nexit 2\n",
    )
    fasta = tmp_path / "contigs.fasta"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ PROFILE_ASSEMBLY_FASTA_FOR_REPORT }} from '{FASTX_MODULE}'

workflow {{
    PROFILE_ASSEMBLY_FASTA_FOR_REPORT(Channel.of(tuple('sample_A', 'illumina', 'single', 'spades', file('{fasta}'))))
    Channel.value('scientific-sentinel').view {{ value -> "SENTINEL: ${{value}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "SENTINEL: scientific-sentinel" in completed.stdout
    assert attempts.read_text(encoding="utf-8").count("attempt") == 3


def test_contig_processing_stops_at_the_first_empty_checkpoint(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_process_contig_fakes(bin_dir)

    assembly_empty = tmp_path / "assembly_empty.fasta"
    assembly_empty.touch()
    length_empty = tmp_path / "length_empty.fasta"
    length_empty.write_text(">contig\nACGT\n", encoding="utf-8")
    retained = tmp_path / "retained.fasta"
    retained.write_text(">contig\nACGT\n", encoding="utf-8")

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ PROCESS_CONTIGS }} from '{PROCESS_CONTIGS}'

params.entropy = 0.7
params.qtrim = 'f'
params.min_consecutive_bases = 1

workflow {{
    PROCESS_CONTIGS(Channel.of(
        tuple('assembly_empty', 'illumina', 'single', 'spades', file('{assembly_empty}')),
        tuple('length_empty', 'illumina', 'single', 'spades', file('{length_empty}')),
        tuple('retained', 'illumina', 'single', 'spades', file('{retained}')),
    ))

    PROCESS_CONTIGS.out.no_contigs.view {{ sample -> "NO_CONTIGS: ${{sample[0]}}" }}
    PROCESS_CONTIGS.out.contigs.view {{ sample -> "CONTIGS: ${{sample[0]}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    assert completed.returncode == 0, diagnostics
    assert completed.stdout.count("NO_CONTIGS: assembly_empty") == 1
    assert completed.stdout.count("NO_CONTIGS: length_empty") == 1
    assert completed.stdout.count("CONTIGS: retained") == 1
    assert "NO_CONTIGS: retained" not in completed.stdout
    assert "nvd.contig_route" not in completed.stdout

    log_text = (tmp_path / ".nextflow.log").read_text(encoding="utf-8")
    assert "sample_id=assembly_empty" in log_text
    assert "stage=collect_contigs" in log_text
    assert "sample_id=length_empty" in log_text
    assert "stage=length_filter" in log_text
    assert "MASK_LOW_COMPLEXITY (assembly_empty)" not in log_text
    assert "FILTER_SHORT_CONTIGS (assembly_empty)" not in log_text


def write_fastq(path: Path, read_name: str) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"@{read_name}\nACGT\n+\n!!!!\n")
    return path


def test_short_read_ineligibility_emits_a_no_contig_sample(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    reads = write_fastq(tmp_path / "reads.fastq.gz", "short")

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ SHORT_READ_DENOVO_ASSEMBLY }} from '{SHORT_READ_ASSEMBLY}'

workflow {{
    SHORT_READ_DENOVO_ASSEMBLY(Channel.of(tuple(
        [id: 'ineligible', platform: 'illumina', read_structure: 'single', sequence_count: 99],
        [[reads: file('{reads}')]],
    )))

    SHORT_READ_DENOVO_ASSEMBLY.out.no_contigs.view {{ sample ->
        "SHORT_INELIGIBLE: ${{sample[0]}}"
    }}
    SHORT_READ_DENOVO_ASSEMBLY.out.contigs.view {{ sample ->
        "UNEXPECTED_CONTIGS: ${{sample[0]}}"
    }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    assert completed.returncode == 0, diagnostics
    assert completed.stdout.count("SHORT_INELIGIBLE: ineligible") == 1
    assert "UNEXPECTED_CONTIGS:" not in completed.stdout


def write_deacon_fake(bin_dir: Path) -> None:
    write_executable(
        bin_dir / "deacon",
        """#!/usr/bin/env python3
import json
import sys
from pathlib import Path

output = Path(sys.argv[sys.argv.index("--output") + 1])
output.touch()
summary = Path(sys.argv[sys.argv.index("--summary") + 1])
summary.write_text(json.dumps({
    "seqs_in": 1,
    "seqs_out": 0,
    "seqs_removed": 1,
    "bp_in": 4,
    "bp_out": 0,
    "bp_removed": 4,
}))
""",
    )


def test_no_contig_samples_route_full_reads_without_mapback(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_deacon_fake(bin_dir)

    lib = tmp_path / "lib"
    lib.mkdir()
    shutil.copy2(ROOT / "lib" / "NvdUtils.groovy", lib / "NvdUtils.groovy")

    assembled = tmp_path / "deacon_empty.short_filtered.fasta"
    assembled.write_text(">contig\nACGT\n", encoding="utf-8")
    query_lookup = tmp_path / "deacon_empty.query_sequences.sqlite"
    with sqlite3.connect(query_lookup):
        pass
    target_index = tmp_path / "target.idx"
    target_index.touch()
    depletion_index = tmp_path / "depletion.idx"
    depletion_index.touch()

    merged = write_fastq(tmp_path / "upstream_empty.merged.fastq.gz", "merged")
    single = write_fastq(tmp_path / "upstream_empty.single.fastq.gz", "single")
    deacon_single = write_fastq(tmp_path / "deacon_empty.single.fastq.gz", "deacon")

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ PREPARE_BLAST_QUERIES }} from '{PREPARE_BLAST_QUERIES}'

params.experimental = true
params.dedup = false
params.dedup_seq = false
params.dedup_pos = false
params.no_enrichment = true
params.min_consecutive_bases = 1
params.virus_abs_threshold = 1
params.virus_rel_threshold = 0.0
params.host_abs_threshold = 1
params.host_rel_threshold = 0.0

workflow {{
    PREPARE_BLAST_QUERIES(
        Channel.of(tuple(
            'deacon_empty',
            'illumina',
            'single',
            file('{assembled}'),
            file('{query_lookup}'),
        )),
        Channel.of(tuple('upstream_empty', 'illumina')),
        Channel.of(tuple(
            'upstream_empty',
            'illumina',
            [file('{merged}')],
            [file('{single}')],
        )),
        Channel.of(tuple(
            'deacon_empty',
            'illumina',
            'single',
            file('{deacon_single}'),
        )),
        Channel.value(file('{target_index}')),
        Channel.value(tuple(false, file('{depletion_index}'))),
    )

    PREPARE_BLAST_QUERIES.out.queries.view {{ query ->
        "QUERY: ${{query[0]}}:${{query[2]}}"
    }}
    PREPARE_BLAST_QUERIES.out.contig_read_counts.view {{ counts ->
        "COUNTS: ${{counts[0]}}:${{counts[1].size()}}"
    }}
    PREPARE_BLAST_QUERIES.out.no_contigs.view {{ sample ->
        "NO_CONTIGS: ${{sample[0]}}"
    }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(
        workflow,
        bin_dir=bin_dir,
        parameters=[
            "--experimental",
            "true",
            "--dedup",
            "false",
            "--dedup_seq",
            "false",
            "--dedup_pos",
            "false",
            "--no_enrichment",
            "true",
            "--min_consecutive_bases",
            "1",
            "--virus_abs_threshold",
            "1",
            "--virus_rel_threshold",
            "0.0",
            "--host_abs_threshold",
            "1",
            "--host_rel_threshold",
            "0.0",
        ],
    )
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    assert completed.returncode == 0, diagnostics
    assert completed.stdout.count("QUERY: upstream_empty:overlap_merged_pair") == 1
    assert completed.stdout.count("QUERY: upstream_empty:single_read") == 1
    assert completed.stdout.count("QUERY: deacon_empty:single_read") == 1
    assert completed.stdout.count("COUNTS: upstream_empty:0") == 1
    assert completed.stdout.count("COUNTS: deacon_empty:0") == 1
    assert completed.stdout.count("NO_CONTIGS: upstream_empty") == 1
    assert completed.stdout.count("NO_CONTIGS: deacon_empty") == 1
    assert "nvd.contig_route" not in completed.stdout

    log_text = (tmp_path / ".nextflow.log").read_text(encoding="utf-8")
    assert "sample_id=deacon_empty" in log_text
    assert "stage=deacon_filter" in log_text
    assert "MAP_SINGLE_READS (deacon_empty)" not in log_text
    assert "SELECT_BLAST_QUERIES (deacon_empty)" not in log_text
    assert list(tmp_path.rglob("*.bam")) == []
    assert list(tmp_path.rglob("*.bam.bai")) == []
    assert list(tmp_path.rglob("*_mapped_counts.txt")) == []
    assert list(tmp_path.rglob("*.crumbs.coverage.tsv")) == []


def write_long_read_profile(path: Path, *, eligible: bool) -> Path:
    threshold_count = 1 if eligible else 0
    path.write_text(
        json.dumps(
            {
                "sample_id": path.stem,
                "sequence_count": 1,
                "total_bases": 1000,
                "length": {"min": 1000, "max": 1000, "mean": 1000},
                "thresholds": [
                    {
                        "name": name,
                        "axis": "length",
                        "value": 1000,
                        "sequence_count_at_or_above": (
                            threshold_count if name == "myloasm_min_read_length" else 0
                        ),
                        "bases_at_or_above": (
                            1000
                            if eligible and name == "myloasm_min_read_length"
                            else 0
                        ),
                    }
                    for name in (
                        "metamdbg_min_read_overlap",
                        "myloasm_min_read_length",
                        "metaflye_min_read_length",
                    )
                ],
            },
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def test_long_read_query_candidates_survive_zero_available_assembler_contigs(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_executable(
        bin_dir / "myloasm",
        """#!/bin/sh
case "$1" in
    *assembled*)
        mkdir -p output
        printf '>contig\nACGT\n' > output/assembly_primary.fa
        ;;
    *)
        exit 1
        ;;
esac
""",
    )
    reads = write_fastq(tmp_path / "reads.fastq.gz", "long")
    assembled_reads = write_fastq(tmp_path / "assembled.fastq.gz", "assembled")
    ineligible = write_long_read_profile(tmp_path / "ineligible.json", eligible=False)
    eligible = write_long_read_profile(tmp_path / "eligible.json", eligible=True)
    assembled = write_long_read_profile(tmp_path / "assembled.json", eligible=True)

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ LONG_READ_DENOVO_ENSEMBLY }} from '{LONG_READ_ENSEMBLE}'

workflow {{
    LONG_READ_DENOVO_ENSEMBLY(Channel.of(
        tuple(
            [id: 'ineligible', platform: 'ont', read_structure: 'single'],
            file('{reads}'),
            file('{ineligible}'),
            file('{ineligible}'),
        ),
        tuple(
            [id: 'eligible', platform: 'ont', read_structure: 'single'],
            file('{reads}'),
            file('{eligible}'),
            file('{eligible}'),
        ),
        tuple(
            [id: 'assembled', platform: 'ont', read_structure: 'single'],
            file('{assembled_reads}'),
            file('{assembled}'),
            file('{assembled}'),
        ),
    ))

    LONG_READ_DENOVO_ENSEMBLY.out.no_contigs.view {{ sample ->
        "READS_ONLY: ${{sample[0]}}"
    }}
    LONG_READ_DENOVO_ENSEMBLY.out.contigs.view {{ sample ->
        "CONTIGS: ${{sample[0]}}"
    }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(
        workflow,
        bin_dir=bin_dir,
        parameters=["--experimental", "true"],
    )
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"

    assert completed.returncode == 0, diagnostics
    assert completed.stdout.count("READS_ONLY: ineligible") == 1
    assert completed.stdout.count("READS_ONLY: eligible") == 1
    assert completed.stdout.count("CONTIGS: assembled") == 1
    assert "READS_ONLY: assembled" not in completed.stdout
