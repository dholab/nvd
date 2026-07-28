from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def write_executable(path: Path, source: str) -> None:
    path.write_text(source, encoding="utf-8")
    path.chmod(0o755)


def write_streaming_tool_fakes(
    bin_dir: Path,
    *,
    fail_nuclease: bool,
) -> tuple[Path, Path, Path]:
    bin_dir.mkdir()
    nuclease_args = bin_dir.parent / "nuclease.args"
    nuclease_calls = bin_dir.parent / "nuclease.calls"
    deacon_args = bin_dir.parent / "deacon.args"
    nuclease_exit = 42 if fail_nuclease else 0
    write_executable(
        bin_dir / "nuclease",
        f"""#!/usr/bin/env bash
set -eu

printf x >> {nuclease_calls}
printf '%s ' "$@" > {nuclease_args}
if (( {nuclease_exit} != 0 )); then
    printf '%s' 'simulated nuclease ENA failure' >&2
    exit {nuclease_exit}
fi
summary=''
while (( $# )); do
    case $1 in
        --summary)
            summary=$2
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
printf '%s' '{{"reads_seen": 2, "reads_emitted": 2}}' > "$summary"
printf '%s' '@read/1
ACGT
+
!!!!
@read/2
TGCA
+
!!!!
'
""",
    )
    write_executable(
        bin_dir / "deacon",
        f"""#!/usr/bin/env bash
set -eu

printf '%s ' "$@" > {deacon_args}
output=''
summary=''
while (( $# )); do
    case $1 in
        --output)
            output=$2
            shift 2
            ;;
        --summary)
            summary=$2
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done
cat > "$output"
printf '%s' '{{"seqs_in": 2}}' > "$summary"
""",
    )
    return nuclease_args, nuclease_calls, deacon_args


def run_streaming_process(
    tmp_path: Path,
    *,
    fail_nuclease: bool = False,
) -> tuple[subprocess.CompletedProcess[str], Path, Path, Path, Path]:
    bin_dir = tmp_path / "bin"
    nuclease_args, nuclease_calls, deacon_args = write_streaming_tool_fakes(
        bin_dir,
        fail_nuclease=fail_nuclease,
    )
    index = tmp_path / "target.idx"
    index.write_text("index", encoding="utf-8")
    output_dir = tmp_path / "outputs"
    workflow = tmp_path / "stream_sra.nf"
    workflow.write_text(
        f"""
nextflow.enable.dsl = 2

include {{ DEACON_ENRICH_SRA_READS }} from '{(ROOT / "modules" / "deacon").as_posix()}'

workflow {{
    ch_inputs = channel
        .of(tuple('sample_A', 'illumina', 'SRR_TEST'))
        .combine(channel.fromPath(params.index))
        .combine(channel.value(true))
    DEACON_ENRICH_SRA_READS(ch_inputs)
}}
""",
        encoding="utf-8",
    )
    config = tmp_path / "nextflow.config"
    config.write_text(
        """
params.max_concurrent_downloads = 1
params.virus_abs_threshold = 1
params.virus_rel_threshold = 0.0
process {
    executor = 'local'
    cpus = 2
    memory = 2.GB
    withName: 'DEACON_ENRICH_SRA_READS' {
        publishDir = [path: params.output_dir, mode: 'copy']
    }
}
""",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{ROOT / 'bin'}{os.pathsep}{env['PATH']}"
    completed = subprocess.run(  # noqa: S603
        [  # noqa: S607
            "nextflow",
            "run",
            str(workflow),
            "-c",
            str(config),
            "-work-dir",
            str(tmp_path / "work"),
            "--index",
            str(index),
            "--output_dir",
            str(output_dir),
        ],
        cwd=tmp_path,
        env=env,
        check=False,
        text=True,
        capture_output=True,
        timeout=60,
    )
    return completed, output_dir, nuclease_args, nuclease_calls, deacon_args


def test_gather_reads_routes_sra_accessions_without_downloading(
    tmp_path: Path,
) -> None:
    local_fastq = tmp_path / "local.fastq"
    local_fastq.write_text("@local\nACGT\n+\n!!!!\n", encoding="utf-8")
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text(
        "sample_id,srr,platform,fastq1,fastq2\n"
        f"local_sample,,illumina,{local_fastq},\n"
        "sra_sample,SRR_TEST,illumina,,\n",
        encoding="utf-8",
    )
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    unexpected_download = tmp_path / "unexpected-sracha-call"
    write_executable(
        bin_dir / "sracha",
        f"""#!/usr/bin/env bash
touch {unexpected_download}
exit 99
""",
    )
    workflow = tmp_path / "gather_streamed_sra.nf"
    workflow.write_text(
        f"""
nextflow.enable.dsl = 2

include {{ GATHER_READS }} from '{(ROOT / "subworkflows" / "gather_reads").as_posix()}'

workflow {{
    GATHER_READS(channel.fromPath(params.samplesheet))
    GATHER_READS.out.reads.view {{ meta, _reads -> "READ:${{meta.id}}:${{meta.source}}" }}
    GATHER_READS.out.sra_accessions.view {{ id, platform, accession -> "SRA:${{id}}:${{platform}}:${{accession}}" }}
}}
""",
        encoding="utf-8",
    )
    config = tmp_path / "gather.config"
    config.write_text(
        """
params.stream_sra = true
params.max_concurrent_downloads = 1
process.executor = 'local'
process.cpus = 1
process.memory = 1.GB
""",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{ROOT / 'bin'}{os.pathsep}{env['PATH']}"

    completed = subprocess.run(  # noqa: S603
        [  # noqa: S607
            "nextflow",
            "run",
            str(workflow),
            "-c",
            str(config),
            "-work-dir",
            str(tmp_path / "work"),
            "--samplesheet",
            str(samplesheet),
        ],
        cwd=tmp_path,
        env=env,
        check=False,
        text=True,
        capture_output=True,
        timeout=60,
    )

    diagnostics = completed.stdout + completed.stderr
    assert completed.returncode == 0, diagnostics
    assert "READ:local_sample:single_file" in diagnostics
    assert "SRA:sra_sample:illumina:SRR_TEST" in diagnostics
    assert not unexpected_download.exists()


@pytest.mark.parametrize(
    ("mode_args", "expected_messages"),
    [
        (
            ["--experimental", "false", "--skip_fastqc", "true"],
            [
                "ENA streaming for SRA run accessions is available only when experimental features are enabled",
                "experimental = true",
                "stream_sra = false",
            ],
        ),
        (
            ["--experimental", "true", "--skip_fastqc", "false"],
            [
                "ENA streaming cannot run with raw FastQC enabled",
                "skip_fastqc = true",
                "stream_sra = false",
            ],
        ),
    ],
)
def test_sra_streaming_rejects_incompatible_modes_with_corrective_message(
    tmp_path: Path,
    mode_args: list[str],
    expected_messages: list[str],
) -> None:
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text(
        "sample_id,srr,platform,fastq1,fastq2\n",
        encoding="utf-8",
    )

    completed = subprocess.run(  # noqa: S603
        [  # noqa: S607
            "nextflow",
            "run",
            ".",
            "-profile",
            "test",
            "--samplesheet",
            str(samplesheet),
            "--stream_sra",
            "true",
            *mode_args,
            "--skip_blast",
            "true",
            "--skip_assembly",
            "true",
        ],
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=True,
        timeout=30,
    )

    diagnostics = completed.stdout + completed.stderr
    assert completed.returncode != 0
    for message in expected_messages:
        assert message in diagnostics


def test_nuclease_passthrough_streams_interleaved_ena_reads_into_deacon(
    tmp_path: Path,
) -> None:
    completed, output_dir, nuclease_args, nuclease_calls, deacon_args = (
        run_streaming_process(tmp_path)
    )

    diagnostics = completed.stdout + completed.stderr
    assert completed.returncode == 0, diagnostics
    assert (output_dir / "sample_A.sra_read_structure.txt").read_text(
        encoding="utf-8",
    ).strip() == "interleaved"
    assert "@read/1" in (output_dir / "sample_A.target_enriched.fastq.gz").read_text(
        encoding="utf-8",
    )
    assert (output_dir / "sample_A.deacon_filter.json").is_file()
    assert (output_dir / "sample_A.nuclease.json").is_file()
    assert nuclease_calls.read_text(encoding="utf-8") == "x"

    argv = nuclease_args.read_text(encoding="utf-8").split()
    assert argv == [
        "--ena",
        "SRR_TEST",
        "--passthrough",
        "--output-encoding",
        "plain",
        "--summary",
        "sample_A.nuclease.json",
    ]
    assert deacon_args.read_text(encoding="utf-8").split()[-2:] == ["-", "-"]
    assert not list((tmp_path / "work").glob("**/*.sra"))


def test_nuclease_ena_failure_retries_then_ignores_sample(tmp_path: Path) -> None:
    completed, output_dir, _nuclease_args, nuclease_calls, _deacon_args = (
        run_streaming_process(tmp_path, fail_nuclease=True)
    )

    diagnostics = completed.stdout + completed.stderr
    assert completed.returncode == 0, diagnostics
    task_errors = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (tmp_path / "work").glob("**/.command.err")
    )
    assert "simulated nuclease ENA failure" in task_errors
    assert nuclease_calls.read_text(encoding="utf-8") == "xxx"
    assert not (output_dir / "sample_A.target_enriched.fastq.gz").exists()
