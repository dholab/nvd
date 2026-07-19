"""Focused Nextflow contracts for the NVD MultiQC reporting spine."""

from __future__ import annotations

import json
import os
import shutil
import stat
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
BUNDLING_SUBWORKFLOW = ROOT / "subworkflows" / "multiqc_bundling"
FASTQC_MODULE = ROOT / "modules" / "fastqc"
MULTIQC_MODULE = ROOT / "modules" / "multiqc"
EXPECTED_EXPANDED_UNITS = 7
EXPECTED_PAIRED_UNITS_PER_END = 2
EXPECTED_RENDERED_FASTQC_UNITS = 4


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
    isolated_config: bool = True,
) -> subprocess.CompletedProcess[str]:
    nextflow = shutil.which("nextflow")
    assert nextflow is not None
    command = [nextflow]
    if isolated_config:
        command.extend(["-C", "/dev/null"])
    command.extend(["run", str(workflow)])
    return subprocess.run(  # noqa: S603
        command,
        cwd=workflow.parent,
        env=workflow_environment(bin_dir),
        text=True,
        capture_output=True,
        check=False,
    )


def copy_reporting_lib(tmp_path: Path) -> None:
    lib = tmp_path / "lib"
    lib.mkdir(exist_ok=True)
    shutil.copy2(ROOT / "lib" / "NvdReporting.groovy", lib / "NvdReporting.groovy")


def copy_results_config(tmp_path: Path) -> None:
    conf = tmp_path / "conf"
    conf.mkdir()
    shutil.copy2(ROOT / "conf" / "results.config", conf / "results.config")
    (tmp_path / "nextflow.config").write_text(
        f"""\
params.results = '{tmp_path / "results"}'
params.experimental = false
params.no_enrichment = true
includeConfig 'conf/results.config'
""",
        encoding="utf-8",
    )


def write_fastq(path: Path) -> Path:
    path.write_text("@read\nACGT\n+\n!!!!\n", encoding="utf-8")
    return path


def write_tool_fakes(
    bin_dir: Path,
    *,
    fastqc_fails: bool = False,
    fake_multiqc: bool = True,
    fastqc_attempt_log: Path | None = None,
    multiqc_fails: bool = False,
    multiqc_attempt_log: Path | None = None,
) -> None:
    if fastqc_fails:
        log_line = (
            f"printf 'fastqc\\n' >> {json.dumps(str(fastqc_attempt_log))}\n"
            if fastqc_attempt_log is not None
            else ""
        )
        write_executable(bin_dir / "fastqc", f"#!/bin/sh\n{log_line}exit 2\n")
    else:
        write_executable(
            bin_dir / "fastqc",
            """#!/usr/bin/env python3
import sys
from pathlib import Path
from zipfile import ZipFile

outdir = Path(sys.argv[sys.argv.index('--outdir') + 1])
outdir.mkdir(exist_ok=True)
input_path = Path(sys.argv[-1])
name = input_path.name
for suffix in ('.gz', '.fastq', '.fq'):
    if name.endswith(suffix):
        name = name[: -len(suffix)]
base = outdir / f'{name}_fastqc'
fastqc_data = f'''##FastQC\t0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Filename\t{input_path.name}
File type\tConventional base calls
Encoding\tSanger / Illumina 1.9
Total Sequences\t1
Sequences flagged as poor quality\t0
Sequence length\t4
%GC\t50
>>END_MODULE
'''
with ZipFile(str(base) + '.zip', 'w') as archive:
    archive.writestr(f'{base.name}/', '')
    archive.writestr(f'{base.name}/fastqc_data.txt', fastqc_data)
    archive.writestr(f'{base.name}/summary.txt', f'PASS\tBasic Statistics\t{input_path.name}\\n')
Path(str(base) + '.html').write_text('<html>fastqc</html>')
""",
        )
    if multiqc_fails:
        log_line = (
            f"printf 'multiqc\\n' >> {json.dumps(str(multiqc_attempt_log))}\n"
            if multiqc_attempt_log is not None
            else ""
        )
        write_executable(bin_dir / "multiqc", f"#!/bin/sh\n{log_line}exit 2\n")
    elif fake_multiqc:
        write_executable(
            bin_dir / "multiqc",
            """#!/usr/bin/env python3
from pathlib import Path

Path('multiqc_report.html').write_text('<html>NVD MultiQC</html>')
data = Path('multiqc_data')
data.mkdir(exist_ok=True)
(data / 'multiqc_data.json').write_text('{}')
""",
        )


def write_roster(
    path: Path,
    *,
    sample_id: str = "sample_A",
    source: str = "paired_files",
    read_structure: str = "paired",
    read_count: int = 1,
) -> Path:
    record = {
        "sample_id": sample_id,
        "platform": "illumina",
        "source": source,
    }
    if read_structure == "single":
        record["read_structure"] = "single"
        record["read_counts"] = {"single": read_count}
        record["reads"] = [
            f"/fake/private/single_{index}.fastq" for index in range(read_count)
        ]
    else:
        record["read_structure"] = "paired"
        record["read_counts"] = {"R1": read_count, "R2": read_count}
        record["r1"] = [
            f"/fake/private/R1_{index}.fastq" for index in range(read_count)
        ]
        record["r2"] = [
            f"/fake/private/R2_{index}.fastq" for index in range(read_count)
        ]
    path.write_text(
        json.dumps(record) + "\n",
        encoding="utf-8",
    )
    return path


def write_query_batch_summary(path: Path) -> Path:
    path.write_text(
        "sample_id\tplatform\tquery_class\tquery_source\tn_query_sequences\tquery_fasta_present\tquery_lookup_present\n"
        "sample_A\tillumina\tshort_assembly_contig\tcontig\t2\ttrue\ttrue\n"
        "sample_A\tillumina\tlong_assembly_contig\tcontig\t0\tfalse\tfalse\n"
        "sample_A\tillumina\toverlap_merged_pair\tread_query\t0\tfalse\tfalse\n"
        "sample_A\tillumina\tsingle_read\tread_query\t0\tfalse\tfalse\n",
        encoding="utf-8",
    )
    return path


def nextflow_file(path: Path) -> str:
    escaped = str(path).replace("\\", "\\\\").replace("'", "\\'")
    return f"file('{escaped}')"


def empty_domain_inputs(config: Path) -> str:
    channels = ["Channel.value(false)" for _ in range(3)]
    channels.extend("Channel.empty()" for _ in range(10))
    channels.append(f"Channel.value(file('{config}'))")
    return ",\n        ".join(channels) + ","


def write_version(path: Path) -> Path:
    path.write_text("version=3.3.0\nrevision=test-rev\n", encoding="utf-8")
    return path


def read_text_tree(root: Path) -> str:
    return "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in sorted(item for item in root.rglob("*") if item.is_file())
    )


def read_native_fastqc_data(data_dir: Path) -> str:
    native_fastqc = data_dir / "multiqc_fastqc.txt"
    if not native_fastqc.is_file():
        raise AssertionError(f"No retained native FastQC data found: {native_fastqc}")
    return native_fastqc.read_text(encoding="utf-8")


def test_fastqc_unit_expansion_for_supported_bundle_shapes(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    single = write_fastq(tmp_path / "single.fastq")
    r1_a = write_fastq(tmp_path / "lane.fastq")
    repeated = tmp_path / "repeated"
    repeated.mkdir()
    r1_b = write_fastq(repeated / "lane.fastq")
    r2_a = write_fastq(tmp_path / "mate.fastq")
    r2_b = write_fastq(repeated / "mate.fastq")
    sra_r1 = write_fastq(tmp_path / "sra_1.fastq")
    sra_r2 = write_fastq(tmp_path / "sra_2.fastq")

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

workflow {{
    Channel.of(
        tuple([id: 'single_sample', platform: 'ont', source: 'single_file', read_mode: 'single', r1_count: 1], [{nextflow_file(single)}]),
        tuple([id: 'paired_sample', platform: 'illumina', source: 'paired_globs', read_mode: 'paired', r1_count: 2], [{nextflow_file(r1_a)}, {nextflow_file(r1_b)}, {nextflow_file(r2_a)}, {nextflow_file(r2_b)}]),
        tuple([id: 'sra_sample', platform: 'illumina', source: 'sra', read_mode: 'paired', r1_count: 1], [{nextflow_file(sra_r1)}, {nextflow_file(sra_r2)}]),
    ).flatMap {{ meta, reads -> NvdReporting.processReadyFastqcTuples(meta, reads) }}
     .view {{ planned -> "UNIT: ${{planned[0].sample_id}}:${{planned[0].source}}:${{planned[0].read_end}}:${{planned[0].input_ordinal}}:${{planned[0].alias}}:${{planned[0].staged_filename}}:${{planned[0].package_name}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    unit_lines = [
        line for line in completed.stdout.splitlines() if line.startswith("UNIT: ")
    ]
    assert len(unit_lines) == EXPECTED_EXPANDED_UNITS
    assert any(
        ":single_file:single:1:single_sample.raw.single.001:single_sample.raw.single.001.fastq:single_sample.raw.single.001.fastqc"
        in line
        for line in unit_lines
    )
    assert (
        sum(":paired_globs:R1:" in line for line in unit_lines)
        == EXPECTED_PAIRED_UNITS_PER_END
    )
    assert (
        sum(":paired_globs:R2:" in line for line in unit_lines)
        == EXPECTED_PAIRED_UNITS_PER_END
    )
    assert sum(":sra:R1:1:" in line for line in unit_lines) == 1
    assert sum(":sra:R2:1:" in line for line in unit_lines) == 1
    aliases = [line.split(":")[-3] for line in unit_lines]
    assert len(aliases) == len(set(aliases))


def test_target_enrichment_package_crosses_explicit_process_boundary(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    write_tool_fakes(bin_dir)
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")
    config = ROOT / "assets" / "multiqc_config.yaml"
    stats = tmp_path / "deacon.json"
    stats.write_text(
        json.dumps(
            {
                "version": "deacon-test",
                "seqs_in": 10,
                "seqs_out": 8,
                "seqs_removed": 2,
                "seqs_out_proportion": 0.8,
                "seqs_removed_proportion": 0.2,
                "bp_in": 1000,
                "bp_out": 800,
                "bp_removed": 200,
                "bp_out_proportion": 0.8,
                "bp_removed_proportion": 0.2,
            },
        ),
        encoding="utf-8",
    )
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'

params.results = '{tmp_path / "results"}'
params.experimental = false
workflow {{
    MULTIQC_BUNDLING(
        Channel.empty(),
        Channel.empty(),
        Channel.value(file('{roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        Channel.value(true),
        Channel.value(false),
        Channel.value(false),
        Channel.of(tuple('sample_A', {nextflow_file(stats)})),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.value(file('{config}')),
    )
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    manifests = sorted(tmp_path.glob("work/**/nvd_report_manifest.json"))
    assert manifests, diagnostics
    manifest = json.loads(manifests[-1].read_text(encoding="utf-8"))
    assert len(manifest["report_packages"]) == 1
    assert "nvd_target_enrichment" in manifest["sections"]
    sections = sorted(tmp_path.glob("work/**/nvd_target_enrichment_mqc.yaml"))
    assert sections, diagnostics
    section_text = sections[-1].read_text(encoding="utf-8")
    assert "status: observed" in section_text
    assert "reads_in: 10" in section_text
    assert "\n  sample_A:\n" in section_text
    assert "row_0001" not in section_text
    assert "sample_id:" not in section_text


def test_prepared_query_summary_crosses_explicit_package_boundary(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")
    summary = write_query_batch_summary(tmp_path / "sample_A.blast_query_batches.tsv")
    config = ROOT / "assets" / "multiqc_config.yaml"
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'

params.results = '{tmp_path / "results"}'
params.experimental = false
workflow {{
    MULTIQC_BUNDLING(
        Channel.empty(),
        Channel.empty(),
        Channel.value(file('{roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        Channel.value(false),
        Channel.value(false),
        Channel.value(true),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.of(tuple('sample_A', file('{summary}'))),
        Channel.value(file('{config}')),
    )
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    manifests = sorted(tmp_path.glob("work/**/nvd_report_manifest.json"))
    assert manifests, diagnostics
    manifest = json.loads(manifests[-1].read_text(encoding="utf-8"))
    assert manifest["report_packages"] == [
        {
            "schema_version": "nvd.report-package/v1",
            "domain": "query_preparation",
            "artifact_type": "prepared_query_batches",
            "sample_id": "sample_A",
            "summary": "summary.payload",
        },
    ]
    sections = sorted(
        tmp_path.glob("work/**/nvd_prepared_blast_query_batches_mqc.yaml"),
    )
    assert sections, diagnostics
    section = yaml.safe_load(sections[-1].read_text(encoding="utf-8"))
    assert section["data"]["sample_A | short_assembly_contig"][
        "query_sequences"
    ] == 2
    assert section["data"]["sample_A | single_read"]["availability"] == (
        "disabled"
    )


def test_assembly_decision_without_profile_is_reported_as_absent(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    write_tool_fakes(bin_dir)
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")
    config = ROOT / "assets" / "multiqc_config.yaml"
    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'

params.results = '{tmp_path / "results"}'
params.experimental = false
workflow {{
    MULTIQC_BUNDLING(
        Channel.empty(),
        Channel.empty(),
        Channel.value(file('{roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        Channel.value(false),
        Channel.value(false),
        Channel.value(true),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.of(tuple(
            [id: 'sample_A', producer: 'spades'],
            'run',
            120,
            100,
            'short_read_minimum_sequence_count',
        )),
        Channel.empty(),
        Channel.empty(),
        Channel.empty(),
        Channel.value(file('{config}')),
    )
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    sections = sorted(tmp_path.glob("work/**/nvd_assembly_mqc.yaml"))
    assert sections, diagnostics
    section_text = sections[-1].read_text(encoding="utf-8")
    assert "status: expected_output_absent" in section_text
    assert "eligibility_decision: run" in section_text
    assert "\n  sample_A | spades:\n" in section_text
    assert "row_0001" not in section_text
    assert "sample_id:" not in section_text


def test_physical_unit_aliases_survive_staging_and_rendering(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    copy_results_config(tmp_path)
    write_tool_fakes(bin_dir, fake_multiqc=False)
    r1_a = write_fastq(tmp_path / "sample_A.raw.R1.001.fastq")
    (tmp_path / "repeated").mkdir()
    r1_b = write_fastq(tmp_path / "repeated" / "lane's.fastq")
    write_fastq(r1_b)
    r2_a = write_fastq(tmp_path / "$(touch command_substitution_side_effect).fastq")
    r2_b = write_fastq(tmp_path / "repeated" / "mate.fastq")
    roster = write_roster(
        tmp_path / "resolved_reads.jsonl",
        source="paired_globs",
        read_count=2,
    )
    version = write_version(tmp_path / "nvd_version.txt")
    config = ROOT / "assets" / "multiqc_config.yaml"
    poison_dir = tmp_path / "results" / "nvd"
    poison_dir.mkdir(parents=True)
    poison_id = "poison_unique_nvd_regression"
    (poison_dir / "poison_mqc.yaml").write_text(
        "\n".join(
            [
                f"id: {poison_id}",
                "section_name: Poison Custom Content",
                "plot_type: table",
                "data:",
                "  poison:",
                "    value: should-not-be-scanned",
                "",
            ],
        ),
        encoding="utf-8",
    )

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'
include {{ GENERATE_MULTIQC_REPORT }} from '{MULTIQC_MODULE}'
include {{ FASTQC_RAW }} from '{FASTQC_MODULE}'

params.results = '{tmp_path / "results"}'
params.experimental = false

workflow {{
    reads = Channel.of(tuple(
        [id: 'sample_A', platform: 'illumina', source: 'paired_globs', read_mode: 'paired', r1_count: 2],
        [{nextflow_file(r1_a)}, {nextflow_file(r1_b)}, {nextflow_file(r2_a)}, {nextflow_file(r2_b)}],
    ))
    fastqc_units = reads.flatMap {{ meta, read_files -> NvdReporting.processReadyFastqcTuples(meta, read_files) }}
    FASTQC_RAW(fastqc_units)
    MULTIQC_BUNDLING(
        FASTQC_RAW.out.packages,
        FASTQC_RAW.out.zips,
        Channel.value(file('{roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        {empty_domain_inputs(config)}
    )
    GENERATE_MULTIQC_REPORT(
        MULTIQC_BUNDLING.out.fastqc_zips,
        MULTIQC_BUNDLING.out.inputs,
        MULTIQC_BUNDLING.out.config,
    )
    GENERATE_MULTIQC_REPORT.out.report.view {{ report -> "REPORT: ${{report.name}}" }}
    GENERATE_MULTIQC_REPORT.out.data.view {{ data -> "DATA: ${{data.name}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir, isolated_config=False)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "REPORT: multiqc_report.html" in completed.stdout
    assert "DATA: multiqc_data" in completed.stdout
    assert (tmp_path / "results" / "nvd" / "multiqc_report.html").is_file()
    assert not (tmp_path / "command_substitution_side_effect").exists()
    assert not sorted(tmp_path.glob("work/**/command_substitution_side_effect"))

    manifests = sorted(tmp_path.glob("work/**/nvd_report_manifest.json"))
    assert manifests, diagnostics
    manifest = json.loads(manifests[-1].read_text(encoding="utf-8"))
    aliases = [unit["alias"] for unit in manifest["raw_fastqc"]]
    assert aliases == sorted(aliases)
    assert len(aliases) == EXPECTED_RENDERED_FASTQC_UNITS
    assert len(set(aliases)) == EXPECTED_RENDERED_FASTQC_UNITS
    assert all(alias.startswith("sample_A.raw.") for alias in aliases)
    assert aliases == [
        "sample_A.raw.R1.001",
        "sample_A.raw.R1.002",
        "sample_A.raw.R2.001",
        "sample_A.raw.R2.002",
    ]
    assert {unit["read_end"] for unit in manifest["raw_fastqc"]} == {"R1", "R2"}
    published_fastqc = (
        tmp_path
        / "results"
        / "nvd"
        / "00_input_preparation"
        / "raw_fastq_qc"
        / "fastqc"
    )
    assert sorted(path.name for path in published_fastqc.glob("*_fastqc.zip")) == [
        f"{alias}_fastqc.zip" for alias in aliases
    ]
    assert sorted(path.name for path in published_fastqc.glob("*_fastqc.html")) == [
        f"{alias}_fastqc.html" for alias in aliases
    ]
    assert poison_id not in json.dumps(manifest)

    data_dirs = sorted(tmp_path.glob("work/**/multiqc_data"))
    assert data_dirs, diagnostics
    data_text = read_text_tree(data_dirs[-1])
    native_fastqc_data = read_native_fastqc_data(data_dirs[-1])
    assert str(r1_a) not in data_text
    assert str(r1_b) not in data_text
    assert str(r2_a) not in data_text
    assert str(r2_b) not in data_text
    for alias in aliases:
        assert alias in native_fastqc_data
    assert (data_dirs[-1] / "nvd_inputs" / "multiqc_config.yaml").is_file()
    assert poison_id not in data_text
    assert data_text.count("id: nvd_sample_roster") == 1
    assert data_text.count("NVD: 3.3.0 (test-rev)") == 1


def test_ignored_fastqc_failure_still_renders_roster_report(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    fastqc_attempt_log = tmp_path / "fastqc_attempts.log"
    write_tool_fakes(bin_dir, fastqc_fails=True, fastqc_attempt_log=fastqc_attempt_log)
    read = write_fastq(tmp_path / "single.fastq")
    roster = write_roster(
        tmp_path / "resolved_reads.jsonl",
        source="single_file",
        read_structure="single",
    )
    version = write_version(tmp_path / "nvd_version.txt")
    config = ROOT / "assets" / "multiqc_config.yaml"

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'
include {{ GENERATE_MULTIQC_REPORT }} from '{MULTIQC_MODULE}'
include {{ FASTQC_RAW }} from '{FASTQC_MODULE}'

params.results = '{tmp_path / "results"}'
params.experimental = false
workflow {{
    reads = Channel.of(tuple(
        [id: 'sample_A', platform: 'illumina', source: 'single_file', read_mode: 'single', r1_count: 1],
        [{nextflow_file(read)}],
    ))
    fastqc_units = reads.flatMap {{ meta, read_files -> NvdReporting.processReadyFastqcTuples(meta, read_files) }}
    FASTQC_RAW(fastqc_units)
    MULTIQC_BUNDLING(
        FASTQC_RAW.out.packages,
        FASTQC_RAW.out.zips,
        Channel.value(file('{roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        {empty_domain_inputs(config)}
    )
    GENERATE_MULTIQC_REPORT(
        MULTIQC_BUNDLING.out.fastqc_zips,
        MULTIQC_BUNDLING.out.inputs,
        MULTIQC_BUNDLING.out.config,
    )
    GENERATE_MULTIQC_REPORT.out.report.view {{ report -> "REPORT: ${{report.name}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "REPORT: multiqc_report.html" in completed.stdout
    manifests = sorted(tmp_path.glob("work/**/nvd_report_manifest.json"))
    assert manifests, diagnostics
    manifest = json.loads(manifests[-1].read_text(encoding="utf-8"))
    assert manifest["samples"] == [
        {
            "sample_id": "sample_A",
            "platform": "illumina",
            "source": "single_file",
            "read_structure": "single",
            "warnings": [],
        },
    ]
    assert manifest["raw_fastqc"] == []
    assert fastqc_attempt_log.read_text(encoding="utf-8").splitlines() == [
        "fastqc",
        "fastqc",
        "fastqc",
    ]


def test_single_read_fastqc_renders_through_real_multiqc(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    write_tool_fakes(bin_dir, fake_multiqc=False)
    read = write_fastq(tmp_path / "single.fastq")
    sample_id = "sample_A.fastq"
    roster = write_roster(
        tmp_path / "resolved_reads.jsonl",
        sample_id=sample_id,
        source="single_file",
        read_structure="single",
    )
    version = write_version(tmp_path / "nvd_version.txt")
    config = ROOT / "assets" / "multiqc_config.yaml"

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'
include {{ GENERATE_MULTIQC_REPORT }} from '{MULTIQC_MODULE}'
include {{ FASTQC_RAW }} from '{FASTQC_MODULE}'

params.results = '{tmp_path / "results"}'
params.experimental = false
workflow {{
    reads = Channel.of(tuple(
        [id: '{sample_id}', platform: 'illumina', source: 'single_file', read_mode: 'single', r1_count: 1],
        [{nextflow_file(read)}],
    ))
    fastqc_units = reads.flatMap {{ meta, read_files -> NvdReporting.processReadyFastqcTuples(meta, read_files) }}
    FASTQC_RAW(fastqc_units)
    MULTIQC_BUNDLING(
        FASTQC_RAW.out.packages,
        FASTQC_RAW.out.zips,
        Channel.value(file('{roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        {empty_domain_inputs(config)}
    )
    GENERATE_MULTIQC_REPORT(
        MULTIQC_BUNDLING.out.fastqc_zips,
        MULTIQC_BUNDLING.out.inputs,
        MULTIQC_BUNDLING.out.config,
    )
    GENERATE_MULTIQC_REPORT.out.report.view {{ report -> "REPORT: ${{report.name}}" }}
    GENERATE_MULTIQC_REPORT.out.data.view {{ data -> "DATA: ${{data.name}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "REPORT: multiqc_report.html" in completed.stdout
    assert "DATA: multiqc_data" in completed.stdout
    manifests = sorted(tmp_path.glob("work/**/nvd_report_manifest.json"))
    assert manifests, diagnostics
    manifest = json.loads(manifests[-1].read_text(encoding="utf-8"))
    assert manifest["raw_fastqc"] == [
        {
            "schema_version": "nvd.fastqc-raw-unit/v1",
            "sample_id": sample_id,
            "platform": "illumina",
            "source": "single_file",
            "read_end": "single",
            "input_ordinal": 1,
            "alias": "sample_A.fastq.raw.single.001",
            "zip_alias": "sample_A.fastq.raw.single.001_fastqc.zip",
            "html_alias": "sample_A.fastq.raw.single.001_fastqc.html",
        },
    ]
    data_dirs = sorted(tmp_path.glob("work/**/multiqc_data"))
    assert data_dirs, diagnostics
    assert "sample_A.fastq.raw.single.001" in read_native_fastqc_data(data_dirs[-1])
    fastqc_rows = (data_dirs[-1] / "multiqc_fastqc.txt").read_text(
        encoding="utf-8",
    ).splitlines()
    assert fastqc_rows[1].split("\t", 1)[0] == "sample_A.fastq.raw.single.001"
    roster_rows = (data_dirs[-1] / "multiqc_nvd_sample_roster.txt").read_text(
        encoding="utf-8",
    ).splitlines()
    assert "sample_id" not in roster_rows[0].split("\t")
    assert roster_rows[1].split("\t", 1)[0] == sample_id


def test_renderer_failure_is_ancillary_to_independent_sentinel(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    multiqc_attempt_log = tmp_path / "multiqc_attempts.log"
    write_tool_fakes(
        bin_dir,
        multiqc_fails=True,
        multiqc_attempt_log=multiqc_attempt_log,
    )
    roster = write_roster(tmp_path / "resolved_reads.jsonl")
    version = write_version(tmp_path / "nvd_version.txt")
    config = ROOT / "assets" / "multiqc_config.yaml"

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'
include {{ GENERATE_MULTIQC_REPORT }} from '{MULTIQC_MODULE}'

params.results = '{tmp_path / "results"}'
params.experimental = false
workflow {{
    sentinel = Channel.value('scientific-complete')
    MULTIQC_BUNDLING(
        Channel.empty(),
        Channel.empty(),
        Channel.value(file('{roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        {empty_domain_inputs(config)}
    )
    GENERATE_MULTIQC_REPORT(
        MULTIQC_BUNDLING.out.fastqc_zips,
        MULTIQC_BUNDLING.out.inputs,
        MULTIQC_BUNDLING.out.config,
    )
    GENERATE_MULTIQC_REPORT.out.report.view {{ report -> "REPORT: ${{report.name}}" }}
    GENERATE_MULTIQC_REPORT.out.data.view {{ data -> "DATA: ${{data.name}}" }}
    sentinel.view {{ value -> "SENTINEL: ${{value}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "SENTINEL: scientific-complete" in completed.stdout
    assert "REPORT: multiqc_report.html" not in completed.stdout
    assert "DATA: multiqc_data" not in completed.stdout
    assert not sorted(tmp_path.glob("work/**/multiqc_report.html"))
    assert not sorted(tmp_path.glob("work/**/multiqc_data"))
    assert multiqc_attempt_log.read_text(encoding="utf-8").splitlines() == [
        "multiqc",
        "multiqc",
        "multiqc",
    ]


def test_compiler_failure_remains_ancillary_to_independent_sentinel(
    tmp_path: Path,
) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    copy_reporting_lib(tmp_path)
    write_tool_fakes(bin_dir)
    compiler_attempt_log = tmp_path / "compiler_attempts.log"
    write_executable(
        bin_dir / "build_nvd_multiqc_inputs.py",
        f"#!/bin/sh\nprintf 'compiler\\n' >> {json.dumps(str(compiler_attempt_log))}\nexit 2\n",
    )
    bad_roster = tmp_path / "resolved_reads.jsonl"
    bad_roster.write_text('{"sample_id": "sample_A"}\n', encoding="utf-8")
    version = write_version(tmp_path / "nvd_version.txt")
    config = ROOT / "assets" / "multiqc_config.yaml"

    workflow = tmp_path / "main.nf"
    workflow.write_text(
        f"""\
nextflow.enable.dsl = 2

include {{ MULTIQC_BUNDLING }} from '{BUNDLING_SUBWORKFLOW}'

params.results = '{tmp_path / "results"}'
params.experimental = false
workflow {{
    sentinel = Channel.value('scientific-complete')
    MULTIQC_BUNDLING(
        Channel.empty(),
        Channel.empty(),
        Channel.value(file('{bad_roster}')),
        Channel.value(file('{version}')),
        Channel.value(false),
        {empty_domain_inputs(config)}
    )
    sentinel.view {{ value -> "SENTINEL: ${{value}}" }}
}}
""",
        encoding="utf-8",
    )

    completed = run_nextflow(workflow, bin_dir=bin_dir)
    diagnostics = f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.returncode == 0, diagnostics
    assert "SENTINEL: scientific-complete" in completed.stdout
    assert compiler_attempt_log.read_text(encoding="utf-8").splitlines() == [
        "compiler",
        "compiler",
        "compiler",
    ]
