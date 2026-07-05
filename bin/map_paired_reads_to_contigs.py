#!/usr/bin/env python3
"""Map paired-derived FASTQs to contigs through one shared BAM writer."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

EVIDENCE_CLASSES = ("overlap_merged_pair", "single_read")
SAFE_TOKEN = re.compile(r"^[A-Za-z0-9_.-]+$")


class MapbackError(RuntimeError):
    """Raised when paired read mapback cannot complete safely."""


@dataclass(frozen=True)
class MapbackConfig:
    """Nextflow task inputs needed to map one sample."""

    sample_id: str
    platform: str
    contigs: Path
    threads: int
    work_dir: Path
    minimap2_bin: str
    samtools_bin: str


@dataclass(frozen=True)
class FastqMappingInput:
    """One manifest row: a FASTQ plus its mapback provenance."""

    evidence_class: str
    read_group_id: str
    fastq: Path


@dataclass(frozen=True)
class MappedBamFifo:
    """Task-local FIFO carrying mapped BAM records for one FASTQ input."""

    input: FastqMappingInput
    path: Path


@dataclass(frozen=True)
class StartedCommand:
    """A named subprocess so failures can report the mapback stage that failed."""

    name: str
    command: list[str]
    process: subprocess.Popen[bytes]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map paired-derived FASTQs to contigs with one shared positional-dedup pass.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--contigs", required=True, type=Path)
    parser.add_argument(
        "--fastq-input",
        action="append",
        nargs=3,
        metavar=("EVIDENCE_CLASS", "READ_GROUP_ID", "FASTQ"),
        required=True,
        help="One FASTQ input for mapback. Repeat for merged and unmerged paired-derived reads.",
    )
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--dedup-pos", action="store_true")
    parser.add_argument("--minimap2-bin", default="minimap2")
    parser.add_argument("--samtools-bin", default="samtools")
    return parser.parse_args(argv)


def config_from_args(args: argparse.Namespace) -> MapbackConfig:
    sample_id = validate_safe_token("sample_id", args.sample_id)
    if args.threads < 1:
        message = "--threads must be at least 1"
        raise MapbackError(message)
    return MapbackConfig(
        sample_id=sample_id,
        platform=args.platform,
        contigs=args.contigs,
        threads=args.threads,
        work_dir=args.work_dir or Path(f"{sample_id}.mapback_work"),
        minimap2_bin=args.minimap2_bin,
        samtools_bin=args.samtools_bin,
    )


def validate_safe_token(name: str, value: str) -> str:
    if not value or not SAFE_TOKEN.fullmatch(value):
        message = f"{name} must contain only letters, numbers, dots, underscores, or hyphens: {value!r}"
        raise MapbackError(message)
    return value


def fastq_mapping_inputs_from_args(
    fastq_inputs: Sequence[Sequence[str]],
) -> tuple[FastqMappingInput, ...]:
    inputs = tuple(
        FastqMappingInput(
            evidence_class=validate_safe_token("evidence_class", evidence_class),
            read_group_id=validate_safe_token("read_group_id", read_group_id),
            fastq=Path(fastq),
        )
        for evidence_class, read_group_id, fastq in fastq_inputs
    )
    validate_fastq_mapping_inputs(inputs)
    return inputs


def validate_fastq_mapping_inputs(inputs: Sequence[FastqMappingInput]) -> None:
    if not inputs:
        message = "provide at least one FASTQ mapping input"
        raise MapbackError(message)
    duplicate_read_groups = duplicates(input_.read_group_id for input_ in inputs)
    if duplicate_read_groups:
        message = "read_group_id values must be unique: " + ", ".join(
            duplicate_read_groups,
        )
        raise MapbackError(message)
    duplicate_classes = duplicates(input_.evidence_class for input_ in inputs)
    if duplicate_classes:
        message = (
            "evidence_class values must be unique because unmapped FASTQs are class-split: "
            + ", ".join(duplicate_classes)
        )
        raise MapbackError(message)
    unknown_classes = sorted(
        {input_.evidence_class for input_ in inputs} - set(EVIDENCE_CLASSES),
    )
    if unknown_classes:
        message = "unsupported evidence_class values: " + ", ".join(unknown_classes)
        raise MapbackError(message)
    missing_fastqs = [input_.fastq for input_ in inputs if not input_.fastq.is_file()]
    if missing_fastqs:
        message = "FASTQ mapping inputs contain missing files:\n" + "\n".join(
            f"  {path}" for path in missing_fastqs
        )
        raise MapbackError(message)


def duplicates(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    repeated: set[str] = set()
    for value in values:
        if value in seen:
            repeated.add(value)
        seen.add(value)
    return sorted(repeated)


def run_mapback(
    config: MapbackConfig,
    inputs: Sequence[FastqMappingInput],
    *,
    dedup_pos: bool,
) -> None:
    fifos = create_mapped_bam_fifos(inputs, config.work_dir)
    commands: list[StartedCommand] = []
    completed = False

    try:
        write_empty_unmapped_outputs(config)
        header = config.work_dir / "combined.header.sam"
        write_combined_sam_header(config, inputs, header)

        commands.extend(
            start_pipeline(
                "sample BAM writer",
                sample_bam_commands(config, fifos, header, dedup_pos=dedup_pos),
            ),
        )
        mapper_threads = threads_per_mapper(config.threads, len(inputs))
        for fifo in fifos:
            commands.extend(start_fastq_mapper(config, fifo, mapper_threads))

        wait_all_or_stop(commands)
        index_bam(config)
        write_empty_unmapped_counts(config)
        completed = True
    finally:
        terminate_running(commands)
        remove_mapped_bam_fifos(fifos)
        if completed:
            shutil.rmtree(config.work_dir, ignore_errors=True)


def create_mapped_bam_fifos(
    inputs: Sequence[FastqMappingInput],
    work_dir: Path,
) -> tuple[MappedBamFifo, ...]:
    work_dir.mkdir(parents=True, exist_ok=True)
    fifos = tuple(
        MappedBamFifo(
            input=input_,
            path=work_dir / f"{input_.read_group_id}.mapped.bam.fifo",
        )
        for input_ in inputs
    )
    for fifo in fifos:
        if fifo.path.exists():
            fifo.path.unlink()
        os.mkfifo(fifo.path)
    return fifos


def remove_mapped_bam_fifos(fifos: Sequence[MappedBamFifo]) -> None:
    for fifo in fifos:
        fifo.path.unlink(missing_ok=True)


def write_empty_unmapped_outputs(config: MapbackConfig) -> None:
    for evidence_class in EVIDENCE_CLASSES:
        with gzip.open(
            unmapped_fastq_path(config.sample_id, evidence_class),
            "wb",
        ) as handle:
            handle.write(b"")


def write_empty_unmapped_counts(config: MapbackConfig) -> None:
    header = "sample_id\tplatform\tevidence_class\tunmapped_reads\n"
    for evidence_class in EVIDENCE_CLASSES:
        row = f"{config.sample_id}\t{config.platform}\t{evidence_class}\t0\n"
        unmapped_counts_path(config.sample_id, evidence_class).write_text(
            header + row,
            encoding="utf-8",
        )


def write_combined_sam_header(
    config: MapbackConfig,
    inputs: Sequence[FastqMappingInput],
    header: Path,
) -> None:
    fai_path = Path(str(config.contigs) + ".fai")
    run_checked([config.samtools_bin, "faidx", str(config.contigs)])
    lines = ["@HD\tVN:1.6\tSO:unsorted\n"]
    for row in fai_path.read_text(encoding="utf-8").splitlines():
        contig_name, contig_length, *_rest = row.split("\t")
        lines.append(f"@SQ\tSN:{contig_name}\tLN:{contig_length}\n")
    lines.extend(
        sam_header_read_group_line(config.sample_id, input_) for input_ in inputs
    )
    header.write_text("".join(lines), encoding="utf-8")


def sample_bam_commands(
    config: MapbackConfig,
    fifos: Sequence[MappedBamFifo],
    header: Path,
    *,
    dedup_pos: bool,
) -> list[list[str]]:
    cat = [
        config.samtools_bin,
        "cat",
        "-h",
        str(header),
        *(str(fifo.path) for fifo in fifos),
    ]
    if dedup_pos:
        return [
            cat,
            [
                config.samtools_bin,
                "collate",
                "-@",
                str(config.threads),
                "-O",
                "-u",
                "-",
            ],
            [
                config.samtools_bin,
                "fixmate",
                "-@",
                str(config.threads),
                "-m",
                "-u",
                "-",
                "-",
            ],
            [config.samtools_bin, "sort", "-@", str(config.threads), "-u", "-"],
            [
                config.samtools_bin,
                "markdup",
                "-@",
                str(config.threads),
                "-s",
                "-r",
                "--duplicate-count",
                "-",
                str(sample_bam_path(config)),
            ],
        ]
    return [
        cat,
        [config.samtools_bin, "view", "-b", "-F", "4", "-"],
        [
            config.samtools_bin,
            "sort",
            "-@",
            str(config.threads),
            "-o",
            str(sample_bam_path(config)),
        ],
    ]


def start_fastq_mapper(
    config: MapbackConfig,
    fifo: MappedBamFifo,
    threads: int,
) -> list[StartedCommand]:
    minimap2 = subprocess.Popen(  # noqa: S603
        minimap2_command(config, fifo.input, threads),
        stdout=subprocess.PIPE,
    )
    mapped = subprocess.Popen(  # noqa: S603
        [
            config.samtools_bin,
            "view",
            "-u",
            "-F",
            "4",
            "-o",
            str(fifo.path),
            "-",
        ],
        stdin=minimap2.stdout,
    )
    if minimap2.stdout is not None:
        minimap2.stdout.close()
    return [
        StartedCommand(
            name=f"minimap2 {fifo.input.read_group_id}",
            command=minimap2_command(config, fifo.input, threads),
            process=minimap2,
        ),
        StartedCommand(
            name=f"samtools mapped BAM {fifo.input.read_group_id}",
            command=[
                config.samtools_bin,
                "view",
                "-u",
                "-F",
                "4",
                "-o",
                str(fifo.path),
                "-",
            ],
            process=mapped,
        ),
    ]


def start_pipeline(name: str, commands: list[list[str]]) -> list[StartedCommand]:
    started: list[StartedCommand] = []
    previous_stdout = None
    for index, command in enumerate(commands):
        is_last = index == len(commands) - 1
        process = subprocess.Popen(  # noqa: S603
            command,
            stdin=previous_stdout,
            stdout=None if is_last else subprocess.PIPE,
        )
        if previous_stdout is not None:
            previous_stdout.close()
        previous_stdout = process.stdout
        started.append(
            StartedCommand(
                name=f"{name}: {command[1] if len(command) > 1 else command[0]}",
                command=command,
                process=process,
            ),
        )
    return started


def wait_all_or_stop(commands: Sequence[StartedCommand]) -> None:
    failures: list[str] = []
    for started in commands:
        return_code = started.process.wait()
        if return_code != 0:
            failures.append(
                f"{started.name} exited {return_code}: {format_command(started.command)}",
            )
            terminate_running(commands)
    if failures:
        raise MapbackError("mapback command failed:\n" + "\n".join(failures))


def terminate_running(commands: Sequence[StartedCommand]) -> None:
    for started in commands:
        if started.process.poll() is None:
            started.process.terminate()


def minimap2_command(
    config: MapbackConfig,
    input_: FastqMappingInput,
    threads: int,
) -> list[str]:
    return [
        config.minimap2_bin,
        "-ax",
        minimap2_preset(config.platform),
        "-R",
        minimap2_read_group_arg(config.sample_id, input_),
        "-t",
        str(threads),
        str(config.contigs),
        str(input_.fastq),
    ]


def minimap2_preset(platform: str) -> str:
    return "map-ont" if platform == "ont" else "sr"


def minimap2_read_group_arg(sample_id: str, input_: FastqMappingInput) -> str:
    return (
        r"@RG\tID:"
        + input_.read_group_id
        + r"\tSM:"
        + sample_id
        + r"\tDS:evidence_class="
        + input_.evidence_class
    )


def sam_header_read_group_line(sample_id: str, input_: FastqMappingInput) -> str:
    return (
        f"@RG\tID:{input_.read_group_id}\tSM:{sample_id}"
        f"\tDS:evidence_class={input_.evidence_class}\n"
    )


def threads_per_mapper(total_threads: int, input_count: int) -> int:
    return max(1, total_threads // (input_count + 1))


def sample_bam_path(config: MapbackConfig) -> Path:
    return Path(f"{config.sample_id}.bam")


def unmapped_fastq_path(sample_id: str, evidence_class: str) -> Path:
    return Path(f"{sample_id}.{evidence_class}.mapback_unmapped.fastq.gz")


def unmapped_counts_path(sample_id: str, evidence_class: str) -> Path:
    return Path(f"{sample_id}.{evidence_class}.mapback_unmapped_counts.tsv")


def index_bam(config: MapbackConfig) -> None:
    run_checked([config.samtools_bin, "index", str(sample_bam_path(config))])


def run_checked(command: list[str]) -> None:
    completed = subprocess.run(command, check=False)  # noqa: S603
    if completed.returncode != 0:
        message = f"command exited {completed.returncode}: {format_command(command)}"
        raise MapbackError(message)


def format_command(command: Sequence[str]) -> str:
    return " ".join(command)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    config = config_from_args(args)
    inputs = fastq_mapping_inputs_from_args(args.fastq_input)
    run_mapback(
        config,
        inputs,
        dedup_pos=bool(args.dedup_pos),
    )


if __name__ == "__main__":
    main()
