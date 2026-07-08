#!/usr/bin/env python3
"""Map paired-derived FASTQs to contigs through one shared BAM writer."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

QUERY_CLASSES = ("overlap_merged_pair", "single_read")
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
    minimap2_bin: str
    samtools_bin: str


@dataclass(frozen=True)
class FastqMappingInput:
    """One manifest row: a FASTQ plus its mapback provenance."""

    query_class: str
    read_group_id: str
    fastq: Path


@dataclass(frozen=True)
class MappedBamFifo:
    """Task-local FIFO carrying mapped BAM records for one FASTQ input."""

    input: FastqMappingInput
    path: Path


@dataclass(frozen=True)
class UnmappedBamFifo:
    """Task-local FIFO carrying unmapped BAM records for one FASTQ input."""

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
        metavar=("QUERY_CLASS", "READ_GROUP_ID", "FASTQ"),
        required=True,
        help="One FASTQ input for mapback. Repeat for merged and unmerged paired-derived reads.",
    )
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--dedup-pos", action="store_true")
    parser.add_argument("--extract-unmapped", action="store_true")
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
            query_class=validate_safe_token("query_class", query_class),
            read_group_id=validate_safe_token("read_group_id", read_group_id),
            fastq=Path(fastq),
        )
        for query_class, read_group_id, fastq in fastq_inputs
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
    duplicate_classes = duplicates(input_.query_class for input_ in inputs)
    if duplicate_classes:
        message = (
            "query_class values must be unique because unmapped FASTQs are class-split: "
            + ", ".join(duplicate_classes)
        )
        raise MapbackError(message)
    unknown_classes = sorted(
        {input_.query_class for input_ in inputs} - set(QUERY_CLASSES),
    )
    if unknown_classes:
        message = "unsupported query_class values: " + ", ".join(unknown_classes)
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
    extract_unmapped: bool,
) -> None:
    mapped_fifos = create_mapped_bam_fifos(inputs)
    unmapped_fifos = create_unmapped_bam_fifos(inputs) if extract_unmapped else ()
    commands: list[StartedCommand] = []

    try:
        write_empty_unmapped_outputs(config)
        header = combined_header_path()
        write_combined_sam_header(config, inputs, header)

        commands.extend(
            start_pipeline(
                "sample BAM writer",
                sample_bam_commands(config, mapped_fifos, header, dedup_pos=dedup_pos),
            ),
        )

        if extract_unmapped:
            for fifo in unmapped_fifos:
                commands.extend(start_unmapped_fastq_writer(config, fifo))

        mapper_threads = threads_per_mapper(config.threads, len(inputs))
        for mapped_fifo in mapped_fifos:
            commands.extend(
                start_fastq_mapper(
                    config,
                    mapped_fifo,
                    matching_unmapped_fifo(mapped_fifo, unmapped_fifos),
                    mapper_threads,
                ),
            )

        wait_all_or_stop(commands)
        index_bam(config)
        write_unmapped_counts(config, inputs)
    finally:
        terminate_running(commands)
        remove_fifos(mapped_fifos)
        remove_fifos(unmapped_fifos)
        combined_header_path().unlink(missing_ok=True)


def create_mapped_bam_fifos(
    inputs: Sequence[FastqMappingInput],
) -> tuple[MappedBamFifo, ...]:
    fifos = tuple(
        MappedBamFifo(
            input=input_,
            path=mapped_bam_fifo_path(input_.read_group_id),
        )
        for input_ in inputs
    )
    create_fifos(fifos)
    return fifos


def create_unmapped_bam_fifos(
    inputs: Sequence[FastqMappingInput],
) -> tuple[UnmappedBamFifo, ...]:
    fifos = tuple(
        UnmappedBamFifo(
            input=input_,
            path=unmapped_bam_fifo_path(input_.read_group_id),
        )
        for input_ in inputs
    )
    create_fifos(fifos)
    return fifos


def create_fifos(fifos: Sequence[MappedBamFifo | UnmappedBamFifo]) -> None:
    for fifo in fifos:
        if fifo.path.exists():
            fifo.path.unlink()
        os.mkfifo(fifo.path)


def remove_fifos(fifos: Sequence[MappedBamFifo | UnmappedBamFifo]) -> None:
    for fifo in fifos:
        fifo.path.unlink(missing_ok=True)


def write_empty_unmapped_outputs(config: MapbackConfig) -> None:
    for query_class in QUERY_CLASSES:
        with gzip.open(
            unmapped_fastq_path(config.sample_id, query_class),
            "wb",
        ) as handle:
            handle.write(b"")


def write_unmapped_counts(
    config: MapbackConfig,
    inputs: Sequence[FastqMappingInput],
) -> None:
    header = "sample_id\tplatform\tquery_class\tunmapped_reads\n"
    for input_ in inputs:
        unmapped_reads = count_fastq_records(
            unmapped_fastq_path(config.sample_id, input_.query_class),
        )
        row = (
            f"{config.sample_id}\t{config.platform}\t{input_.query_class}"
            f"\t{unmapped_reads}\n"
        )
        unmapped_counts_path(config.sample_id, input_.query_class).write_text(
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
    unmapped_fifo: UnmappedBamFifo | None,
    threads: int,
) -> list[StartedCommand]:
    view_command = [
        config.samtools_bin,
        "view",
        "-u",
        "-F",
        "4",
        "-o",
        str(fifo.path),
    ]
    if unmapped_fifo is not None:
        view_command += ["-U", str(unmapped_fifo.path)]
    view_command += ["-"]

    minimap2 = subprocess.Popen(  # noqa: S603
        minimap2_command(config, fifo.input, threads),
        stdout=subprocess.PIPE,
    )
    mapped = subprocess.Popen(  # noqa: S603
        view_command,
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
            command=view_command,
            process=mapped,
        ),
    ]


def start_unmapped_fastq_writer(
    config: MapbackConfig,
    fifo: UnmappedBamFifo,
) -> list[StartedCommand]:
    fastq_command = [config.samtools_bin, "fastq", str(fifo.path)]
    gzip_command = ["gzip", "-c"]
    gzip_output = unmapped_fastq_path(config.sample_id, fifo.input.query_class)

    fastq = subprocess.Popen(  # noqa: S603
        fastq_command,
        stdout=subprocess.PIPE,
    )
    gzip_handle = gzip_output.open("wb")
    gzip_process = subprocess.Popen(  # noqa: S603
        gzip_command,
        stdin=fastq.stdout,
        stdout=gzip_handle,
    )
    gzip_handle.close()
    if fastq.stdout is not None:
        fastq.stdout.close()
    return [
        StartedCommand(
            name=f"samtools unmapped FASTQ {fifo.input.read_group_id}",
            command=fastq_command,
            process=fastq,
        ),
        StartedCommand(
            name=f"gzip unmapped FASTQ {fifo.input.read_group_id}",
            command=gzip_command,
            process=gzip_process,
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
        + r"\tDS:query_class="
        + input_.query_class
    )


def sam_header_read_group_line(sample_id: str, input_: FastqMappingInput) -> str:
    return (
        f"@RG\tID:{input_.read_group_id}\tSM:{sample_id}"
        f"\tDS:query_class={input_.query_class}\n"
    )


def threads_per_mapper(total_threads: int, input_count: int) -> int:
    return max(1, total_threads // (input_count + 1))


def matching_unmapped_fifo(
    mapped_fifo: MappedBamFifo,
    unmapped_fifos: Sequence[UnmappedBamFifo],
) -> UnmappedBamFifo | None:
    for unmapped_fifo in unmapped_fifos:
        if unmapped_fifo.input.read_group_id == mapped_fifo.input.read_group_id:
            return unmapped_fifo
    return None


def combined_header_path() -> Path:
    return Path(".mapback.combined.header.sam")


def mapped_bam_fifo_path(read_group_id: str) -> Path:
    return Path(f".mapback.{read_group_id}.mapped.bam.fifo")


def unmapped_bam_fifo_path(read_group_id: str) -> Path:
    return Path(f".mapback.{read_group_id}.unmapped.bam.fifo")


def sample_bam_path(config: MapbackConfig) -> Path:
    return Path(f"{config.sample_id}.bam")


def unmapped_fastq_path(sample_id: str, query_class: str) -> Path:
    return Path(f"{sample_id}.{query_class}.mapback_unmapped.fastq.gz")


def unmapped_counts_path(sample_id: str, query_class: str) -> Path:
    return Path(f"{sample_id}.{query_class}.mapback_unmapped_counts.tsv")


def index_bam(config: MapbackConfig) -> None:
    run_checked([config.samtools_bin, "index", str(sample_bam_path(config))])


def count_fastq_records(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for _line in handle) // 4


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
        extract_unmapped=bool(args.extract_unmapped),
    )


if __name__ == "__main__":
    main()
