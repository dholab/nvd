"""Tests for the paired read mapback helper."""

from __future__ import annotations

import gzip
import stat
from pathlib import Path

import pytest
from map_paired_reads_to_contigs import (
    FastqMappingInput,
    MapbackConfig,
    MapbackError,
    MappedBamFifo,
    UnmappedBamFifo,
    create_mapped_bam_fifos,
    fastq_mapping_inputs_from_args,
    mapped_bam_fifo_path,
    minimap2_command,
    minimap2_read_group_arg,
    remove_fifos,
    sam_header_read_group_line,
    sample_bam_commands,
    start_fastq_mapper,
    threads_per_mapper,
    write_empty_unmapped_outputs,
    write_unmapped_counts,
)

ONE_INPUT_MAPPER_THREADS = 2


def config(tmp_path: Path) -> MapbackConfig:
    return MapbackConfig(
        sample_id="sample1",
        platform="illumina",
        contigs=tmp_path / "contigs.fasta",
        threads=4,
        minimap2_bin="minimap2",
        samtools_bin="samtools",
    )


def write_fastq_gz(path: Path, record_count: int) -> Path:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for index in range(record_count):
            handle.write(f"@read{index}\nACGT\n+\n!!!!\n")
    return path


def test_read_group_strings_use_minimap2_escaped_tabs_and_header_real_tabs() -> None:
    input_ = FastqMappingInput(
        query_class="single_read",
        read_group_id="unmerged_reads",
        fastq=Path("reads.fastq.gz"),
    )

    assert minimap2_read_group_arg("sample1", input_) == (
        r"@RG\tID:unmerged_reads\tSM:sample1\tDS:query_class=single_read"
    )
    assert sam_header_read_group_line("sample1", input_) == (
        "@RG\tID:unmerged_reads\tSM:sample1\tDS:query_class=single_read\n"
    )


def test_minimap2_command_keeps_read_group_as_one_argument(tmp_path: Path) -> None:
    cfg = config(tmp_path)
    input_ = FastqMappingInput(
        query_class="overlap_merged_pair",
        read_group_id="overlap_merged_pair",
        fastq=Path("merged.fastq.gz"),
    )

    command = minimap2_command(cfg, input_, threads=2)

    read_group_arg = command[command.index("-R") + 1]
    assert read_group_arg == (
        r"@RG\tID:overlap_merged_pair\tSM:sample1\tDS:query_class=overlap_merged_pair"
    )
    assert "\t" not in read_group_arg


def test_fastq_input_args_reject_duplicate_query_classes(tmp_path: Path) -> None:
    first = write_fastq_gz(tmp_path / "first.fastq.gz", 1)
    second = write_fastq_gz(tmp_path / "second.fastq.gz", 1)

    with pytest.raises(MapbackError, match="query_class values must be unique"):
        fastq_mapping_inputs_from_args(
            [
                ["single_read", "unmerged_reads", str(first)],
                ["single_read", "other_unmerged_reads", str(second)],
            ],
        )


def test_mapped_bam_fifos_are_task_local_hidden_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    input_ = FastqMappingInput(
        query_class="single_read",
        read_group_id="unmerged_reads",
        fastq=Path("reads.fastq.gz"),
    )

    fifos = create_mapped_bam_fifos([input_])
    try:
        assert fifos[0].path == Path(".mapback.unmerged_reads.mapped.bam.fifo")
        assert stat.S_ISFIFO(fifos[0].path.stat().st_mode)
    finally:
        remove_fifos(fifos)

    assert not fifos[0].path.exists()


def test_sample_bam_commands_name_the_shared_dedup_pass(tmp_path: Path) -> None:
    cfg = config(tmp_path)
    input_ = FastqMappingInput(
        query_class="single_read",
        read_group_id="unmerged_reads",
        fastq=Path("reads.fastq.gz"),
    )
    fifos = [MappedBamFifo(input=input_, path=mapped_bam_fifo_path("unmerged_reads"))]
    commands = sample_bam_commands(
        cfg,
        fifos,
        Path(".mapback.combined.header.sam"),
        dedup_pos=True,
    )

    flattened = [item for command in commands for item in command]
    assert "markdup" in flattened
    assert ".mapback.unmerged_reads.mapped.bam.fifo" in flattened


def test_mapper_command_can_split_unmapped_bam_with_samtools_u(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = config(tmp_path)
    input_ = FastqMappingInput(
        query_class="single_read",
        read_group_id="single_reads",
        fastq=Path("reads.fastq.gz"),
    )
    mapped_fifo = MappedBamFifo(
        input=input_,
        path=Path(".mapback.single_reads.mapped.bam.fifo"),
    )
    unmapped_fifo = UnmappedBamFifo(
        input=input_,
        path=Path(".mapback.single_reads.unmapped.bam.fifo"),
    )
    started_commands = []

    class FakePopen:
        def __init__(self, command: list[str], **_kwargs: object) -> None:
            self.command = command
            self.stdout = None
            started_commands.append(command)

    monkeypatch.setattr("map_paired_reads_to_contigs.subprocess.Popen", FakePopen)

    commands = start_fastq_mapper(cfg, mapped_fifo, unmapped_fifo, threads=2)

    view_command = commands[1].command
    assert view_command == [
        "samtools",
        "view",
        "-u",
        "-F",
        "4",
        "-o",
        ".mapback.single_reads.mapped.bam.fifo",
        "-U",
        ".mapback.single_reads.unmapped.bam.fifo",
        "-",
    ]
    assert started_commands[1] == view_command


def test_unmapped_counts_are_written_per_external_query_class(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    cfg = config(tmp_path)

    write_empty_unmapped_outputs(cfg)
    expected_counts = {"overlap_merged_pair": 2, "single_read": 1}
    for query_class, record_count in expected_counts.items():
        write_fastq_gz(
            Path(f"sample1.{query_class}.mapback_unmapped.fastq.gz"),
            record_count,
        )
    inputs = [
        FastqMappingInput(
            query_class="overlap_merged_pair",
            read_group_id="overlap_merged_pairs",
            fastq=Path("merged.fastq.gz"),
        ),
        FastqMappingInput(
            query_class="single_read",
            read_group_id="single_reads",
            fastq=Path("single.fastq.gz"),
        ),
    ]
    write_unmapped_counts(cfg, inputs)

    assert not Path("sample1.mapback_unmapped_counts.tsv").exists()
    for query_class, expected_count in expected_counts.items():
        assert Path(
            f"sample1.{query_class}.mapback_unmapped_counts.tsv",
        ).read_text(encoding="utf-8") == (
            "sample_id\tplatform\tquery_class\tunmapped_reads\n"
            f"sample1\tillumina\t{query_class}\t{expected_count}\n"
        )


def test_unmapped_counts_include_absent_query_classes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    cfg = config(tmp_path)
    write_empty_unmapped_outputs(cfg)
    write_fastq_gz(Path("sample1.single_read.mapback_unmapped.fastq.gz"), 1)

    write_unmapped_counts(
        cfg,
        [
            FastqMappingInput(
                query_class="single_read",
                read_group_id="single_reads",
                fastq=Path("single.fastq.gz"),
            ),
        ],
    )

    assert Path(
        "sample1.overlap_merged_pair.mapback_unmapped_counts.tsv",
    ).read_text(encoding="utf-8") == (
        "sample_id\tplatform\tquery_class\tunmapped_reads\n"
        "sample1\tillumina\toverlap_merged_pair\t0\n"
    )
    assert Path("sample1.single_read.mapback_unmapped_counts.tsv").exists()


def test_threads_per_mapper_leaves_capacity_for_the_shared_bam_writer() -> None:
    assert (
        threads_per_mapper(total_threads=4, input_count=1) == ONE_INPUT_MAPPER_THREADS
    )
    assert threads_per_mapper(total_threads=4, input_count=2) == 1
