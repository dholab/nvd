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
    create_mapped_bam_fifos,
    fastq_mapping_inputs_from_args,
    minimap2_command,
    minimap2_read_group_arg,
    remove_mapped_bam_fifos,
    sam_header_read_group_line,
    sample_bam_commands,
    threads_per_mapper,
    write_empty_unmapped_counts,
    write_empty_unmapped_outputs,
)

ONE_INPUT_MAPPER_THREADS = 2


def config(tmp_path: Path) -> MapbackConfig:
    return MapbackConfig(
        sample_id="sample1",
        platform="illumina",
        contigs=tmp_path / "contigs.fasta",
        threads=4,
        work_dir=tmp_path / "sample1.mapback_work",
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
        evidence_class="single_read",
        read_group_id="unmerged_reads",
        fastq=Path("reads.fastq.gz"),
    )

    assert minimap2_read_group_arg("sample1", input_) == (
        r"@RG\tID:unmerged_reads\tSM:sample1\tDS:evidence_class=single_read"
    )
    assert sam_header_read_group_line("sample1", input_) == (
        "@RG\tID:unmerged_reads\tSM:sample1\tDS:evidence_class=single_read\n"
    )


def test_minimap2_command_keeps_read_group_as_one_argument(tmp_path: Path) -> None:
    cfg = config(tmp_path)
    input_ = FastqMappingInput(
        evidence_class="overlap_merged_pair",
        read_group_id="overlap_merged_pair",
        fastq=Path("merged.fastq.gz"),
    )

    command = minimap2_command(cfg, input_, threads=2)

    read_group_arg = command[command.index("-R") + 1]
    assert read_group_arg == (
        r"@RG\tID:overlap_merged_pair\tSM:sample1\tDS:evidence_class=overlap_merged_pair"
    )
    assert "\t" not in read_group_arg


def test_fastq_input_args_reject_duplicate_evidence_classes(tmp_path: Path) -> None:
    first = write_fastq_gz(tmp_path / "first.fastq.gz", 1)
    second = write_fastq_gz(tmp_path / "second.fastq.gz", 1)

    with pytest.raises(MapbackError, match="evidence_class values must be unique"):
        fastq_mapping_inputs_from_args(
            [
                ["single_read", "unmerged_reads", str(first)],
                ["single_read", "other_unmerged_reads", str(second)],
            ],
        )


def test_mapped_bam_fifos_are_task_local_and_removed(tmp_path: Path) -> None:
    input_ = FastqMappingInput(
        evidence_class="single_read",
        read_group_id="unmerged_reads",
        fastq=Path("reads.fastq.gz"),
    )
    work_dir = tmp_path / "sample1.mapback_work"

    fifos = create_mapped_bam_fifos([input_], work_dir)
    try:
        assert fifos[0].path.parent == work_dir
        assert stat.S_ISFIFO(fifos[0].path.stat().st_mode)
    finally:
        remove_mapped_bam_fifos(fifos)

    assert not fifos[0].path.exists()


def test_sample_bam_commands_name_the_shared_dedup_pass(tmp_path: Path) -> None:
    cfg = config(tmp_path)
    input_ = FastqMappingInput(
        evidence_class="single_read",
        read_group_id="unmerged_reads",
        fastq=Path("reads.fastq.gz"),
    )
    fifos = create_mapped_bam_fifos([input_], cfg.work_dir)
    try:
        commands = sample_bam_commands(
            cfg,
            fifos,
            cfg.work_dir / "combined.header.sam",
            dedup_pos=True,
        )
    finally:
        remove_mapped_bam_fifos(fifos)

    flattened = [item for command in commands for item in command]
    assert "markdup" in flattened
    assert str(cfg.work_dir / "unmerged_reads.mapped.bam.fifo") in flattened


def test_empty_unmapped_outputs_keep_external_evidence_class_header(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)
    cfg = config(tmp_path)

    write_empty_unmapped_outputs(cfg)
    write_empty_unmapped_counts(cfg)

    for evidence_class in ("overlap_merged_pair", "single_read"):
        assert Path(
            f"sample1.{evidence_class}.mapback_unmapped_counts.tsv",
        ).read_text(encoding="utf-8") == (
            "sample_id\tplatform\tevidence_class\tunmapped_reads\n"
            f"sample1\tillumina\t{evidence_class}\t0\n"
        )
        with gzip.open(
            Path(f"sample1.{evidence_class}.mapback_unmapped.fastq.gz"),
            "rb",
        ) as handle:
            assert handle.read() == b""


def test_threads_per_mapper_leaves_capacity_for_the_shared_bam_writer() -> None:
    assert (
        threads_per_mapper(total_threads=4, input_count=1) == ONE_INPUT_MAPPER_THREADS
    )
    assert threads_per_mapper(total_threads=4, input_count=2) == 1
