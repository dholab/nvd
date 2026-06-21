#!/usr/bin/env python3
"""Stream ordered FASTQ bundles into Deacon through named pipes."""

from __future__ import annotations

import argparse
import errno
import lzma
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from threading import Event
from typing import BinaryIO

from py_nvd.read_filenames import Compression, compression_kind

DeaconRunner = Callable[[list[str]], subprocess.CompletedProcess[str]]


class StreamError(RuntimeError):
    """Raised when FASTQ streaming or Deacon execution fails."""


@dataclass(frozen=True)
class DeaconStreamConfig:
    sample_id: str
    index: Path
    output: Path
    summary: Path
    reads: tuple[Path, ...]
    r1: tuple[Path, ...]
    r2: tuple[Path, ...]
    threads: int
    abs_threshold: int
    rel_threshold: float
    deacon_bin: str


@dataclass(frozen=True)
class ProducerResult:
    name: str
    files: tuple[Path, ...]
    bytes_written: int


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run deacon filter on one or more ordered FASTQ inputs without concatenating them on disk.",
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--index", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--reads-list", type=Path)
    parser.add_argument("--r1-list", type=Path)
    parser.add_argument("--r2-list", type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--abs-threshold", type=int, default=1)
    parser.add_argument("--rel-threshold", type=float, default=0.0)
    parser.add_argument("--deacon-bin", default="deacon")
    return parser.parse_args(argv)


def read_path_list(path: Path | None) -> tuple[Path, ...]:
    if path is None:
        return ()
    return tuple(
        Path(line.strip())
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )


def config_from_args(args: argparse.Namespace) -> DeaconStreamConfig:
    return DeaconStreamConfig(
        sample_id=args.sample_id,
        index=args.index,
        output=args.output,
        summary=args.summary,
        reads=read_path_list(args.reads_list),
        r1=read_path_list(args.r1_list),
        r2=read_path_list(args.r2_list),
        threads=args.threads,
        abs_threshold=args.abs_threshold,
        rel_threshold=args.rel_threshold,
        deacon_bin=args.deacon_bin,
    )


def one_compression(files: tuple[Path, ...]) -> Compression:
    compressions = {compression_kind(path) for path in files}
    if len(compressions) == 1:
        return next(iter(compressions))
    raise StreamError(
        "FASTQ bundles streamed to Deacon must use one compression type; found "
        + ", ".join(sorted(compressions)),
    )


def validate_files(name: str, files: tuple[Path, ...]) -> None:
    if not files:
        message = f"{name} file list is empty"
        raise StreamError(message)
    missing = [path for path in files if not path.is_file()]
    if missing:
        raise StreamError(
            f"{name} file list contains missing files:\n"
            + "\n".join(f"  {path}" for path in missing),
        )
    one_compression(files)


def validate_sample_compression(files: tuple[Path, ...]) -> None:
    one_compression(files)


def validate_config(
    config: DeaconStreamConfig,
    *,
    require_deacon_executable: bool = True,
) -> None:
    if require_deacon_executable:
        deacon_path = shutil.which(config.deacon_bin)
        explicit_path = Path(config.deacon_bin)
        explicit_executable = explicit_path.is_file() and os.access(
            explicit_path,
            os.X_OK,
        )
        if deacon_path is None and not explicit_executable:
            message = f"deacon executable was not found: {config.deacon_bin}"
            raise StreamError(message)
    single = bool(config.reads)
    paired = bool(config.r1 or config.r2)
    if single == paired:
        message = "Provide either --reads-list for single-end input or both --r1-list and --r2-list for paired input"
        raise StreamError(message)
    if single:
        validate_files("reads", config.reads)
        validate_sample_compression(config.reads)
        return
    validate_files("r1", config.r1)
    validate_files("r2", config.r2)
    validate_sample_compression(config.r1 + config.r2)
    if len(config.r1) != len(config.r2):
        message = (
            f"R1/R2 file list length mismatch: {len(config.r1)} != {len(config.r2)}"
        )
        raise StreamError(message)


def fifo_suffix(compression: Compression) -> str:
    match compression:
        case "gz" | "zst":
            return f".fastq.{compression}"
        case "none" | "xz":
            return ".fastq"


def open_deacon_payload(path: Path) -> BinaryIO:
    match compression_kind(path):
        case "none" | "gz" | "zst":
            return path.open("rb")
        case "xz":
            return lzma.open(path, "rb")


def open_fifo_for_write(fifo: Path, stop: Event) -> BinaryIO:
    while not stop.is_set():
        try:
            fd = os.open(fifo, os.O_WRONLY | os.O_NONBLOCK)
            os.set_blocking(fd, True)
            return os.fdopen(fd, "wb")
        except OSError as error:
            if error.errno != errno.ENXIO:
                raise
            time.sleep(0.01)
    message = f"no reader opened FIFO before Deacon exited: {fifo}"
    raise StreamError(message)


def stream_files_to_fifo(
    name: str,
    files: tuple[Path, ...],
    fifo: Path,
    stop: Event,
) -> ProducerResult:
    bytes_written = 0
    with open_fifo_for_write(fifo, stop) as output:
        for path in files:
            with open_deacon_payload(path) as source:
                shutil.copyfileobj(source, output)
                bytes_written += path.stat().st_size
    return ProducerResult(name=name, files=files, bytes_written=bytes_written)


def deacon_command(config: DeaconStreamConfig, inputs: tuple[Path, ...]) -> list[str]:
    return [
        config.deacon_bin,
        "filter",
        "--threads",
        str(config.threads),
        "--abs-threshold",
        str(config.abs_threshold),
        "--rel-threshold",
        str(config.rel_threshold),
        "--summary",
        str(config.summary),
        "--output",
        str(config.output),
        str(config.index),
        *(str(path) for path in inputs),
    ]


def subprocess_deacon_runner(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=False, text=True)  # noqa: S603


def await_producers(
    futures: tuple[Future[ProducerResult], ...],
) -> tuple[ProducerResult, ...]:
    results: list[ProducerResult] = []
    errors: list[str] = []
    for future in futures:
        try:
            results.append(future.result())
        except Exception as error:  # noqa: BLE001
            errors.append(str(error))
    if errors:
        raise StreamError("FASTQ stream producer failed:\n" + "\n".join(errors))
    return tuple(results)


def run_deacon_stream(
    config: DeaconStreamConfig,
    runner: DeaconRunner | None = None,
) -> None:
    if runner is None:
        runner = subprocess_deacon_runner
    validate_config(
        config,
        require_deacon_executable=runner is subprocess_deacon_runner,
    )
    config.output.parent.mkdir(parents=True, exist_ok=True)
    config.summary.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"{config.sample_id}.deacon-input.") as tmp:
        tmpdir = Path(tmp)
        if config.reads:
            suffix = fifo_suffix(compression_kind(config.reads[0]))
            fifo_specs = (("reads", config.reads, tmpdir / f"reads{suffix}"),)
        else:
            suffix = fifo_suffix(compression_kind(config.r1[0]))
            fifo_specs = (
                ("r1", config.r1, tmpdir / f"r1{suffix}"),
                ("r2", config.r2, tmpdir / f"r2{suffix}"),
            )

        command_inputs = tuple(fifo for _name, _files, fifo in fifo_specs)

        for _name, _files, fifo in fifo_specs:
            os.mkfifo(fifo)

        stop = Event()
        with ThreadPoolExecutor(max_workers=len(fifo_specs)) as executor:
            futures = tuple(
                executor.submit(stream_files_to_fifo, name, files, fifo, stop)
                for name, files, fifo in fifo_specs
            )
            command = deacon_command(config, command_inputs)
            completed = runner(command)
            stop.set()
            try:
                await_producers(futures)
            except StreamError:
                if completed.returncode != 0:
                    message = (
                        f"deacon filter failed with exit code {completed.returncode}"
                    )
                    raise StreamError(message) from None
                raise
            if completed.returncode != 0:
                message = f"deacon filter failed with exit code {completed.returncode}"
                raise StreamError(message)


def main() -> None:
    try:
        run_deacon_stream(config_from_args(parse_args()))
    except StreamError as error:
        sys.stderr.write(f"{error}\n")
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
