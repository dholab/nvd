"""Stack per-sample Big Table TSVs after validating their schemas."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from collections.abc import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stack per-sample Big Table TSVs into one all-sample TSV.",
    )
    parser.add_argument("--input", type=Path, nargs="*", default=[])
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def scan_tsv(path: Path) -> pl.LazyFrame:
    return pl.scan_csv(path, separator="\t", infer_schema=False)


def schema_columns(path: Path) -> list[str]:
    return scan_tsv(path).collect_schema().names()


def schema_mismatch_message(
    path: Path,
    expected: Sequence[str],
    observed: Sequence[str],
) -> str:
    missing = [column for column in expected if column not in observed]
    extra = [column for column in observed if column not in expected]
    return "\n".join(
        [
            f"{path}: schema does not match first input",
            f"expected columns: {', '.join(expected)}",
            f"observed columns: {', '.join(observed)}",
            f"missing columns: {', '.join(missing) if missing else 'none'}",
            f"extra columns: {', '.join(extra) if extra else 'none'}",
        ],
    )


def validate_matching_schemas(paths: Sequence[Path]) -> None:
    if not paths:
        return
    expected = schema_columns(paths[0])
    for path in paths[1:]:
        observed = schema_columns(path)
        if observed != expected:
            raise ValueError(schema_mismatch_message(path, expected, observed))


def stack_big_tables(paths: Sequence[Path]) -> pl.LazyFrame:
    sorted_paths = sorted(paths)
    validate_matching_schemas(sorted_paths)
    return pl.concat([scan_tsv(path) for path in sorted_paths], how="vertical")


def input_paths(explicit_inputs: Sequence[Path], input_dir: Path | None) -> list[Path]:
    directory_inputs = sorted(input_dir.glob("*.tsv")) if input_dir else []
    return sorted([*explicit_inputs, *directory_inputs])


def write_empty(path: Path) -> None:
    path.write_text("", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    paths = input_paths(args.input, args.input_dir)
    if not paths:
        write_empty(args.output)
        return
    stack_big_tables(paths).collect().write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
