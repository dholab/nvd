#!/usr/bin/env python3
"""Normalize taxonomic profile summary display names for visualization tools.

The current input contract is a CSV with at least ``rank`` and ``lineage``
columns, matching sourmash ``csv_summary`` output. The normalizer preserves all
columns and row order, changing only lineage display names when a visualization
tool would otherwise see duplicate taxon node names in different tree branches.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path

STANDARD_RANKS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
)


@dataclass(frozen=True)
class Node:
    prefix: tuple[str, ...]
    rank: str


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Normalize a taxonomic profile summary CSV so taxon display names "
            "are globally unique for visualization tools."
        ),
    )
    subcommands = parser.add_subparsers(dest="command")

    normalize_parser = subcommands.add_parser(
        "normalize",
        help="Build a normalization map from one summary and apply it to that summary.",
    )
    normalize_parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input taxonomic profile summary CSV.",
    )
    normalize_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output CSV with normalized lineage display names.",
    )

    build_map_parser = subcommands.add_parser(
        "build-map",
        help="Build one normalization map from one or more summary CSVs.",
    )
    build_map_parser.add_argument(
        "--input",
        required=True,
        type=Path,
        action="append",
        help="Input taxonomic profile summary CSV. May be supplied more than once.",
    )
    build_map_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output normalization map JSON.",
    )

    apply_map_parser = subcommands.add_parser(
        "apply-map",
        help="Apply an existing normalization map to one summary CSV.",
    )
    apply_map_parser.add_argument(
        "--map",
        required=True,
        type=Path,
        help="Normalization map JSON from build-map.",
    )
    apply_map_parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input taxonomic profile summary CSV.",
    )
    apply_map_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output CSV with normalized lineage display names.",
    )
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        parser.error("a subcommand is required")
    if argv and argv[0].startswith("-"):
        argv = ["normalize", *argv]
    return parser.parse_args(argv)


def taxburst_ignores(row: dict[str, str]) -> bool:
    """Return whether taxburst skips this bare unclassified row."""
    return row.get("lineage") == "unclassified" and row.get("rank") != "superkingdom"


def split_lineage(lineage: str) -> tuple[str, ...]:
    if not lineage:
        return ()
    return tuple(lineage.split(";"))


def rank_for_depth(depth: int) -> str:
    rank_index = depth - 1
    if 0 <= rank_index < len(STANDARD_RANKS):
        return STANDARD_RANKS[rank_index]
    return f"rank {depth}"


def lineage_prefixes(lineage: str) -> list[tuple[str, ...]]:
    parts = split_lineage(lineage)
    return [parts[:depth] for depth in range(1, len(parts) + 1)]


def collect_nodes(rows: list[dict[str, str]]) -> OrderedDict[tuple[str, ...], Node]:
    nodes: OrderedDict[tuple[str, ...], Node] = OrderedDict()
    for row in rows:
        if taxburst_ignores(row):
            continue

        for prefix in lineage_prefixes(row.get("lineage", "")):
            if prefix not in nodes:
                nodes[prefix] = Node(prefix=prefix, rank=rank_for_depth(len(prefix)))

    return nodes


def short_hash(prefix: tuple[str, ...]) -> str:
    digest = hashlib.sha1(";".join(prefix).encode("utf-8"), usedforsecurity=False)
    return digest.hexdigest()[:8]


def original_parent_context(prefix: tuple[str, ...], context_depth: int) -> str:
    parent = prefix[:-1]
    if not parent:
        return "root"
    return " > ".join(parent[-context_depth:])


def candidate_name(node: Node, context_depth: int) -> str:
    name = node.prefix[-1]
    context = original_parent_context(node.prefix, context_depth)
    if context == "root":
        return f"{name} [{node.rank} at root]"
    return f"{name} [{node.rank} under {context}]"


def choose_duplicate_names(
    duplicate_nodes: list[Node],
    reserved_names: set[str],
) -> dict[tuple[str, ...], str]:
    """Choose short, deterministic display names for one duplicate-name group."""
    context_depths = {node.prefix: 1 for node in duplicate_nodes}

    while True:
        names_by_prefix = {
            node.prefix: candidate_name(node, context_depths[node.prefix])
            for node in duplicate_nodes
        }
        counts = Counter(names_by_prefix.values())
        colliding_prefixes = {
            prefix
            for prefix, name in names_by_prefix.items()
            if counts[name] > 1 or name in reserved_names
        }
        if not colliding_prefixes:
            return names_by_prefix

        can_expand = False
        for node in duplicate_nodes:
            if node.prefix not in colliding_prefixes:
                continue
            max_context_depth = max(1, len(node.prefix) - 1)
            if context_depths[node.prefix] < max_context_depth:
                context_depths[node.prefix] += 1
                can_expand = True

        if not can_expand:
            return {
                prefix: f"{name} #{short_hash(prefix)}"
                if counts[name] > 1 or name in reserved_names
                else name
                for prefix, name in names_by_prefix.items()
            }


def build_normalized_prefixes(
    nodes: OrderedDict[tuple[str, ...], Node],
) -> dict[tuple[str, ...], tuple[str, ...]]:
    nodes_by_name: defaultdict[str, list[Node]] = defaultdict(list)
    for node in nodes.values():
        nodes_by_name[node.prefix[-1]].append(node)

    duplicate_groups = {
        name: group for name, group in nodes_by_name.items() if len(group) > 1
    }
    reserved_names = {
        name for name, group in nodes_by_name.items() if name not in duplicate_groups
    }

    replacement_names: dict[tuple[str, ...], str] = {}
    for group in duplicate_groups.values():
        chosen_names = choose_duplicate_names(group, reserved_names)
        replacement_names.update(chosen_names)
        reserved_names.update(chosen_names.values())

    normalized: dict[tuple[str, ...], tuple[str, ...]] = {}
    for prefix in nodes:
        parent = prefix[:-1]
        normalized_parent = normalized.get(parent, ())
        normalized_name = replacement_names.get(prefix, prefix[-1])
        normalized[prefix] = (*normalized_parent, normalized_name)

    return normalized


def build_lineage_replacements(rows: list[dict[str, str]]) -> dict[str, str]:
    nodes = collect_nodes(rows)
    normalized_prefixes = build_normalized_prefixes(nodes)
    replacements: dict[str, str] = {}

    for row in rows:
        if taxburst_ignores(row):
            continue
        original = split_lineage(row.get("lineage", ""))
        normalized = normalized_prefixes.get(original)
        if normalized is not None and normalized != original:
            replacements[";".join(original)] = ";".join(normalized)

    return replacements


def validate_columns(fieldnames: list[str] | None) -> list[str]:
    if fieldnames is None:
        msg = "Input CSV is missing a header row."
        raise ValueError(msg)

    required = {"rank", "lineage"}
    missing = required.difference(fieldnames)
    if missing:
        joined = ", ".join(sorted(missing))
        msg = f"Input CSV is missing required column(s): {joined}"
        raise ValueError(msg)

    return fieldnames


def read_summary_rows(input_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with input_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = validate_columns(reader.fieldnames)
        rows = list(reader)
    return fieldnames, rows


def write_normalization_map(*, input_paths: list[Path], output_path: Path) -> None:
    all_rows: list[dict[str, str]] = []
    for input_path in input_paths:
        _fieldnames, rows = read_summary_rows(input_path)
        all_rows.extend(rows)

    replacements = build_lineage_replacements(all_rows)
    output_path.write_text(
        json.dumps(
            {
                "format": "nvd.taxonomic_profile_summary.normalization_map.v1",
                "replacements": replacements,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def read_normalization_map(map_path: Path) -> dict[str, str]:
    payload = json.loads(map_path.read_text(encoding="utf-8"))
    if payload.get("format") != "nvd.taxonomic_profile_summary.normalization_map.v1":
        msg = f"Unsupported normalization map format in {map_path}"
        raise ValueError(msg)
    replacements = payload.get("replacements")
    if not isinstance(replacements, dict):
        msg = f"Normalization map is missing replacements object in {map_path}"
        raise TypeError(msg)
    return {
        str(original): str(normalized) for original, normalized in replacements.items()
    }


def apply_lineage_replacements(
    *,
    input_path: Path,
    output_path: Path,
    replacements: dict[str, str],
) -> None:
    fieldnames, rows = read_summary_rows(input_path)
    if not replacements:
        shutil.copyfile(input_path, output_path)
        return

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            lineage = row.get("lineage", "")
            output_row = row
            if lineage in replacements:
                output_row = {**row, "lineage": replacements[lineage]}
            writer.writerow(output_row)


def normalize_taxonomic_profile_summary(*, input_path: Path, output_path: Path) -> None:
    _fieldnames, rows = read_summary_rows(input_path)
    replacements = build_lineage_replacements(rows)
    apply_lineage_replacements(
        input_path=input_path,
        output_path=output_path,
        replacements=replacements,
    )


def main() -> None:
    args = parse_args()
    if args.command in (None, "normalize"):
        normalize_taxonomic_profile_summary(
            input_path=args.input,
            output_path=args.output,
        )
    elif args.command == "build-map":
        write_normalization_map(input_paths=args.input, output_path=args.output)
    elif args.command == "apply-map":
        replacements = read_normalization_map(args.map)
        apply_lineage_replacements(
            input_path=args.input,
            output_path=args.output,
            replacements=replacements,
        )
    else:
        msg = f"Unknown command: {args.command}"
        raise ValueError(msg)


if __name__ == "__main__":
    main()
