#!/usr/bin/env python3
"""Render a merged multi-sample taxonomic profile sunburst HTML report."""

from __future__ import annotations

import argparse
import csv
import html
from dataclasses import dataclass, field
from pathlib import Path

from taxburst.output import env

RANKS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
)


@dataclass
class TreeNode:
    name: str
    rank: str
    counts: list[float]
    children: dict[str, TreeNode] = field(default_factory=dict)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render one multi-sample Krona-style sunburst from summary CSVs.",
    )
    parser.add_argument(
        "--summary",
        action="append",
        required=True,
        help="Sample summary as SAMPLE_ID=PATH. May be supplied more than once.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output HTML path.",
    )
    return parser.parse_args()


def parse_summary_arg(value: str) -> tuple[str, Path]:
    sample_id, separator, path = value.partition("=")
    if not separator or not sample_id or not path:
        msg = f"Summary must be formatted as SAMPLE_ID=PATH, got: {value}"
        raise ValueError(msg)
    return sample_id, Path(path)


def taxburst_ignores(row: dict[str, str]) -> bool:
    return row.get("lineage") == "unclassified" and row.get("rank") != "superkingdom"


def rank_for_depth(depth: int) -> str:
    rank_index = depth - 1
    if 0 <= rank_index < len(RANKS):
        return RANKS[rank_index]
    return f"rank {depth}"


def row_count(row: dict[str, str]) -> float:
    return 1000 * float(row["f_weighted_at_rank"])


def add_lineage(
    root: TreeNode,
    *,
    lineage: str,
    rank: str,
    count: float,
    sample_index: int,
) -> None:
    node = root
    parts = lineage.split(";")
    for depth, part in enumerate(parts, start=1):
        node = node.children.setdefault(
            part,
            TreeNode(
                name=part,
                rank=rank_for_depth(depth),
                counts=[0.0] * len(root.counts),
            ),
        )
    node.rank = rank
    node.counts[sample_index] = count


def build_tree(samples: list[tuple[str, Path]]) -> TreeNode:
    root = TreeNode(name="all samples", rank="root", counts=[0.0] * len(samples))
    for sample_index, (_sample_id, summary_path) in enumerate(samples):
        with summary_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if taxburst_ignores(row):
                    continue
                add_lineage(
                    root,
                    lineage=row["lineage"],
                    rank=row["rank"],
                    count=row_count(row),
                    sample_index=sample_index,
                )
    root.counts = [
        sum(child.counts[sample_index] for child in root.children.values())
        for sample_index in range(len(samples))
    ]
    return root


def vals(values: list[float]) -> str:
    return "".join(f"<val>{value:.1f}</val>" for value in values)


def empty_vals(count: int) -> str:
    return "".join("<val></val>" for _ in range(count))


def render_node(node: TreeNode, *, indent: int = 2) -> str:
    space = "  " * indent
    pieces = [
        f'{space}<node name="{html.escape(node.name, quote=True)}">',
        f"{space}  <members>{empty_vals(len(node.counts))}</members>",
        f"{space}  <rank><val>{html.escape(node.rank)}</val></rank>",
        f"{space}  <count>{vals(node.counts)}</count>",
    ]
    pieces.extend(
        render_node(child, indent=indent + 1)
        for child in sorted(node.children.values(), key=lambda child: child.name)
    )
    pieces.append(f"{space}</node>")
    return "\n".join(pieces)


def render_krona_body(root: TreeNode, sample_ids: list[str]) -> str:
    dataset_xml = "\n".join(
        f"    <dataset>{html.escape(sample_id)}</dataset>" for sample_id in sample_ids
    )
    return f"""<krona collapse="true" key="false">
  <attributes magnitude="count">
    <data>members</data>
    <attribute display="Count" dataAll="members">count</attribute>
    <attribute display="Unassigned" dataNode="members">unassigned</attribute>
    <attribute display="Rank" mono="true">rank</attribute>
  </attributes>
  <datasets>
{dataset_xml}
  </datasets>
  <color attribute="count" hueStart="0" hueEnd="120" valueStart="0.413" valueEnd="0.932367571822115" default="false" ></color>
{render_node(root, indent=1)}
</krona>
"""


def render_html(samples: list[tuple[str, Path]]) -> str:
    root = build_tree(samples)
    header = env.get_template("krona-header.html").render()
    footer = env.get_template("krona-footer.html").render()
    return (
        header
        + "\n<style>#lastDataset { font: 11px sans-serif !important; }</style>\n"
        + "\n"
        + render_krona_body(root, [sample[0] for sample in samples])
        + footer
    )


def main() -> None:
    args = parse_args()
    samples = [parse_summary_arg(value) for value in args.summary]
    args.output.write_text(render_html(samples), encoding="utf-8")


if __name__ == "__main__":
    main()
