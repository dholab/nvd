#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "lxml",
#     "snakemake",
# ]
# ///

import argparse
import collections
import logging
import sys
from typing import Literal

from lxml.builder import E
from py_nvd import __version__, taxonomy

# snakemake setup
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"

# Configuration
PROGRESS_GRANULARITY = 100000

# Set up logging
logger = logging.getLogger("tax_analysis")


def parse_input(file, tax, args):
    """Parse input file and count taxids."""
    counter = collections.Counter()
    counter[1] = 0  # explicitly add root
    for tax_id in args.include_tax_id:
        counter[tax_id] = 0

    if args.wgs_mode:
        parse_wgs_mode(file, counter)
    else:
        parse_standard_mode(file, counter, tax, args.compact, args.collated)

    return counter


def parse_wgs_mode(file, counter):
    """Parse WGS mode input format."""
    for line in file:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            logger.warning(f"Skipping malformed line: {line}")
            continue
        hits = parts[1:]
        if not hits:
            logger.warning(f"No hits for spot: {parts[0]}")
            continue
        for hit in hits:
            if "x" in hit:
                tax_id, count = map(int, hit.split("x"))
            else:
                tax_id, count = int(hit), 1
            counter[tax_id] += count


def parse_standard_mode(file, counter, tax, compact, collated):
    """Parse standard mode input format."""
    iterator = iterate_merged_spots_compact if compact else iterate_merged_spots
    for hits, copies in iterator(file, collated):
        if not hits:
            logger.warning("Empty hits set encountered")
            continue
        tax_id = deduce_tax_id(hits, tax)
        if tax_id is None:
            logger.warning(f"Could not deduce tax_id for hits: {hits}")
            continue
        counter[tax_id] += copies


def iterate_merged_spots(file, collated):
    """Iterate over merged spots in standard format."""
    last_spot = None
    last_hits = None
    for line in file:
        parts = line.split("\t")
        if len(parts) < 2:
            logger.warning(f"Skipping malformed line: {line}")
            continue
        spot = parts[0]
        try:
            hits = set(int(p.split("x")[0]) for p in parts[1:])
        except ValueError as e:
            logger.warning(f"Error parsing hits in line: {line}. Error: {e}")
            continue
        if not hits:
            logger.warning(f"No hits for spot: {spot}")
            continue
        if collated:
            yield hits, 1
        elif spot == last_spot:
            last_hits |= hits
        else:
            if last_spot:
                yield last_hits, 1
            last_spot = spot
            last_hits = hits
    if last_spot:
        yield last_hits, 1


def iterate_merged_spots_compact(file, collated):
    """Iterate over merged spots in compact format."""
    last_spot = None
    last_hits = None
    for line in file:
        if not line.strip():
            continue
        if line[0] != "\t":
            if last_spot:
                last_spot = None
                yield last_hits, 1
            parts = line.split("\t")
            if len(parts) < 2:
                logger.warning(f"Skipping malformed line: {line}")
                continue
            try:
                copies = int(parts[0])
                hits = set(int(p.split("x")[0]) for p in parts[1:])
            except ValueError as e:
                logger.warning(f"Error parsing line: {line}. Error: {e}")
                continue
            yield hits, copies
        else:
            line = line[1:]
            parts = line.split("\t")
            if len(parts) < 2:
                logger.warning(f"Skipping malformed line: {line}")
                continue
            spot = parts[0]
            try:
                hits = set(int(p.split("x")[0]) for p in parts[1:])
            except ValueError as e:
                logger.warning(f"Error parsing hits in line: {line}. Error: {e}")
                continue
            if not hits:
                logger.warning(f"No hits for spot: {spot}")
                continue
            if collated:
                yield hits, 1
            elif spot == last_spot:
                last_hits |= hits
            else:
                if last_spot:
                    yield last_hits, 1
                last_spot = spot
                last_hits = hits
    if last_spot:
        yield last_hits, 1


def deduce_tax_id(hits, tax):
    """
    Deduce the consensus taxid for a set of hits using LCA.

    Args:
        hits: Set of taxids from BLAST hits
        tax: TaxonomyDB instance

    Returns:
        The LCA taxid, or None if no valid taxids
    """
    if len(hits) == 1:
        return next(iter(hits))

    # Use taxonomy module's find_lca for consistency with annotate_blast_lca.py
    return tax.find_lca(list(hits))


def build_tree(counter, tax):
    """Build XML tree from taxid counts."""
    logger.info("Building XML tree")
    nodes = {}
    while counter:
        tax_id = next(iter(counter))
        get_or_add_node(nodes, counter, tax, tax_id)

    root = nodes[1]
    logger.info("Calculating totals")
    calculate_total_counts(root)
    return root


def get_or_add_node(nodes, counter, tax, tax_id):
    """Get or create a node in the tree."""
    if tax_id in nodes:
        return nodes[tax_id]

    count = counter.pop(tax_id, 0)
    node = E.taxon(tax_id=str(tax_id), self_count=str(count))

    tax_info = tax.get_taxon(tax_id)
    if tax_info:
        if tax_info.rank:
            node.attrib["rank"] = tax_info.rank
        node.attrib["name"] = tax_info.scientific_name
        parent_tax_id = tax_info.parent_tax_id
    else:
        node.attrib["name"] = "unknown"
        parent_tax_id = 1

    nodes[tax_id] = node
    if tax_id != 1:
        parent_node = get_or_add_node(nodes, counter, tax, parent_tax_id)
        parent_node.append(node)
    return node


def calculate_total_counts(node):
    """Calculate total counts for each node in the tree."""
    total_count = int(node.attrib["self_count"])
    for child in node:
        calculate_total_counts(child)
        total_count += int(child.attrib["total_count"])
    node[:] = sorted(
        node,
        key=lambda child: int(child.attrib["total_count"]),
        reverse=True,
    )
    node.attrib["total_count"] = str(total_count)


def format_tax_tree(tree, args):
    """Format the taxonomy tree for output."""
    if len(tree) == 0:
        return "Tree is empty"
    root = tree[0]
    total_count = int(root.attrib["total_count"])

    res = [
        format_node(node, total_count, "", args)
        for node in sorted(
            root,
            key=lambda x: int(x.attrib["total_count"]),
            reverse=True,
        )
    ]
    res = list(flatten(res))
    if args.no_padding:
        res = [
            f"{name}{args.separator}{stats}"
            for name, stats in (
                line.split("\t", 1) for line in res if isinstance(line, str)
            )
        ]
    else:
        res = list(pad_tree(res, args.separator))
    return "\n".join(res)


def format_node(node, grand_total, offset, args):
    """Format a single node for output."""
    self_count = int(node.attrib["self_count"])
    if len(node) == 1 and check_cutoff(self_count, grand_total, args):
        if args.skip == "none":
            skip = False
        elif args.skip == "unranked":
            skip = node.attrib.get("rank") is None
        elif args.skip == "all":
            skip = True
        else:
            assert False
        if skip:
            for subnode in sorted(
                node,
                key=lambda x: int(x.attrib["total_count"]),
                reverse=True,
            ):
                yield from format_node(subnode, grand_total, offset, args)
            return

    name = node.attrib["name"]
    total_count = int(node.attrib["total_count"])

    if check_cutoff(total_count, grand_total, args):
        return

    rate = float(total_count) / grand_total
    percent = rate * 100

    if args.no_padding:
        percent_precision = f".{args.precision}"
        hits_precision = ""
    else:
        percent_precision = f"{args.precision + 3}.{args.precision}"
        hits_precision = str(len(str(grand_total)))

    pattern = f"%s%s\t%{percent_precision}f%%  (%{hits_precision}d hits)"
    yield pattern % (offset, name, percent, total_count)

    for subnode in sorted(
        node,
        key=lambda x: int(x.attrib["total_count"]),
        reverse=True,
    ):
        yield from format_node(subnode, grand_total, offset + args.indent, args)


def check_cutoff(value, total, args):
    """Check if a value is below the cutoff threshold."""
    if value < args.cutoff_hit_count:
        return True
    rate = float(value) / total
    percent = rate < args.cutoff_percent
    return percent


def flatten(iterable):
    """Flatten nested iterables."""
    for item in iterable:
        if isinstance(item, (list, tuple)) or (
            hasattr(item, "__iter__") and not isinstance(item, str)
        ):
            yield from flatten(item)
        else:
            yield item


def pad_tree(lines: list[str], separator: str):
    """Pad tree lines for aligned output."""
    lines = list(lines)
    split_lines = [line.split("\t", 1) for line in lines if isinstance(line, str)]
    if not split_lines:
        return []
    max_name_len = max(len(name) for name, _ in split_lines)
    for name, stats in split_lines:
        yield name.ljust(max_name_len) + separator + stats


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Taxonomic analysis tool")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose mode",
    )
    # Keep --sqlite-cache for backward compatibility with Nextflow,
    # but it's no longer used (taxonomy.open() handles DB location)
    parser.add_argument(
        "-c",
        "--sqlite-cache",
        default=None,
        help="(Deprecated) Path to SQLite cache file - now handled automatically",
    )
    parser.add_argument(
        "-i",
        "--include-tax-id",
        type=int,
        action="append",
        default=[],
        help="Include taxon into tax tree even if it has no hits",
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Use compact input format",
    )
    parser.add_argument(
        "--wgs-mode",
        action="store_true",
        help="Use WGS mode for parsing",
    )
    parser.add_argument("--collated", action="store_true", help="The input is collated")
    parser.add_argument(
        "--indent",
        default="  ",
        help="Indentation string, default is two spaces",
    )
    parser.add_argument(
        "--separator",
        default="    ",
        help="Name/stats separator string, default is four spaces",
    )
    parser.add_argument(
        "--no-padding",
        action="store_true",
        help="Disable tree padding",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=2,
        help="Number of digits after decimal point, default is 2",
    )
    parser.add_argument(
        "--cutoff-percent",
        type=float,
        default=0.01,
        help="Cutoff percent, default is 0.01",
    )
    parser.add_argument(
        "--cutoff-hit-count",
        type=int,
        default=0,
        help="Cutoff hit count, disabled by default",
    )
    parser.add_argument(
        "--skip",
        choices=["all", "unranked", "none"],
        default="all",
        help="Skip nodes with only one child and less than cutoff exact hits",
    )
    parser.add_argument("input_file", nargs="?", help="Path to the input file")
    parser.add_argument("output_file", nargs="?", help="Path to the output file")
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_arguments()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)

    if "snakemake" in globals() and MODE == "snakemake":
        # Running within Snakemake
        # args.sqlite_cache from snakemake is ignored - taxonomy.open() handles it
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
    else:
        # Running from the command line
        if not args.input_file or not args.output_file:
            logger.error(
                "Input and output files must be specified when running from the command line",
            )
            sys.exit(1)
        input_file = args.input_file
        output_file = args.output_file

    try:
        with taxonomy.open() as tax:
            logger.info(f"Reading {input_file}")
            with open(input_file) as f:
                counter = parse_input(f, tax, args)

            if not counter:
                logger.warning("No valid data found in input. Creating empty output.")
                with open(output_file, "w") as f:
                    f.write("")
                return

            xml_tree = build_tree(counter, tax)
            if xml_tree is None:
                logger.warning("No data to build tree. Creating empty output.")
                with open(output_file, "w") as f:
                    f.write("")
                return

            xml_root = E.taxon_tree(xml_tree, parser_version=__version__)
            formatted_tree = format_tax_tree(xml_root, args)

            logger.info(f"Writing output to {output_file}")
            with open(output_file, "w") as f:
                f.write(formatted_tree)

    except Exception:
        logger.exception("An unexpected error occurred. Detailed traceback:")
        sys.exit(1)


if __name__ == "__main__":
    main()
