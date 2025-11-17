#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "lxml",
#     "snakemake",
#     "sqlite3",
# ]
# ///

import argparse
import collections
import logging
import os
import sqlite3
import sys
import urllib.request
from typing import Literal, Self

from lxml.builder import E

__version__ = "1.1"

# snakemake setup
MODE: Literal["snakemake", "commandline"]
if "snakemake" in globals():
    from snakemake.script import snakemake

    MODE = "snakemake"
else:
    MODE = "commandline"

# Configuration
TAXONOMY_URL = "https://sra-download.ncbi.nlm.nih.gov/traces/sra_references/tax_analysis/gettax.sqlite"
PROGRESS_GRANULARITY = 100000

# Set up logging
logger = logging.getLogger("tax_analysis")

# Data structures
Taxon = collections.namedtuple(
    "Taxon", ["tax_id", "parent_tax_id", "rank", "scientific_name"]
)  # noqa: PYI024


class TaxonomyDatabase:
    def __init__(self, sqlite_cache: str) -> None:
        self.sqlite_cache = sqlite_cache
        self.conn = None

    def connect(self) -> Self | None:
        try:
            if not os.path.exists(self.sqlite_cache):
                logger.error(f"SQLite file does not exist: {self.sqlite_cache}")
                return None
            if not os.access(self.sqlite_cache, os.R_OK):
                logger.error(f"SQLite file is not readable: {self.sqlite_cache}")
                return None
            self.conn = sqlite3.connect(self.sqlite_cache)
        except sqlite3.Error:
            logger.exception("Error connecting to SQLite database")
            return None
        else:
            return self

    def close(self) -> None:
        if self.conn:
            self.conn.close()

    def __enter__(self) -> Self | None:
        return self.connect()

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def get_taxon(self, tax_id) -> Taxon | None:
        if not self.conn:
            logger.error("Database connection is not established")
            return None
        cur = self.conn.cursor()
        cur.execute("SELECT * FROM taxons WHERE tax_id = ?", [int(tax_id)])
        row = cur.fetchone()
        cur.close()
        return Taxon(*row) if row else None


def ensure_taxonomy_file(sqlite_cache: str | None) -> str | None:
    try:
        sqlite_cache = sqlite_cache if sqlite_cache else "gettax.sqlite"
        if not os.path.exists(sqlite_cache):
            logger.info(f"Downloading taxonomy file from {TAXONOMY_URL}")
            urllib.request.urlretrieve(TAXONOMY_URL, sqlite_cache)  # noqa: S310
            logger.info(f"Downloaded taxonomy file to {sqlite_cache}")

        if not os.access(sqlite_cache, os.R_OK | os.W_OK):
            logger.error(f"SQLite file is not readable/writable: {sqlite_cache}")
            return None

        # Ensure the directory exists
        dir_path = os.path.dirname(sqlite_cache)
        if dir_path and not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)
            logger.info(f"Created directory: {dir_path}")

    except Exception:
        logger.exception("Error ensuring taxonomy file")
        return None
    else:
        return sqlite_cache


def parse_input(file, tax_db, args):
    counter = collections.Counter()
    counter[1] = 0  # explicitly add root
    for tax_id in args.include_tax_id:
        counter[tax_id] = 0

    lineage_cache = {}

    if args.wgs_mode:
        parse_wgs_mode(file, counter)
    else:
        parse_standard_mode(
            file, counter, tax_db, lineage_cache, args.compact, args.collated
        )

    return counter


def parse_wgs_mode(file, counter):
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


def parse_standard_mode(file, counter, tax_db, lineage_cache, compact, collated):
    iterator = iterate_merged_spots_compact if compact else iterate_merged_spots
    for hits, copies in iterator(file, collated):
        if not hits:
            logger.warning("Empty hits set encountered")
            continue
        tax_id = deduce_tax_id(hits, lineage_cache, tax_db)
        if tax_id is None:
            logger.warning(f"Could not deduce tax_id for hits: {hits}")
            continue
        counter[tax_id] += copies


def iterate_merged_spots(file, collated):
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


def deduce_tax_id(hits, lineage_cache, tax_db):
    if len(hits) == 1:
        return next(iter(hits))

    lineages = [get_lineage(hit, lineage_cache, tax_db) for hit in hits]
    last_matching_tax_id = None
    for nodes in zip(*lineages, strict=False):
        if len(set(nodes)) > 1:
            return last_matching_tax_id
        last_matching_tax_id = nodes[0]
    return last_matching_tax_id


def get_lineage(tax_id, cache, tax_db):
    if tax_id in cache:
        return cache[tax_id]

    if tax_id in (0, 1):
        lineage = [1]
    else:
        tax_info = tax_db.get_taxon(tax_id)
        parent_tax_id = tax_info.parent_tax_id if tax_info else 0
        parent_lineage = get_lineage(parent_tax_id, cache, tax_db)
        lineage = parent_lineage + [tax_id]

    cache[tax_id] = lineage
    return lineage


def build_tree(counter, tax_db):
    logger.info("Building XML tree")
    nodes = {}
    while counter:
        tax_id = next(iter(counter))
        get_or_add_node(nodes, counter, tax_db, tax_id)

    root = nodes[1]
    logger.info("Calculating totals")
    calculate_total_counts(root)
    return root


def get_or_add_node(nodes, counter, tax_db, tax_id):
    if tax_id in nodes:
        return nodes[tax_id]

    count = counter.pop(tax_id, 0)
    node = E.taxon(tax_id=str(tax_id), self_count=str(count))

    tax_info = tax_db.get_taxon(tax_id)
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
        parent_node = get_or_add_node(nodes, counter, tax_db, parent_tax_id)
        parent_node.append(node)
    return node


def calculate_total_counts(node):
    total_count = int(node.attrib["self_count"])
    for child in node:
        calculate_total_counts(child)
        total_count += int(child.attrib["total_count"])
    node[:] = sorted(
        node, key=lambda child: int(child.attrib["total_count"]), reverse=True
    )
    node.attrib["total_count"] = str(total_count)


def format_tax_tree(tree, args):
    if len(tree) == 0:
        return "Tree is empty"
    root = tree[0]
    total_count = int(root.attrib["total_count"])

    res = [
        format_node(node, total_count, "", args)
        for node in sorted(
            root, key=lambda x: int(x.attrib["total_count"]), reverse=True
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
                node, key=lambda x: int(x.attrib["total_count"]), reverse=True
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
        node, key=lambda x: int(x.attrib["total_count"]), reverse=True
    ):
        yield from format_node(subnode, grand_total, offset + args.indent, args)


def check_cutoff(value, total, args):
    if value < args.cutoff_hit_count:
        return True
    rate = float(value) / total
    percent = rate < args.cutoff_percent
    return percent


def flatten(iterable):
    for item in iterable:
        if isinstance(item, (list, tuple)) or (
            hasattr(item, "__iter__") and not isinstance(item, str)
        ):
            yield from flatten(item)
        else:
            yield item


def pad_tree(lines: list[str], separator: str):
    lines = list(lines)
    split_lines = [line.split("\t", 1) for line in lines if isinstance(line, str)]
    if not split_lines:
        return []
    max_name_len = max(len(name) for name, _ in split_lines)
    for name, stats in split_lines:
        yield name.ljust(max_name_len) + separator + stats


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Taxonomic analysis tool")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose mode"
    )
    parser.add_argument(
        "-c", "--sqlite-cache", default=None, help="Path to the SQLite cache file"
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
        "--compact", action="store_true", help="Use compact input format"
    )
    parser.add_argument(
        "--wgs-mode", action="store_true", help="Use WGS mode for parsing"
    )
    parser.add_argument("--collated", action="store_true", help="The input is collated")
    parser.add_argument(
        "--indent", default="  ", help="Indentation string, default is two spaces"
    )
    parser.add_argument(
        "--separator",
        default="    ",
        help="Name/stats separator string, default is four spaces",
    )
    parser.add_argument(
        "--no-padding", action="store_true", help="Disable tree padding"
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
    args = parse_arguments()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)

    if "snakemake" in globals() and MODE == "snakemake":
        # Running within Snakemake
        args.sqlite_cache = snakemake.params.gettax_sqlite_path
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
        # Ensure gettax SQLite cache is set up
        args.sqlite_cache = ensure_taxonomy_file(args.sqlite_cache)
        if args.sqlite_cache is None:
            logger.error("Failed to ensure taxonomy file. Exiting.")
            sys.exit(1)

        logger.info(f"Working with SQLite file: {args.sqlite_cache}")
        logger.info(f"Current working directory: {os.getcwd()}")

        with TaxonomyDatabase(args.sqlite_cache) as tax_db:
            if tax_db is None:
                logger.error("Failed to connect to the taxonomy database. Exiting.")
                sys.exit(1)

            logger.info(f"Reading {input_file}")
            with open(input_file) as f:
                counter = parse_input(f, tax_db, args)

            if not counter:
                logger.error("No valid data found in input")
                sys.exit(1)

            xml_tree = build_tree(counter, tax_db)
            if xml_tree is None:
                logger.error("Failed to generate XML tree. No data processed.")
                sys.exit(1)

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
