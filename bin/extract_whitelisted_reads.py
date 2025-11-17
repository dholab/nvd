#!/usr/bin/env python3

import argparse
from collections import defaultdict
import re
import os


def load_nodes_dmp(nodes_path):
    """Load taxonomy tree from nodes.dmp file"""
    child_to_parent = {}
    parent_to_children = defaultdict(list)
    rank_map = {}

    with open(nodes_path) as f:
        for line in f:
            parts = [p.strip() for p in line.strip().split("|")]
            if len(parts) < 3:
                continue
            taxid, parent, rank = int(parts[0]), int(parts[1]), parts[2]
            child_to_parent[taxid] = parent
            parent_to_children[parent].append(taxid)
            rank_map[taxid] = rank

    return parent_to_children, rank_map, child_to_parent


def get_descendants(taxid, tree, rank_map, target_ranks=None):
    """Get all descendants of a taxid, optionally filtered by rank"""
    descendants = set()
    stack = [taxid]

    while stack:
        current = stack.pop()
        if target_ranks is None or rank_map.get(current) in target_ranks:
            descendants.add(current)
        stack.extend(tree.get(current, []))

    return descendants


def extract_taxid(header):
    """Extract taxid from FASTA header"""
    match = re.search(r"TAXID=(\d+)", header)
    return int(match.group(1)) if match else None


def parse_fasta(fasta_path):
    """Parse FASTA file and yield (header, sequence, taxid) tuples"""
    with open(fasta_path) as f:
        current_header = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    taxid = extract_taxid(current_header)
                    yield (current_header, "".join(current_seq), taxid)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_header:
            taxid = extract_taxid(current_header)
            yield (current_header, "".join(current_seq), taxid)


def group_reads_by_taxid(fasta_path, allowed_taxids):
    """Group reads by their taxid"""
    reads_by_taxid = defaultdict(list)

    for header, sequence, taxid in parse_fasta(fasta_path):
        if taxid in allowed_taxids:
            reads_by_taxid[taxid].append((header, sequence))

    return reads_by_taxid


def write_strain_fastas(reads_by_taxid, sample_id, output_dir="."):
    """Write separate FASTA files for each strain taxid"""
    written_files = []

    for taxid, reads in reads_by_taxid.items():
        if not reads:
            continue

        output_file = os.path.join(
            output_dir, f"{sample_id}-{taxid}.taxon_specific.fasta"
        )

        with open(output_file, "w") as out:
            for header, sequence in reads:
                out.write(f"{header}\n")
                # Write sequence in 80-character lines
                for i in range(0, len(sequence), 80):
                    out.write(f"{sequence[i : i + 80]}\n")

        written_files.append((taxid, output_file))
        print(f"Wrote {len(reads)} reads for taxid {taxid} to {output_file}")

    return written_files


def main():
    parser = argparse.ArgumentParser(
        description="Extract reads matching strains under a parent taxid and create per-strain FASTAs"
    )
    parser.add_argument(
        "--taxid", type=int, required=True, help="Parent taxid to expand"
    )
    parser.add_argument("--nodes", required=True, help="Path to nodes.dmp")
    parser.add_argument("--fasta", required=True, help="Path to extracted FASTA")
    parser.add_argument(
        "--output", required=True, help="Output FASTA with all matching reads"
    )
    parser.add_argument(
        "--taxid-list",
        required=True,
        help="Path to write list of descendant strain taxids",
    )
    parser.add_argument(
        "--sample-id", required=True, help="Sample ID for naming output files"
    )
    parser.add_argument(
        "--output-dir", default=".", help="Directory for strain-specific FASTA files"
    )

    args = parser.parse_args()

    # Load taxonomy tree
    print(f"Loading taxonomy from {args.nodes}")
    tree, rank_map, child_to_parent = load_nodes_dmp(args.nodes)

    # Get all descendants of the parent taxid
    print(f"Finding descendants of taxid {args.taxid}")
    all_descendants = get_descendants(args.taxid, tree, rank_map)

    # First, try to find strains (including 'no rank' which is common for bacterial/viral strains)
    strain_ranks = {"strain", "no rank", "isolate", "serotype", "serovar"}
    strain_taxids = [
        tid for tid in all_descendants if rank_map.get(tid) in strain_ranks
    ]

    # If no strains found, fall back to species level
    if not strain_taxids:
        species_ranks = {"species", "subspecies"}
        species_taxids = [
            tid for tid in all_descendants if rank_map.get(tid) in species_ranks
        ]
        strain_species_taxids = sorted(species_taxids)
        print(
            f"No strains found under taxid {args.taxid}, using {len(species_taxids)} species-level taxids"
        )
    else:
        strain_species_taxids = sorted(strain_taxids)
        print(f"Found {len(strain_taxids)} strain-level taxids under {args.taxid}")

    print(f"Total taxids to process: {len(strain_species_taxids)}")

    # Write taxid list
    with open(args.taxid_list, "w") as tlist:
        for tid in strain_species_taxids:
            rank = rank_map.get(tid, "unknown")
            tlist.write(f"{tid}\t{rank}\n")

    # Group reads by taxid
    print(f"Extracting reads from {args.fasta}")
    reads_by_taxid = group_reads_by_taxid(args.fasta, set(strain_species_taxids))

    # Write combined output (all matching reads)
    all_reads = []
    for taxid_reads in reads_by_taxid.values():
        all_reads.extend(taxid_reads)

    with open(args.output, "w") as out:
        for header, sequence in all_reads:
            out.write(f"{header}\n")
            for i in range(0, len(sequence), 80):
                out.write(f"{sequence[i : i + 80]}\n")

    print(f"Wrote {len(all_reads)} total reads to {args.output}")

    # Write per-strain FASTA files
    written_files = write_strain_fastas(reads_by_taxid, args.sample_id, args.output_dir)

    if not written_files:
        print("No strain-specific reads found")
    else:
        print(f"Created {len(written_files)} strain-specific FASTA files")


if __name__ == "__main__":
    main()
