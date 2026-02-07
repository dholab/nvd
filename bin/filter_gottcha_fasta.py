#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict


def filter_fasta_strain_or_species(input_fasta: str) -> str:
    lines = input_fasta.strip().splitlines()
    record_map = defaultdict(dict)  # {read_id: {level: [(header, seq_lines)]}}

    current_header = None
    current_seq = []

    for line in lines:
        if line.startswith(">"):
            if current_header:
                read_id = current_header.split("|")[0]  # everything up to first |
                level_match = re.search(r"LEVEL=(\w+)", current_header)
                if level_match:
                    level = level_match.group(1).lower()
                    record_map[read_id].setdefault(level, []).append(
                        (current_header, current_seq),
                    )
            current_header = line.strip()
            current_seq = []
        else:
            current_seq.append(line.strip())

    if current_header:
        read_id = current_header.split("|")[0]
        level_match = re.search(r"LEVEL=(\w+)", current_header)
        if level_match:
            level = level_match.group(1).lower()
            record_map[read_id].setdefault(level, []).append(
                (current_header, current_seq),
            )

    output_lines = []
    for levels in record_map.values():
        if "strain" in levels:
            selected = levels["strain"]
        elif "species" in levels:
            selected = levels["species"]
        else:
            continue
        for header, seq_lines in selected:
            output_lines.append(header)
            output_lines.extend(seq_lines)

    return "\n".join(output_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GOTTCHA2-style FASTA to strain/species only.",
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output FASTA file (filtered)",
    )

    args = parser.parse_args()

    with open(args.input) as infile:
        fasta_data = infile.read()

    filtered = filter_fasta_strain_or_species(fasta_data)

    with open(args.output, "w") as outfile:
        outfile.write(filtered)


if __name__ == "__main__":
    main()
