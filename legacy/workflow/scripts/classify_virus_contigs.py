#!/usr/bin/env python3

import argparse
import logging
import sys
from collections import Counter
from hits_to_report import TaxonomyDatabase, ensure_taxonomy_file, get_lineage

logger = logging.getLogger('classify_virus_contigs')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Classify virus contigs at family level')
    parser.add_argument('-i', '--input', required=True, help='Input file from aligns_to.3.1.1')
    parser.add_argument('-o', '--output', required=True, help='Output file for family classifications')
    parser.add_argument('-c', '--sqlite-cache', default='./gettax.sqlite', help='Path to the SQLite cache file')
    parser.add_argument('-s', '--stringency', type=float, default=0.7, help='Stringency threshold for classification')
    return parser.parse_args()

def get_virus_families():
    return {
        "Adenoviridae", "Anelloviridae", "Arenaviridae", "Arteriviridae", "Astroviridae",
        "Bornaviridae", "Peribunyaviridae", "Caliciviridae", "Coronaviridae", "Filoviridae",
        "Flaviviridae", "Hepadnaviridae", "Hepeviridae", "Herpesviridae", "Orthomyxoviridae",
        "Papillomaviridae", "Paramyxoviridae", "Parvoviridae", "Picobirnaviridae", "Picornaviridae",
        "Pneumoviridae", "Polyomaviridae", "Poxviridae", "Sedoreoviridae", "Retroviridae",
        "Rhabdoviridae", "Togaviridae"
    }

def classify_contig(tax_counts, virus_families, tax_db, stringency):
    total_hits = sum(tax_counts.values())
    family_hits = Counter()

    for tax_id, count in tax_counts.items():
        lineage = get_lineage(tax_id, {}, tax_db)
        for node in lineage:
            taxon = tax_db.get_taxon(node)
            if taxon and taxon.rank == 'family' and taxon.scientific_name in virus_families:
                family_hits[taxon.scientific_name] += count
                break

    if family_hits:
        top_family, top_count = family_hits.most_common(1)[0]
        if top_count / total_hits >= stringency:
            return top_family

    return "Unclassified"

def main():
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO)

    virus_families = get_virus_families()
    args.sqlite_cache = ensure_taxonomy_file(args.sqlite_cache)

    with TaxonomyDatabase(args.sqlite_cache) as tax_db:
        with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                contig_id = parts[0]
                tax_counts = Counter()
                
                for hit in parts[1:]:
                    if 'x' in hit:
                        tax_id, count = hit.split('x')
                        tax_counts[int(tax_id)] += int(count)
                    else:
                        tax_counts[int(hit)] += 1
                
                family = classify_contig(tax_counts, virus_families, tax_db, args.stringency)
                outfile.write(f"{contig_id}\t{family}\n")

    logger.info("Classification complete")

if __name__ == '__main__':
    main()