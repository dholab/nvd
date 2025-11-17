#!/usr/bin/env python3

import argparse
import csv
import os
import subprocess
import sys
import tarfile
import threading
import time
import urllib.request
from collections import defaultdict

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
TAXDB_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"


def format_duration(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    if h:
        return f"{int(h)}h{int(m):02d}m{int(s):02d}s"
    if m:
        return f"{int(m)}m{int(s):02d}s"
    return f"{s:.2f}s"


def spin(done_event, prefix=""):
    spinner = "|/-\\"
    idx = 0
    while not done_event.is_set():
        sys.stdout.write(f"\r{prefix} {spinner[idx % len(spinner)]}")
        sys.stdout.flush()
        idx += 1
        time.sleep(0.1)
    sys.stdout.write("\r" + " " * (len(prefix) + 2) + "\r")
    sys.stdout.flush()


def fetch_taxdump(taxdir):
    nodes = os.path.join(taxdir, "nodes.dmp")
    names = os.path.join(taxdir, "names.dmp")
    if os.path.exists(nodes) and os.path.exists(names):
        return
    archive = os.path.join(taxdir, "taxdump.tar.gz")
    urllib.request.urlretrieve(TAXDUMP_URL, archive)
    with tarfile.open(archive, "r:gz") as tar:
        for member in ("nodes.dmp", "names.dmp"):
            tar.extract(member, path=taxdir)
    os.remove(archive)


def fetch_taxdb(taxdir):
    needed = {"taxdb.bti", "taxdb.btd"}
    present = set(os.listdir(taxdir)) & needed
    if present == needed:
        return
    archive = os.path.join(taxdir, "taxdb.tar.gz")
    urllib.request.urlretrieve(TAXDB_URL, archive)
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            name = os.path.basename(member.name)
            if name in needed:
                tar.extract(member, path=taxdir)
    os.remove(archive)


def load_descendants(nodes_file, target_taxid):
    children = defaultdict(list)
    with open(nodes_file) as fh:
        for line in fh:
            cols = [c.strip() for c in line.split("|")]
            if len(cols) < 2:
                continue
            child, parent = cols[0], cols[1]
            children[parent].append(child)

    desc = {target_taxid}
    stack = [target_taxid]
    while stack:
        node = stack.pop()
        for c in children.get(node, []):
            if c not in desc:
                desc.add(c)
                stack.append(c)
    return desc


def run_blast(query, db, threads, word_size, max_seqs, out_file):
    cmd = [
        "blastn",
        "-task",
        "megablast",
        "-query",
        query,
        "-db",
        db,
        "-num_threads",
        str(threads),
        "-word_size",
        str(word_size),
    ]
    # only add max_target_seqs if user specified
    if max_seqs is not None:
        cmd.extend(["-max_target_seqs", str(max_seqs)])
    cmd.extend(
        [
            "-outfmt",
            "6 qseqid sseqid stitle pident qcovhsp evalue bitscore staxids sscinames",
        ],
    )
    done = threading.Event()
    spinner_thread = threading.Thread(target=spin, args=(done, "BLASTing…"))
    spinner_thread.start()

    t0 = time.time()
    with open(out_file, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)
    duration = time.time() - t0

    done.set()
    spinner_thread.join()
    print(f"[+] BLAST finished in {format_duration(duration)}")
    return duration


def classify_hits(hits_file, descendants, class_out):
    data = defaultdict(list)
    with open(hits_file) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            qid, e, taxstr = row[0], float(row[5]), row[7]
            tids = taxstr.split(";")
            data[qid].append((e, tids))

    unambig = ambig = 0
    t0 = time.time()
    with open(class_out, "w", newline="") as outfh:
        writer = csv.writer(outfh, delimiter="\t")
        writer.writerow(["qseqid", "status", "best_evalue", "best_taxids"])
        for qid, hits in data.items():
            min_ev = min(e for e, _ in hits)
            best = {t for e, tids in hits if e == min_ev for t in tids}
            if best <= descendants:
                status, unambig = "unambiguous", unambig + 1
            else:
                status, ambig = "ambiguous", ambig + 1
            writer.writerow([qid, status, min_ev, ";".join(sorted(best))])
    elapsed = time.time() - t0
    print(f"[+] Classification done in {format_duration(elapsed)}")
    return unambig, ambig


def main():
    p = argparse.ArgumentParser(
        description="Megablast + taxonomic unambiguous/ambiguous classification",
    )
    p.add_argument("--query", required=True, help="FASTA input for BLAST")
    p.add_argument(
        "--db", default="ref/core_nt", help="BLAST DB prefix (e.g. ref/core_nt)"
    )
    p.add_argument("--threads", type=int, default=8, help="BLAST threads")
    p.add_argument("--word_size", type=int, default=64, help="BLAST word size")
    p.add_argument(
        "--max_target_seqs",
        type=int,
        default=None,
        help="Maximum target sequences to retrieve (omit for BLAST default)",
    )
    p.add_argument(
        "--target-taxid",
        required=True,
        help="Parent taxID to include (plus all descendants)",
    )
    p.add_argument("--hits-out", default="hits.tsv", help="Intermediate BLAST output")
    p.add_argument(
        "--class-out",
        default="classification.tsv",
        help="Per-query classification output",
    )
    args = p.parse_args()

    taxdir = os.path.dirname(args.db) or "."

    print("Fetching taxonomy dump…")
    t0 = time.time()
    fetch_taxdump(taxdir)
    print(f"[+] Taxdump ready in {format_duration(time.time() - t0)}")

    print("Fetching BLAST taxdb…")
    t0 = time.time()
    fetch_taxdb(taxdir)
    print(f"[+] BLAST taxdb ready in {format_duration(time.time() - t0)}")

    nodes_f = os.path.join(taxdir, "nodes.dmp")
    print("Building descendant list…")
    t0 = time.time()
    descendants = load_descendants(nodes_f, args.target_taxid)
    print(f"[+] Descendants loaded in {format_duration(time.time() - t0)}")

    # BLAST with spinner + timer
    run_blast(
        args.query,
        args.db,
        args.threads,
        args.word_size,
        args.max_target_seqs,
        args.hits_out,
    )

    # classification with timer
    unambig, ambig = classify_hits(args.hits_out, descendants, args.class_out)

    total = unambig + ambig
    print(f"[+] All done. {total} queries: {unambig} unambiguous, {ambig} ambiguous")


if __name__ == "__main__":
    main()
