# Virus Enrichment Index

NVD v3 uses a compact minimizer index to enrich reads and contigs for human-infecting viruses before the more expensive assembly and BLAST verification stages. The default index is derived from NCBI's SRA Taxonomic Analysis Tool (STAT) dense database, filtered to NVD's curated human-virus taxid list, and converted into [Deacon](https://github.com/bede/deacon) index format.

The enrichment index is a fast screening artifact, not the final taxonomic classifier. It helps NVD retain reads and contigs that are plausible human-virus signal while discarding much of the irrelevant metagenomic background. Candidate sequences are still assembled, searched with BLAST, and filtered with taxonomy-aware logic downstream.

```text
raw reads
  -> deacon enrichment using a STAT-derived index
  -> enriched reads
  -> assembly
  -> BLAST verification
  -> final virus calls
```

## Default release index

Most users should use the release index we distribute via our LabKey server, [human_infecting_viruses.k31w1.idx](https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v3.0.0/v3.0.0/human_infecting_viruses.k31w1.idx), rather than rebuilding it themselves. The v3 release index was generated on May 15, 2026. LabKey reports the artifact size at roughly 138 MB.

The md5 checksum is:

```text
a3e5332ab0acb3c95ac144cb246ba695
```

Use it in a params file with:

```yaml
virus_index: /path/to/human_infecting_viruses.k31w1.idx
virus_index_version: v3.0
virus_kmer_size: 31
virus_window_size: 1
virus_abs_threshold: 1
virus_rel_threshold: 0.0
```

Or pass the index path directly with the CLI:

```bash
nvd run \
  --params-file run.yaml \
  --virus-index /path/to/human_infecting_viruses.k31w1.idx
```

The index path must be visible to the machine or worker process that executes the pipeline. On distributed systems such as CHTC, prefer a shared filesystem path or a shared preset that records the correct cluster-visible path.

## Source data

The upstream source is the SRA Taxonomic Analysis Tool dense database. Note! The following links point directly to downloadable files and may start downloads if opened in a browser. The dense database is large at around 250 gigabytes.

- Dense database: [tree_filter.dbss](https://sra-download.ncbi.nlm.nih.gov/traces/sra_references/tax_analysis/tree_filter.dbss)
- Matching annotation file: [tree_filter.dbss.annotation](https://sra-download.ncbi.nlm.nih.gov/traces/sra_references/tax_analysis/tree_filter.dbss.annotation)

The STAT dense database contains minimizers selected from 64-base windows across the tree of life. NVD does not use the full tree-of-life database for enrichment. During conversion, it keeps only entries assigned to the curated human-virus taxid list committed at [assets/human_viruses_taxlist.txt.zst](../assets/human_viruses_taxlist.txt.zst). As of May 15, 2026, this list contains 179,281 taxids.

The source taxid list is zstd-compressed in the repository. The conversion tool expects a plain text file with one NCBI taxid per line, so the asset must be decompressed before conversion.

## Conversion workflow

The conversion from STAT to Deacon is performed by the [`rust-script`](https://rust-script.org/) program [bin/stat_to_deacon.rs](../bin/stat_to_deacon.rs). Because it is a rust-script, the file is directly executable when `rust-script` and Rust are available. Its lightweight Cargo dependency metadata lives in the file header, so users do not need to create a separate Cargo project just to reproduce the conversion.

The conversion process reads the [tree_filter.dbss](https://sra-download.ncbi.nlm.nih.gov/traces/sra_references/tax_analysis/tree_filter.dbss) and [tree_filter.dbss.annotation](https://sra-download.ncbi.nlm.nih.gov/traces/sra_references/tax_analysis/tree_filter.dbss.annotation) pair. The annotation file acts as a table of contents: for taxids outside the curated list, the converter skips the corresponding k-mer block; for taxids in the list, it reads, decodes, converts, canonicalizes, deduplicates, and writes the k-mers into a Deacon v3 index file.

A reproducible build should look roughly like this, with exact local paths adjusted as needed:

```bash
zstd -dc assets/human_viruses_taxlist.txt.zst > human_viruses_taxlist.txt

rust-script bin/stat_to_deacon.rs \
  --dbss tree_filter.dbss \
  --annotation tree_filter.dbss.annotation \
  --taxids human_viruses_taxlist.txt \
  --target-k 31 \
  --window-size 1 \
  --output human_infecting_viruses.k31w1.idx
```

The conversion source should be understood as the commit where [bin/stat_to_deacon.rs](../bin/stat_to_deacon.rs) was last modified.

## Why the release index is `k31w1`

The release filename includes `k31w1` because the Deacon-formatted artifact is stored and queried with k=31 and w=1. This does not mean the original STAT database was built from 31-base, 1-base-window minimizers. STAT selected its source minimizers from 64-base windows.

The `w=1` setting is a query-time compatibility choice. During NVD runs, Deacon cannot reproduce STAT's exact minimizer selection process from query reads. Setting Deacon's query window to 1 forces every query k-mer to be checked against the STAT-derived index, rather than applying Deacon's own minimizer selection scheme. In other words, `w=1` is about maximizing sensitivity when querying a STAT-derived artifact; it is not a statement about how STAT originally selected the source minimizers.

## Rebuilding or overriding the index

We repeat: most users should not need to rebuild the default index. Rebuilding may make sense when you need to test a new virus taxid list, target a narrower or broader host range, update to a newer STAT database, or tune Deacon parameters for an unusual sequencing strategy.

If you rebuild the index, record enough provenance for another user to repeat the build: source STAT database URL and checksum, annotation URL and checksum, source taxid list commit, [bin/stat_to_deacon.rs](../bin/stat_to_deacon.rs) commit, Rust and rust-script versions, full command line, output filename, output checksum, and build date.

## Interpretation caveats

The enrichment index trades breadth, speed, and sensitivity before the final BLAST verification stage. Users should treat it as a recall-oriented prefilter for the viruses represented in the source list, not as a standalone classifier or an exhaustive screen for all viruses.

Negative results should therefore be interpreted as negative under the configured enrichment and BLAST strategy, not as evidence that every possible virus is absent. Viruses that are highly divergent from represented references or outside the curated human-virus taxid list may be less likely to pass the enrichment step.
