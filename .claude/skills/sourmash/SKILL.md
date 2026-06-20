---
name: sourmash
description: Use for NVD sourmash and sourmash_plugin_branchwater sketch-based metagenomic FASTQ workflows: sketching, gather/fastgather, tax metagenome, signature/database compatibility, branchwater commands, and result interpretation.
---

# sourmash and branchwater for NVD

Use this skill when adding, reviewing, or debugging sourmash-based metagenomic analysis in NVD. In this project, sourmash work is experimental and should stay behind `params.experimental == true` until the release branch deliberately pulls it in. Think of the input reads as a query metagenome: sourmash sketches the query’s k-mers and compares that sketch to reference sketches; it does not classify individual reads.

The expected project versions are `sourmash >=4.9.4,<5` and `sourmash_plugin_branchwater >=0.9.15,<0.10`. The branchwater plugin adds Rust-backed commands under `sourmash scripts ...`. Its commands are much faster for large jobs, but the plugin docs are less strict about semver and output-column stability than sourmash core, so parse CSVs by header and check installed help when exact spelling matters.

## First checks

Prefer the project’s locked environment and recipes:

```bash
pixi run sourmash --version
pixi run sourmash scripts singlesketch --help
pixi run sourmash scripts fastgather --help
pixi run sourmash tax metagenome --help
```

If pipeline params or config are touched, use the repo recipes rather than ad hoc commands:

```bash
just schema-check
just config-check
```

## The smallest NVD slice

NVD post-preprocessing reads are already represented as one FASTQ per sample:

```text
tuple(sample_id, platform, read_structure, fastq)
```

Sketch that file directly. Pairing semantics do not matter to sourmash sketching; it hashes sequences.

```nextflow
process SOURMASH_SKETCH_QUERY_METAGENOME {
    tag "${sample_id}"
    label "low"

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.sourmash.query_metagenome.k31.scaled1000.abund.sig.gz"), emit: query_sketches

    script:
    """
    sourmash scripts singlesketch \\
        ${reads} \\
        --name "${sample_id}" \\
        -p dna,k=31,scaled=1000,abund \\
        -o "${sample_id}.sourmash.query_metagenome.k31.scaled1000.abund.sig.gz"
    """
}
```

Gate sourmash at the workflow callsite, not in every process:

```nextflow
if (params.experimental == true) {
    METAGENOME_PROFILING(PREPROCESS_READS.out.reads)
}
```

The subworkflow should use domain/product vocabulary:

```nextflow
workflow METAGENOME_PROFILING {
    emit:
    query_sketches = SOURMASH_SKETCH_QUERY_METAGENOME.out.query_sketches
}
```

Keep `SOURMASH_` on processes that directly invoke sourmash. Do not force the prefix onto future Python reporting/QC processes that merely consume sourmash sketches.

## Canonical query metagenome to profile recipe

Core sourmash commands:

```bash
sourmash sketch dna \
  -p k=31,scaled=1000,abund \
  -o sample.sig.zip \
  sample_R1.fastq.gz sample_R2.fastq.gz

sourmash gather \
  sample.sig.zip \
  database.zip \
  -k 31 \
  --threshold-bp 50000 \
  -o sample.gather.csv

sourmash tax metagenome \
  --gather-csv sample.gather.csv \
  --taxonomy taxonomy.csv \
  --use-abund \
  --output-format human csv_summary krona kreport \
  --rank species \
  --output-base sample.tax
```

Branchwater fast path:

```bash
sourmash scripts singlesketch \
  input.fastq.gz \
  --name sample_id \
  -p dna,k=31,scaled=1000,abund \
  -o sample.sig.gz

sourmash scripts fastgather \
  sample.sig.gz \
  database.zip \
  --cores 4 \
  -o sample.fastgather.csv
```

For many samples, consider `manysketch` and `fastmultigather`, but only after the workflow benefits from batching enough to justify the extra manifest management:

```bash
sourmash scripts manysketch samples.csv \
  -p k=31,scaled=1000,abund \
  -o samples.zip

sourmash scripts fastmultigather \
  queries.manifest.csv \
  database.zip \
  --cores 8 \
  --save-matches \
  -o all-samples.fastmultigather.csv
```

Be careful with `fastmultigather`: duplicate signature names can overwrite sidecar outputs. Always set unique sample names when sketching.

## Compatibility checklist

Before comparing, gathering, indexing, or summarizing, verify:

```text
query sketch    database sketch/index    taxonomy CSV
     |                 |                      |
     +-- same ksize ---+
     +-- same moltype -+
     +-- scaled compatible or lower query scaled
     +-- abundance present if using abundance-aware tax metagenome
                       +-- taxonomy identifiers match database identifiers
```

Default metagenome profiling should usually start with `dna,k=31,scaled=1000,abund` because sourmash tax docs recommend `k=31` for species-level metagenome work, and abundance is needed for abundance-weighted summaries. Higher `scaled` values are faster and smaller but less sensitive to small overlaps. Branchwater docs suggest `scaled=10000` can be useful when only overlaps larger than roughly 50 kbp matter.

Do not mix taxonomy files and databases from different releases. A plausible-looking `tax metagenome` report can still be wrong if the gather database and taxonomy CSV do not describe the same identifiers.

## Databases and indexes

Prefer sourmash zip collections, standalone manifests, or RocksDB indexes. Avoid pathlists for new branchwater work.

Useful commands:

```bash
# Inspect a collection.
sourmash sig summarize database.zip

# Combine sketches into a zip collection.
sourmash sig cat *.sig.gz -o sketches.zip

# Build a manifest from a path list when needed.
sourmash sig collect pathlist.txt -o manifest.csv -F csv

# Build a sourmash RocksDB index, available in sourmash v4.9+.
sourmash index database.rocksdb sketches.zip -F rocksdb

# Build a branchwater RocksDB index.
sourmash scripts index sketches.zip -o database.rocksdb
```

RocksDB is useful for repeated gather operations against large databases. All sketches in an index must be compatible in k-mer size, scaled value, and molecule type.

## Gather and tax interpretation

`sourmash gather` decomposes one metagenome-like query against one or more reference collections. Its detailed CSV includes fields such as `f_unique_to_query`, `f_unique_weighted`, `average_abund`, `median_abund`, `ksize`, `moltype`, `scaled`, and ANI estimates.

Prefer sum-safe fields like `f_unique_to_query`/`p_query` when aggregating. Raw overlap can double-count because different matches may share query hashes.

In sourmash v4.9.4, `tax metagenome --use-abund` is recommended for abundance-aware results. In sourmash v5, abundance behavior becomes stricter, so missing abundances should be treated as a useful failure rather than silently ignored.

## Search and prefetch adjuncts

Use `search` for quick similarity checks and `prefetch` to reduce very large databases before gather:

```bash
sourmash search query.sig database.zip -o search.csv
sourmash search --containment query.sig database.zip -o containment.csv

sourmash prefetch \
  sample.sig.zip database.zip \
  --threshold-bp 0 \
  --save-matching-hashes matching-hashes.sig \
  --save-matches db-matches.sig \
  -o prefetch.csv

sourmash gather matching-hashes.sig db-matches.sig -o gather.csv
```

Remember that containment wording is easy to invert: `search --containment Q DB` reports query-in-match containment. For metagenome profiling, `gather`/`prefetch` usually use the metagenome sketch as the query.

## Sample similarity QC

NVD also uses the query metagenome sketches for experimental sample-similarity QC. This is not taxonomic profiling and should live in its own subworkflow:

```nextflow
workflow SAMPLE_SIMILARITY_QC
```

Use sourmash names only for processes that directly invoke sourmash:

```nextflow
process SOURMASH_COLLECT_QUERY_SKETCHES
process SOURMASH_COMPARE_QUERY_SKETCHES
```

Use domain names for Python reporting and plotting processes that consume sourmash outputs:

```nextflow
process COMPUTE_SAMPLE_SKETCH_DISTANCES
process COMPUTE_SAMPLE_ORDINATION
process REPORT_POSSIBLE_SAMPLE_MIXUPS
process PLOT_SAMPLE_ORDINATION
```

The basic flow mirrors nt-terroir, but without assuming nt-terroir-specific location/date metadata:

```text
query sketches
  -> sourmash sig cat
  -> sourmash compare (abund and noabund)
  -> pairwise distance TSVs
  -> nearest-neighbor / possible-mix-up TSVs
  -> PCoA coordinates and ordination plots
```

Be conservative about “possible mix-up” language. High sketch similarity between distinct samples is a candidate signal, not proof. Prefer reporting nearest neighbors and thresholded candidates over making hard pass/fail claims until the project has validated thresholds.

## Branchwater caveats

Branchwater is the right default for NVD’s future high-throughput branch, but do not paper over these differences:

- `fastgather` output mostly matches `sourmash gather`, but uses fields like `match_name`, `match_md5`, and `match_filename`, and omits some core sourmash columns such as `potential_false_negative`.
- Branchwater commands prefer zip collections, standalone manifests, and RocksDB indexes.
- RocksDB indexes do not support abundance estimation for containment queries.
- Branchwater docs warn that column order may change. Parse by header.

## Sources

- sourmash classifying signatures: `https://sourmash.readthedocs.io/en/latest/classifying-signatures.html`
- sourmash databases: `https://sourmash.readthedocs.io/en/latest/databases.html`
- sourmash command-line docs: `https://sourmash.readthedocs.io/en/latest/command-line.html`
- branchwater plugin docs: `https://github.com/sourmash-bio/sourmash_plugin_branchwater/blob/main/doc/README.md`
