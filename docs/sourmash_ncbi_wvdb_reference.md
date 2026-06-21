# Building a combined NCBI Virus + WVDB sourmash reference

Alongside its core assembly-then-BLAST workflow, NVD can also provide kmer-based abundance estimates for non-novel taxa in the sample. While this feature is currently experimental, our typical usage of it has been against the 2025.01 release of NCBI Virus, which was pre-sketched by the sourmash maintainers.

NCBI Virus is a good start, but it leaves a fair amount of “dark matter” when checked against a wastewater sample. To address this gap, our collaborator Rose Kantor has been assembling and curating a database called the Wastewater Virus Database (WVDB) using large wastewater samples from the CASPER coalition. The latest version of WVDB is published on Zenodo at <https://zenodo.org/records/20276352>.

To make novel or poorly characterized taxa in WVDB “visible,” we use a procedure where NCBI Virus and its taxonomy are merged with WVDB and its taxonomy. The workflow for performing this merger and generating a new reference artifact is reproduced below. We perform this workflow periodically ourselves and is not part of normal NVD responsibilities.

The runtime artifacts consumed by NVD are only these two files:

```text
ncbi-viruses-2025.01-plus-wvdb-v2.dna.k31.scaled50.sig.zip
ncbi-viruses-2025.01-plus-wvdb-v2.lineages.csv
```

The workflow also produces as reference manifest and build directory, which should be thought of as provenance/debugging artifacts.

## A quick note on organization

Typically, the Nextflow convention is to place scripts in a `bin/` directory when they are executed inside Nextflow processes. While the NCBI Virus + WVDB build uses scripts, these scripts are not executed by Nextflow/NVD. They are thus more accurately described as maintainer/build helpers for creating reference artifacts *before* running NVD. As such, we’ve decided to keep them in a `scripts/` directory instead, which is a common convention for utility scripts that are not part of the core application. See below:

```text
bin/       Nextflow runtime executables
scripts/   maintainer, build, and developer utilities
```

## WVDB lineage policy

At a basic level, the NCBI Virus + WVDB build involves 1) downloading/generating sourmash FracMinHash signatures for NCBI Virus and WVDB, and 2) concatenating them. Of course, this is not the hard part. The hard part is producing a lineage CSV whose identifiers match the WVDB signatures and whose ranks are honest enough for annotating abundance estimates. This is challenging in part because many of the contigs in WVDB are novel/not characterized in other databases in taxonomies. A further challenge is that many contigs do not have a species- or strain-level classification.

Our first stab at projecting WVDB contigs into sourmash’s standard lineage columns is as as follows:

| sourmash column | WVDB source policy |
| --- | --- |
| `ident` | `WVDB|{contig}` |
| `superkingdom` | `Viruses` |
| `phylum` | `Phylum_RdRp`, then `Phylum_gNd` |
| `class` | `Class_RdRp`, then `Class_gNd` |
| `order` | `Order_consensus`, then `Order_RdRp`, then `Order_gNd` |
| `family` | `Family_consensus`, then `Family_RdRp`, then `Family_gNd`, then `Family_ICTV` |
| `genus` | `Genus_RdRp` when present |
| `species` | `Species_RdRp` when present |

Contigs without genus or species assignments are left blank at those ranks. We do not invent species names to satisfy a species-rank output. Blank lower ranks are preferable to false precision.

The FASTA identifiers are normalized with the same `WVDB|` prefix before sketching, so the sourmash signature names and lineage identifiers match.

## Build command

Because the commands used to run the build are intricate and must be run in a specific order with specific directories, we use a `justfile` recipe to provide shorthands. These shorthands require explicit user confirmation before downloading big files:

```bash
just build-sourmash-ncbi-wvdb yes
```

By default, outputs are written under:

```text
build/sourmash/ncbi-virus-wvdb-v2/
```

Use a different build directory when needed:

```bash
just build-sourmash-ncbi-wvdb \
  yes \
  /path/to/reference-build/ncbi-virus-wvdb-v2
```

The workflow is decomposed into smaller recipes so maintainers can rebuild one side without repeating all downloads:

```bash
just fetch-sourmash-ncbi-virus yes
just build-sourmash-wvdb yes
just build-sourmash-ncbi-wvdb yes
```

The top-level `build-sourmash-ncbi-wvdb` recipe uses Just dependencies to run the NCBI and WVDB recipes before combining their outputs. The NCBI recipe downloads the NCBI Virus sourmash sketch and matching lineage CSV first, then validates them with sourmash before any WVDB work. The WVDB recipe downloads WVDB FASTA and annotations, removes derived WVDB outputs from any previous partial build, runs `scripts/prepare_wvdb_sourmash_inputs.py` to write a normalized WVDB FASTA and WVDB lineage CSV, sketches WVDB with `dna,k=31,scaled=50`, writes a sourmash manifest, checks WVDB manifest coverage, and validates the WVDB taxonomy with `sourmash tax prepare`.

The combined recipe removes stale combined outputs, concatenates the NCBI and WVDB sketches, prepares the combined lineage CSV, writes a sourmash manifest, checks WVDB manifest coverage against the combined lineages, and summarizes the combined sketch.

NVD’s justfile contains URLs where it expects to find the prebuild NCBI Virus signature and lineage at UCSD’s data farm. If the farm URLs are temporarily unavailable, use a mirror by overriding the Just variables rather than editing the recipe:

```bash
just \
  --set NCBI_VIRUS_SOURMASH_SIG_URL https://example.org/ncbi-viruses-2025.01.dna.k31.scaled50.sig.zip \
  --set NCBI_VIRUS_SOURMASH_LINEAGES_URL https://example.org/ncbi-viruses-2025.01.lineages.csv \
  wvdb yes
```

## Registering the reference as an NVD preset

As mentioned, our approach is to build this reference periodically and store it in a publicly accessible location (usually out LabKey LIMS server). From there, we use NVD presets to plug it into NVD’s sourmash taxonomic decomposition workflows, which are still experimental and require `params.experimental == true`.

An example YAML-formatted params file for such a preset would look like this:

```yaml
experimental: true

sourmash_ref_url: "https://example.org/nvd/sourmash/ncbi-viruses-2025.01-plus-wvdb-v2.dna.k31.scaled50.sig.zip"
sourmash_lineages_url: "https://example.org/nvd/sourmash/ncbi-viruses-2025.01-plus-wvdb-v2.lineages.csv"
sourmash_ksize: 31
sourmash_scaled: 50
sourmash_threshold_bp: 50
```

Then register and use it:

```bash
nvd preset register \
  sourmash-ncbi-virus-wvdb-v2 \
  --from-file nvd-sourmash-ncbi-virus-wvdb-v2.yaml \
  --description "Experimental sourmash profiling against NCBI Virus 2025.01 plus WVDB v2"

nvd run --preset sourmash-ncbi-virus-wvdb-v2 --params-file run-specific.yaml
```

NVD will download the two sourmash reference files, sketch query metagenomes with abundance tracking, run taxonomic decomposition against the new reference, and summarize results in a variety of formats.

## Caveats

WVDB was assembled from reads in the same broad sequencing program that NVD processes. That is useful for recovering wastewater viral dark matter, but abundance estimates should be described as reference containment/decomposition rather than independent discovery.

Do not mix sketches and lineage CSVs from different builds. The sketch and lineage file are one reference release and should be published, named, and replaced together.
