# NVD2: Identify Human Viruses (and more) in your Metagenomic Samples

[![Release](https://github.com/dhoconno/nvd/actions/workflows/release.yml/badge.svg)](https://github.com/dhoconno/nvd/actions/workflows/release.yml)
[![Latest Release](https://img.shields.io/github/v/release/dhoconno/nvd)](https://github.com/dhoconno/nvd/releases/latest)
[![License](https://img.shields.io/github/license/dhoconno/nvd)](LICENSE)

## Overview

NVD2 is a Nextflow pipeline focused primarily though not exclusively on finding
human viruses in metagenomic samples. It leverages the battle-tested NCBI
toolchain, including
[STAT](https://www.ncbi.nlm.nih.gov/sra/docs/sra-taxonomy-analysis-tool/) and
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), to identify human viruses
while minimizing false positives _and_ false negatives (the STAT+BLAST workflow,
historically called "NVD"). For more general purpose classification, NVD2
implements a [GOTTCHA2](https://github.com/poeli/GOTTCHA2) subworkflow, which
uses a carefully curated database of taxonomically diagnostic reference
sequences spanning well-characterized taxa across the tree of life.

NVD2 was designed from the ground up to handle enormous datasets and performs
particularly well with complex Illumina deep sequencing datasets like those from
wastewater sewersheds. To perform well with these kinds of datasets, it must:

1. Handle highly fragmented genome recovery.
2. Be resilient to wild fluctuations in depth-of-coverage.
3. Resolve ambiguities between closely related organisms with high sequence
   identity.

Many pipelines for classifying mixtures of organisms exist, but none satisfied
these criteria to a satisfactory degree for human viruses--hence, NVD2 was born!

In addition to wastewater sequenced with Illumina, NVD2 also supports datasets
sequenced on Oxford Nanopore instruments and has been tested on Nanopore
libraries generated from air samples.

## Get Started

NVD2 set-up is a multi-phase process, including dependency setup, reference
database setup, sample data setup, and run command construction.

### Dependency setup

#### Option 1: Interactive Setup Script

NVD2 includes an interactive setup script to help configure database paths and
create a configuration file:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash
```

Or download and inspect first:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh -o install.sh
chmod +x install.sh
./install.sh
```

**Prerequisites** (must be installed separately):

- Java 11 or newer
- Nextflow
- Docker, Apptainer/Singularity, or Pixi (for containerized execution)

The script will help you:

- Check that prerequisites are installed
- Detect available execution environments (Docker, Apptainer)
- Configure database paths
- Optionally download reference databases
- Create a configuration file at `~/.nvd2/config/user.config`

See the [Installation Guide](./docs/INSTALLATION.md) for more details.

#### Option 2: Manual Setup

**Required dependencies:**

- Java 11 or newer (required by Nextflow)
- [Nextflow](https://nextflow.io/)
- Container runtime: [Docker](https://www.docker.com/),
  [Apptainer/Singularity](https://apptainer.org/), or [Pixi](https://pixi.sh)
  (for local execution)

With Nextflow and a container runtime installed, the pipeline can run using
containerized dependencies. No additional software installation needed.

**For Conda-based dependencies** (optional): The pipeline includes a
`pyproject.toml` and `pixi.lock` for reproducible Conda environments via
[Pixi](https://pixi.sh/latest/). Note that Pixi is only needed for managing
Conda dependencies - containerized execution (Docker/Apptainer) is generally
preferred and doesn't require Pixi.

```bash
# If using Pixi for local execution
pixi shell --frozen
```

> [!WARNING]\
> The provided Pixi environment will _not_ include NCBI's STAT tool, as it's not
> distributed through Conda or PyPI and must instead be built from source. As
> such, to use the `nvd` subworkflow that implements STAT with two phases of
> BLAST verification, users must use the pre-built Docker or Apptainer
> containers for this project. Alternatively, the bundled
> [Containerfile](./Containerfile) can be built with Docker or Podman to provide
> all dependencies, including STAT as well as those managed by Pixi.

### Reference database setup

The hardest/slowest part of NVD2 setup is getting the necessary reference
datasets organized, though ideally, you'll only need to do it once.

Currently, NVD2 uses three datasets:

1. A BLAST database built from the NCBI core_nt database. This database is used
   to triage putative taxonomic classifications as they make their way through
   the pipeline.
2. An NCBI taxonomic classification database used by STAT.
3. A GOTTCHA2 curated database.

Rather than bundling all three together, we've made each available individually
so that users can choose which subworkflows they'd like to run. If they'd like
to run the STAT+BLAST subworkflow, named nvd after NVD2's predecessor, users
will need to download databases 1 and 2 above. To run GOTTCHA2, they will need
the third.

All three are publicly available via `wget` or `curl` from the O'Connor
Laboratory's LabKey server, like so:

```bash
# download the STAT database
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/stat_db.tar.gz

# download the BLAST database
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/blast_db.tar.gz

# download the GOTTCHA2 database
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/gottcha2.tar.gz
```

(`curl -fSsL` can be substituted for `wget` in the above commands if desired)

> [!IMPORTANT]\
> We strongly recommend users verify their database downloads with the md5
> hashes available in `checksum.txt`, which can be downloaded with
> `wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/checksum.txt`.
> Updates to the reference databases will also be reflected in `CHANGELOG.md`,
> available at the same endpoint as the databases and checksum text file.

Also at that endpoint, if desired, is a pre-built Apptainer image file for use
on HPC cluster or other linux environments:

```bash
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/nvd2.sif
```

All TAR-archived reference databases must be extracted into directories with
`tar xvf` before use with the NVD2 pipeline.

### Sample data setup

With the environment and source code set up, next you'll need to organize the
file paths to your input FASTQ-formatted sequencing read files into a simple,
CSV-formatted samplesheet. It must look like this:

```csv
sample_id,srr,platform,fastq1,fastq2
nanopore_test,,ont,nanopore.fastq.gz,
illumina_test,,illumina,illumina_R1.fastq.gz,illumina_R1.fastq.gz
sra_test,SRR33296246,,
```

Note that this example samplesheet is provided in the repo's [assets](./assets)
directory for convenience.

### Run command setup

With that, you're ready to run the pipeline!

#### Option 1: Using the NVD2 CLI Wrapper

If you have set up `~/.nvd2/config/user.config`, you can use the streamlined CLI
wrapper:

```bash
# Clone the repository
git clone https://github.com/dhoconno/nvd.git
cd nvd

# Install Python dependencies
uv sync
# or for all Conda dependencies
pixi shell --frozen

# Run with minimal flags (prepend with pixi run or uv run)
nvd run --samplesheet samples.csv --experiment-id exp001

# Specify tools and profile
nvd run -s samples.csv -e exp002 --tools gottcha -p docker

# Override database paths if needed
nvd run -s samples.csv -e exp003 \
    --blast-db /custom/blast_db \
    --gottcha2-db /custom/gottcha2.fna

# Dry-run to preview the command
nvd run -s samples.csv -e exp004 --dry-run
```

The CLI wrapper automatically loads configuration from
`~/.nvd2/config/user.config`, auto-detects your execution profile
(Docker/Apptainer/local), and provides a simplified interface.

**Note**: Prepend commands with `pixi run` (for full toolchain access) or
`uv run` (Python-only scripts) when not in an active environment shell.

See the [CLI Wrapper Guide](./docs/cli_wrapper_guide.md) for complete examples
and usage.

#### Option 2: Direct Nextflow Execution

Before you construct your run command, first answer the following questions:

- Are you interested only in human viruses? If so, include `--tools stat_blast`
  (or `nvd`, `stat`, `blast`, `stast` - all equivalent) in your `nextflow run`
  command (more on this below).
- Are you interested in whatever's in your sample and not in human viruses in
  particular? If so, use `--tools gottcha`.
- Are you interested in both? If so, use `--tools stat_blast,gottcha` (or
  `nvd,gottcha` for backward compatibility)
- Do you want the kitchen sink? We use `--tools all` for that.

Beyond the answers to these questions, most of the default run command will
simply be devoted to configuring paths to the required reference database files.
Note that we plan to use presets to simply normal use-case run commands, but for
now, an example run command to use all subworkflows will look like:

```bash
nextflow run dhoconno/nvd \
	--tools all \
	--samplesheet $YOUR_SAMPLESHEET \
	--gottcha2_db $YOUR_REFERENCE_PATH/gottcha2/gottcha_db.species.fna \
	--blast_db $YOUR_REFERENCE_PATH/blast_db \
	--blast_db_prefix core_nt \
	--stat_index $YOUR_REFERENCE_PATH/STAT_db/tree_index.dbs \
	--stat_dbss $YOUR_REFERENCE_PATH/STAT_db/tree_filter.dbss \
	--stat_annotation $YOUR_REFERENCE_PATH/STAT_db/tree_filter.dbss.annotation \
	--human_virus_taxlist $YOUR_REFERENCE_PATH/STAT_db/human_viruses_taxlist.txt \
	--experiment_id github_readme_test
```

(This command assumes you have set `YOUR_SAMPLESHEET` and `YOUR_REFERENCE_PATH`
to the paths to your samplesheet CSV and the parent directory of your extracted
reference databases, respectively; you can replace them with whatever valid path
you have used for these files.)

## Further Documentation

- **[CLI Wrapper Guide](./docs/cli_wrapper_guide.md)** - Comprehensive guide to
  the `nvd` command-line tool (recommended)
- **[Installation Guide](./docs/INSTALLATION.md)** - Detailed installation
  instructions using `install.sh`
- **[Direct Nextflow Examples](./docs/example_commands.md)** - Traditional
  Nextflow command examples
- **[Contributor Guide](./docs/contributor_guide.md)** - Development guidelines
  and best practices

## Grab bag of features

- **Multi-platform sequencing support**: Seamlessly processes both Illumina and
  Oxford Nanopore data with platform-specific optimizations
- **Dual classification engines**: Combines NCBI STAT+BLAST for human virus
  detection with GOTTCHA2 for comprehensive taxonomic profiling
- **Industrial-scale data processing**: Built from the ground up to handle
  massive wastewater datasets with complex read mixtures and variable coverage
  depths
- **Smart contig assembly**: Automatically assembles reads with SPAdes and
  filters contigs for optimal classification accuracy
- **Two-phase BLAST verification**: Uses both megablast and blastn with
  intelligent filtering to minimize false positives
- **Least Common Ancestor (LCA) resolution**: Resolves ambiguous BLAST hits by
  computing taxonomic consensus—either using dominant taxid assignment when one
  organism has strong support (>80% bitscore weight), or calculating the LCA for
  near-tie cases to avoid over-specificity when multiple closely-scoring hits
  disagree at the species level
- **Advanced taxonomic filtering**: Sophisticated lineage-based filtering with
  adjustable stringency for precise organism identification
- **Human read scrubbing**: Built-in capability to remove human sequences for
  privacy-compliant public data sharing
- **Automated data deduplication**: CLUMPIFY workflow removes PCR duplicates and
  optical duplicates to improve analysis accuracy and reduce disk usage for big
  datasets
- **Enterprise data integration**: Native LabKey LIMS integration with WebDAV
  file uploads and structured metadata management
- **Comprehensive quality control**: Read counting, contig metrics, and BLAST
  hit validation throughout the pipeline
- **Flexible workflow orchestration**: Mix-and-match subworkflows (nvd, gottcha,
  clumpify) based on research needs
- **Production-ready deployment**: Docker/Apptainer containerization with Pixi
  environment management for reproducible execution
- **Intelligent error handling**: Robust retry logic and graceful failure modes
  for reliable high-throughput processing
- **SRA integration**: Direct processing of NCBI SRA datasets alongside local
  FASTQ files
- **Real-time validation**: Pre-flight checks for database integrity, API
  connectivity, and experiment ID uniqueness
- **Multi-format output**: Generates taxonomic reports, FASTA sequences, and
  structured CSV files for downstream analysis

## Citation

Coming soon!

## Acknowledgments

- [Marc Johnson](https://medicine.missouri.edu/faculty/marc-johnson-phd) and
  [Shelby O'Connor](https://slo.pathology.wisc.edu), our partners in innovative
  pathogen monitoring from environmental samples.
- Kenneth Katz, NCBI, for developing NCBI STAT, maintaining pre-built databases
  for STAT, and helpful discusssions
- [C. Titus Brown](https://www.vetmed.ucdavis.edu/faculty/c-titus-brown), for
  helpful discussions of using kmer classifiers as part of metagenomic workflows
- Development funded by [Inkfish](https://ink.fish)

## License

See `LICENSE` for more information. NVD2's predecessor used the
[copyleft](https://en.wikipedia.org/wiki/Copyleft) GPLv3 license, which means we
have to as well. This means you’re welcome to use it, share it, and modify it,
**as long as any changes you distribute are also shared under the same
license**.

By contributing to this project, you are thus locked into agreeing that your
code will be released under GPLv3. This means:

- **Share alike** – if you share modified versions of this project, they _must_
  also be under GPLv3.
- **Freedom to use** – anyone can use the code for personal, educational, or
  commercial purposes, as long as they respect the license.
- **Source availability** – if you distribute the software (original or
  modified), you must also make the source code available under GPLv3.
- **Community contributions** – your pull requests and patches automatically
  become part of the GPLv3-licensed project.
- **Commercial use** - You’re free to use NVD2 code in commercial contexts, but
  you can’t combine it with proprietary software without open-sourcing that
  software under GPLv3 too.
- **No Warranty and Liability** - The software is provided as-is, without any
  warranty or liability.
