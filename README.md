# NVD: Identify Low-Abundance Viruses in Big Environmental Metagenomes

[![Release](https://github.com/dholab/nvd/actions/workflows/release.yml/badge.svg)](https://github.com/dholab/nvd/actions/workflows/release.yml)
[![Latest Release](https://img.shields.io/github/v/release/dholab/nvd)](https://github.com/dholab/nvd/releases/latest)
[![License](https://img.shields.io/github/license/dholab/nvd)](LICENSE)

## Overview

NVD is a Nextflow pipeline focused on finding human-infecting viruses in environmental metagenome samples. It leverages a fast _in silico_ enrichment approach together with de novo assembly and BLAST verification to maximize recall and precision. It supports Illumina short-read as well as Nanopore long-read inputs and offers a configurable deduplication approach tuned for shotgun sequencing datasets.

NVD was designed from the ground up to handle enormous datasets and performs particularly well with complex Illumina deep sequencing datasets like those from wastewater sewersheds. To perform well with these kinds of datasets, it must:

1. Handle highly fragmented genome recovery.
2. Be resilient to wild fluctuations in depth-of-coverage.
3. Resolve ambiguities between closely related organisms with high sequence identity.

Many pipelines for classifying mixtures of organisms exist, but none satisfied these criteria to a satisfactory degree for human viruses--hence, NVD was born!

This pipeline in its 3rd major version, which brings with it a helpful CLI and pipeline control system, vast performance improvements over its predecessor, and a tighter focus on BLAST-searching viral hits. This means better and faster results for viruses, but worse results for other taxa like bacteria. We recommend users try [jhuapl-bio/taxtriage](https://github.com/jhuapl-bio/taxtriage) or [nf-core/mag](https://nf-co.re/mag/5.4.2) if they're interested in more than viruses.

## Get Started

NVD set-up has a few phases, including dependency setup, reference database setup, sample data setup, and run command construction. These phases can be handled manually, but we recommend users use our installer script to get started.

### Dependency setup

#### Option 1: Interactive Installer Script

NVD includes an interactive install script to help configure database paths and establish user settings:

```bash
curl -fsSL https://raw.githubusercontent.com/dholab/nvd/main/install.sh | bash
```

Or download and inspect first:

```bash
curl -fsSL https://raw.githubusercontent.com/dholab/nvd/main/install.sh -o install.sh
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
- Create a configuration file at `~/.nvd/user.config`

See the [Installation Guide](./docs/INSTALLATION.md) for more details.

#### Option 2: Manual Setup

**Required dependencies:**

- Java 11 or newer (required by Nextflow)
- [Nextflow](https://nextflow.io/)
- Container runtime: [Docker](https://www.docker.com/),
  [Apptainer/Singularity](https://apptainer.org/), or [Pixi](https://pixi.sh)
  (for local execution)

With Nextflow and a container runtime installed, the pipeline can run using containerized dependencies. No additional software installation needed.

**For Conda-based dependencies** (optional): The pipeline includes a `pyproject.toml` and `pixi.lock` for reproducible Conda environments via [Pixi](https://pixi.sh/latest/). Note that Pixi is only needed for managing Conda dependencies - containerized execution (Docker/Apptainer) is generally preferred and doesn't require Pixi.

```bash
# If using Pixi for local execution
pixi shell --frozen
```

Note also that new versions of our container image are automatically built and pushed to [docker hub](https://hub.docker.com/r/nrminor/nvd/tags) on each NVD release.

### Reference database setup

NVD is predominantly designed to be run offline, which means users must download the required reference files ahead of time. Since version 3.0, the list of required files is significantly decreased; users just need our [Deacon](https://github.com/bede/deacon) index of human infecting viruses, which is less than half a gigabyte, and the BLAST core_nt database, which is about 230 gigabytes. We provide both via our LabKey server as follows.

> [!IMPORTANT]\
> Prefer the installer for reference setup. The installer can download, verify, arrange, and place the reference artifacts for you, including extracting the BLAST database archive and checking that your system has enough available space for it. Manual download is still supported, but it is easier to get paths, checksums, or extraction layout wrong by hand.

Both of these reference datasets are publicly available via `wget` or `curl` from the O'Connor Laboratory's LabKey server, like so:

```bash
# download the BLAST database archive
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v3.0.0/v3.0.0/blast_db_v3_0.tar.gz

# download the deacon human-infecting virus index
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v3.0.0/v3.0.0/human_infecting_viruses.k31w1.idx
```

(`curl -fSsL` can be substituted for `wget` in the above commands if desired)

If you download manually, verify the files against the checksum manifest:

```bash
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v3.0.0/v3.0.0/checksums_v3_0.txt
```

Also at that endpoint, if desired, is a pre-built Apptainer image file for use on HPC cluster or other linux environments:

```bash
wget https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v3.0.0/v3.0.0/nvd-v3.0.0.sif
```

The installer extracts the BLAST tarball for you. If you download manually, extract `blast_db_v3_0.tar.gz` before use and point `blast_db` at the extracted directory. The deacon index is already a single `.idx` file and does not need extraction.

#### The virus enrichment index

To reduce sequence search space before assembly and BLAST verification, NVD uses a the aforementioned minimizer index representing human-infecting viruses. This index is derived from NCBI's STAT dense tree-of-life minimizer database, filtered to NVD's curated human-virus taxid list, and converted into deacon's `.idx` format. It is used as a fast screening step, not as the final taxonomic classifier.

You can pass it through a params file as `virus_index`, or with the CLI:

```bash
nvd run --params-file run.yaml --virus-index /path/to/human_infecting_viruses.k31w1.idx
```

For the rationale, source STAT database, curated taxid list, conversion script, build parameters, and reproducibility notes, see the [Virus Enrichment Index Guide](./docs/virus_enrichment_index.md).

### NVD setup management

After the initial installer has cloned NVD and installed the Pixi environment, you can re-run the setup portion directly with:

```bash
nvd setup
```

`nvd setup` handles the following:

- detects whether you are on an O'Connor Lab [CHTC](https://chtc.cs.wisc.edu/) access point or a generic system
- detects your shell and can install the `nvd setup shell-hook` line into your shell startup file
- writes or refreshes `~/.nvd/setup.conf`, including the stable repo path used by the wrapper
- on CHTC, records the default profile, shared taxonomy location, and shared preset store location
- on CHTC, creates or refreshes `~/.nvd/user.config` from the CHTC template if you approve overwriting it
- installs or refreshes the lightweight `~/.local/bin/nvd` wrapper script if you approve overwriting it
- generates shell completions
- on CHTC, offers to download or re-download the release SIF into `~/.nvd`

Those conveniences are intentionally separate from run-specific configuration. BLAST database paths, the deacon virus index path, samplesheets, output locations, preprocessing choices, and LabKey settings should still be provided through params files, presets, explicit CLI flags, or direct Nextflow config when you really need it.

Useful setup variants include:

```bash
nvd setup --force              # overwrite existing setup-managed files without repeated prompts
nvd setup --skip-container     # skip the CHTC SIF download prompt
nvd setup --skip-shell-hook    # leave shell startup files untouched
nvd setup --config-dir ~/.nvd  # choose the NVD config directory explicitly
```

### Sample data setup

With the environment and source code set up, next you'll need to organize your input reads into a simple CSV samplesheet. It must look like this:

```csv
sample_id,srr,platform,fastq1,fastq2
nanopore_test,,ont,nanopore.fastq.gz,
illumina_test,,illumina,illumina_R1.fastq.gz,illumina_R1.fastq.gz
sra_test,SRR33296246,,
```

Note that this example samplesheet is provided in the repo's [assets](./assets) directory for convenience.

Again, while manual setup is supported, it's also more error-prone. We instead recommend users take advantage of the `nvd samplesheet` CLI. It can scan a directory of FASTQ files, infer paired-end samples by filename, preview the generated table, write the CSV, and validate the result before you run the pipeline.

For local FASTQ files:

```bash
nvd samplesheet generate --from-dir ./fastqs --platform illumina --output samplesheet.csv
```

For Nanopore/ONT files, use `--platform ont`:

```bash
nvd samplesheet generate --from-dir ./nanopore-fastqs --platform ont --output samplesheet.csv
```

For public SRA inputs, put one accession per line in a text file:

```text
SRR33296246
SRR33296247
```

Then generate a samplesheet from that accession list:

```bash
nvd samplesheet generate --from-sra accessions.txt --platform illumina --output samplesheet.csv
```

If you want to inspect what NVD would write before touching the filesystem, use dry-run mode:

```bash
nvd samplesheet generate --from-dir ./fastqs --platform illumina --dry-run
```

You can also validate any existing samplesheet explicitly:

```bash
nvd samplesheet validate samplesheet.csv
```

Like other NVD subcommands, the samplesheet commands also have shorter aliases for command line wizards:

```bash
nvd samplesheet gen -d ./fastqs -p illumina -o samplesheet.csv
nvd samplesheet val samplesheet.csv
nvd ss gen -d ./fastqs -p illumina -o samplesheet.csv
nvd ss val samplesheet.csv
```

### Running NVD

If you have reference databases downloaded and a samplesheet generated, you're ready to run NVD!

#### Option 1: Using `nvd run`

The impetus for the `nvd` command line interface was to make complex Nextflow runs easier to launch, validate, resume, and repeat. The `nvd run` subcommand is the recommended entry point for normal use. It accepts the same kinds of runtime values you would otherwise pass to Nextflow, but it also handles params-file merging, preset lookup, config discovery, CHTC setup integration, taxonomy environment handling, and friendlier dry-run behavior.

For reproducibility, we recommend putting run settings in a YAML or JSON params file rather than building one long shell command. The generated params files include schema references, so editors like VS Code and Neovim can offer autocomplete and inline validation when YAML or JSON language support is installed.

Start by generating a params file and editing it for your environment:

```bash
nvd params init run.yaml
```

At minimum, a first run usually needs a samplesheet, an experiment ID, BLAST database settings, a deacon enrichment index, and a taxonomy directory. Those can all live in `run.yaml`; paths should point to locations visible to the machine or cluster worker that will execute the pipeline.

Before launching, validate the params file:

```bash
nvd params check run.yaml
```

If you are preparing the file on a laptop but the paths only exist on CHTC or another cluster, skip local path checks:

```bash
nvd params check run.yaml --no-check-paths
```

Then launch the run:

```bash
nvd run --params-file run.yaml --profile docker
```

You can preview the generated Nextflow command without starting the pipeline:

```bash
nvd run --params-file run.yaml --profile docker --dry-run
```

If you'd rather spell out your params in the command line, that's supported too:

```bash
nvd run \
  --samplesheet samplesheet.csv \
  --experiment-id exp001 \
  --blast-db /path/to/blast_db \
  --blast-db-prefix core_nt \
  --virus-index /path/to/human_infecting_viruses.k31w1.idx \
  --taxonomy-dir /path/to/taxdump \
  --profile docker
```

If your group has a shared preset, use it for common reference paths and defaults, then keep run-specific values in your params file:

```bash
nvd run --preset chtc-defaults --params-file run.yaml
```

CLI flags override both presets and params files, which is useful for small one-off changes:

```bash
nvd run --preset chtc-defaults --params-file run.yaml --results ./results/exp001
```

Resume with the same inputs by adding `--resume`, or ask NVD to resume the last cached command:

```bash
nvd run --params-file run.yaml --resume
nvd resume
```

For more on authoring params, managing presets, samplesheet helpers, secrets, and taxonomy setup, see the [NVD CLI Guide](./docs/nvd_cli_guide.md).

**Note**: Prepend commands with `pixi run` when not in an active environment shell, for example `pixi run nvd run --params-file run.yaml`.

#### Option 2: Direct Nextflow Execution

Direct Nextflow execution is still available, but we recommend treating it as a lower-level fallback for debugging workflow behavior. The `nvd run` CLI handles parameter merging, preset lookup, config discovery, CHTC setup integration, taxonomy environment handling, and friendlier validation before it launches Nextflow.

If you do need to run Nextflow directly, pass the same explicit runtime values that `nvd run` would pass for you. In v3 there is no `--tools` selector as in v2; the main workflow runs the current NVD pipeline, and optional behavior is controlled by params.

```bash
nextflow run . \
  -profile docker \
  --samplesheet $YOUR_SAMPLESHEET \
  --experiment_id github_readme_test \
  --blast_db $YOUR_REFERENCE_PATH/blast_db \
  --blast_db_prefix core_nt \
  --virus_index $YOUR_REFERENCE_PATH/human_infecting_viruses.k31w1.idx \
  --taxonomy_dir $YOUR_TAXDUMP_PATH
```

This command assumes `YOUR_SAMPLESHEET` points to your samplesheet CSV, `YOUR_REFERENCE_PATH` points to the directory where the BLAST database directory and deacon index live, and `YOUR_TAXDUMP_PATH` points to an NCBI taxdump directory visible to the machine or worker running the pipeline.

If you are using a generated params file, direct Nextflow can consume it too:

```bash
nextflow run . -profile docker -params-file run.yaml
```

That said, prefer this equivalent CLI form unless you specifically need to bypass the wrapper:

```bash
nvd run --params-file run.yaml --profile docker
```

## Important Notes

### Network Requirements

NVD is designed so heavyweight references are prepared before pipeline execution. The BLAST database, deacon virus index, and NCBI taxonomy dump should be available before launching a serious run, especially on CHTC or other distributed systems where worker jobs may not have the same filesystem view or outbound network access as the login node.

Use `nvd taxonomy ensure` or the installer/setup flow to prepare taxonomy, and use `nvd taxonomy status` to inspect it. On CHTC, `nvd setup` records the shared taxonomy location in `NVD_TAXONOMY_DB`, and `nvd run` passes that through to the pipeline.

### Reference paths are explicit

The installer can download and arrange reference artifacts, but it does not silently persist BLAST or deacon paths as hidden runtime state. Put those paths in a params file, preset, CLI flags, or explicit Nextflow config.

### Taxonomy data

NVD uses NCBI taxonomy files for BLAST annotation and LCA resolution. The taxonomy directory should contain the standard taxdump files:

- `nodes.dmp`, `names.dmp`, `merged.dmp` - Raw NCBI taxonomy files
- `taxonomy.sqlite` - Indexed database for fast lookups

For local CLI taxonomy commands, NVD can manage a local taxonomy directory. For pipeline runs, especially distributed runs, the taxonomy directory should be explicit and visible to worker jobs.

Useful commands:

```bash
nvd taxonomy ensure --taxonomy-dir /path/to/taxdump
nvd taxonomy status --taxonomy-dir /path/to/taxdump
```

You can also set:

```bash
export NVD_TAXONOMY_DB=/shared/path/taxdump
```

When using `nvd run`, this environment variable is materialized into the pipeline's `--taxonomy_dir` parameter if you do not pass `--taxonomy-dir` directly.

### Offline Mode

For offline or restricted-network environments, pre-populate the BLAST database, deacon virus index, and taxonomy directory before launching the pipeline. To prevent taxonomy helpers from attempting network refreshes, set:

```bash
export NVD_TAXONOMY_OFFLINE=1
```

Offline mode requires the taxonomy directory to already contain the expected NCBI taxdump files.

## Further Documentation

- **[NVD CLI Guide](./docs/nvd_cli_guide.md)** - Comprehensive guide to
  the `nvd` command-line tool, params files, and presets (recommended)
- **[Installation Guide](./docs/INSTALLATION.md)** - Detailed installation
  instructions using `install.sh`
- **[Direct Nextflow Examples](./docs/example_commands.md)** - Traditional
  Nextflow command examples
- **[Virus Enrichment Index Guide](./docs/virus_enrichment_index.md)** - Notes on
  the deacon human-infecting virus enrichment index
- **[Contributor Guide](./docs/contributor_guide.md)** - Development guidelines
  and best practices

## Grab bag of features

- **Multi-platform sequencing support**: Seamlessly processes both Illumina and Oxford Nanopore data with platform-specific optimizations
- **Smart contig assembly**: Automatically assembles reads with SPAdes and filters contigs for optimal classification accuracy
- **Two-phase BLAST verification**: Uses both megablast and blastn with intelligent filtering to minimize false positives
- **Least Common Ancestor (LCA) resolution**: Resolves ambiguous BLAST hits by computing taxonomic consensus—either using dominant taxid assignment when one organism has strong support (>80% bitscore weight), or calculating the LCA for near-tie cases to avoid over-specificity when multiple closely-scoring hits disagree at the species level
- **Advanced taxonomic filtering**: Sophisticated lineage-based filtering with adjustable stringency for precise organism identification
- **Human read scrubbing**: Built-in capability to remove human sequences for privacy-compliant public data sharing
- **LIMS data integration**: Native LabKey LIMS integration with WebDAV file uploads and structured metadata management
- **Comprehensive quality control**: Read counting, contig metrics, and BLAST hit validation throughout the pipeline
- **Flexible workflow orchestration**: Mix-and-match subworkflows (nvd, gottcha, clumpify) based on research needs
- **Production-ready deployment**: Docker/Apptainer containerization with Pixi environment management for reproducible execution
- **Intelligent error handling**: Robust retry logic and graceful failure modes for reliable high-throughput processing
- **SRA integration**: Direct processing of NCBI SRA datasets alongside local FASTQ files
- **Real-time validation**: Pre-flight checks for database integrity, API connectivity, and experiment ID uniqueness
- **Multi-format output**: Generates taxonomic reports, FASTA sequences, and structured CSV files for downstream analysis

## Citation

Coming soon!

## Acknowledgments

- [Marc Johnson](https://medicine.missouri.edu/faculty/marc-johnson-phd) and
  [Shelby O'Connor](https://slo.pathology.wisc.edu), our partners in innovative pathogen monitoring from environmental samples.
- Kenneth Katz, NCBI, for developing NCBI STAT, maintaining pre-built databases
  for STAT, and helpful discusssions
- [C. Titus Brown](https://www.vetmed.ucdavis.edu/faculty/c-titus-brown), for helpful discussions of using kmer classifiers as part of metagenomic workflows
- Development funded by [Inkfish](https://ink.fish)

## License

See `LICENSE` for more information. NVD2's predecessor used the [copyleft](https://en.wikipedia.org/wiki/Copyleft) GPLv3 license, which means we have to as well. This means you’re welcome to use it, share it, and modify it, **as long as any changes you distribute are also shared under the same
license**.

By contributing to this project, you are thus locked into agreeing that your code will be released under GPLv3. This means:

- **Share alike** – if you share modified versions of this project, they _must_ also be under GPLv3.
- **Freedom to use** – anyone can use the code for personal, educational, or commercial purposes, as long as they respect the license.
- **Source availability** – if you distribute the software (original or modified), you must also make the source code available under GPLv3.
- **Community contributions** – your pull requests and patches automatically become part of the GPLv3-licensed project.
- **Commercial use** - You’re free to use NVD2 code in commercial contexts, but you can’t combine it with proprietary software without open-sourcing that software under GPLv3 too.
- **No Warranty and Liability** - The software is provided as-is, without any warranty or liability.
