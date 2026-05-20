# NVD Installation Guide

This guide describes the supported v3 installation path. The recommended interface is the `nvd` CLI: use it to set up the local wrapper, author parameter files, validate configuration, and launch the pipeline. Direct `nextflow run` execution remains available for debugging and advanced fallback use, but it is not the primary user workflow.

## Quick start

Run the installer from a terminal:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash
```

To review the script first:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh -o install.sh
less install.sh
bash install.sh
```

To pass installer flags while using the `curl | bash` form, use `bash -s --`:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash -s -- --dry-run
```

## Prerequisites

NVD expects these tools to already be available:

- Java 11 or newer, required by Nextflow
- Nextflow
- Pixi, used by the checked-out NVD CLI environment
- one container/runtime option: Docker, Podman, Apptainer, or Singularity
- Git and tar

On CHTC, Java, Nextflow, Pixi, and Apptainer are typically provided on the access point. On a workstation, Docker or Podman is usually the easiest runtime.

## What the installer does

The installer is intentionally a bootstrapper, not a hidden state manager. It checks prerequisites, clones or updates the NVD repository under `~/.nvd`, installs the locked Pixi environment, and runs `nvd setup`.

On CHTC, setup also writes `~/.nvd/setup.conf`, configures the default profile, records the shared taxonomy location, records the shared preset store location, and can help place the release SIF under `~/.nvd`.

The optional reference wizard can download the BLAST database archive and deacon vertebrate-virus index. If you choose a custom reference path, that path is printed for you to put in a params file, preset, CLI flags, or explicit Nextflow config. The installer does not persist BLAST/deacon reference paths in another hidden state file.

## Installer modes

### Interactive mode

Interactive mode is the default:

```bash
bash install.sh
```

It asks before updating an existing checkout, runs `nvd setup`, and asks whether to enter the reference download wizard.

### Dry-run mode

Dry-run mode walks the installer and reference wizard prompts but avoids side effects such as cloning, pulling updates, downloading references, extracting archives, or removing downloaded archives.

```bash
bash install.sh --dry-run
```

or, through the piped installer form:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash -s -- --dry-run
```

Use this to verify the setup path and see which reference actions would be taken without downloading large artifacts.

### Non-interactive mode

Non-interactive mode checks prerequisites and exits. It is meant for CI and readiness checks, not for configuring a user environment.

```bash
bash install.sh --non-interactive
```

## CHTC setup

On an O'Connor Lab CHTC access point, `nvd setup` detects the host and uses CHTC defaults. The generated `~/.nvd/setup.conf` should include values like:

```bash
NVD_REPO=/home/you/.nvd/latest
NVD_CONFIG_DIR=/home/you/.nvd
NVD_TAXONOMY_DB=/staging/groups/oconnor_group/nvd/taxdump
NVD_PRESET_STORE=/staging/groups/oconnor_group/nvd/presets.sqlite
NVD_DEFAULT_PROFILE=chtc_htc
```

The taxonomy and preset-store paths are setup/runtime integration points. BLAST and deacon reference paths are still runtime parameters and should be supplied through params files, presets, or CLI flags.

The shell hook installed by setup loads this configuration in new shells. After setup, either start a new shell or run:

```bash
source ~/.bashrc
```

## Reference artifacts

NVD v3 uses these main reference artifacts:

- a BLAST database directory, usually from the v3 BLAST database archive
- a BLAST database prefix, currently `core_nt` for the provided database
- a deacon vertebrate-virus index, usually `human_infecting_viruses.k31w1.idx`
- an NCBI taxonomy taxdump directory for taxonomy resolution and LCA annotation

The installer reference wizard can download the BLAST archive and deacon index from the v3 release endpoint and verify checksums when the manifest is available. If the checksum manifest or a specific checksum entry is unavailable, the installer warns and continues.

After downloading references, put the printed paths in a params file, preset, or command-line flags. A params file is the recommended route:

```yaml
# yaml-language-server: $schema=https://raw.githubusercontent.com/dhoconno/nvd/main/schemas/nvd-params.latest.schema.json
samplesheet: /path/to/samplesheet.csv
experiment_id: experiment-001
results: /path/to/results
blast_db: /path/to/references/blast_db
blast_db_prefix: core_nt
virus_index: /path/to/references/human_infecting_viruses.k31w1.idx
taxonomy_dir: /path/to/taxdump
```

For CHTC users, `taxonomy_dir` can usually be omitted from the params file if `NVD_TAXONOMY_DB` was written by setup and the `nvd run` wrapper is used.

## Parameter-file workflow

Generate a schema-backed params file:

```bash
nvd params init run.yaml
```

Open it in an editor with YAML or JSON language-server support. VS Code users should install a YAML extension so the schema comment enables autocomplete and inline validation.

Validate before launching:

```bash
nvd params check run.yaml
nvd params check run.yaml --no-check-paths   # useful when paths exist only on CHTC
```

Run with the params file:

```bash
nvd run --params-file run.yaml
```

Values are resolved in this order:

```text
CLI flags > params file > preset > pipeline defaults
```

## Running NVD

The preferred launch path is the CLI:

```bash
nvd run \
  --samplesheet samples.csv \
  --experiment-id experiment-001 \
  --blast-db /path/to/blast_db \
  --blast-db-prefix core_nt \
  --virus-index /path/to/human_infecting_viruses.k31w1.idx \
  --profile docker
```

Or use a params file:

```bash
nvd run --params-file run.yaml --profile docker
```

On CHTC, setup may provide a default profile, so the profile flag can often be omitted:

```bash
nvd run --params-file run.yaml
```

Direct Nextflow execution is available when debugging the workflow layer directly:

```bash
nextflow run . \
  -profile docker \
  --samplesheet samples.csv \
  --experiment_id experiment-001 \
  --blast_db /path/to/blast_db \
  --blast_db_prefix core_nt \
  --virus_index /path/to/human_infecting_viruses.k31w1.idx \
  --taxonomy_dir /path/to/taxdump
```

Prefer `nvd run` for normal use because it applies the CLI params model, preset handling, config discovery, taxonomy environment handling, and friendlier validation before launching Nextflow.

## Updating an existing install

Re-run the installer. If `~/.nvd/latest` already exists, the installer asks before pulling updates. Declining that prompt must leave the checkout untouched.

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash
```

After updates, setup is run again so shell hooks and `setup.conf` stay current.

## Troubleshooting

### Installer prompts do not appear interactive

Run the installer from a real terminal. The `curl | bash` path uses the script pipe as standard input, so the installer reads prompts from the controlling terminal. If no terminal is available, it should fail closed instead of choosing defaults.

For a no-side-effect check of the prompt flow, use dry-run:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash -s -- --dry-run
```

### Docker is installed but not running

Start Docker Desktop on macOS, or start the Docker daemon on Linux:

```bash
sudo systemctl start docker
```

Alternatively use Apptainer, Singularity, Podman, or the CHTC profile where appropriate.

### Insufficient disk space

Choose a reference path on a filesystem with enough space, download fewer artifacts, or skip the reference wizard and configure paths later. The BLAST database archive is large; do not place reference databases on short-lived scratch storage.

### Java or Nextflow is missing

Install Java 11 or newer and Nextflow, or load the appropriate modules on an HPC system:

```bash
module load java/17
module load nextflow
module load apptainer
```

Then check readiness:

```bash
bash install.sh --non-interactive
```

### Params validate locally but paths fail on CHTC

Use `--no-check-paths` when validating params on a machine that cannot see the cluster filesystem:

```bash
nvd params check run.yaml --no-check-paths
```

The cluster run still needs paths that are visible to worker jobs.

## Security and permissions

The installer does not need root. It writes under `~/.nvd`, shell startup files when you approve shell-hook installation, and reference locations you choose. It should ask before mutating an existing checkout and should fail rather than silently defaulting when an interactive answer cannot be read.

Reference downloads are verified when a checksum manifest and matching entry are available. If verification cannot be performed because the manifest or entry is missing, the installer warns and continues.

## Getting help

Useful commands:

```bash
nvd --help
nvd setup --help
nvd params --help
nvd params check run.yaml --no-check-paths
nvd taxonomy status --taxonomy-dir /path/to/taxdump
```

Project documentation and issues are available at <https://github.com/dhoconno/nvd>. For day-to-day usage after installation, see the [NVD CLI Guide](./nvd_cli_guide.md).

## License

NVD is licensed under GPLv3. See the repository `LICENSE` file for details.
