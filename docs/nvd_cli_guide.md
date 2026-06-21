# NVD CLI Guide

The `nvd` command is the recommended interface for NVD v3. It wraps the Nextflow pipeline with parameter validation, params-file authoring, presets, samplesheet helpers, setup integration, taxonomy preparation, and friendlier run/resume behavior. Direct `nextflow run` remains available, but users will find that experience comes with more friction.

## Mental model

NVD has a few layers:

```text
install.sh
  bootstraps ~/.nvd/latest and runs nvd setup

nvd setup
  wires the CLI into your shell and environment

nvd params / nvd preset / nvd samplesheet / nvd taxonomy
  prepare explicit run inputs

nvd run
  merges params, validates inputs, and launches Nextflow
```

The CLI is not meant to hide the pipeline from you. It is meant to keep important choices explicit while removing repetitive error-prone plumbing.

## Typical workflow

```bash
nvd setup # run on install; does not need to be run for each nvd run

nvd samplesheet generate --from-dir ./fastqs --platform illumina --sanitize --output samplesheet.csv
nvd params init run.yaml # then edit run.yaml in an editor with YAML language support
nvd params check run.yaml
nvd run --params-file run.yaml
```

For CHTC, `nvd setup` can record the shared taxonomy directory, shared preset store, and default profile in `~/.nvd/setup.conf`. BLAST and deacon reference paths remain run parameters; put them in a params file, preset, or explicit CLI flags.

## Schema-backed params files

Params files are the best way to author run configuration. They are easier to review than long commands and work well with editor schema support.

Generate a YAML template:

```bash
nvd params init run.yaml
```

Generate JSON instead:

```bash
nvd params init run.json
```

The YAML template starts with a `yaml-language-server` schema comment, and the JSON template includes `$schema`. Install YAML or JSON language-server support in your editor, e.g. VS Code, so you get autocomplete, inline validation, and typo detection while editing params. We cannot recommend this part of the `nvd` setup enough!

Validate before launching:

```bash
nvd params check run.yaml
```

If you are validating on a laptop but the paths only exist on a cluster, skip local path checks:

```bash
nvd params check run.yaml --no-check-paths
```

Merge layered params files when you want a base config plus a run-specific override:

```bash
nvd params merge base.yaml run-overrides.yaml --output merged.yaml
```

## Presets

Presets are named collections of params. They are a core quality-of-life feature for shared research workflows: one person can register a known-to-be-good combination of reference paths and analysis defaults, and other runs/users can reuse it without pasting large config blocks. Critically, unlike editing nextflow params files, presets are always validated; if you try to register an invalid `nvd` preset, it will be rejected.

On CHTC, `nvd setup` records a shared preset store in `NVD_PRESET_STORE`, currently intended for O'Connor Lab shared presets. Outside CHTC, presets default to a local SQLite file under your NVD config directory unless you set `NVD_PRESET_STORE` yourself (we recommend you do!).

Register a preset from a schema-backed params file:

```bash
nvd preset register chtc-defaults --from-file chtc-defaults.yaml --description "Shared CHTC reference paths and defaults"
```

List and inspect presets:

```bash
nvd preset list
nvd preset show chtc-defaults # pretty-prints the preset's validated params
```

Use a preset when checking a params file:

```bash
nvd params check run.yaml --preset chtc-defaults --no-check-paths
```

Use a preset when running:

```bash
nvd run --preset chtc-defaults --params-file run.yaml
```

CLI flags override both the params file and the preset, so this is safe for one-off changes:

```bash
nvd run --preset chtc-defaults --params-file run.yaml --results /staging/groups/oconnor_group/my-run/results
```

Compare two presets before launching:

```bash
nvd preset diff chtc-defaults stringent-blast
```

Merge presets when you want a new named combination:

```bash
nvd preset merge chtc-defaults stringent-blast --name chtc-stringent --description "CHTC defaults plus stricter BLAST settings"
```

Export/import presets for review or migration:

```bash
nvd preset export chtc-defaults --output chtc-defaults.yaml
nvd preset import chtc-defaults.yaml --name chtc-defaults
```

Delete a preset when it is no longer valid:

```bash
nvd preset delete old-defaults
```

Aliases are available for common interactive use:

```bash
nvd preset ls
nvd preset reg quick-test --from-file quick-test.yaml
nvd preset rm quick-test
```

## Samplesheets

Generate a samplesheet from a directory of FASTQ files:

```bash
nvd samplesheet generate --from-dir ./fastqs --platform illumina --output samplesheet.csv
```

For Illumina/CASAVA filenames, add `--sanitize` when generated sample IDs should drop the sample-number and lane suffixes. For example, `patient-001_S7_L003_R1_001.fastq.gz` and `patient-001_S7_L003_R2_001.fastq.gz` become sample ID `patient-001` instead of `patient-001_S7_L003`:

```bash
nvd samplesheet generate --from-dir ./fastqs --platform illumina --sanitize --output samplesheet.csv
```

`--sanitize` only changes generated sample IDs. It does not group multiple lanes by itself; multi-lane Illumina samples still produce one row per discovered read pair.

For Illumina/CASAVA directories where each biological sample may have reads split across multiple lanes, add `--group-lanes`. This writes one row per CASAVA sample prefix and leaves the FASTQ files untouched. NVD expands the grouped lane declarations when the run starts, so generation fails early if a lane is missing its R1 or R2 mate or if filenames do not follow the expected CASAVA pattern.

```bash
nvd samplesheet generate --from-dir ./fastqs --platform illumina --group-lanes --sanitize --output samplesheet.csv
```

The generated CSV uses the lower-level `fastq1_glob` and `fastq2_glob` columns only when grouped lanes are requested; ordinary exact-path samplesheets can keep the shorter five-column shape. Grouped-lane glob patterns should be considered “live”: NVD resolves them during validation and at run startup. If new matching lane files appear later on a `nextflow run ... -resume`, say, those files will become part of the run. See `assets/grouped_lane_samplesheet.csv` for an illustrative expanded form. Its `/data/run42` paths are placeholders and are expected to fail filesystem-sensitive validation unless you replace them with paths that exist in your run environment.

Samplesheet validation is includes filesystem checks. Local FASTQ paths and glob patterns must thus be absolute and accessible when validation runs. NVD preserves user-facing absolute symlink paths rather than rewriting them through `realpath`; on clusters, use the stable absolute namespace that nodes can also see.

For Nanopore/ONT reads:

```bash
nvd samplesheet generate --from-dir ./nanopore-fastqs --platform ont --output samplesheet.csv
```

Generate from SRA accessions, one accession per line:

```bash
nvd samplesheet generate --from-sra accessions.txt --platform illumina --output samplesheet.csv
```

Preview without writing:

```bash
nvd samplesheet generate --from-dir ./fastqs --platform illumina --dry-run
```

Validate an existing samplesheet:

```bash
nvd samplesheet validate samplesheet.csv
```

Short aliases exist for interactive use:

```bash
nvd samplesheet gen -d ./fastqs -p illumina -o samplesheet.csv
nvd samplesheet val samplesheet.csv
nvd ss gen -d ./fastqs -p illumina -o samplesheet.csv
nvd ss val samplesheet.csv
```

## Taxonomy

NVD uses NCBI taxonomy for BLAST annotation and LCA resolution. For local work, you can ask the CLI to prepare a taxonomy directory:

```bash
nvd taxonomy ensure --taxonomy-dir /path/to/taxdump
```

Existing taxonomy data is reused even when old, while missing required files still trigger download/build. This is the safest behavior for shared HPC references because worker jobs do not mutate a shared taxonomy directory merely because `nodes.dmp` is older than a freshness window.


Inspect the current state:

```bash
nvd taxonomy status --taxonomy-dir /path/to/taxdump
```

On CHTC, setup records the shared taxonomy location in `NVD_TAXONOMY_DB`, and `nvd run` materializes that environment variable into the pipeline's `--taxonomy_dir` parameter when you do not pass `--taxonomy-dir` explicitly.

For offline environments, pre-populate taxonomy and set:

```bash
export NVD_TAXONOMY_OFFLINE=1
```

## Setup and configuration

Run setup after installation, or re-run it later when you change machines, shells, or CHTC defaults:

```bash
nvd setup
```

Useful setup options:

```bash
nvd setup --force
nvd setup --skip-container
nvd setup --skip-shell-hook
nvd setup --config-dir ~/.nvd
```

Setup manages shell integration, the lightweight `~/.local/bin/nvd` wrapper, completions, `~/.nvd/setup.conf`, and CHTC-specific defaults. It does not silently persist BLAST/deacon reference paths.

Inspect configuration discovery:

```bash
nvd config show
nvd config path
```

Open the active config in your editor:

```bash
nvd config edit
```

Config aliases:

```bash
nvd cfg s
nvd cfg p
nvd cfg e
```

## Running the pipeline

Run with a params file:

```bash
nvd run --params-file run.yaml
```

Run with a preset plus run-specific params:

```bash
nvd run --preset chtc-defaults --params-file run.yaml
```

Run with explicit flags:

```bash
nvd run \
  --samplesheet samplesheet.csv \
  --experiment-id experiment-001 \
  --blast-db /path/to/blast_db \
  --blast-db-prefix core_nt \
  --virus-index /path/to/human_infecting_viruses.k31w1.idx \
  --taxonomy-dir /path/to/taxdump \
  --profile docker
```

Dry-run to see the generated Nextflow command without executing it:

```bash
nvd run --params-file run.yaml --dry-run
```

Pass Nextflow-only flags after `--`:

```bash
nvd run --params-file run.yaml -- -with-report -with-trace
```

Resume using the same CLI inputs:

```bash
nvd run --params-file run.yaml --resume
```

Resume the exact previous command cached by `nvd run`:

```bash
nvd resume
nvd resume --interactive
```

The short alias for `run` is available as `nvd r`, but examples use `nvd run` for clarity.

## Validation and diagnostics

Check core dependencies:

```bash
nvd validate deps
```

Validate a samplesheet through the validation namespace:

```bash
nvd validate samplesheet samplesheet.csv
```

Validate reference paths from config/params where applicable:

```bash
nvd validate databases --config ~/.nvd/user.config
```

Run the validation checks exposed by the validation namespace:

```bash
nvd validate all --samplesheet samplesheet.csv
```

For params files, prefer the params-native command because it understands presets and source precedence:

```bash
nvd params check run.yaml --preset chtc-defaults --no-check-paths
```

Show the installed CLI and dependency versions:

```bash
nvd version
nvd v
```

## Secrets

The secrets namespace wraps Nextflow secrets for values that should not be placed in params files or configs, such as tokens.

```bash
nvd secrets list
nvd secrets set SLACK_BOT_TOKEN
nvd secrets check
nvd secrets delete SLACK_BOT_TOKEN
```

Aliases are available for listing and removal:

```bash
nvd secrets ls
nvd secrets rm SLACK_BOT_TOKEN
```

## Command surface summary

```text
nvd run / nvd r                         launch the pipeline
nvd resume                              rerun the previous cached command with resume
nvd setup                               configure wrapper, shell hook, CHTC defaults, completions, SIF
nvd setup shell-hook                    print shell integration code
nvd params init                         create schema-backed YAML/JSON params files
nvd params check                        validate params and show merged sources
nvd params merge                        merge params files
nvd preset list / ls                    list presets
nvd preset show                         inspect one preset
nvd preset register / reg               create or update a preset
nvd preset merge                        combine presets into a new preset
nvd preset diff                         compare two presets
nvd preset export                       write a preset to YAML/JSON
nvd preset import                       import a preset from YAML/JSON
nvd preset delete / rm                  remove a preset
nvd samplesheet generate / gen          generate a samplesheet
nvd samplesheet validate / val          validate a samplesheet
nvd ss gen / nvd ss val                 short samplesheet aliases
nvd taxonomy ensure                     prepare taxonomy data
nvd taxonomy status                     inspect taxonomy data
nvd config show / path / edit           inspect or edit discovered config
nvd validate deps / samplesheet / databases / all / params
nvd secrets list / set / get / delete / check
nvd version / nvd v                     show version information
```

## Direct Nextflow fallback

Use direct Nextflow only when you are deliberately debugging the workflow layer or reproducing a lower-level issue. A direct command should mirror what `nvd run --dry-run` prints:

```bash
nextflow run . \
  -profile docker \
  -params-file run.yaml
```

For ordinary use, prefer:

```bash
nvd run --params-file run.yaml --profile docker
```
