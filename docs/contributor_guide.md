# NVD Contributor Guide

Welcome to the NVD development guide. This document is for people changing the pipeline, CLI, setup flow, tests, or documentation. It is intentionally focused on the current v3 architecture: the `nvd` CLI is the recommended user interface, and direct `nextflow run` is available mostly as an escape hatch for debugging.

## Quick start

For full-stack development, use Pixi. It installs the Python package, Nextflow, and the bioinformatics tools used by the workflow.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
git clone https://github.com/dhoconno/nvd.git
cd nvd
pixi install --frozen
pixi shell
```

For Python-only work, `uv` is faster and usually enough for CLI, params, preset, setup, taxonomy, and helper-script changes.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
git clone https://github.com/dhoconno/nvd.git
cd nvd
uv sync
```

The Nix/direnv setup remains available for contributors who want a more portable system shell:

```bash
curl -L https://nixos.org/nix/install | sh
curl -sfL https://direnv.net/install.sh | bash
git clone https://github.com/dhoconno/nvd.git
cd nvd
direnv allow
```

## Current usage model

The `nvd` CLI is the preferred interface for installing, configuring, authoring parameters, and launching the pipeline. Contributors should treat direct Nextflow commands as a lower-level fallback, not the main user path.

The common user flow should look like this:

```text
nvd setup
   ↓
nvd params init run.yaml
   ↓
edit run.yaml with schema-backed completion
   ↓
nvd params check run.yaml --no-check-paths
   ↓
nvd run --params-file run.yaml
```

On CHTC, `nvd setup` also writes `~/.nvd/setup.conf`, configures the shared preset store location, writes the taxonomy location, and can help place the SIF under `~/.nvd`. Reference artifacts such as the BLAST database and deacon virus index are still runtime parameters; the installer can download them and print the paths, but it should not create another hidden state layer for those paths.

## Parameter authoring

The params file workflow is central to v3. Prefer documenting and testing new runtime options through the params model rather than through hand-written Nextflow-only examples.

Generate a schema-backed template:

```bash
nvd params init run.yaml
```

The generated YAML starts with a `yaml-language-server` schema comment. Editors such as VS Code can use that schema for completions, type checking, and inline validation if a YAML language server is installed. JSON templates include the same schema through `$schema`.

For VS Code, contributors and users should install a YAML or JSON language server extension. The goal is for common mistakes such as misspelled fields, wrong boolean values, or invalid enum-like strings to be caught while editing, before the user launches a cluster job.

Validate params before running:

```bash
nvd params check run.yaml
nvd params check run.yaml --no-check-paths   # useful when paths exist only on CHTC
```

Multiple params files can be merged when a workflow needs a base config plus a small run-specific override:

```bash
nvd params merge base.yaml run-overrides.yaml --output merged.yaml
```

The runtime precedence is:

```text
CLI flags > params file > preset > pipeline defaults
```

That precedence matters when debugging. If a value appears surprising, start by checking whether it came from a CLI flag, params file, preset, generated config, or default.

## Running NVD during development

Use the CLI for ordinary local or cluster launches:

```bash
nvd run \
  --params-file run.yaml \
  --profile docker
```

or, for a small explicit command:

```bash
nvd run \
  --samplesheet assets/example_samplesheet.csv \
  --experiment-id dev-smoke \
  --blast-db /path/to/blast_db \
  --blast-db-prefix core_nt \
  --virus-index /path/to/human_infecting_viruses.k31w1.idx \
  --taxonomy-dir /path/to/taxdump \
  --profile docker
```

Direct Nextflow execution is still useful for isolating workflow-level issues, but it should usually mirror what `nvd run` would produce:

```bash
nextflow run . \
  -profile docker \
  --samplesheet assets/example_samplesheet.csv \
  --experiment_id dev-smoke \
  --blast_db /path/to/blast_db \
  --blast_db_prefix core_nt \
  --virus_index /path/to/human_infecting_viruses.k31w1.idx \
  --taxonomy_dir /path/to/taxdump
```

The repository does not yet include tiny BLAST, deacon, and taxonomy fixtures for a full local pipeline smoke test. Until those exist, config parsing, unit tests, focused script tests, and CHTC smoke runs carry more of the validation load.

## Project structure

```text
nvd/
├── bin/                    # Python helper scripts run by Nextflow processes
├── conf/                   # Nextflow configuration snippets
├── docs/                   # User and contributor documentation
├── lib/
│   ├── *.groovy            # Nextflow/Groovy helpers
│   └── py_nvd/             # Python CLI, params, presets, taxonomy helpers
├── modules/                # Individual Nextflow process modules
├── schemas/                # JSON schemas for params files
├── subworkflows/           # Reusable Nextflow subworkflow modules
├── workflows/              # Main workflow orchestration
├── .github/scripts/        # CI validation scripts
├── install.sh              # User installer and CHTC setup entrypoint
├── main.nf                 # Nextflow entrypoint
├── nextflow.config         # Pipeline defaults and profile includes
├── pixi.lock               # Locked Pixi environment
├── pyproject.toml          # Python package, Pixi deps, and tool config
└── uv.lock                 # Locked Python dependencies
```

Important files to know:

- `lib/py_nvd/models.py` defines the Pydantic params model.
- `schemas/nvd-params.v3.0.0.schema.json` is the public params schema.
- `lib/py_nvd/cli/commands/run.py` maps CLI flags into params and launches Nextflow.
- `lib/py_nvd/cli/commands/params.py` owns params templates, validation, and merging.
- `lib/py_nvd/cli/commands/setup.py` owns `nvd setup`, CHTC defaults, shell hooks, taxonomy, preset-store, and SIF setup behavior.
- `install.sh` bootstraps the repo and then delegates setup to the CLI.

## Development workflow

When adding or changing a runtime parameter, update the whole parameter surface together:

1. Add or change the field in `lib/py_nvd/models.py`.
2. Update `nextflow.config` if the pipeline needs a default.
3. Update the CLI mapping in `lib/py_nvd/cli/commands/run.py` if users should be able to set it directly.
4. Update `schemas/nvd-params.v3.0.0.schema.json` so params files remain schema-backed.
5. Run `.github/scripts/validate_schema_completeness.py`.
6. Update documentation and examples.

When adding a process or changing process names, check `conf/results.config` so publish rules still point at real process names. CI has a strict process-name check for this because stale selectors silently degrade outputs.

When changing installer/setup behavior, test both normal shell execution and the `curl ... | bash` shape. The checked-in pseudo-terminal regression script exists because ordinary `printf ... | bash install.sh` does not reproduce piped installer prompts correctly.

## Testing

Useful checks before handing off a change:

```bash
uv run pytest . -m "not slow"
uv run python .github/scripts/validate_schema_completeness.py
uv run python .github/scripts/test_installer_interactive.py
bash -n install.sh && shellcheck install.sh
```

Parse the supported profiles after Nextflow changes:

```bash
nextflow config -profile standard
nextflow config -profile docker
nextflow config -profile apptainer
nextflow config -profile chtc_hpc
nextflow config -profile local
```

For focused Python changes, run the relevant tests directly with `uv run pytest` and use `ruff` on the files you touched. A broad `ruff check bin lib/py_nvd` may surface older unrelated issues, so do not treat a focused cleanup as an excuse to mix in a style migration across the whole tree.

For installer changes, the important test is:

```bash
uv run python .github/scripts/test_installer_interactive.py
```

That test verifies the standard piped installer shape, including declining a checkout update and entering the reference wizard in dry-run mode without downloading large artifacts.

## Code style

Python code uses Ruff for linting and formatting. The project currently selects Ruff's full rule set with targeted ignores for this codebase, so expect new code to be explicit about types, paths, subprocess use, and error handling.

Nextflow code should keep orchestration in subworkflows where possible. Prefer stable DAGs controlled by `Channel.empty()` and `Channel.value()` over large top-level `if/else` blocks or process-level `when:` guards. Use variable assignment rather than postfix `.set {}`. Process names should be descriptive and stable because result publishing and CI checks depend on them.

Avoid hiding dependencies in local imports unless there is a clear reason. Most imports should live at module top level so readers can see what a file depends on without reading every function body.

## Gitignore strategy

NVD uses an inverted `.gitignore`: everything is ignored by default, and tracked paths are explicitly allowed. This prevents accidental commits of large data, work directories, downloaded references, SIFs, logs, and analysis outputs.

If you add a new tracked file type or directory, update `.gitignore` in the same commit. Common tracked paths already include:

- root project files such as `README.md`, `install.sh`, `pyproject.toml`, `pixi.lock`, `uv.lock`, `Containerfile`, and `.dockerignore`
- Nextflow files under `workflows/`, `subworkflows/`, and `modules/`
- Python scripts under `bin/` and Python package files under `lib/py_nvd/`
- schemas under `schemas/`
- CI workflows and helper scripts under `.github/`
- selected example assets under `assets/`

When a file appears to be missing from version control, check whether `.gitignore` needs an allow rule before forcing it in.

## Troubleshooting development environments

Reset Pixi if the full development environment is confused:

```bash
rm -rf .pixi
pixi install --frozen
pixi shell
```

Reset uv if the Python environment is confused:

```bash
rm -rf .venv
uv sync
```

Clear local Nextflow state when debugging profile/config issues:

```bash
rm -rf .nextflow* work/
nextflow config -profile docker
```

On CHTC, be explicit about which paths live in user home, shared staging, and per-run work/results directories. `NVD_CONFIG_DIR`, `NVD_TAXONOMY_DB`, and `NVD_PRESET_STORE` are setup/runtime integration points; BLAST and deacon reference paths should still be provided through params, presets, or CLI flags.

## Before submitting

Before asking for review, make sure the change has the right scope and the tests match the kind of code touched. At minimum, run focused tests for your change and the schema/profile checks if you touched params, config, or workflow structure.

Update docs whenever user-facing CLI flags, params, setup behavior, reference requirements, or CHTC behavior changes. The recommended docs path should remain CLI-first, schema-backed params second, and direct Nextflow only as a fallback.

This project supports public-health and environmental monitoring workflows. Be especially careful with changes that can mutate existing installs, submit cluster jobs, upload to LabKey, or write shared configuration.
