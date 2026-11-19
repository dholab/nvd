set shell := ["bash", "-cu"]

PYTHON_CI_PATHS := "lib/"
REPORTING_TYPECHECK_PATHS := "lib/py_nvd/read_inputs.py lib/py_nvd/multiqc_report.py lib/py_nvd/multiqc_packages.py lib/py_nvd/multiqc_fastx.py lib/py_nvd/multiqc_assembly.py lib/py_nvd/multiqc_query_preparation.py lib/py_nvd/multiqc_blast.py lib/py_nvd/multiqc_taxonomy.py lib/py_nvd/multiqc_domains.py bin/write_nvd_fastqc_receipt.py bin/package_nvd_report.py bin/assess_long_read_assembly.py"

# Default recipe: show available commands
default:
    @just --list

# choose recipes interactively
choose:
    @just --choose

# === Project Setup ===

# install the locked Pixi environment and sync Python dependencies
setup:
    pixi install --frozen
    uv sync

# install the locked Pixi environment with explicit development dependencies
setup-dev:
    pixi install -e dev --frozen

# enter the locked Pixi development shell
shell:
    pixi shell --frozen

# enter the locked Pixi development shell
shell-dev:
    pixi shell -e dev --frozen

# === Day-to-day Checks ===

# run the usual local checks, excluding the slow networked integration test
check: fmt-check lint typecheck-reporting schema-check test config-check
    @echo "Project checks passed"

# check Python formatting for paths covered by CI
fmt-check:
    uv run ruff format --check {{ PYTHON_CI_PATHS }}

# format Python paths covered by CI
fmt:
    uv run ruff format {{ PYTHON_CI_PATHS }}

# lint Python paths covered by CI
lint:
    uv run ruff check {{ PYTHON_CI_PATHS }}

# lint and apply safe Python fixes for paths covered by CI
lint-fix:
    uv run ruff check {{ PYTHON_CI_PATHS }} --fix

# statically check the typed MultiQC reporting boundary
typecheck-reporting:
    pixi run -e dev ty check {{ REPORTING_TYPECHECK_PATHS }}

# validate that the params schema covers pipeline params
schema-check:
    uv run python .github/scripts/validate_schema_completeness.py

# run the fast pytest suite only
test:
    pixi run -e dev fast-test

# run one pytest file or node, e.g. just test-one lib/py_nvd/test_models.py
test-one path:
    pixi run -e dev pytest "{{ path }}"

# === Nextflow Development ===

# show the Pixi-provided Nextflow version
nextflow-version:
    pixi run nextflow -version

# render a single Nextflow config profile
config profile="test":
    pixi run nextflow config -profile "{{ profile }}"

# validate all repo-defined Nextflow config profiles render
config-check:
    @for profile in standard docker apptainer chtc_hpc local test; do \
        echo "Checking Nextflow config profile: ${profile}"; \
        pixi run nextflow config -profile "${profile}" > /dev/null; \
    done

# === Reference Builds ===

# === End-to-end Integration ===

# run the slow mini SRA end-to-end test with progress output
e2e profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test with experimental features enabled
e2e-experimental profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" NVD_INTEGRATION_EXPERIMENTAL=1 pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test with contigs as the only BLAST query class
e2e-skip-unassembled-read-queries profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" NVD_INTEGRATION_SKIP_UNASSEMBLED_READ_QUERIES=1 pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test without scheduling SPAdes assembly
e2e-skip-assembly profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" NVD_INTEGRATION_SKIP_ASSEMBLY=1 pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test with experimental features but no SPAdes assembly
e2e-experimental-skip-assembly profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" NVD_INTEGRATION_EXPERIMENTAL=1 NVD_INTEGRATION_SKIP_ASSEMBLY=1 pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test as CI does
e2e-ci profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" pixi run -e dev e2e-test-ci

# run the opt-in real-LabKey e2e (needs the labkey-e2e preset; prompts for LABKEY_API_KEY if unset)
e2e-labkey profile="test":
    @if ! pixi run nvd secrets check > /dev/null 2>&1; then \
        echo "LABKEY_API_KEY secret not set; prompting (input hidden)..."; \
        pixi run nvd secrets set LABKEY_API_KEY; \
    fi
    NVD_INTEGRATION_PROFILE="{{ profile }}" pixi run -e dev e2e-labkey-test

# print the latest end-to-end run directory
e2e-latest:
    @if [ -f .e2e/latest.txt ]; then \
        cat .e2e/latest.txt; \
    else \
        echo "No .e2e/latest.txt found; run 'just e2e' first."; \
        exit 1; \
    fi

# remove local end-to-end test output
clean-e2e:
    rm -rf .e2e

# === Container Development ===

# build the NVD container image locally
docker-build tag="nvd:test":
    docker build -f Containerfile -t "{{ tag }}" .

# smoke-test a locally built NVD container image
docker-test tag="nvd:test":
    docker run --rm "{{ tag }}" bash -c 'nvd --help && nvd version'

alias t := test
alias c := check
alias e := e2e
alias e2e-test := e2e
alias test-drive := e2e
alias ci-e2e := e2e-ci
alias e2e-exp := e2e-experimental
alias e2e-lite := e2e-skip-assembly
alias e2e-exp-lite := e2e-experimental-skip-assembly
alias nf-config := config
alias config-test := config
alias clean := clean-e2e
alias build := docker-build
