set shell := ["bash", "-cu"]

PYTHON_CI_PATHS := "lib/"

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

# enter the locked Pixi development shell
shell:
    pixi shell --frozen

# === Day-to-day Checks ===

# run the usual local checks, excluding the slow networked integration test
check: fmt-check lint schema-check test config-check
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

# validate that the params schema covers pipeline params
schema-check:
    uv run python .github/scripts/validate_schema_completeness.py

# run the fast pytest suite only
test:
    uv run pytest -m "not slow and not network"

# run one pytest file or node, e.g. just test-one lib/py_nvd/test_models.py
test-one path:
    uv run pytest "{{ path }}"

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

# === End-to-end Integration ===

# run the slow mini SRA end-to-end test with progress output
e2e profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" pixi run e2e-test

# run the slow mini SRA end-to-end test as CI does
e2e-ci profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" pixi run e2e-test-ci

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
alias nf-config := config
alias config-test := config
alias clean := clean-e2e
alias build := docker-build
