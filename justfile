set shell := ["bash", "-cu"]

PYTHON_CI_PATHS := "lib/ scripts/"
REPORTING_TYPECHECK_PATHS := "lib/py_nvd/read_inputs.py lib/py_nvd/multiqc_report.py lib/py_nvd/multiqc_packages.py lib/py_nvd/multiqc_fastx.py lib/py_nvd/multiqc_assembly.py lib/py_nvd/multiqc_domains.py bin/write_nvd_fastqc_receipt.py bin/package_nvd_report.py bin/assess_long_read_assembly.py"
NCBI_VIRUS_SOURMASH_SIG_URL := "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/ncbi-viruses-2025.01/ncbi-viruses-2025.01.dna.k=31.sig.zip"
NCBI_VIRUS_SOURMASH_LINEAGES_URL := "https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/ncbi-viruses-2025.01/ncbi-viruses.2025.01.lineages.csv"
WVDB_FASTA_URL := "https://zenodo.org/records/20276352/files/WVDB_v1.0.fasta?download=1"
WVDB_ANNOTATIONS_URL := "https://zenodo.org/records/20276352/files/WVDB_v1.0_annotations.tsv?download=1"

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
    pixi run -e dev pytest -m "not slow and not network"

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

# fetch and validate the NCBI Virus sourmash reference; downloads large inputs
fetch-sourmash-ncbi-virus approval="" build_dir="build/sourmash/ncbi-virus-wvdb-v1.0":
    @if [ "{{ approval }}" != "yes" ]; then echo "This recipe downloads large reference files. Re-run with: just fetch-sourmash-ncbi-virus yes"; exit 2; fi
    mkdir -p "{{ build_dir }}/.work/ncbi"
    rm -f "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.manifest.csv"
    curl -fL "{{ NCBI_VIRUS_SOURMASH_SIG_URL }}" -o "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.dna.k31.scaled50.sig.zip"
    curl -fL "{{ NCBI_VIRUS_SOURMASH_LINEAGES_URL }}" -o "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.lineages.csv"
    @if file "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.dna.k31.scaled50.sig.zip" | grep -qi 'HTML'; then echo "Downloaded NCBI Virus sourmash sketch is HTML, not a sourmash zip. The upstream Farm URL may be serving a maintenance page."; exit 1; fi
    @if file "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.lineages.csv" | grep -qi 'HTML'; then echo "Downloaded NCBI Virus lineages file is HTML, not CSV. The upstream Farm URL may be serving a maintenance page."; exit 1; fi
    pixi run sourmash sig summarize "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.dna.k31.scaled50.sig.zip"
    pixi run sourmash sig manifest "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.dna.k31.scaled50.sig.zip" -o "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.manifest.csv"
    pixi run sourmash tax prepare --taxonomy-csv "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.lineages.csv" --keep-identifier-versions -F csv -o "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.lineages.validated.csv"

# build and validate the WVDB sourmash reference side; downloads large WVDB inputs
build-sourmash-wvdb approval="" build_dir="build/sourmash/ncbi-virus-wvdb-v1.0": (fetch-sourmash-ncbi-virus approval build_dir)
    @if [ "{{ approval }}" != "yes" ]; then echo "This recipe downloads large WVDB reference files. Re-run with: just build-sourmash-wvdb yes"; exit 2; fi
    mkdir -p "{{ build_dir }}/.work/wvdb"
    curl -fL "{{ WVDB_FASTA_URL }}" -o "{{ build_dir }}/.work/wvdb/WVDB_v1.0.fasta"
    curl -fL "{{ WVDB_ANNOTATIONS_URL }}" -o "{{ build_dir }}/.work/wvdb/WVDB_v1.0_annotations.tsv"
    rm -f "{{ build_dir }}/.work/wvdb/WVDB_v1.0.normalized.fasta" "{{ build_dir }}/.work/wvdb/wvdb-v1.0.sourmash.lineages.csv" "{{ build_dir }}/.work/wvdb/wvdb-v1.0.dna.k31.scaled50.sig.zip" "{{ build_dir }}/.work/wvdb/wvdb-v1.0.manifest.csv" "{{ build_dir }}/.work/wvdb/wvdb-v1.0.sourmash.lineages.validated.csv"
    uv run scripts/prepare_wvdb_sourmash_inputs.py prepare-inputs --fasta "{{ build_dir }}/.work/wvdb/WVDB_v1.0.fasta" --annotations-tsv "{{ build_dir }}/.work/wvdb/WVDB_v1.0_annotations.tsv" --reference-lineages-csv "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.lineages.csv" --normalized-fasta "{{ build_dir }}/.work/wvdb/WVDB_v1.0.normalized.fasta" --lineages-csv "{{ build_dir }}/.work/wvdb/wvdb-v1.0.sourmash.lineages.csv"
    pixi run sourmash sketch dna "{{ build_dir }}/.work/wvdb/WVDB_v1.0.normalized.fasta" --singleton -p dna,k=31,scaled=50 -o "{{ build_dir }}/.work/wvdb/wvdb-v1.0.dna.k31.scaled50.sig.zip"
    pixi run sourmash sig summarize "{{ build_dir }}/.work/wvdb/wvdb-v1.0.dna.k31.scaled50.sig.zip"
    pixi run sourmash sig manifest "{{ build_dir }}/.work/wvdb/wvdb-v1.0.dna.k31.scaled50.sig.zip" -o "{{ build_dir }}/.work/wvdb/wvdb-v1.0.manifest.csv"
    uv run scripts/prepare_wvdb_sourmash_inputs.py check-manifest-coverage --manifest-csv "{{ build_dir }}/.work/wvdb/wvdb-v1.0.manifest.csv" --lineages-csv "{{ build_dir }}/.work/wvdb/wvdb-v1.0.sourmash.lineages.csv"
    pixi run sourmash tax prepare --taxonomy-csv "{{ build_dir }}/.work/wvdb/wvdb-v1.0.sourmash.lineages.csv" --keep-identifier-versions -F csv -o "{{ build_dir }}/.work/wvdb/wvdb-v1.0.sourmash.lineages.validated.csv"

# build a combined NCBI Virus 2025.01 + WVDB v1.0 sourmash reference; downloads large inputs
build-sourmash-ncbi-wvdb approval="" build_dir="build/sourmash/ncbi-virus-wvdb-v1.0": (build-sourmash-wvdb approval build_dir)
    mkdir -p "{{ build_dir }}/.work/combined"
    rm -rf "{{ build_dir }}/ncbi-viruses-2025.01.dna.k31.scaled50.sig.zip" "{{ build_dir }}/ncbi-viruses-2025.01.lineages.csv" "{{ build_dir }}/ncbi-viruses-2025.01.lineages.validated.csv" "{{ build_dir }}/WVDB_v1.0.fasta" "{{ build_dir }}/WVDB_v1.0_annotations.tsv" "{{ build_dir }}/WVDB_v1.0.normalized.fasta" "{{ build_dir }}/wvdb-v1.0.dna.k31.scaled50.sig.zip" "{{ build_dir }}/wvdb-v1.0.manifest.csv" "{{ build_dir }}/wvdb-v1.0.sourmash.lineages.csv" "{{ build_dir }}/wvdb-v1.0.sourmash.lineages.validated.csv" "{{ build_dir }}/wvdb-cat-smoke.dna.k31.scaled50.sig.zip"
    rm -rf "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.sig" "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.manifest.csv" "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.taxonomy-conflicts.tsv" "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.candidate.sig.zip" "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.candidate.csv" "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.validated.csv" "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.manifest.csv"
    uv run scripts/prepare_wvdb_sourmash_inputs.py curate-lineages --reference-manifest-csv "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.manifest.csv" --wvdb-manifest-csv "{{ build_dir }}/.work/wvdb/wvdb-v1.0.manifest.csv" --reference-lineages-csv "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.lineages.csv" --wvdb-lineages-csv "{{ build_dir }}/.work/wvdb/wvdb-v1.0.sourmash.lineages.csv" --lineages-csv "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.candidate.csv" --diagnostics-tsv "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.taxonomy-conflicts.tsv"
    pixi run sourmash sig cat "{{ build_dir }}/.work/ncbi/ncbi-viruses-2025.01.dna.k31.scaled50.sig.zip" "{{ build_dir }}/.work/wvdb/wvdb-v1.0.dna.k31.scaled50.sig.zip" --dna -k 31 --unique -o "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.candidate.sig.zip"
    pixi run sourmash sig manifest "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.candidate.sig.zip" -o "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.manifest.csv"
    uv run scripts/prepare_wvdb_sourmash_inputs.py validate-lineages --manifest-csv "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.manifest.csv" --lineages-csv "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.candidate.csv"
    pixi run sourmash tax prepare --taxonomy-csv "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.candidate.csv" --keep-identifier-versions -F csv -o "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.validated.csv"
    pixi run sourmash sig summarize "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.candidate.sig.zip"
    rm -f "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.sha256"
    mv "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.candidate.csv" "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.csv"
    mv "{{ build_dir }}/.work/combined/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.candidate.sig.zip" "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.sig.zip"
    uv run scripts/prepare_wvdb_sourmash_inputs.py write-checksums --output "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.sha256" "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.dna.k31.scaled50.sig.zip" "{{ build_dir }}/ncbi-viruses-2025.01-plus-wvdb-v1.0.lineages.csv"

# === End-to-end Integration ===

# run the slow mini SRA end-to-end test with progress output
e2e profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test with experimental features enabled
e2e-experimental profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" NVD_INTEGRATION_EXPERIMENTAL=1 pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test without scheduling SPAdes assembly
e2e-skip-assembly profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" NVD_INTEGRATION_SKIP_ASSEMBLY=1 pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test with experimental features but no SPAdes assembly
e2e-experimental-skip-assembly profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" NVD_INTEGRATION_EXPERIMENTAL=1 NVD_INTEGRATION_SKIP_ASSEMBLY=1 pixi run -e dev e2e-test

# run the slow mini SRA end-to-end test as CI does
e2e-ci profile="test":
    NVD_INTEGRATION_PROFILE="{{ profile }}" pixi run -e dev e2e-test-ci

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
alias wvdb := build-sourmash-ncbi-wvdb
