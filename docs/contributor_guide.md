# NVD2 Contributor Guide

Welcome to the NVD2 development community! This guide will help you set up a
development environment and understand our project structure, tooling, and
contribution workflow.

## Table of Contents

- [Quick Start](#quick-start)
- [Development Environment Setup](#development-environment-setup)
- [Project Structure](#project-structure)
- [Development Workflow](#development-workflow)
- [Testing](#testing)
- [Code Style and Quality](#code-style-and-quality)
- [Understanding the Gitignore Logic](#understanding-the-gitignore-logic)
- [Contributing Guidelines](#contributing-guidelines)

## Quick Start

**The fastest way to start contributing:**

1. **Install Pixi** (recommended for full development):
   ```bash
   curl -fsSL https://pixi.sh/install.sh | bash
   ```

2. **Clone and enter the project**:
   ```bash
   git clone https://github.com/dhoconno/nvd.git
   cd nvd
   pixi install --frozen
   pixi shell
   ```

3. **Start developing**:
   ```bash
   python -c "import nvd2; print('Ready to contribute!')"
   ```

**Alternative setups:**

- **For Python-only development** (working on `bin/` scripts):
  ```bash
  # Install uv (ultra-fast Python package manager)
  curl -LsSf https://astral.sh/uv/install.sh | sh

  # Clone and setup
  git clone https://github.com/dhoconno/nvd.git
  cd nvd
  uv sync
  source .venv/bin/activate
  python -c "import nvd2; print('Ready for Python development!')"
  ```

- **For system-level reproducibility** (Nix alternative):
  ```bash
  # Install Nix and direnv
  curl -L https://nixos.org/nix/install | sh
  curl -sfL https://direnv.net/install.sh | bash

  # Clone and auto-setup
  git clone https://github.com/dhoconno/nvd.git
  cd nvd
  direnv allow  # Automatically sets up environment
  ```

## Development Environment Setup

NVD2 uses a modern, reproducible development environment with multiple tool
options:

### Primary Development Stack

- **[Pixi](https://pixi.sh/)**: Recommended package manager for reproducible
  environments (conda-forge + PyPI)
- **[uv](https://docs.astral.sh/uv/)**: Ultra-fast Python package manager for
  Python-only development
- **[Nix](https://nixos.org/)**: System-level reproducibility and development
  shell (alternative)
- **[direnv](https://direnv.net/)**: Automatic environment activation (with Nix)
- **Python 3.11+**: Core language requirement
- **Nextflow**: Workflow orchestration engine

### Setup Options

#### Option 1: Pixi (Recommended for Full Development)

Best for developing Nextflow workflows, Python scripts, and working with
bioinformatics tools:

```bash
# Install Pixi
curl -fsSL https://pixi.sh/install.sh | bash

# Setup project
git clone https://github.com/dhoconno/nvd.git
cd nvd
pixi install --frozen
pixi shell
```

**What gets installed:**

- All bioinformatics tools (samtools, BLAST, GOTTCHA2, SPAdes, etc.)
- Nextflow workflow engine
- Python dependencies and development tools
- Container runtimes (Apptainer on Linux)

#### Option 2: uv (Recommended for Python-Only Development)

Perfect when you only need to work on Python scripts in the `bin/` directory:

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Setup project
git clone https://github.com/dhoconno/nvd.git
cd nvd
uv sync  # Creates virtual environment and installs dependencies
```

**Working with uv:**

```bash
# Activate the virtual environment
source .venv/bin/activate

# Now run commands normally (no uv run prefix needed)
python bin/script_name.py
pytest
ruff format .

# Install new dependencies (from outside or inside venv)
uv add package_name

# Deactivate when done
deactivate
```

**When to use uv:**

- You're only modifying Python scripts in `bin/`
- You don't need bioinformatics tools installed locally
- You want the fastest possible Python dependency management
- You're doing pure Python development work

#### Option 3: Nix + direnv (Alternative for System-Level Reproducibility)

For maximum reproducibility across different systems:

```bash
# Install prerequisites
curl -L https://nixos.org/nix/install | sh
curl -sfL https://direnv.net/install.sh | bash

# Setup (automatic environment activation)
git clone https://github.com/dhoconno/nvd.git
cd nvd
direnv allow  # Automatically sets up Nix shell + Pixi environment
```

The `.envrc` file automatically:

- Enters a Nix development shell with system dependencies
- Installs Pixi and runs `pixi install --frozen`
- Sets up the Python environment with all bioinformatics tools

#### Option 4: Manual Installation (Not Recommended)

For compatibility, you can use standard Python tooling, but you'll need to
manually install bioinformatics tools:

```bash
git clone https://github.com/dhoconno/nvd.git
cd nvd
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
# Manually install: samtools, minimap2, seqkit, bbmap, blast, gottcha2, spades, nextflow, etc.
```

### What Gets Installed

The development environment includes:

**Bioinformatics Tools:**

- samtools, minimap2, seqkit, bbmap, vsearch
- BLAST, GOTTCHA2, SPAdes
- Nextflow (workflow engine)

**Python Dependencies:**

- Core: biopython, pandas, polars, pysam
- Analysis: lxml, more-itertools, pyyaml
- Integration: labkey, snakemake
- CLI: typer, loguru

**Development Tools:**

- Testing: pytest, tox
- Code Quality: ruff (linting/formatting), basedpyright (type checking)
- Notebooks: jupyter, marimo
- Containers: apptainer (Linux only)

### Building NCBI STAT (Advanced)

**Note**: NCBI STAT is pre-built in the Docker/Apptainer containers. You only
need to build it manually if:

- You're using Option 2 (uv) or Option 4 (manual) setup
- You want to use the `nvd` workflow locally without containers
- You're contributing to STAT-related functionality

**Prerequisites:**

```bash
# Install build dependencies
sudo apt-get update
sudo apt-get install -y curl wget git gcc g++ cmake util-linux

# Or on macOS
brew install cmake gcc wget git
```

**Build process:**

```bash
# Create build directory
mkdir /tmp/ncbi-build && cd /tmp/ncbi-build

# Clone required repositories
git config --global http.sslVerify false  # If needed for firewalls
git clone -b 3.2.0 https://github.com/ncbi/ngs-tools.git
git clone -b 3.2.1 https://github.com/ncbi/ncbi-vdb.git
git clone -b 3.2.1 https://github.com/ncbi/sra-tools.git

# Fix CMake version requirement in ngs-tools
sed -i 's/cmake_minimum_required[[:space:]]*([[:space:]]*VERSION[[:space:]]*2\.8\.12[[:space:]]*)/cmake_minimum_required(VERSION 3.5)/' ngs-tools/CMakeLists.txt

# Build ncbi-vdb
cd ncbi-vdb
./configure --relative-build-out-dir
make -j$(nproc)

# Build sra-tools
cd ../sra-tools
./configure --relative-build-out-dir
make -j$(nproc)

# Build ngs-tools (contains STAT)
cd ../ngs-tools
./configure --relative-build-out-dir
make -j$(nproc)

# Install binaries to your PATH
sudo find OUTDIR -type f -executable -exec cp {} /usr/local/bin/ \;

# Clean up
cd / && rm -rf /tmp/ncbi-build
```

**Verify installation:**

```bash
# Check that STAT tools are available
which aligns_to
which tax_analysis
```

**NCBI configuration:**

```bash
# Create NCBI configuration directory
mkdir -p ~/.ncbi

# Copy the provided configuration (if using the repo version)
cp conf/user-settings.mkfg ~/.ncbi/user-settings.mkfg
```

**When STAT build fails:**

- Check that all dependencies are installed
- Ensure you have sufficient disk space (build requires ~2GB)
- On macOS, you may need to set `CC=gcc-13` if using Homebrew GCC
- For compilation errors, check that CMake version >= 3.5

## Project Structure

```
nvd/
├── bin/                    # Python CLI scripts and utilities
├── workflows/              # Main Nextflow workflow definitions
├── subworkflows/           # Reusable Nextflow subworkflow modules
├── modules/                # Individual Nextflow process modules
├── conf/                   # Configuration files
├── assets/                 # Example data and static resources
├── docs/                   # Documentation (you are here!)
├── scripts/                # Utility scripts (shell, awk, perl, python, lua)
├── legacy/                 # Original Snakemake implementation
├── pyproject.toml          # Python package configuration and dependencies
├── pixi.lock              # Locked dependency versions for reproducibility
├── flake.nix              # Nix development environment specification
├── .envrc                 # direnv automatic environment activation
├── nextflow.config        # Nextflow pipeline configuration
└── .gitignore             # Inverted gitignore (see section below)
```

### Key Files

- **`main.nf`**: Main Nextflow pipeline entry point
- **`pyproject.toml`**: Defines Python dependencies, development tools, and
  project metadata
- **`pixi.lock`**: Locked versions of all conda/PyPI dependencies for
  reproducibility
- **`flake.nix`**: Nix flake for system-level reproducible development
  environment
- **`.envrc`**: direnv configuration that automatically activates the Nix dev
  shell

## Development Workflow

### Making Changes

1. **Create a branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** to the relevant files:
   - **Python scripts**: `bin/` directory
   - **Nextflow workflows**: `workflows/`, `subworkflows/`, `modules/`
   - **Configuration**: `conf/`, `nextflow.config`
   - **Documentation**: `docs/`, `README.md`

3. **Test your changes**:

   **With Pixi (full development):**
   ```bash
   # In pixi shell, run commands directly
   pytest

   # Test Nextflow syntax
   nextflow run . --help

   # Run a small test dataset
   nextflow run . --tools nvd --samplesheet assets/example_samplesheet.csv [other params]
   ```

   **With uv (Python-only development):**
   ```bash
   # Activate environment first
   source .venv/bin/activate

   # Run Python tests
   pytest

   # Test Python scripts directly
   python bin/your_script.py --help

   # Run specific tests
   pytest tests/test_your_module.py
   ```

4. **Format and lint your code**:

   **With Pixi:**
   ```bash
   # In pixi shell, run commands directly
   ruff format .
   ruff check .
   basedpyright .  # optional
   ```

   **With uv:**
   ```bash
   # Activate environment first
   source .venv/bin/activate

   # Format and lint Python code
   ruff format .
   ruff check .
   basedpyright .  # optional
   ```

### Adding New Dependencies

#### Python Dependencies

Add to `pyproject.toml`:

```toml
dependencies = [
  # ... existing deps
  "new-package>=1.0.0",
]
```

#### Bioinformatics Tools

Add to `pyproject.toml` under `[tool.pixi.dependencies]`:

```toml
[tool.pixi.dependencies]
# ... existing tools
new-tool = ">=2.0.0,<3"
```

#### Development Dependencies

Add to `pyproject.toml` under `[dependency-groups]`:

```toml
[dependency-groups]
dev = [
  # ... existing dev deps
  "new-dev-tool>=1.0.0",
]
```

**After adding dependencies:**

**With Pixi:**

```bash
pixi install  # Update the environment
```

**With uv:**

```bash
uv add package_name  # Add Python dependencies
uv sync              # Update environment after manual pyproject.toml edits
# No need to reactivate - changes apply immediately to active environment
```

## Testing

### Running Tests

**With Pixi (recommended for full testing):**

```bash
# In pixi shell, run commands directly
pytest

# Run specific test file
pytest tests/test_specific.py

# Run with coverage
pytest --cov=bin

# Run tests across multiple Python versions
tox
```

**With uv (Python-only development):**

```bash
# Activate environment first
source .venv/bin/activate

# Run all Python tests
pytest

# Run specific test file  
pytest tests/test_specific.py

# Run with coverage
pytest --cov=bin

# Run tests with verbose output
pytest -v
```

### Testing Nextflow Workflows

```bash
# In pixi shell, run commands directly
nextflow run . --help

# Test with example data
nextflow run . \
  --tools nvd \
  --samplesheet assets/example_samplesheet.csv \
  --blast_db db \
  --blast_db_prefix PP819512-nt \
  --stat_index db/tree_index.dense.dbs \
  --stat_dbss db/tree_index.dense.dbss \
  --stat_annotation db/tree_index.dense.dbss.annotation \
  --human_virus_taxlist db/human_viruses_taxlist.txt \
  --experiment_id 99999

# Validate workflow DAG
nextflow run . -preview
```

### Adding Tests

Create test files in a `tests/` directory:

```python
# tests/test_your_module.py
import pytest
from bin.your_module import your_function

def test_your_function():
    result = your_function("input")
    assert result == "expected_output"
```

## Code Style and Quality

### Python Code Style

We use **Ruff** for both formatting and linting:

**With Pixi:**

```bash
# In pixi shell, run commands directly
ruff format .
ruff check .
ruff check --fix .
```

**With uv:**

```bash
# Activate environment first
source .venv/bin/activate

# Format and lint Python code
ruff format .
ruff check .
ruff check --fix .
```

### Configuration

Our `pyproject.toml` includes minimal Ruff configuration:

```toml
[tool.ruff.lint]
select = []  # Use default rules
ignore = ["G004", "PTH"]  # Ignore specific rules

[tool.basedpyright]
typeCheckingMode = "off"  # Type checking currently disabled
```

### Nextflow Style

- Use meaningful process names: `EXTRACT_HUMAN_VIRUS_READS` not `process1`
- Include helpful comments and documentation
- Use consistent indentation (4 spaces)
- Group related processes in subworkflows

## Understanding the Gitignore Logic

**Important**: NVD2 uses an **inverted gitignore strategy**. This is crucial to
understand!

### How It Works

Instead of listing files to ignore, we:

1. **Ignore everything** with `*` at the top
2. **Explicitly allow** needed files with `!pattern`

```gitignore
# Block everything
*

# Allow specific files/directories
!.gitignore
!README.md
!main.nf
!nextflow.config
!pyproject.toml

# Allow entire directories
!/workflows
!/workflows/*.nf

!/subworkflows
!/subworkflows/*.nf

!/bin
!/bin/*.py

!/docs
!/docs/*.sh
!/docs/*.md
```

### Why This Approach?

1. **Prevents accidental commits** of large datasets, intermediate files, etc.
2. **Forces intentional decisions** about what should be version controlled
3. **Keeps the repository clean** as the project grows
4. **Avoids the Sisyphean task** of continuously adding new ignore patterns

### Adding New Files to Git

When you create new files that should be tracked:

1. **Check if your file type is already allowed**:
   ```bash
   # These are already allowed
   touch workflows/new_workflow.nf      # *.nf in workflows/
   touch bin/new_script.py             # *.py in bin/
   touch docs/new_doc.md               # *.md in docs/
   ```

2. **If you need a new file type or location**, update `.gitignore`:
   ```gitignore
   # Add to .gitignore
   !/conf
   !/conf/*.config
   ```

3. **For completely new directories**, add both the directory and its contents:
   ```gitignore
   !/tests
   !/tests/*.py
   ```

### Common Patterns Already Allowed

- **Python**: `!/bin/*.py`
- **Nextflow**: `!/workflows/*.nf`, `!/subworkflows/*.nf`, `!/modules/*.nf`
- **Configuration**: `!/conf/*.config`
- **Documentation**: `!/docs/*.md`, `!/docs/*.sh`
- **Scripts**: `!/scripts/*.{sh,awk,pl,py,lua}`
- **Assets**: `!/assets/*.csv`

## Contributing Guidelines

### Before Submitting

1. **Ensure tests pass**: `pixi run pytest`
2. **Format code**: `pixi run ruff format .`
3. **Check style**: `pixi run ruff check .`
4. **Test with real data** if possible
5. **Update documentation** if needed

### Pull Request Process

1. **Create a descriptive branch name**: `feature/add-nanopore-support` or
   `fix/memory-leak-in-blast`
2. **Write a clear PR description** explaining:
   - What the change does
   - Why it's needed
   - How to test it
3. **Reference any related issues**
4. **Keep PRs focused** - one feature or fix per PR
5. **Update relevant documentation**

### Code Review

- PRs require at least one approval
- Address all review comments
- Ensure CI passes
- Be responsive to feedback

## Getting Help

### Documentation

- **Main README**: Project overview and basic usage
- **Example Commands**: `docs/example_commands.md` - comprehensive usage
  examples
- **This Guide**: Development environment and contribution workflow

### Community

- **Issues**: Use GitHub issues for bugs and feature requests
- **Discussions**: Use GitHub discussions for questions and ideas
- **Email**: Contact maintainers for urgent issues

### Troubleshooting

#### Environment Issues

**With Pixi:**

```bash
# Reset Pixi environment
rm -rf .pixi
pixi install --frozen
pixi shell

# Check tool availability (in pixi shell)
which nextflow
which samtools
python -c "import biopython; print('OK')"
```

**With uv:**

```bash
# Reset uv environment
rm -rf .venv
uv sync
source .venv/bin/activate

# Check Python dependencies (in activated venv)
python -c "import biopython; print('OK')"
python -c "import nvd2; print('Package installed')"
```

**With Nix:**

```bash
# Reset Nix environment (if using Nix)
direnv reload

# Rebuild Nix environment if needed
nix-collect-garbage
direnv allow
```

#### Nextflow Issues

```bash
# Clear Nextflow cache
rm -rf .nextflow*
rm -rf work/

# Test basic functionality (in pixi shell)
nextflow info
nextflow run hello
```

#### Git Issues with Inverted Gitignore

```bash
# See what files are ignored
git status --ignored

# Force add a file that should be tracked
git add -f path/to/file

# Check if .gitignore needs updating
git check-ignore -v path/to/file
```

## Development Tips

1. **Environment management:**
   - **Pixi**: Use `pixi shell` to enter the environment, then run commands
     directly
   - **uv**: Use `source .venv/bin/activate` to activate, then run commands
     normally
   - **Nix**: `direnv` automatically activates/deactivates environment when you
     `cd`

2. **Testing workflow:**
   - **Test small changes** with the example samplesheet first
   - **Use `nextflow -resume`** to restart failed pipeline runs
   - **Keep work directories** during development (set `cleanup = false`)

3. **Python development:**
   - **Use uv for speed** when only working on Python scripts
   - **Switch to Pixi** when you need to test full workflows
   - **Always activate your environment** before running commands

4. **General best practices:**
   - **Monitor resource usage** during testing
   - **Use meaningful commit messages**
   - **Document complex algorithms** in code comments
   - **Add tests for new functionality**

---

**Happy Contributing!** 🧬

This project builds important tools for public health surveillance and
environmental monitoring. Your contributions help researchers worldwide detect
and track pathogens in complex metagenomic samples.

For questions or issues with this guide, please open a GitHub issue or contact
the maintainers.
