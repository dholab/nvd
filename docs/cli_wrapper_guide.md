# NVD2 CLI Wrapper Guide

**The modern, user-friendly interface for NVD2**

This guide shows you how to use the **NVD2 CLI wrapper** (`nvd` command), 
a streamlined Python interface that simplifies running the NVD2 Nextflow pipeline.

The CLI wrapper integrates seamlessly with the `install.sh` script and provides
validation tools, auto-configuration, and a clean command-line experience.

> **Note**: For traditional direct Nextflow execution, see [example_commands.md](./example_commands.md).

## Table of Contents

- [Why Use the CLI Wrapper?](#why-use-the-cli-wrapper)
- [Installation and Setup](#installation-and-setup)
- [Quick Start](#quick-start)
- [Configuration Management](#configuration-management)
- [Running the Pipeline](#running-the-pipeline)
- [Validation Commands](#validation-commands)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)

## Why Use the CLI Wrapper?

The `nvd` CLI wrapper provides several advantages over direct Nextflow execution:

### **Simplified Interface**
```bash
# CLI wrapper - clean and intuitive
nvd run -s samples.csv -e exp001

# Direct Nextflow - verbose and error-prone
nextflow run dhoconno/nvd \
    --samplesheet samples.csv \
    --experiment_id exp001 \
    --blast_db /long/path/to/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /long/path/to/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /long/path/to/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /long/path/to/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /long/path/to/STAT_db/human_viruses_taxlist.txt \
    --gottcha2_db /long/path/to/gottcha2.fna
```

### **Smart Defaults**
- **Auto-detects execution profile** (Docker, Apptainer, or local)
- **Loads configuration** from `~/.nvd/user.config` automatically
- **Validates inputs** before execution (saves time catching errors early)
- **100% parameter coverage** - all Nextflow parameters available as CLI options

### **Integration with install.sh**
The CLI wrapper works seamlessly with the `install.sh` script:

1. **install.sh** sets up databases and creates `~/.nvd/user.config`
2. **nvd CLI** reads that config automatically - no need to specify paths every time
3. **Validation tools** verify everything is configured correctly

This tight integration means you configure once and run everywhere!

## Installation and Setup

### Option 1: Using install.sh

The interactive setup script helps configure database paths:

```bash
# Download and run
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash

# Prerequisites (must be installed separately):
# - Java 11+
# - Nextflow  
# - Docker or Apptainer/Singularity

# The script will:
# 1. Check that prerequisites are installed
# 2. Help you configure database paths
# 3. Create ~/.nvd/user.config
# 4. Optionally download reference databases
```

After setup, install dependencies and verify:

```bash
cd nvd

# Install Python dependencies only
uv sync

# OR install all Conda dependencies (for full toolchain)
pixi shell --frozen

# Verify everything works
nvd version              # Check CLI version
nvd validate deps        # Verify dependencies
nvd validate databases   # Check database paths
```

**Note**: When not in an active shell (`pixi shell`), prepend commands with `pixi run` 
for full toolchain access or `uv run` for Python-only scripts.

### Option 2: Manual Setup

If you prefer manual configuration:

```bash
# 1. Ensure prerequisites are installed
# - Java 11+ (required by Nextflow)
# - Nextflow (https://nextflow.io)
# - Docker (https://docker.com) or Apptainer (https://apptainer.org)

# 2. Clone the repository
git clone https://github.com/dhoconno/nvd.git
cd nvd

# 3. Install dependencies
uv sync                  # Python dependencies only
# or
pixi shell --frozen      # All Conda dependencies

# 4. Create config file manually
mkdir -p ~/.nvd
cat > ~/.nvd/user.config <<'EOF'
params {
    // STAT database paths
    stat_index              = "/path/to/STAT_db/tree_index.20240830.dbs"
    stat_dbss               = "/path/to/STAT_db/tree_filter.20240830.dbss"
    stat_annotation         = "/path/to/STAT_db/tree_filter.20240830.dbss.annotation"
    human_virus_taxlist     = "/path/to/STAT_db/human_viruses_taxlist.txt"
    
    // BLAST database paths
    blast_db                = "/path/to/blast_db"
    blast_db_prefix         = "core_nt"
    
    // GOTTCHA2 database
    gottcha2_db             = "/path/to/gottcha2.fna"
}
EOF

# 5. Verify installation
nvd validate all
```

## Quick Start

### Basic Usage Pattern

```bash
nvd run --samplesheet SAMPLES.csv --experiment-id UNIQUE_ID [OPTIONS]
```

### Minimal Examples

**Human viruses only (STAT+BLAST workflow):**
```bash
nvd run -s samples.csv -e exp001 --tools stat_blast
```

**General metagenomics (GOTTCHA2):**
```bash
nvd run -s samples.csv -e exp002 --tools gottcha
```

**Everything (comprehensive analysis):**
```bash
nvd run -s samples.csv -e exp003 --tools all
```

### Preview Before Running

Use `--dry-run` to see the exact Nextflow command that will be executed:

```bash
nvd run -s samples.csv -e exp004 --tools all --dry-run
```

Output:
```
✓ Auto-detected execution profile: docker
ℹ Using config: /Users/you/.nvd/user.config

Executing command:
nextflow run dhoconno/nvd -profile docker -c /Users/you/.nvd/user.config \
    --samplesheet samples.csv --experiment_id exp004 --tools all

✓ Dry-run mode: command shown above but not executed
```

## Configuration Management

### View Current Configuration

```bash
# Show all configured parameters
nvd config show

# Output:
# Configuration Parameters
# ┏━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
# ┃ Parameter                ┃ Value                          ┃
# ┡━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
# │ blast_db                 │ /data/blast_db                 │
# │ blast_db_prefix          │ core_nt                        │
# │ gottcha2_db              │ /data/gottcha2.fna             │
# │ stat_index               │ /data/STAT_db/tree_index.20240830.dbs │
# └──────────────────────────┴────────────────────────────────┘

# Show config file location
nvd config path
```

### Using Custom Configuration

```bash
# Use a different config file
nvd run -s samples.csv -e exp005 --config /path/to/custom.config
```

### Overriding Configuration Parameters

Even with a config file, you can override specific parameters:

```bash
# Use config defaults but override one database
nvd run -s samples.csv -e exp006 \
    --gottcha2-db /experimental/gottcha2_v3.fna

# Override multiple parameters
nvd run -s samples.csv -e exp007 \
    --blast-db /custom/blast_db \
    --blast-db-prefix custom_nt \
    --cutoff-percent 0.01
```

## Running the Pipeline

### Tool Selection

The `--tools` flag determines which workflows run:

```bash
# Human virus detection only
nvd run -s samples.csv -e exp008 --tools stat_blast
# Aliases: nvd, stat, blast, stast also work

# General taxonomic profiling
nvd run -s samples.csv -e exp009 --tools gottcha

# Both workflows
nvd run -s samples.csv -e exp010 --tools stat_blast,gottcha

# Everything (includes deduplication)
nvd run -s samples.csv -e exp011 --tools all
```

### Profile Selection

The CLI auto-detects your execution environment, but you can override:

```bash
# Auto-detect (default behavior)
nvd run -s samples.csv -e exp012

# Force Docker
nvd run -s samples.csv -e exp013 --profile docker

# Use Apptainer (HPC environments)
nvd run -s samples.csv -e exp014 --profile apptainer

# Local execution (no containers)
nvd run -s samples.csv -e exp015 --profile local
```

### Output Management

```bash
# Custom results directory
nvd run -s samples.csv -e exp016 --results /project/results_2024

# Clean up intermediate files after success
nvd run -s samples.csv -e exp017 --cleanup

# Custom work directory (useful for fast local storage)
nvd run -s samples.csv -e exp018 --work-dir /scratch/nvd_work
```

### Resuming Failed Runs

```bash
# Resume from last checkpoint
nvd run -s samples.csv -e exp019 --resume
```

The `--resume` flag is passed to Nextflow, allowing you to restart from where the pipeline left off.

## Validation Commands

The CLI includes powerful validation tools to catch problems before wasting compute time.

### Validate Dependencies

```bash
nvd validate deps

# Output:
# Checking Dependencies
#
# ✓ Java version "17.0.1" 2021-10-19 LTS
# ✓ Nextflow version 23.10.0 build 5889
# ✓ Docker (running)
#
# ✓ All critical dependencies are installed
```

### Validate Samplesheet

```bash
nvd validate samplesheet samples.csv

# Output:
# Validating Samplesheet: samples.csv
#
# ✓ Header valid: sample_id, srr, platform, fastq1, fastq2
# ✓ Found 24 valid samples
#
# ✓ Samplesheet is valid
```

The validator checks:
- Required columns present
- Platform values are valid (`illumina`, `ont`, or empty)
- Either SRR or FASTQ files provided
- No duplicate sample IDs
- No empty sample_id fields

### Validate Database Paths

```bash
nvd validate databases

# Output:
# Validating Database Paths
#
# ℹ Using config: /Users/you/.nvd/user.config
#
# ✓ STAT index file             - /data/STAT_db/tree_index.20240830.dbs
# ✓ STAT dbss file               - /data/STAT_db/tree_filter.20240830.dbss
# ✓ STAT annotation file         - /data/STAT_db/tree_filter.20240830.dbss.annotation
# ✓ Human virus taxlist          - /data/STAT_db/human_viruses_taxlist.txt
# ✓ BLAST database directory     - /data/blast_db
# ✓ GOTTCHA2 database file       - /data/gottcha2.fna
#
# ✓ All database paths are valid
```

### Validate Parameter Coverage (CI/Development)

This validates that the CLI wrapper covers all configurable Nextflow parameters:

```bash
nvd validate params

# Output:
# Validating Parameter Coverage
#
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# COVERAGE ANALYSIS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# ✓ Covered Parameters (23):
#   • samplesheet → --samplesheet, -s
#   • experiment_id → --experiment-id, -e
#   • tools → --tools, -t
#   ... and 20 more
#
# Coverage: 23/23 (100.0%)
# ✓ Status: Full coverage
```

Useful for developers to ensure the CLI stays in sync with Nextflow config.

### Run All Validations

```bash
# Run all validations at once
nvd validate all

# Optionally include samplesheet
nvd validate all --samplesheet samples.csv
```

## Advanced Usage

### Analysis Parameter Tuning

#### Sensitivity Adjustment

```bash
# High sensitivity (detect rare organisms)
nvd run -s samples.csv -e exp020 \
    --cutoff-percent 0.0001 \
    --tax-stringency 0.5 \
    --entropy 0.7 \
    --min-gottcha-reads 50

# High specificity (minimize false positives)
nvd run -s samples.csv -e exp021 \
    --cutoff-percent 0.01 \
    --tax-stringency 0.9 \
    --entropy 0.95 \
    --min-gottcha-reads 1000
```

#### BLAST Configuration

```bash
# More exhaustive BLAST analysis
nvd run -s samples.csv -e exp022 \
    --max-blast-targets 200 \
    --blast-retention-count 10

# Faster, less exhaustive BLAST
nvd run -s samples.csv -e exp023 \
    --max-blast-targets 50 \
    --blast-retention-count 3
```

#### Read Quality Control

```bash
# Illumina optimization
nvd run -s samples.csv -e exp024 \
    --min-consecutive-bases 150 \
    --qtrim 't'

# Nanopore optimization
nvd run -s samples.csv -e exp025 \
    --min-consecutive-bases 300 \
    --qtrim 'f'
```

### SRA Dataset Processing

```bash
# Limit concurrent downloads (respect NCBI rate limits)
nvd run -s sra_samples.csv -e exp026 \
    --max-concurrent-downloads 3

# Faster download for large batches (use with caution)
nvd run -s sra_samples.csv -e exp027 \
    --max-concurrent-downloads 8
```

### Taxonomic Classification Options

```bash
# Include child taxa in classification
nvd run -s samples.csv -e exp028 \
    --include-children

# Exclude child taxa (stricter classification)
nvd run -s samples.csv -e exp029 \
    --no-include-children
```

### Privacy and Data Sharing

```bash
# Scrub human reads before sharing
nvd run -s samples.csv -e exp030 \
    --human-read-scrub /path/to/human_filter_db \
    --tools clumpify
```

### LabKey Integration

```bash
# Enable LabKey LIMS upload
nvd run -s samples.csv -e exp031 \
    --labkey

# LabKey configuration is in ~/.nvd/user.config
# or can be specified via Nextflow config file
```

## Practical Examples

### Example 1: Wastewater Surveillance

**Scenario**: Weekly wastewater monitoring for human viruses

```bash
# Week 1
nvd run \
    --samplesheet wastewater_week1.csv \
    --experiment-id WW_2024_W01 \
    --tools stat_blast,gottcha \
    --profile docker \
    --results /project/wastewater_2024 \
    --cleanup

# Week 2 (resume if interrupted)
nvd run \
    -s wastewater_week2.csv \
    -e WW_2024_W02 \
    --tools stat_blast,gottcha \
    -p docker \
    --results /project/wastewater_2024 \
    --resume
```

### Example 2: Clinical Sample Screening

**Scenario**: High-specificity viral detection in patient samples

```bash
nvd run \
    --samplesheet clinical_batch_5.csv \
    --experiment-id CLINIC_2024_B5 \
    --tools stat_blast \
    --profile docker \
    --cutoff-percent 0.01 \
    --tax-stringency 0.9 \
    --entropy 0.95 \
    --min-consecutive-bases 300 \
    --results /secure/clinical_results \
    --cleanup
```

### Example 3: Mixed Platform Dataset

**Scenario**: Combining Illumina and Nanopore data from air samples

```bash
# First validate the samplesheet
nvd validate samplesheet mixed_platform_samples.csv

# Then run with balanced parameters
nvd run \
    -s mixed_platform_samples.csv \
    -e AIR_SAMPLING_2024_Q1 \
    --tools all \
    -p apptainer \
    --min-gottcha-reads 250 \
    --min-consecutive-bases 200 \
    --results /project/air_surveillance
```

### Example 4: Rapid Outbreak Investigation

**Scenario**: Emergency viral detection with maximum sensitivity

```bash
# High sensitivity, preserve all intermediate files
nvd run \
    --samplesheet outbreak_samples.csv \
    --experiment-id OUTBREAK_2024_001 \
    --tools stat_blast \
    --profile docker \
    --cutoff-percent 0.0001 \
    --tax-stringency 0.5 \
    --max-blast-targets 200 \
    --max-concurrent-downloads 8 \
    --results /urgent/outbreak_analysis \
    # Don't use --cleanup, keep all intermediate files for review
```

### Example 5: HPC Cluster Execution

**Scenario**: Large-scale analysis on institutional HPC

```bash
# Submit as batch job
#!/bin/bash
#SBATCH --job-name=nvd_analysis
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G

module load java/17
module load nextflow

nvd run \
    --samplesheet /shared/data/large_dataset.csv \
    --experiment-id HPC_2024_BATCH_03 \
    --tools all \
    --profile apptainer \
    --work-dir /scratch/$USER/nvd_work \
    --results /shared/results/batch_03 \
    --cleanup
```

### Example 6: Development and Testing

**Scenario**: Test pipeline with single sample before batch processing

```bash
# Create test samplesheet with one sample
head -2 full_dataset.csv > test_sample.csv

# Validate
nvd validate all -s test_sample.csv

# Test run with dry-run
nvd run -s test_sample.csv -e TEST_001 --dry-run

# Actual test run
nvd run -s test_sample.csv -e TEST_001

# If successful, run full dataset
nvd run -s full_dataset.csv -e PROD_2024_001 --resume
```

## Troubleshooting

### Common Issues

#### "Nextflow not found"

```bash
# Check if Nextflow is installed
nvd validate deps

# Install Nextflow if missing
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

#### "No config file found"

```bash
# Check config location
nvd config path

# Output shows:
# Default config location: /Users/you/.nvd/user.config
# ⚠ Config file does not exist

# Solution: Run install.sh or create manually
./install.sh
# or
mkdir -p ~/.nvd
# ... create user.config manually
```

#### "Docker daemon not running"

```bash
# Check dependencies
nvd validate deps

# Output:
# ⚠ Docker (daemon not running)

# macOS: Start Docker Desktop
# Linux: sudo systemctl start docker
```

#### "Database path not found"

```bash
# Check what's configured
nvd validate databases

# Output:
# ✗ BLAST database directory - NOT FOUND: /wrong/path/blast_db

# Fix in config file
nvd config path  # Shows location
# Edit ~/.nvd/user.config with correct paths

# Or override on command line
nvd run -s samples.csv -e exp032 --blast-db /correct/path/blast_db
```

### Validation Workflow

Before running a large batch, always validate:

```bash
# 1. Check system setup
nvd validate deps

# 2. Check databases
nvd validate databases

# 3. Check samplesheet
nvd validate samplesheet my_samples.csv

# 4. Dry-run to preview command
nvd run -s my_samples.csv -e TEST --dry-run

# 5. Test with single sample
head -2 my_samples.csv > test.csv
nvd run -s test.csv -e TEST_SINGLE

# 6. Run full dataset if test passes
nvd run -s my_samples.csv -e PRODUCTION_RUN
```

### Getting Help

```bash
# General help
nvd --help

# Command-specific help
nvd run --help
nvd validate --help
nvd config --help

# Version information
nvd version
```

## Complete Parameter Reference

### Required Parameters

| Parameter | Short | Description |
|-----------|-------|-------------|
| `--samplesheet` | `-s` | CSV file with sample information |
| `--experiment-id` | `-e` | Unique identifier for this run |

### Workflow Selection

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--tools` | `-t` | Workflows to run: `stat_blast`, `gottcha`, `all` | `all` |

### Execution Options

| Parameter | Short | Description | Default |
|-----------|-------|-------------|---------|
| `--profile` | `-p` | Execution profile: `docker`, `apptainer`, `local` | Auto-detect |
| `--config` | `-c` | Custom config file | `~/.nvd/user.config` |
| `--results` | `-r` | Results output directory | `./results` |
| `--work-dir` | `-w` | Nextflow work directory | `./work` |
| `--resume` | | Resume from checkpoint | `false` |
| `--cleanup` | | Remove work dir after success | `false` |
| `--dry-run` | | Show command without executing | `false` |

### Database Paths

| Parameter | Description |
|-----------|-------------|
| `--gottcha2-db` | GOTTCHA2 database file |
| `--blast-db` | BLAST database directory |
| `--blast-db-prefix` | BLAST database prefix name |
| `--stat-index` | STAT index file |
| `--stat-dbss` | STAT dbss file |
| `--stat-annotation` | STAT annotation file |
| `--human-virus-taxlist` | Human virus taxonomy list |

### Analysis Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--cutoff-percent` | Abundance cutoff threshold | `0.001` |
| `--tax-stringency` | Taxonomy confidence level | `0.7` |
| `--entropy` | Sequence entropy threshold | `0.9` |
| `--min-consecutive-bases` | Minimum consecutive bases | `200` |
| `--min-gottcha-reads` | Minimum reads for GOTTCHA2 | `250` |
| `--max-blast-targets` | Maximum BLAST targets to consider | `100` |
| `--blast-retention-count` | Top BLAST hits to retain | `5` |
| `--qtrim` | Quality trimming mode | `'t'` |
| `--include-children` / `--no-include-children` | Include taxonomic children | `true` |

### Data Management

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--max-concurrent-downloads` | Max concurrent SRA downloads | `3` |
| `--human-read-scrub` | Human read scrubbing database path | None |

### Integration

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--labkey` | Enable LabKey integration | `false` |

## See Also

- [Installation Guide](./INSTALLATION.md) - Complete installation instructions with install.sh
- [Direct Nextflow Examples](./example_commands.md) - Traditional Nextflow command examples
- [Contributor Guide](./contributor_guide.md) - Development guidelines
- [Main README](../README.md) - Project overview

## Tips and Best Practices

1. **Always validate before large runs** - Use `nvd validate all` to catch config issues early
2. **Use dry-run** - Preview commands with `--dry-run` before committing to long analyses
3. **Start with test samples** - Validate pipeline behavior on small datasets first
4. **Leverage config files** - Store database paths in `~/.nvd/user.config` to avoid repetition
5. **Enable cleanup for production** - Use `--cleanup` to save disk space on routine runs
6. **Respect rate limits** - Keep `--max-concurrent-downloads` low (≤3) for SRA downloads
7. **Use resume** - Always add `--resume` when restarting failed runs to save time
8. **Match sensitivity to purpose** - High sensitivity for surveillance, high specificity for diagnostics
9. **Document your parameters** - Use dry-run output to record exact parameters for reproducibility
10. **Check parameter coverage** - Developers: Run `nvd validate params` to ensure CLI stays updated
