# NVD2 Example Run Commands

This guide provides comprehensive examples of how to run NVD2 for different use
cases. Each example focuses on specific scenarios that researchers commonly
encounter, with explanations of when and why you'd use each configuration.

For basic setup instructions, see the main [README.md](../README.md) file.

## Tool Selection Note

**STAT+BLAST Workflow Aliases**: The human virus detection workflow (STAT + two-phase BLAST) 
can be invoked using any of these tool names: `stat_blast`, `nvd`, `stat`, `blast`, or `stast`. 
They all trigger the same workflow. The legacy name `nvd` is maintained for backward compatibility.

Examples in this guide use various aliases to demonstrate flexibility, but they're all equivalent.

## Table of Contents

- [Quick Start Examples](#quick-start-examples)
- [Deployment Profiles](#deployment-profiles)
- [Platform-Specific Configurations](#platform-specific-configurations)
- [Quality and Sensitivity Adjustments](#quality-and-sensitivity-adjustments)
- [Data Sources and Input Types](#data-sources-and-input-types)
- [Output and Results Management](#output-and-results-management)
- [Privacy and Compliance](#privacy-and-compliance)
- [Enterprise Integration (LabKey LIMS)](#enterprise-integration-labkey-lims)
- [Specialized Use Cases](#specialized-use-cases)
- [Troubleshooting and Debugging](#troubleshooting-and-debugging)
- [Configuration Files and Parameters](#configuration-files-and-parameters)
- [Useful Nextflow Options](#useful-nextflow-options)

## Quick Start Examples

### Minimal STAT+BLAST run (human viruses only)

Use this when you're specifically interested in detecting human viruses and want
the fastest, most focused analysis:

**Note**: Tool aliases `stat_blast`, `nvd`, `stat`, `blast`, and `stast` all invoke 
the same workflow. Use whichever name is most intuitive for your team.

```bash
nextflow run dhoconno/nvd \
    --tools stat_blast \
    --samplesheet assets/example_samplesheet.csv \
    --gottcha2_db db/PP819512_gottcha2.fasta \
    --blast_db db \
    --blast_db_prefix PP819512-nt \
    --stat_index db/tree_index.dense.dbs \
    --stat_dbss db/tree_index.dense.dbss \
    --stat_annotation db/tree_index.dense.dbss.annotation \
    --human_virus_taxlist db/human_viruses_taxlist.txt \
    --experiment_id 12345
```

**When to use**: Wastewater surveillance, clinical samples where you
specifically need human virus detection, outbreak investigations focusing on
viral pathogens.

### Minimal GOTTCHA2 run (general metagenomics)

Use this for broad taxonomic profiling when you want to see everything in your
sample:

```bash
nextflow run dhoconno/nvd \
    --tools gottcha \
    --samplesheet assets/example_samplesheet.csv \
    --gottcha2_db db/PP819512_gottcha2.fasta \
    --experiment_id 12346
```

**When to use**: Environmental samples, microbiome studies, exploratory analysis
of unknown samples, quality control to see overall sample composition.

### Kitchen sink (all workflows)

Run everything for comprehensive analysis - human virus detection, general
taxonomic profiling, and data deduplication:

```bash
nextflow run dhoconno/nvd \
    --tools all \
    --samplesheet assets/example_samplesheet.csv \
    --gottcha2_db db/PP819512_gottcha2.fasta \
    --blast_db db \
    --blast_db_prefix PP819512-nt \
    --stat_index db/tree_index.dense.dbs \
    --stat_dbss db/tree_index.dense.dbss \
    --stat_annotation db/tree_index.dense.dbss.annotation \
    --human_virus_taxlist db/human_viruses_taxlist.txt \
    --experiment_id 12347
```

**When to use**: Research projects where you want comprehensive analysis, when
compute resources aren't limiting, for samples where you're unsure what you
might find.

## Deployment Profiles

### Local execution (no containers)

Best for development or when you have all tools installed locally:

```bash
nextflow run dhoconno/nvd \
    -profile local \
    --tools nvd \
    --samplesheet my_samples.csv \
    --blast_db /path/to/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /path/to/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /path/to/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /path/to/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /path/to/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12348
```

**When to use**: Development environments, when you've manually installed all
dependencies, debugging tool-specific issues.

### Docker execution (recommended for local development)

Provides consistent environment across different systems:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools all \
    --samplesheet my_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12349
```

**When to use**: Local workstations, development environments, when you want
reproducible results without manual tool installation.

### Apptainer/Singularity execution (HPC environments)

Ideal for high-performance computing clusters that don't support Docker:

```bash
nextflow run dhoconno/nvd \
    -profile apptainer \
    --tools nvd,gottcha \
    --samplesheet large_dataset.csv \
    --gottcha2_db /shared/gottcha2/gottcha_db.species.fna \
    --blast_db /shared/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /shared/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /shared/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /shared/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /shared/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12350
```

**When to use**: University HPC clusters, government computing facilities, any
environment where Docker is restricted but Singularity/Apptainer is available.

### CHTC HPC cluster execution

Configured specifically for Center for High Throughput Computing environments:

```bash
nextflow run dhoconno/nvd \
    -profile chtc_hpc \
    --tools all \
    --samplesheet wastewater_samples.csv \
    --gottcha2_db /staging/gottcha2/gottcha_db.species.fna \
    --blast_db /staging/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /staging/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /staging/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /staging/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /staging/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12351 \
    --max_concurrent_downloads 5
```

**When to use**: CHTC infrastructure, SLURM-based clusters with shared storage,
when you need high parallelization for large datasets.

## Platform-Specific Configurations

### Illumina-heavy dataset optimization

Optimized parameters for high-quality short-read Illumina data:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd,gottcha \
    --samplesheet illumina_wastewater.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --min_gottcha_reads 500 \
    --min_consecutive_bases 150 \
    --experiment_id 12352
```

**Key parameters**:

- `--min_gottcha_reads 500`: Higher threshold for reliable classification on
  short reads
- `--min_consecutive_bases 150`: Optimized for typical Illumina read lengths

**When to use**: Wastewater surveillance with Illumina sequencing, clinical
samples, any high-quality short-read dataset.

### Oxford Nanopore dataset optimization

Optimized for longer, more error-prone Nanopore reads:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd,gottcha \
    --samplesheet nanopore_air_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --min_gottcha_reads 100 \
    --min_consecutive_bases 300 \
    --qtrim 'f' \
    --experiment_id 12353
```

**Key parameters**:

- `--min_gottcha_reads 100`: Lower threshold appropriate for longer reads
- `--min_consecutive_bases 300`: Takes advantage of longer Nanopore reads
- `--qtrim 'f'`: Disables quality trimming (less critical for Nanopore)

**When to use**: Air sampling studies, environmental surveillance with portable
sequencers, when you need longer reads for better assembly.

### Mixed platform dataset

Balanced parameters for datasets containing both Illumina and Nanopore data:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools all \
    --samplesheet mixed_platform_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --min_gottcha_reads 250 \
    --experiment_id 12354
```

**When to use**: Comparative studies using multiple sequencing platforms, when
combining historical and new data from different instruments.

## Quality and Sensitivity Adjustments

### High sensitivity mode (lower thresholds)

Maximize detection of low-abundance organisms, accepting higher false positive
risk:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet low_biomass_samples.csv \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --cutoff_percent 0.0001 \
    --tax_stringency 0.5 \
    --entropy 0.7 \
    --min_consecutive_bases 100 \
    --min_gottcha_reads 50 \
    --experiment_id 12355
```

**Key parameters**:

- `--cutoff_percent 0.0001`: Very low threshold for rare organisms
- `--tax_stringency 0.5`: Lower confidence requirements
- `--entropy 0.7`: Allow more repetitive sequences
- `--min_consecutive_bases 100`: Shorter minimum sequences

**When to use**: Low-biomass samples, early outbreak detection, environmental
surveillance where you don't want to miss anything.

### High specificity mode (higher thresholds)

Minimize false positives for high-confidence results:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd,gottcha \
    --samplesheet high_biomass_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --cutoff_percent 0.01 \
    --tax_stringency 0.9 \
    --entropy 0.95 \
    --min_consecutive_bases 500 \
    --min_gottcha_reads 1000 \
    --experiment_id 12356
```

**Key parameters**:

- `--cutoff_percent 0.01`: Higher threshold for confident detection
- `--tax_stringency 0.9`: Require high confidence for classification
- `--entropy 0.95`: Exclude repetitive/low-complexity regions
- `--min_consecutive_bases 500`: Require longer, more reliable sequences

**When to use**: Clinical diagnostics, regulatory reporting, when false
positives have serious consequences.

### Custom virus family targeting

Focus analysis on specific virus families of interest:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet respiratory_samples.csv \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12357 \
    -params-file custom_virus_families.json
```

Create `custom_virus_families.json` for respiratory virus surveillance:

```json
{
  "human_virus_families": [
    "Coronaviridae",
    "Orthomyxoviridae",
    "Paramyxoviridae",
    "Pneumoviridae",
    "Picornaviridae"
  ]
}
```

**When to use**: Targeted surveillance (respiratory, enteric, etc.), when you
want to focus computational resources on specific virus types.

## Data Sources and Input Types

### Local FASTQ files only

Process locally stored sequencing data:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools gottcha \
    --samplesheet local_fastq_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --experiment_id 12358
```

**Sample `local_fastq_samples.csv`:**

```csv
sample_id,srr,platform,fastq1,fastq2
sample1,,illumina,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,,ont,/path/to/sample2.fastq.gz,
```

**When to use**: Processing local sequencing runs, private data that cannot be
shared publicly, when you have fast local storage.

### SRA datasets only

Process public datasets from NCBI's Sequence Read Archive:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet sra_samples.csv \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --max_concurrent_downloads 2 \
    --experiment_id 12359
```

**Sample `sra_samples.csv`:**

```csv
sample_id,srr,platform,fastq1,fastq2
wastewater_1,SRR12345678,,,
wastewater_2,SRR12345679,,,
wastewater_3,SRR12345680,,,
```

**Key parameters**:

- `--max_concurrent_downloads 2`: Respect NCBI rate limits

**When to use**: Meta-analyses of public data, comparative studies, when
reproducing published results.

### Mixed local and SRA datasets

Combine your own data with public datasets for comparison:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools all \
    --samplesheet mixed_sources.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --max_concurrent_downloads 3 \
    --experiment_id 12360
```

**Sample `mixed_sources.csv`:**

```csv
sample_id,srr,platform,fastq1,fastq2
local_sample1,,illumina,/data/local1_R1.fastq.gz,/data/local1_R2.fastq.gz
sra_sample1,SRR12345678,,,
local_sample2,,ont,/data/local2.fastq.gz,
sra_sample2,SRR12345679,,,
```

**When to use**: Benchmarking your samples against public data, temporal studies
combining historical and new data.

## Output and Results Management

### Custom results directory

Organize results in a specific location:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet my_samples.csv \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --results /project/nvd_analysis_2024 \
    --experiment_id 12361
```

**When to use**: Project-specific organization, shared storage systems, when you
need results in a particular location for downstream analysis.

### Automatic cleanup enabled

Save disk space by removing intermediate files after successful completion:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools gottcha \
    --samplesheet my_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --cleanup true \
    --experiment_id 12362
```

**When to use**: Production runs with limited disk space, when you only need
final results, routine surveillance processing.

### Resume failed run with work directory preservation

Restart from the last successful checkpoint:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools all \
    --samplesheet my_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --cleanup false \
    --experiment_id 12363 \
    -resume
```

**When to use**: After system failures, when debugging pipeline issues, for
long-running analyses that might be interrupted.

## Privacy and Compliance

### Human read scrubbing for public sharing

Remove human genetic material before sharing data publicly:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools clumpify \
    --samplesheet sensitive_samples.csv \
    --human_read_scrub true \
    --human_filter_db /data/human_filter_db \
    --experiment_id 12364
```

**When to use**: Before uploading to SRA, sharing with external collaborators,
when working with clinical samples that may contain human DNA.

### Deduplication and compression for storage

Optimize data for long-term storage while preserving analysis capability:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools clumpify,gottcha \
    --samplesheet large_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --human_read_scrub false \
    --experiment_id 12365
```

**When to use**: Archiving large datasets, when storage costs are a concern, for
datasets you may want to reanalyze later.

## Enterprise Integration (LabKey LIMS)

### Basic LabKey integration

Upload results directly to LabKey Laboratory Information Management System:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet my_samples.csv \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12366 \
    -params-file labkey_config.json
```

**Example `labkey_config.json`:**

```json
{
  "labkey": true,
  "labkey_server": "https://your-labkey-server.org",
  "labkey_project_name": "YourProject",
  "labkey_schema": "lists",
  "labkey_webdav": "https://your-labkey-server.org/_webdav/YourProject/@files/",
  "labkey_blast_meta_hits_list": "nvd_blast_hits",
  "labkey_blast_fasta_list": "nvd_fasta_sequences"
}
```

**When to use**: Laboratory workflows, regulatory compliance, when you need
structured data management and audit trails.

### Full enterprise integration

Complete integration with all workflows and LabKey systems:

```bash
nextflow run dhoconno/nvd \
    -profile apptainer \
    --tools all \
    --samplesheet enterprise_samples.csv \
    --gottcha2_db /shared/gottcha2/gottcha_db.species.fna \
    --blast_db /shared/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /shared/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /shared/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /shared/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /shared/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12368 \
    -params-file full_enterprise_config.json
```

**When to use**: Production surveillance systems, institutional deployments,
when you need comprehensive data management.

## Specialized Use Cases

### Wastewater surveillance pipeline

Optimized for routine wastewater monitoring programs:

```bash
nextflow run dhoconno/nvd \
    -profile chtc_hpc \
    --tools nvd,gottcha \
    --samplesheet wastewater_surveillance.csv \
    --gottcha2_db /staging/gottcha2/gottcha_db.species.fna \
    --blast_db /staging/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /staging/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /staging/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /staging/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /staging/STAT_db/human_viruses_taxlist.txt \
    --cutoff_percent 0.001 \
    --tax_stringency 0.7 \
    --min_gottcha_reads 100 \
    --max_concurrent_downloads 10 \
    --experiment_id 20240001 \
    -params-file wastewater_labkey.json
```

**Key features**: Balanced sensitivity, high throughput, enterprise integration
for public health reporting.

### Clinical sample screening

High-specificity analysis for diagnostic applications:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet clinical_samples.csv \
    --blast_db /data/clinical_blast_db \
    --blast_db_prefix clinical_nt \
    --stat_index /data/clinical_STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/clinical_STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/clinical_STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/clinical_STAT_db/human_viruses_taxlist.txt \
    --cutoff_percent 0.01 \
    --tax_stringency 0.9 \
    --entropy 0.95 \
    --min_consecutive_bases 300 \
    --experiment_id 20240002
```

**Key features**: High specificity parameters, clinical-grade databases, minimal
false positives.

### Rapid outbreak investigation

Fast turnaround for emergency response:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet outbreak_samples.csv \
    --blast_db /data/outbreak_blast_db \
    --blast_db_prefix outbreak_nt \
    --stat_index /data/outbreak_STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/outbreak_STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/outbreak_STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/outbreak_STAT_db/human_viruses_taxlist.txt \
    --cutoff_percent 0.0001 \
    --tax_stringency 0.5 \
    --max_concurrent_downloads 8 \
    --cleanup false \
    --experiment_id 20240004 \
    -params-file outbreak_labkey.json
```

**Key features**: High sensitivity, fast processing, preserved intermediate
files for follow-up analysis.

## Troubleshooting and Debugging

### Debug mode with comprehensive logging

Generate detailed reports for troubleshooting:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet debug_samples.csv \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --cleanup false \
    --experiment_id 99999 \
    -with-report debug_report.html \
    -with-trace debug_trace.txt \
    -with-timeline debug_timeline.html \
    -with-dag debug_dag.png
```

**Generates**:

- `debug_report.html`: Execution summary and resource usage
- `debug_trace.txt`: Detailed process execution log
- `debug_timeline.html`: Visual timeline of pipeline execution
- `debug_dag.png`: Workflow diagram

**When to use**: When troubleshooting failures, optimizing performance,
understanding resource requirements.

### Single sample test run

Test pipeline with minimal data for quick validation:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools nvd \
    --samplesheet single_sample.csv \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 99998
```

**Single sample CSV example:**

```csv
sample_id,srr,platform,fastq1,fastq2
test_sample,,illumina,test_R1.fastq.gz,test_R2.fastq.gz
```

**When to use**: Testing new configurations, validating database setups, quick
pipeline checks.

### Dry run (validate configuration without execution)

Check configuration and workflow structure without running analysis:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools all \
    --samplesheet my_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 99997 \
    -preview
```

**When to use**: Validating complex configurations, checking parameter
compatibility, planning resource allocation.

## Configuration Files and Parameters

### Using a custom nextflow.config file

Create reusable configurations for your environment:

```bash
nextflow run dhoconno/nvd \
    -c custom_nvd.config \
    --tools all \
    --samplesheet my_samples.csv \
    --experiment_id 12369
```

**Example `custom_nvd.config`:**

```groovy
params {
    gottcha2_db = "/shared/databases/gottcha2/gottcha_db.species.fna"
    blast_db = "/shared/databases/blast_db"
    blast_db_prefix = "core_nt"
    stat_index = "/shared/databases/STAT_db/tree_index.20240830.dbs"
    stat_dbss = "/shared/databases/STAT_db/tree_filter.20240830.dbss"
    stat_annotation = "/shared/databases/STAT_db/tree_filter.20240830.dbss.annotation"
    human_virus_taxlist = "/shared/databases/STAT_db/human_viruses_taxlist.txt"
    results = "/project/results"
    labkey = true
    labkey_server = "https://company-labkey.org"
    labkey_project_name = "ProjectName"
}

profiles {
    company_hpc {
        apptainer.enabled = true
        process.executor = 'slurm'
        process.queue = 'high-memory'
        process.cpus = 32
        process.memory = 128.GB
    }
}
```

### Using parameters file

Store all parameters in a JSON file for reproducibility:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    -params-file complete_params.json
```

**Example `complete_params.json`:**

```json
{
  "tools": "all",
  "samplesheet": "comprehensive_samples.csv",
  "experiment_id": 12370,
  "gottcha2_db": "/data/gottcha2/gottcha_db.species.fna",
  "blast_db": "/data/blast_db",
  "blast_db_prefix": "core_nt",
  "stat_index": "/data/STAT_db/tree_index.20240830.dbs",
  "stat_dbss": "/data/STAT_db/tree_filter.20240830.dbss",
  "stat_annotation": "/data/STAT_db/tree_filter.20240830.dbss.annotation",
  "human_virus_taxlist": "/data/STAT_db/human_viruses_taxlist.txt",
  "cutoff_percent": 0.001,
  "tax_stringency": 0.8,
  "entropy": 0.9,
  "min_consecutive_bases": 200,
  "min_gottcha_reads": 250,
  "max_concurrent_downloads": 5,
  "results": "/project/nvd_results",
  "cleanup": true,
  "labkey": true,
  "labkey_server": "https://labkey.example.org",
  "labkey_project_name": "MetagenomicsSurveillance"
}
```

**When to use**: Reproducible research, sharing configurations, complex
enterprise setups.

## Useful Nextflow Options

### Limit resource usage

Control computational resource consumption:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools gottcha \
    --samplesheet small_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --experiment_id 12371 \
    -process.cpus 4 \
    -process.memory 8.GB
```

**When to use**: Shared systems, resource-constrained environments, development
work.

### Use work directory on different filesystem

Improve performance by using fast storage for temporary files:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools all \
    --samplesheet my_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12372 \
    -work-dir /scratch/nvd_work
```

**When to use**: HPC environments with fast local storage, when your home
directory has limited space.

### Generate comprehensive execution reports

Create detailed analytics for pipeline optimization:

```bash
nextflow run dhoconno/nvd \
    -profile docker \
    --tools all \
    --samplesheet my_samples.csv \
    --gottcha2_db /data/gottcha2/gottcha_db.species.fna \
    --blast_db /data/blast_db \
    --blast_db_prefix core_nt \
    --stat_index /data/STAT_db/tree_index.20240830.dbs \
    --stat_dbss /data/STAT_db/tree_filter.20240830.dbss \
    --stat_annotation /data/STAT_db/tree_filter.20240830.dbss.annotation \
    --human_virus_taxlist /data/STAT_db/human_viruses_taxlist.txt \
    --experiment_id 12373 \
    -with-report execution_report.html \
    -with-trace execution_trace.txt \
    -with-timeline execution_timeline.html \
    -with-dag workflow_dag.png
```

**When to use**: Performance optimization, resource planning, documentation for
publications.

## Important Notes

1. **Unique experiment IDs**: Always ensure your `experiment_id` is unique for
   each run
2. **Absolute paths**: Database paths should be absolute paths to avoid
   confusion
3. **Samplesheet format**: The CSV format is critical - see
   `assets/example_samplesheet.csv`
4. **LabKey integration**: Ensure the `nvd2` secret is properly configured
   before using LabKey features
5. **Resume capability**: Use `-resume` to restart failed workflows from the
   last successful checkpoint
6. **Container considerations**: Use `-profile apptainer` for HPC environments
   without Docker
7. **Disk space**: Consider using `--cleanup true` when processing large
   datasets
8. **NCBI limits**: Respect rate limits with appropriate
   `max_concurrent_downloads` for SRA data

For additional help and troubleshooting, see the main [README.md](../README.md)
file and [contributor guide](contributor_guide.md).
