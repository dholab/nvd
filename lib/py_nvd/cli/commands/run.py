"""
Run command for the NVD CLI.

Commands:
    nvd run  - Run the NVD pipeline
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import typer

from py_nvd.cli.provenance import nextflow_environment
from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    PANEL_ANALYSIS,
    PANEL_CORE,
    PANEL_DATABASES,
    PANEL_LABKEY,
    PANEL_NOTIFICATIONS,
    PANEL_PREPROCESSING,
    RESUME_FILE,
    auto_detect_profile,
    check_command_exists,
    console,
    error,
    find_config_file,
    format_command_for_display,
    get_default_profile,
    get_pipeline_root,
    info,
    success,
    warning,
)
from py_nvd.models import NvdParams
from py_nvd.params import load_params_file
from py_nvd.paths import ENV_VAR_TAXONOMY
from py_nvd.presets import get_preset_store


def run(
    # Complexity/boolean/B008 warnings are acceptable for CLI with many options
    ctx: typer.Context,
    # -------------------------------------------------------------------------
    # Core Options
    # -------------------------------------------------------------------------
    samplesheet: Path | None = typer.Option(
        None,
        "--samplesheet",
        "-s",
        help="Path to samplesheet CSV (required; can come from --preset or --params-file)",
        file_okay=True,
        dir_okay=False,
        rich_help_panel=PANEL_CORE,
    ),
    experiment_id: str | None = typer.Option(
        None,
        "--experiment-id",
        "-e",
        help="Unique experiment identifier (required for LabKey uploads)",
        rich_help_panel=PANEL_CORE,
    ),
    results: Path | None = typer.Option(
        None,
        "--results",
        "-r",
        help="Results directory (default: ./results)",
        rich_help_panel=PANEL_CORE,
    ),
    profile: str | None = typer.Option(
        None,
        "--profile",
        "-p",
        help="Execution profile: docker, apptainer, local (auto-detect if not set)",
        rich_help_panel=PANEL_CORE,
    ),
    config: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Custom config file (default: NVD_CONFIG_DIR/user.config or NVD_CONFIG)",
        exists=True,
        rich_help_panel=PANEL_CORE,
    ),
    resume: bool | None = typer.Option(
        None,
        "--resume/--no-resume",
        help="Resume from checkpoint",
        rich_help_panel=PANEL_CORE,
    ),
    cleanup: bool | None = typer.Option(
        None,
        "--cleanup/--no-cleanup",
        help="Clean work directory after success",
        rich_help_panel=PANEL_CORE,
    ),
    work_dir: Path | None = typer.Option(
        None,
        "--work-dir",
        "-w",
        help="Nextflow work directory",
        rich_help_panel=PANEL_CORE,
    ),
    taxonomy_dir: Path | None = typer.Option(
        None,
        "--taxonomy-dir",
        help="Explicit taxonomy database directory",
        rich_help_panel=PANEL_CORE,
    ),
    taxonomy_mode: str | None = typer.Option(
        None,
        "--taxonomy-mode",
        help="Pipeline taxonomy mode: read_only or missing",
        rich_help_panel=PANEL_CORE,
    ),
    taxonomy_refresh: str | None = typer.Option(
        None,
        "--taxonomy-refresh",
        help="Taxonomy preflight refresh policy: missing, stale, or force",
        rich_help_panel=PANEL_CORE,
    ),
    taxonomy_max_age_days: int | None = typer.Option(
        None,
        "--taxonomy-max-age-days",
        help="Freshness threshold for taxonomy refreshes and warnings",
        rich_help_panel=PANEL_CORE,
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        "-P",
        help="Use a registered preset (CLI args override preset values)",
        rich_help_panel=PANEL_CORE,
    ),
    params_file: Path | None = typer.Option(
        None,
        "--params-file",
        "-f",
        help="JSON/YAML params file (merged with CLI args; CLI takes precedence)",
        exists=True,
        file_okay=True,
        dir_okay=False,
        rich_help_panel=PANEL_CORE,
    ),
    experimental: bool | None = typer.Option(
        None,
        "--experimental",
        help="Enable experimental release-candidate features",
        rich_help_panel=PANEL_CORE,
    ),
    skip_assembly: bool | None = typer.Option(
        None,
        "--skip-assembly",
        help="Skip SPAdes assembly and downstream contig classification",
        rich_help_panel=PANEL_CORE,
    ),
    skip_blast: bool | None = typer.Option(
        None,
        "--skip-blast",
        help="Skip MEGABLAST and BLASTN contig search.",
        rich_help_panel=PANEL_CORE,
    ),
    # -------------------------------------------------------------------------
    # Reference Paths
    # -------------------------------------------------------------------------
    blast_db: Path | None = typer.Option(
        None,
        "--blast-db",
        help="Override BLAST database directory",
        rich_help_panel=PANEL_DATABASES,
    ),
    blast_db_prefix: str | None = typer.Option(
        None,
        "--blast-db-prefix",
        help="Override BLAST database prefix",
        rich_help_panel=PANEL_DATABASES,
    ),
    virus_index: Path | None = typer.Option(
        None,
        "--virus-index",
        help="Path to prebuilt vertebrate-infecting virus deacon index (.idx file)",
        rich_help_panel=PANEL_DATABASES,
    ),
    virus_index_url: str | None = typer.Option(
        None,
        "--virus-index-url",
        help="URL to download a prebuilt vertebrate-infecting virus deacon index",
        rich_help_panel=PANEL_DATABASES,
    ),
    virus_reference_fasta: Path | None = typer.Option(
        None,
        "--virus-reference-fasta",
        help="Custom vertebrate-infecting virus FASTA for building an enrichment index",
        rich_help_panel=PANEL_DATABASES,
    ),
    no_enrichment: bool | None = typer.Option(
        None,
        "--no-enrichment",
        help="Disable target enrichment even when an index source is configured",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_ref_path: Path | None = typer.Option(
        None,
        "--sourmash-ref-path",
        help="Path to a prebuilt sourmash reference sketch database",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_ref_url: str | None = typer.Option(
        None,
        "--sourmash-ref-url",
        help="URL to download a prebuilt sourmash reference sketch database",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_ref_fasta: Path | None = typer.Option(
        None,
        "--sourmash-ref-fasta",
        help="Local FASTA to sketch as an experimental sourmash reference database",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_lineages_path: Path | None = typer.Option(
        None,
        "--sourmash-lineages-path",
        help="Path to a sourmash taxonomy lineages CSV matching the reference sketch database",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_lineages_url: str | None = typer.Option(
        None,
        "--sourmash-lineages-url",
        help="URL to download a sourmash taxonomy lineages CSV matching the reference sketch database",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_ksize: int | None = typer.Option(
        None,
        "--sourmash-ksize",
        help="K-mer size for experimental sourmash sketching (default: 31)",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_scaled: int | None = typer.Option(
        None,
        "--sourmash-scaled",
        help="Scaled value for experimental sourmash sketching (default: 50)",
        rich_help_panel=PANEL_DATABASES,
    ),
    sourmash_threshold_bp: int | None = typer.Option(
        None,
        "--sourmash-threshold-bp",
        help="Minimum estimated base-pair overlap for experimental sourmash gather (default: 50)",
        rich_help_panel=PANEL_DATABASES,
    ),
    virus_kmer_size: int | None = typer.Option(
        None,
        "--virus-kmer-size",
        help="K-mer size for building a custom virus enrichment index (default: 31)",
        rich_help_panel=PANEL_DATABASES,
    ),
    virus_window_size: int | None = typer.Option(
        None,
        "--virus-window-size",
        help="Minimizer window size for building a custom virus enrichment index (default: 1)",
        rich_help_panel=PANEL_DATABASES,
    ),
    virus_abs_threshold: int | None = typer.Option(
        None,
        "--virus-abs-threshold",
        help="Minimum absolute minimizer hits for virus read enrichment (default: 1)",
        rich_help_panel=PANEL_DATABASES,
    ),
    virus_rel_threshold: float | None = typer.Option(
        None,
        "--virus-rel-threshold",
        help="Minimum relative minimizer proportion for virus read enrichment (default: 0.0)",
        rich_help_panel=PANEL_DATABASES,
    ),
    # -------------------------------------------------------------------------
    # Database Versions
    # -------------------------------------------------------------------------
    blast_db_version: str | None = typer.Option(
        None,
        "--blast-db-version",
        help="BLAST database version (auto-resolved from registry if path registered)",
        rich_help_panel=PANEL_DATABASES,
    ),
    # -------------------------------------------------------------------------
    # Analysis Parameters
    # -------------------------------------------------------------------------
    cutoff_percent: float | None = typer.Option(
        None,
        "--cutoff-percent",
        help="Cutoff percentage (default: 0.001)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    tax_stringency: float | None = typer.Option(
        None,
        "--tax-stringency",
        help="Taxonomy stringency (default: 0.7)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    entropy: float | None = typer.Option(
        None,
        "--entropy",
        help="Entropy threshold (default: 0.9)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    max_blast_targets: int | None = typer.Option(
        None,
        "--max-blast-targets",
        help="Maximum BLAST targets (default: 100)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    blast_retention_count: int | None = typer.Option(
        None,
        "--blast-retention-count",
        help="Retain top X BLAST hits (default: 5)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    min_consecutive_bases: int | None = typer.Option(
        None,
        "--min-consecutive-bases",
        help="Minimum consecutive bases (default: 200)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    qtrim: str | None = typer.Option(
        None,
        "--qtrim",
        help="Quality trimming mode (default: 't')",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    include_children: bool | None = typer.Option(
        None,
        "--include-children/--no-include-children",
        help="Include children in taxonomy (default: true)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    max_concurrent_downloads: int | None = typer.Option(
        None,
        "--max-concurrent-downloads",
        help="Maximum concurrent SRA downloads (default: 3)",
        rich_help_panel=PANEL_ANALYSIS,
    ),
    # -------------------------------------------------------------------------
    # Read Preprocessing
    # -------------------------------------------------------------------------
    preprocess: bool | None = typer.Option(
        None,
        "--preprocess/--no-preprocess",
        help="Enable default preprocessing steps",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    dedup: bool | None = typer.Option(
        None,
        "--dedup",
        help="Deduplicate reads (umbrella: enables both --dedup-seq and --dedup-pos)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    dedup_seq: bool | None = typer.Option(
        None,
        "--dedup-seq",
        help="Enable sequence-based deduplication with clumpify",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    dedup_pos: bool | None = typer.Option(
        None,
        "--dedup-pos",
        help="Enable positional deduplication with samtools markdup",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    trim_adapters: bool | None = typer.Option(
        None,
        "--trim-adapters/--no-trim-adapters",
        help="Trim Illumina adapters (default: follows --preprocess)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    merge_pairs: bool | None = typer.Option(
        None,
        "--merge-pairs/--no-merge-pairs",
        help="Merge overlapping paired-end reads before contig mapback",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    host_index: Path | None = typer.Option(
        None,
        "--host-index",
        help="Path to prebuilt host/contaminant index (.idx file)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    host_index_url: str | None = typer.Option(
        None,
        "--host-index-url",
        help="URL to download a prebuilt host/contaminant index",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    host_contaminants_fasta: Path | None = typer.Option(
        None,
        "--host-contaminants-fasta",
        help="Custom contaminant FASTA to build and union with other host indexes",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    filter_reads: bool | None = typer.Option(
        None,
        "--filter-reads/--no-filter-reads",
        help="Filter reads by quality/length (default: follows --preprocess)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    filter_low_complexity_reads: bool | None = typer.Option(
        None,
        "--filter-low-complexity-reads/--no-filter-low-complexity-reads",
        help="Filter reads with low normalized 5-mer entropy (default: false)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    min_read_quality_illumina: int | None = typer.Option(
        None,
        "--min-read-quality-illumina",
        help="Minimum average quality for Illumina reads (default: 20)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    min_read_quality_nanopore: int | None = typer.Option(
        None,
        "--min-read-quality-nanopore",
        help="Minimum average quality for Nanopore reads (default: 12)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    min_read_length: int | None = typer.Option(
        None,
        "--min-read-length",
        help="Minimum read length (default: 50)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    max_read_length: int | None = typer.Option(
        None,
        "--max-read-length",
        help="Maximum read length (default: no limit)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    min_read_entropy: float | None = typer.Option(
        None,
        "--min-read-entropy",
        help="Minimum normalized 5-mer entropy over 50-base windows (default: 0.9)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    # -------------------------------------------------------------------------
    # LabKey Integration
    # -------------------------------------------------------------------------
    labkey: bool | None = typer.Option(
        None,
        "--labkey/--no-labkey",
        help="Enable LabKey integration (requires all labkey-* params to be set)",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_server: str | None = typer.Option(
        None,
        "--labkey-server",
        help="LabKey server URL (e.g., 'dholk.primate.wisc.edu')",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_project_name: str | None = typer.Option(
        None,
        "--labkey-project-name",
        help="LabKey project name/path (e.g., 'dho/projects/lungfish/InfinitePath')",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_webdav: str | None = typer.Option(
        None,
        "--labkey-webdav",
        help="LabKey WebDAV URL for file uploads",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_schema: str | None = typer.Option(
        None,
        "--labkey-schema",
        help="LabKey database schema name (e.g., 'lists')",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_blast_meta_hits_list: str | None = typer.Option(
        None,
        "--labkey-blast-meta-hits-list",
        help="LabKey list name for BLAST metagenomic hits",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_blast_fasta_list: str | None = typer.Option(
        None,
        "--labkey-blast-fasta-list",
        help="LabKey list name for BLAST FASTA results",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_exp_id_guard_list: str | None = typer.Option(
        None,
        "--labkey-exp-id-guard-list",
        help="LabKey list name for experiment ID guard (prevents duplicate uploads)",
        rich_help_panel=PANEL_LABKEY,
    ),
    # -------------------------------------------------------------------------
    # Notifications
    # -------------------------------------------------------------------------
    slack_channel: str | None = typer.Option(
        None,
        "--slack-channel",
        help="Slack channel ID for notifications (e.g., 'C0123456789')",
        rich_help_panel=PANEL_NOTIFICATIONS,
    ),
    no_slack: bool = typer.Option(
        False,
        "--no-slack",
        help="Disable Slack notifications for this run",
        rich_help_panel=PANEL_NOTIFICATIONS,
    ),
    # -------------------------------------------------------------------------
    # Execution Control
    # -------------------------------------------------------------------------
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Show command without executing",
        rich_help_panel=PANEL_CORE,
    ),
) -> None:
    """
    Run the NVD pipeline.

    This command wraps 'nextflow run' with a simpler interface, running the
    pipeline from the local installation. Reference paths and settings are
    loaded from NVD_CONFIG_DIR/user.config unless overridden with command-line
    options or the NVD_CONFIG file override.

    The command is saved to .nfresume for easy resumption with 'nvd resume'.

    Parameter precedence (highest to lowest):
        1. CLI arguments (--blast-db, --dedup, etc.)
        2. Params file (--params-file params.yaml)
        3. Preset values (--preset production)
        4. Pipeline defaults (nextflow.config)

    The samplesheet can come from any of these sources. For example, you can
    define it in a params file or preset and omit --samplesheet from the CLI.

    Any arguments after '--' are passed directly to Nextflow. This allows
    using Nextflow-native options like -with-tower, -with-trace, etc.

    Examples:

        # Simple run with auto-detected profile
        nvd run -s samples.csv -e exp001

        # Specify profile and results directory
        nvd run -s samples.csv -e exp002 -p docker -r ./my_results

        # Resume a failed run (two ways)
        nvd run -s samples.csv -e exp003 --resume
        nvd resume  # Uses saved command from .nfresume

        # Use a preset (CLI args override preset values)
        nvd run -s samples.csv -e exp004 --preset production

        # Use a preset with overrides
        nvd run -s samples.csv -e exp005 --preset production --cutoff-percent 0.01

        # Load params from a YAML file (can include samplesheet)
        nvd run --params-file run-config.yaml

        # Combine params file with CLI overrides
        nvd run -f run-config.yaml --dedup --cutoff-percent 0.01

        # Pass other Nextflow options via '--' separator
        nvd run -s samples.csv -- -with-tower -with-trace
    """

    # =========================================================================
    # STEP 1: Load preset params (if specified)
    # =========================================================================
    preset_params: dict[str, Any] = {}
    if preset:
        preset_store = get_preset_store()
        preset_obj = preset_store.get(preset)
        if not preset_obj:
            available = preset_store.list()
            console.print(f"[red]✗ Preset not found: {preset}[/red]")
            if available:
                console.print("\n[cyan]Available presets:[/cyan]")
                for p in available:
                    desc = f" - {p.description}" if p.description else ""
                    console.print(f"  • {p.name}{desc}")
            else:
                console.print("\n[dim]No presets registered.[/dim]")
                console.print(
                    "[dim]Create one with: nvd preset register <name> --from-file params.yaml[/dim]",
                )
            raise typer.Exit(1)

        preset_params = preset_obj.params
        info(f"Using preset '{preset}' ({len(preset_params)} parameters)")

    # =========================================================================
    # STEP 2: Build CLI args dict and merge with precedence
    #
    # Precedence (highest to lowest): CLI args > params_file > preset > defaults
    # NvdParams.merge() handles None filtering and validation.
    # Nextflow-native options (profile, config, resume) are handled separately.
    # =========================================================================

    # All pipeline params from CLI (None values are filtered by merge)
    # NOTE: profile, config, resume are Nextflow-native, not pipeline params
    resolved_taxonomy_dir = taxonomy_dir
    if resolved_taxonomy_dir is None and ENV_VAR_TAXONOMY in os.environ:
        resolved_taxonomy_dir = Path(os.environ[ENV_VAR_TAXONOMY])

    cli_args: dict[str, Any] = {
        # Core
        "samplesheet": samplesheet,
        "experiment_id": experiment_id,
        "results": results,
        "cleanup": cleanup,
        "work_dir": work_dir,
        "taxonomy_dir": resolved_taxonomy_dir,
        "taxonomy_mode": taxonomy_mode,
        "taxonomy_refresh": taxonomy_refresh,
        "taxonomy_max_age_days": taxonomy_max_age_days,
        "experimental": experimental,
        "skip_assembly": skip_assembly,
        "skip_blast": skip_blast,
        # Reference paths
        "blast_db": blast_db,
        "blast_db_prefix": blast_db_prefix,
        "virus_index": virus_index,
        "virus_index_url": virus_index_url,
        "virus_reference_fasta": virus_reference_fasta,
        "no_enrichment": no_enrichment,
        "virus_kmer_size": virus_kmer_size,
        "virus_window_size": virus_window_size,
        "virus_abs_threshold": virus_abs_threshold,
        "virus_rel_threshold": virus_rel_threshold,
        "sourmash_ref_path": sourmash_ref_path,
        "sourmash_ref_url": sourmash_ref_url,
        "sourmash_ref_fasta": sourmash_ref_fasta,
        "sourmash_lineages_path": sourmash_lineages_path,
        "sourmash_lineages_url": sourmash_lineages_url,
        "sourmash_ksize": sourmash_ksize,
        "sourmash_scaled": sourmash_scaled,
        "sourmash_threshold_bp": sourmash_threshold_bp,
        # Reference versions
        "blast_db_version": blast_db_version,
        # Analysis
        "cutoff_percent": cutoff_percent,
        "tax_stringency": tax_stringency,
        "entropy": entropy,
        "max_blast_targets": max_blast_targets,
        "blast_retention_count": blast_retention_count,
        "min_consecutive_bases": min_consecutive_bases,
        "qtrim": qtrim,
        "include_children": include_children,
        "max_concurrent_downloads": max_concurrent_downloads,
        # Preprocessing
        "preprocess": preprocess,
        "dedup": dedup,
        "dedup_seq": dedup_seq,
        "dedup_pos": dedup_pos,
        "trim_adapters": trim_adapters,
        "merge_pairs": merge_pairs,
        "host_index": host_index,
        "host_index_url": host_index_url,
        "host_contaminants_fasta": host_contaminants_fasta,
        "filter_reads": filter_reads,
        "filter_low_complexity_reads": filter_low_complexity_reads,
        "min_read_quality_illumina": min_read_quality_illumina,
        "min_read_quality_nanopore": min_read_quality_nanopore,
        "min_read_length": min_read_length,
        "max_read_length": max_read_length,
        "min_read_entropy": min_read_entropy,
        # LabKey
        "labkey": labkey,
        "labkey_server": labkey_server,
        "labkey_project_name": labkey_project_name,
        "labkey_webdav": labkey_webdav,
        "labkey_schema": labkey_schema,
        "labkey_blast_meta_hits_list": labkey_blast_meta_hits_list,
        "labkey_blast_fasta_list": labkey_blast_fasta_list,
        "labkey_exp_id_guard_list": labkey_exp_id_guard_list,
        # Notifications
        "slack_enabled": False if no_slack else None,  # Only override if --no-slack
        "slack_channel": slack_channel,
    }

    # Load params file if provided (we merge it ourselves, not Nextflow)
    params_file_dict = load_params_file(params_file) if params_file else None

    # Merge with precedence: preset < params_file < CLI
    # NvdParams.merge() filters None values and applies validation
    params = NvdParams.merge(preset_params, params_file_dict, cli_args)

    # =========================================================================
    # STEP 3: Post-merge validation
    # =========================================================================

    # Validate samplesheet is present (required, but can come from any source)
    if params.samplesheet is None:
        error(
            "Samplesheet is required. Provide via --samplesheet, --preset, or --params-file",
        )
        raise typer.Exit(1)

    if not params.samplesheet.exists():
        error(f"Samplesheet not found: {params.samplesheet}")
        raise typer.Exit(1)

    # =========================================================================
    # STEP 4: Handle Nextflow-native options (not pipeline params)
    # =========================================================================

    # Determine execution profile:
    # 1. Command-line --profile takes precedence
    # 2. NVD_DEFAULT_PROFILE from NVD_CONFIG_DIR/setup.conf
    # 3. Auto-detect based on available container runtime
    effective_profile = profile
    if effective_profile is None:
        default_profile = get_default_profile()
        if default_profile:
            effective_profile = default_profile
            info(f"Using default profile from setup.conf: {effective_profile}")
        else:
            effective_profile = auto_detect_profile()
            info(f"Auto-detected execution profile: {effective_profile}")

    # Find config file (user.config auto-discovery)
    effective_config = find_config_file(config)
    if effective_config:
        info(f"Using config: {effective_config}")
    else:
        warning("No config file found. Using command-line parameters only.")
        warning(f"Consider creating a config file at: {DEFAULT_CONFIG}")

    # Check Nextflow is installed before real execution. Dry-run mode only needs
    # to render the command, which is useful on systems where Nextflow is not
    # installed yet.
    if not dry_run and not check_command_exists("nextflow"):
        error(
            "Nextflow not found. Please install Nextflow:\n"
            "  curl -s https://get.nextflow.io | bash\n"
            "  Or visit: https://www.nextflow.io/docs/latest/getstarted.html",
        )

    # =========================================================================
    # STEP 5: Build and execute command
    # =========================================================================

    # Build command with pipeline params
    pipeline_root = get_pipeline_root()
    cmd = params.to_nextflow_args(pipeline_root)

    # Add Nextflow-native options (not pipeline params)
    if effective_profile:
        cmd.extend(["-profile", effective_profile])
    if effective_config:
        cmd.extend(["-c", str(effective_config)])
    if resume:
        cmd.append("-resume")

    # Append any extra args passed through to Nextflow (e.g., -with-tower)
    if ctx.args:
        cmd.extend(ctx.args)

    # Format for display (pretty-printed with continuations)
    cmd_display = format_command_for_display(cmd)

    # Single-line version for .nfresume file
    cmd_str = " ".join(cmd)

    # Show command
    console.print("\n[bold]Executing command:[/bold]")
    console.print(f"[dim]{cmd_display}[/dim]\n")

    if dry_run:
        success("Dry-run mode: command shown above but not executed")
        return

    # Save resume-enabled version for 'nvd resume' command
    resume_cmd = cmd_str if "-resume" in cmd_str else f"{cmd_str} -resume"
    RESUME_FILE.write_text(resume_cmd, encoding="utf-8")
    info(f"Saved resume command to {RESUME_FILE}")

    # Execute nextflow
    try:
        result = subprocess.run(  # noqa: S603
            cmd,
            check=False,
            env=nextflow_environment(pipeline_root),
        )
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        console.print("\n[yellow]Pipeline interrupted by user[/yellow]")
        sys.exit(130)
    except OSError as e:
        error(f"Failed to execute Nextflow: {e}")
