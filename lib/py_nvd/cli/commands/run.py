"""
Run command for the NVD CLI.

Commands:
    nvd run  - Run the NVD2 pipeline
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path  # noqa: TC003
from typing import Any

import typer

from py_nvd import state
from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    PANEL_ANALYSIS,
    PANEL_CORE,
    PANEL_DATABASES,
    PANEL_LABKEY,
    PANEL_NOTIFICATIONS,
    PANEL_PREPROCESSING,
    PANEL_SRA,
    PIPELINE_ROOT,
    RESUME_FILE,
    auto_detect_profile,
    check_command_exists,
    console,
    error,
    find_config_file,
    get_default_profile,
    info,
    success,
    warning,
)
from py_nvd.models import NvdParams
from py_nvd.params import load_params_file
from py_nvd.state import resolve_database_versions

_MIN_CMD_PARTS_FOR_CONTINUATION = 3


def _format_command_for_display(cmd: list[str]) -> str:
    """
    Format a command list for readable display with line continuations.

    Each argument pair (--flag value) gets its own line, indented and
    with shell continuation characters for copy-paste compatibility.
    """
    if len(cmd) < _MIN_CMD_PARTS_FOR_CONTINUATION:
        return " ".join(cmd)

    lines = []
    # First line: nextflow run <pipeline_root>
    lines.append(f"{cmd[0]} {cmd[1]} {cmd[2]} \\")

    # Process remaining args in pairs where possible
    i = 3
    while i < len(cmd):
        arg = cmd[i]

        # Check if this is a flag that takes a value (not a standalone flag like -resume)
        if i + 1 < len(cmd) and not cmd[i + 1].startswith("-"):
            # Flag with value: --param value
            value = cmd[i + 1]
            if i + 2 < len(cmd):
                lines.append(f"    {arg} {value} \\")
            else:
                lines.append(f"    {arg} {value}")
            i += 2
        else:
            # Standalone flag like -resume
            if i + 1 < len(cmd):
                lines.append(f"    {arg} \\")
            else:
                lines.append(f"    {arg}")
            i += 1

    return "\n".join(lines)


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
    # NOTE: Default is None, not "all" - this lets us detect if user provided it
    # vs using preset value. Default "all" is applied after preset merge.
    tools: str | None = typer.Option(
        None,
        "--tools",
        "-t",
        help="Workflow(s) to run: all, stat_blast, gottcha, blast, clumpify. "
        "Combine with commas, e.g. 'stat_blast,clumpify' (default: all)",
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
        help="Custom config file (default: ~/.nvd/user.config or NVD_CONFIG)",
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
    state_dir: Path | None = typer.Option(
        None,
        "--state-dir",
        help="State directory for run tracking (default: ~/.nvd/ or NVD_STATE_DIR)",
        rich_help_panel=PANEL_CORE,
    ),
    stateless: bool | None = typer.Option(
        None,
        "--stateless/--no-stateless",
        help="Run without state management (disables run tracking, LabKey, Slack)",
        rich_help_panel=PANEL_CORE,
    ),
    taxonomy_dir: Path | None = typer.Option(
        None,
        "--taxonomy-dir",
        help="Taxonomy database directory (required when --stateless)",
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
    # -------------------------------------------------------------------------
    # Database Paths
    # -------------------------------------------------------------------------
    gottcha2_db: Path | None = typer.Option(
        None,
        "--gottcha2-db",
        help="Override GOTTCHA2 database path",
        rich_help_panel=PANEL_DATABASES,
    ),
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
    stat_index: Path | None = typer.Option(
        None,
        "--stat-index",
        help="Override STAT index file",
        rich_help_panel=PANEL_DATABASES,
    ),
    stat_dbss: Path | None = typer.Option(
        None,
        "--stat-dbss",
        help="Override STAT dbss file",
        rich_help_panel=PANEL_DATABASES,
    ),
    stat_annotation: Path | None = typer.Option(
        None,
        "--stat-annotation",
        help="Override STAT annotation file",
        rich_help_panel=PANEL_DATABASES,
    ),
    human_virus_taxlist: Path | None = typer.Option(
        None,
        "--human-virus-taxlist",
        help="Override human virus taxlist file",
        rich_help_panel=PANEL_DATABASES,
    ),
    # -------------------------------------------------------------------------
    # Database Versions
    # -------------------------------------------------------------------------
    gottcha2_db_version: str | None = typer.Option(
        None,
        "--gottcha2-db-version",
        help="GOTTCHA2 database version (auto-resolved from registry if path registered)",
        rich_help_panel=PANEL_DATABASES,
    ),
    blast_db_version: str | None = typer.Option(
        None,
        "--blast-db-version",
        help="BLAST database version (auto-resolved from registry if path registered)",
        rich_help_panel=PANEL_DATABASES,
    ),
    stat_db_version: str | None = typer.Option(
        None,
        "--stat-db-version",
        help="STAT database version (auto-resolved from registry if path registered)",
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
    min_gottcha_reads: int | None = typer.Option(
        None,
        "--min-gottcha-reads",
        help="Minimum reads for GOTTCHA2 (default: 250)",
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
        help="Enable all preprocessing steps (dedup, trim, scrub, filter)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    merge_pairs: bool | None = typer.Option(
        None,
        "--merge-pairs/--no-merge-pairs",
        help="Merge paired read mates based on overlaps",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    dedup: bool | None = typer.Option(
        None,
        "--dedup/--no-dedup",
        help="Deduplicate reads (default: follows --preprocess)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    trim_adapters: bool | None = typer.Option(
        None,
        "--trim-adapters/--no-trim-adapters",
        help="Trim Illumina adapters (default: follows --preprocess)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    scrub_host_reads: bool | None = typer.Option(
        None,
        "--scrub-host-reads/--no-scrub-host-reads",
        help="Remove host reads with STAT (requires --sra-human-db; default: follows --preprocess)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    filter_reads: bool | None = typer.Option(
        None,
        "--filter-reads/--no-filter-reads",
        help="Filter reads by quality/length (default: follows --preprocess)",
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
    # -------------------------------------------------------------------------
    # SRA Submission
    # -------------------------------------------------------------------------
    sra_human_db: Path | None = typer.Option(
        None,
        "--sra-human-db",
        help="Path to human reads database for SRA submission scrubbing",
        rich_help_panel=PANEL_SRA,
    ),
    # DEPRECATED: Use --sra-human-db instead
    human_read_scrub: Path | None = typer.Option(
        None,
        "--human-read-scrub",
        help="[DEPRECATED] Use --sra-human-db instead",
        hidden=True,
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
    labkey_gottcha_fasta_list: str | None = typer.Option(
        None,
        "--labkey-gottcha-fasta-list",
        help="LabKey list name for GOTTCHA2 FASTA results",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_gottcha_full_list: str | None = typer.Option(
        None,
        "--labkey-gottcha-full-list",
        help="LabKey list name for full GOTTCHA2 results",
        rich_help_panel=PANEL_LABKEY,
    ),
    labkey_gottcha_blast_verified_full_list: str | None = typer.Option(
        None,
        "--labkey-gottcha-blast-verified-full-list",
        help="LabKey list name for BLAST-verified GOTTCHA2 results",
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
    Run the NVD2 pipeline.

    This command wraps 'nextflow run' with a simpler interface, running the
    pipeline from the local installation. Database paths and settings are
    loaded from ~/.nvd/user.config unless overridden with command-line
    options or NVD_CONFIG env var.

    The command is saved to .nfresume for easy resumption with 'nvd resume'.

    Parameter precedence (highest to lowest):
        1. CLI arguments (--tools, --blast-db, etc.)
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

        # Run only GOTTCHA2 workflow
        nvd run -s samples.csv -e exp003 -t gottcha

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
        nvd run -f run-config.yaml --tools blast --cutoff-percent 0.01

        # Pass other Nextflow options via '--' separator
        nvd run -s samples.csv -- -with-tower -with-trace
    """

    # =========================================================================
    # STEP 1: Load preset params (if specified)
    # =========================================================================
    preset_params: dict[str, Any] = {}
    if preset:
        preset_obj = state.get_preset(preset)
        if not preset_obj:
            available = state.list_presets()
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

    # Handle deprecation warning for human_read_scrub -> sra_human_db
    effective_sra_human_db = sra_human_db
    if human_read_scrub is not None and sra_human_db is None:
        warning(
            "DEPRECATION: --human-read-scrub is deprecated. Please use --sra-human-db instead.",
        )
        effective_sra_human_db = human_read_scrub

    # All pipeline params from CLI (None values are filtered by merge)
    # NOTE: profile, config, resume are Nextflow-native, not pipeline params
    cli_args: dict[str, Any] = {
        # Core
        "samplesheet": samplesheet,
        "experiment_id": experiment_id,
        "tools": tools,
        "results": results,
        "cleanup": cleanup,
        "work_dir": work_dir,
        "state_dir": state_dir,
        "stateless": stateless,
        "taxonomy_dir": taxonomy_dir,
        # Database paths
        "gottcha2_db": gottcha2_db,
        "blast_db": blast_db,
        "blast_db_prefix": blast_db_prefix,
        "stat_index": stat_index,
        "stat_dbss": stat_dbss,
        "stat_annotation": stat_annotation,
        "human_virus_taxlist": human_virus_taxlist,
        # Database versions
        "gottcha2_db_version": gottcha2_db_version,
        "blast_db_version": blast_db_version,
        "stat_db_version": stat_db_version,
        # Analysis
        "cutoff_percent": cutoff_percent,
        "tax_stringency": tax_stringency,
        "entropy": entropy,
        "min_gottcha_reads": min_gottcha_reads,
        "max_blast_targets": max_blast_targets,
        "blast_retention_count": blast_retention_count,
        "min_consecutive_bases": min_consecutive_bases,
        "qtrim": qtrim,
        "include_children": include_children,
        "max_concurrent_downloads": max_concurrent_downloads,
        # Preprocessing
        "preprocess": preprocess,
        "merge_pairs": merge_pairs,
        "dedup": dedup,
        "trim_adapters": trim_adapters,
        "scrub_host_reads": scrub_host_reads,
        "filter_reads": filter_reads,
        "min_read_quality_illumina": min_read_quality_illumina,
        "min_read_quality_nanopore": min_read_quality_nanopore,
        "min_read_length": min_read_length,
        "max_read_length": max_read_length,
        # SRA (with deprecation handling)
        "sra_human_db": effective_sra_human_db,
        # LabKey
        "labkey": labkey,
        "labkey_server": labkey_server,
        "labkey_project_name": labkey_project_name,
        "labkey_webdav": labkey_webdav,
        "labkey_schema": labkey_schema,
        "labkey_gottcha_fasta_list": labkey_gottcha_fasta_list,
        "labkey_gottcha_full_list": labkey_gottcha_full_list,
        "labkey_gottcha_blast_verified_full_list": labkey_gottcha_blast_verified_full_list,
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

    # Tools validation is handled by NvdParams, just show info
    if params.tools and params.tools != "all":
        info(f"Using tools: {params.tools}")

    # =========================================================================
    # STEP 3b: Resolve database versions from registry
    # =========================================================================
    resolution = resolve_database_versions(
        blast_db=params.blast_db,
        blast_db_version=params.blast_db_version,
        gottcha2_db=params.gottcha2_db,
        gottcha2_db_version=params.gottcha2_db_version,
        stat_index=params.stat_index,
        stat_db_version=params.stat_db_version,
        state_dir=state_dir,
    )

    # Display warnings (unregistered paths, version mismatches)
    for warn in resolution.warnings:
        warning(warn)

    # Display info for auto-registrations
    for db_type, version, path in resolution.auto_registered:
        info(f"Registered {db_type} database: {version} at {path}")

    # Apply resolved versions to params
    resolved_updates: dict[str, Any] = {}
    if resolution.blast_db_version is not None:
        resolved_updates["blast_db_version"] = resolution.blast_db_version
    if resolution.gottcha2_db_version is not None:
        resolved_updates["gottcha2_db_version"] = resolution.gottcha2_db_version
    if resolution.stat_db_version is not None:
        resolved_updates["stat_db_version"] = resolution.stat_db_version

    if resolved_updates:
        params = NvdParams.merge(params, resolved_updates)

    # =========================================================================
    # STEP 4: Handle Nextflow-native options (not pipeline params)
    # =========================================================================

    # Determine execution profile:
    # 1. Command-line --profile takes precedence
    # 2. NVD_DEFAULT_PROFILE from ~/.nvd/setup.conf
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

    # Check Nextflow is installed
    if not check_command_exists("nextflow"):
        error(
            "Nextflow not found. Please install Nextflow:\n"
            "  curl -s https://get.nextflow.io | bash\n"
            "  Or visit: https://www.nextflow.io/docs/latest/getstarted.html",
        )

    # =========================================================================
    # STEP 5: Build and execute command
    # =========================================================================

    # Build command with pipeline params
    cmd = params.to_nextflow_args(PIPELINE_ROOT)

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
    cmd_display = _format_command_for_display(cmd)

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
        result = subprocess.run(cmd, check=False)  # noqa: S603
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        console.print("\n[yellow]Pipeline interrupted by user[/yellow]")
        sys.exit(130)
    except OSError as e:
        error(f"Failed to execute Nextflow: {e}")
