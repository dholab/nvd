# ruff: noqa: B008, FBT001, FBT003
"""
Run command for the NVD CLI.

Commands:
    nvd run  - Run the NVD2 pipeline
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any

import typer

from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    PANEL_ANALYSIS,
    PANEL_CORE,
    PANEL_DATABASES,
    PANEL_LABKEY,
    PANEL_PREPROCESSING,
    PANEL_SRA,
    PIPELINE_ROOT,
    RESUME_FILE,
    VALID_PROFILES,
    VALID_TOOLS,
    auto_detect_profile,
    check_command_exists,
    console,
    error,
    find_config_file,
    info,
    success,
    warning,
)

# Nextflow-native options that use special syntax (not --param value)
NEXTFLOW_NATIVE_OPTIONS = frozenset({"profile", "config", "work_dir", "resume"})


def build_nextflow_command(params: dict[str, Any]) -> list[str]:
    """
    Construct the nextflow run command from a unified params dict.

    All parameters flow through this single dict. Special handling:
    - profile → -profile value
    - config → -c value
    - work_dir → -work-dir value
    - resume → -resume (flag)
    - All others → --param-name value

    Booleans are converted to "true"/"false" strings for Nextflow.
    Underscores in param names are converted to hyphens.
    """
    cmd = ["nextflow", "run", str(PIPELINE_ROOT)]

    # Handle Nextflow-native options first
    if params.get("profile"):
        cmd.extend(["-profile", params["profile"]])

    if params.get("config"):
        cmd.extend(["-c", str(params["config"])])

    if params.get("work_dir"):
        cmd.extend(["-work-dir", str(params["work_dir"])])

    if params.get("resume"):
        cmd.append("-resume")

    # Handle all pipeline parameters (--param-name value)
    for key, value in params.items():
        # Skip Nextflow-native options (already handled above)
        if key in NEXTFLOW_NATIVE_OPTIONS:
            continue

        # Skip None values
        if value is None:
            continue

        # Convert param name: Python snake_case → Nextflow kebab-case
        param_name = key.replace("_", "-")

        # Convert value to string, handling booleans specially
        if isinstance(value, bool):
            str_value = "true" if value else "false"
        elif isinstance(value, Path):
            str_value = str(value)
        else:
            str_value = str(value)

        cmd.extend([f"--{param_name}", str_value])

    return cmd


def _validate_tools(tools: str) -> None:
    """Validate tools option (supports comma-separated values)."""
    tools_list = [t.strip() for t in tools.split(",")]
    invalid_tools = [t for t in tools_list if t not in VALID_TOOLS]
    if invalid_tools:
        error(
            f"Invalid tools option(s): {', '.join(invalid_tools)}. "
            f"Must be one of: {', '.join(VALID_TOOLS)}",
        )


def _validate_profile(profile: str) -> None:
    """Validate profile option."""
    if profile not in VALID_PROFILES:
        error(
            f"Invalid profile: {profile}. Must be one of: {', '.join(VALID_PROFILES)}"
        )


def _format_command_for_display(cmd: list[str]) -> str:
    """
    Format a command list for readable display with line continuations.

    Each argument pair (--flag value) gets its own line, indented and
    with shell continuation characters for copy-paste compatibility.
    """
    if len(cmd) < 3:
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


def run(  # noqa: PLR0913, PLR0912, PLR0915, C901
    # Complexity/boolean/B008 warnings are acceptable for CLI with many options
    ctx: typer.Context,
    # -------------------------------------------------------------------------
    # Core Options
    # -------------------------------------------------------------------------
    samplesheet: Path = typer.Option(
        ...,
        "--samplesheet",
        "-s",
        help="Path to samplesheet CSV",
        exists=True,
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
        help="JSON/YAML params file (passed to Nextflow, its precedence rules apply)",
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
        help="Remove host reads with hostile (default: follows --preprocess)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    hostile_index: Path | None = typer.Option(
        None,
        "--hostile-index",
        help="Path to local hostile index (for offline use)",
        rich_help_panel=PANEL_PREPROCESSING,
    ),
    hostile_index_name: str | None = typer.Option(
        None,
        "--hostile-index-name",
        help="Standard hostile index name (default: human-t2t-hla.rs-viral-202401_ml-phage-202401)",
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
        2. Preset values (--preset production)
        3. User config file (~/.nvd/user.config)
        4. Pipeline defaults (nextflow.config)

    Any arguments after '--' are passed directly to Nextflow. This allows
    using Nextflow-native options like -params-file, -with-tower, -with-trace,
    etc. Nextflow handles precedence for -params-file contents.

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

        # Load additional params from a JSON/YAML file
        nvd run -s samples.csv --params-file extra-params.json

        # Pass other Nextflow options via '--' separator
        nvd run -s samples.csv -- -with-tower -with-trace
    """

    # =========================================================================
    # STEP 1: Load preset params (if specified)
    # =========================================================================
    preset_params: dict[str, Any] = {}
    if preset:
        from py_nvd import state

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
                    "[dim]Create one with: nvd preset register <name> --from-file params.yaml[/dim]"
                )
            raise typer.Exit(1)

        preset_params = preset_obj.params
        info(f"Using preset '{preset}' ({len(preset_params)} parameters)")

    # =========================================================================
    # STEP 2: Build unified params dict with proper precedence
    #
    # Precedence: CLI args > preset params > defaults
    # We start with preset, then overlay CLI args (if provided, i.e., not None)
    # =========================================================================
    params: dict[str, Any] = {}

    # Layer 1: Start with preset params
    for key, value in preset_params.items():
        if value is not None:
            params[key] = value

    # Layer 2: CLI args override (only if actually provided, not None)
    # Required params (always provided)
    params["samplesheet"] = samplesheet

    # Optional params - only override if CLI provided a value
    if experiment_id is not None:
        params["experiment_id"] = experiment_id
    if tools is not None:
        params["tools"] = tools
    if results is not None:
        params["results"] = results
    if profile is not None:
        params["profile"] = profile
    if config is not None:
        params["config"] = config
    if resume is not None:
        params["resume"] = resume
    if cleanup is not None:
        params["cleanup"] = cleanup
    if work_dir is not None:
        params["work_dir"] = work_dir
    if state_dir is not None:
        params["state_dir"] = state_dir

    # Database paths
    if gottcha2_db is not None:
        params["gottcha2_db"] = gottcha2_db
    if blast_db is not None:
        params["blast_db"] = blast_db
    if blast_db_prefix is not None:
        params["blast_db_prefix"] = blast_db_prefix
    if stat_index is not None:
        params["stat_index"] = stat_index
    if stat_dbss is not None:
        params["stat_dbss"] = stat_dbss
    if stat_annotation is not None:
        params["stat_annotation"] = stat_annotation
    if human_virus_taxlist is not None:
        params["human_virus_taxlist"] = human_virus_taxlist

    # Analysis parameters
    if cutoff_percent is not None:
        params["cutoff_percent"] = cutoff_percent
    if tax_stringency is not None:
        params["tax_stringency"] = tax_stringency
    if entropy is not None:
        params["entropy"] = entropy
    if min_gottcha_reads is not None:
        params["min_gottcha_reads"] = min_gottcha_reads
    if max_blast_targets is not None:
        params["max_blast_targets"] = max_blast_targets
    if blast_retention_count is not None:
        params["blast_retention_count"] = blast_retention_count
    if min_consecutive_bases is not None:
        params["min_consecutive_bases"] = min_consecutive_bases
    if qtrim is not None:
        params["qtrim"] = qtrim
    if include_children is not None:
        params["include_children"] = include_children
    if max_concurrent_downloads is not None:
        params["max_concurrent_downloads"] = max_concurrent_downloads

    # Preprocessing options
    if preprocess is not None:
        params["preprocess"] = preprocess
    if merge_pairs is not None:
        params["merge_pairs"] = merge_pairs
    if dedup is not None:
        params["dedup"] = dedup
    if trim_adapters is not None:
        params["trim_adapters"] = trim_adapters
    if scrub_host_reads is not None:
        params["scrub_host_reads"] = scrub_host_reads
    if hostile_index is not None:
        params["hostile_index"] = hostile_index
    if hostile_index_name is not None:
        params["hostile_index_name"] = hostile_index_name
    if filter_reads is not None:
        params["filter_reads"] = filter_reads
    if min_read_quality_illumina is not None:
        params["min_read_quality_illumina"] = min_read_quality_illumina
    if min_read_quality_nanopore is not None:
        params["min_read_quality_nanopore"] = min_read_quality_nanopore
    if min_read_length is not None:
        params["min_read_length"] = min_read_length
    if max_read_length is not None:
        params["max_read_length"] = max_read_length

    # SRA options (handle deprecation)
    if human_read_scrub is not None and sra_human_db is None:
        warning(
            "DEPRECATION: --human-read-scrub is deprecated. Please use --sra-human-db instead.",
        )
        params["sra_human_db"] = human_read_scrub
    elif sra_human_db is not None:
        params["sra_human_db"] = sra_human_db

    # LabKey options
    if labkey is not None:
        params["labkey"] = labkey
    if labkey_server is not None:
        params["labkey_server"] = labkey_server
    if labkey_project_name is not None:
        params["labkey_project_name"] = labkey_project_name
    if labkey_webdav is not None:
        params["labkey_webdav"] = labkey_webdav
    if labkey_schema is not None:
        params["labkey_schema"] = labkey_schema
    if labkey_gottcha_fasta_list is not None:
        params["labkey_gottcha_fasta_list"] = labkey_gottcha_fasta_list
    if labkey_gottcha_full_list is not None:
        params["labkey_gottcha_full_list"] = labkey_gottcha_full_list
    if labkey_gottcha_blast_verified_full_list is not None:
        params["labkey_gottcha_blast_verified_full_list"] = (
            labkey_gottcha_blast_verified_full_list
        )
    if labkey_blast_meta_hits_list is not None:
        params["labkey_blast_meta_hits_list"] = labkey_blast_meta_hits_list
    if labkey_blast_fasta_list is not None:
        params["labkey_blast_fasta_list"] = labkey_blast_fasta_list
    if labkey_exp_id_guard_list is not None:
        params["labkey_exp_id_guard_list"] = labkey_exp_id_guard_list

    # =========================================================================
    # STEP 3: Apply defaults for required params not yet set
    # =========================================================================
    params.setdefault("tools", "all")
    params.setdefault("resume", False)
    params.setdefault("cleanup", False)

    # =========================================================================
    # STEP 4: Auto-detect and validate
    # =========================================================================

    # Auto-detect profile if not specified
    if params.get("profile") is None:
        params["profile"] = auto_detect_profile()
        info(f"Auto-detected execution profile: {params['profile']}")
    else:
        _validate_profile(params["profile"])

    # Validate tools
    _validate_tools(params["tools"])
    if params["tools"] != "all":
        info(f"Using tools: {params['tools']}")

    # Find config file (user.config auto-discovery)
    config_file = find_config_file(params.get("config"))
    if config_file:
        params["config"] = config_file
        info(f"Using config: {config_file}")
    else:
        # Remove config key if no file found (don't pass None to Nextflow)
        params.pop("config", None)
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
    cmd = build_nextflow_command(params)

    # Append params file if provided (Nextflow handles precedence)
    if params_file is not None:
        cmd.extend(["-params-file", str(params_file)])

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
