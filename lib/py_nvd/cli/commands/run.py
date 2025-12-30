# ruff: noqa: B008, FBT001, FBT002, FBT003
"""
Run command for the NVD CLI.

Commands:
    nvd run  - Run the NVD2 pipeline
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import typer

from py_nvd.cli.utils import (
    DEFAULT_CONFIG,
    PANEL_ANALYSIS,
    PANEL_CORE,
    PANEL_DATABASES,
    PANEL_LABKEY,
    PANEL_PREPROCESSING,
    PANEL_SRA,
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


def build_nextflow_command(  # noqa: PLR0913
    samplesheet: Path,
    experiment_id: str,
    tools: str,
    profile: str,
    config: Path | None,
    results: Path | None = None,
    resume: bool = False,
    cleanup: bool = False,
    work_dir: Path | None = None,
    **extra_params: str | Path | float | bool,
) -> list[str]:
    """Construct the nextflow run command."""

    cmd = ["nextflow", "run", "dhoconno/nvd"]

    # Profile
    cmd.extend(["-profile", profile])

    # Config file
    if config:
        cmd.extend(["-c", str(config)])

    # Work directory
    if work_dir:
        cmd.extend(["-work-dir", str(work_dir)])

    # Resume flag
    if resume:
        cmd.append("-resume")

    # Required pipeline parameters
    cmd.extend(
        [
            "--samplesheet",
            str(samplesheet),
            "--experiment_id",
            experiment_id,
            "--tools",
            tools,
        ],
    )

    # Optional pipeline parameters
    if results:
        cmd.extend(["--results", str(results)])

    if cleanup:
        cmd.append("--cleanup")

    # Add any extra parameters passed as database overrides
    for key, value in extra_params.items():
        if value is not None:
            # Convert Python arg names (gottcha2_db) to Nextflow params (--gottcha2-db)
            param_name = key.replace("_", "-")
            cmd.extend([f"--{param_name}", str(value)])

    return cmd


def run(  # noqa: PLR0913, PLR0912, PLR0915, C901
    # Complexity/boolean/B008 warnings are acceptable for CLI with many options
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
    experiment_id: str = typer.Option(
        ...,
        "--experiment-id",
        "-e",
        help="Unique experiment identifier",
        rich_help_panel=PANEL_CORE,
    ),
    tools: str = typer.Option(
        "all",
        "--tools",
        "-t",
        help="Tools to run: stat_blast, gottcha, or all",
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
        help="Custom config file (default: ~/.nvd2/config/user.config)",
        exists=True,
        rich_help_panel=PANEL_CORE,
    ),
    resume: bool = typer.Option(
        False,
        "--resume",
        help="Resume from checkpoint",
        rich_help_panel=PANEL_CORE,
    ),
    cleanup: bool = typer.Option(
        False,
        "--cleanup",
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
        help="State directory for run tracking and upload deduplication (default: ~/.cache/nvd or NVD_STATE_DIR)",
        rich_help_panel=PANEL_CORE,
    ),
    preset: str | None = typer.Option(
        None,
        "--preset",
        "-P",
        help="Use a registered preset (CLI args override preset values)",
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
    preprocess: bool = typer.Option(
        False,
        "--preprocess",
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
    labkey: bool = typer.Option(
        False,
        "--labkey",
        help="Enable LabKey integration",
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

    This command wraps 'nextflow run dhoconno/nvd' with a simpler interface.
    Database paths and settings are loaded from ~/.nvd2/config/user.config
    unless overridden with command-line options.

    Examples:

        # Simple run with auto-detected profile
        nvd run -s samples.csv -e exp001

        # Specify profile and results directory
        nvd run -s samples.csv -e exp002 -p docker -r ./my_results

        # Run only GOTTCHA2 workflow
        nvd run -s samples.csv -e exp003 -t gottcha

        # Resume a failed run
        nvd run -s samples.csv -e exp003 --resume

        # Use a preset (CLI args override preset values)
        nvd run -s samples.csv -e exp004 --preset production

        # Use a preset with overrides
        nvd run -s samples.csv -e exp005 --preset production --cutoff-percent 0.01
    """

    # Load preset if specified
    preset_params: dict[str, str | int | float | bool] = {}
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

    # Validate tools option
    if tools not in VALID_TOOLS:
        error(
            f"Invalid tools option: {tools}. Must be one of: {', '.join(VALID_TOOLS)}",
        )

    # Auto-detect profile if not specified
    if not profile:
        profile = auto_detect_profile()
        info(f"Auto-detected execution profile: {profile}")
    elif profile not in VALID_PROFILES:
        valid_profiles = ", ".join(VALID_PROFILES)
        error(f"Invalid profile: {profile}. Must be one of: {valid_profiles}")

    # Find config file
    config_file = find_config_file(config)
    if config_file:
        info(f"Using config: {config_file}")
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

    # Start with preset params (if any), CLI values will override
    extra_params: dict[str, str | Path | float | bool | int] = {}
    for key, value in preset_params.items():
        # Convert booleans to Nextflow string format
        if isinstance(value, bool):
            extra_params[key] = "true" if value else "false"
        else:
            extra_params[key] = value

    # CLI overrides (these take precedence over preset values)
    if gottcha2_db:
        extra_params["gottcha2_db"] = gottcha2_db
    if blast_db:
        extra_params["blast_db"] = blast_db
    if blast_db_prefix:
        extra_params["blast_db_prefix"] = blast_db_prefix
    if stat_index:
        extra_params["stat_index"] = stat_index
    if stat_dbss:
        extra_params["stat_dbss"] = stat_dbss
    if stat_annotation:
        extra_params["stat_annotation"] = stat_annotation
    if human_virus_taxlist:
        extra_params["human_virus_taxlist"] = human_virus_taxlist
    if cutoff_percent is not None:
        extra_params["cutoff_percent"] = cutoff_percent
    if tax_stringency is not None:
        extra_params["tax_stringency"] = tax_stringency
    if entropy is not None:
        extra_params["entropy"] = entropy
    if min_gottcha_reads is not None:
        extra_params["min_gottcha_reads"] = min_gottcha_reads
    if max_blast_targets is not None:
        extra_params["max_blast_targets"] = max_blast_targets
    if blast_retention_count is not None:
        extra_params["blast_retention_count"] = blast_retention_count
    if min_consecutive_bases is not None:
        extra_params["min_consecutive_bases"] = min_consecutive_bases
    if qtrim is not None:
        extra_params["qtrim"] = qtrim
    if include_children is not None:
        extra_params["include_children"] = "true" if include_children else "false"
    if max_concurrent_downloads is not None:
        extra_params["max_concurrent_downloads"] = max_concurrent_downloads
    # Preprocessing options
    if preprocess:
        extra_params["preprocess"] = "true"
    if merge_pairs is not None:
        extra_params["merge_pairs"] = "true" if merge_pairs else "false"
    if dedup is not None:
        extra_params["dedup"] = "true" if dedup else "false"
    if trim_adapters is not None:
        extra_params["trim_adapters"] = "true" if trim_adapters else "false"
    if scrub_host_reads is not None:
        extra_params["scrub_host_reads"] = "true" if scrub_host_reads else "false"
    if hostile_index:
        extra_params["hostile_index"] = hostile_index
    if hostile_index_name:
        extra_params["hostile_index_name"] = hostile_index_name
    if filter_reads is not None:
        extra_params["filter_reads"] = "true" if filter_reads else "false"
    if min_read_quality_illumina is not None:
        extra_params["min_read_quality_illumina"] = min_read_quality_illumina
    if min_read_quality_nanopore is not None:
        extra_params["min_read_quality_nanopore"] = min_read_quality_nanopore
    if min_read_length is not None:
        extra_params["min_read_length"] = min_read_length
    if max_read_length is not None:
        extra_params["max_read_length"] = max_read_length
    # SRA submission options
    # Handle deprecated --human-read-scrub option
    if human_read_scrub and not sra_human_db:
        warning(
            "DEPRECATION: --human-read-scrub is deprecated. "
            "Please use --sra-human-db instead.",
        )
        extra_params["sra_human_db"] = human_read_scrub
    elif sra_human_db:
        extra_params["sra_human_db"] = sra_human_db
    if labkey:
        extra_params["labkey"] = "true"
    # State directory for run tracking and upload deduplication
    if state_dir:
        extra_params["state_dir"] = state_dir

    # Build command
    cmd = build_nextflow_command(
        samplesheet=samplesheet,
        experiment_id=experiment_id,
        tools=tools,
        profile=profile,
        config=config_file,
        results=results,
        resume=resume,
        cleanup=cleanup,
        work_dir=work_dir,
        **extra_params,
    )

    # Show command
    console.print("\n[bold]Executing command:[/bold]")
    console.print(f"[dim]{' '.join(cmd)}[/dim]\n")

    if dry_run:
        success("Dry-run mode: command shown above but not executed")
        return

    # Execute nextflow
    try:
        result = subprocess.run(cmd, check=False)  # noqa: S603
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        console.print("\n[yellow]Pipeline interrupted by user[/yellow]")
        sys.exit(130)
    except OSError as e:
        error(f"Failed to execute Nextflow: {e}")
