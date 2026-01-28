"""
Resume command for the NVD CLI.

Commands:
    nvd resume  - Resume a previous pipeline run from .nfresume
"""

from __future__ import annotations

import shlex
import subprocess
import sys
import time
from typing import Annotated

import typer

from py_nvd.cli.prompts import (
    PromptTimeout,
    confirm_with_timeout,
    is_interactive,
)
from py_nvd.cli.utils import (
    RESUME_FILE,
    console,
    error,
    get_editor,
    info,
    warning,
)

# If editor exits faster than this, warn about GUI editors needing --wait
EDITOR_MIN_DURATION_SECONDS = 1.0


def _read_resume_command() -> str:
    """
    Read the resume command from .nfresume file.

    Returns:
        The command string from the resume file.

    Raises:
        SystemExit: If no resume file exists.
    """
    if not RESUME_FILE.exists():
        error(
            "No previous run detected (.nfresume not found).\n"
            "Make sure you start with 'nvd run' before using 'nvd resume'.",
        )
    return RESUME_FILE.read_text(encoding="utf-8").strip()


def _execute_command(command: str) -> None:
    """
    Execute a shell command string.

    Args:
        command: The command string to execute.

    Raises:
        SystemExit: With the command's return code, or 130 on interrupt.
    """
    try:
        split_command = shlex.split(command)
        result = subprocess.run(split_command, check=False)  # noqa: S603
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        console.print("\n[yellow]Pipeline interrupted by user[/yellow]")
        sys.exit(130)
    except OSError as e:
        error(f"Failed to execute Nextflow: {e}")


def _open_editor(filepath: str) -> bool:
    """
    Open a file in the user's preferred editor.

    Uses $VISUAL, then $EDITOR, then falls back to vi.

    Note: GUI editors (VS Code, Sublime, etc.) may return immediately
    unless configured with a --wait flag. If the editor exits in less
    than 1 second, we warn the user about this.

    Args:
        filepath: Path to the file to edit.

    Returns:
        True if editor exited successfully, False otherwise.
    """
    editor = get_editor()
    info(f"Opening {filepath} in {editor}...")

    try:
        start_time = time.monotonic()
        result = subprocess.run([editor, filepath], check=False)  # noqa: S603
        elapsed = time.monotonic() - start_time

        # Warn if editor exited suspiciously fast (likely a GUI editor without --wait)
        if elapsed < EDITOR_MIN_DURATION_SECONDS:
            warning(
                f"Editor exited in {elapsed:.1f}s - did it open in the background?\n"
                "   GUI editors like VS Code need: export EDITOR='code --wait'\n"
                "   Or use a terminal editor like vim, nano, or helix.",
            )

        return result.returncode == 0
    except FileNotFoundError:
        warning(f"Editor '{editor}' not found. Set $EDITOR to your preferred editor.")
        return False
    except OSError as e:
        warning(f"Failed to open editor: {e}")
        return False


def _interactive_resume(command: str) -> None:
    """
    Handle the interactive resume flow with confirmation and optional editing.

    Flow:
    1. Show command, prompt to run
    2. If no: prompt to edit
    3. If edit: open editor, then prompt to run (no second edit chance)
    4. Execute or exit based on responses

    Args:
        command: The initial command from .nfresume

    Raises:
        SystemExit: On timeout, user decline, or after execution.
    """
    console.print("\n[bold]Resume command:[/bold]")
    console.print(f"[dim]{command}[/dim]\n")

    try:
        # First confirmation: run as-is?
        if confirm_with_timeout(console, "Run this command?", default=True):
            _execute_command(command)
            return  # _execute_command calls sys.exit, but be explicit

        # User said no - offer to edit
        if not confirm_with_timeout(console, "Edit the command?", default=False):
            info("Resume cancelled.")
            raise typer.Exit(0)

        # Open editor
        if not _open_editor(str(RESUME_FILE)):
            error("Editor failed. Resume cancelled.")

        # Re-read the (possibly modified) command
        edited_command = _read_resume_command()

        console.print("\n[bold]Edited command:[/bold]")
        console.print(f"[dim]{edited_command}[/dim]\n")

        # Final confirmation (no more edit chances)
        if confirm_with_timeout(console, "Run this command?", default=True):
            _execute_command(edited_command)
        else:
            info("Resume cancelled.")
            raise typer.Exit(0)

    except PromptTimeout:
        console.print("\n")
        warning("Prompt timed out after 5 minutes. Resume cancelled.")
        raise typer.Exit(1) from None


def resume(
    interactive: Annotated[
        bool,
        typer.Option(
            "--interactive",
            "-i",
            help="Prompt for confirmation before running; optionally edit the command",
        ),
    ] = False,
) -> None:
    """
    Resume a previous pipeline run.

    This command reads the cached Nextflow command from the .nfresume file
    (created during a previous 'nvd run') and re-executes it with the
    -resume flag to continue from where it left off.

    [dim]Note:[/dim] This is different from 'nvd run --resume', which
    requires you to re-specify all parameters. The 'resume' subcommand
    uses the exact same parameters from your previous run.

    [bold]Editor Configuration:[/bold]

    In interactive mode (-i), you can edit the command using your preferred
    editor ($VISUAL or $EDITOR). GUI editors like VS Code or Sublime Text
    need a --wait flag to block until the file is closed:

        export EDITOR='code --wait'      # VS Code
        export EDITOR='subl --wait'      # Sublime Text

    Terminal editors (vim, nano, helix, emacs -nw) work without modification.

    [dim]Examples:[/dim]

        # First, run the pipeline
        nvd run -s samples.csv -e exp001 -p docker

        # If it fails or is interrupted, resume with:
        nvd resume

        # Or review/edit the command before running:
        nvd resume -i
    """
    command = _read_resume_command()

    if interactive:
        if not is_interactive():
            error(
                "Interactive mode requires a terminal.\n"
                "Use 'nvd resume' without -i for non-interactive execution.",
            )
        _interactive_resume(command)
    else:
        info(f"Resuming: {command}")
        _execute_command(command)
