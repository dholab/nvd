"""
Interactive prompts with timeout support for NVD CLI.

Provides Rich-styled prompts with configurable timeout, suitable for
samplesheet generation and other interactive workflows.

The timeout mechanism uses select.select() on Unix systems. On Windows
or non-TTY environments, prompts fall back to blocking input.

Usage:
    from py_nvd.cli.prompts import prompt_with_timeout, confirm_with_timeout
    from py_nvd.cli.utils import console

    # String prompt with choices
    platform = prompt_with_timeout(
        console,
        "Sequencing platform",
        choices=["illumina", "ont"],
        default="illumina",
    )

    # Confirmation prompt
    if confirm_with_timeout(console, "Write samplesheet?", default=True):
        write_file()
"""

from __future__ import annotations

import select
import sys
from typing import TYPE_CHECKING

from rich.text import Text

if TYPE_CHECKING:
    from rich.console import Console

# Default timeout for all interactive prompts (5 minutes)
PROMPT_TIMEOUT_SECONDS = 300


class PromptTimeout(Exception):
    """
    Raised when a prompt times out waiting for user input.

    Attributes:
        timeout_seconds: The timeout duration that was exceeded
    """

    def __init__(self, timeout_seconds: int):
        self.timeout_seconds = timeout_seconds
        super().__init__(f"Prompt timed out after {timeout_seconds} seconds")


def _read_with_timeout(timeout_seconds: int) -> str:
    """
    Read a line from stdin with timeout.

    On Windows or non-TTY, falls back to blocking read.

    Args:
        timeout_seconds: Maximum seconds to wait for input

    Returns:
        User input string (without trailing newline)

    Raises:
        PromptTimeout: If no input received within timeout
        EOFError: If stdin is closed
    """
    # Windows doesn't support select on stdin; non-TTY can't timeout
    if sys.platform == "win32" or not sys.stdin.isatty():
        return input()

    ready, _, _ = select.select([sys.stdin], [], [], timeout_seconds)
    if not ready:
        raise PromptTimeout(timeout_seconds)
    return input()


def prompt_with_timeout(
    console: Console,
    prompt: str,
    *,
    choices: list[str] | None = None,
    default: str | None = None,
    timeout_seconds: int = PROMPT_TIMEOUT_SECONDS,
) -> str:
    """
    Display a Rich-styled prompt with timeout.

    Prompts the user for input with optional choice validation. If the user
    doesn't respond within the timeout, raises PromptTimeout.

    Args:
        console: Rich Console instance for styled output
        prompt: Prompt message to display
        choices: Optional list of valid choices (case-insensitive matching)
        default: Default value if user presses Enter without input
        timeout_seconds: Maximum seconds to wait (default: 300)

    Returns:
        User input string, or default if Enter pressed

    Raises:
        PromptTimeout: If no input received within timeout
        EOFError: If stdin is closed

    Example:
        >>> platform = prompt_with_timeout(
        ...     console,
        ...     "Sequencing platform",
        ...     choices=["illumina", "ont"],
        ...     default="illumina",
        ... )
    """
    while True:
        # Build styled prompt text
        prompt_text = Text(prompt, style="bold")

        if choices:
            prompt_text.append(" ")
            prompt_text.append(f"[{'/'.join(choices)}]", style="cyan")

        if default is not None:
            prompt_text.append(" ")
            prompt_text.append(f"({default})", style="dim")

        prompt_text.append(": ")

        console.print(prompt_text, end="")

        try:
            response = _read_with_timeout(timeout_seconds).strip()
        except EOFError:
            console.print()
            raise

        # Handle empty response with default
        if not response and default is not None:
            return default

        # Handle empty response without default
        if not response:
            console.print("[red]Input required.[/red]")
            continue

        # Validate against choices if provided
        if choices:
            # Case-insensitive matching, return canonical form
            for choice in choices:
                if choice.lower() == response.lower():
                    return choice
            console.print(
                f"[red]Invalid choice.[/red] Please select from: {', '.join(choices)}"
            )
            continue

        return response


def confirm_with_timeout(
    console: Console,
    prompt: str,
    *,
    default: bool = False,
    timeout_seconds: int = PROMPT_TIMEOUT_SECONDS,
) -> bool:
    """
    Display a Rich-styled confirmation prompt with timeout.

    Prompts the user for yes/no confirmation. If the user doesn't respond
    within the timeout, raises PromptTimeout.

    Args:
        console: Rich Console instance for styled output
        prompt: Confirmation message to display
        default: Default value if user presses Enter (True=yes, False=no)
        timeout_seconds: Maximum seconds to wait (default: 300)

    Returns:
        True if user confirmed, False otherwise

    Raises:
        PromptTimeout: If no input received within timeout
        EOFError: If stdin is closed

    Example:
        >>> if confirm_with_timeout(console, "Write samplesheet?", default=True):
        ...     write_file()
    """
    while True:
        # Build styled prompt text
        prompt_text = Text(prompt, style="bold")
        prompt_text.append(" ")

        if default:
            prompt_text.append("[Y/n]", style="cyan")
        else:
            prompt_text.append("[y/N]", style="cyan")

        prompt_text.append(": ")

        console.print(prompt_text, end="")

        try:
            response = _read_with_timeout(timeout_seconds).strip().lower()
        except EOFError:
            console.print()
            return False

        # Handle empty response with default
        if not response:
            return default

        # Parse yes/no
        if response in ("y", "yes"):
            return True
        if response in ("n", "no"):
            return False

        console.print("[red]Please enter 'y' or 'n'[/red]")


def is_interactive() -> bool:
    """
    Check if the current environment supports interactive prompts.

    Returns True if stdin is a TTY (terminal), False otherwise.
    Use this to decide whether to prompt or require CLI flags.

    Returns:
        True if interactive prompts are possible
    """
    return sys.stdin.isatty()
