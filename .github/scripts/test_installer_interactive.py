#!/usr/bin/env python3
"""Regression tests for piped interactive installer prompts.

These tests intentionally run install.sh the same way users do with
`curl ... | bash`: the script body is piped into bash, while answers are typed
through a controlling pseudo-terminal. Plain `printf ... | bash install.sh` is
not an equivalent test because it puts answers on stdin, which is already used
for the piped script body in the real install path.
"""

from __future__ import annotations

import errno
import os
import select
import shutil
import stat
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INSTALLER = ROOT / "install.sh"


class PtySession:
    def __init__(
        self, command: str, env: dict[str, str], timeout: float = 20.0
    ) -> None:
        self.command = command
        self.env = env
        self.timeout = timeout
        self.output = bytearray()
        self.pid: int | None = None
        self.fd: int | None = None

    def __enter__(self) -> PtySession:
        pid, fd = os.forkpty()
        if pid == 0:
            os.execvpe("bash", ["bash", "-c", self.command], self.env)
        self.pid = pid
        self.fd = fd
        os.set_blocking(fd, False)
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        if self.fd is not None:
            try:
                os.close(self.fd)
            except OSError:
                pass
        if self.pid is not None:
            try:
                waited_pid, status = os.waitpid(self.pid, os.WNOHANG)
                if waited_pid == 0:
                    os.kill(self.pid, 15)
                    os.waitpid(self.pid, 0)
                elif status != 0 and exc_type is None:
                    raise AssertionError(
                        f"command exited non-zero: {status}\n{self.text}"
                    )
            except ChildProcessError:
                pass

    @property
    def text(self) -> str:
        return self.output.decode("utf-8", errors="replace")

    def _read_available(self) -> None:
        assert self.fd is not None
        while True:
            try:
                chunk = os.read(self.fd, 4096)
            except BlockingIOError:
                return
            except OSError as e:
                if e.errno == errno.EIO:
                    return
                raise
            if not chunk:
                return
            self.output.extend(chunk)

    def expect(self, needle: str) -> None:
        deadline = time.monotonic() + self.timeout
        encoded = needle.encode()
        while time.monotonic() < deadline:
            self._read_available()
            if encoded in self.output:
                return
            time.sleep(0.05)
        raise AssertionError(f"timed out waiting for {needle!r}\n{self.text}")

    def send(self, text: str) -> None:
        assert self.fd is not None
        os.write(self.fd, text.encode())

    def wait(self) -> int:
        assert self.pid is not None
        deadline = time.monotonic() + self.timeout
        while time.monotonic() < deadline:
            self._read_available()
            waited_pid, status = os.waitpid(self.pid, os.WNOHANG)
            if waited_pid != 0:
                self.pid = None
                return status
            time.sleep(0.05)
        raise AssertionError(f"timed out waiting for command exit\n{self.text}")


def write_executable(path: Path, content: str) -> None:
    path.write_text(textwrap.dedent(content).lstrip())
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def make_fake_install(home: Path) -> Path:
    fake_bin = home / "fake-bin"
    fake_bin.mkdir(parents=True)

    install_root = home / ".nvd"
    versioned_repo = install_root / "v3.0.0"
    nvd_bin_dir = versioned_repo / ".pixi" / "envs" / "default" / "bin"
    nvd_bin_dir.mkdir(parents=True)
    (install_root / "latest").symlink_to("v3.0.0")

    write_executable(
        fake_bin / "git",
        """
        #!/usr/bin/env bash
        if [[ "$1" == "--version" ]]; then
          echo "git version 2.52.0"
          exit 0
        fi
        echo "UNEXPECTED_GIT_COMMAND:$*" >&2
        exit 70
        """,
    )
    write_executable(
        fake_bin / "pixi",
        """
        #!/usr/bin/env bash
        echo "pixi 0.68.1"
        """,
    )
    write_executable(
        fake_bin / "java",
        """
        #!/usr/bin/env bash
        echo 'openjdk version "17.0.0"' >&2
        """,
    )
    write_executable(
        fake_bin / "nextflow",
        """
        #!/usr/bin/env bash
        echo "Nextflow version 25.10.4"
        """,
    )
    write_executable(
        fake_bin / "apptainer",
        """
        #!/usr/bin/env bash
        echo "apptainer version 1.5.0"
        """,
    )
    write_executable(
        nvd_bin_dir / "nvd",
        """
        #!/usr/bin/env bash
        if [[ "$1" == "setup" && "${2:-}" == "--help" ]]; then
          exit 0
        fi
        if [[ "$1" == "setup" ]]; then
          if [[ -t 0 ]]; then
            echo "FAKE_NVD_STDIN_TTY=yes"
          else
            echo "FAKE_NVD_STDIN_TTY=no"
            exit 43
          fi
          printf 'Config directory [/tmp/default]: '
          if ! IFS= read -r answer; then
            echo "FAKE_NVD_EOF"
            exit 42
          fi
          echo "FAKE_NVD_READ=${answer}"
          echo "✓ Setup complete!"
          exit 0
        fi
        echo "unexpected nvd command: $*" >&2
        exit 64
        """,
    )

    return fake_bin


def base_env(home: Path, fake_bin: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.update(
        {
            "HOME": str(home),
            "PATH": f"{fake_bin}{os.pathsep}{env['PATH']}",
            "SHELL": "/bin/bash",
            "NVD_INSTALLER_DEBUG": "1",
        },
    )
    return env


def run_piped_installer(env: dict[str, str], *args: str) -> PtySession:
    quoted_args = " ".join(subprocess.list2cmdline([arg]) for arg in args)
    command = (
        "export NVD_INSTALLER_DEBUG=1; "
        f"/bin/cat {subprocess.list2cmdline([str(INSTALLER)])} | bash -s -- {quoted_args}"
    )
    return PtySession(command, env)


def test_existing_install_can_decline_update_and_setup_reads_tty() -> None:
    tmp = Path(tempfile.mkdtemp(prefix="nvd-installer-"))
    try:
        fake_bin = make_fake_install(tmp)
        with run_piped_installer(base_env(tmp, fake_bin)) as session:
            session.expect("Pull latest updates? [Y/n]:")
            session.send("n\r")
            session.expect("[debug] read_prompt: bytes=1 value=n")
            session.expect("[debug] prompt_yes_no: no")
            session.expect("FAKE_NVD_STDIN_TTY=yes")
            session.expect("Config directory [/tmp/default]:")
            session.send(f"{tmp}/config\r")
            session.expect(f"FAKE_NVD_READ={tmp}/config")
            session.expect("Download reference databases? [y/N]:")
            session.send("n\r")
            session.expect("Installation Complete")
            status = session.wait()
        output = session.text
        assert status == 0, output
        assert "Pulling updates" not in output
        assert "UNEXPECTED_GIT_COMMAND" not in output
        assert "FAKE_NVD_EOF" not in output
    finally:
        shutil.rmtree(tmp)


def test_dry_run_database_wizard_walks_planned_downloads() -> None:
    tmp = Path(tempfile.mkdtemp(prefix="nvd-installer-"))
    try:
        fake_bin = make_fake_install(tmp)
        with run_piped_installer(base_env(tmp, fake_bin), "--dry-run") as session:
            session.expect("Pull latest updates? [Y/n]:")
            session.send("n\r")
            session.expect("[debug] prompt_yes_no: no")
            session.expect("Download reference databases? [y/N]:")
            session.send("y\r")
            session.expect("[debug] prompt_yes_no: yes")
            session.expect("Database Setup")
            session.expect("Choice [1-4]:")
            session.send("3\r")
            session.expect("Where should references be stored?")
            session.send("\r")
            session.expect("[DRY RUN] Would check for approximately 501GB free")
            session.expect("[DRY RUN] Would download checksum manifest")
            session.expect("[DRY RUN] Would download BLAST database archive")
            session.expect("[DRY RUN] Would extract BLAST database")
            session.expect("Remove BLAST archive to save disk space? [Y/n]:")
            session.send("n\r")
            session.expect("[DRY RUN] Would download deacon vertebrate-virus index")
            session.expect("Reference Paths")
            session.expect("Installation Complete")
            status = session.wait()
        output = session.text
        assert status == 0, output
        assert "Downloading BLAST database archive" not in output
        assert "[DRY RUN] Would remove" not in output
        assert not (tmp / ".nvd" / "references" / "blast_db_v3_0.tar.gz").exists()
        assert not (
            tmp / ".nvd" / "references" / "human_infecting_viruses.k31w1.idx"
        ).exists()
    finally:
        shutil.rmtree(tmp)


def main() -> None:
    test_existing_install_can_decline_update_and_setup_reads_tty()
    test_dry_run_database_wizard_walks_planned_downloads()
    print("installer interactive prompt tests passed")


if __name__ == "__main__":
    main()
