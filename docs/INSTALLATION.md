# NVD2 Installation Guide

## Quick Start

NVD2 includes an interactive setup script to help configure database paths:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash
```

Or download and inspect first:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh -o install.sh
chmod +x install.sh
./install.sh
```

**Prerequisites** (must be installed separately):
- Java 11 or newer
- Nextflow
- Docker, Apptainer/Singularity, or Pixi

## Installation Modes

### Interactive Mode (Default)

Guides you through the complete setup process:

```bash
./install.sh
```

**What it does:**
- Checks that prerequisites are installed
- Detects available execution environments (Docker, Apptainer)
- Helps configure reference database paths
- Optionally downloads databases (100s of GB)
- Creates configuration file at `~/.nvd2/config/user.config`

**Time:** 5-10 minutes + database downloads (if selected)

### Dry-Run Mode

Test the installer without making any changes:

```bash
./install.sh --dry-run
```

**Use cases:**
- Preview what the installer will do
- Test on a new system before committing
- Verify the installation flow

### Non-Interactive Mode (CI/CD)

Check prerequisites only, no prompts or modifications:

```bash
./install.sh --non-interactive
```

**Exit codes:**
- `0` - All prerequisites met
- `1` - Java missing or wrong version
- `2` - Nextflow missing
- `3` - No execution environment available

### Verify Mode

Check an existing installation:

```bash
./install.sh --verify
```

**What it checks:**
- Dependencies are installed and correct versions
- Configuration file exists and is valid
- Database paths are correct and files exist
- Reports any issues found

### Uninstall Mode

Remove NVD2 from your system:

```bash
./install.sh --uninstall
```

**What it removes:**
- Configuration files (`~/.nvd2/`)
- Cached pipeline (`~/.nextflow/assets/dhoconno/nvd`)
- Optionally: databases (asks for confirmation)
- Shows commands for removing dependencies

## System Requirements

### Required Dependencies

1. **Java 11 or newer**
   - OpenJDK or Oracle JDK
   - Required by Nextflow

2. **Nextflow**
   - Latest version recommended
   - Pipeline execution engine

3. **Execution Environment** (choose one):
   - **Docker** (recommended for workstations)
   - **Podman** (rootless Docker alternative)
   - **Apptainer/Singularity** (HPC environments)
   - **Pixi** (local execution, GOTTCHA2 only)

### Disk Space Requirements

**For databases:**
- STAT database: ~500GB
- BLAST database: ~500GB
- GOTTCHA2 database: ~500GB

**Plus 20% safety buffer**

You can configure databases on separate volumes or skip download entirely.

## Manual Configuration

If you prefer to configure NVD2 manually without using `install.sh`:

### 1. Ensure Prerequisites are Installed

You must have these installed:
- **Java 11+**: OpenJDK or Oracle JDK (required by Nextflow)
- **Nextflow**: https://www.nextflow.io/docs/latest/getstarted.html
- **Container runtime** (choose one):
  - Docker: https://docs.docker.com/get-docker/
  - Apptainer/Singularity: https://apptainer.org/docs/admin/main/installation.html
  - Pixi: https://pixi.sh (for local Conda-based execution)

### 2. Configure Databases

Create a configuration file at `~/.nvd2/config/user.config`:

```groovy
params {
    // STAT database
    stat_index      = "/path/to/databases/stat_db/tree_index.dbs"
    stat_dbss       = "/path/to/databases/stat_db/tree_filter.dbss"
    stat_annotation = "/path/to/databases/stat_db/tree_filter.dbss.annotation"
    human_virus_taxlist = "/path/to/databases/stat_db/human_viruses_taxlist.txt"
    
    // BLAST database
    blast_db        = "/path/to/databases/blast_db"
    blast_db_prefix = "core_nt"
    
    // GOTTCHA2 database
    gottcha2_db     = "/path/to/databases/gottcha2_db/gottcha_db.species.fna"
}
```

## Database Setup

### Option 1: Download via Installer

The installer can download and extract databases for you:

```bash
./install.sh
# Choose "Help me set up database downloads" when prompted
```

### Option 2: Manual Download

Databases are available from:
```
https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/@files/release-v2.0.0/
```

Files:
- `stat_db.tar.gz` (MD5: 2641a1754f6986eedbe9f38e188f2b0c)
- `blast_db.tar.gz` (MD5: cebc4c59ea572c09fb93aa4e3594bf3e)
- `gottcha2.tar.gz` (MD5: d33fb5d1b2d22f7a174239f1dfc142cb)

Extract to your desired location and update the configuration file.

## Running NVD2

After installation, you can run NVD2 using either the CLI wrapper or direct Nextflow execution.

### Option 1: CLI Wrapper

The `nvd` CLI wrapper provides a simplified interface:

```bash
# Clone repository and install dependencies
git clone https://github.com/dhoconno/nvd.git
cd nvd
uv sync                 # Python dependencies only
# or
pixi shell --frozen     # All Conda dependencies

# Run with auto-configuration
nvd run --samplesheet samples.csv --experiment-id exp001

# Specify tools and profile
nvd run -s samples.csv -e exp002 --tools gottcha -p docker

# Get help
nvd --help
```

**Note**: When not in an active shell, prepend commands with `pixi run` (for full toolchain) 
or `uv run` (Python-only scripts).

The CLI wrapper automatically:
- Loads configuration from `~/.nvd2/config/user.config`
- Auto-detects your execution profile (Docker/Apptainer/local)
- Validates inputs before running

See the **[CLI Wrapper Guide](./cli_wrapper_guide.md)** for comprehensive examples.

### Option 2: Direct Nextflow Execution

You can also run NVD2 with direct Nextflow commands:

```bash
nextflow run dhoconno/nvd \
  -c ~/.nvd2/config/user.config \
  -profile docker \
  --samplesheet samples.csv \
  --experiment_id exp001 \
  --tools all
```

**Execution profiles:**
- `-profile docker` - Use Docker containers (local workstations)
- `-profile apptainer` - Use Apptainer/Singularity (HPC clusters)
- `-profile local` - Local execution without containers

**Tool options:**
- `--tools stat_blast` - STAT + BLAST workflow (human virus detection)
- `--tools gottcha` - GOTTCHA2 workflow (general classification)
- `--tools all` - Run both workflows

See **[Direct Nextflow Examples](./example_commands.md)** for more traditional usage patterns.

## Troubleshooting

### Docker Not Running

**Error:** "Docker is installed but not running"

**Solution:**
- macOS: Open Docker Desktop and wait for it to start
- Linux: `sudo systemctl start docker`

### Insufficient Disk Space

**Error:** "Insufficient disk space"

**Solutions:**
1. Choose a different location with more space
2. Download fewer databases
3. Skip database download and configure manually later

### Java Version Too Old

**Error:** "Java version too old"

**Solution:** Install Java 11 or newer. See https://adoptium.net/ or your system's package manager.

### Network Issues During Download

Downloads are resumable. If interrupted:
1. Re-run the installer
2. Choose the same database location
3. Download will resume from where it left off

### Permission Errors

**Error:** "Cannot write to /path"

**Solutions:**
- Choose a different path you have write access to
- Create the directory first: `mkdir -p /path && chmod 755 /path`
- Use sudo (not recommended): `sudo ./install.sh`

## HPC-Specific Notes

### Module Systems

If your HPC uses modules:

```bash
module load java/17
module load nextflow
module load apptainer
```

Run installer after loading modules:
```bash
./install.sh --non-interactive  # Check if ready
```

### Shared Installations

For system-wide installation:
1. Admin installs to shared location (e.g., `/opt/nvd2/`)
2. Users create personal config pointing to shared databases
3. Each user runs: `./install.sh` and configures paths

### Scratch vs Permanent Storage

**Databases:**
- Install to permanent storage (not /scratch)
- Databases are 1+ TB total

**Working directory:**
- Can use scratch space for pipeline execution
- Set with `-work-dir` flag when running Nextflow

## Advanced Options

### Custom Installation Path

```bash
# Pipeline source
./install.sh
# Choose option 4: "Clone to custom directory"

# Databases
# Specify custom paths when prompted
```

### Offline Installation

1. Download installer and databases on internet-connected machine
2. Transfer to offline system
3. Run: `./install.sh`
4. Choose "I already have the databases"
5. Point to transferred database locations

### Multiple Versions

Install multiple pipeline versions:

```bash
# Version 1
nextflow pull dhoconno/nvd -r v1.0.0

# Version 2  
nextflow pull dhoconno/nvd -r v2.0.0

# Run specific version
nextflow run dhoconno/nvd -r v1.0.0 ...
```

## Getting Help

**CLI Wrapper Guide:** [cli_wrapper_guide.md](./cli_wrapper_guide.md) - Modern CLI interface

**Documentation:** https://github.com/dhoconno/nvd

**Issues:** https://github.com/dhoconno/nvd/issues

**Validation Commands:**
```bash
./install.sh --verify        # Verify installation
nvd validate all             # Validate complete setup with CLI
nvd --help                   # CLI help
```

## Platform-Specific Notes

### macOS

- Apple Silicon (M1/M2/M3): Fully supported
- Intel: Fully supported
- Docker Desktop required (no native Docker daemon on macOS)

### Linux

- Ubuntu 22.04+: Fully tested
- Debian 11+: Compatible
- RHEL/CentOS 8+: Compatible
- Arch: Compatible (install java-openjdk)

### Windows

- WSL2 required
- Follow Linux instructions within WSL2
- Docker Desktop for Windows with WSL2 backend

## Security Considerations

### curl | bash Pattern

The quick start uses `curl | bash`. To review first:

```bash
curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh -o install.sh
less install.sh  # Review the script
bash install.sh
```

### Permissions

The installer:
- Never requires root/sudo
- Only writes to `~/.nvd2/` and user-specified database locations
- Does not modify system files
- Does not install dependencies automatically (shows instructions only)

### Database Integrity

All database downloads are verified with MD5 checksums.

## FAQ

**Q: How long does installation take?**  
A: 5-10 minutes for setup. Database downloads: 2-6 hours depending on connection.

**Q: Can I install without downloading databases?**  
A: Yes! Choose "Skip" when prompted. Configure database paths manually later.

**Q: Do I need Docker?**  
A: Docker is recommended but not required. You can use Podman, Apptainer, or Pixi.

**Q: Can I install databases on a different drive?**  
A: Yes! Specify any path when configuring database locations.

**Q: How do I update NVD2?**  
A: `nextflow pull dhoconno/nvd` updates the pipeline. Databases are versioned separately.

**Q: How do I uninstall?**  
A: Run `./install.sh --uninstall` and follow the prompts.

**Q: What if installation fails?**  
A: Run `./install.sh --verify` to see what's missing. Check the error log at `~/.nvd2/install-error.log`.

## License

NVD2 installer is licensed under GPLv3. See LICENSE file for details.
