#!/usr/bin/env bash
#
# NVD Installer
#
# Bootstrap script that:
#   1. Checks dependencies (pixi, Java, Nextflow, container engine)
#   2. Clones the NVD repo to a versioned directory
#   3. Runs `nvd setup` for configuration
#   4. Optionally downloads reference databases
#
# Usage:
#   curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash
#   OR
#   ./install.sh [OPTIONS]
#
# License: GPLv3
# Project: https://github.com/dhoconno/nvd
#

set -euo pipefail

# Installer script version (semver)
# MAJOR: Breaking changes to CLI flags, directory structure, or exit codes
# MINOR: New features, new flags, new modes
# PATCH: Bug fixes, cross-platform improvements, better error messages
readonly VERSION="2.0.0"

# Cleanup handler for interrupts - removes temporary files
cleanup() {
	local exit_code=$?
	# Remove temp clone directory if it exists (don't remove partial downloads - they can resume)
	if [[ -n "${NVD_DIR:-}" ]] && [[ -d "${NVD_DIR}/.clone-temp" ]]; then
		rm -rf "${NVD_DIR}/.clone-temp"
	fi
	exit $exit_code
}
trap cleanup EXIT INT TERM

# =============================================================================
# Colors and Symbols
# =============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
DIM='\033[2m'
RESET='\033[0m'

CHECKMARK="✓"
CROSS="✗"
ARROW="→"

# ASCII fallbacks for non-UTF-8 terminals
if [[ ! "${LC_ALL:-${LC_CTYPE:-${LANG:-}}}" =~ UTF-8 ]]; then
	CHECKMARK="[OK]"
	CROSS="[X]"
	ARROW="->"
fi

# =============================================================================
# Print Utilities
# =============================================================================

print_header() {
	local message="$1"
	echo
	echo -e "${BOLD}${CYAN}${message}${RESET}"
	echo -e "${DIM}$(printf '%.0s─' {1..60})${RESET}"
}

print_success() {
	echo -e "  ${GREEN}${CHECKMARK}${RESET} $1"
}

print_error() {
	echo -e "  ${RED}${CROSS}${RESET} $1" >&2
}

print_warning() {
	echo -e "  ${YELLOW}!${RESET} $1"
}

print_info() {
	echo -e "  ${CYAN}${ARROW}${RESET} $1"
}

# =============================================================================
# Prompt Utilities
# =============================================================================

prompt_yes_no() {
	local prompt="$1"
	local default="${2:-y}"

	if [[ "$DRY_RUN" == "true" ]] || [[ "$NON_INTERACTIVE" == "true" ]]; then
		[[ "$default" == "y" ]]
		return
	fi

	local yn_hint
	if [[ "$default" == "y" ]]; then
		yn_hint="[Y/n]"
	else
		yn_hint="[y/N]"
	fi

	echo -en "  ${prompt} ${yn_hint}: "
	local response
	read -r response
	response="${response:-$default}"

	case "$response" in
	[Yy]*) return 0 ;;
	*) return 1 ;;
	esac
}

prompt_choice() {
	local prompt="$1"
	shift
	local options=("$@")

	# Display goes to stderr so it doesn't interfere with captured output
	echo >&2
	echo -e "  ${BOLD}${prompt}${RESET}" >&2
	echo >&2

	local i=1
	for opt in "${options[@]}"; do
		echo -e "    ${CYAN}[$i]${RESET} $opt" >&2
		((i++))
	done

	echo >&2

	if [[ "$DRY_RUN" == "true" ]] || [[ "$NON_INTERACTIVE" == "true" ]]; then
		echo "  Choice [1-${#options[@]}]: 1 (auto-selected)" >&2
		echo 1
		return
	fi

	while true; do
		echo -en "  Choice [1-${#options[@]}]: " >&2
		local choice
		read -r choice

		if [[ "$choice" =~ ^[0-9]+$ ]] && [[ "$choice" -ge 1 ]] && [[ "$choice" -le "${#options[@]}" ]]; then
			echo "$choice"
			return
		fi
		print_error "Invalid choice. Enter a number between 1 and ${#options[@]}"
	done
}

prompt_path() {
	local prompt="$1"
	local default="$2"

	if [[ "$DRY_RUN" == "true" ]] || [[ "$NON_INTERACTIVE" == "true" ]]; then
		echo "$default"
		return
	fi

	echo -en "  ${prompt} [${default}]: " >&2
	local response
	read -r response
	response="${response:-$default}"

	# Expand ~ and $HOME
	response="${response/#\~/$HOME}"
	response="${response/\$HOME/$HOME}"

	echo "$response"
}

# =============================================================================
# Dependency Checks
# =============================================================================

check_pixi() {
	if ! command -v pixi &>/dev/null; then
		print_error "pixi not found"
		echo
		echo "    NVD requires pixi for environment management."
		echo "    Install it with:"
		echo
		echo "      curl -fsSL https://pixi.sh/install.sh | bash"
		echo
		echo "    Then run this installer again."
		echo
		return 1
	fi
	local version
	version=$(pixi --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
	print_success "pixi: ${version}"
}

check_java() {
	if ! command -v java &>/dev/null; then
		print_error "Java not found"
		echo
		echo "    NVD requires Java 11 or newer for Nextflow."
		echo "    Install it with:"
		echo
		if [[ "$OSTYPE" == "darwin"* ]]; then
			echo "      brew install openjdk@17"
		elif command -v apt &>/dev/null; then
			echo "      sudo apt install openjdk-17-jdk"
		elif command -v yum &>/dev/null; then
			echo "      sudo yum install java-17-openjdk"
		else
			echo "      # Install OpenJDK 17 via your package manager"
		fi
		echo
		echo "    Or use SDKMAN: https://sdkman.io"
		echo
		return 1
	fi

	# Extract version - handles both "17.0.1" and "1.8.0_292" formats
	local version major
	version=$(java -version 2>&1 | head -1 | sed -n 's/.*version "\([0-9.]*\).*/\1/p')
	major=$(echo "$version" | cut -d. -f1)

	# Handle old 1.x versioning (1.8 = Java 8)
	if [[ "$major" == "1" ]]; then
		major=$(echo "$version" | cut -d. -f2)
	fi

	# Guard against empty/non-numeric major version
	if [[ -z "$major" ]] || [[ ! "$major" =~ ^[0-9]+$ ]]; then
		print_error "Could not parse Java version"
		echo
		echo "    java -version output: $(java -version 2>&1 | head -1)"
		echo
		return 1
	fi

	if [[ "$major" -lt 11 ]]; then
		print_error "Java version too old: ${version} (need 11+)"
		echo
		echo "    Please upgrade Java to version 11 or newer."
		echo
		return 1
	fi
	print_success "Java: ${version}"
}

check_nextflow() {
	if ! command -v nextflow &>/dev/null; then
		print_error "Nextflow not found"
		echo
		echo "    NVD requires Nextflow to run pipelines."
		echo "    Install it with:"
		echo
		echo "      curl -s https://get.nextflow.io | bash"
		echo "      sudo mv nextflow /usr/local/bin/"
		echo
		return 1
	fi
	local version
	version=$(nextflow -version 2>&1 | grep -i version | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
	print_success "Nextflow: ${version}"
}

check_git() {
	if ! command -v git &>/dev/null; then
		print_error "git not found"
		echo
		echo "    NVD requires git to clone the repository."
		echo "    Install it via your package manager:"
		echo
		if [[ "$OSTYPE" == "darwin"* ]]; then
			echo "      xcode-select --install"
		elif command -v apt &>/dev/null; then
			echo "      sudo apt install git"
		elif command -v yum &>/dev/null; then
			echo "      sudo yum install git"
		else
			echo "      # Install git via your package manager"
		fi
		echo
		return 1
	fi
	local version
	version=$(git --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
	print_success "git: ${version}"
}

check_tar() {
	if ! command -v tar &>/dev/null; then
		print_error "tar not found"
		echo
		echo "    NVD requires tar to extract archives."
		echo "    This should be pre-installed on most systems."
		echo
		return 1
	fi
	# tar --version output varies wildly between implementations, just confirm it exists
	print_success "tar: available"
}

check_container_engine() {
	local found=false

	# Check Docker (with timeout to avoid hanging if daemon is unresponsive)
	if command -v docker &>/dev/null; then
		# Portable timeout: try GNU timeout, then gtimeout (Homebrew), then fallback
		local docker_ok=false
		if command -v timeout &>/dev/null; then
			timeout 5 docker info &>/dev/null 2>&1 && docker_ok=true
		elif command -v gtimeout &>/dev/null; then
			gtimeout 5 docker info &>/dev/null 2>&1 && docker_ok=true
		else
			# No timeout available - try docker info directly (may hang on unresponsive daemon)
			docker info &>/dev/null 2>&1 && docker_ok=true
		fi

		if [[ "$docker_ok" == "true" ]]; then
			local version
			version=$(docker --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
			print_success "Docker: ${version}"
			found=true
		else
			print_warning "Docker installed but not running"
		fi
	fi

	# Check Apptainer/Singularity
	if command -v apptainer &>/dev/null; then
		local version
		version=$(apptainer --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
		print_success "Apptainer: ${version}"
		found=true
	elif command -v singularity &>/dev/null; then
		local version
		version=$(singularity --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
		print_success "Singularity: ${version}"
		found=true
	fi

	# Check Podman
	if command -v podman &>/dev/null; then
		local version
		version=$(podman --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
		print_success "Podman: ${version}"
		found=true
	fi

	if [[ "$found" == "false" ]]; then
		print_error "No container engine found"
		echo
		echo "    NVD requires Docker, Apptainer, or Podman."
		echo
		echo "    Install Docker:    https://docs.docker.com/get-docker/"
		echo "    Install Apptainer: https://apptainer.org/docs/admin/main/installation.html"
		echo
		return 1
	fi
}

run_preflight_checks() {
	local failed=false

	check_git || failed=true
	check_tar || failed=true
	check_pixi || failed=true
	check_java || failed=true
	check_nextflow || failed=true
	check_container_engine || failed=true

	if [[ "$failed" == "true" ]]; then
		echo
		print_error "Some dependencies are missing. Please install them and try again."
		exit 1
	fi

	echo
	print_success "All dependencies satisfied"
}

# =============================================================================
# Versioned Repo Setup
# =============================================================================

NVD_REPO_URL="https://github.com/dhoconno/nvd.git"
NVD_DIR="${HOME}/.nvd"
LATEST_LINK="${NVD_DIR}/latest"

# Set by setup_repo, used by main
NVD_REPO=""

get_repo_version() {
	local repo_dir="$1"
	local pyproject="${repo_dir}/pyproject.toml"

	[[ -n "$repo_dir" ]] || {
		print_error "INTERNAL: repo_dir is empty"
		return 1
	}

	if [[ ! -f "$pyproject" ]]; then
		print_error "pyproject.toml not found in ${repo_dir}"
		return 1
	fi

	# Use [[:space:]] instead of \s for POSIX compatibility
	# Handle both single and double quotes around version string
	local version
	version=$(grep -E '^[[:space:]]*version[[:space:]]*=' "$pyproject" | head -1 | sed -E "s/.*['\"]([^'\"]+)['\"].*/\1/")

	if [[ -z "$version" ]]; then
		print_error "Could not parse version from pyproject.toml"
		return 1
	fi

	echo "$version"
}

# Helper: create symlink atomically to avoid race conditions
create_version_symlink() {
	local version="$1"
	local temp_link="${LATEST_LINK}.tmp.$$"
	ln -sfn "v${version}" "$temp_link"
	mv -f "$temp_link" "$LATEST_LINK"
}

# Helper: safely remove a directory with sanity checks
safe_rmrf() {
	local dir="$1"
	# Sanity check: never rm -rf empty path, root, or home
	if [[ -z "$dir" ]] || [[ "$dir" == "/" ]] || [[ "$dir" == "$HOME" ]]; then
		print_error "INTERNAL ERROR: refusing to rm -rf '${dir}'"
		return 1
	fi
	rm -rf "$dir"
}

setup_repo() {
	mkdir -p "$NVD_DIR" || {
		print_error "Failed to create directory: ${NVD_DIR}"
		return 1
	}

	if [[ -L "$LATEST_LINK" ]] && [[ -d "$LATEST_LINK" ]]; then
		# Existing install - offer to update
		print_info "NVD is already installed at ${LATEST_LINK}"

		if prompt_yes_no "Pull latest updates?" "y"; then
			if [[ "$DRY_RUN" == "true" ]]; then
				print_info "[DRY RUN] Would pull updates"
				print_info "[DRY RUN] Would run pixi install"
			else
				# Check for uncommitted changes before pulling
				if ! git -C "$LATEST_LINK" diff-index --quiet HEAD -- 2>/dev/null; then
					print_warning "Local modifications detected in ${LATEST_LINK}"
					print_info "Please commit or stash changes before updating"
					print_info "Try: cd ${LATEST_LINK} && git status"
					return 1
				fi

				print_info "Pulling updates..."
				if ! git -C "$LATEST_LINK" pull --ff-only; then
					print_error "Failed to pull updates"
					print_info "You may need to resolve conflicts manually"
					print_info "Try: cd ${LATEST_LINK} && git status"
					return 1
				fi

				# Check if version changed after pull
				local new_version old_dir new_dir
				new_version=$(get_repo_version "$LATEST_LINK") || return 1
				old_dir=$(readlink "$LATEST_LINK")
				new_dir="${NVD_DIR}/v${new_version}"

				# Handle relative vs absolute symlink (readlink returns target as-is)
				if [[ "$old_dir" != /* ]]; then
					old_dir="${NVD_DIR}/${old_dir}"
				fi

				# Validate symlink target exists
				if [[ ! -d "$old_dir" ]]; then
					print_error "Symlink target does not exist: ${old_dir}"
					return 1
				fi

				# Rename directory if version changed
				if [[ "$old_dir" != "$new_dir" ]] && [[ ! -d "$new_dir" ]]; then
					print_info "Version changed, updating directory..."
					mv "$old_dir" "$new_dir" || {
						print_error "Failed to rename directory"
						return 1
					}
					create_version_symlink "$new_version"
					print_success "Updated to v${new_version}"
				else
					print_success "Already at v${new_version}"
				fi

				# Re-run pixi install in case dependencies changed
				print_info "Updating pixi environment..."
				if ! (cd "$LATEST_LINK" && pixi install); then
					print_error "pixi install failed"
					return 1
				fi
			fi
		fi
	else
		# Fresh install - clone to temp, then move to versioned dir
		if [[ "$DRY_RUN" == "true" ]]; then
			print_info "[DRY RUN] Would clone ${NVD_REPO_URL}"
			print_info "[DRY RUN] Would create versioned directory and symlink"
			print_info "[DRY RUN] Would run pixi install"
		else
			local temp_dir="${NVD_DIR}/.clone-temp"
			safe_rmrf "$temp_dir" || return 1

			print_info "Cloning NVD repository..."
			if ! git clone "$NVD_REPO_URL" "$temp_dir"; then
				print_error "Failed to clone repository"
				return 1
			fi

			local version
			version=$(get_repo_version "$temp_dir") || {
				safe_rmrf "$temp_dir"
				return 1
			}

			local versioned_dir="${NVD_DIR}/v${version}"

			if [[ -d "$versioned_dir" ]]; then
				print_info "Version ${version} already exists, using existing"
				safe_rmrf "$temp_dir" || return 1
			else
				mv "$temp_dir" "$versioned_dir" || {
					print_error "Failed to move cloned repo to ${versioned_dir}"
					return 1
				}
			fi

			# Create symlink atomically
			create_version_symlink "$version"

			print_info "Installing pixi environment..."
			if ! (cd "$LATEST_LINK" && pixi install); then
				print_error "pixi install failed"
				return 1
			fi
			print_success "Installed v${version}"
		fi
	fi

	NVD_REPO="$LATEST_LINK"
}

# =============================================================================
# Database Wizard
# =============================================================================

# Database URLs and checksums (update these for new releases)
STAT_DB_URL="https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/stat_db.tar.gz"
BLAST_DB_URL="https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/blast_db.tar.gz"
GOTTCHA2_DB_URL="https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/gottcha2.tar.gz"

STAT_DB_MD5="2641a1754f6986eedbe9f38e188f2b0c"
BLAST_DB_MD5="cebc4c59ea572c09fb93aa4e3594bf3e"
GOTTCHA2_DB_MD5="d33fb5d1b2d22f7a174239f1dfc142cb"

# Approximate size per database in GB (for disk space checks)
DB_SIZE_GB=500

check_disk_space() {
	local path="$1"
	local required_gb="$2"

	[[ -n "$path" ]] || {
		print_error "INTERNAL: path is empty"
		return 1
	}
	[[ -n "$required_gb" ]] || {
		print_error "INTERNAL: required_gb is empty"
		return 1
	}

	# Find existing parent directory for disk space check
	local check_path="$path"
	while [[ ! -d "$check_path" ]] && [[ "$check_path" != "/" ]]; do
		check_path=$(dirname "$check_path")
	done

	# Use POSIX-compliant df -P -k for maximum portability
	# -P forces POSIX output format (consistent column layout)
	# -k uses 1024-byte blocks (universally supported)
	local available_kb available_gb
	available_kb=$(df -P -k "$check_path" 2>/dev/null | awk 'NR==2 {print $4}')

	if [[ -z "$available_kb" ]] || [[ ! "$available_kb" =~ ^[0-9]+$ ]]; then
		print_warning "Could not determine available disk space"
		return 0
	fi

	available_gb=$((available_kb / 1024 / 1024))

	if [[ "$available_gb" -lt "$required_gb" ]]; then
		print_error "Insufficient disk space at ${path}"
		echo "      Required:  ${required_gb}GB"
		echo "      Available: ${available_gb}GB"
		return 1
	fi

	print_success "Disk space OK: ${available_gb}GB available (need ${required_gb}GB)"
	return 0
}

download_database() {
	local url="$1"
	local dest="$2"
	local name="$3"
	local expected_md5="$4"

	# Assertions
	[[ -n "$url" ]] || {
		print_error "INTERNAL: url is empty"
		return 1
	}
	[[ -n "$dest" ]] || {
		print_error "INTERNAL: dest is empty"
		return 1
	}
	[[ -n "$name" ]] || {
		print_error "INTERNAL: name is empty"
		return 1
	}
	[[ -n "$expected_md5" ]] || {
		print_error "INTERNAL: expected_md5 is empty"
		return 1
	}

	local archive="${dest}.tar.gz"

	print_info "Downloading ${name} database..."
	echo "      Source: ${url}"
	echo "      Target: ${archive}"
	echo

	if [[ "$DRY_RUN" == "true" ]]; then
		print_info "[DRY RUN] Would download ${name} database"
		print_info "[DRY RUN] Would verify MD5 checksum"
		print_info "[DRY RUN] Would extract to ${dest}"
		return 0
	fi

	mkdir -p "$(dirname "$archive")" || {
		print_error "Failed to create directory for archive"
		return 1
	}

	# Download with resume support
	if command -v wget &>/dev/null; then
		wget --continue --progress=bar:force -O "$archive" "$url" || {
			print_error "Download failed"
			return 1
		}
	elif command -v curl &>/dev/null; then
		curl -fSL -C - -o "$archive" "$url" || {
			print_error "Download failed"
			return 1
		}
	else
		print_error "Neither wget nor curl found"
		return 1
	fi

	# Verify MD5 checksum
	print_info "Verifying checksum..."
	local actual_md5
	if command -v md5sum &>/dev/null; then
		actual_md5=$(md5sum "$archive" | awk '{print $1}')
	elif command -v md5 &>/dev/null; then
		actual_md5=$(md5 -q "$archive")
	else
		print_warning "No MD5 tool found, skipping verification"
		actual_md5="$expected_md5"
	fi

	if [[ "$actual_md5" != "$expected_md5" ]]; then
		print_error "Checksum mismatch!"
		echo "      Expected: ${expected_md5}"
		echo "      Got:      ${actual_md5}"
		return 1
	fi
	print_success "Checksum verified"

	# Extract
	print_info "Extracting to ${dest}..."
	mkdir -p "$dest" || {
		print_error "Failed to create extraction directory"
		return 1
	}
	tar -xzf "$archive" -C "$dest" || {
		print_error "Extraction failed"
		return 1
	}
	print_success "Extracted successfully"

	# Offer to remove archive to save space
	if prompt_yes_no "Remove archive to save disk space?" "y"; then
		rm -f "$archive"
		print_success "Archive removed"
	fi

	return 0
}

database_wizard() {
	print_header "Database Setup"

	echo "  NVD uses large reference databases for analysis."
	echo "  Each database is approximately ${DB_SIZE_GB}GB."
	echo

	local choice
	choice=$(prompt_choice "Which databases do you need?" \
		"STAT+BLAST (human virus detection)" \
		"GOTTCHA2 (general taxonomic classification)" \
		"Both (recommended for full functionality)" \
		"Skip database download")

	if [[ "$choice" == "4" ]]; then
		print_info "Skipping database download"
		echo
		echo "  You can download databases later by running this installer again"
		echo "  or by downloading manually from:"
		echo "    https://github.com/dhoconno/nvd#databases"
		echo
		return 0
	fi

	# Get database path
	local default_db_path="${HOME}/.nvd/databases"
	echo
	local db_path
	db_path=$(prompt_path "Where should databases be stored?" "$default_db_path")
	echo

	local stat_path="" blast_path="" gottcha_path=""
	local required_space=0

	# Determine which databases to download
	if [[ "$choice" == "1" ]] || [[ "$choice" == "3" ]]; then
		stat_path="${db_path}/stat_db"
		blast_path="${db_path}/blast_db"
		required_space=$((required_space + DB_SIZE_GB * 2))
	fi

	if [[ "$choice" == "2" ]] || [[ "$choice" == "3" ]]; then
		gottcha_path="${db_path}/gottcha2_db"
		required_space=$((required_space + DB_SIZE_GB))
	fi

	# Check disk space
	if ! check_disk_space "$db_path" "$required_space"; then
		echo
		if ! prompt_yes_no "Continue anyway?" "n"; then
			print_info "Database download cancelled"
			return 1
		fi
	fi

	echo

	# Download databases
	if [[ -n "$stat_path" ]]; then
		download_database "$STAT_DB_URL" "$stat_path" "STAT" "$STAT_DB_MD5" || {
			print_error "STAT database download failed"
			return 1
		}
		echo
	fi

	if [[ -n "$blast_path" ]]; then
		download_database "$BLAST_DB_URL" "$blast_path" "BLAST" "$BLAST_DB_MD5" || {
			print_error "BLAST database download failed"
			return 1
		}
		echo
	fi

	if [[ -n "$gottcha_path" ]]; then
		download_database "$GOTTCHA2_DB_URL" "$gottcha_path" "GOTTCHA2" "$GOTTCHA2_DB_MD5" || {
			print_error "GOTTCHA2 database download failed"
			return 1
		}
		echo
	fi

	# Print summary of paths for user to add to config
	print_header "Database Paths"
	echo
	echo "  Add these paths to your Nextflow config (~/.nvd/user.config):"
	echo
	if [[ -n "$stat_path" ]]; then
		echo "    params {"
		echo "        stat_index          = \"${stat_path}/tree_index.dbs\""
		echo "        stat_dbss           = \"${stat_path}/tree_filter.dbss\""
		echo "        stat_annotation     = \"${stat_path}/tree_filter.dbss.annotation\""
		echo "        human_virus_taxlist = \"${stat_path}/human_viruses_taxlist.txt\""
		echo "    }"
		echo
	fi
	if [[ -n "$blast_path" ]]; then
		echo "    params {"
		echo "        blast_db        = \"${blast_path}\""
		echo "        blast_db_prefix = \"core_nt\""
		echo "    }"
		echo
	fi
	if [[ -n "$gottcha_path" ]]; then
		echo "    params {"
		echo "        gottcha2_db = \"${gottcha_path}/gottcha_db.species.fna\""
		echo "    }"
		echo
	fi
}

# =============================================================================
# Argument Parsing and Main
# =============================================================================

DRY_RUN=false
NON_INTERACTIVE=false

show_help() {
	cat <<'EOF'
NVD Installer

Usage: install.sh [OPTIONS]

Options:
    (none)              Interactive installation
    --dry-run           Simulate without making changes
    --non-interactive   Check dependencies only (for CI/CD)
    -h, --help          Show this help
    -v, --version       Show version

Examples:
    ./install.sh                    # Full interactive install
    ./install.sh --dry-run          # See what would happen
    ./install.sh --non-interactive  # CI/CD dependency check

For more information: https://github.com/dhoconno/nvd
EOF
}

main() {
	# Parse arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--dry-run)
			DRY_RUN=true
			shift
			;;
		--non-interactive)
			NON_INTERACTIVE=true
			shift
			;;
		-h | --help)
			show_help
			exit 0
			;;
		-v | --version)
			echo "NVD Installer v${VERSION}"
			exit 0
			;;
		*)
			echo "Unknown option: $1" >&2
			show_help
			exit 1
			;;
		esac
	done

	if [[ "$DRY_RUN" == "true" ]]; then
		echo
		echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
		echo -e "${YELLOW}  DRY RUN MODE - No changes will be made${RESET}"
		echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
	fi

	print_header "NVD Installer v${VERSION}"

	# Phase 1: Preflight checks
	print_info "Checking dependencies..."
	echo
	run_preflight_checks
	echo

	# Exit early for non-interactive mode (CI/CD dependency check only)
	if [[ "$NON_INTERACTIVE" == "true" ]]; then
		echo
		echo "Non-interactive mode: dependency check complete."
		exit 0
	fi

	# Phase 2: Repo setup
	print_info "Setting up NVD repository..."
	setup_repo
	echo

	# Phase 3: Run nvd setup
	print_info "Running nvd setup..."
	local nvd_bin="${NVD_REPO}/.pixi/envs/default/bin/nvd"
	if [[ "$DRY_RUN" == "true" ]]; then
		print_info "[DRY RUN] Would run: ${nvd_bin} setup"
	elif "${nvd_bin}" setup --help &>/dev/null; then
		"${nvd_bin}" setup
	else
		print_warning "nvd setup command not available in this version"
		print_info "You may need to run 'nvd setup' manually after updating"
	fi
	echo

	# Phase 4: Database wizard
	if prompt_yes_no "Download reference databases?" "n"; then
		database_wizard
	fi

	# Done
	print_header "Installation Complete"
	echo
	echo "  Next steps:"
	echo "    1. Start a new shell (or run: source ~/.bashrc)"
	echo "    2. Run: nvd --help"
	echo "    3. See documentation: https://github.com/dhoconno/nvd"
	echo
}

main "$@"
