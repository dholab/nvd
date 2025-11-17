#!/usr/bin/env bash
#
# NVD2 Interactive Installer
#
# This script helps set up NVD2 (a Nextflow pipeline for metagenomic analysis)
# by checking dependencies, optionally downloading databases, and generating
# configuration files.
#
# Usage:
#   curl -fsSL https://raw.githubusercontent.com/dhoconno/nvd/main/install.sh | bash
#   OR
#   ./install.sh [OPTIONS]
#
# License: GPLv3
# Project: https://github.com/dhoconno/nvd
#

# ShellCheck notes:
# - SC2155 disabled globally: we don't check return codes on local var assignments
# - SC2034 (unused variables): Some are intentional for readability or future use
# - SC2086 (unquoted variables): Remaining cases are safe in context
# shellcheck disable=SC2155

set -euo pipefail

# Script version
readonly VERSION="0.1.0"

# Error handling trap
handle_error() {
	local exit_code=$?
	local line_number=$1

	echo >&2
	echo -e "\033[0;31m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\033[0m" >&2
	echo -e "\033[0;31m  ERROR: Installation script failed\033[0m" >&2
	echo -e "\033[0;31m━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\033[0m" >&2
	echo >&2
	echo "An unexpected error occurred at line $line_number (exit code: $exit_code)" >&2
	echo >&2
	echo "This might be due to:" >&2
	echo "  • Network connectivity issues" >&2
	echo "  • Insufficient permissions" >&2
	echo "  • Missing system dependencies" >&2
	echo "  • Disk space issues" >&2
	echo >&2

	# Save error log
	local log_dir="${HOME}/.nvd2"
	local log_file="${log_dir}/install-error.log"
	mkdir -p "$log_dir" 2>/dev/null || true

	{
		echo "NVD2 Installation Error Log"
		echo "Timestamp: $(date)"
		echo "Script Version: $VERSION"
		echo "Exit Code: $exit_code"
		echo "Line Number: $line_number"
		echo "OS: $OSTYPE"
		echo ""
		echo "Error occurred during NVD2 installation"
	} >"$log_file" 2>/dev/null || true

	if [[ -f "$log_file" ]]; then
		echo "Error details saved to: $log_file" >&2
	fi

	echo >&2
	echo "To get help:" >&2
	echo "  • Check the documentation: https://github.com/dhoconno/nvd" >&2
	echo "  • Report issues: https://github.com/dhoconno/nvd/issues" >&2
	echo >&2

	exit $exit_code
}

# Set up error trap (skip in non-interactive mode to avoid breaking CI)
if [[ "${NON_INTERACTIVE:-false}" != "true" ]]; then
	trap 'handle_error $LINENO' ERR
fi

# Global flags (modified by command-line arguments)
DRY_RUN=false
NON_INTERACTIVE=false
QUIET=false

# Global variables for dependency versions (populated by check functions)
JAVA_VERSION=""
NEXTFLOW_VERSION=""
DOCKER_VERSION=""
PODMAN_VERSION=""
APPTAINER_VERSION=""
APPTAINER_CMD=""
PIXI_VERSION=""

# ============================================================================
# ANSI Colors and Symbols
# ============================================================================

# Colors
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly BLUE='\033[0;34m'
readonly YELLOW='\033[0;33m'
readonly CYAN='\033[0;36m'
readonly DIM='\033[2m'
readonly BOLD='\033[1m'
readonly RESET='\033[0m'

# Semantic colors
readonly COLOR_HEADER="${BOLD}${BLUE}"
readonly COLOR_SUCCESS="${GREEN}"
readonly COLOR_ERROR="${RED}"
readonly COLOR_WARNING="${YELLOW}"
readonly COLOR_INFO="${CYAN}"
readonly COLOR_DRY_RUN="${YELLOW}"

# Check for Unicode support
if [[ "${LANG}" =~ UTF-8 ]]; then
	readonly CHECKMARK="✓"
	readonly CROSS="✗"
	readonly ARROW="→"
	readonly INFO="ℹ️ "
	readonly WARNING="⚠️ "
	readonly PACKAGE="📦"
	readonly DOWNLOAD="⬇️ "
	readonly WRENCH="🔧"
	readonly DNA="🧬"

	# Box drawing characters for headers
	readonly BOX_TL="┏"
	readonly BOX_TR="┓"
	readonly BOX_BL="┗"
	readonly BOX_BR="┛"
	readonly BOX_H="━"
	readonly BOX_V="┃"
else
	# ASCII fallbacks
	readonly CHECKMARK="[OK]"
	readonly CROSS="[X]"
	readonly ARROW=">"
	readonly INFO="[i]"
	readonly WARNING="[!]"
	readonly PACKAGE="[PKG]"
	readonly DOWNLOAD="[DOWN]"
	readonly WRENCH="[TOOL]"
	readonly DNA="[DNA]"

	# ASCII box drawing
	readonly BOX_TL="+"
	readonly BOX_TR="+"
	readonly BOX_BL="+"
	readonly BOX_BR="+"
	readonly BOX_H="="
	readonly BOX_V="|"
fi

# ============================================================================
# Utility Functions
# ============================================================================

print_header() {
	local message="$1"
	local width=70

	# Calculate padding: width - message length - 1 (for leading space)
	local padding=$((width - ${#message} - 1))

	echo
	echo -e "${COLOR_HEADER}${BOX_TL}$(printf "${BOX_H}%.0s" $(seq 1 $width))${BOX_TR}${RESET}"
	echo -e "${COLOR_HEADER}${BOX_V} ${message}$(printf ' %.0s' $(seq 1 $padding))${BOX_V}${RESET}"
	echo -e "${COLOR_HEADER}${BOX_BL}$(printf "${BOX_H}%.0s" $(seq 1 $width))${BOX_BR}${RESET}"
	echo
}

print_success() {
	local message="$1"
	echo -e "${COLOR_SUCCESS}  ${CHECKMARK} ${message}${RESET}"
}

print_error() {
	local message="$1"
	echo -e "${COLOR_ERROR}  ${CROSS} ${message}${RESET}" >&2
}

print_warning() {
	local message="$1"
	echo -e "${COLOR_WARNING}  ${WARNING}${message}${RESET}"
}

print_info() {
	local message="$1"
	echo -e "${COLOR_INFO}  ${INFO}${message}${RESET}"
}

print_dry_run() {
	local message="$1"
	echo -e "${COLOR_DRY_RUN}  ${INFO}[DRY RUN] ${message}${RESET}"
}

expand_path() {
	local path="$1"

	# Expand tilde
	path="${path/#\~/$HOME}"

	# Expand $HOME
	path="${path/\$HOME/$HOME}"

	# Expand environment variables
	path=$(eval echo "$path")

	# Get absolute path
	if [[ -d "$path" ]]; then
		path=$(cd "$path" && pwd)
	else
		# Directory doesn't exist yet, make absolute based on current dir
		if [[ ! "$path" =~ ^/ ]]; then
			path="$(pwd)/$path"
		fi
	fi

	echo "$path"
}

validate_path_writeable() {
	local path="$1"

	# Try to create directory
	if mkdir -p "$path" 2>/dev/null; then
		# Check if we can write to it
		if touch "$path/.nvd2_test" 2>/dev/null; then
			rm "$path/.nvd2_test"
			return 0
		else
			print_error "Cannot write to $path (permission denied)"
			return 1
		fi
	else
		print_error "Cannot create directory $path"
		return 1
	fi
}

# ============================================================================
# Error Handling Functions
# ============================================================================

handle_docker_not_running() {
	print_header "Docker Not Running"

	echo "Docker is installed but the Docker daemon is not currently running."
	echo

	# OS-specific instructions
	if [[ "$OSTYPE" == "darwin"* ]]; then
		echo "To start Docker on macOS:"
		echo "  1. Open Docker Desktop from Applications"
		echo "  2. Wait for the Docker icon to appear in the menu bar"
		echo "  3. Verify the icon shows 'Docker Desktop is running'"
	else
		echo "To start Docker on Linux:"
		echo "  1. Run: sudo systemctl start docker"
		echo "  2. Or: sudo service docker start"
		echo
		echo "To enable Docker on boot:"
		echo "  sudo systemctl enable docker"
	fi

	echo
	echo -e "${DIM}After starting Docker, you can:"
	echo "  • Press ENTER to retry the check"
	echo "  • Press Ctrl+C to exit and re-run this installer later${RESET}"
	echo

	# Retry loop
	while true; do
		prompt_continue "Press ENTER to retry Docker check..."

		echo "Checking Docker status..."
		if docker info &>/dev/null; then
			echo
			print_success "Docker is now running!"
			DOCKER_VERSION=$(docker --version 2>&1 | sed -n 's/.*version \([0-9][0-9.]*\).*/\1/p')
			return 0
		else
			echo
			print_warning "Docker is still not running"
			echo

			if ! prompt_yes_no "Try again?" "y"; then
				return 1
			fi
		fi
	done
}

handle_insufficient_disk_space() {
	local path="$1"
	local required_gb="$2"
	local available_gb="$3"
	local shortage_gb=$((required_gb - available_gb))

	print_header "Insufficient Disk Space"

	echo "The selected location does not have enough free space."
	echo
	echo -e "  ${BOLD}Location:${RESET}    $path"
	echo -e "  ${BOLD}Required:${RESET}    ${required_gb}GB"
	echo -e "  ${BOLD}Available:${RESET}   ${available_gb}GB"
	echo -e "  ${BOLD}Shortage:${RESET}    ${RED}${shortage_gb}GB${RESET}"
	echo

	echo "What would you like to do?"
	echo

	local options=(
		"Choose a different location with more space"
		"Configure fewer databases (reduce space requirement)"
		"Skip database setup (configure manually later)"
		"Cancel installation"
	)

	local choice
	choice=$(prompt_choice "Select an option:" "${options[@]}")

	return $((choice - 1)) # Return 0-3 for the selected option
}

handle_network_failure() {
	local db_name="$1"
	local partial_size="$2"

	print_header "Download Interrupted"

	echo "Network connection was lost while downloading ${db_name} database."
	echo

	if [[ -n "$partial_size" ]]; then
		echo -e "  ${BOLD}Downloaded so far:${RESET} $partial_size"
		echo
	fi

	echo "The partial download has been preserved and can be resumed."
	echo
	echo "What would you like to do?"
	echo

	local options=(
		"Retry download now"
		"Skip this database (continue with installation)"
		"Exit (resume installation later)"
	)

	local choice
	choice=$(prompt_choice "Select an option:" "${options[@]}")

	return $((choice - 1)) # Return 0-2 for the selected option
}

handle_old_java_version() {
	local found_version="$1"
	local required_version="11"

	print_header "Java Version Too Old"

	echo "Your Java version is too old for Nextflow."
	echo
	echo -e "  ${BOLD}Found:${RESET}     Java $found_version"
	echo -e "  ${BOLD}Required:${RESET}  Java ${required_version}+"
	echo
	echo "Nextflow requires Java 11 or newer to run."
	echo

	echo "Installation options:"
	echo

	# OS-specific instructions
	if [[ "$OSTYPE" == "darwin"* ]]; then
		echo "  ${BOLD}macOS:${RESET}"
		echo "    brew install openjdk@17"
		echo
	elif command -v apt &>/dev/null; then
		echo "  ${BOLD}Debian/Ubuntu:${RESET}"
		echo "    sudo apt update"
		echo "    sudo apt install openjdk-17-jdk"
		echo
	elif command -v yum &>/dev/null; then
		echo "  ${BOLD}RHEL/CentOS:${RESET}"
		echo "    sudo yum install java-17-openjdk"
		echo
	fi

	echo "  ${BOLD}SDKMAN (any OS):${RESET}"
	echo "    curl -s https://get.sdkman.io | bash"
	echo "    sdk install java 17.0.9-tem"
	echo

	echo -e "${DIM}After installing Java, re-run this installer.${RESET}"
	echo

	return 1
}

# ============================================================================
# Interactive Prompts
# ============================================================================

prompt_yes_no() {
	local prompt="$1"
	local default="${2:-y}"
	local response

	# In dry-run mode, simulate user input
	if [[ "$DRY_RUN" == "true" ]]; then
		if [[ "$default" == "y" ]]; then
			echo -e "${BOLD}${prompt}${RESET} ${DIM}[Y/n]:${RESET} ${YELLOW}(dry-run: defaulting to yes)${RESET}"
		else
			echo -e "${BOLD}${prompt}${RESET} ${DIM}[y/N]:${RESET} ${YELLOW}(dry-run: defaulting to no)${RESET}"
		fi
		echo
		[[ "$default" == "y" ]] && return 0 || return 1
	fi

	# Show prompt with appropriate default indicator (to stderr to avoid capture)
	if [[ "$default" == "y" ]]; then
		echo -en "${BOLD}${prompt}${RESET} ${DIM}[Y/n]:${RESET} " >&2
	else
		echo -en "${BOLD}${prompt}${RESET} ${DIM}[y/N]:${RESET} " >&2
	fi

	read -r response

	# Use default if empty response
	response=${response:-$default}

	# Normalize to lowercase
	response=$(echo "$response" | tr '[:upper:]' '[:lower:]')

	# Evaluate response
	case "$response" in
	y | yes)
		return 0
		;;
	n | no)
		return 1
		;;
	*)
		# Invalid input - treat as default
		[[ "$default" == "y" ]] && return 0 || return 1
		;;
	esac
}

prompt_continue() {
	local message="${1:-Press ENTER to continue, or Ctrl+C to cancel}"

	if [[ "$DRY_RUN" == "true" ]]; then
		echo -e "${DIM}${message}${RESET} ${YELLOW}(dry-run: auto-continuing)${RESET}"
		echo
		return 0
	fi

	echo -e "${DIM}${message}${RESET}"
	read -r
	echo
}

prompt_input() {
	local prompt="$1"
	local default="$2"
	local response

	# In dry-run mode, simulate user input with default
	if [[ "$DRY_RUN" == "true" ]]; then
		if [[ -n "$default" ]]; then
			echo -e "${BOLD}${prompt}${RESET}" >&2
			echo -e "${DIM}Default: ${CYAN}${default}${RESET}" >&2
			echo -e "${YELLOW}(dry-run: using default)${RESET}" >&2
			echo >&2
			echo "$default"
			return 0
		else
			echo -e "${BOLD}${prompt}${RESET}" >&2
			echo -e "${YELLOW}(dry-run: skipping input)${RESET}" >&2
			echo >&2
			return 1
		fi
	fi

	# Show prompt (to stderr so it doesn't get captured)
	echo -e "${BOLD}${prompt}${RESET}" >&2

	# Show default if provided
	if [[ -n "$default" ]]; then
		echo -e "${DIM}Default: ${CYAN}${default}${RESET}" >&2
		echo -en "${DIM}> ${RESET}" >&2
	else
		echo -en "${DIM}> ${RESET}" >&2
	fi

	read -r response

	# Use default if empty response
	if [[ -z "$response" ]] && [[ -n "$default" ]]; then
		response="$default"
	fi

	echo >&2

	# Return the response via stdout (only this gets captured)
	echo "$response"
}

prompt_path() {
	local prompt="$1"
	local default="$2"
	local path

	# In dry-run mode, use default
	if [[ "$DRY_RUN" == "true" ]]; then
		if [[ -n "$default" ]]; then
			echo -e "${BOLD}${prompt}${RESET}" >&2
			echo -e "${DIM}Default: ${CYAN}${default}${RESET}" >&2
			echo -e "${YELLOW}(dry-run: using default)${RESET}" >&2
			echo >&2
			# Expand the default path
			local expanded="${default/#\~/$HOME}"
			expanded="${expanded//\$HOME/$HOME}"
			expanded=$(eval echo "$expanded")
			if [[ ! "$expanded" =~ ^/ ]]; then
				expanded="$(pwd)/$expanded"
			fi
			echo "$expanded"
			return 0
		else
			echo -e "${BOLD}${prompt}${RESET}" >&2
			echo -e "${YELLOW}(dry-run: skipping path input)${RESET}" >&2
			echo >&2
			return 1
		fi
	fi

	while true; do
		# Get input
		path=$(prompt_input "$prompt" "$default")

		# Expand tilde and environment variables
		path="${path/#\~/$HOME}"
		path="${path//\$HOME/$HOME}"
		path=$(eval echo "$path")

		# If path is empty and no default, ask again
		if [[ -z "$path" ]]; then
			print_error "Path cannot be empty" >&2
			echo >&2
			continue
		fi

		# Convert to absolute path if relative
		if [[ ! "$path" =~ ^/ ]]; then
			path="$(pwd)/$path"
		fi

		# Show the expanded path (to stderr)
		echo -e "${DIM}Expanded path: ${CYAN}${path}${RESET}" >&2
		echo >&2

		# Return the path (to stdout)
		echo "$path"
		return 0
	done
}

prompt_choice() {
	local prompt="$1"
	shift
	local options=("$@")
	local choice

	# Validate we have options
	if [[ ${#options[@]} -eq 0 ]]; then
		echo "ERROR: prompt_choice called with no options" >&2
		return 1
	fi

	# Display prompt and options (same for both dry-run and normal mode)
	echo -e "${BOLD}${prompt}${RESET}" >&2
	echo >&2

	local i=1
	for opt in "${options[@]}"; do
		if [[ "$DRY_RUN" == "true" ]]; then
			echo -e "  ${DIM}[$i]${RESET} ${opt}" >&2
		else
			echo -e "  ${CYAN}[$i]${RESET} ${opt}" >&2
		fi
		((i++))
	done

	echo >&2

	# In dry-run mode, automatically select option 1
	if [[ "$DRY_RUN" == "true" ]]; then
		echo -e "${DIM}Enter choice [1-${#options[@]}]:${RESET} ${YELLOW}(dry-run: auto-selecting option 1)${RESET}" >&2
		echo >&2
		echo "1"
		return 0
	fi

	# Get valid choice (interactive in normal mode)
	while true; do
		echo -en "${DIM}Enter choice [1-${#options[@]}]:${RESET} " >&2
		read -r choice

		# Check if input is a number
		if [[ "$choice" =~ ^[0-9]+$ ]]; then
			# Check if in valid range
			if [[ "$choice" -ge 1 ]] && [[ "$choice" -le "${#options[@]}" ]]; then
				echo >&2
				echo "$choice"
				return 0
			fi
		fi

		# Invalid input
		print_error "Invalid choice. Please enter a number between 1 and ${#options[@]}" >&2
		echo >&2
	done
}

# ============================================================================
# Help Text
# ============================================================================

show_help() {
	echo -e "
${BOLD}${BLUE}NVD2 Installation Script${RESET}

This script helps set up NVD2, a Nextflow pipeline for metagenomic analysis
of human viruses and general taxonomic classification.

${BOLD}USAGE:${RESET}
    ${CYAN}./install.sh${RESET} ${DIM}[OPTIONS]${RESET}

${BOLD}OPTIONS:${RESET}
    ${DIM}(none)${RESET}                      Run interactive installation wizard
    ${CYAN}-d${RESET}, ${CYAN}--dry-run${RESET}, ${CYAN}--test${RESET}       Simulate installation without making changes
    ${CYAN}-y${RESET}, ${CYAN}--yes${RESET}, ${CYAN}--non-interactive${RESET} Check prerequisites only (for CI/CD)
    ${CYAN}--verify${RESET}                    Verify existing NVD2 installation
    ${CYAN}--uninstall${RESET}                 Remove NVD2 configuration and databases
    ${CYAN}--quiet${RESET}                     Minimal output
    ${CYAN}-h${RESET}, ${CYAN}--help${RESET}                  Show this help message
    ${CYAN}-v${RESET}, ${CYAN}--version${RESET}               Show version information

${BOLD}EXAMPLES:${RESET}
    ${DIM}# Interactive installation with full wizard${RESET}
    ${CYAN}./install.sh${RESET}

    ${DIM}# Test the installer without making any changes${RESET}
    ${CYAN}./install.sh --dry-run${RESET}

    ${DIM}# Check prerequisites in CI/CD environment${RESET}
    ${CYAN}./install.sh --non-interactive${RESET}

    ${DIM}# Verify existing installation${RESET}
    ${CYAN}./install.sh --verify${RESET}

    ${DIM}# Remove NVD2 from system${RESET}
    ${CYAN}./install.sh --uninstall${RESET}

${BOLD}MODES:${RESET}
    ${GREEN}Interactive${RESET} ${DIM}(default)${RESET}
        - Guides you through dependency installation
        - Offers to download reference databases
        - Generates configuration files
        - Provides detailed explanations

    ${YELLOW}Dry-run${RESET} ${DIM}(--dry-run)${RESET}
        - Simulates the entire installation process
        - Shows what would happen without making changes
        - Useful for testing and preview

    ${BLUE}Non-interactive${RESET} ${DIM}(--non-interactive)${RESET}
        - Only checks if prerequisites are met
        - No prompts, no installations, no modifications
        - Returns exit code 0 if ready, non-zero otherwise
        - Suitable for automated environments

    ${CHECKMARK} ${GREEN}Verify${RESET} ${DIM}(--verify)${RESET}
        - Checks dependencies, configuration, and databases
        - Validates existing installation
        - Reports any issues found

    ${WARNING} ${RED}Uninstall${RESET} ${DIM}(--uninstall)${RESET}
        - Removes NVD2 configuration files
        - Optionally removes databases
        - Shows commands for removing dependencies

${BOLD}MORE INFORMATION:${RESET}
    Documentation: ${CYAN}https://github.com/dhoconno/nvd${RESET}
    Report issues: ${CYAN}https://github.com/dhoconno/nvd/issues${RESET}
"
}

show_version() {
	echo "NVD2 Installer version ${VERSION}"
}

# ============================================================================
# Disk Space Utilities
# ============================================================================

format_size_gb() {
	local size_gb="$1"

	# Convert to TB if >= 1000GB
	if [[ "$size_gb" -ge 1000 ]]; then
		local size_tb=$(awk "BEGIN {printf \"%.1f\", $size_gb / 1000}")
		echo "${size_tb}TB"
	else
		echo "${size_gb}GB"
	fi
}

get_available_space_gb() {
	local path="$1"
	local available_gb

	# Create parent directory if it doesn't exist (for checking mount point)
	local check_path="$path"
	if [[ ! -d "$check_path" ]]; then
		# Find first existing parent directory
		while [[ ! -d "$check_path" ]] && [[ "$check_path" != "/" ]]; do
			check_path=$(dirname "$check_path")
		done
	fi

	# Get available space based on OS
	if [[ "$OSTYPE" == "darwin"* ]]; then
		# macOS - use df -h and parse human-readable output
		local avail_raw
		avail_raw=$(df -h "$check_path" 2>/dev/null | awk 'NR==2 {print $4}')

		# Parse the human-readable size (could be 83G, 1.5T, etc.)
		if [[ "$avail_raw" =~ ^([0-9.]+)([KMGT])i?$ ]]; then
			local number="${BASH_REMATCH[1]}"
			local unit="${BASH_REMATCH[2]}"

			case "$unit" in
			K) available_gb=$(awk "BEGIN {printf \"%d\", $number / 1024 / 1024}") ;;
			M) available_gb=$(awk "BEGIN {printf \"%d\", $number / 1024}") ;;
			G) available_gb=$(awk "BEGIN {printf \"%d\", $number}") ;;
			T) available_gb=$(awk "BEGIN {printf \"%d\", $number * 1024}") ;;
			esac
		fi
	else
		# Linux - use df with block size
		available_gb=$(df -BG "$check_path" 2>/dev/null | awk 'NR==2 {print $4}' | sed 's/G//')
	fi

	# Return the value (or 0 if failed)
	if [[ -n "$available_gb" ]] && [[ "$available_gb" =~ ^[0-9]+$ ]]; then
		echo "$available_gb"
	else
		echo "0"
	fi
}

check_disk_space() {
	local path="$1"
	local required_gb="$2"
	local show_output="${3:-true}" # Optional: set to false for silent check

	# Calculate safety buffer (20% of required, minimum 50GB)
	local buffer_percent=20
	local buffer_gb=$((required_gb * buffer_percent / 100))
	if [[ "$buffer_gb" -lt 50 ]]; then
		buffer_gb=50
	fi

	local required_with_buffer=$((required_gb + buffer_gb))

	# Get available space
	local available_gb
	available_gb=$(get_available_space_gb "$path")

	if [[ "$available_gb" -eq 0 ]]; then
		if [[ "$show_output" == "true" ]]; then
			print_error "Could not determine available disk space at $path"
			echo -e "  ${DIM}Please verify the path is valid and accessible${RESET}"
			echo
		fi
		return 1
	fi

	# Check if sufficient space
	if [[ "$available_gb" -lt "$required_with_buffer" ]]; then
		if [[ "$show_output" == "true" ]]; then
			print_error "Insufficient disk space"
			echo
			echo -e "  ${BOLD}Location:${RESET}       $path"
			echo -e "  ${BOLD}Required:${RESET}       $(format_size_gb $required_gb)"
			echo -e "  ${BOLD}Safety buffer:${RESET}  $(format_size_gb $buffer_gb) ${DIM}(${buffer_percent}%)${RESET}"
			echo -e "  ${BOLD}Total needed:${RESET}   ${YELLOW}$(format_size_gb $required_with_buffer)${RESET}"
			echo -e "  ${BOLD}Available:${RESET}      ${RED}$(format_size_gb $available_gb)${RESET}"
			echo -e "  ${BOLD}Shortage:${RESET}       ${RED}$(format_size_gb $((required_with_buffer - available_gb)))${RESET}"
			echo
		fi
		return 1
	else
		# Sufficient space
		local margin=$((available_gb - required_with_buffer))

		if [[ "$show_output" == "true" ]]; then
			print_success "Sufficient disk space available"
			echo
			echo -e "  ${BOLD}Location:${RESET}       $path"
			echo -e "  ${BOLD}Required:${RESET}       $(format_size_gb $required_gb)"
			echo -e "  ${BOLD}Safety buffer:${RESET}  $(format_size_gb $buffer_gb) ${DIM}(${buffer_percent}%)${RESET}"
			echo -e "  ${BOLD}Total needed:${RESET}   $(format_size_gb $required_with_buffer)"
			echo -e "  ${BOLD}Available:${RESET}      ${GREEN}$(format_size_gb $available_gb)${RESET}"
			echo -e "  ${BOLD}Extra margin:${RESET}   ${GREEN}$(format_size_gb $margin)${RESET}"
			echo
		fi
		return 0
	fi
}

# ============================================================================
# Database Configuration
# ============================================================================

# Database sizes in GB (for disk space checking)
readonly STAT_DB_SIZE=500
readonly BLAST_DB_SIZE=500
readonly GOTTCHA2_DB_SIZE=500

configure_database_path() {
	local db_name="$1"
	local db_size_gb="$2"
	local default_path="$3"
	local require_disk_check="$4"

	# In dry-run mode, use fake path to avoid complexity
	if [[ "$DRY_RUN" == "true" ]]; then
		echo -e "${BOLD}📦 ${db_name} Database${RESET}" >&2
		echo >&2
		echo -e "  ${DIM}Required space: ~${db_size_gb}GB${RESET}" >&2
		echo >&2
		echo -e "${YELLOW}(dry-run: auto-configuring with default path)${RESET}" >&2
		echo >&2

		local fake_path
		if [[ -n "$default_path" ]]; then
			fake_path="$default_path"
		else
			fake_path="/data/nvd2/databases/${db_name,,}_db"
		fi

		print_success "Path configured: ${CYAN}${fake_path}${RESET}" >&2
		print_dry_run "Would create directory: ${fake_path}" >&2
		echo >&2

		echo "$fake_path"
		return 0
	fi

	echo -e "${BOLD}📦 ${db_name} Database${RESET}" >&2
	echo >&2
	echo -e "  ${DIM}Required space: ~${db_size_gb}GB${RESET}" >&2
	echo >&2

	# Ask if user wants to configure this database
	local configure_it
	if prompt_yes_no "Would you like to configure the ${db_name} database path?" "y"; then
		configure_it="yes"
	else
		configure_it="no"
	fi

	if [[ "$configure_it" == "no" ]]; then
		print_info "${db_name} database skipped - will not be configured" >&2
		echo >&2
		return 1
	fi

	local path
	local disk_check_passed=false

	while true; do
		# Get path from user
		path=$(prompt_path "Where should the ${db_name} database be located?" "$default_path")

		# If disk check required, verify space
		if [[ "$require_disk_check" == "true" ]]; then
			echo -e "${DIM}Checking disk space...${RESET}" >&2
			echo >&2

			if check_disk_space "$path" "$db_size_gb"; then
				disk_check_passed=true
				break
			else
				# Insufficient space - use handler
				local available_gb
				available_gb=$(get_available_space_gb "$path")

				# Calculate required with buffer
				local buffer_gb=$((db_size_gb * 20 / 100))
				if [[ "$buffer_gb" -lt 50 ]]; then
					buffer_gb=50
				fi
				local required_with_buffer=$((db_size_gb + buffer_gb))

				echo >&2
				handle_insufficient_disk_space "$path" "$required_with_buffer" "$available_gb"
				local handler_choice=$?

				case $handler_choice in
				0)
					# Choose different location - continue loop
					continue
					;;
				1)
					# Configure fewer databases - return to let caller handle
					print_info "Returning to database selection..." >&2
					echo >&2
					return 2 # Special return code
					;;
				2)
					# Skip database setup
					print_info "${db_name} database skipped" >&2
					echo >&2
					return 1
					;;
				3)
					# Cancel installation
					print_info "Installation cancelled" >&2
					exit 0
					;;
				esac
			fi
		else
			# No disk check required (databases already exist)
			print_success "Path configured: ${CYAN}${path}${RESET}" >&2
			echo >&2
			break
		fi
	done

	# Create directory if it doesn't exist (in real mode)
	if [[ "$DRY_RUN" != "true" ]]; then
		if [[ ! -d "$path" ]]; then
			if mkdir -p "$path" 2>/dev/null; then
				print_success "Created directory: ${path}" >&2
			else
				print_warning "Could not create directory: ${path}" >&2
				echo -e "  ${DIM}You may need to create it manually with appropriate permissions${RESET}" >&2
			fi
			echo >&2
		fi
	else
		print_dry_run "Would create directory: ${path}" >&2
		echo >&2
	fi

	# Return path via stdout
	echo "$path"
	return 0
}

configure_all_databases() {
	local mode="$1" # "have" or "download"

	print_header "Database Path Configuration" >&2

	if [[ "$mode" == "have" ]]; then
		echo "Configuring paths for existing databases..." >&2
		echo >&2
		echo -e "${DIM}Note: No disk space checks will be performed${RESET}" >&2
		echo -e "${DIM}(assuming databases already exist at specified locations)${RESET}" >&2
		echo >&2
	else
		echo "Configuring paths for database downloads..." >&2
		echo >&2
		echo -e "${YELLOW}${BOLD}Important:${RESET}${YELLOW} Disk space will be checked at each location${RESET}" >&2
		echo -e "${DIM}Each database requires 400-500GB + 20% safety buffer${RESET}" >&2
		echo >&2
	fi

	local require_checks="false"
	if [[ "$mode" == "download" ]]; then
		require_checks="true"
	fi

	# Configure STAT database (no default on first)
	local stat_path=""
	stat_path=$(configure_database_path "STAT" "$STAT_DB_SIZE" "" "$require_checks") || true

	# Configure BLAST database (default to same parent as STAT if STAT was configured)
	local blast_path=""
	local blast_default=""
	if [[ -n "$stat_path" ]]; then
		blast_default="$(dirname "$stat_path")/blast_db"
	fi
	blast_path=$(configure_database_path "BLAST" "$BLAST_DB_SIZE" "$blast_default" "$require_checks") || true

	# Configure GOTTCHA2 database (default to same parent as previous)
	local gottcha_path=""
	local gottcha_default=""
	if [[ -n "$stat_path" ]]; then
		gottcha_default="$(dirname "$stat_path")/gottcha2_db"
	elif [[ -n "$blast_path" ]]; then
		gottcha_default="$(dirname "$blast_path")/gottcha2_db"
	fi
	gottcha_path=$(configure_database_path "GOTTCHA2" "$GOTTCHA2_DB_SIZE" "$gottcha_default" "$require_checks") || true

	# Return paths as space-separated string
	echo "$stat_path|$blast_path|$gottcha_path"
}

# ============================================================================
# Database Download & Extraction
# ============================================================================

# Database URLs and checksums
readonly STAT_DB_URL="https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/stat_db.tar.gz"
readonly BLAST_DB_URL="https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/blast_db.tar.gz"
readonly GOTTCHA2_DB_URL="https://dholk.primate.wisc.edu/_webdav/dho/projects/lungfish/InfinitePath/public/%40files/release-v2.0.0/gottcha2.tar.gz"

readonly STAT_DB_MD5="2641a1754f6986eedbe9f38e188f2b0c"
readonly BLAST_DB_MD5="cebc4c59ea572c09fb93aa4e3594bf3e"
readonly GOTTCHA2_DB_MD5="d33fb5d1b2d22f7a174239f1dfc142cb"

download_database() {
	local url="$1"
	local dest_path="$2"
	local db_name="$3"
	local expected_md5="$4"

	local dest_dir="$(dirname "$dest_path")"
	local archive_file="${dest_path}.tar.gz"

	print_header "Downloading ${db_name} Database" >&2

	echo -e "  ${BOLD}Source:${RESET} ${url}" >&2
	echo -e "  ${BOLD}Target:${RESET} ${archive_file}" >&2
	echo >&2

	# Dry-run mode: simulate download
	if [[ "$DRY_RUN" == "true" ]]; then
		print_dry_run "Would download ${db_name} database" >&2
		echo >&2
		print_dry_run "Simulating download progress..." >&2
		echo >&2

		# Simulate realistic download progress (design doc style)
		local total_steps=30
		local bar_width=30
		for i in $(seq 0 $total_steps); do
			local percent=$((i * 100 / total_steps))
			local filled=$((bar_width * i / total_steps))
			local empty=$((bar_width - filled))

			# Calculate simulated download stats
			local downloaded=$(awk "BEGIN {printf \"%.1f\", 500 * $i / $total_steps}")
			local speed=$(awk "BEGIN {printf \"%.1f\", 15.2}")
			local remaining=$((total_steps - i))

			# Print progress bar: ████████████░░░░░  67%
			printf "\r  " >&2
			printf "%${filled}s" | tr ' ' '█' >&2
			printf "%${empty}s" | tr ' ' '░' >&2
			printf "  %3d%%" "$percent" >&2

			sleep 0.05
		done
		echo >&2
		echo >&2

		print_success "Download simulated successfully" >&2
		echo >&2
		return 0
	fi

	# Check if archive already exists
	if [[ -f "$archive_file" ]]; then
		print_info "Archive already exists: ${archive_file}" >&2
		echo -e "  ${DIM}Verifying existing file...${RESET}" >&2
		echo >&2

		# Verify MD5 if file exists
		if command -v md5sum &>/dev/null; then
			local actual_md5=$(md5sum "$archive_file" | awk '{print $1}')
			if [[ "$actual_md5" == "$expected_md5" ]]; then
				print_success "Existing file verified (MD5 matches)" >&2
				echo >&2
				return 0
			else
				print_warning "Existing file MD5 mismatch - will re-download" >&2
				echo >&2
				rm -f "$archive_file"
			fi
		elif command -v md5 &>/dev/null; then
			# macOS md5 command
			local actual_md5=$(md5 -q "$archive_file")
			if [[ "$actual_md5" == "$expected_md5" ]]; then
				print_success "Existing file verified (MD5 matches)" >&2
				echo >&2
				return 0
			else
				print_warning "Existing file MD5 mismatch - will re-download" >&2
				echo >&2
				rm -f "$archive_file"
			fi
		else
			print_warning "Cannot verify MD5 (no md5sum/md5 command)" >&2
			echo -e "  ${DIM}Will attempt to use existing file${RESET}" >&2
			echo >&2
			return 0
		fi
	fi

	# Ensure destination directory exists
	mkdir -p "$dest_dir"

	echo -e "${YELLOW}Large download in progress (hundreds of GB)${RESET}" >&2
	echo -e "${DIM}Download can be interrupted (Ctrl+C) and will resume from the same point${RESET}" >&2
	echo >&2

	# Try wget first (better resume support), then curl
	if command -v wget &>/dev/null; then
		print_info "Using wget for download (with resume support)" >&2
		echo >&2

		wget \
			--continue \
			--progress=bar:force:noscroll \
			--tries=5 \
			--timeout=60 \
			--read-timeout=60 \
			--retry-connrefused \
			--waitretry=3 \
			-O "$archive_file" \
			"$url"

		local wget_result=$?

		if [[ $wget_result -eq 0 ]]; then
			echo >&2
			print_success "Download complete" >&2
		else
			echo >&2

			# Get partial download size if it exists
			local partial_size=""
			if [[ -f "$archive_file" ]]; then
				partial_size=$(du -h "$archive_file" | cut -f1)
			fi

			# Use error handler
			handle_network_failure "$db_name" "$partial_size"
			local handler_choice=$?

			case $handler_choice in
			0)
				# Retry - recursive call
				echo >&2
				print_info "Retrying download..." >&2
				echo >&2
				download_database "$url" "$dest_path" "$db_name" "$expected_md5"
				return $?
				;;
			1)
				# Skip this database
				print_info "Skipping ${db_name} database" >&2
				return 1
				;;
			2)
				# Exit
				print_info "Installation will resume on next run" >&2
				exit 0
				;;
			esac
		fi

	elif command -v curl &>/dev/null; then
		print_info "Using curl for download (with resume support)" >&2
		echo >&2

		curl \
			-fSL \
			-C - \
			--connect-timeout 60 \
			--max-time 0 \
			--retry 5 \
			--retry-delay 3 \
			-o "$archive_file" \
			"$url"

		local curl_result=$?

		if [[ $curl_result -eq 0 ]]; then
			echo >&2
			print_success "Download complete" >&2
		else
			echo >&2

			# Get partial download size if it exists
			local partial_size=""
			if [[ -f "$archive_file" ]]; then
				partial_size=$(du -h "$archive_file" | cut -f1)
			fi

			# Use error handler
			handle_network_failure "$db_name" "$partial_size"
			local handler_choice=$?

			case $handler_choice in
			0)
				# Retry - recursive call
				echo >&2
				print_info "Retrying download..." >&2
				echo >&2
				download_database "$url" "$dest_path" "$db_name" "$expected_md5"
				return $?
				;;
			1)
				# Skip this database
				print_info "Skipping ${db_name} database" >&2
				return 1
				;;
			2)
				# Exit
				print_info "Installation will resume on next run" >&2
				exit 0
				;;
			esac
		fi
	else
		print_error "Neither wget nor curl is available" >&2
		echo -e "  ${DIM}Please install wget or curl and try again${RESET}" >&2
		return 1
	fi

	# Verify MD5 checksum
	echo >&2
	echo -e "${DIM}Verifying download integrity...${RESET}" >&2

	if command -v md5sum &>/dev/null; then
		local actual_md5=$(md5sum "$archive_file" | awk '{print $1}')
		if [[ "$actual_md5" == "$expected_md5" ]]; then
			print_success "MD5 checksum verified" >&2
		else
			print_error "MD5 checksum mismatch!" >&2
			echo -e "  ${BOLD}Expected:${RESET} $expected_md5" >&2
			echo -e "  ${BOLD}Actual:${RESET}   $actual_md5" >&2
			echo >&2
			echo -e "  ${YELLOW}File may be corrupted. Try downloading again.${RESET}" >&2
			return 1
		fi
	elif command -v md5 &>/dev/null; then
		local actual_md5=$(md5 -q "$archive_file")
		if [[ "$actual_md5" == "$expected_md5" ]]; then
			print_success "MD5 checksum verified" >&2
		else
			print_error "MD5 checksum mismatch!" >&2
			echo -e "  ${BOLD}Expected:${RESET} $expected_md5" >&2
			echo -e "  ${BOLD}Actual:${RESET}   $actual_md5" >&2
			echo >&2
			echo -e "  ${YELLOW}File may be corrupted. Try downloading again.${RESET}" >&2
			return 1
		fi
	else
		print_warning "Cannot verify MD5 checksum (no md5sum/md5 command)" >&2
		echo -e "  ${DIM}Proceeding anyway, but file integrity not guaranteed${RESET}" >&2
	fi

	echo >&2
	return 0
}

extract_database() {
	local archive_file="$1"
	local dest_path="$2"
	local db_name="$3"

	print_header "Extracting ${db_name} Database" >&2

	echo -e "  ${BOLD}Archive:${RESET} ${archive_file}" >&2
	echo -e "  ${BOLD}Target:${RESET}  ${dest_path}" >&2
	echo >&2

	# Dry-run mode: simulate extraction
	if [[ "$DRY_RUN" == "true" ]]; then
		print_dry_run "Would extract ${db_name} database" >&2
		echo >&2
		print_dry_run "Simulating extraction..." >&2
		echo >&2

		# Simulate extraction progress (design doc style)
		local total_steps=30
		local bar_width=30
		for i in $(seq 0 $total_steps); do
			local percent=$((i * 100 / total_steps))
			local filled=$((bar_width * i / total_steps))
			local empty=$((bar_width - filled))

			# Print progress bar: ████████████████████████████████  100%
			printf "\r  " >&2
			printf "%${filled}s" | tr ' ' '█' >&2
			printf "%${empty}s" | tr ' ' '░' >&2
			printf "  %3d%%" "$percent" >&2

			sleep 0.05
		done
		echo >&2
		echo >&2

		print_success "Extraction simulated successfully" >&2
		echo >&2
		return 0
	fi

	# Check if archive exists
	if [[ ! -f "$archive_file" ]]; then
		print_error "Archive file not found: $archive_file" >&2
		return 1
	fi

	# Create destination directory
	mkdir -p "$dest_path"

	# Extract with tar
	echo -e "${DIM}Extracting (this may take several minutes)...${RESET}" >&2
	echo >&2

	if tar -xzf "$archive_file" -C "$dest_path" 2>&1 | sed 's/^/  /' >&2; then
		echo >&2
		print_success "Extraction complete" >&2
		echo >&2

		# Ask if user wants to remove archive
		if prompt_yes_no "Remove archive file to save space? (archive: $(du -h "$archive_file" | cut -f1))" "n"; then
			rm -f "$archive_file"
			print_success "Archive removed" >&2
		else
			print_info "Archive kept: ${archive_file}" >&2
		fi
		echo >&2

		return 0
	else
		echo >&2
		print_error "Extraction failed" >&2
		echo >&2
		return 1
	fi
}

download_and_extract_all() {
	local stat_path="$1"
	local blast_path="$2"
	local gottcha_path="$3"

	print_header "Database Download & Extraction" >&2

	echo "Preparing to download reference databases..." >&2
	echo >&2
	echo -e "${YELLOW}${BOLD}Note:${RESET}${YELLOW} Each database is 400-500GB and will take significant time${RESET}" >&2
	echo -e "${DIM}Stable internet connection required. Downloads are resumable if interrupted.${RESET}" >&2
	echo >&2

	local download_failed=false

	# Download and extract STAT database
	if [[ -n "$stat_path" ]]; then
		if download_database "$STAT_DB_URL" "$stat_path" "STAT" "$STAT_DB_MD5"; then
			extract_database "${stat_path}.tar.gz" "$stat_path" "STAT" || download_failed=true
		else
			download_failed=true
		fi
	fi

	# Download and extract BLAST database
	if [[ -n "$blast_path" ]]; then
		if download_database "$BLAST_DB_URL" "$blast_path" "BLAST" "$BLAST_DB_MD5"; then
			extract_database "${blast_path}.tar.gz" "$blast_path" "BLAST" || download_failed=true
		else
			download_failed=true
		fi
	fi

	# Download and extract GOTTCHA2 database
	if [[ -n "$gottcha_path" ]]; then
		if download_database "$GOTTCHA2_DB_URL" "$gottcha_path" "GOTTCHA2" "$GOTTCHA2_DB_MD5"; then
			extract_database "${gottcha_path}.tar.gz" "$gottcha_path" "GOTTCHA2" || download_failed=true
		else
			download_failed=true
		fi
	fi

	echo >&2
	if [[ "$download_failed" == "true" ]]; then
		print_warning "Some downloads or extractions failed" >&2
		echo >&2
		echo -e "  ${DIM}You can re-run this installer to retry failed downloads${RESET}" >&2
		echo -e "  ${DIM}Downloads will resume from where they left off${RESET}" >&2
		echo >&2
		return 1
	else
		print_success "All databases downloaded and extracted successfully!" >&2
		echo >&2
		return 0
	fi
}

# ============================================================================
# Pipeline Source Code Setup
# ============================================================================

setup_pipeline_source() {
	print_header "Pipeline Source Code Setup" >&2

	echo "NVD2 pipeline code can be obtained in different ways:" >&2
	echo >&2
	echo -e "${DIM}The pipeline source (~few MB) is separate from databases (~1TB).${RESET}" >&2
	echo >&2

	local source_options=(
		"${GREEN}Auto-download on first run${RESET} ${DIM}(recommended)${RESET}
        ${DIM}Nextflow fetches from GitHub automatically: dhoconno/nvd${RESET}"

		"${CYAN}Pre-cache with nextflow pull${RESET} ${DIM}(offline-ready)${RESET}
        ${DIM}Downloads now to: ~/.nextflow/assets/dhoconno/nvd${RESET}"

		"${YELLOW}Clone to current directory${RESET} ${DIM}(for development)${RESET}
        ${DIM}git clone to: $(pwd)/nvd${RESET}"

		"${YELLOW}Clone to custom directory${RESET} ${DIM}(specify location)${RESET}
        ${DIM}You choose where to clone the repository${RESET}"
	)

	local source_choice
	source_choice=$(prompt_choice "How would you like to set up the pipeline source code?" "${source_options[@]}")

	local clone_path=""

	case "$source_choice" in
	1)
		# Auto-download - do nothing
		print_success "Pipeline will auto-download on first run" >&2
		echo >&2
		echo -e "${DIM}Nextflow will automatically fetch the pipeline when you run:${RESET}" >&2
		echo -e "${CYAN}nextflow run dhoconno/nvd ...${RESET}" >&2
		echo >&2
		echo "$clone_path" # Return empty string
		;;
	2)
		# Pre-cache with nextflow pull
		echo
		print_info "Pre-caching pipeline with nextflow pull..."
		echo

		if [[ "$DRY_RUN" == "true" ]]; then
			print_dry_run "Would run: nextflow pull dhoconno/nvd"
			echo
			print_success "Pipeline would be cached to ~/.nextflow/assets/dhoconno/nvd"
		else
			if nextflow pull dhoconno/nvd; then
				echo >&2
				print_success "Pipeline cached successfully!" >&2
				echo >&2
				echo -e "${GREEN}✓${RESET} Location: ${CYAN}~/.nextflow/assets/dhoconno/nvd${RESET}" >&2
			else
				echo >&2
				print_warning "Failed to cache pipeline" >&2
				echo >&2
				echo -e "${DIM}Not critical - Nextflow will download on first run${RESET}" >&2
			fi
		fi
		echo >&2
		echo "$clone_path" # Return empty string (use dhoconno/nvd)
		;;
	3)
		# Clone to current directory
		clone_path="$(pwd)/nvd"
		echo >&2
		echo -e "Clone target: ${CYAN}${clone_path}${RESET}" >&2
		echo >&2

		if [[ "$DRY_RUN" == "true" ]]; then
			print_dry_run "Would run: git clone https://github.com/dhoconno/nvd.git ${clone_path}"
			echo
			print_success "Would clone to ${clone_path}"
		else
			if [[ -d "$clone_path" ]]; then
				print_warning "Directory already exists: ${clone_path}"
				echo
				if prompt_yes_no "Pull latest changes instead?" "Y"; then
					if git -C "$clone_path" pull; then
						print_success "Updated to latest version"
					else
						print_error "Failed to update repository"
						clone_path=""
					fi
				fi
			else
				if git clone https://github.com/dhoconno/nvd.git "$clone_path"; then
					echo >&2
					print_success "Cloned successfully!" >&2
					echo >&2
					echo -e "${GREEN}✓${RESET} Location: ${CYAN}${clone_path}${RESET}" >&2
				else
					echo >&2
					print_error "Failed to clone repository" >&2
					clone_path=""
				fi
			fi
		fi
		echo >&2
		echo "$clone_path"
		;;
	4)
		# Clone to custom directory
		echo >&2
		echo "Where should the pipeline be cloned?" >&2
		echo -e "${DIM}(default: current directory)${RESET}" >&2
		echo -n "> " >&2
		read -r custom_path

		if [[ -z "$custom_path" ]]; then
			custom_path="."
		fi

		# Expand and validate path
		custom_path=$(expand_path "$custom_path")
		clone_path="${custom_path}/nvd"

		echo >&2
		echo -e "Clone target: ${CYAN}${clone_path}${RESET}" >&2
		echo >&2

		if [[ "$DRY_RUN" == "true" ]]; then
			print_dry_run "Would run: git clone https://github.com/dhoconno/nvd.git ${clone_path}"
			echo
			print_success "Would clone to ${clone_path}"
		else
			# Create parent directory if needed
			mkdir -p "$custom_path" 2>/dev/null

			if [[ -d "$clone_path" ]]; then
				print_warning "Directory already exists: ${clone_path}"
				echo
				if prompt_yes_no "Pull latest changes instead?" "Y"; then
					if git -C "$clone_path" pull; then
						print_success "Updated to latest version"
					else
						print_error "Failed to update repository"
						clone_path=""
					fi
				fi
			else
				if git clone https://github.com/dhoconno/nvd.git "$clone_path"; then
					echo >&2
					print_success "Cloned successfully!" >&2
					echo >&2
					echo -e "${GREEN}✓${RESET} Location: ${CYAN}${clone_path}${RESET}" >&2
				else
					echo >&2
					print_error "Failed to clone repository" >&2
					clone_path=""
				fi
			fi
		fi
		echo >&2
		echo "$clone_path"
		;;
	esac
}

# ============================================================================
# Configuration File Generation
# ============================================================================

generate_user_config() {
	local stat_path="$1"
	local blast_path="$2"
	local gottcha_path="$3"
	local exec_env="$4" # docker, podman, apptainer, or pixi

	print_header "Generating Configuration File"

	# Define config directory and file
	local config_dir="${HOME}/.nvd2/config"
	local config_file="${config_dir}/user.config"

	# Create config directory if it doesn't exist
	mkdir -p "$config_dir"

	echo "Creating user configuration at:"
	echo -e "  ${CYAN}${config_file}${RESET}"
	echo

	# Build the config content
	local config_content="/*
 * NVD2 User Configuration
 * 
 * Generated by: install.sh
 * Date: $(date)
 * 
 * This configuration file contains paths to reference databases and
 * execution settings for the NVD2 pipeline.
 * 
 * Usage:
 *   nextflow run dhoconno/nvd -c ${config_file} --samplesheet samples.csv
 * 
 * To modify database paths or settings, edit this file directly.
 */

params {
    // Execution Environment
    // ----------------------------------------"

	# Add execution environment comment
	case "$exec_env" in
	docker)
		config_content+="
    // Running with Docker containers"
		;;
	podman)
		config_content+="
    // Running with Podman containers"
		;;
	apptainer)
		config_content+="
    // Running with Apptainer containers"
		;;
	pixi)
		config_content+="
    // Running with Pixi environment (local execution)"
		;;
	esac

	config_content+="

    // Database Paths
    // ----------------------------------------"

	# Add STAT database paths (only if configured)
	if [[ -n "$stat_path" ]]; then
		config_content+="
    // STAT database files
    stat_index      = \"${stat_path}/tree_index.dbs\"
    stat_dbss       = \"${stat_path}/tree_filter.dbss\"
    stat_annotation = \"${stat_path}/tree_filter.dbss.annotation\"
    human_virus_taxlist = \"${stat_path}/human_viruses_taxlist.txt\""
	fi

	# Add BLAST database paths (only if configured)
	if [[ -n "$blast_path" ]]; then
		config_content+="
    
    // BLAST database files
    blast_db        = \"${blast_path}\"
    blast_db_prefix = \"core_nt\""
	fi

	# Add GOTTCHA2 database paths (only if configured)
	if [[ -n "$gottcha_path" ]]; then
		config_content+="
    
    // GOTTCHA2 database files
    gottcha2_db     = \"${gottcha_path}/gottcha_db.species.fna\""
	fi

	config_content+="
}
"

	# Show preview of config
	echo -e "${DIM}Preview of configuration file:${RESET}"
	echo -e "${DIM}────────────────────────────────────────${RESET}"
	echo "$config_content" | sed 's/^/  /'
	echo -e "${DIM}────────────────────────────────────────${RESET}"
	echo

	# Ask for confirmation
	if prompt_yes_no "Write this configuration file?" "Y"; then
		echo "$config_content" >"$config_file"

		if [[ $? -eq 0 ]]; then
			print_success "Configuration file created!"
			echo
			echo -e "${GREEN}✓${RESET} Config saved to: ${CYAN}${config_file}${RESET}"
			echo
			echo -e "  ${DIM}You can edit this file anytime to update paths or settings${RESET}"
			echo
			return 0
		else
			print_error "Failed to write configuration file"
			echo
			return 1
		fi
	else
		print_info "Configuration file not created"
		echo
		echo -e "  ${DIM}You'll need to manually specify database paths when running NVD2${RESET}"
		echo
		return 1
	fi
}

# ============================================================================
# Dependency Checking
# ============================================================================

check_java() {
	# Silent check - stores result in JAVA_VERSION and returns status code
	if ! command -v java &>/dev/null; then
		JAVA_VERSION=""
		return 1
	fi

	# Get Java version
	local version_output
	version_output=$(java -version 2>&1 | head -n 1)
	local version
	# Extract version using sed (works on both macOS and Linux)
	version=$(echo "$version_output" | sed -n 's/.*version "\?\([0-9][0-9.]*\).*/\1/p')

	if [[ -z "$version" ]]; then
		version="unknown"
	fi

	JAVA_VERSION="$version"

	# Parse major version (handle both old 1.8.0 and new 11.0.0 formats)
	local major
	if [[ "$version" =~ ^1\. ]]; then
		# Old format: 1.8.0 -> major is 8
		major=$(echo "$version" | cut -d. -f2)
	else
		# New format: 11.0.0 -> major is 11
		major=$(echo "$version" | cut -d. -f1)
	fi

	# Check if major version is sufficient
	if [[ "$major" =~ ^[0-9]+$ ]] && [[ "$major" -ge 11 ]]; then
		return 0
	else
		return 2 # Version too old
	fi
}

check_nextflow() {
	# Silent check - stores result in NEXTFLOW_VERSION and returns status code
	if ! command -v nextflow &>/dev/null; then
		NEXTFLOW_VERSION=""
		return 1
	fi

	# Get Nextflow version
	local version
	version=$(nextflow -version 2>&1 | grep -i "version" | head -n 1 | sed -n 's/.*version \([0-9][0-9.]*\).*/\1/p')

	if [[ -z "$version" ]]; then
		version="unknown"
	fi

	NEXTFLOW_VERSION="$version"
	return 0
}

check_docker() {
	# Silent check - stores result in DOCKER_VERSION and returns status code
	# Return codes: 0=running, 1=not found, 2=installed but not running
	if ! command -v docker &>/dev/null; then
		DOCKER_VERSION=""
		return 1
	fi

	# Check if Docker daemon is running (with timeout to avoid hanging)
	if timeout 2 docker info &>/dev/null; then
		local version
		version=$(docker --version 2>&1 | sed -n 's/.*version \([0-9][0-9.]*\).*/\1/p')

		if [[ -z "$version" ]]; then
			version="unknown"
		fi

		DOCKER_VERSION="$version"
		return 0
	else
		DOCKER_VERSION="installed"
		return 2 # Installed but not running
	fi
}

check_podman() {
	# Silent check - stores result in PODMAN_VERSION and returns status code
	if ! command -v podman &>/dev/null; then
		PODMAN_VERSION=""
		return 1
	fi

	local version
	version=$(podman --version 2>&1 | sed -n 's/.*version \([0-9][0-9.]*\).*/\1/p')

	if [[ -z "$version" ]]; then
		version="unknown"
	fi

	PODMAN_VERSION="$version"
	return 0
}

check_apptainer() {
	# Silent check - stores result in APPTAINER_VERSION and returns status code
	# Check for both apptainer and singularity (legacy name)
	local cmd=""
	if command -v apptainer &>/dev/null; then
		cmd="apptainer"
	elif command -v singularity &>/dev/null; then
		cmd="singularity"
	else
		APPTAINER_VERSION=""
		APPTAINER_CMD=""
		return 1
	fi

	local version
	version=$($cmd --version 2>&1 | sed -n 's/.*version \([0-9][0-9.]*\).*/\1/p')

	if [[ -z "$version" ]]; then
		version="unknown"
	fi

	APPTAINER_VERSION="$version"
	APPTAINER_CMD="$cmd" # Store which command is available
	return 0
}

check_pixi() {
	# Silent check - stores result in PIXI_VERSION and returns status code
	if ! command -v pixi &>/dev/null; then
		PIXI_VERSION=""
		return 1
	fi

	local version
	version=$(pixi --version 2>&1 | sed -n 's/.*pixi \([0-9][0-9.]*\).*/\1/p')

	if [[ -z "$version" ]]; then
		version="unknown"
	fi

	PIXI_VERSION="$version"
	return 0
}

check_execution_environment() {
	echo "Checking execution environments..."
	echo

	local docker_status=0
	local podman_status=0
	local apptainer_status=0
	local pixi_status=0

	check_docker && docker_status=1 || docker_status=$?
	check_podman && podman_status=1
	check_apptainer && apptainer_status=1
	check_pixi && pixi_status=1

	# Return 0 if at least one is available
	if [[ $docker_status -eq 1 ]] || [[ $docker_status -eq 2 ]] ||
		[[ $podman_status -eq 1 ]] || [[ $apptainer_status -eq 1 ]] ||
		[[ $pixi_status -eq 1 ]]; then
		return 0
	else
		return 1
	fi
}

# ============================================================================
# Dependency Summary
# ============================================================================

print_compact_dependency_summary() {
	# Compact, single-screen dependency summary (no redundancy)
	print_header "System Check"

	echo "Dependencies:"

	# Java
	if [[ -n "$JAVA_VERSION" ]]; then
		local java_result
		check_java
		local java_check=$?
		if [[ $java_check -eq 0 ]]; then
			print_success "Java $JAVA_VERSION"
		elif [[ $java_check -eq 2 ]]; then
			print_error "Java $JAVA_VERSION ${DIM}(version too old, need 11+)${RESET}"
		fi
	else
		print_error "Java not found ${DIM}(required)${RESET}"
	fi

	# Nextflow
	if [[ -n "$NEXTFLOW_VERSION" ]]; then
		print_success "Nextflow $NEXTFLOW_VERSION"
	else
		print_error "Nextflow not found ${DIM}(required)${RESET}"
	fi

	# Execution environments
	local has_execution=false

	# Docker
	if [[ -n "$DOCKER_VERSION" ]]; then
		if [[ "$DOCKER_VERSION" == "installed" ]]; then
			print_warning "Docker installed but daemon not running"
		else
			print_success "Docker $DOCKER_VERSION"
			has_execution=true
		fi
	fi

	# Podman
	if [[ -n "$PODMAN_VERSION" ]]; then
		print_success "Podman $PODMAN_VERSION"
		has_execution=true
	fi

	# Apptainer
	if [[ -n "$APPTAINER_VERSION" ]]; then
		if [[ "$APPTAINER_CMD" == "singularity" ]]; then
			print_success "Singularity $APPTAINER_VERSION ${DIM}(legacy Apptainer)${RESET}"
		else
			print_success "Apptainer $APPTAINER_VERSION"
		fi
		has_execution=true
	fi

	# Pixi
	if [[ -n "$PIXI_VERSION" ]]; then
		print_success "Pixi $PIXI_VERSION ${DIM}(GOTTCHA2 only)${RESET}"
		has_execution=true
	fi

	# Show what execution environment will be used or warn if none
	echo
	if [[ "$has_execution" == "true" ]]; then
		# Determine which will be selected
		local selected_env=""
		if [[ -n "$DOCKER_VERSION" ]] && [[ "$DOCKER_VERSION" != "installed" ]]; then
			selected_env="Docker"
		elif [[ -n "$PODMAN_VERSION" ]]; then
			selected_env="Podman"
		elif [[ -n "$APPTAINER_VERSION" ]]; then
			selected_env="Apptainer"
		elif [[ -n "$PIXI_VERSION" ]]; then
			selected_env="Pixi"
		fi

		if [[ -n "$selected_env" ]]; then
			echo -e "${BOLD}Execution Environment:${RESET} ${GREEN}$selected_env${RESET}"
		fi
	else
		print_error "No execution environment available"
		echo -e "${DIM}Need at least one of: Docker, Podman, Apptainer, or Pixi${RESET}"
	fi
}

# ============================================================================
# Mode: Interactive Installation
# ============================================================================

interactive_mode() {
	print_header "NVD2 Interactive Installation"

	echo "This installer will guide you through setting up NVD2, a Nextflow pipeline"
	echo "for metagenomic analysis and viral detection."
	echo
	echo "Steps:"
	echo "  1. Check system dependencies (Java, Nextflow, container runtime)"
	echo "  2. Configure execution environment"
	echo "  3. Set up reference databases (optional)"
	echo "  4. Generate configuration file"
	echo
	echo "Estimated time: 5-10 minutes"
	echo -e "${DIM}Note: Database downloads are large (100s of GB) and optional${RESET}"
	echo

	if [[ "$DRY_RUN" == "true" ]]; then
		print_dry_run "Running in dry-run mode - no files will be modified"
		echo
	fi

	prompt_continue

	# Check dependencies silently (results stored in global variables)
	check_java || true
	check_nextflow || true
	check_docker || true
	check_podman || true
	check_apptainer || true
	check_pixi || true

	# Print compact summary (single screen, no redundancy)
	print_compact_dependency_summary

	echo

	# Handle critical dependency errors before proceeding
	local java_check_result
	check_java
	java_check_result=$?

	if [[ $java_check_result -eq 2 ]]; then
		# Java version too old
		echo
		if [[ "$DRY_RUN" != "true" ]]; then
			handle_old_java_version "$JAVA_VERSION"
			exit 1
		else
			print_dry_run "Would show Java version error and exit"
		fi
	elif [[ $java_check_result -eq 1 ]]; then
		# Java not found
		if [[ "$DRY_RUN" != "true" ]]; then
			echo
			print_error "Java is required but not found"
			echo
			echo "Please install Java 11 or newer and re-run this installer."
			echo
			handle_old_java_version "not installed"
			exit 1
		fi
	fi

	# Ask if user wants to proceed
	if prompt_yes_no "Would you like to proceed with the installation wizard?" "y"; then
		print_success "Great! Let's continue..."
		echo
	else
		print_info "Installation cancelled by user"
		echo
		echo -e "  ${DIM}You can run this script again anytime to continue.${RESET}"
		echo
		exit 0
	fi

	# Phase 2g: Execution environment selection
	print_header "Execution Environment Selection"

	echo "NVD2 needs an execution environment to provide the bioinformatics tools."
	echo

	echo "Based on what we detected on your system, please choose which to use:"
	echo

	# Build options list dynamically based on what's available
	local env_options=()
	local env_values=()

	if [[ -n "$DOCKER_VERSION" ]] && [[ "$DOCKER_VERSION" != "installed" ]]; then
		env_options+=("${GREEN}Docker${RESET} ${DIM}(detected, ready to use)${RESET}
        ${DIM}Works with all workflows (STAT+BLAST, GOTTCHA2)${RESET}")
		env_values+=("docker")
	elif [[ -n "$DOCKER_VERSION" ]] && [[ "$DOCKER_VERSION" == "installed" ]]; then
		env_options+=("${YELLOW}Docker${RESET} ${DIM}(detected but not running)${RESET}
        ${DIM}Start Docker Desktop first, then re-run this installer${RESET}")
		env_values+=("docker_not_running")
	fi

	if [[ -n "$PODMAN_VERSION" ]]; then
		env_options+=("${GREEN}Podman${RESET} ${DIM}(detected, ready to use)${RESET}
        ${DIM}Works with all workflows (STAT+BLAST, GOTTCHA2)${RESET}")
		env_values+=("podman")
	fi

	if [[ -n "$APPTAINER_VERSION" ]]; then
		env_options+=("${GREEN}Apptainer${RESET} ${DIM}(detected, ready to use)${RESET}
        ${DIM}Works with all workflows (STAT+BLAST, GOTTCHA2)${RESET}")
		env_values+=("apptainer")
	fi

	if [[ -n "$PIXI_VERSION" ]]; then
		env_options+=("${GREEN}Pixi${RESET} ${DIM}(detected, ready to use)${RESET}
        ${DIM}⚠️  GOTTCHA2 only - STAT tool not available${RESET}")
		env_values+=("pixi")
	fi

	# Always offer skip option
	env_options+=("${DIM}Skip - I'll configure this manually later${RESET}")
	env_values+=("skip")

	# Show menu and get choice
	local choice
	choice=$(prompt_choice "Which execution environment would you like to use?" "${env_options[@]}")

	# Convert choice number to value
	local selected_env="${env_values[$((choice - 1))]}"

	echo -e "${BOLD}Selected:${RESET} ${env_options[$((choice - 1))]}" >&2
	echo >&2

	# Handle the selection
	case "$selected_env" in
	docker)
		print_success "Using Docker"
		echo -e "  ${DIM}Configuration will be set up for Docker${RESET}"
		echo
		;;
	docker_not_running)
		echo
		if [[ "$DRY_RUN" == "true" ]]; then
			print_dry_run "Would prompt to start Docker"
			print_warning "Docker not running - skipping in dry-run mode"
			selected_env="docker"
		else
			# Offer to help user start Docker
			if handle_docker_not_running; then
				print_success "Using Docker"
				selected_env="docker"
			else
				print_info "Continuing without Docker"
				echo -e "  ${DIM}You can configure Docker later or use a different execution environment${RESET}"
				echo
				# Let them continue anyway - they might configure manually
				selected_env="docker"
			fi
		fi
		echo
		;;
	podman)
		print_success "Using Podman"
		echo -e "  ${DIM}Configuration will be set up for Podman${RESET}"
		echo
		;;
	apptainer)
		print_success "Using Apptainer"
		echo -e "  ${DIM}Configuration will be set up for Apptainer${RESET}"
		echo
		;;
	pixi)
		print_success "Using Pixi"
		echo
		print_warning "Note: STAT+BLAST workflow will not be available"
		echo -e "  ${DIM}Only GOTTCHA2 workflow will work${RESET}"
		echo
		;;
	skip)
		print_info "Execution environment configuration skipped"
		echo -e "${DIM}You'll need to configure this manually in your Nextflow config${RESET}"
		echo
		;;
	esac

	# Phase 3a: Pipeline Source Code Setup
	local pipeline_path
	pipeline_path=$(setup_pipeline_source)

	# Phase 3b: Database Setup
	print_header "Database Setup"

	echo "NVD2 uses reference databases for analysis. These are large files (400-500GB each)."
	echo
	echo -e "${BOLD}Important:${RESET} Databases are selected at runtime with ${CYAN}--tools${RESET} flag:"
	echo -e "${DIM}• ${CYAN}--tools stat_blast${RESET}${DIM} uses STAT + BLAST databases${RESET}"
	echo -e "${DIM}• ${CYAN}--tools gottcha${RESET}${DIM} uses GOTTCHA2 database${RESET}"
	echo -e "${DIM}• ${CYAN}--tools all${RESET}${DIM} uses all configured databases${RESET}"
	echo

	# Ask: do they have databases already, want to download, or skip?
	local db_mode_options=(
		"${GREEN}I already have the databases${RESET} ${DIM}(just configure paths)${RESET}"
		"${CYAN}Help me set up database downloads${RESET} ${DIM}(with disk space checks)${RESET}"
		"${DIM}Skip - I'll configure this manually later${RESET}"
	)

	local db_mode_choice
	db_mode_choice=$(prompt_choice "Do you have the NVD2 databases already?" "${db_mode_options[@]}")

	local stat_path=""
	local blast_path=""
	local gottcha_path=""

	case "$db_mode_choice" in
	1)
		# User has databases - configure paths without disk checks
		echo
		print_success "Configuring paths for existing databases"
		echo
		local db_paths
		db_paths=$(configure_all_databases "have")
		IFS='|' read -r stat_path blast_path gottcha_path <<<"$db_paths"
		;;
	2)
		# User wants to download - configure paths WITH disk checks
		echo
		print_success "Setting up database download locations"
		echo
		local db_paths
		db_paths=$(configure_all_databases "download")
		IFS='|' read -r stat_path blast_path gottcha_path <<<"$db_paths"

		# Now download and extract the databases
		if [[ -n "$stat_path" ]] || [[ -n "$blast_path" ]] || [[ -n "$gottcha_path" ]]; then
			echo
			download_and_extract_all "$stat_path" "$blast_path" "$gottcha_path"
		fi
		;;
	3)
		# Skip
		echo
		print_info "Database configuration skipped"
		echo
		echo -e "${DIM}You'll need to manually configure database paths in:${RESET}"
		echo -e "${CYAN}~/.nvd2/config/user.config${RESET}"
		echo
		;;
	esac

	# Summary of what was configured
	if [[ -n "$stat_path" ]] || [[ -n "$blast_path" ]] || [[ -n "$gottcha_path" ]]; then
		echo
		print_header "Database Configuration Summary"

		if [[ -n "$stat_path" ]]; then
			echo -e "${GREEN}✓${RESET} STAT database:     ${CYAN}${stat_path}${RESET}"
		else
			echo -e "  ${DIM}○ STAT database:     Not configured${RESET}"
		fi

		if [[ -n "$blast_path" ]]; then
			echo -e "${GREEN}✓${RESET} BLAST database:    ${CYAN}${blast_path}${RESET}"
		else
			echo -e "  ${DIM}○ BLAST database:    Not configured${RESET}"
		fi

		if [[ -n "$gottcha_path" ]]; then
			echo -e "${GREEN}✓${RESET} GOTTCHA2 database: ${CYAN}${gottcha_path}${RESET}"
		else
			echo -e "  ${DIM}○ GOTTCHA2 database: Not configured${RESET}"
		fi

		echo
	fi

	# Phase 3e: Generate Configuration File
	if [[ -n "$stat_path" ]] || [[ -n "$blast_path" ]] || [[ -n "$gottcha_path" ]]; then
		if [[ "$selected_env" != "skip" ]]; then
			echo
			generate_user_config "$stat_path" "$blast_path" "$gottcha_path" "$selected_env"
		fi
	fi

	# Phase 3f: Final Summary and Next Steps
	if [[ "$DRY_RUN" != "true" ]]; then
		print_header "Installation Complete"

		echo "NVD2 has been configured and is ready to use."
		echo

		echo -e "${BOLD}Usage Example:${RESET}"
		echo

		# Determine pipeline path
		local run_path="dhoconno/nvd"
		if [[ -n "$pipeline_path" ]]; then
			run_path="$pipeline_path"
		fi

		# Determine profile
		local profile_flag=""
		case "$selected_env" in
		docker)
			profile_flag="-profile docker"
			;;
		podman)
			profile_flag="-profile docker" # podman uses docker profile
			;;
		apptainer)
			profile_flag="-profile apptainer"
			;;
		pixi)
			profile_flag="" # pixi uses standard profile
			;;
		esac

		# Show ONE simple example using config file
		echo -e "  ${CYAN}nextflow run ${run_path} \\"
		if [[ -n "$profile_flag" ]]; then
			echo -e "    ${profile_flag} \\"
		fi
		echo -e "    -c ~/.nvd2/config/user.config \\"
		echo -e "    --samplesheet your_samples.csv \\"
		echo -e "    --experiment_id my_experiment \\"

		# Suggest appropriate --tools flag
		if [[ -n "$stat_path" ]] && [[ -n "$blast_path" ]] && [[ -n "$gottcha_path" ]]; then
			echo -e "    --tools all${RESET}  ${DIM}# or stat_blast or gottcha${RESET}"
		elif [[ -n "$stat_path" ]] && [[ -n "$blast_path" ]]; then
			echo -e "    --tools stat_blast${RESET}"
		elif [[ -n "$gottcha_path" ]]; then
			echo -e "    --tools gottcha${RESET}"
		fi
		echo
		echo

		echo -e "${BOLD}Next Steps:${RESET}"
		echo
		echo -e "  ${GREEN}1.${RESET} Create a samplesheet: ${CYAN}assets/example_samplesheet.csv${RESET}"
		echo -e "  ${GREEN}2.${RESET} Run the command above with your samplesheet"
		echo -e "  ${GREEN}3.${RESET} Results will be in: ${CYAN}./results${RESET}"
		echo

		echo -e "${DIM}Documentation: ${CYAN}https://github.com/dhoconno/nvd${RESET}"
		echo
	else
		print_dry_run "Would generate configuration file for: ${selected_env}"
		echo
		print_success "Dry-run complete! No changes were made."
	fi
}

# ============================================================================
# Mode: Verify Installation
# ============================================================================

verify_mode() {
	print_header "Verify NVD2 Installation"

	echo "Verifying NVD2 installation and configuration..."
	echo

	local issues_found=0

	# Check dependencies
	echo -e "${BOLD}Checking Dependencies:${RESET}"
	echo

	check_java
	local java_status=$?
	if [[ $java_status -eq 0 ]]; then
		print_success "Java $JAVA_VERSION"
	elif [[ $java_status -eq 2 ]]; then
		print_error "Java $JAVA_VERSION (too old, need 11+)"
		((issues_found++))
	else
		print_error "Java not found"
		((issues_found++))
	fi

	check_nextflow
	if [[ $? -eq 0 ]]; then
		print_success "Nextflow $NEXTFLOW_VERSION"
	else
		print_error "Nextflow not found"
		((issues_found++))
	fi

	# Check execution environment
	check_docker
	local docker_status=$?
	check_podman
	local podman_status=$?
	check_apptainer
	local apptainer_status=$?
	check_pixi
	local pixi_status=$?

	local has_runtime=false
	if [[ $docker_status -eq 0 ]]; then
		print_success "Docker $DOCKER_VERSION"
		has_runtime=true
	elif [[ $docker_status -eq 2 ]]; then
		print_warning "Docker installed but not running"
	fi

	if [[ $podman_status -eq 0 ]]; then
		print_success "Podman $PODMAN_VERSION"
		has_runtime=true
	fi

	if [[ $apptainer_status -eq 0 ]]; then
		print_success "Apptainer $APPTAINER_VERSION"
		has_runtime=true
	fi

	if [[ $pixi_status -eq 0 ]]; then
		print_success "Pixi $PIXI_VERSION"
		has_runtime=true
	fi

	if [[ "$has_runtime" == "false" ]]; then
		print_error "No execution environment available"
		((issues_found++))
	fi

	echo

	# Check configuration
	echo -e "${BOLD}Checking Configuration:${RESET}"
	echo

	local config_file="${HOME}/.nvd2/config/user.config"
	if [[ -f "$config_file" ]]; then
		print_success "Config file exists: $config_file"

		# Try to validate it's readable
		if grep -q "params" "$config_file" 2>/dev/null; then
			print_success "Config file format looks valid"
		else
			print_warning "Config file may be malformed"
			((issues_found++))
		fi
	else
		print_warning "Config file not found: $config_file"
		echo -e "  ${DIM}Run installer to generate configuration${RESET}"
	fi

	echo

	# Check databases (if config exists)
	if [[ -f "$config_file" ]]; then
		echo -e "${BOLD}Checking Databases:${RESET}"
		echo

		# Extract database paths from config
		local stat_index=$(grep "stat_index" "$config_file" 2>/dev/null | sed 's/.*= *"\([^"]*\)".*/\1/')
		local blast_db=$(grep "blast_db[^_]" "$config_file" 2>/dev/null | sed 's/.*= *"\([^"]*\)".*/\1/')
		local gottcha_db=$(grep "gottcha2_db" "$config_file" 2>/dev/null | sed 's/.*= *"\([^"]*\)".*/\1/')

		if [[ -n "$stat_index" ]]; then
			if [[ -f "$stat_index" ]]; then
				print_success "STAT database found: $(dirname "$stat_index")"
			else
				print_error "STAT database not found: $stat_index"
				((issues_found++))
			fi
		fi

		if [[ -n "$blast_db" ]]; then
			if [[ -d "$blast_db" ]]; then
				print_success "BLAST database found: $blast_db"
			else
				print_error "BLAST database not found: $blast_db"
				((issues_found++))
			fi
		fi

		if [[ -n "$gottcha_db" ]]; then
			if [[ -f "$gottcha_db" ]]; then
				print_success "GOTTCHA2 database found: $gottcha_db"
			else
				print_error "GOTTCHA2 database not found: $gottcha_db"
				((issues_found++))
			fi
		fi

		if [[ -z "$stat_index" ]] && [[ -z "$blast_db" ]] && [[ -z "$gottcha_db" ]]; then
			print_warning "No database paths configured"
		fi

		echo
	fi

	# Summary
	print_header "Verification Summary"

	if [[ $issues_found -eq 0 ]]; then
		print_success "No issues found - NVD2 is properly configured"
		echo
		echo "You're ready to run NVD2!"
		echo
		exit 0
	else
		print_warning "Found $issues_found issue(s) that need attention"
		echo
		echo "Recommendations:"
		if ! check_java &>/dev/null; then
			echo "  • Install Java 11 or newer"
		fi
		if ! check_nextflow &>/dev/null; then
			echo "  • Install Nextflow"
		fi
		if [[ "$has_runtime" == "false" ]]; then
			echo "  • Install Docker, Podman, Apptainer, or Pixi"
		fi
		echo
		exit 1
	fi
}

# ============================================================================
# Mode: Uninstall
# ============================================================================

uninstall_mode() {
	print_header "Uninstall NVD2"

	echo "This will remove NVD2 configuration and optionally databases."
	echo
	echo -e "${YELLOW}${BOLD}Warning:${RESET}${YELLOW} This action cannot be undone.${RESET}"
	echo

	# Show what will be removed
	echo "The following will be checked for removal:"
	echo

	local config_dir="${HOME}/.nvd2"
	local nextflow_assets="${HOME}/.nextflow/assets/dhoconno/nvd"

	local items_to_remove=()
	local items_found=false

	# Check config directory
	if [[ -d "$config_dir" ]]; then
		echo -e "  ${CHECKMARK} Config directory: ${CYAN}${config_dir}${RESET}"
		items_to_remove+=("$config_dir")
		items_found=true
	else
		echo -e "  ${DIM}• Config directory: $config_dir (not found)${RESET}"
	fi

	# Check cached pipeline
	if [[ -d "$nextflow_assets" ]]; then
		echo -e "  ${CHECKMARK} Cached pipeline: ${CYAN}${nextflow_assets}${RESET}"
		items_to_remove+=("$nextflow_assets")
		items_found=true
	else
		echo -e "  ${DIM}• Cached pipeline: $nextflow_assets (not found)${RESET}"
	fi

	# Check for database paths in config
	local config_file="${config_dir}/config/user.config"
	local db_paths=()

	if [[ -f "$config_file" ]]; then
		echo
		echo "Database locations found in config:"

		# Extract database parent directories
		local stat_db=$(grep "stat_index" "$config_file" 2>/dev/null | sed 's/.*= *"\([^"]*\)".*/\1/' | xargs dirname 2>/dev/null)
		local blast_db=$(grep "blast_db[^_]" "$config_file" 2>/dev/null | sed 's/.*= *"\([^"]*\)".*/\1/')
		local gottcha_db=$(grep "gottcha2_db" "$config_file" 2>/dev/null | sed 's/.*= *"\([^"]*\)".*/\1/' | xargs dirname 2>/dev/null)

		[[ -n "$stat_db" ]] && [[ -d "$stat_db" ]] && echo -e "  • STAT: ${CYAN}${stat_db}${RESET}" && db_paths+=("$stat_db")
		[[ -n "$blast_db" ]] && [[ -d "$blast_db" ]] && echo -e "  • BLAST: ${CYAN}${blast_db}${RESET}" && db_paths+=("$blast_db")
		[[ -n "$gottcha_db" ]] && [[ -d "$gottcha_db" ]] && echo -e "  • GOTTCHA2: ${CYAN}${gottcha_db}${RESET}" && db_paths+=("$gottcha_db")
	fi

	echo

	if [[ "$items_found" == "false" ]] && [[ ${#db_paths[@]} -eq 0 ]]; then
		print_info "No NVD2 files found on this system"
		echo
		exit 0
	fi

	# Confirm removal of config and cached pipeline
	if [[ ${#items_to_remove[@]} -gt 0 ]]; then
		echo -e "${BOLD}Remove config and cached pipeline?${RESET}"
		if ! prompt_yes_no "Remove ${#items_to_remove[@]} item(s)?" "n"; then
			print_info "Uninstall cancelled"
			exit 0
		fi

		# Remove items
		for item in "${items_to_remove[@]}"; do
			echo "Removing: $item"
			if rm -rf "$item" 2>/dev/null; then
				print_success "Removed: $item"
			else
				print_error "Failed to remove: $item"
				echo -e "  ${DIM}You may need to remove it manually with: sudo rm -rf $item${RESET}"
			fi
		done

		echo
	fi

	# Ask about databases
	if [[ ${#db_paths[@]} -gt 0 ]]; then
		echo -e "${BOLD}Database Removal:${RESET}"
		echo
		echo "Found ${#db_paths[@]} database location(s)."
		echo -e "${YELLOW}Databases are large (100s of GB) and may be shared with other tools.${RESET}"
		echo

		if prompt_yes_no "Remove databases?" "n"; then
			for db_path in "${db_paths[@]}"; do
				echo
				echo -e "Remove: ${CYAN}${db_path}${RESET}"
				echo -e "${DIM}Size: $(du -sh "$db_path" 2>/dev/null | cut -f1)${RESET}"

				if prompt_yes_no "Confirm removal of this database?" "n"; then
					echo "Removing: $db_path"
					if rm -rf "$db_path" 2>/dev/null; then
						print_success "Removed: $db_path"
					else
						print_error "Failed to remove: $db_path"
						echo -e "  ${DIM}You may need to remove it manually with: sudo rm -rf $db_path${RESET}"
					fi
				else
					print_info "Skipped: $db_path"
				fi
			done
		else
			print_info "Databases were not removed"
			echo
			echo "To remove databases manually:"
			for db_path in "${db_paths[@]}"; do
				echo "  rm -rf $db_path"
			done
		fi

		echo
	fi

	# Show commands for removing dependencies
	print_header "Removing Dependencies (Optional)"

	echo "NVD2 dependencies can be removed if not needed for other projects:"
	echo

	echo -e "${BOLD}Nextflow:${RESET}"
	if command -v nextflow &>/dev/null; then
		local nf_location=$(which nextflow)
		echo "  rm $nf_location"
		echo "  rm -rf ~/.nextflow"
	else
		echo "  (not installed)"
	fi

	echo
	echo -e "${BOLD}Java:${RESET}"
	if command -v java &>/dev/null; then
		if [[ "$OSTYPE" == "darwin"* ]]; then
			echo "  brew uninstall openjdk@17  # if installed via Homebrew"
		else
			echo "  sudo apt remove openjdk-17-jdk  # Ubuntu/Debian"
			echo "  sudo yum remove java-17-openjdk  # RHEL/CentOS"
		fi
	else
		echo "  (not installed)"
	fi

	echo
	echo -e "${BOLD}Docker:${RESET}"
	if command -v docker &>/dev/null; then
		if [[ "$OSTYPE" == "darwin"* ]]; then
			echo "  Uninstall Docker Desktop from Applications folder"
		else
			echo "  sudo apt remove docker-ce  # Ubuntu/Debian"
			echo "  sudo yum remove docker-ce  # RHEL/CentOS"
		fi
	else
		echo "  (not installed)"
	fi

	echo
	print_success "Uninstall complete"
}

# ============================================================================
# Mode: Non-Interactive (CI/CD)
# ============================================================================

non_interactive_mode() {
	print_header "NVD2 Prerequisite Check (Non-Interactive)"

	echo "Checking prerequisites for NVD2..."
	echo

	local exit_code=0

	# Check Java
	if ! check_java; then
		exit_code=1
	fi

	# Check Nextflow
	if ! check_nextflow; then
		exit_code=2
	fi

	# Check execution environment
	if ! check_execution_environment; then
		exit_code=3
	fi

	echo

	if [[ $exit_code -eq 0 ]]; then
		print_success "All prerequisites met"
		echo
		print_info "System is ready to run NVD2"
		echo
	else
		print_error "Some prerequisites are missing or incorrect"
		echo
		echo "Install missing dependencies and try again."
		echo
		echo "Required:"
		echo "- Java 11+"
		echo "- Nextflow"
		echo "- At least one of: Docker, Podman, Apptainer, or Pixi"
		echo
	fi

	exit $exit_code
}

# ============================================================================
# Argument Parsing
# ============================================================================

parse_arguments() {
	# Track which mode to run
	local mode="interactive"

	while [[ $# -gt 0 ]]; do
		case $1 in
		--dry-run | --test | -d)
			DRY_RUN=true
			echo
			echo -e "${COLOR_WARNING}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
			echo -e "${COLOR_WARNING}  DRY RUN MODE - No actual changes will be made${RESET}"
			echo -e "${COLOR_WARNING}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${RESET}"
			shift
			;;
		--non-interactive | --ci | -y | --yes)
			NON_INTERACTIVE=true
			mode="non-interactive"
			shift
			;;
		--verify)
			mode="verify"
			shift
			;;
		--uninstall)
			mode="uninstall"
			shift
			;;
		--quiet)
			QUIET=true
			shift
			;;
		-h | --help)
			show_help
			exit 0
			;;
		-v | --version)
			show_version
			exit 0
			;;
		*)
			echo -e "${COLOR_ERROR}Error: Unknown option: $1${RESET}" >&2
			echo
			show_help
			exit 1
			;;
		esac
	done

	# Store mode for main function
	echo "$mode"
}

# ============================================================================
# Main Entry Point
# ============================================================================

main() {
	# Parse command-line arguments and get mode
	local mode
	mode=$(parse_arguments "$@") || mode=""

	# Route to appropriate mode
	case "$mode" in
	interactive)
		interactive_mode
		;;
	non-interactive)
		non_interactive_mode
		;;
	verify)
		verify_mode
		;;
	uninstall)
		uninstall_mode
		;;
	"")
		# Empty mode means --help or --version was called (already handled)
		exit 0
		;;
	*)
		if [[ -n "$mode" ]]; then
			echo "Error: Unknown mode: $mode" >&2
		fi
		exit 1
		;;
	esac
}

# Run main function with all arguments
main "$@"
