#!/usr/bin/env python
# Note that this script is designed to be run on CHTC, so the paths are specific to that environment.
# Use the master branch for the NVD repo, not the nvd-chtc branch, for typical use

import os
import csv
import re
import sys
import subprocess
import argparse
from pathlib import Path

def parse_args():
    """
    Parse command-line arguments and return an object with
    the folder_path, apptainer_path, and profile.
    """
    parser = argparse.ArgumentParser(
        description="Prepare and run a CHTC HTC job for NVD workflows."
    )
    parser.add_argument(
        "folder_path",
        type=str,
        help="Path to the folder containing FASTQ files (e.g. /staging/groups/oconnor_group/30906)."
    )
    parser.add_argument(
        "apptainer_path",
        type=str,
        help="Path to the Apptainer image (e.g. /home/dhoconno/apptainer_images/nvd.30572.sif)."
    )
    parser.add_argument(
        "--profile",
        type=str,
        choices=["normal", "testing"],
        default="normal",
        help="Profile for configuration. Choices are 'normal' or 'testing'. Default: 'normal'."
    )
    return parser.parse_args()

def sanitize_filename(filename):
    """
    Replace any character that is not a letter, digit, underscore, or dot with an underscore.
    """
    sanitized = re.sub(r'[^a-zA-Z0-9_.]', '_', filename)
    if sanitized != filename:
        print(f"[DEBUG] Sanitized filename: '{filename}' -> '{sanitized}'")
    return sanitized

def ensure_logs_directory():
    """
    Ensure that a 'logs' directory exists in the current working directory.
    """
    logs_dir = Path("logs")
    if not logs_dir.exists():
        logs_dir.mkdir(parents=True)
        print(f"[INFO] Created logs directory: {logs_dir}")
    else:
        print(f"[INFO] Logs directory already exists: {logs_dir}")

def extract_sample_info(folder_path_input):
    """
    Extract the experiment number and sample information from the given folder path.

    Args:
        folder_path_input (Path): Path object to the folder containing FASTQ files.

    Returns:
        tuple: (experiment_number (int), sample_dict (dict))
    """
    path_resolved = folder_path_input.resolve()
    experiment_number = None

    print(f"[DEBUG] Resolving experiment number from base path: {path_resolved}")
    
    # Traverse the path to find the first numeric directory component
    for part in reversed(path_resolved.parts):
        if part.isdigit():
            experiment_number = int(part)
            break

    if experiment_number is None:
        msg = f"[ERROR] No numeric experiment number found in the path: {folder_path_input}"
        print(msg)
        raise ValueError(msg)

    print(f"[INFO] Extracted experiment_number: {experiment_number}")

    # Updated pattern to match both naming conventions
    fastq_pattern = re.compile(r"(.+)_R(\d)\.fastq\.gz|(.+)\.R(\d)\.fq\.gz")
    sample_dict = {}

    for root, _, files in os.walk(path_resolved):
        root_path = Path(root)
        for file in files:
            sanitized_file = sanitize_filename(file)
            original_path = root_path / file
            new_path = root_path / sanitized_file

            if original_path != new_path:
                try:
                    original_path.rename(new_path)
                    print(f"[INFO] Renamed '{original_path}' to '{new_path}'")
                except OSError as e:
                    print(f"[WARN] Error renaming '{original_path}' to '{new_path}': {e}")
                    continue  # Skip to the next file

            match = fastq_pattern.match(sanitized_file)
            if match:
                # Extracting sample name and read number from the matched groups
                sample_name = match.group(1) or match.group(3)  # Handle both patterns
                read_number = match.group(2) or match.group(4)

                if not read_number:
                    print(f"[WARN] Read number not found in file: '{sanitized_file}'")
                    continue  # Skip files without a read number

                # Construct the user-visible path by joining the base path with the relative path
                try:
                    rel_path = new_path.relative_to(path_resolved)
                    abs_path_user_visible = folder_path_input / rel_path
                    abs_path_user_visible_str = str(abs_path_user_visible)
                except ValueError:
                    # In case new_path is not under path_resolved, which shouldn't happen
                    print(f"[ERROR] File '{new_path}' is not under the base path '{path_resolved}'")
                    continue

                if sample_name not in sample_dict:
                    sample_dict[sample_name] = {'R1': [], 'R2': []}

                if read_number in ['1', '2']:
                    sample_dict[sample_name][f'R{read_number}'].append(abs_path_user_visible_str)
                    print(f"[DEBUG] Found file: '{abs_path_user_visible_str}' as R{read_number} for sample '{sample_name}'")
                else:
                    print(f"[WARN] Unexpected read number '{read_number}' in file: '{sanitized_file}'")

    return experiment_number, sample_dict

def write_csv_file(experiment_number, sample_dict):
    """
    Write the sample information to a CSV file without a header.

    Args:
        experiment_number (int): The experiment number.
        sample_dict (dict): Dictionary containing sample information.
    """
    csv_file = f"{experiment_number}.csv"
    print(f"[INFO] Writing CSV file (no header): '{csv_file}'...")
    try:
        with open(csv_file, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # No header row (as requested)
            for sample, paths in sample_dict.items():
                r1 = paths['R1'][0] if paths['R1'] else ''
                r2 = paths['R2'][0] if paths['R2'] else ''
                writer.writerow([sample, r1, r2])
        print(f"[INFO] CSV file '{csv_file}' written successfully (no header).")
    except Exception as e:
        msg = f"[ERROR] Error writing CSV file '{csv_file}': {e}"
        print(msg)
        sys.exit(1)

def write_config_file(experiment_number, profile):
    """
    Write the configuration YAML file based on the provided profile.

    Args:
        experiment_number (int): The experiment number.
        profile (str): 'normal' or 'testing'.
        Testing uses a virus-only database for faster testing. See exp 30906.
    """
    config_file = f"{experiment_number}.config.yaml"
    
    print(f"[DEBUG] Writing config file for experiment_number={experiment_number} with profile={profile}...")

    if profile == 'testing':
        # Testing configuration content
        config_content = f"""global:
    experiment: '{experiment_number}'
    out_dir: '/staging/groups/oconnor_group/{experiment_number}/output'
    blast_db_name: 'ref_viruses_rep_genomes'
    blast_db_path: 'resources/ref_viruses_rep_genomes'
    stat_dbss: 'resources/viruses.filtered.20240625.dbss'
    stat_index: 'resources/viruses.filtered.20240625.dbs'
    stat_annotation: 'resources/viruses.filtered.20240625.dbss.annotation'
    resources_tar_zst: '/staging/groups/oconnor_group/nvd/30906-nvd-testing-resources.tar.zst'
    human_virus_taxlist: 'resources/human_viruses_taxlist.txt'
    gettax_sqlite_path: 'resources/gettax.sqlite'
    labkey_server: 'XXX'
    labkey_project_name: 'XXX'
    labkey_api_key: 'XXX'
    labkey_username: 'XXX'
    labkey_password: 'XXX'
    webdav_url: 'XXX'
    """
    else:
        # Normal configuration content
        config_content = f"""global:
  experiment: '{experiment_number}'
  out_dir: '/staging/groups/oconnor_group/{experiment_number}/output'
  blast_db_name: 'core_nt.20240803'
  blast_db_path: 'resources/core_nt'
  stat_dbss: 'resources/tree_filter.20240830.dbss'
  stat_index: 'resources/tree_index.20240830.dbs'
  stat_annotation: 'resources/tree_filter.20240830.dbss.annotation'
  resources_tar_zst: '/staging/groups/oconnor_group/nvd/resources.20241204.tar.zst'
  human_virus_taxlist: 'resources/human_viruses_taxlist.txt'
  gettax_sqlite_path: 'resources/gettax.sqlite'
  labkey_server: 'XXX'
  labkey_project_name: 'XXX'
  labkey_api_key: 'XXX'
  labkey_username: 'XXX'
  labkey_password: 'XXX'
  webdav_url: 'XXX'
"""

    print("[DEBUG] Final config content:\n", config_content)

    print(f"[INFO] Writing configuration file: '{config_file}'...")
    try:
        with open(config_file, 'w') as f:
            f.write(config_content)
        print(f"[INFO] Config file '{config_file}' written successfully.")
    except Exception as e:
        msg = f"[ERROR] Error writing config file '{config_file}': {e}"
        print(msg)
        sys.exit(1)

def write_run_script(experiment_number):
    """
    Write the run script shell file.

    Args:
        experiment_number (int): The experiment number.
    """
    run_script_file = f"{experiment_number}.run.sh"
    run_script_content = f"""#! /bin/bash

# Initialize variables
name=""
r1_fastq=""
r2_fastq=""
run_id=""
config_file="{experiment_number}.config.yaml"

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --name) name="$2"; shift ;;
        --r1_fastq) r1_fastq="$2"; shift ;;
        --r2_fastq) r2_fastq="$2"; shift ;;
        --run_id) run_id="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Print received values
echo "Sample Name: $name"
echo "Run ID: $run_id"

# Append name and files to config.yaml
yaml_entry="- name: $name"
[[ -n "$r1_fastq" ]] && yaml_entry="$yaml_entry"$'\\n'"  r1_fastq: $r1_fastq"
[[ -n "$r2_fastq" ]] && yaml_entry="$yaml_entry"$'\\n'"  r2_fastq: $r2_fastq"

if grep -q "^samples:" "$config_file"; then
    echo "$yaml_entry" >> "$config_file"
else
    echo -e "\\nsamples:" >> "$config_file"
    echo "$yaml_entry" >> "$config_file"
fi

echo "New sample added to $config_file:"
echo "$yaml_entry"

# Clone GitHub NVD repo workflow files to execute node
# Use branch nvd-chtc created specifically to work with this script
git init
git remote add origin https://github.com/dhoconno/nvd.git
git sparse-checkout init --cone
git sparse-checkout set workflow
git fetch origin
git checkout -b nvd-chtc origin/nvd-chtc

# Create config directory and copy config file
mkdir -p config
cp "$config_file" config/config.yaml

snakemake --cores 4 --debug-dag --verbose --config run_id=$run_id
"""

    print(f"[INFO] Writing run script: '{run_script_file}'...")
    try:
        with open(run_script_file, 'w') as f:
            f.write(run_script_content)
        # Make the run script executable
        os.chmod(run_script_file, 0o755)
        print(f"[INFO] Run script '{run_script_file}' written and made executable.")
    except Exception as e:
        msg = f"[ERROR] Error writing run script '{run_script_file}': {e}"
        print(msg)
        sys.exit(1)

def create_submit_file(experiment_number, apptainer_path):
    """
    Create and submit the Condor submit file.

    Args:
        experiment_number (int): The experiment number.
        apptainer_path (str): Path to the Apptainer image.
    """
    ensure_logs_directory()
    # Create a .gitignore file and add 'logs' to it
    with open(".gitignore", "w") as f:
        f.write("logs\n")
    submit_file = f"{experiment_number}.sub"
    csv_file = f"{experiment_number}.csv"
    run_script_file = f"{experiment_number}.run.sh"
    config_file = f"{experiment_number}.config.yaml"

    submit_content = f"""# Execute node type
universe = container
container_image = {apptainer_path}

+AccountingGroup = "Pathology_OConnor"
requirements = (Target.HasCHTCStaging == True)

output = logs/$(db_type)_out.$(Cluster).$(Process).out
error = logs/$(db_type)_err.$(Cluster).$(Process).err
log = logs/$(db_type)_log.$(Cluster).$(Process).log

request_cpus = 8
request_memory = 64GB
request_disk = 1000GB

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {config_file},{run_script_file}
transfer_output_files = logs

Arguments = bash {run_script_file} --name $(NAME) --r1_fastq $(R1) --r2_fastq $(R2) --run_id $(Cluster)

queue NAME, R1, R2 from {csv_file}
"""

    print(f"[INFO] Creating submit file: '{submit_file}'...")
    try:
        with open(submit_file, 'w') as f:
            f.write(submit_content)
        print(f"[INFO] Submit file '{submit_file}' created successfully.")
    except Exception as e:
        msg = f"[ERROR] Error creating submit file '{submit_file}': {e}"
        print(msg)
        sys.exit(1)

    # Submit the job to Condor
    print(f"[INFO] Submitting job using submit file '{submit_file}'...")
    try:
        result = subprocess.run(["condor_submit", submit_file], check=True, capture_output=True, text=True)
        print(f"[INFO] Job submitted successfully:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Error submitting job with Condor:\n{e.stderr}")
        sys.exit(1)

def main():
    # Parse command-line arguments
    args = parse_args()
    folder_path_input = Path(args.folder_path)
    apptainer_path = Path(args.apptainer_path).resolve()
    profile = args.profile

    print(f"[INFO] Received arguments:")
    print(f"       folder_path = {folder_path_input}")
    print(f"       apptainer_path = {apptainer_path}")
    print(f"       profile = {profile}")

    # Basic validity checks
    if not folder_path_input.is_dir():
        msg = f"[ERROR] The folder_path '{folder_path_input}' does not exist or is not a directory."
        print(msg)
        sys.exit(1)
    if not apptainer_path.is_file():
        msg = f"[ERROR] The apptainer_path '{apptainer_path}' does not exist or is not a file."
        print(msg)
        sys.exit(1)

    # Proceed with logic
    try:
        experiment_number, sample_dict = extract_sample_info(folder_path_input)
        print(f"[INFO] Writing CSV and config with experiment_number={experiment_number}, profile={profile}")
        write_csv_file(experiment_number, sample_dict)
        write_config_file(experiment_number, profile)
        write_run_script(experiment_number)
        create_submit_file(experiment_number, str(apptainer_path))
    except Exception as e:
        print(f"[ERROR] An error occurred during processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
