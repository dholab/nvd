import os
import glob
import shutil

def cleanup(local_output_dir, samples, dirs_to_remove):
    # Define patterns for files to delete using {sample} wildcard
    patterns_to_delete = [
        f"{local_output_dir}/01-prepare-input/{{sample}}.fq.gz",
        f"{local_output_dir}/08-de-novo-assemble-reads/{{sample}}_cleanup.done",
        f"{local_output_dir}/09-mask-low-complexity-contigs/{{sample}}.fa",
        f"{local_output_dir}/10-exclude-short-contigs/{{sample}}.fa",
        # Add more patterns as needed
    ]

    # Add any non-sample-specific files
    patterns_to_delete.extend([
        f"{local_output_dir}/11-list-contig-files/datasets.list",
    ])

    # Add directories to delete
    directories_to_delete = [
        f"{local_output_dir}/08-de-novo-assemble-reads/{{sample}}",
        # Add more directories as needed
    ]

    for sample in samples:
        sample_name = sample['name'].replace('-', '_')
        
        # Delete files
        for pattern in patterns_to_delete:
            sample_pattern = pattern.format(sample=sample_name)
            for file_path in glob.glob(sample_pattern):
                if os.path.exists(file_path):
                    os.remove(file_path)
        
        # Delete directories
        for dir_pattern in directories_to_delete:
            dir_path = dir_pattern.format(sample=sample_name)
            if os.path.exists(dir_path) and os.path.isdir(dir_path):
                shutil.rmtree(dir_path)

    # Remove any directories defined in dirs_to_remove
    for dir_to_remove in dirs_to_remove:
        if os.path.exists(dir_to_remove):
            shutil.rmtree(dir_to_remove)

# Snakemake parameters
local_output_dir = snakemake.params.local_output_dir
samples = snakemake.config['samples']
dirs_to_remove = snakemake.params.dirs_to_remove

# Perform cleanup
cleanup(local_output_dir, samples, dirs_to_remove)
