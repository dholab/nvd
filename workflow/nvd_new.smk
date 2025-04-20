#### filepath: nvd_new.smk
import os
import glob
import pandas as pd
import re
from snakemake.utils import min_version

# --- Minimum Snakemake version ---
min_version("6.0")

# --- Configuration ---
configfile: "config.yaml"

# --- Import Utility Functions ---
# Assuming workflow_functions.py is in a 'scripts' subdirectory
try:
    from scripts.workflow_functions import (
        get_input_files,
        get_input_type,
        get_sra_accession,
        # Add other functions from workflow_functions.py if needed
    )
except ImportError:
    raise ImportError("Could not import helper functions from scripts/workflow_functions.py. Ensure the file exists.")


# --- Sample Handling ---
# Normalize sample names from config
# Ensure 'samples' is a list of dictionaries, each with at least a 'name' key
if not isinstance(config.get('samples'), list):
     raise ValueError("config['samples'] must be a list of dictionaries.")

for sample_info in config['samples']:
    if 'name' not in sample_info:
        raise ValueError("Each entry in config['samples'] must have a 'name' key.")
    # Use original name for lookup, store normalized name for file paths
    sample_info['original_name'] = sample_info['name']
    sample_info['name'] = re.sub(r'[\s\-]+', '_', str(sample_info['name'])) # Ensure name is string and normalize

SAMPLES = {s["name"]: s for s in config["samples"]}

# --- Global Config Variables ---
# Define variables used across multiple rules for clarity
EXPERIMENT = config['global']['experiment']
RESOURCES_TAR_ZST = config['global'].get("resources_tar_zst") # Optional archive
RESOURCES_DIR_CONFIG = config['global'].get("resources_dir") # Optional existing dir
BLAST_DB_PATH = config['global']['blast_db_path'] # e.g., core_nt
STAT_DBSS = config['global']['stat_dbss'] # e.g., resources/ncbi_stat_db/database.fasta.dbss
STAT_INDEX = config['global']['stat_index'] # e.g., resources/ncbi_stat_db/database.fasta.index
STAT_ANNOTATION = config['global']['stat_annotation'] # e.g., resources/ncbi_stat_db/database.fasta.annotation.gz
HUMAN_VIRUS_TAXLIST = config['global']['human_virus_taxlist'] # e.g., resources/human_virus_taxids.txt
GOTTCHA2_DB_SPECIES = config['global']['gottcha2_db_species'] # e.g., resources/gottcha_db/20250221.gottcha_db.species.fna
BLAST_VERIFY_DB = config['blast_verify']['db_prefix'] # e.g., resources/ref/core_nt
GETTAX_SQLITE_PATH = config["global"]["gettax_sqlite_path"] # e.g., resources/gettax.sqlite

# LabKey Config (handle missing keys gracefully for non-LabKey runs)
LABKEY_SERVER = config["global"].get("labkey_server")
LABKEY_PROJECT_NAME = config["global"].get("labkey_project_name")
LABKEY_API_KEY = config["global"].get("labkey_api_key")
LABKEY_USERNAME = config["global"].get("labkey_username")
LABKEY_PASSWORD = config["global"].get("labkey_password")
WEBDAV_URL = config["global"].get("webdav_url")

# Output Directories
RESULTS_DIR = "results"
LOG_DIR = "logs"
FINAL_OUTPUT_PATH = f"{config['global']['nvd_out_dir']}/{EXPERIMENT}" # Final destination for LabKey uploads/CSVs

# --- Validation Taxa ---
VALIDATION_TAXA = config.get("validation_taxa", []) # Expects a list of taxids

# --- Resource Management ---
RESOURCES_DIR = "resources" # Local working resources directory

# --- Target Files ---
# Define all final expected outputs using the numbered structure
ALL_TARGETS = []

# Clumpify (optional, depends on input type) - Assuming paired-end Illumina default output naming
# Adjust if single-end is the primary target
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/00_clumpify/{{sample}}_R1.clumped.fastq.gz", sample=SAMPLES))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/00_clumpify/{{sample}}_R2.clumped.fastq.gz", sample=SAMPLES))
# Add single-end output if needed based on config analysis

# NVD Pipeline Outputs
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/01_prepare_input/{{sample}}.fq.gz", sample=SAMPLES))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/09_extract_human_virus_contigs/{{sample}}.fa", sample=SAMPLES))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/21_extract_unclassified_contigs/{{sample}}.unclassified.fa", sample=SAMPLES))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/22_map_reads_to_contigs/{{sample}}.bam", sample=SAMPLES))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/23_count_mapped_reads/{{sample}}_mapped_counts.txt", sample=SAMPLES))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/20_merge_annotated_blast_results/{{sample}}.txt", sample=SAMPLES))

# GOTTCHA2 Outputs
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/10_gottcha2/{{sample}}/{{sample}}.full.tsv", sample=SAMPLES))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/10_gottcha2/{{sample}}/{{sample}}.gottcha_strain.sam", sample=SAMPLES))

# GOTTCHA2 Validation Outputs
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.classification.tsv", sample=SAMPLES, taxid=VALIDATION_TAXA))
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.extract.fasta", sample=SAMPLES, taxid=VALIDATION_TAXA))

# GOTTCHA2 20 Reads Extraction Output
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.extract_20.fasta", sample=SAMPLES))

# Archive Output
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/24_create_tarball/{{sample}}.tar.zst", sample=SAMPLES))

# LabKey Upload Placeholders/Tokens
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/25_labkey/{{sample}}/nvd_blast_upload.tkn", sample=SAMPLES)) # NVD BLAST results
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/25_labkey/{{sample}}/nvd_fasta_upload.tkn", sample=SAMPLES)) # NVD contigs
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/25_labkey/{{sample}}/gottcha_upload.tkn", sample=SAMPLES)) # GOTTCHA2 results + validation summary
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/25_labkey/{{sample}}/validation_fasta_upload_{{taxid}}.tkn", sample=SAMPLES, taxid=VALIDATION_TAXA)) # GOTTCHA2 validation FASTAs
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/25_labkey/{{sample}}/extract_20_fasta_upload.tkn", sample=SAMPLES)) # GOTTCHA2 20-reads FASTAs
ALL_TARGETS.extend(expand(f"{RESULTS_DIR}/25_labkey/{{sample}}/archive_upload.tkn", sample=SAMPLES)) # Final archive upload confirmation

# --- Main Rule ---
rule all:
    input: ALL_TARGETS
    run:
        print("Workflow finished. All target files generated:")
        for target in ALL_TARGETS:
             print(f"- {target}")


# --- Resource Handling Rule ---
rule manage_resources:
    output:
        touch(os.path.join(RESOURCES_DIR, ".resources_ready")),
        stat_dbss_local=os.path.join(RESOURCES_DIR, os.path.basename(STAT_DBSS)),
        stat_index_local=os.path.join(RESOURCES_DIR, os.path.basename(STAT_INDEX)),
        stat_annotation_local=os.path.join(RESOURCES_DIR, os.path.basename(STAT_ANNOTATION)),
        human_virus_taxlist_local=os.path.join(RESOURCES_DIR, os.path.basename(HUMAN_VIRUS_TAXLIST)),
        gottcha_species_local=os.path.join(RESOURCES_DIR, os.path.basename(GOTTCHA2_DB_SPECIES) + ".mmi"),
    params:
        resources_archive=RESOURCES_TAR_ZST,
        resources_dir_config=RESOURCES_DIR_CONFIG,
        target_dir=RESOURCES_DIR
    log:
        os.path.join(LOG_DIR, "00_manage_resources.log")
    shell:
        """
        echo "Managing resources..." > {log}
        TARGET_DIR="{params.target_dir}"
        SOURCE_ARCHIVE="{params.resources_archive}"
        SOURCE_DIR="{params.resources_dir_config}"
        
        if [ -n "$SOURCE_ARCHIVE" ] && [ -f "$SOURCE_ARCHIVE" ]; then
            mkdir -p "$TARGET_DIR"
            echo "Unpacking resources from archive: $SOURCE_ARCHIVE to $TARGET_DIR" >> {log}
            cp "$SOURCE_ARCHIVE" "$TARGET_DIR/$(basename $SOURCE_ARCHIVE)" >> {log} 2>&1
            tar -I pzstd -xf "$TARGET_DIR/$(basename $SOURCE_ARCHIVE)" -C "$TARGET_DIR" >> {log} 2>&1
            rm "$TARGET_DIR/$(basename $SOURCE_ARCHIVE)" >> {log} 2>&1
            echo "Unpacking complete." >> {log}
            RES_DIR=$(realpath "$TARGET_DIR")
        elif [ -n "$SOURCE_DIR" ] && [ -d "$SOURCE_DIR" ]; then
            RES_DIR=$(realpath "$SOURCE_DIR")
            # Ensure that TARGET_DIR is a symlink to RES_DIR.
            if [ -L "$TARGET_DIR" ]; then
                echo "Target resource directory '$TARGET_DIR' is already a symlink." >> {log}
            elif [ -e "$TARGET_DIR" ]; then
                echo "Target resource directory '$TARGET_DIR' exists but is not a symlink." >> {log}
                echo "Removing it and creating a symlink to the source resources." >> {log}
                rm -rf "$TARGET_DIR"
                ln -s "$RES_DIR" "$TARGET_DIR" >> {log} 2>&1
            else
                echo "Creating symlink from $RES_DIR to $TARGET_DIR" >> {log}
                ln -s "$RES_DIR" "$TARGET_DIR" >> {log} 2>&1
            fi
        else
            echo "Error: Configuration must provide either a valid resources archive or an existing resources directory." >> {log}
            exit 1
        fi

        # Use RES_DIR for verification. This ensures that the essential files are looked up
        # in the actual source directory now accessible via the symlink TARGET_DIR.
        for f in "$RES_DIR/$(basename {output.stat_dbss_local})" "$RES_DIR/$(basename {output.stat_index_local})" "$RES_DIR/$(basename {output.stat_annotation_local})" "$RES_DIR/$(basename {output.human_virus_taxlist_local})" "$RES_DIR/$(basename {output.gottcha_species_local})"; do
            if [ ! -e "$f" ]; then
                echo "Error: Essential resource file $f not found after setup." >> {log}
                exit 1
            fi
        done
        echo "Resource setup finished successfully." >> {log}
        """

# --- Clumpify Rule (Low Priority) ---
rule clumpify:
    priority: 1
    input:
        r1=lambda w: get_input_files(w, config)[0], # Assumes get_input_files returns list [R1] or [R1, R2]
        r2=lambda w: get_input_files(w, config)[1] if len(get_input_files(w, config)) > 1 else None,
    output:
        # Define outputs based on input type (ONT/Illumina)
        # Using expand-like logic within output definition
        r1_clumped=f"{RESULTS_DIR}/00_clumpify/{{sample}}_R1.clumped.fastq.gz",
        r2_clumped=f"{RESULTS_DIR}/00_clumpify/{{sample}}_R2.clumped.fastq.gz",
        ont_clumped=f"{RESULTS_DIR}/00_clumpify/{{sample}}.clumped.fastq.gz"
    params:
        sample="{sample}",
        output_dir=f"{RESULTS_DIR}/00_clumpify",
        input_type=lambda w: get_input_type(w, config),
        extra=config["clumpify"].get("extra_params", "dedupe=t optical=t distance=10"), # Example params
        mem_gb=config["clumpify"].get("mem_gb", 100) # Memory in GB
    threads: 16 # Runs later, can use more threads
    log:
        f"{LOG_DIR}/00_clumpify/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/00_clumpify/{{sample}}.benchmark.txt"
    run:
        # Create output directory
        shell("mkdir -p {params.output_dir}")

        # Construct command based on input type
        if params.input_type in ['ont', 'sra']: # Treat SRA as single-end for clumpify unless specified otherwise
            # Ensure only the single-end output is generated
            shell(f"touch {output.r1_clumped} {output.r2_clumped}") # Touch unused outputs
            shell(
                "clumpify.sh "
                "in={input.r1} "
                "out={output.ont_clumped} "
                "zl=9 pigz fastawrap=100000 "
                "threads={threads} "
                "-Xmx{params.mem_gb}g "
                "{params.extra} "
                "> {log} 2>&1"
            )
        elif params.input_type == 'illumina' and input.r2:
             # Ensure only the paired-end outputs are generated
            shell(f"touch {output.ont_clumped}") # Touch unused output
            shell(
                "clumpify.sh "
                "in={input.r1} in2={input.r2} "
                "out={output.r1_clumped} out2={output.r2_clumped} "
                "zl=9 pigz fastawrap=100000 "
                "threads={threads} "
                "-Xmx{params.mem_gb}g "
                "{params.extra} "
                "> {log} 2>&1"
            )
        else: # Assume single-end Illumina if r2 is missing
            shell(f"touch {output.r1_clumped} {output.r2_clumped}") # Touch unused outputs
            shell(
                "clumpify.sh "
                "in={input.r1} "
                "out={output.ont_clumped} " # Output to single-end file name
                "zl=9 pigz fastawrap=100000 "
                "threads={threads} "
                "-Xmx{params.mem_gb}g "
                "{params.extra} "
                "> {log} 2>&1"
            )


# --- NVD Pipeline Rules (Adapted from nvd.smk) ---

rule prepare_input:
    priority: 9 # Run early
    input:
        # Input files determined dynamically by the script based on config
        files = lambda wildcards: get_input_files(wildcards, config)
    output:
        fq_gz=f"{RESULTS_DIR}/01_prepare_input/{{sample}}.fq.gz"
    params:
        input_type = lambda wildcards: get_input_type(wildcards, config),
        sra_accession = lambda wildcards: get_sra_accession(wildcards, config) if get_input_type(wildcards, config) == "sra" else None,
        temp_dir = lambda wildcards: f"{wildcards.sample}.temp", # Script likely creates and removes this
        sample = "{sample}",
        out_dir = f"{RESULTS_DIR}/01_prepare_input" # Pass output dir to script
    threads: 1 # Script might download or cat, usually single-threaded
    log:
        f"{LOG_DIR}/01_prepare_input/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/01_prepare_input/{{sample}}.benchmark.txt"
    script:
        "scripts/prepare_input.py" # This script needs access to config and wildcards

rule extract_potential_human_virus_family_reads:
    priority: 9
    input:
        reads=rules.prepare_input.output.fq_gz,
        human_virus_taxlist=rules.manage_resources.output.human_virus_taxlist_local,
        stat_dbss=rules.manage_resources.output.stat_dbss_local
    output:
        fq_gz=f"{RESULTS_DIR}/02_extract_human_virus_reads/{{sample}}.fq.gz"
    threads: 4
    log:
        f"{LOG_DIR}/02_extract_human_virus_reads/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/02_extract_human_virus_reads/{{sample}}.benchmark.txt"
    shell:
        """
        seqkit fq2fa --threads 1 {input.reads} | \
        aligns_to \
            -dbss {input.stat_dbss} \
            -num_threads {threads} \
            -tax_list {input.human_virus_taxlist} \
            stdin | \
        cut -f1 | \
        seqkit grep --threads {threads} -f - {input.reads} -o {output.fq_gz} > {log} 2>&1
        """

rule run_spades:
    priority: 8
    input:
        reads=rules.extract_potential_human_virus_family_reads.output.fq_gz
    output:
        # Use touch on a marker file, as SPAdes creates a directory
        done=touch(f"{RESULTS_DIR}/03_de_novo_assemble/{{sample}}/spades.done")
    params:
        outdir=f"{RESULTS_DIR}/03_de_novo_assemble/{{sample}}/assembly",
        input_type=lambda w: get_input_type(w, config),
        # Choose SPAdes command based on input type
        spades_cmd=lambda w: "spades.py --sewage --only-assembler -s" if get_input_type(w, config) in ['ont', 'sra'] else "spades.py --sewage --12"
    threads: 8 # SPAdes can use multiple cores
    log:
        f"{LOG_DIR}/03_de_novo_assemble/{{sample}}_spades.log"
    benchmark:
        f"{LOG_DIR}/03_de_novo_assemble/{{sample}}_spades.benchmark.txt"
    shell:
        """
        mkdir -p {params.outdir}
        ({params.spades_cmd} {input.reads} -t {threads} -o {params.outdir} || echo "SPAdes finished, possibly with warnings/errors") > {log} 2>&1
        # Check if contigs file exists, even if SPAdes had non-zero exit code
        if [ ! -f "{params.outdir}/contigs.fasta" ]; then
            echo "SPAdes failed to produce contigs.fasta" >> {log}
            # Optionally create an empty contigs file to allow downstream steps to proceed gracefully
            # touch "{params.outdir}/contigs.fasta"
            # exit 1 # Or exit if contigs are essential
        fi
        """

rule move_spades_output:
    priority: 8
    input:
        # Depend on the marker file from run_spades
        spades_done=rules.run_spades.output.done,
    output:
        contigs_fa=f"{RESULTS_DIR}/03_de_novo_assemble/{{sample}}/{{sample}}.fa"
    params:
        assembly_dir=f"{RESULTS_DIR}/03_de_novo_assemble/{{sample}}/assembly", # Directory to remove,
        contigs_raw=f"{RESULTS_DIR}/03_de_novo_assemble/{{sample}}/assembly/contigs.fasta"
 
    log:
        f"{LOG_DIR}/03_de_novo_assemble/{{sample}}_move.log"
    shell:
        """
        # Check if the raw contigs file exists and has size > 0
        if [ -s {params.contigs_raw} ]; then
            echo "Moving SPAdes contigs." > {log}
            mv "{params.contigs_raw}" "{output.contigs_fa}" >> {log} 2>&1
            # Clean up the SPAdes output directory
            rm -rf "{params.assembly_dir}" >> {log} 2>&1
        else
            echo "SPAdes output contigs.fasta is missing or empty. Creating empty output file." > {log}
            touch "{output.contigs_fa}"
            # Optionally clean up assembly dir even if contigs are empty/missing
            rm -rf "{params.assembly_dir}" >> {log} 2>&1
        fi
        """

rule mask_low_complexity_contigs:
    priority: 8
    input:
        contigs=rules.move_spades_output.output.contigs_fa
    output:
        masked_contigs=f"{RESULTS_DIR}/04_mask_low_complexity_contigs/{{sample}}.fa"
    params:
        entropy=0.9,
        mem_mb=8192
    threads: 1
    log:
        f"{LOG_DIR}/04_mask_low_complexity_contigs/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/04_mask_low_complexity_contigs/{{sample}}.benchmark.txt"
    shell:
        """
        if [ -s {input.contigs} ]; then
            bbmask.sh -Xmx{params.mem_mb}m \
                in={input.contigs} \
                out={output.masked_contigs} \
                entropy={params.entropy} > {log} 2>&1
        else
            echo "Input contigs file is empty. Creating empty output file." > {log}
            touch {output.masked_contigs}
        fi
        """

rule exclude_short_contigs:
    priority: 8
    input:
        masked_contigs=rules.mask_low_complexity_contigs.output.masked_contigs
    output:
        filtered_contigs=f"{RESULTS_DIR}/05_exclude_short_contigs/{{sample}}.fa"
    params:
        min_consecutive_bases=200,
        qtrim="t",
        mem_mb=8192
    threads: 1
    log:
        f"{LOG_DIR}/05_exclude_short_contigs/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/05_exclude_short_contigs/{{sample}}.benchmark.txt"
    shell:
        """
        if [ -s {input.masked_contigs} ]; then
            reformat.sh -Xmx{params.mem_mb}m \
                in={input.masked_contigs} \
                out={output.filtered_contigs} \
                qtrim={params.qtrim} \
                minconsecutivebases={params.min_consecutive_bases} > {log} 2>&1
        else
            echo "Input masked contigs file is empty. Creating empty output file." > {log}
            touch {output.filtered_contigs}
        fi
        """

rule classify_contigs_first_pass:
    priority: 8
    input:
        fasta_file=rules.exclude_short_contigs.output.filtered_contigs,
        stat_index=rules.manage_resources.output.stat_index_local
    output:
        first_pass_stat_file=f"{RESULTS_DIR}/06_classify_contigs/{{sample}}.firstpass.txt"
    threads: 8 # Adjust based on server capacity
    log:
        f"{LOG_DIR}/06_classify_contigs/{{sample}}.classify_contigs_first_pass.log"
    benchmark:
        f"{LOG_DIR}/06_classify_contigs/{{sample}}.classify_contigs_first_pass.benchmark.txt"
    shell:
        """
        if [ -s {input.fasta_file} ]; then
            aligns_to -dbs {input.stat_index} \
                -num_threads {threads} {input.fasta_file} > {output.first_pass_stat_file} 2> {log}
        else
            echo "Input filtered contigs file is empty. Creating empty output file." > {log}
            touch {output.first_pass_stat_file}
        fi
        """

rule generate_contigs_tax_list:
    priority: 8
    input:
        first_pass_stat_file=rules.classify_contigs_first_pass.output.first_pass_stat_file
    output:
        tax_list=f"{RESULTS_DIR}/06_classify_contigs/tax_list/{{sample}}.tax_list"
    log:
        f"{LOG_DIR}/06_classify_contigs/{{sample}}.generate_tax_list.log"
    benchmark:
        f"{LOG_DIR}/06_classify_contigs/{{sample}}.generate_tax_list.benchmark.txt"
    params:
        out_dir=f"{RESULTS_DIR}/06_classify_contigs/tax_list" # Pass dir to script
    script:
        "scripts/generate_tax_list.py" # Script needs input/output paths

rule classify_contigs_second_pass:
    priority: 8
    input:
        tax_list=rules.generate_contigs_tax_list.output.tax_list,
        fasta_file=rules.exclude_short_contigs.output.filtered_contigs,
        stat_dbss=rules.manage_resources.output.stat_dbss_local
    output:
        second_pass_stat_file=f"{RESULTS_DIR}/06_classify_contigs/{{sample}}.secondpass.txt"
    threads: 8 # Adjust based on server capacity
    log:
        f"{LOG_DIR}/06_classify_contigs/{{sample}}.classify_contigs_second_pass.log"
    benchmark:
        f"{LOG_DIR}/06_classify_contigs/{{sample}}.classify_contigs_second_pass.benchmark.txt"
    shell:
        """
        if [ -s {input.fasta_file} ] && [ -s {input.tax_list} ]; then
            aligns_to \
                -tax_list {input.tax_list} \
                -dbss {input.stat_dbss} \
                -num_threads {threads} \
                {input.fasta_file} > {output.second_pass_stat_file} 2> {log}
        else
            echo "Input fasta or tax_list file is empty. Creating empty output file." > {log}
            touch {output.second_pass_stat_file}
        fi
        """

rule generate_stat_contig_report:
    priority: 8
    input:
        hits_file=rules.classify_contigs_second_pass.output.second_pass_stat_file
    output:
        report=f"{RESULTS_DIR}/07_generate_stat_contig_report/{{sample}}.report"
    params:
        cutoff_percent=0.001,
        gettax_sqlite_path=GETTAX_SQLITE_PATH # Path to the SQLite DB
    threads: 1
    log:
        f"{LOG_DIR}/07_generate_stat_contig_report/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/07_generate_stat_contig_report/{{sample}}.benchmark.txt"
    script:
        "scripts/hits_to_report.py" # Script needs input/output paths and params

rule identify_human_virus_family_contigs:
    priority: 8
    input:
        hits=rules.classify_contigs_second_pass.output.second_pass_stat_file
    output:
        filtered=f"{RESULTS_DIR}/08_identify_human_virus_family_contigs/{{sample}}.txt"
    params:
        # List of target taxa families
        taxa=["Adenoviridae", "Anelloviridae", "Arenaviridae", "Arteriviridae",
              "Astroviridae", "Bornaviridae", "Peribunyaviridae", "Caliciviridae",
              "Coronaviridae", "Filoviridae", "Flaviviridae", "Hepadnaviridae",
              "Hepeviridae", "Orthoherpesviridae", "Orthomyxoviridae", "Papillomaviridae",
              "Paramyxoviridae", "Parvoviridae", "Picobirnaviridae", "Picornaviridae",
              "Pneumoviridae", "Polyomaviridae", "Poxviridae", "Sedoreoviridae",
              "Retroviridae", "Rhabdoviridae", "Togaviridae", "Kolmioviridae"],
        stringency=0.7,
        include_children=True
    threads: 1 # Script is likely single-threaded
    log:
        f"{LOG_DIR}/08_identify_human_virus_family_contigs/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/08_identify_human_virus_family_contigs/{{sample}}.benchmark.txt"
    script:
        "scripts/extract_taxa_spots.py" # Script needs input/output paths and params

rule extract_human_virus_contigs:
    priority: 8
    input:
        human_virus_family_hits=rules.identify_human_virus_family_contigs.output.filtered,
        contigs=rules.exclude_short_contigs.output.filtered_contigs # Use the filtered contigs
    output:
        extracted_fa=f"{RESULTS_DIR}/09_extract_human_virus_contigs/{{sample}}.fa"
    threads: 1
    log:
        f"{LOG_DIR}/09_extract_human_virus_contigs/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/09_extract_human_virus_contigs/{{sample}}.benchmark.txt"
    shell:
        """
        # Check if the hits file exists and is not empty
        if [ -s {input.human_virus_family_hits} ]; then
            seqkit grep --threads {threads} -f {input.human_virus_family_hits} {input.contigs} -o {output.extracted_fa} > {log} 2>&1
        else
            echo "No human virus family hits found. Creating empty output FASTA." > {log}
            touch {output.extracted_fa}
        fi
        """

# --- GOTTCHA2 Rules (High Priority) ---
rule run_gottcha2:
    priority: 10
    input:
        # Use the prepared input reads
        fastq=rules.prepare_input.output.fq_gz,
        # Gottcha needs the actual DB files, not just the marker
        db_species=rules.manage_resources.output.gottcha_species_local,
    output:
        # Define outputs within the sample-specific subdirectory
        tsv_full=f"{RESULTS_DIR}/10_gottcha2/{{sample}}/{{sample}}.full.tsv",
        sam_strain=f"{RESULTS_DIR}/10_gottcha2/{{sample}}/{{sample}}.gottcha_strain.sam",
        log_out=f"{RESULTS_DIR}/10_gottcha2/{{sample}}/{{sample}}.gottcha_strain.log" # Combined log
    params:
        outdir=f"{RESULTS_DIR}/10_gottcha2/{{sample}}",
        prefix="{sample}.gottcha",
        db_species_path=lambda w, input: os.path.splitext(input.db_species)[0], # Path without .mmi
    threads: 12 # High priority, allocate significant threads
    log:
        # Use the output log file directly
        f"{LOG_DIR}/10_gottcha2/{{sample}}/{{sample}}.gottcha.log"
    benchmark:
        f"{LOG_DIR}/10_gottcha2/{{sample}}.benchmark.txt"
    shell:
        """
        mkdir -p {params.outdir}
        echo "Running GOTTCHA2 Species Level for {wildcards.sample}" > {log}
        gottcha2.py \
            -i {input.fastq} \
            -d {params.db_species_path} \
            -t {threads} \
            -o {params.outdir} \
            -l strain \
            --noCutoff \
            >> {log} 2>&1
        """

# --- GOTTCHA2 Validation Rules ---
rule gottcha2_extract_taxon:
    priority: 10
    input:
        sam=rules.run_gottcha2.output.sam_strain,
        db_species=rules.manage_resources.output.gottcha_species_local # Use the actual DB file
    output:
        # Note: gottcha2.py -e creates output based on prefix, need to rename
        fastq_tmp=temp(f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.gottcha_strain.extract.fastq"),
        fastq=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.extract.fastq",
        log=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.extract.log"
    params:
        taxid="{taxid}",
        outdir=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}",
        prefix="{sample}_{taxid}", # Prefix used by gottcha2.py
        level="strain",
        db_species_path=lambda w, input: os.path.splitext(input.db_species)[0] # Path without .fna
    threads: 4
    log:
         f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.extract.log" # Use output log
    benchmark:
        f"{LOG_DIR}/11_gottcha2_verify/{{sample}}_{{taxid}}.extract.benchmark.txt"
    shell:
        """
        mkdir -p {params.outdir}
        echo "Extracting reads for taxid {params.taxid} from {input.sam}" > {log}
        gottcha2.py -e {params.taxid} --noCutoff \
            -s {input.sam} \
            -d {params.db_species_path} \
            -t {threads} \
            -o {params.outdir}/{params.prefix} \
            -l {params.level} >> {log} 2>&1

        # Check and rename the output file created by gottcha2.py
        GOTTCHA_OUT_FQ="{params.outdir}/{params.prefix}.{params.level}.extract.fastq"
        if [ -f "$GOTTCHA_OUT_FQ" ]; then
            mv "$GOTTCHA_OUT_FQ" "{output.fastq}" >> {log} 2>&1
            echo "Extraction complete." >> {log}
        else
            echo "Error: Expected output file $GOTTCHA_OUT_FQ not found. Creating empty file." >> {log}
            touch "{output.fastq}"
            # exit 1 # Optionally exit
        fi
        """

rule reformat_extracted_fastq:
    priority: 10
    input:
        fastq=rules.gottcha2_extract_taxon.output.fastq
    output:
        fasta=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.extract.fasta",
        log=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.reformat.log"
    params:
        max_reads=1000
    threads: 1
    log:
        f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.reformat.log" # Use output log
    benchmark:
        f"{LOG_DIR}/11_gottcha2_verify/{{sample}}_{{taxid}}.reformat.benchmark.txt"
    shell:
        """
        if [ -s {input.fastq} ]; then
            reformat.sh \
                in={input.fastq} \
                out={output.fasta} \
                nullifybrokenquality=t \
                samplereadstarget={params.max_reads} > {log} 2>&1
        else
            echo "Input FASTQ file is empty. Creating empty FASTA file." > {log}
            touch {output.fasta}
        fi
        """

rule gottcha_blast_verify:
    priority: 10
    input:
        fasta=rules.reformat_extracted_fastq.output.fasta,
        # Ensure resources are ready, blast_verify.py might download taxdump/taxdb
        ready=rules.manage_resources.output[0]
    output:
        class_tsv=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.classification.tsv",
        hits_tsv=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.hits.tsv",
        log=f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.verify.log"
    params:
        taxid="{taxid}",
        # Assuming blast_verify.py uses a DB prefix relative to RESOURCES_DIR
        blast_db=os.path.join(RESOURCES_DIR, os.path.basename(BLAST_VERIFY_DB)),
        script=config["scripts"]["blast_verify"], # Path to blast_verify.py
        word_size=config["blast_verify"].get("word_size", 64),
        max_target_seqs=config["blast_verify"].get("max_target_seqs", None) # Optional
    threads: 8 # BLAST verification can use threads
    log:
        f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.verify.log" # Use output log
    benchmark:
        f"{LOG_DIR}/11_gottcha2_verify/{{sample}}_{{taxid}}.verify.benchmark.txt"
    shell:
        """
        if [ -s {input.fasta} ]; then
            python {params.script} \
                --query {input.fasta} \
                --target-taxid {params.taxid} \
                --db {params.blast_db} \
                --threads {threads} \
                --word_size {params.word_size} \
                $(if [ ! -z "{params.max_target_seqs}" ]; then echo '--max_target_seqs {params.max_target_seqs}'; fi) \
                --hits-out {output.hits_tsv} \
                --class-out {output.class_tsv} > {log} 2>&1
        else
            echo "Input FASTA file is empty. Skipping verification, creating empty output files." > {log}
            touch {output.class_tsv} {output.hits_tsv}
        fi
        """

# --- GOTTCHA2 Extract 20 Reads Per Taxon Rules ---
rule gottcha2_extract_all_taxa:
    priority: 10
    input:
        sam=rules.run_gottcha2.output.sam_strain,
        db_species=rules.manage_resources.output.gottcha_species_local
    output:
        # Expecting gottcha2.py -ef to create {prefix}.{level}.extract_filtered.fastq
        fastq_tmp=temp(f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.gottcha_strain.extract_filtered.fastq"),
        fastq=f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.extract_20.fastq",
        log=f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.extract_20.log"
    params:
        outdir=f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}",
        prefix="{sample}", # Prefix used by gottcha2.py
        level="strain",
        db_species_path=lambda w, input: os.path.splitext(input.db_species)[0] # Path without .fna
    threads: 4
    log:
        f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.extract_20.log" # Use output log
    benchmark:
        f"{LOG_DIR}/12_gottcha2_extract_all/{{sample}}.extract.benchmark.txt"
    shell:
        """
        mkdir -p {params.outdir}
        echo "Extracting top 20 reads per taxon from {input.sam}" > {log}
        gottcha2.py -ef --noCutoff \
            -s {input.sam} \
            -d {params.db_species_path} \
            -t {threads} \
            -o {params.outdir}/{params.prefix} \
            -l {params.level} >> {log} 2>&1

        # Check and rename the output file created by gottcha2.py -ef
        GOTTCHA_OUT_FQ="{params.outdir}/{params.prefix}.{params.level}.extract_filtered.fastq"
        if [ -f "$GOTTCHA_OUT_FQ" ]; then
            mv "$GOTTCHA_OUT_FQ" "{output.fastq}" >> {log} 2>&1
            echo "Extraction complete." >> {log}
        else
            echo "Error: Expected output file $GOTTCHA_OUT_FQ not found. Creating empty file." >> {log}
            touch "{output.fastq}"
            # exit 1 # Optionally exit
        fi
        """

rule reformat_all_extracted_fastq:
    priority: 10
    input:
        fastq=rules.gottcha2_extract_all_taxa.output.fastq
    output:
        fasta=f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.extract_20.fasta",
        log=f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.reformat_20.log"
    threads: 1
    log:
        f"{RESULTS_DIR}/12_gottcha2_extract_all/{{sample}}/{{sample}}.reformat_20.log" # Use output log
    benchmark:
        f"{LOG_DIR}/12_gottcha2_extract_all/{{sample}}.reformat.benchmark.txt"
    shell:
        """
        if [ -s {input.fastq} ]; then
            reformat.sh \
                in={input.fastq} \
                out={output.fasta} \
                nullifybrokenquality=t > {log} 2>&1
        else
            echo "Input extracted FASTQ file is empty. Creating empty FASTA file." > {log}
            touch {output.fasta}
        fi
        """

# --- NVD BLAST and Post-processing Rules (Renumbered) ---

rule megablast:
    priority: 7
    input:
        contigs=rules.extract_human_virus_contigs.output.extracted_fa,
        ready=rules.manage_resources.output[0] # Ensure resources (like BLASTDB) are set up
    output:
        blast_results=f"{RESULTS_DIR}/13_megablast/{{sample}}.txt"
    params:
        db=os.path.basename(BLAST_DB_PATH), # Use basename, BLASTDB env var points to dir
        outfmt="6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids",
        max_target_seqs=5,
        task="megablast",
        blastdb=RESOURCES_DIR # Point to the local resources directory
    threads: 4
    log:
        f"{LOG_DIR}/13_megablast/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/13_megablast/{{sample}}.benchmark.txt"
    shell:
        """
        export BLASTDB={params.blastdb}
        if [ -s {input.contigs} ]; then
            blastn -task {params.task} \
                -db {params.db} \
                -query {input.contigs} \
                -num_threads {threads} \
                -outfmt "{params.outfmt}" \
                -max_target_seqs {params.max_target_seqs} \
                -out {output.blast_results} > {log} 2>&1
        else
            echo "Input contigs file is empty. Creating empty BLAST output." > {log}
            touch {output.blast_results}
        fi
        """

rule annotate_megablast_results:
    priority: 7
    input:
        blast_results=rules.megablast.output.blast_results
    output:
        annotated_results=f"{RESULTS_DIR}/14_annotate_megablast_results/{{sample}}.txt"
    params:
        sample="{sample}",
        gettax_sqlite_path=GETTAX_SQLITE_PATH,
        task="megablast"
    threads: 1
    log:
        f"{LOG_DIR}/14_annotate_megablast_results/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/14_annotate_megablast_results/{{sample}}.benchmark.txt"
    script:
        "scripts/annotate_blast_results.py" # Script needs input/output/params

rule filter_non_virus_megablast_nodes:
    priority: 7
    input:
        annotated_results=rules.annotate_megablast_results.output.annotated_results
    output:
        filtered_results=f"{RESULTS_DIR}/15_filter_non_virus_megablast_nodes/{{sample}}.txt"
    threads: 1
    log:
        f"{LOG_DIR}/15_filter_non_virus_megablast_nodes/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/15_filter_non_virus_megablast_nodes/{{sample}}.benchmark.txt"
    script:
        "scripts/filter_non_virus_megablast_nodes.py" # Script needs input/output

rule remove_megablast_mapped_contigs:
    priority: 7
    input:
        megablast_results=rules.megablast.output.blast_results, # Use original BLAST output for IDs
        contigs_fasta=rules.extract_human_virus_contigs.output.extracted_fa
    output:
        classified_contigs=f"{RESULTS_DIR}/16_remove_megablast_mapped_contigs/{{sample}}.classified.txt",
        pruned_contigs=f"{RESULTS_DIR}/16_remove_megablast_mapped_contigs/{{sample}}.pruned.fa"
    threads: 1
    log:
        f"{LOG_DIR}/16_remove_megablast_mapped_contigs/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/16_remove_megablast_mapped_contigs/{{sample}}.benchmark.txt"
    script:
        "scripts/remove_megablast_mapped_contigs.py" # Script needs input/output

rule blastn_classify:
    priority: 6 # Lower priority than megablast
    input:
        contigs=rules.remove_megablast_mapped_contigs.output.pruned_contigs,
        ready=rules.manage_resources.output[0]
    output:
        blast_results=f"{RESULTS_DIR}/17_blastn_classify/{{sample}}.txt"
    params:
        db=os.path.basename(BLAST_DB_PATH),
        outfmt="6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids",
        max_target_seqs=5,
        task="blastn",
        blastdb=RESOURCES_DIR
    threads: 4
    log:
        f"{LOG_DIR}/17_blastn_classify/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/17_blastn_classify/{{sample}}.benchmark.txt"
    shell:
        """
        export BLASTDB={params.blastdb}
        if [ -s {input.contigs} ]; then
            blastn -task {params.task} \
                -db {params.db} \
                -query {input.contigs} \
                -num_threads {threads} \
                -outfmt "{params.outfmt}" \
                -max_target_seqs {params.max_target_seqs} \
                -out {output.blast_results} > {log} 2>&1
        else
            echo "Input pruned contigs file is empty. Creating empty BLAST output." > {log}
            touch {output.blast_results}
        fi
        """

rule annotate_blastn_results:
    priority: 6
    input:
        blast_results=rules.blastn_classify.output.blast_results
    output:
        annotated_results=f"{RESULTS_DIR}/18_annotate_blastn_results/{{sample}}.txt"
    params:
        sample="{sample}",
        gettax_sqlite_path=GETTAX_SQLITE_PATH,
        task="blastn"
    threads: 1
    log:
        f"{LOG_DIR}/18_annotate_blastn_results/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/18_annotate_blastn_results/{{sample}}.benchmark.txt"
    script:
        "scripts/annotate_blast_results.py" # Script needs input/output/params

rule filter_non_virus_blastn_nodes:
    priority: 6
    input:
        annotated_results=rules.annotate_blastn_results.output.annotated_results
    output:
        filtered_results=f"{RESULTS_DIR}/19_filter_non_virus_blastn_nodes/{{sample}}.txt"
    threads: 1
    log:
        f"{LOG_DIR}/19_filter_non_virus_blastn_nodes/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/19_filter_non_virus_blastn_nodes/{{sample}}.benchmark.txt"
    script:
        "scripts/filter_non_virus_megablast_nodes.py" # Re-use script

rule merge_annotated_blast_results:
    priority: 5
    input:
        megablast=rules.filter_non_virus_megablast_nodes.output.filtered_results,
        blastn=rules.filter_non_virus_blastn_nodes.output.filtered_results
    output:
        merged=f"{RESULTS_DIR}/20_merge_annotated_blast_results/{{sample}}.txt"
    threads: 1
    log:
        f"{LOG_DIR}/20_merge_annotated_blast_results/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/20_merge_annotated_blast_results/{{sample}}.benchmark.txt"
    shell:
        """
        # Concatenate files, ensuring output exists even if inputs are empty
        cat {input.megablast} {input.blastn} > {output.merged} 2> {log}
        # Verify output was created
        if [ ! -f {output.merged} ]; then
            echo "Failed to create merged BLAST results file." >> {log}
            touch {output.merged} # Create empty file to satisfy downstream rules
        fi
        """

rule extract_unclassified_contigs:
    priority: 5
    input:
        contigs=rules.extract_human_virus_contigs.output.extracted_fa,
        megablast_results=rules.megablast.output.blast_results, # Original megablast output
        blastn_results=rules.blastn_classify.output.blast_results # Original blastn output
    output:
        unclassified=f"{RESULTS_DIR}/21_extract_unclassified_contigs/{{sample}}.unclassified.fa"
    threads: 1
    log:
        f"{LOG_DIR}/21_extract_unclassified_contigs/{{sample}}.log"
    benchmark:
        f"{LOG_DIR}/21_extract_unclassified_contigs/{{sample}}.benchmark.txt"
    script:
        "scripts/extract_unclassified_contigs.py" # Script needs input/output

rule map_reads_to_contigs:
    priority: 4
    input:
        contigs=rules.extract_human_virus_contigs.output.extracted_fa,
        reads=rules.extract_potential_human_virus_family_reads.output.fq_gz # Use the filtered reads
    output:
        bam=f"{RESULTS_DIR}/22_map_reads_to_contigs/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/22_map_reads_to_contigs/{{sample}}.bam.bai"
    params:
        mapper=lambda w: "map-ont" if get_input_type(w, config) in ['ont', 'sra'] else "sr"
    threads: 4
    log:
        f"{LOG_DIR}/22_map_reads_to_contigs/{{sample}}.map.log"
    benchmark:
        f"{LOG_DIR}/22_map_reads_to_contigs/{{sample}}.map.benchmark.txt"
    shell:
        """
        # Check if contigs file is empty before running minimap2
        if [ ! -s {input.contigs} ]; then
            echo "Contigs file is empty, skipping read mapping." > {log}
            # Create empty BAM and BAI files to satisfy downstream rules
            samtools view -b -o {output.bam} /dev/null
            touch {output.bai}
        else
            (minimap2 -ax {params.mapper} -t {threads} {input.contigs} {input.reads} | \
            samtools view -b -F 4 | \
            samtools sort -@ {threads} -o {output.bam} && \
            samtools index {output.bam}) > {log} 2>&1
        fi
        """

rule count_mapped_reads:
    priority: 4
    input:
        bam=rules.map_reads_to_contigs.output.bam,
        bai=rules.map_reads_to_contigs.output.bai # Ensure index is built first
    output:
        counts=f"{RESULTS_DIR}/23_count_mapped_reads/{{sample}}_mapped_counts.txt"
    log:
        f"{LOG_DIR}/23_count_mapped_reads/{{sample}}.count.log"
    benchmark:
        f"{LOG_DIR}/23_count_mapped_reads/{{sample}}.count.benchmark.txt"
    shell:
        """
        # Check if BAM file exists and is not empty before counting
        if [ -s {input.bam} ]; then
            # samtools idxstats requires a non-empty indexed BAM
            samtools idxstats {input.bam} | awk '{{print $1 "\\t" $3}}' > {output.counts} 2> {log}
        else
            echo "Input BAM file is empty or missing. Creating empty counts file." > {log}
            touch {output.counts}
        fi
        """

# --- Aggregation and Archiving ---
rule aggregate_validation_results:
    """Aggregates ambiguous/unambiguous counts from GOTTCHA2 verification TSVs."""
    priority: 3
    input:
        # Use expand with sample wildcard fixed
        expand(f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.classification.tsv", taxid=VALIDATION_TAXA, sample=SAMPLES)
    output:
        summary=f"{RESULTS_DIR}/25_labkey/{{sample}}/gottcha_validation_summary.tsv"
    params:
        sample="{sample}"
    log:
        f"{LOG_DIR}/25_labkey/{{sample}}.aggregate_validation.log"
    run:
        all_results = []
        print(f"Aggregating validation results for sample {params.sample}", file=open(log.log, 'w'))
        print(f"Input files: {input}", file=open(log.log, 'a'))
        for f in input:
            if os.path.exists(f) and os.path.getsize(f) > 0:
                try:
                    df = pd.read_csv(f, sep='\t')
                    # Extract taxid from filename path more robustly
                    taxid = os.path.basename(os.path.dirname(f))
                    unambig_count = df[df['status'] == 'unambiguous'].shape[0]
                    ambig_count = df[df['status'] == 'ambiguous'].shape[0]
                    all_results.append({
                        "sample_id": params.sample,
                        "validation_taxid": taxid,
                        "unambiguous_reads": unambig_count,
                        "ambiguous_reads": ambig_count
                    })
                    print(f"Processed {f}: Unambig={unambig_count}, Ambig={ambig_count}", file=open(log.log, 'a'))
                except Exception as e:
                    print(f"Warning: Could not process file {f}: {e}", file=open(log.log, 'a'))
            else:
                 print(f"Warning: Skipping missing or empty file {f}", file=open(log.log, 'a'))

        if all_results:
            summary_df = pd.DataFrame(all_results)
            summary_df.to_csv(output.summary, sep='\t', index=False)
            print(f"Aggregation summary saved to {output.summary}", file=open(log.log, 'a'))
        else:
            # Create empty file with header if no results
            print(f"No valid validation results found for sample {params.sample}. Creating empty summary file.", file=open(log.log, 'a'))
            with open(output.summary, 'w') as outfile:
                outfile.write("sample_id\tvalidation_taxid\tunambiguous_reads\tambiguous_reads\n")


rule create_tarball:
    priority: 2
    input:
        # NVD outputs
        virus_contigs=rules.extract_human_virus_contigs.output.extracted_fa,
        unclassified_contigs=rules.extract_unclassified_contigs.output.unclassified,
        mapped_reads=rules.map_reads_to_contigs.output.bam,
        mapped_reads_index=rules.map_reads_to_contigs.output.bai,
        # GOTTCHA2 outputs
        gottcha_full_tsv=rules.run_gottcha2.output.tsv_full,
        # GOTTCHA2 validation FASTAs (use expand with fixed sample)
        validation_fastas=expand(f"{RESULTS_DIR}/11_gottcha2_verify/{{sample}}/{{taxid}}/{{sample}}_{{taxid}}.extract.fasta", taxid=VALIDATION_TAXA, sample=SAMPLES),
        # GOTTCHA2 20-reads FASTA
        extract_20_fasta=rules.reformat_all_extracted_fastq.output.fasta,
    output:
        tarball=f"{RESULTS_DIR}/24_create_tarball/{{sample}}.tar.zst"
    params:
        # Define files relative to RESULTS_DIR for archiving
        files_to_archive=lambda w, input: [
            os.path.relpath(f, RESULTS_DIR) for f in [
                input.virus_contigs, input.unclassified_contigs,
                input.mapped_reads, input.mapped_reads_index,
                input.gottcha_full_tsv, input.extract_20_fasta
            ] + input.validation_fastas if os.path.exists(f) and os.path.getsize(f) > 0 # Only include existing, non-empty files
        ]
    threads: 2
    log:
        f"{LOG_DIR}/24_create_tarball/{{sample}}.log"
    shell:
        """
        echo "Creating archive for {wildcards.sample}" > {log}
        # Change to results directory to archive with relative paths
        cd {RESULTS_DIR}

        # Check if there are any files to archive
        files=({params.files_to_archive})
        if [ ${{#files[@]}} -eq 0 ]; then
            echo "No valid files found to archive for sample {wildcards.sample}. Creating empty archive." >> ../{log}
            # Create empty tarball outside results dir
            touch ../{output.tarball}
        else
            echo "Archiving files:" >> ../{log}
            printf '%s\n' "${{files[@]}}" >> ../{log}
            # Create archive outside results dir, using relative paths from within results dir
            tar -I pzstd --warning=no-file-changed -cf ../{output.tarball} "${{files[@]}}" >> ../{log} 2>&1
            echo "Archive creation complete." >> ../{log}
        fi
        """


# --- LabKey Upload Rules (Low Priority) ---
# Using the script-based approach from nvd.smk, adapted for new uploads

rule labkey_upload_nvd_blast:
    priority: 0
    input:
        blast_results=rules.merge_annotated_blast_results.output.merged,
        mapped_counts=rules.count_mapped_reads.output.counts,
        # Need original prepared input for context if script uses it
        prepared_input=rules.prepare_input.output.fq_gz
    output:
        token=touch(f"{RESULTS_DIR}/25_labkey/{{sample}}/nvd_blast_upload.tkn"),
        # Also output CSV as fallback or primary output if LabKey isn't configured
        csv_file=f"{FINAL_OUTPUT_PATH}/{{sample}}.nvd_blast_results.csv"
    params:
        sample="{sample}",
        experiment=EXPERIMENT,
        blast_db_name=os.path.basename(BLAST_DB_PATH),
        snakemake_run_id=config['run_id'], # Use run_id from config
        labkey_server=LABKEY_SERVER,
        project_name=LABKEY_PROJECT_NAME,
        api_key=LABKEY_API_KEY,
        stat_db_version=os.path.basename(STAT_DBSS),
        out_dir=FINAL_OUTPUT_PATH, # Directory for CSV output
        upload_type="blast" # Tell script what to upload
    threads: 1
    log:
        f"{LOG_DIR}/25_labkey/{{sample}}.nvd_blast_upload.log"
    benchmark:
        f"{LOG_DIR}/25_labkey/{{sample}}.nvd_blast_upload.benchmark.txt"
    script:
        "scripts/labkey_upload.py" # Assumes this script handles both LabKey and CSV output

rule labkey_upload_nvd_fasta:
    priority: 0
    input:
        fasta=rules.extract_human_virus_contigs.output.extracted_fa
    output:
        token=touch(f"{RESULTS_DIR}/25_labkey/{{sample}}/nvd_fasta_upload.tkn")
    params:
        sample="{sample}",
        experiment=EXPERIMENT,
        labkey_server=LABKEY_SERVER,
        project_name=LABKEY_PROJECT_NAME,
        api_key=LABKEY_API_KEY,
        snakemake_run_id=config['run_id'],
        singleline_fasta_path=temp(f"{RESULTS_DIR}/25_labkey/{{sample}}/nvd_contigs.tmp.fasta"),
        upload_type="nvd_fasta" # Specify upload type
    threads: 1
    log:
        f"{LOG_DIR}/25_labkey/{{sample}}.nvd_fasta_upload.log"
    benchmark:
        f"{LOG_DIR}/25_labkey/{{sample}}.nvd_fasta_upload.benchmark.txt"
    script:
        "scripts/labkey_upload_fasta.py" # Assumes this script handles the upload

rule labkey_upload_gottcha:
    priority: 0
    input:
        # Main GOTTCHA results
        tsv_full=rules.run_gottcha2.output.tsv_full,
        # Aggregated validation summary
        validation_summary=rules.aggregate_validation_results.output.summary
    output:
        token=touch(f"{RESULTS_DIR}/25_labkey/{{sample}}/gottcha_upload.tkn")
    params:
        sample="{sample}",
        experiment=EXPERIMENT,
        labkey_server=LABKEY_SERVER,
        project_name=LABKEY_PROJECT_NAME,
        api_key=LABKEY_API_KEY,
        gottcha2_db_version=os.path.basename(GOTTCHA2_DB_SPECIES), # Use DB name as version
        upload_type="gottcha" # Specify upload type
    threads: 1
    log:
        f"{LOG_DIR}/25_labkey/{{sample}}.gottcha_upload.log"
    benchmark:
        f"{LOG_DIR}/25_labkey/{{sample}}.gottcha_upload.benchmark.txt"
    script:
        "scripts/labkey_upload_gottcha.py" # Needs a dedicated script to handle TSV + summary

rule labkey_upload_validation_fasta:
    priority: 0
    input:
        fasta=rules.reformat_extracted_fastq.output.fasta
    output:
        token=touch(f"{RESULTS_DIR}/25_labkey/{{sample}}/validation_fasta_upload_{{taxid}}.tkn")
    params:
        sample="{sample}",
        taxid="{taxid}",
        experiment=EXPERIMENT,
        labkey_server=LABKEY_SERVER,
        project_name=LABKEY_PROJECT_NAME,
        api_key=LABKEY_API_KEY,
        snakemake_run_id=config['run_id'],
        upload_type="validation_fasta" # Specify upload type
    threads: 1
    log:
        f"{LOG_DIR}/25_labkey/{{sample}}.validation_fasta_upload_{{taxid}}.log"
    benchmark:
        f"{LOG_DIR}/25_labkey/{{sample}}.validation_fasta_upload_{{taxid}}.benchmark.txt"
    script:
        "scripts/labkey_upload_fasta.py" # Reuse or adapt fasta upload script

rule labkey_upload_extract_20_fasta:
    priority: 0
    input:
        fasta=rules.reformat_all_extracted_fastq.output.fasta
    output:
        token=touch(f"{RESULTS_DIR}/25_labkey/{{sample}}/extract_20_fasta_upload.tkn")
    params:
        sample="{sample}",
        experiment=EXPERIMENT,
        labkey_server=LABKEY_SERVER,
        project_name=LABKEY_PROJECT_NAME,
        api_key=LABKEY_API_KEY,
        snakemake_run_id=config['run_id'],
        upload_type="extract_20_fasta" # Specify upload type
    threads: 1
    log:
        f"{LOG_DIR}/25_labkey/{{sample}}.extract_20_fasta_upload.log"
    benchmark:
        f"{LOG_DIR}/25_labkey/{{sample}}.extract_20_fasta_upload.benchmark.txt"
    script:
        "scripts/labkey_upload_fasta.py" # Reuse or adapt fasta upload script

rule upload_archive_to_labkey:
    priority: 0
    input:
        tarball=rules.create_tarball.output.tarball
    output:
        token=touch(f"{RESULTS_DIR}/25_labkey/{{sample}}/archive_upload.tkn")
    params:
        sample="{sample}",
        experiment=EXPERIMENT,
        labkey_server=LABKEY_SERVER,
        project_name=LABKEY_PROJECT_NAME,
        webdav_url=WEBDAV_URL, # WebDAV URL for file uploads
        username=LABKEY_USERNAME,
        password=LABKEY_PASSWORD,
        snakemake_run_id=config['run_id'],
        final_output_path=FINAL_OUTPUT_PATH, # Base path on server
        upload_type="archive" # Specify upload type
    threads: 1
    log:
        f"{LOG_DIR}/25_labkey/{{sample}}.archive_upload.log"
    benchmark:
        f"{LOG_DIR}/25_labkey/{{sample}}.archive_upload.benchmark.txt"
    script:
        "scripts/upload_files_to_labkey.py" # Assumes this script handles archive uploads via WebDAV

# --- Cleanup Rule ---
rule clean:
    """Remove all generated files and directories."""
    shell:
        """
        echo "Cleaning up results, logs, and temporary files..."
        rm -rf {RESULTS_DIR} {LOG_DIR} resources *.temp* decompressed.txt resources.tar.zst
        # Add specific temp files if needed: rm -f fasterq.tmp* *.temp
        echo "Cleanup complete."
        """
