import os
import glob
import re

# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Retrieve configuration parameters
# From command line
sample = config["sample_name"]
r1_fastq = config["r1_fastq"]
r2_fastq = config["r2_fastq"]

# From parameters file
clumped_output = config['global']["clumped_dir"]
gottcha2_output = config['global']["gottcha2_dir"]

##########################################
# CONFIGURATION FROM THE COMMAND LINE
##########################################
# Base folder where final outputs are copied (includes experiment number)
BASE_FOLDER = gottcha2_output

# Comma-separated lists of full paths to the original FASTQ files.
R1_FILES = [f.strip() for f in config["r1_files"].split(",")]
R2_FILES = [f.strip() for f in config["r2_files"].split(",")]

# Debug flag: when true, intermediate files are preserved.
DEBUG = config.get("debug_intermediates", False)

# Final output directory: BASE_FOLDER/out
OUT_DIR = os.path.join(BASE_FOLDER, "out")

##########################################
# BUILD SAMPLE-TO-FILE MAPPINGS
##########################################
# We extract the sample name from the basename by assuming that everything before "_R1" or "_R2" is the sample name.
sample_to_r1 = {}
sample_to_r2 = {}
for f in R1_FILES:
    base = os.path.basename(f)
    m = re.match(r"(.*)_R1.*\.fastq\.gz", base)
    if m:
        sample = m.group(1)
        sample_to_r1[sample] = f

for f in R2_FILES:
    base = os.path.basename(f)
    m = re.match(r"(.*)_R2.*\.fastq\.gz", base)
    if m:
        sample = m.group(1)
        sample_to_r2[sample] = f

# Only use samples that have both R1 and R2.
samples = sorted(set(sample_to_r1.keys()) & set(sample_to_r2.keys()))

##########################################
# INTERMEDIATE FILE NAMING (optionally temporary)
##########################################
if DEBUG:
    FULL_OUT    = "out/{sample}_gottcha2.out.full.tsv"
    LINEAGE_OUT = "out/{sample}_gottcha2.out.lineage.tsv"
    TSV_OUT     = "out/{sample}_gottcha2.out.tsv"
else:
    FULL_OUT    = temp("out/{sample}_gottcha2.out.full.tsv")
    LINEAGE_OUT = temp("out/{sample}_gottcha2.out.lineage.tsv")
    TSV_OUT     = temp("out/{sample}_gottcha2.out.tsv")

##########################################
# FINAL TARGET: Copy final results to BASE_FOLDER/out/<sample>/
##########################################
rule all:
    input:
        expand(os.path.join(OUT_DIR, "{sample}", "{sample}_gottcha2.out.full.tsv"), sample=samples),
        expand(os.path.join(OUT_DIR, "{sample}", "{sample}_gottcha2.out.lineage.tsv"), sample=samples),
        expand(os.path.join(OUT_DIR, "{sample}", "{sample}_gottcha2.out.tsv"), sample=samples),

##########################################
# REFERENCE FILES (ref/)
##########################################
rule copy_gottcha_db:
    output:
        mmi   = "ref/gottcha_db.species.fna.mmi",
        stats = "ref/gottcha_db.species.fna.stats",
        tax   = "ref/gottcha_db.species.fna.tax.tsv"
    input:
        mmi   = "/staging/groups/oconnor_group/gottcha2/gottcha_db.species.fna.mmi",
        stats = "/staging/groups/oconnor_group/gottcha2/gottcha_db.species.fna.stats",
        tax   = "/staging/groups/oconnor_group/gottcha2/gottcha_db.species.fna.tax.tsv"
    shell:
        """
        mkdir -p ref
        cp {input.mmi} {output.mmi}
        cp {input.stats} {output.stats}
        cp {input.tax} {output.tax}
        """

##########################################
# STAGE INPUT FILES: Copy the FASTQ files into the local "in/" folder.
##########################################
rule copy_input_R1:
    output:
        "in/{sample}_R1.fastq.gz"
    input:
        lambda wc: sample_to_r1[wc.sample]
    shell:
        """
        mkdir -p in
        cp {input} {output}
        """

rule copy_input_R2:
    output:
        "in/{sample}_R2.fastq.gz"
    input:
        lambda wc: sample_to_r2[wc.sample]
    shell:
        """
        mkdir -p in
        cp {input} {output}
        """

##########################################
# GOTTCHA2 ANALYSIS
##########################################
rule gottcha2:
    input:
        r1 = "in/{sample}_R1.fastq.gz",
        r2 = "in/{sample}_R2.fastq.gz",
        db_mmi   = "ref/gottcha_db.species.fna.mmi",
        db_stats = "ref/gottcha_db.species.fna.stats",
        db_tax   = "ref/gottcha_db.species.fna.tax.tsv"
    output:
        full    = FULL_OUT,
        lineage = LINEAGE_OUT,
        tsv     = TSV_OUT,
        sam     = SAM_OUT
    params:
        output_path = "out/{sample}_gottcha2.out"
    threads: 16
    shell:
        """
        gottcha2.py \
          --database "ref/gottcha_db.species.fna" \
          --prefix {params.output_path} \
          -t {threads} \
          -i {input.r1} {input.r2}
        """

##########################################
# COPY FINAL RESULTS TO THE OUTPUT FOLDER
##########################################
rule copy_results:
    input:
        full    = FULL_OUT,
        lineage = LINEAGE_OUT,
        tsv     = TSV_OUT,
    output:
        full    = os.path.join(OUT_DIR, "{sample}", "{sample}_gottcha2.out.full.tsv"),
        lineage = os.path.join(OUT_DIR, "{sample}", "{sample}_gottcha2.out.lineage.tsv"),
        tsv     = os.path.join(OUT_DIR, "{sample}", "{sample}_gottcha2.out.tsv"),
    shell:
        """
        mkdir -p $(dirname {output.full})
        cp {input.full} {output.full}
        cp {input.lineage} {output.lineage}
        cp {input.tsv} {output.tsv}
        """
