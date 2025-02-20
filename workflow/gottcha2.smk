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

# Define output paths using Python string formatting
FULL_OUT    = f"out/{sample}_gottcha2.out.full.tsv"
LINEAGE_OUT = f"out/{sample}_gottcha2.out.lineage.tsv"
TSV_OUT     = f"out/{sample}_gottcha2.out.tsv"

rule all:
    input:
        os.path.join(gottcha2_output, sample, f"{sample}_gottcha2.out.full.tsv"),
        os.path.join(gottcha2_output, sample, f"{sample}_gottcha2.out.lineage.tsv"),
        os.path.join(gottcha2_output, sample, f"{sample}_gottcha2.out.tsv")

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
rule copy_input:
    input:
        r1_fastq,
        r2_fastq
    output:
        temp(f"{sample}_R1.fastq.gz"),
        temp(f"{sample}_R2.fastq.gz"),
    shell:
        """
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """

##########################################
# GOTTCHA2 ANALYSIS
##########################################
rule gottcha2:
    input:
        r1 = f"{sample}_R1.fastq.gz",
        r2 = f"{sample}_R2.fastq.gz",
        db_mmi   = "ref/gottcha_db.species.fna.mmi",
        db_stats = "ref/gottcha_db.species.fna.stats",
        db_tax   = "ref/gottcha_db.species.fna.tax.tsv"
    output:
        full    = FULL_OUT,
        lineage = LINEAGE_OUT,
        tsv     = TSV_OUT,
    params:
        output_path = f"out/{sample}",
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
        full    = os.path.join(gottcha2_output, sample, f"{sample}_gottcha2.out.full.tsv"),
        lineage = os.path.join(gottcha2_output, sample, f"{sample}_gottcha2.out.lineage.tsv"),
        tsv     = os.path.join(gottcha2_output, sample, f"{sample}_gottcha2.out.tsv"),
    shell:
        """
        mkdir -p $(dirname {output.full})
        cp {input.full} {output.full}
        cp {input.lineage} {output.lineage}
        cp {input.tsv} {output.tsv}
        """
