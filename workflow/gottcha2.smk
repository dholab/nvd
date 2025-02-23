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

# From parameters file - updated paths
gottcha2_out_dir = config['global']["gottcha2_out_dir"]

# Define output paths using Python string formatting
FULL_OUT    = f"{sample}.full.tsv"
LINEAGE_OUT = f"{sample}.lineage.tsv"
TSV_OUT     = f"{sample}.tsv"

rule all:
    input:
        os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.full.tsv"),
        os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.lineage.tsv"),
        os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.summary.tsv"),
        os.path.join(gottcha2_out_dir, sample, "gottcha2_full.upload.done")

##########################################
# REFERENCE FILES (ref/)
##########################################
rule copy_gottcha_db:
    output:
        mmi   = "ref/gottcha_db.species.fna.mmi",
        stats = "ref/gottcha_db.species.fna.stats",
        tax   = "ref/gottcha_db.species.fna.tax.tsv"
    input:
        mmi   = config['global']['gottcha2_mmi'],
        stats = config['global']['gottcha2_stats'],
        tax   = config['global']['gottcha2_tax']
    shell:
        """
        mkdir -p ref
        cp {input.mmi} {output.mmi}
        cp {input.stats} {output.stats}
        cp {input.tax} {output.tax}
        """

##########################################
# STAGE INPUT FILES: Copy the FASTQ files into local working directory
##########################################
rule copy_input:
    input:
        r1_fastq,
        r2_fastq
    output:
        temp(f"{sample}_R1.fastq.gz"),
        temp(f"{sample}_R2.fastq.gz")
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
        tsv     = TSV_OUT
    params:
        prefix = sample
    threads: 16
    shell:
        """
        gottcha2.py \
          --database "ref/gottcha_db.species.fna" \
          --prefix {params.prefix} \
          -t {threads} \
          -i {input.r1} {input.r2}
        """

##########################################
# COPY FINAL RESULTS TO OUTPUT FOLDER
##########################################
rule copy_results:
    input:
        full    = FULL_OUT,
        lineage = LINEAGE_OUT,
        tsv     = TSV_OUT
    output:
        full    = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.full.tsv"),
        lineage = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.lineage.tsv"),
        tsv     = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.summary.tsv")
    shell:
        """
        mkdir -p $(dirname {output.full})
        cp {input.full} {output.full}
        cp {input.lineage} {output.lineage}
        cp {input.tsv} {output.tsv}
        """

...existing code...

rule upload_gottcha2_full:
    input:
        tsv = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.full.tsv")
    output:
        token = touch(os.path.join(gottcha2_out_dir, sample, "gottcha2_full.upload.done"))
    params:
        experiment = config["global"]["experiment"],
        sample = config["sample_name"],
        labkey_server = config["global"]["labkey_server"],
        project_name = config["global"]["labkey_project_name"],
        api_key = config["global"]["labkey_api_key"],
        db_version = config["global"]["gottcha2_db_version"]
    run:
        import pandas as pd
        from labkey.api_wrapper import APIWrapper

        # Read TSV and add required columns
        df = pd.read_csv(input.tsv, sep='\t')
        df['SAMPLE'] = params.sample
        df['EXPERIMENT'] = params.experiment
        df['GOTTCHA2_DB_VERSION'] = params.db_version

        # Initialize LabKey connection
        api = APIWrapper(
            params.labkey_server,
            params.project_name,
            api_key=params.api_key,
            use_ssl=True
        )

        # Upload in chunks of 1000 rows
        for chunk in [df[i:i+1000] for i in range(0, len(df), 1000)]:
            api.query.insert_rows(
                schema_name='lists',
                query_name='gottcha2_full',
                rows=chunk.to_dict('records')
            )