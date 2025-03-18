import os
import glob
import re

# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Retrieve configuration parameters
# From command line
sample = config["sample_name"]
r1_fastq = config["r1_fastq"]
r2_fastq = config.get("r2_fastq", None)  # Optional for ONT
is_nanopore = config.get("is_nanopore", False)  # New parameter to identify ONT data

# From parameters file - updated paths
gottcha2_out_dir = config['global']["gottcha2_out_dir"]

# Define output paths using Python string formatting
FULL_OUT    = f"{sample}.full.tsv"
LINEAGE_OUT = f"{sample}.lineage.tsv"
TSV_OUT     = f"{sample}.tsv"
EXTRACT_FASTA = f"{sample}.extract.fasta"
PARSED_FASTA_TSV = f"{sample}.extract.tsv"

rule all:
    input:
        os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.full.tsv"),
        os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.lineage.tsv"),
        os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.summary.tsv"),
        os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.extracted_sequences.tsv"),
        os.path.join(gottcha2_out_dir, sample, "gottcha2_full.upload.done"),
        os.path.join(gottcha2_out_dir, sample, "gottcha_fasta.upload.done")

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
        r1 = r1_fastq,
        r2 = r2_fastq if not is_nanopore else []
    output:
        r1 = temp(f"{sample}_R1.fastq.gz"),
        r2 = temp(f"{sample}_R2.fastq.gz") if not is_nanopore else []
    run:
        shell(f"cp {input.r1} {output.r1}")
        if not is_nanopore:
            shell(f"cp {input.r2} {output.r2}")

##########################################
# GOTTCHA2 ANALYSIS
##########################################
rule gottcha2:
    input:
        r1 = f"{sample}_R1.fastq.gz",
        r2 = f"{sample}_R2.fastq.gz" if not is_nanopore else [],
        db_mmi   = "ref/gottcha_db.species.fna.mmi",
        db_stats = "ref/gottcha_db.species.fna.stats",
        db_tax   = "ref/gottcha_db.species.fna.tax.tsv"
    output:
        full    = FULL_OUT,
        lineage = LINEAGE_OUT,
        tsv     = TSV_OUT,
        extract = EXTRACT_FASTA
    params:
        prefix = sample,
        nanopore_flag = "--nanopore" if is_nanopore else ""
    threads: 16
    run:
        if is_nanopore:
            shell("""
            gottcha2.py \
              --database "ref/gottcha_db.species.fna" \
              --prefix {params.prefix} \
              -t {threads} \
              -i {input.r1} \
              {params.nanopore_flag} \
              -ef 20
            """)
        else:
            shell("""
            gottcha2.py \
              --database "ref/gottcha_db.species.fna" \
              --prefix {params.prefix} \
              -t {threads} \
              -i {input.r1} {input.r2} \
              -ef 20
            """)

##########################################
# PARSE EXTRACTED FASTA TO TSV
##########################################
rule parse_extract_fasta:
    input:
        extract = EXTRACT_FASTA
    output:
        tsv = PARSED_FASTA_TSV
    params:
        experiment = config["global"]["experiment"],
        sample = config["sample_name"]
    run:
        import pandas as pd
        
        # Parse FASTA file
        headers = []
        sequences = []
        taxids = []
        
        current_header = None
        current_seq = ""
        
        with open(input.extract, 'r') as fasta_file:
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    if current_header is not None:
                        headers.append(current_header)
                        sequences.append(current_seq)
                        # Extract taxid from header
                        taxid = None
                        if "TAXID=" in current_header:
                            taxid = current_header.split("TAXID=")[1].split()[0]
                        taxids.append(taxid)
                    
                    current_header = line[1:]  # Remove '>' character
                    current_seq = ""
                else:
                    current_seq += line
            
            # Don't forget the last sequence
            if current_header is not None:
                headers.append(current_header)
                sequences.append(current_seq)
                taxid = None
                if "TAXID=" in current_header:
                    taxid = current_header.split("TAXID=")[1].split()[0]
                taxids.append(taxid)
        
        # Create DataFrame
        df = pd.DataFrame({
            'experiment': params.experiment,
            'sample_id': params.sample,
            'fasta_header': headers,
            'taxid': taxids,
            'fasta_sequence': sequences
        })
        
        # Save to TSV
        df.to_csv(output.tsv, sep='\t', index=False)

##########################################
# COPY FINAL RESULTS TO OUTPUT FOLDER
##########################################
rule copy_results:
    input:
        full      = FULL_OUT,
        lineage   = LINEAGE_OUT,
        tsv       = TSV_OUT,
        fasta_tsv = PARSED_FASTA_TSV
    output:
        full      = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.full.tsv"),
        lineage   = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.lineage.tsv"),
        tsv       = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.summary.tsv"),
        fasta_tsv = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.extracted_sequences.tsv")
    shell:
        """
        mkdir -p $(dirname {output.full})
        cp {input.full} {output.full}
        cp {input.lineage} {output.lineage}
        cp {input.tsv} {output.tsv}
        cp {input.fasta_tsv} {output.fasta_tsv}
        """

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

##########################################
# UPLOAD EXTRACTED FASTA SEQUENCES TO LABKEY
##########################################
rule upload_gottcha_fasta:
    input:
        tsv = os.path.join(gottcha2_out_dir, sample, f"{sample}_gottcha2.extracted_sequences.tsv")
    output:
        token = touch(os.path.join(gottcha2_out_dir, sample, "gottcha_fasta.upload.done"))
    params:
        labkey_server = config["global"]["labkey_server"],
        project_name = config["global"]["labkey_project_name"],
        api_key = config["global"]["labkey_api_key"]
    run:
        import pandas as pd
        from labkey.api_wrapper import APIWrapper

        # Read TSV
        df = pd.read_csv(input.tsv, sep='\t')

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
                query_name='gottcha_fasta',
                rows=chunk.to_dict('records')
            )