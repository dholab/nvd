## Example commands: ##
# On dhogal2 and other machines when using Apptainer
# Uses per-rule Apptainer to prevent dependency collisions
# conda activate snakemake
# snakemake --cores 40 --config run_id=`date +%s` --software-deployment-method apptainer
#
# On laptops and other machines when using monolithic Docker image
# docker load -i docker/nvd.tar
# docker run -it --user $(id -u):$(id -g) -v $(pwd)/:/scratch -w /scratch nvd:30575
# snakemake --cores 10 --config run_id=`date +%s`

# Import utility functions
from scripts.workflow_functions import (
    get_input_files,
    get_input_type,
    get_sra_accession,
    create_directory,
    upload_file
)

import re

# Use samples from configfile: "config.yaml"
configfile: "config/config.yaml"

# Normalize sample names by replacing dashes, spaces, and special characters with underscores
for sample in config['samples']:
    sample['name'] = re.sub(r'[\s\-]', '_', sample['name'])

# Timestamp set when running starting snakemake run
snakemake_run_id = str(config['run_id'])

# From config.yaml globals
experiment = config['global']['experiment']
blast_db_path = config['global']['blast_db_path'] # downlaoded from NCBI
blast_db_name = config['global']['blast_db_name'] # basename of BLAST db
stat_dbss = config['global']['stat_dbss'] # provided by Kenneth Katz, NCBI
stat_index = config['global']['stat_index'] # provided by Kenneth Katz, NCBI
stat_annotation = config['global']['stat_annotation'] # provided by Kenneth Katz, NCBI
final_output_path = f"{config['global']['out_dir']}/{experiment}" # for local tarball of results
human_virus_taxlist = config['global']['human_virus_taxlist'] # human virus family taxa + arteriviridae

# Path to results directory
local_output_dir = 'results'
log_dir = 'logs'

rule all:
    input:
        expand([
            f"{{local_output_dir}}/09_extract_human_virus_contigs/{{sample}}.fa",
            f"{{local_output_dir}}/18_extract_unclassified_contigs/{{sample}}.unclassified.fa",
            f"{{local_output_dir}}/19_labkey/{{sample}}.blast.tkn",
            f"{{local_output_dir}}/19_labkey/{{sample}}.fasta.tkn",
            f"{{local_output_dir}}/20_map_reads_to_contigs/{{sample}}.bam",
            f"{{local_output_dir}}/21_create_tarball/{{sample}}.tar.zst",
            f"{{local_output_dir}}/19_labkey/{{sample}}.upload.done"
        ],
        local_output_dir=local_output_dir,
        sample=[s['name'] for s in config['samples']]
    )

rule decompress_resources:
    """
    Copies resources.tar.zst from staging (exact file specified in config.yaml)
    into the local 'resources/' directory and decompresses it.
    """
    input:
        tar_file = config["global"]["resources_tar_zst"]
    output:
        "decompressed.txt",
        stat_dbss,
        stat_index,
        stat_annotation,
        human_virus_taxlist
    log:
        f"{log_dir}/00_decompress_resources.log"
    shell:
        """
        mkdir -p resources \
        && cp {input.tar_file} resources.tar.zst \
        && tar -I zstd -xvf resources.tar.zst \
        && rm -f resources.tar.zst \
        && touch "decompressed.txt"
        """

rule prepare_input:
    """
    This workflow expects FASTQ sequences
    These can be from paired-end Illumina reads, ONT sequences, or SRA accessions
    If they are SRA accessions, they will be treated as ONT sequences
    """
    input:
        files = lambda wildcards: get_input_files(wildcards, config)
    output:
        f"{local_output_dir}/01_prepare_input/{{sample}}.fq.gz"
    params:
        input_type = lambda wildcards: get_input_type(wildcards, config),
        sra_accession = lambda wildcards: get_sra_accession(wildcards, config) if get_input_type(wildcards, config) == "sra" else None,
        temp_dir = lambda wildcards: f"{wildcards.sample}.temp",
        sample = "{sample}"
    threads: 1
    log:
        f"{log_dir}/01_prepare_input/{{sample}}.log"
    benchmark:
        f"{log_dir}/01_prepare_input/{{sample}}.txt"
    script:
        "scripts/prepare_input.py"

rule extract_potential_human_virus_family_reads:
    """
    Filter to extract reads that are classified as human-infecting virus families by STAT
    There will be false positives that will need to be dealt with later
    However, this step dramatically reduces the number of reads that need to be examined.
    In testing, 88m reads were reduced to 82.5k reads
    This also dramatically reduces memory requirements which is important for laptop-scale analyses
    And saves a lot of time even on larger servers.

    - Use seqkit to convert FASTQ to FASTA
    - Use aligns_to to classify reads
    - Extract first column with read name
    - Use seqkit to extract reads that align to human-infecting virus families
    - Save output as FASTQ files that can be used as SPAdes input
    """
    input:
        reads = f"{local_output_dir}/01_prepare_input/{{sample}}.fq.gz",
        human_virus_taxlist = human_virus_taxlist,
        stat_dbss = stat_dbss,
        decompressed = "decompressed.txt"
    output:
        f"{local_output_dir}/02_extract_human_virus_reads/{{sample}}.fq.gz"
    log: f"{log_dir}/02_extract_human_virus_reads/{{sample}}.log"
    benchmark: f"{log_dir}/02_extract_human_virus_reads/{{sample}}.txt"
    threads: 4
    shell:
        """
        {{
        seqkit fq2fa --threads 1 {input.reads} \
        | aligns_to.3.1.1 \
        -dbss {input.stat_dbss} \
        -num_threads 2 \
        -tax_list {input.human_virus_taxlist} \
        stdin \
        | cut -f1 \
        | seqkit grep -f - \
        {input.reads} \
        -o {output}
        }} > {log} 2>&1
        """

rule run_spades:
    """
    Assemble putative human virus family reads.
    Sewage mode gives excellent results with ONT and Illumina reads.
    """
    input:
        reads = f"{local_output_dir}/02_extract_human_virus_reads/{{sample}}.fq.gz"
    output:
        touch(f"{local_output_dir}/03_de_novo_assemble/{{sample}}/assembly/contigs.fasta")
    threads: 4
    log:
        f"{log_dir}/03_de_novo_assemble/{{sample}}_spades.log"
    benchmark:
        f"{log_dir}/03_de_novo_assemble/{{sample}}_spades.txt"
    params:
        outdir = f"{local_output_dir}/03_de_novo_assemble/{{sample}}/assembly",
        spades_cmd = lambda wildcards, input: (
            "spades.py --sewage --only-assembler -s"
            if get_input_type(wildcards, config) in ['ont', 'sra'] else
            "spades.py --sewage --12"
        )
    shell:
        """
        {{
        mkdir -p {params.outdir}
        ({params.spades_cmd} {input.reads} -t {threads} -o {params.outdir} || true) 2>&1
        }} > {log} 2>&1
        """

rule move_spades_output:
    """
    Moves SPAdes contigs to reliably accessible location or create an empty file if assembly failed.
    """
    input:
        contigs = f"{local_output_dir}/03_de_novo_assemble/{{sample}}/assembly/contigs.fasta"
    output:
        final_contigs = f"{local_output_dir}/03_de_novo_assemble/{{sample}}.fa"
    params:
        assembly_dir = f"{local_output_dir}/03_de_novo_assemble/{{sample}}"
    log:
        f"{log_dir}/03_de_novo_assemble/{{sample}}_move.log"
    shell:
        """
        {{
        if [ -s {input.contigs} ]; then
            mv {input.contigs} {output.final_contigs}
            rm -rf {params.assembly_dir}
        else
            touch {output.final_contigs}
            echo "No contigs were assembled or SPAdes encountered an error. Created an empty file."
        fi
        }} > {log} 2>&1
        """

rule mask_low_complexity_contigs:
    """
    Use bbmask to N-mask low complexity sequences that interfere with classification.
    """
    input:
        contigs = f"{local_output_dir}/03_de_novo_assemble/{{sample}}.fa"
    output:
        masked_contigs = f"{local_output_dir}/04_mask_low_complexity_contigs/{{sample}}.fa"
    threads: 1
    resources:
        mem_mb = 8192  # 8GB of memory
    log:
        f"{log_dir}/04_mask_low_complexity_contigs/{{sample}}.log"
    benchmark:
        f"{log_dir}/04_mask_low_complexity_contigs/{{sample}}.txt"
    params:
        entropy = 0.9
    shell:
        """
        {{
        if [ -s {input.contigs} ]; then
            bbmask.sh -Xmx{resources.mem_mb}m \
                in={input.contigs} \
                out={output.masked_contigs} \
                entropy={params.entropy}
        else
            touch {output.masked_contigs}
            echo "Input file is empty. Created an empty output file."
        fi
        }} > {log} 2>&1
        """

rule exclude_short_contigs:
    """
    Remove N from ends of sequences and require at least 200bp of called bases after trimming.
    This step is necessary because otherwise BLAST can hang.
    """
    input:
        masked_contigs = f"{local_output_dir}/04_mask_low_complexity_contigs/{{sample}}.fa"
    output:
        filtered_contigs = touch(f"{local_output_dir}/05_exclude_short_contigs/{{sample}}.fa")
    threads: 1
    resources:
        mem_mb = 8192  # 8GB of memory
    log:
        f"{log_dir}/05_exclude_short_contigs/{{sample}}.log"
    benchmark:
        f"{log_dir}/05_exclude_short_contigs/{{sample}}.txt"
    params:
        min_consecutive_bases = 200,
        qtrim = "t"
    shell:
        """
        {{       
        if [ -s {input.masked_contigs} ]; then
            reformat.sh -Xmx{resources.mem_mb}m \
                in={input.masked_contigs} \
                out={output.filtered_contigs} \
                qtrim={params.qtrim} \
                minconsecutivebases={params.min_consecutive_bases} 
        else
            touch {output.filtered_contigs}
            echo "Input file is empty. Created an empty output file." 
        fi
        }} > {log} 2>&1
        """

rule classify_contigs_first_pass:
    """
    Classify contigs using NCBI STAT tree_index database
    "Coarse" classification that will be refined in the second pass
    Classification helps identify "viral" contigs that could be better assigned to other taxa
    """
    input:
        fasta_file = f"{local_output_dir}/05_exclude_short_contigs/{{sample}}.fa",
        tree_filter = stat_index,
        decompressed = "decompressed.txt"
    output:
        first_pass_stat_file = f"{local_output_dir}/06_classify_contigs/{{sample}}.firstpass.txt"
    threads: workflow.cores
    log:
        f"{log_dir}/06_classify_contigs/{{sample}}.classify_contigs_first_pass.log"
    benchmark:
        f"{log_dir}/06_classify_contigs/{{sample}}.classify_contigs_first_pass.txt"
    shell:
        "aligns_to.3.1.1 -dbs {input.tree_filter}  "
        "-num_threads {threads} {input.fasta_file} > {output.first_pass_stat_file} 2> {log}"

rule generate_contigs_tax_list:
    """
    Create list of used taxa based on STAT classification of contigs
    """
    input:
        first_pass_stat_file = f"{local_output_dir}/06_classify_contigs/{{sample}}.firstpass.txt"
    output:
        tax_list = f"{local_output_dir}/06_classify_contigs/tax_list/{{sample}}.tax_list"
    log:
        f"{log_dir}/06_classify_contigs/{{sample}}.generate_tax_list.log"
    benchmark:
        f"{log_dir}/06_classify_contigs/{{sample}}.generate_tax_list.txt"
    script:
        "scripts/generate_tax_list.py"

rule classify_contigs_second_pass:
    """
    Create a second pass classification of contigs using the tax_list generated in the first pass
    By using -dbss, only the taxa that are present in the contigs are further analyzed
    """
    input:
        tax_list = f"{local_output_dir}/06_classify_contigs/tax_list/{{sample}}.tax_list",
        fasta_file = f"{local_output_dir}/05_exclude_short_contigs/{{sample}}.fa",
        tree_filter = stat_dbss,
        decompressed = "decompressed.txt"
    output:
        second_pass_stat_file = f"{local_output_dir}/06_classify_contigs/{{sample}}.secondpass.txt"
    threads: workflow.cores
    log:
        f"{log_dir}/06_classify_contigs/{{sample}}.classify_reads_second_pass.log"
    benchmark:
        f"{log_dir}/06_classify_contigs/{{sample}}.classify_reads_second_pass.txt"
    shell:
        """
        {{
        aligns_to.3.1.1 \
            -tax_list {input.tax_list} \
            -dbss {input.tree_filter} \
            -num_threads 1 \
            {input.fasta_file} \
            > {output.second_pass_stat_file}
        }} > {log} 2>&1
        """

rule generate_stat_contig_report:
    """
    generate NCBI STAT report that classifies contig hits
    modify default threshold to >0.001 to report hits more sensitively
    """
    input:
        f"{local_output_dir}/06_classify_contigs/{{sample}}.secondpass.txt"
    output:
        f"{local_output_dir}/07_generate_stat_contig_report/{{sample}}.report"
    threads: 1
    log:
        f"{log_dir}/07_generate_stat_contig_report/{{sample}}.log"
    benchmark:
        f"{log_dir}/07_generate_stat_contig_report/{{sample}}.txt"
    params:
        cutoff_percent = 0.001,
        gettax_sqlite_path = config["global"]["gettax_sqlite_path"]
    script:
        "scripts/hits_to_report.py"

rule identify_human_virus_family_contigs:
    """
    Identifies contigs potentially originating from human viruses.
    This rule filters the classified contigs to focus on those matching known human virus families.
    
    Uses same container as hits_to_report.py since it shares dependencies
    """
    input:
        hits = f"{local_output_dir}/06_classify_contigs/{{sample}}.secondpass.txt"
    output:
        filtered = touch(f"{local_output_dir}/08_identify_human_virus_family_contigs/{{sample}}.txt")
    threads: workflow.cores
    log: f"{log_dir}/08_identify_human_virus_family_contigs/{{sample}}.log"
    benchmark: f"{log_dir}/08_identify_human_virus_family_contigs/{{sample}}.txt"
    params:
        taxa = ["Adenoviridae", "Anelloviridae", "Arenaviridae", "Arteriviridae", 
                "Astroviridae", "Bornaviridae", "Peribunyaviridae", "Caliciviridae", 
                "Coronaviridae", "Filoviridae", "Flaviviridae", "Hepadnaviridae", 
                "Hepeviridae", "Orthoherpesviridae", "Orthomyxoviridae", "Papillomaviridae", 
                "Paramyxoviridae", "Parvoviridae", "Picobirnaviridae", "Picornaviridae", 
                "Pneumoviridae", "Polyomaviridae", "Poxviridae", "Sedoreoviridae", 
                "Retroviridae", "Rhabdoviridae", "Togaviridae", "Kolmioviridae"],
        stringency = 0.7,
        include_children = True
    script:
        "scripts/extract_taxa_spots.py"

rule extract_human_virus_contigs:
    """
    use seqkit to extract contigs from original data that have kmers matching human virus families
    """
    input:
        human_virus_family_hits=f"{local_output_dir}/08_identify_human_virus_family_contigs/{{sample}}.txt",
        fastq=f"{local_output_dir}/05_exclude_short_contigs/{{sample}}.fa"
    output:
        touch(f"{local_output_dir}/09_extract_human_virus_contigs/{{sample}}.fa")
    threads: 1
    log: f"{log_dir}/09_extract_human_virus_contigs/{{sample}}.log"
    benchmark: f"{log_dir}/09_extract_human_virus_contigs/{{sample}}.txt"
    shell:
        """
        seqkit grep -f {input.human_virus_family_hits} {input.fastq} -o {output} 2> {log}
        """

## NCBI BLAST classification of contigs
# For putative human virus family contigs, use two-stage BLAST analysis
# First, analyze contigs with megablast to find close matches to existing targets in core-nt
# Then inspect these results, filtering out those that are phages or are non-viral false positives from STAT
# Next, identify contigs that were not classified by megablast.
# Attempt to classify these contigs with blastn, performing similar filtering on the output 
# Merge the megablast and blastn output into a single file per sample. 

rule megablast:
    """
    Perform rapid sequence similarity search using MEGABLAST.

    This rule aligns the assembled contigs against a comprehensive nucleotide database
    to identify similar known sequences. It's crucial for initial, broad_scale
    identification of viral sequences.
    """
    input:
        contigs = f"{local_output_dir}/09_extract_human_virus_contigs/{{sample}}.fa",
        decompressed = "decompressed.txt"
    output:
        blast_results = touch(f"{local_output_dir}/10_megablast/{{sample}}.txt")
    params:
        db = blast_db_path,
        outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids",
        max_target_seqs = 5,
        task = "megablast",
        blastdb = os.environ.get("BLASTDB", "resources/") 
    threads: 4
    log:
        f"{log_dir}/10_megablast/{{sample}}.log"
    benchmark:
        f"{log_dir}/10_megablast/{{sample}}.benchmark.txt"
    shell:
        """
        {{
        export BLASTDB={params.blastdb}
        if [ -s {input.contigs} ]; then
            blastn -task {params.task} \
                -db {params.db} \
                -query {input.contigs} \
                -num_threads {threads} \
                -outfmt "{params.outfmt}" \
                -max_target_seqs {params.max_target_seqs} \
                -out {output.blast_results}
        else
            touch {output.blast_results}
            echo "Input file is empty. Created an empty output file."
        fi
        }} > {log} 2>&1
        """

rule annotate_megablast_results:
    """
    Create file with megablast results annotated with megablast task and full taxonomic rank.
    
    Note:
        gettax.sqlite is created by the hits_to_report.py script and is used to extract full rank
    """
    input:
        blast_results=f"{local_output_dir}/10_megablast/{{sample}}.txt"
    output:
        annotated_results=touch(f"{local_output_dir}/11_annotate_megablast_results/{{sample}}.txt")
    params:
        sample="{sample}",
        gettax_sqlite_path=config["global"]["gettax_sqlite_path"],
        task="megablast"
    threads: 1
    resources:
        mem_mb=1000
    log:
        f"{log_dir}/11_annotate_megablast_results/{{sample}}.log"
    benchmark:
        f"{log_dir}/11_annotate_megablast_results/{{sample}}.benchmark.txt"
    script:
        "scripts/annotate_blast_results.py"

rule filter_non_virus_megablast_nodes:
    """
    Remove any contigs that do not have at least one hit corresponding to viruses.
    
    This rule handles cases where the input file is empty and filters out
    groups where none of the entries are from superkingdom: Viruses
    or where the only viral hits are phages.
    """
    input:
        annotated_results=f"{local_output_dir}/11_annotate_megablast_results/{{sample}}.txt"
    output:
        filtered_results=touch(f"{local_output_dir}/12_filter_non_virus_megablast_nodes/{{sample}}.txt")
    log:
        f"{log_dir}/12_filter_non_virus_megablast_nodes/{{sample}}.log"
    benchmark:
        f"{log_dir}/12_filter_non_virus_megablast_nodes/{{sample}}.benchmark.txt"
    script:
        "scripts/filter_non_virus_megablast_nodes.py"

rule remove_megablast_mapped_contigs:
    """
    Create a list of contigs that are not classified by megablast.
    
    This rule generates two outputs:
    1. A list of classified contigs
    2. A pruned FASTA file containing only unclassified contigs
    
    The pruned file is used as input for more sensitive blastn searching.
    """
    input:
        megablast_results=f"{local_output_dir}/10_megablast/{{sample}}.txt",
        contigs_fasta=f"{local_output_dir}/09_extract_human_virus_contigs/{{sample}}.fa"
    output:
        classified_contigs=f"{local_output_dir}/13_remove_megablast_mapped_contigs/{{sample}}.classified.txt",
        pruned_contigs=f"{local_output_dir}/13_remove_megablast_mapped_contigs/{{sample}}.pruned.fa"
    threads: 1
    log:
        f"{log_dir}/13_remove_megablast_mapped_contigs/{{sample}}.log"
    benchmark:
        f"{log_dir}/13_remove_megablast_mapped_contigs/{{sample}}.benchmark.txt"
    script:
        "scripts/remove_megablast_mapped_contigs.py"

rule blastn_classify:
    """
    Perform sequence similarity search using BLASTN.

    This rule aligns the assembled contigs against a comprehensive nucleotide database
    to identify similar known sequences with BLASTN
    """
    input:
        contigs = f"{local_output_dir}/13_remove_megablast_mapped_contigs/{{sample}}.pruned.fa",
        decompressed = "decompressed.txt"
    output:
        blast_results = touch(f"{local_output_dir}/14_blastn_classify/{{sample}}.txt")
    params:
        db = blast_db_path,
        outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids",
        max_target_seqs = 5,
        task = "blastn",
        blastdb = os.environ.get("BLASTDB", "resources/") 
    threads: 4
    log:
        f"{log_dir}/14_blastn_classify/{{sample}}.log"
    benchmark:
        f"{log_dir}/14_blastn_classify/{{sample}}.benchmark.txt"
    shell:
        """
        {{
        export BLASTDB={params.blastdb}
        if [ -s {input.contigs} ]; then
            blastn -task {params.task} \
                -db {params.db} \
                -query {input.contigs} \
                -num_threads {threads} \
                -outfmt "{params.outfmt}" \
                -max_target_seqs {params.max_target_seqs} \
                -out {output.blast_results}
        else
            touch {output.blast_results}
            echo "Input file is empty. Created an empty output file."
        fi
        }} > {log} 2>&1
        """

rule annotate_blastn_results:
    """
    Create file with blastn results annotated with blastn task and full taxonomic rank.
    
    Note:
        gettax.sqlite is created by the hits_to_report.py script.
    """
    input:
        blast_results = f"{local_output_dir}/14_blastn_classify/{{sample}}.txt"
    output:
        annotated_results = touch(f"{local_output_dir}/15_annotate_blastn_results/{{sample}}.txt")
    params:
        sample = "{sample}",
        gettax_sqlite_path = config["global"]["gettax_sqlite_path"],
        task = "blastn"
    threads: 1
    log:
        f"{log_dir}/15_annotate_blastn_results/{{sample}}.log"
    benchmark:
        f"{log_dir}/15_annotate_blastn_results/{{sample}}.benchmark.txt"
    script:
        "scripts/annotate_blast_results.py"

rule filter_non_virus_blastn_nodes:
    """
    Remove any contigs that do not have at least one hit corresponding to viruses.
    
    This rule handles cases where the input file is empty and filters out
    groups where none of the entries are from superkingdom: Viruses
    or where the only viral hits are phages.
    """
    input:
        annotated_results = f"{local_output_dir}/15_annotate_blastn_results/{{sample}}.txt"
    output:
        filtered_results = touch(f"{local_output_dir}/16_filter_non_virus_blastn_nodes/{{sample}}.txt")
    log:
        f"{log_dir}/16_filter_non_virus_blastn_nodes/{{sample}}.log"
    threads: 1
    benchmark:
        f"{log_dir}/16_filter_non_virus_blastn_nodes/{{sample}}.benchmark.txt"
    script:
        "scripts/filter_non_virus_megablast_nodes.py"

rule merge_annotated_blast_results:
    """
    Create file with merged megablast and blastn results.
    """
    input:
        megablast = f"{local_output_dir}/12_filter_non_virus_megablast_nodes/{{sample}}.txt",
        blastn = f"{local_output_dir}/16_filter_non_virus_blastn_nodes/{{sample}}.txt"
    output:
        merged = f"{local_output_dir}/17_merge_annotated_blast_results/{{sample}}.txt"
    threads: 1
    log:
        f"{log_dir}/17_merge_annotated_blast_results/{{sample}}.log"
    benchmark:
        f"{log_dir}/17_merge_annotated_blast_results/{{sample}}.benchmark.txt"
    shell:
        """
        {{
        if [ -s {input.megablast} ] || [ -s {input.blastn} ]; then
            cat {input.megablast} {input.blastn} > {output.merged}
        else
            touch {output.merged}
            echo "Both input files are empty. Created an empty output file."
        fi
        }} > {log} 2>&1
        """

rule extract_unclassified_contigs:
    """
    Extract unclassified contigs that do not have BLAST hits.
    These could be highly divergent viruses, but are most often junk.
    """
    input:
        contigs = f"{local_output_dir}/09_extract_human_virus_contigs/{{sample}}.fa",
        megablast_results = f"{local_output_dir}/10_megablast/{{sample}}.txt",
        blastn_results = f"{local_output_dir}/14_blastn_classify/{{sample}}.txt"
    output:
        unclassified = touch(f"{local_output_dir}/18_extract_unclassified_contigs/{{sample}}.unclassified.fa")
    threads: 1
    resources:
        mem_mb = 4000
    log:
        f"{log_dir}/18_extract_unclassified_contigs/{{sample}}.log"
    benchmark:
        f"{log_dir}/18_extract_unclassified_contigs/{{sample}}.benchmark.txt"
    script:
        "scripts/extract_unclassified_contigs.py"

## Load results into LabKey ##
# Import the BLAST output into a LabKey table for exploration
# Import contig FASTA sequences into a LabKey table for easy access

rule labkey_upload_blast:
    """
    Uploads BLAST results to LabKey or saves to CSV if LabKey API is missing.
    """
    input:
        blast_results = f"{local_output_dir}/17_merge_annotated_blast_results/{{sample}}.txt",
        mapped_reads = f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}_mapped_counts.txt",
        fastq_file = f"{local_output_dir}/01_prepare_input/{{sample}}.fq.gz"
    output:
        token = touch(f"{local_output_dir}/19_labkey/{{sample}}.blast.tkn"),
        csv_file = f"{config['global']['out_dir']}/{{sample}}.csv"  # Ensures the CSV is always an output
    threads: 1
    resources:
        mem_mb = 4000
    params:
        experiment = config["global"]["experiment"],
        blast_db_name = config["global"]["blast_db_name"],
        snakemake_run_id = snakemake_run_id,
        labkey_server = config["global"].get("labkey_server", ""),
        project_name = config["global"].get("labkey_project_name", ""),
        api_key = config["global"].get("labkey_api_key", ""),
        stat_db_version = config["global"]["stat_dbss"],
        out_dir = config["global"]["out_dir"]
    log:
        f"{log_dir}/19_labkey/{{sample}}.blast.log"
    benchmark:
        f"{log_dir}/19_labkey/{{sample}}.blast.benchmark.txt"
    script:
        "scripts/labkey_upload.py"

rule labkey_upload_fasta:
    """
    Uploads FASTA results to LabKey data explorer
    """
    input:
        fasta = f"{local_output_dir}/09_extract_human_virus_contigs/{{sample}}.fa"
    output:
        token = touch(f"{local_output_dir}/19_labkey/{{sample}}.fasta.tkn")
    threads: 1
    log:
        f"{log_dir}/19_labkey{{sample}}.fasta.log"
    benchmark:
        f"{log_dir}/19_labkey/{{sample}}.fasta.txt"
    params:
        experiment = config["global"]["experiment"],
        labkey_server = config["global"].get("labkey_server", ""),
        project_name = config["global"].get("labkey_project_name", ""),        
        api_key = config["global"].get("labkey_api_key", ""),
        snakemake_run_id = snakemake_run_id,
        singleline_fasta_path = f"{local_output_dir}/19_labkey/{{sample}}.fasta"
    script:
        "scripts/labkey_upload_fasta.py"

rule map_reads_to_contigs:
    """
    Maps potential human virus reads to identified SPAdes contigs.
    It's important for validating the assembled viral sequences and assessing their coverage.
    """
    input:
        contigs = f"{local_output_dir}/09_extract_human_virus_contigs/{{sample}}.fa",
        reads = f"{local_output_dir}/02_extract_human_virus_reads/{{sample}}.fq.gz"
    output:
        bam = f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}.bam",
        bai = f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}.bam.bai"
    threads: 4
    log:
        f"{log_dir}/20_map_reads_to_contigs/{{sample}}.map.log"
    benchmark:
        f"{log_dir}/20_map_reads_to_contigs/{{sample}}.map.benchmark.txt"
    params:
        mapper = lambda wildcards, input: "map-ont" if get_input_type(wildcards, config) in ['ont', 'sra'] else "sr"
    shell:
        """
        {{
        (minimap2 -ax {params.mapper} -t {threads} {input.contigs} {input.reads} | \
        samtools view -b -F 4 | \
        samtools sort -@ {threads} -o {output.bam} && \
        samtools index {output.bam})
        }} > {log} 2>&1
        """

rule count_mapped_reads:
    """
    Counts the number of reads mapped to each contig in the BAM file.
    Outputs a two_column file with contig name and the number of mapped reads.
    samtools flag 2304 to retain only primary alignments in counting
    """
    input:
        bam = f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}.bam"
    output:
        filtered_bam=temp(f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}.filtered.bam"),
        filtered_bam_index=temp(f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}.filtered.bam.bai"),
        counts = f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}_mapped_counts.txt"
    log:
        f"{log_dir}/20_map_reads_to_contigs/{{sample}}.count.log"
    benchmark:
        f"{log_dir}/20_map_reads_to_contigs/{{sample}}.count.benchmark.txt"
    shell:
        """
        {{
        samtools view -F 2304 -b {input.bam} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        samtools idxstats {output.filtered_bam} | awk '{{print $1, $3}}' > {output.counts}
        }} > {log} 2>&1
        """

rule create_tarball:
    input:
        virus_contigs = f"{local_output_dir}/09_extract_human_virus_contigs/{{sample}}.fa",
        unclassified_contigs = f"{local_output_dir}/18_extract_unclassified_contigs/{{sample}}.unclassified.fa",
        mapped_reads = f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}.bam",
        mapped_reads_index = f"{local_output_dir}/20_map_reads_to_contigs/{{sample}}.bam.bai"
    output:
        tarball = f"{local_output_dir}/21_create_tarball/{{sample}}.tar.zst"
    log:
        f"{log_dir}/21_create_tarball/{{sample}}.log"
    shell:
        """
        {{
            # Initialize an empty list to hold the valid files
            valid_files=""
            
            # Check each file for existence and non-zero size
            for file in {input.virus_contigs} {input.unclassified_contigs} {input.mapped_reads} {input.mapped_reads_index}; do
                if [ -s "$file" ]; then
                    valid_files="$valid_files $file"
                else
                    echo "Skipping missing or empty file: $file"
                fi
            done
            
            # Only create the tarball if there are valid files
            if [ -n "$valid_files" ]; then
                tar -I zstd -cvf {output.tarball} $valid_files
            else
                echo "No valid files to include in the tarball."
            fi
        }} > {log} 2>&1
        """
        
rule upload_files_to_labkey:
    input:
        tarball = rules.create_tarball.output.tarball,
        config_file = "config/config.yaml",
    output:
        done = f"{local_output_dir}/19_labkey/{{sample}}.upload.done"
    params:
        labkey_server = config["global"].get("labkey_server", ""),
        project_name = config["global"].get("labkey_project_name", ""),   
        webdav_url = config["global"].get("webdav_url", ""),
        username = config["global"].get("labkey_username", ""),
        password = config["global"].get("labkey_password", ""),
        experiment = config["global"]["experiment"],
        snakemake_run_id = snakemake_run_id,
        final_output_path = config["global"]["out_dir"]
    log:
        f"{log_dir}/19_labkey/{{sample}}.files.log"
    script:
        "scripts/upload_files_to_labkey.py"

## Cleanup ##
# Utility rule to run in isolation to remove existing output files.

rule clean:
    """
    Remove all generated files and directories.
    """
    params:
        local_output_dir = local_output_dir,
        dirs_to_remove = [
            "results",
            "fasterq.tmp*",
            "*.temp"
        ]
    script:
        "scripts/clean.py"
