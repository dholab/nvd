process MEGABLAST {
    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(human_virus_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def task = "megablast"
    """
    export BLASTDB=resources/
    blastn -task ${task} \
    -db ${blast_db} \
    -query ${human_virus_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    -out .
    """
}

process ANNOTATE_MEGABLAST_RESULTS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(blast_txt), path(sqlite_file)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_megablast.txt"), emit: hits
    path "*.sqlite", emit: taxa_sqlite

    script:
    """
    annotate_blast_results.py \
    --sqlite_cache ${sqlite_file} \
    --input_file ${blast_txt} \
    --output_file ${sample_id}.annotated_megablast.txt \
    --task 'megablast'
    """
}

process FILTER_NON_VIRUS_MEGABLAST_NODES {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.virus_only.txt")

    script:
    """
    filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.virus_only.txt
    """
}

process REMOVE_MEGABLAST_MAPPED_CONTIGS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(filtered_megablast_nodes), path(human_virus_contigs)

    output:
    tuple val(sample_id), path("${sample_id}.classified.txt"), path("${sample_id}.pruned.fasta")

    script:
    """
    remove_megablast_mapped_contigs.py \
    --megablast_results ${filtered_megablast_nodes} \
    --contigs_fasta ${human_virus_contigs} \
    --pruned_contigs ${sample_id}.pruned.fasta
    """
}

process BLASTN_CLASSIFY {
    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(classified), path(pruned_contigs), path("resources/*")

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def task = "blastn"
    """
    export BLASTDB=resources/
    blastn -task ${task} \
    -db resources/${params.blast_db_prefix} \
    -query ${pruned_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    -out .
    """
}

process ANNOTATE_BLASTN_RESULTS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(blastn_txt), path(taxa_sqlite)

    output:
    tuple val(sample_id), path("${sample_id}.txt")

    script:
    """
    annotate_blast_results.py \
    --sqlite_cache ${taxa_sqlite} \
    --input_file ${blastn_txt} \
    --output_file ${sample_id}.annotated_megablast.txt \
    --task 'blastn'
    """
}

process FILTER_NON_VIRUS_BLASTN_NODES {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.virus_only.txt")

    script:
    """
    filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.virus_only.txt
    """
}

process EXTRACT_UNCLASSIFIED_CONTIGS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(contigs), path("${sample_id}_megablast.txt"), path("${sample_id}_blastn.txt")

    output:
    tuple val(sample_id), path("${sample_id}.unclassified.fasta")

    script:
    """
    extract_unclassified_contigs.py \
    --contigs ${contigs} \
    --megablast_results ${sample_id}_megablast.txt \
    --blastn_results ${sample_id}_blastn.txt \
    --unclassified ${sample_id}.unclassified.fasta
    """
}

