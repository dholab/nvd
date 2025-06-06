/*
Perform rapid sequence similarity search using MEGABLAST.

This rule aligns the assembled contigs against a comprehensive nucleotide database
to identify similar known sequences. It's crucial for initial, broad_scale
identification of viral sequences.
*/
process MEGABLAST {
    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(human_virus_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.megablast.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def blast_task = "megablast"
    """
    export BLASTDB=${blast_db}
    blastn -task ${blast_task} \
    -db ${blast_db}/${params.blast_db_prefix} \
    -query ${human_virus_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    > ${sample_id}.megablast.txt
    """
}

//Create file with megablast results annotated with megablast task and full taxonomic rank.
process ANNOTATE_MEGABLAST_RESULTS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(blast_txt), path(sqlite_file)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_megablast.txt"), emit: hits

    script:
    """
    annotate_blast_results.py \
    --sqlite_cache ${sqlite_file} \
    --input_file ${blast_txt} \
    --output_file ${sample_id}.annotated_megablast.txt \
    --sample_name ${sample_id} \
    --task 'megablast'
    """
}

/*
Remove any contigs that do not have at least one hit corresponding to viruses.

This rule handles cases where the input file is empty and filters out
groups where none of the entries are from superkingdom: Viruses
or where the only viral hits are phages.
*/
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

/*
Create a list of contigs that are not classified by megablast.

This rule generates two outputs:
1. A list of classified contigs
2. A pruned FASTA file containing only unclassified contigs

The pruned file is used as input for more sensitive blastn searching.
*/
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
    --pruned_contigs ${sample_id}.pruned.fasta \
    --classified_contigs ${sample_id}.classified.txt
    """
}

/*
Perform sequence similarity search using BLASTN.

This rule aligns the assembled contigs against a comprehensive nucleotide database
to identify similar known sequences with BLASTN
*/
process BLASTN_CLASSIFY {
    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), path(classified), path(pruned_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.blastn.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def blast_task = "blastn"
    """
    export BLASTDb=${blast_db}/
    blastn -task ${blast_task} \
    -db ${blast_db}/${params.blast_db_prefix} \
    -query ${pruned_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    > ${sample_id}.blastn.txt
    """
}

// Create file with blastn results annotated with blastn task and full taxonomic rank.
process ANNOTATE_BLASTN_RESULTS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(blastn_txt), path(taxa_sqlite)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_megablast.txt")

    script:
    """
    annotate_blast_results.py \
    --sample_name ${sample_id} \
    --sqlite_cache ${taxa_sqlite} \
    --input_file ${blastn_txt} \
    --output_file ${sample_id}.annotated_megablast.txt \
    --task 'blastn'
    """
}

/*
Remove any contigs that do not have at least one hit corresponding to viruses.

This rule handles cases where the input file is empty and filters out
groups where none of the entries are from superkingdom: Viruses
or where the only viral hits are phages.
*/
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

// Create file with merged megablast and blastn results.
process MERGE_ANNOTATED_BLAST_RESULTS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(megablast_viruses)
    tuple val(sample_id), path(blastn_viruses)

    output:
    tuple val(sample_id), path("${sample_id}_blast.merged.txt")

    script:
    """
        {{
        if [ -s ${megablast_viruses} ] || [ -s ${blastn_viruses} ]; then
            cat ${megablast_viruses} ${blastn_viruses} > ${sample_id}_blast.merged.txt
        else
            touch ${sample_id}_blast.merged.txt
            echo "Both input files are empty. Created an empty output file."
        fi
        }} > {log} 2>&1
    """
}

// Extract unclassified contigs that do not have BLAST hits.
// These could be highly divergent viruses, but are most often junk.
process EXTRACT_UNCLASSIFIED_CONTIGS {

    tag "${sample_id}"

    cpus 1

    input:
    tuple val(sample_id), path(classified_list), path(contigs), path("${sample_id}_megablast.txt"), path("${sample_id}_blastn.txt")

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
