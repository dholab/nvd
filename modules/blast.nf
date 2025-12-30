/*
Perform rapid sequence similarity search using MEGABLAST.

This rule aligns the assembled contigs against a comprehensive nucleotide database
to identify similar known sequences. It's crucial for initial, broad_scale
identification of viral sequences.
*/
process MEGABLAST {
    tag "${sample_id}"
	label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(human_virus_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.megablast.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def blast_task = "megablast"
    """
    export BLASTDB=${blast_db}
    blastn -task ${blast_task} \
    -db ${blast_db}/${params.blast_db_prefix} \
    -query ${human_virus_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${params.max_blast_targets} \
    > ${sample_id}.megablast.txt
    """
}

//Create file with megablast results annotated with megablast task and full taxonomic rank.
process ANNOTATE_MEGABLAST_RESULTS {

    tag "${sample_id}"
	label "medium"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(blast_txt)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_megablast.txt"), emit: hits

    script:
    """
    annotate_blast_results.py \
    --input_file ${blast_txt} \
    --output_file ${sample_id}.annotated_megablast.txt \
    --sample_name ${sample_id} \
    --task 'megablast'
    """
}

process SELECT_TOP_BLAST_HITS {

    tag "${sample_id}"
	label "low"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(blast_txt)

    output:
    tuple val(sample_id), path("${sample_id}.top${params.blast_retention_count}_hits.txt")

    script:
    """
    select_top_blast_hits.py \\
    --input-file ${blast_txt} \\
    --output-file ${sample_id}.top${params.blast_retention_count}_hits.txt \\
    --blast-retention-count ${params.blast_retention_count}
    """
}/*

Remove any contigs that do not have at least one hit corresponding to viruses.

This rule handles cases where the input file is empty and filters out
blast query read groups where none of the entries are from superkingdom: Viruses
or where the only viral hits are phages. This will leave in any query seq that has things that are 
non virus as long as one of the qseqid reads has a viral hit.
*/
process FILTER_NON_VIRUS_MEGABLAST_NODES {

    tag "${sample_id}"
	label "low"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.mb_virus_only.txt")

    script:
    """
    filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.mb_virus_only.txt
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
    label "low"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

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
    label "high"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(classified), path(pruned_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.blastn.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def blast_task = "blastn"
    """
    export BLASTDb=${blast_db}/
    blastn -task ${blast_task} \
    -db ${blast_db}/${params.blast_db_prefix} \
    -query ${pruned_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${params.max_blast_targets} \
    > ${sample_id}.blastn.txt
    """
}

// Create file with blastn results annotated with blastn task and full taxonomic rank.
process ANNOTATE_BLASTN_RESULTS {

    tag "${sample_id}"
    label "low"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(blastn_txt)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_blastn.txt")

    script:
    """
    annotate_blast_results.py \
    --sample_name ${sample_id} \
    --input_file ${blastn_txt} \
    --output_file ${sample_id}.annotated_blastn.txt \
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
    label "low"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.nt_virus_only.txt")

    script:
    """
    filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.nt_virus_only.txt
    """
}

// Create file with merged megablast and blastn results.
process MERGE_FILTERED_BLAST_RESULTS {

    tag "${sample_id}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(megablast_viruses), path(blastn_input)

    output:
    tuple val(sample_id), path("${sample_id}_blast.merged.tsv")

    script:
    """
    concat_multiblast_hits.py \\
        --megablast-hits ${megablast_viruses} \\
        --blastn-hits ${blastn_input} \\
        --output-file ${sample_id}_blast.merged.tsv
    """
}

// Extract unclassified contigs that do not have BLAST hits.
// These could be highly divergent viruses, but are most often junk.
// Not used currently, removed from logic of nvd-2 by WG
process EXTRACT_UNCLASSIFIED_CONTIGS {

    tag "${sample_id}"
    label "low"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

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
