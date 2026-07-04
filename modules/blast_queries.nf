process SELECT_BLAST_QUERIES {
    /* Select screened records admitted to BLAST by class-aware rules. */

    tag "${sample_id}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(fasta), path(lookup)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.blast_queries.fasta"), path(lookup), emit: queries

    script:
    """
    select_blast_queries.py \
        --input-fasta ${fasta} \
        --contig-lookup ${lookup} \
        --output-fasta ${sample_id}.blast_queries.fasta
    """
}
