process COLLECT_CONTIGS {
    /* Collect producer-emitted contigs under stable NVD query IDs. */

    tag "${sample_id}"
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), val(producer), path(fasta)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.contigs.fasta"), path("${sample_id}.query_sequences.sqlite"), emit: contigs

    script:
    """
    collect_contigs.py \
        --sample-id '${sample_id}' \
        --input-fasta ${fasta} \
        --output-fasta ${sample_id}.contigs.fasta \
        --query-lookup ${sample_id}.query_sequences.sqlite \
        --producer '${producer}' \
        --long-contig-min-length 10000
    """
}
