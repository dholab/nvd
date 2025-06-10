process EXTRACT_HUMAN_VIRUS_CONTIGS {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(human_virus_family_hits), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.human_virus.fasta")

    script:
    """
    seqkit grep -f ${human_virus_family_hits} ${fasta} -o ${sample_id}.human_virus.fasta
    """
}

