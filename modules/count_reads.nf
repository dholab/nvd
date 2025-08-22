process COUNT_READS {
    tag "${meta}"
    label 'process_low'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(meta), val(library), path(fastq)

    output:
    tuple val(meta), env(total_reads), emit: counts

    when:
    params.tools && (params.tools.contains("nvd") || params.tools.contains("all"))

    script:
    """
    if [[ "${fastq}" == *.gz ]]; then
        total_reads=\$(zcat ${fastq} | wc -l | awk '{print \$1/4}')
    else
        total_reads=\$(cat ${fastq} | wc -l | awk '{print \$1/4}')
    fi
    """
}
