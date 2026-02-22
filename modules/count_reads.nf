process COUNT_READS {
    tag "${sample_id}"
    label 'low'

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(fastq)

    output:
    tuple val(sample_id), env(total_reads), emit: counts

    script:
    """
    if [[ "${fastq}" == *.gz ]]; then
        total_reads=\$(zcat ${fastq} | wc -l | awk -v OFMT="%.0f" '{print \$1/4}')
    else
        total_reads=\$(cat ${fastq} | wc -l | awk -v OFMT="%.0f" '{print \$1/4}')
    fi

    """
}

// Lightweight threshold check: reads only min_reads records instead of the
// entire file. For large gzipped FASTQs this is O(threshold) vs O(file_size).
process CHECK_MIN_READS {
    tag "${sample_id}"
    label 'low'

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(fastq)

    output:
    tuple val(sample_id), env(read_count), emit: counts

    script:
    """
    read_count=\$(seqkit head -n ${params.min_gottcha_reads} ${fastq} | seqkit stats -T | tail -1 | cut -f4)
    """
}
