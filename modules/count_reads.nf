process COUNT_READS {
    /*
     * Count total reads across R1 and R2 (paired) or a single file (single-end).
     * Accepts the same R1/R2-or-sentinel tuple pattern as DEACON_FILTER_HUMAN_VIRUS_READS:
     *   Paired:    (sample_id, platform, R1, R2)        — counts both files
     *   Single:    (sample_id, platform, fastq, NO_R2)  — counts one file
     */

    tag "${sample_id}"
    label 'low'

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(platform), path(reads), path(reads2)

    output:
    tuple val(sample_id), env(total_reads), emit: counts

    script:
    def is_paired = reads2.name != "NO_R2"
    if (is_paired)
        """
        r1_count=\$(zcat ${reads} | wc -l | awk -v OFMT="%.0f" '{print \$1/4}')
        r2_count=\$(zcat ${reads2} | wc -l | awk -v OFMT="%.0f" '{print \$1/4}')
        total_reads=\$((r1_count + r2_count))
        """
    else
        """
        if [[ "${reads}" == *.gz ]]; then
            total_reads=\$(zcat ${reads} | wc -l | awk -v OFMT="%.0f" '{print \$1/4}')
        else
            total_reads=\$(cat ${reads} | wc -l | awk -v OFMT="%.0f" '{print \$1/4}')
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
