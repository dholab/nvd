// idxstats: The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments.
process COUNT_MAPPED_READS {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai"), emit: filtered_bam
    path "${sample_id}_mapped_counts.txt", emit: mapped_counts

    script:
    """
    samtools view -F 2304 -b ${bam} > ${sample_id}.filtered.bam
    samtools index ${sample_id}.filtered.bam
    samtools idxstats ${sample_id}.filtered.bam | cut -f1,3 > ${sample_id}_mapped_counts.txt 
    """
    
}
