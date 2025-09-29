process MAP_READS_TO_CONTIGS {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), val(platform), path(reads)
    tuple val(_sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    def preset = platform == 'ont' || platform == 'sra'
        ? "map-ont"
        : "sr"
    """
    minimap2 -ax ${preset} -t ${task.cpus} ${contigs} ${reads} \
    | samtools view -b -F ${task.cpus} \
    | samtools sort -@ ${task.cpus} -o ${sample_id}.bam && \
    samtools index ${sample_id}.bam
    """
    
}
