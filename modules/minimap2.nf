/*
 * Map reads to contigs using minimap2, with optional positional duplicate marking.
 *
 * When positional dedup is enabled (dedup_pos or dedup), the pipeline includes
 * samtools collate/fixmate/markdup to identify and remove PCR/optical duplicates.
 * This is recommended for amplicon or high-duplication libraries but adds
 * computational overhead.
 *
 * When positional dedup is disabled, reads are simply filtered (unmapped removed)
 * and coordinate-sorted, which is sufficient for many viral metagenomics applications.
 *
 * Output is always a coordinate-sorted, indexed BAM file.
 */
process MAP_READS_TO_CONTIGS {

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), val(platform), path(reads), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    def preset = platform == 'ont' || platform == 'sra'
        ? "map-ont"
        : "sr"
    def should_dedup_pos = params.dedup_pos ?: params.dedup ?: params.preprocess
    if (should_dedup_pos) {
        """
        minimap2 -ax ${preset} -t ${task.cpus} ${contigs} ${reads} \\
        | samtools view -b -F 4 \\
        | samtools collate -@ ${task.cpus} -O -u - \\
        | samtools fixmate -@ ${task.cpus} -m -u - - \\
        | samtools sort -@ ${task.cpus} -u - \\
        | samtools markdup -@ ${task.cpus} -s -r --duplicate-count - ${sample_id}.bam

        samtools index ${sample_id}.bam
        """
    } else {
        """
        minimap2 -ax ${preset} -t ${task.cpus} ${contigs} ${reads} \\
        | samtools view -b -F 4 \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.bam

        samtools index ${sample_id}.bam
        """
    }
}
