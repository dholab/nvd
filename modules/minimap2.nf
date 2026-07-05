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
process MAP_SINGLE_READS {

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    // Platform describes sequencing chemistry; SRA is only an input source.
    def preset = platform == 'ont'
        ? "map-ont"
        : "sr"
    def should_dedup_pos = params.dedup || params.dedup_pos
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

process MAP_PAIRED_READS {

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    // Platform describes sequencing chemistry; SRA is only an input source.
    def preset = platform == 'ont'
        ? "map-ont"
        : "sr"
    def should_dedup_pos = params.dedup || params.dedup_pos
    def read_group = "@RG\\tID:unmerged_reads\\tSM:${sample_id}\\tDS:evidence_class=single_read"
    if (should_dedup_pos) {
        """
        minimap2 -ax ${preset} -R '${read_group}' -t ${task.cpus} ${contigs} ${reads} \
        | samtools view -b -F 4 \
        | samtools collate -@ ${task.cpus} -O -u - \
        | samtools fixmate -@ ${task.cpus} -m -u - - \
        | samtools sort -@ ${task.cpus} -u - \
        | samtools markdup -@ ${task.cpus} -s -r --duplicate-count - ${sample_id}.bam

        samtools index ${sample_id}.bam
        """
    } else {
        """
        minimap2 -ax ${preset} -R '${read_group}' -t ${task.cpus} ${contigs} ${reads} \
        | samtools view -b -F 4 \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.bam

        samtools index ${sample_id}.bam
        """
    }
}

process EXTRACT_UNMAPPED_READS {
    /* Extract reads that do not map back to screened contigs. */

    tag "${sample_id}"
    label "medium"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads), path(contigs)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.mapback_unmapped.fastq.gz"), emit: reads
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.mapback_unmapped_counts.tsv"), emit: counts

    script:
    // Platform describes sequencing chemistry; SRA is only an input source.
    def preset = platform == 'ont'
        ? "map-ont"
        : "sr"
    """
    minimap2 -ax ${preset} -t ${task.cpus} ${contigs} ${reads} \
    | samtools view -b -o ${sample_id}.mapback_all.bam

    unmapped_reads=\$(samtools view -c -f 4 ${sample_id}.mapback_all.bam)
    samtools fastq -f 4 ${sample_id}.mapback_all.bam \
    | gzip -c > ${sample_id}.mapback_unmapped.fastq.gz

    printf 'sample_id\tplatform\tread_structure\tunmapped_reads\n' > ${sample_id}.mapback_unmapped_counts.tsv
    printf '${sample_id}\t${platform}\t${read_structure}\t%s\n' "\${unmapped_reads}" >> ${sample_id}.mapback_unmapped_counts.tsv
    printf 'nvd.mapback_unmapped_reads sample_id=${sample_id} platform=${platform} read_structure=${read_structure} unmapped_reads=%s\n' "\${unmapped_reads}" >&2
    """
}
