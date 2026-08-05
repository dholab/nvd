process CONCAT_READS_AS_FASTA {

    tag "${sample_id}"
    label "medium"

    cpus 4

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.reads.fasta.gz")

    script:
    """
    seqkit scat \
        --find-only \
        --in-format fastq \
        --out-format fasta \
        --threads ${task.cpus} \
        --out-file ${sample_id}.reads.fasta.gz \
        .
    """
}
