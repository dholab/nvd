process RUN_SPADES {

    tag "${sample_id}"

    cpus 4

    input:
    tuple val(sample_id), val(sample_type), path(reads)

    output:
    tuple val(sample_id), val(sample_type), path("contigs.fasta")

    script:
    def spades_cmd = sample_type == 'ont' || sample_type == 'sra'
        ? "spades.py --sewage --only-assembler -s"
        : "spades.py --sewage --12"
    """
    (${spades_cmd} ${reads} -t ${task.cpus} -o . || true)
    """
}

