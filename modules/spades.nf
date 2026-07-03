process RUN_SPADES {

    tag "${sample_id}"
    label "high"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.spades.contigs.fasta")

    script:
    // Use interleaved mode (--12) only for actual interleaved paired-end reads.
    // The non-Illumina fallback intentionally preserves the old transitional
    // behavior by treating long reads as unpaired input. SPAdes --nanopore is
    // an auxiliary long-read library option and is invalid without a primary
    // short-read library.
    def short_read_spades_cmd = read_structure == "interleaved"
        ? "spades.py --sewage --12"
        : "spades.py --sewage --only-assembler -s"
    def spades_cmd = platform == "illumina"
        ? short_read_spades_cmd
        : "spades.py --sewage --only-assembler -s"
    """
    mkdir -p output
    ${spades_cmd} ${reads} -t ${task.cpus} -o output

    # SPAdes may exit 0 but not produce contigs for low-coverage samples.
    # Create an empty FASTA so downstream processes can handle gracefully.
    if [ -s output/contigs.fasta ]; then
        cp output/contigs.fasta ${sample_id}.spades.contigs.fasta
    else
        touch ${sample_id}.spades.contigs.fasta
    fi
    """
}
