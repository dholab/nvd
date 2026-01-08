process RUN_SPADES {

    tag "${sample_id}"
    label "ludicrous"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("output/contigs.fasta")

    script:
    // Use interleaved mode (--12) only for actual interleaved paired-end reads;
    // all other cases (merged, single, ONT) use single-end mode (-s)
    def spades_cmd = read_structure == "interleaved"
        ? "spades.py --sewage --12"
        : "spades.py --sewage --only-assembler -s"
    """
    mkdir -p output
    ${spades_cmd} ${reads} -t ${task.cpus} -o output

    # SPAdes may exit 0 but not produce contigs for low-coverage samples.
    # Create empty contigs.fasta so downstream processes can handle gracefully.
    if [ ! -f output/contigs.fasta ]; then
        touch output/contigs.fasta
    fi
    """
}

