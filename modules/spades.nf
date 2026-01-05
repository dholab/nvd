process RUN_SPADES {

    tag "${sample_id}"
    label "ludicrous"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), val(sample_type), path(reads)

    output:
    tuple val(sample_id), val(sample_type), path("output/contigs.fasta")

    script:
    def spades_cmd = sample_type == 'ont' || sample_type == 'sra'
        ? "spades.py --sewage --only-assembler -s"
        : "spades.py --sewage --12"
    """
    mkdir -p output
    ${spades_cmd} ${reads} -t ${task.cpus} -o output || true
    
    # Ensure contigs.fasta exists (SPAdes may not create it for low-read samples)
    if [ ! -f output/contigs.fasta ]; then
        touch output/contigs.fasta
    fi
    """
}

