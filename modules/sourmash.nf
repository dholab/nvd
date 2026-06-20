/*
 * sourmash: sketch-based metagenomic analysis
 *
 * Experimental branchwater-backed commands live here while the sourmash branch
 * is still being shaped. Keep the interface aligned with post-preprocessing NVD
 * read tuples so later gather/tax steps can extend the same subworkflow.
 */

process SOURMASH_SKETCH_QUERY_METAGENOME {

    tag "${sample_id}"
    label "low"

    input:
    tuple val(sample_id), val(platform), val(read_structure), path(reads)

    output:
    tuple val(sample_id), val(platform), val(read_structure), path("${sample_id}.sourmash.query_metagenome.k31.scaled1000.abund.sig.gz"), emit: query_sketches

    script:
    """
    sourmash scripts singlesketch \
        ${reads} \
        --name "${sample_id}" \
        -p dna,k=31,scaled=1000,abund \
        -o "${sample_id}.sourmash.query_metagenome.k31.scaled1000.abund.sig.gz"
    """
}
