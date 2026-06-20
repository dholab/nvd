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

process SOURMASH_COLLECT_QUERY_SKETCHES {

    tag "query_metagenomes"
    label "low"

    input:
    path(query_sketches)

    output:
    path("query_metagenomes.k31.scaled1000.zip"), emit: collection

    script:
    def sketch_list = query_sketches.collect { sketch -> "\"${sketch}\"" }.join(" ")
    """
    sourmash sig cat \
        ${sketch_list} \
        -o query_metagenomes.k31.scaled1000.zip
    """
}

process SOURMASH_COMPARE_QUERY_SKETCHES {

    tag "${metric}"
    label "medium"

    input:
    tuple val(metric), path(query_sketch_collection)

    output:
    tuple val(metric), path("query_sketch_similarity.${metric}.npy"), path("query_sketch_similarity.${metric}.csv"), emit: matrices

    script:
    def ignore_abundance = metric == "noabund" ? "--ignore-abundance" : ""
    """
    sourmash compare \
        ${query_sketch_collection} \
        ${ignore_abundance} \
        --output query_sketch_similarity.${metric}.npy \
        --csv query_sketch_similarity.${metric}.csv
    """
}
