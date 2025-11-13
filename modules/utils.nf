process RETRIEVE_GETTAX {

    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2
    
    output:
    path "gettax.sqlite"

    when:
    params.tools && (params.tools.contains("nvd") || params.tools.contains("stat") || params.tools.contains("blast") || params.tools.contains("stat_blast") || params.tools.contains("stast") || params.tools.contains("all"))

    script:
    """
    retrieve_gettax.py
    """
}

process ANNOTATE_LEAST_COMMON_ANCESTORS {

    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2
    
    input:
    tuple val(sample_id), path(all_blast_hits)

    output:
    tuple val(sample_id), path("${sample_id}_blast.merged_with_lca.tsv")

    script:
    """
    annotate_blast_lca.py -i ${all_blast_hits} -o ${sample_id}_blast.merged_with_lca.tsv
    """

}
