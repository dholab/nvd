process RETRIEVE_GETTAX {

    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2
    
    output:
    path "gettax.sqlite"

    when:
    params.tools && (params.tools.contains("nvd") || params.tools.contains("all"))

    script:
    """
    retrieve_gettax.py
    """
}
