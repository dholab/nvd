process RETRIEVE_GETTAX {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2
    
    output:
    path "gettax.sqlite"

    script:
    """
    retrieve_gettax.py
    """
}
