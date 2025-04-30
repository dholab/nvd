process RETRIEVE_GETTAX {
    
    output:
    path "gettax.sqlite"

    script:
    """
    retrieve_gettax.py
    """
}
