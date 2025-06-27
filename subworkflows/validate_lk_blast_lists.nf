workflow VALIDATE_LK_BLAST {
    VALIDATE_BLAST_HITS_LIST()
    VALIDATE_BLAST_FASTA_LIST()
}

process VALIDATE_BLAST_HITS_LIST {
    label 'process_low'
    secret 'nvd2'

    output:
    file 'blast_validation_report.txt'

    script:
    """
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_blast_meta_hits_list}' \
        --api_key \$nvd2 \
        --type blast > blast_validation_report.txt 2>&1
    """
}

process VALIDATE_BLAST_FASTA_LIST {
    label 'process_low'
    secret 'nvd2'

    output:
    file 'fasta_validation_report.txt'

    script:
    """
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_blast_fasta_list}' \
        --api_key \$nvd2 \
        --type fasta > fasta_validation_report.txt 2>&1
    """
}
