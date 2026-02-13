workflow VALIDATE_LK_GOTTCHA2 {

    VALIDATE_GOTTCHA2_FULL_LIST()
    VALIDATE_GOTTCHA2_FASTA_LIST()

}

/*
This workflow is basically a duplicate of validate_lk_blast_lists.nf.
Refer to that workflow for more notes.
*/

process VALIDATE_GOTTCHA2_FULL_LIST {
    label 'low'
    secret 'LABKEY_API_KEY'

    output:
    path 'gottcha2_full_validation_report.txt'

    script:
    """
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_gottcha_full_list}' \
        --api_key \$LABKEY_API_KEY \
        --experiment_id ${params.experiment_id} \
        --type gottcha2_full > gottcha2_full_validation_report.txt 2>&1
    """
}

process VALIDATE_GOTTCHA2_FASTA_LIST {
    label 'low'
    secret 'LABKEY_API_KEY'

    output:
    path 'gottcha2_fasta_validation_report.txt'

    script:
    """
    validate_labkey.py \\
        --server '${params.labkey_server}' \\
        --container '${params.labkey_project_name}' \\
        --list '${params.labkey_gottcha_fasta_list}' \\
        --api_key \$LABKEY_API_KEY \\
        --experiment_id ${params.experiment_id} \\
        --type gottcha2_fasta > gottcha2_fasta_validation_report.txt 2>&1
    """
}
