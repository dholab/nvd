workflow VALIDATE_LK_GOTTCHA2 {

    VALIDATE_GOTTCHA2_FULL_LIST()
    VALIDATE_GOTTCHA2_BLAST_VERIFIED_LIST()
    VALIDATE_GOTTCHA2_FASTA_LIST()

}

/*
This workflow is basically a duplicate of validate_lk_blast_lists.nf.
Refer to that workflow for more notes.
*/

process VALIDATE_GOTTCHA2_FULL_LIST {
    label 'low'
    secret 'nvd2'

    output:
    path 'gottcha2_full_validation_report.txt'

    when:
    params.tools && (params.tools.contains("gottcha") || params.tools.contains("all"))

    script:
    """
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_gottcha_full_list}' \
        --api_key \$nvd2 \
        --experiment_id ${params.experiment_id} \
        --type gottcha2_full > gottcha2_full_validation_report.txt 2>&1
    """
}

// We can reuse the schema and everything from VALIDATE_GOTTCHA2_FULL_LIST as the blast verified list
// is just a duplicate of the full list with less hits uploaded after filtering. We just need to ensure
// that the correct different list is passed through.
process VALIDATE_GOTTCHA2_BLAST_VERIFIED_LIST {
    label 'low'
    secret 'nvd2'

    output:
    path 'gottcha2_blast_verified_validation_report.txt'

    when:
    params.tools && (params.tools.contains("gottcha") || params.tools.contains("all"))

    script:
    """
    validate_labkey.py \
        --server '${params.labkey_server}' \
        --container '${params.labkey_project_name}' \
        --list '${params.labkey_gottcha_blast_verified_full_list}' \
        --api_key \$nvd2 \
        --experiment_id ${params.experiment_id} \
        --type gottcha2_full > gottcha2_blast_verified_validation_report.txt 2>&1
    """
}

process VALIDATE_GOTTCHA2_FASTA_LIST {
    label 'low'
    secret 'nvd2'

    output:
    path 'gottcha2_fasta_validation_report.txt'

    when:
    params.tools && (params.tools.contains("gottcha") || params.tools.contains("all"))

    script:
    """
    validate_labkey.py \\
        --server '${params.labkey_server}' \\
        --container '${params.labkey_project_name}' \\
        --list '${params.labkey_gottcha_fasta_list}' \\
        --api_key \$nvd2 \\
        --experiment_id ${params.experiment_id} \\
        --type gottcha2_fasta > gottcha2_fasta_validation_report.txt 2>&1
    """
}
