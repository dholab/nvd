// Pre-check: Ensure experiment ID hasn't been uploaded before
process VALIDATE_LK_EXP_FRESH {

    cache false
    secret 'nvd2'
    tag "exp_${params.experiment_id}"
    label "low"

    input:
    val trigger  // Any value to trigger execution

    output:
    val true, emit: validated

    when:
    params.labkey && params.experiment_id

    script:
    """
    labkey_check_guard_list.py \\
        --mode check \\
        --server ${params.labkey_server} \\
        --container ${params.labkey_project_name} \\
        --guard_list ${params.labkey_exp_id_guard_list} \\
        --api_key \$nvd2 \\
        --experiment_id ${params.experiment_id}
    """
}

// Post-upload: Register experiment ID to prevent future duplicates
process REGISTER_LK_EXPERIMENT {

    cache false
    secret 'nvd2'
    tag "exp_${params.experiment_id}"
    label 'low'

    input:
    val upload_complete  // Depends on uploads finishing

    output:
    val true, emit: registered

    when:
    params.labkey && params.experiment_id

    script:
    """
    labkey_check_guard_list.py \\
        --mode register \\
        --server ${params.labkey_server} \\
        --container ${params.labkey_project_name} \\
        --guard_list ${params.labkey_exp_id_guard_list} \\
        --api_key \$nvd2 \\
        --experiment_id ${params.experiment_id}
    """
}
