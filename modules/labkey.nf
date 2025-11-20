// Process adds an additional safeguard to large NVD2 runs where the experiment_id
// and htcondor cluster are used to confirm a new NVD2 run is not accidentally
// uploading results to labkey for samples from the wrong experiment.
process VALIDATE_LK_EXP_TO_CLUSTER {

    secret 'nvd2'

    tag "exp_${params.experiment_id}" // FIXME: params.condor_cluster needs to be renamed and have a default
    label "low"

    output:
    path "lk_validation_log.txt", emit: log
    stdout emit: validation_result

    when:
    params.labkey && params.experiment_id

    script:

    def uuid_to_use = params.uuid ?: workflow.sessionId
    """
    echo "Starting guard list validation for experiment ${params.experiment_id}"
    echo "UUID: ${uuid_to_use}"
    echo "Guard List: ${params.labkey_exp_id_guard_list}"
    echo "Server: ${params.labkey_server}"
    echo "Container: ${params.labkey_project_name}"
    echo "Timestamp: \$(date)"
    echo "----------------------------------------"

    labkey_check_guard_list.py \\
        --server ${params.labkey_server} \\
        --container ${params.labkey_project_name} \\
        --guard_list ${params.labkey_exp_id_guard_list} \\
        --api_key \$nvd2 \\
        --unique_id ${uuid_to_use} \\
        --experiment_id ${params.experiment_id} \\
        2>&1 | tee lk_validation_log.txt

    # Check if the Python script succeeded
    if [ \${PIPESTATUS[0]} -eq 0 ]; then
        echo "✅ Guard list validation PASSED"
        echo "VALIDATION_SUCCESS"
    else
        echo "❌ Guard list validation FAILED"
        echo "VALIDATION_FAILED"
        echo "Does your experiment ID appear in the guard list?"
        echo "If so, either delete the experiment_id row from the guard list or set the uuid parameter to match the previous workflow.sessionId or uuid used for that experiment."
        exit 1
    fi
    """
}
