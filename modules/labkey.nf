// Process adds an additional safeguard to large NVD2 runs where the experiment_id
// and htcondor cluster are used to confirm a new NVD2 run is not accidentally
// uploading results to labkey for samples from the wrong experiment.
process VALIDATE_LK_EXP_TO_CLUSTER {

    secret 'nvd2'

    tag "exp_${params.experiment_id}_cluster_${params.condor_cluster}"

    // Option 1: Use a dummy input channel if you need to trigger from workflow
    // input:
    // val ready from Channel.value(true)

    // Option 2: No input at all - process runs immediately when called
    input:
    // Leave empty for no inputs

    output:
    path("lk_validation_log.txt"), emit: log
    stdout emit: validation_result

    when:
    params.condor_cluster && params.labkey && params.experiment_id

    script:
    """
    echo "Starting guard list validation for experiment ${params.experiment_id}"
    echo "Cluster ID: ${params.condor_cluster}"
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
        --cluster_id ${params.condor_cluster} \\
        --experiment_id ${params.experiment_id} \\
        2>&1 | tee lk_validation_log.txt

    # Check if the Python script succeeded
    if [ \${PIPESTATUS[0]} -eq 0 ]; then
        echo "✅ Guard list validation PASSED"
        echo "VALIDATION_SUCCESS"
    else
        echo "❌ Guard list validation FAILED"
        echo "VALIDATION_FAILED"
        exit 1
    fi
    """
}