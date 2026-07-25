process NOTIFY_RAPID_SCREENING_SLACK {
    tag "${sample_id}"
    label "low"
    maxForks 1

    errorStrategy 'ignore'

    secret 'SLACK_BOT_TOKEN'

    input:
    tuple val(sample_id), path(risk_summary), path(enrichment_summary)
    val target_enrichment_enabled
    val run_id

    output:
    tuple val(sample_id), val(true), emit: done

    when:
    params.experimental && params.slack_enabled && params.slack_channel

    script:
    def enrichment_arg = target_enrichment_enabled ? "--target-enrichment-enabled" : ""
    """
    notify_rapid_screening_slack.py \
        --sample-id '${sample_id}' \
        --run-id '${run_id}' \
        --channel '${params.slack_channel}' \
        --risk-summary ${risk_summary} \
        --enrichment ${enrichment_summary} \
        --minimum-risk-group RG3 \
        ${enrichment_arg} \
        -v
    """
}


process NOTIFY_RUN_COMPLETION_SLACK {
    tag "${run_id}"
    label "low"

    errorStrategy 'ignore'

    secret 'SLACK_BOT_TOKEN'

    input:
    path experiment_results
    val sample_set_id
    val run_id

    output:
    val true, emit: done

    when:
    params.slack_enabled && params.slack_channel

    script:
    def encoded_experiment_id = params.experiment_id != null
        ? params.experiment_id.toString().bytes.encodeBase64().toString()
        : null
    def experiment_arg = encoded_experiment_id
        ? "--experiment-id-base64 ${encoded_experiment_id}"
        : ""
    def labkey_config = groovy.json.JsonOutput.toJson([
        server: params.labkey_server,
        project: params.labkey_project_name,
        results_list: params.labkey_blast_meta_hits_list,
    ])
    def labkey_args = params.labkey
        ? "--labkey-config-base64 ${labkey_config.bytes.encodeBase64().toString()}"
        : ""
    """
    notify_slack.py \
        --run-id '${run_id}' \
        --channel '${params.slack_channel}' \
        --sample-set-id '${sample_set_id}' \
        --experiment-results ${experiment_results} \
        -v ${experiment_arg} ${labkey_args}
    """
}
