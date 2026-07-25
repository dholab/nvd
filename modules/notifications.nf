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
