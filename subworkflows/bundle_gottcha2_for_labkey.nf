workflow BUNDLE_GOTTCHA2_FOR_LABKEY {
    take:
    gottcha2_results           // queue channel: [sample_id, full_tsv, ref_mmi, stats, tax_tsv]
    gottcha2_extracted_fastas  // queue channel: [sample_id, fasta, full_tsv]
    experiment_id              // value channel: experiment ID
    run_id                     // value channel: workflow.runName
    validation_complete        // value channel: gate ensuring validation passed
    run_context                // value channel: [sample_set_id, state_dir] from CHECK_RUN_STATE

    main:

    LABKEY_UPLOAD_GOTTCHA2_FULL(
        gottcha2_results,
        experiment_id,
        run_id,
        run_context,
        validation_complete
    )

    LABKEY_UPLOAD_GOTTCHA2_FASTA(
        gottcha2_extracted_fastas,
        experiment_id,
        run_id,
        run_context,
        validation_complete
    )

    LABKEY_WEBDAV_UPLOAD_FILES(
        gottcha2_extracted_fastas,
        experiment_id,
        validation_complete
    )

    emit:
    upload_log = LABKEY_UPLOAD_GOTTCHA2_FULL.out.log
        .mix(LABKEY_UPLOAD_GOTTCHA2_FASTA.out.log)
}

process LABKEY_UPLOAD_GOTTCHA2_FULL {
    label 'low'

    tag "${sample_id}"

    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(full_tsv), path(ref_mmi), path(stats), path(tax_tsv)
    val experiment_id
    val run_id
    tuple val(sample_set_id), val(state_dir)  // run_context from CHECK_RUN_STATE
    val _validation_complete  // Gate: ensures validation passed

    output:
    path "fasta_labkey_upload.log", emit: log

    script:
    def sample_set_arg = sample_set_id ? "--sample-set-id '${sample_set_id}'" : ""
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    def run_id_arg = run_id ? "--run-id '${run_id}'" : ""
    """
    labkey_upload_gottcha2_full.py \
    --input-tsv ${full_tsv} \
    --sample ${sample_id} \
    --experiment '${experiment_id}' \
    ${sample_set_arg} \
    ${state_dir_arg} \
    ${run_id_arg} \
    --labkey-server '${params.labkey_server}' \
    --labkey-project-name '${params.labkey_project_name}' \
    --labkey-api-key \$LABKEY_API_KEY \
    --db-version ${params.gottcha2_db_version} \
    --labkey-schema '${params.labkey_schema}' \
    --table-name '${params.labkey_gottcha_full_list}' \
    """
}

process LABKEY_UPLOAD_GOTTCHA2_FASTA {

    label 'low'
    tag "${sample_id}"

    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(fasta), path(full_tsv)
    val experiment_id
    val run_id
    tuple val(sample_set_id), val(state_dir)  // run_context from CHECK_RUN_STATE
    val _validation_complete  // Gate: ensures validation passed

    output:
    path("${sample_id}_df.tsv"), emit: log

    script:
    def sample_set_arg = sample_set_id ? "--sample_set_id '${sample_set_id}'" : ""
    def state_dir_arg = state_dir ? "--state_dir '${state_dir}'" : ""
    """
    labkey_upload_gottcha2_fasta.py \
        --fasta ${fasta} \
        --sample_id ${sample_id} \
        --output_tsv ${sample_id}_df.tsv \
        --server ${params.labkey_server} \
        --experiment_id ${experiment_id} \
        --container ${params.labkey_project_name} \
        --list ${params.labkey_gottcha_fasta_list} \
        --api_key \$LABKEY_API_KEY \
        --batch_size 10000 \
        --notes "NVD2 upload" \
        --run_id '${run_id}' \
        ${sample_set_arg} \
        ${state_dir_arg} \
    """
}

// Process to upload Gottcha2 report and extracted fastas
process LABKEY_WEBDAV_UPLOAD_FILES {

    tag "${sample_id}"

    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(fasta), path(full_tsv)
    val experiment_id
    val _validation_complete  // Gate: ensures validation passed

    output:
    tuple val(sample_id), path("${fasta}.gz"), path("${full_tsv}.gz")

    script:
    """
    gzip -c ${fasta} > ${fasta}.gz
    gzip -c ${full_tsv} > ${full_tsv}.gz

    # Upload tabular Gottcha2 report
    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${full_tsv}.gz ${experiment_id}/${sample_id}/gottcha2/${full_tsv}.gz

    # Upload fasta file
    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${fasta}.gz ${experiment_id}/${sample_id}/gottcha2/${fasta}.gz
    """
}
