workflow BUNDLE_GOTTCHA2_FOR_LABKEY {
    take:
    gottcha2_results     // channel: [ meta, csv ]
    gottcha2_extracted_fastas  // channel: [ meta, fasta ]
    // Note: sample_set_id is not currently passed from gottcha2_workflow.nf
    // The upload processes support it but it won't be used until CHECK_RUN_STATE
    // is added to the GOTTCHA2 workflow

    main:

    // DSL2 automatically forks channels when used by multiple consumers.
    // No explicit tap() needed - just use the channels directly.

    LABKEY_UPLOAD_GOTTCHA2_FULL(gottcha2_results)

    LABKEY_UPLOAD_GOTTCHA2_FASTA(gottcha2_extracted_fastas)

    LABKEY_WEBDAV_UPLOAD_FILES(gottcha2_extracted_fastas)
}

process LABKEY_UPLOAD_GOTTCHA2_FULL {
    label 'low'

    tag "${sample_id}"

    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(full_tsv), path(ref_mmi), path(stats), path(tax_tsv)

    output:
    path "fasta_labkey_upload.log", emit: log

    script:
    """
    labkey_upload_gottcha2_full.py \
    --input-tsv ${full_tsv} \
    --sample ${sample_id} \
    --experiment '${params.experiment_id}' \
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

    output:
    path("${sample_id}_df.tsv")

    script:
    """
    labkey_upload_gottcha2_fasta.py \
        --fasta ${fasta} \
        --sample_id ${sample_id} \
        --output_tsv ${sample_id}_df.tsv \
        --server ${params.labkey_server} \
        --experiment_id  ${params.experiment_id} \
        --container ${params.labkey_project_name} \
        --list ${params.labkey_gottcha_fasta_list} \
        --api_key \$LABKEY_API_KEY \
        --batch_size 10000 \
        --notes "NVD2 upload" \
        --run_id ${params.condor_cluster} \
    """
}

// Process to upload Gottcha2 report and extracted fastas
process LABKEY_WEBDAV_UPLOAD_FILES {

    tag "${sample_id}"

    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(fasta), path(full_tsv)

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
        upload ${full_tsv}.gz ${params.experiment_id}/${sample_id}/gottcha2/${full_tsv}.gz

    # Upload fasta file
    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${fasta}.gz ${params.experiment_id}/${sample_id}/gottcha2/${fasta}.gz
    """
}
