process BUILD_MULTIQC_INPUTS {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path resolved_reads
    path nvd_version
    path "fastqc_packages/*"
    path "target_enrichment_packages/*"
    path "depletion_packages/*"
    path "fastx_packages/*"
    path "assembly_packages/*"
    path "query_preparation_packages/*"
    path "blast_packages/*"
    path "taxonomy_packages/*"
    val experimental_enabled
    val target_enrichment_enabled
    val depletion_enabled
    val assembly_enabled
    val blast_enabled

    output:
    path "nvd_inputs", emit: inputs

    script:
    """
    mkdir -p fastqc_packages
    mkdir -p target_enrichment_packages depletion_packages fastx_packages assembly_packages query_preparation_packages blast_packages taxonomy_packages
    build_nvd_multiqc_inputs.py \
        --resolved-reads '${resolved_reads}' \
        --nvd-version '${nvd_version}' \
        --fastqc-root fastqc_packages \
        --target-enrichment-root target_enrichment_packages \
        --depletion-root depletion_packages \
        --fastx-root fastx_packages \
        --assembly-root assembly_packages \
        --query-preparation-root query_preparation_packages \
        --blast-root blast_packages \
        --taxonomy-root taxonomy_packages \
        --experimental-enabled '${experimental_enabled}' \
        --target-enrichment-enabled '${target_enrichment_enabled}' \
        --depletion-enabled '${depletion_enabled}' \
        --assembly-enabled '${assembly_enabled}' \
        --blast-enabled '${blast_enabled}' \
        --output-dir nvd_inputs
    """
}

process PACKAGE_NVD_DEACON_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(domain), val(sample_id), val(query_class), path(stats, stageAs: "stats.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py deacon \
        --domain='${domain}' \
        --sample-id='${sample_id}' \
        --query-class='${query_class}' \
        --stats=stats.input
    """
}

process PACKAGE_NVD_PROFILE_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(domain), val(sample_id), val(stage), val(query_class), val(producer), path(profile, stageAs: "profile.input"), path(length_histogram, stageAs: "length_histogram.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py profile \
        --domain='${domain}' \
        --sample-id='${sample_id}' \
        --stage='${stage}' \
        --query-class='${query_class}' \
        --producer='${producer}' \
        --profile=profile.input \
        --length-histogram=length_histogram.input
    """
}

process PACKAGE_NVD_PROFILE_WITH_QUALITY_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(domain), val(sample_id), val(stage), val(query_class), val(producer), path(profile, stageAs: "profile.input"), path(length_histogram, stageAs: "length_histogram.input"), path(quality_histogram, stageAs: "quality_histogram.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py profile \
        --domain='${domain}' \
        --sample-id='${sample_id}' \
        --stage='${stage}' \
        --query-class='${query_class}' \
        --producer='${producer}' \
        --profile=profile.input \
        --length-histogram=length_histogram.input \
        --quality-histogram=quality_histogram.input
    """
}

process PACKAGE_NVD_ELIGIBILITY_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(producer), val(decision), val(sequence_count), val(threshold), val(reason)

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py eligibility \
        --sample-id='${sample_id}' \
        --producer='${producer}' \
        --decision='${decision}' \
        --sequence-count='${sequence_count}' \
        --threshold='${threshold}' \
        --reason='${reason}'
    """
}

process PACKAGE_NVD_LONG_READ_ELIGIBILITY_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(summary, stageAs: "summary.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py long-read-eligibility \
        --sample-id='${sample_id}' \
        --summary=summary.input
    """
}

process PACKAGE_NVD_UNION_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(summary, stageAs: "summary.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py union \
        --sample-id='${sample_id}' \
        --summary=summary.input
    """
}

process PACKAGE_NVD_PREPARED_QUERY_BATCHES_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(summary, stageAs: "summary.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py prepared-query-batches \
        --sample-id='${sample_id}' \
        --summary=summary.input
    """
}

process PACKAGE_NVD_MEGABLAST_QUERY_PARTITION_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), val(query_class), path(summary, stageAs: "summary.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py megablast-query-partition \
        --sample-id='${sample_id}' \
        --query-class='${query_class}' \
        --summary=summary.input
    """
}

process PACKAGE_NVD_TAXON_BIG_TABLE_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(sample_id), path(table, stageAs: "table.input")

    output:
    path "*.pkg", arity: '1', emit: packages

    script:
    """
    package_nvd_report.py taxon-big-table \
        --sample-id='${sample_id}' \
        --table=table.input
    """
}

process GENERATE_MULTIQC_REPORT {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path "fastqc_zips/*"
    path nvd_inputs
    path multiqc_config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    mkdir -p fastqc_zips
    multiqc nvd_inputs fastqc_zips --config '${multiqc_config}' --outdir . --filename multiqc_report.html --force
    mkdir -p multiqc_data/nvd_inputs
    if [ -d multiqc_report_data ]; then
        cp -R multiqc_report_data/. multiqc_data/
    fi
    cp -R nvd_inputs/. multiqc_data/nvd_inputs/
    cp '${multiqc_config}' multiqc_data/nvd_inputs/multiqc_config.yaml
    """
}
