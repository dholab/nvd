process BUILD_NVD_MULTIQC_INPUTS {
    label "low"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path resolved_reads
    path nvd_version
    path "fastqc_packages/*"
    val experimental_enabled

    output:
    path "nvd_inputs", emit: inputs

    script:
    """
    mkdir -p fastqc_packages
    build_nvd_multiqc_inputs.py \
        --resolved-reads '${resolved_reads}' \
        --nvd-version '${nvd_version}' \
        --fastqc-root fastqc_packages \
        --experimental-enabled '${experimental_enabled}' \
        --output-dir nvd_inputs
    """
}

process MULTIQC_NVD {
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
