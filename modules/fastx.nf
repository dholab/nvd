process PROFILE_FASTX {

    tag { meta.query_class ? "${meta.id}, ${meta.query_class}" : "${meta.id}" }
    label "low"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastx_profile.json"), path("*.length_histogram.tsv"), path("*.sequence_count.txt"), emit: profiled
    tuple val(meta), path("*.quality_histogram.tsv"), optional: true, emit: quality_histogram
    tuple val(meta), path("publish/*"), optional: true, emit: public_artifacts

    script:
    def output_prefix_arg = meta.profile_prefix ? "--output-prefix '${meta.profile_prefix}'" : ""
    def threshold_args = meta.thresholds
        ? meta.thresholds.collect { threshold -> "--threshold '${threshold.name}:${threshold.axis}:${threshold.value}'" }.join(" ")
        : ""
    """
    profile_fastx.py \
        --sample-id '${meta.id}' \
        --stage '${meta.profile_stage}' \
        ${output_prefix_arg} \
        ${threshold_args} \
        --input ${reads}

    sequence_count=\$(tr -d '[:space:]' < *.sequence_count.txt)
    if [ "\${sequence_count}" -gt 0 ]; then
        mkdir -p publish
        ln -s ../${reads} publish/\$(basename ${reads})
        for artifact in *.fastx_profile.json *.length_histogram.tsv *.quality_histogram.tsv; do
            [ -e "\${artifact}" ] || continue
            ln -s ../"\${artifact}" publish/"\${artifact}"
        done
    fi
    """
}

process PLOT_FASTX_LENGTH_PROFILE {

    tag { meta.query_class ? "${meta.id}, ${meta.query_class}" : "${meta.id}" }
    label "low"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(profile_json), path(length_histogram)

    output:
    tuple val(meta), path("*.length_histogram.png"), emit: plot

    script:
    """
    plot_fastx_profile.py \
        --profile-json ${profile_json} \
        --length-histogram ${length_histogram}
    """
}

process PLOT_FASTX_QUALITY_PROFILE {

    tag { meta.query_class ? "${meta.id}, ${meta.query_class}" : "${meta.id}" }
    label "low"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(profile_json), path(quality_histogram)

    output:
    tuple val(meta), path("*.quality_histogram.png"), emit: plot

    script:
    """
    plot_fastx_profile.py \
        --profile-json ${profile_json} \
        --quality-histogram ${quality_histogram}
    """
}
