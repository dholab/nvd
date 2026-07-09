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
