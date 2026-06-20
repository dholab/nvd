process COMPUTE_SAMPLE_SKETCH_DISTANCES {

    tag "${metric}"
    label "low"

    input:
    tuple val(metric), path(_matrix_npy), path(matrix_csv)

    output:
    tuple val(metric), path("sample_sketch_distances.${metric}.tsv"), emit: pairwise_distances

    script:
    """
    sample_sketch_distances.py \
        --metric ${metric} \
        --matrix ${matrix_csv} \
        --output sample_sketch_distances.${metric}.tsv
    """
}

process COMPUTE_SAMPLE_ORDINATION {

    tag "${metric}"
    label "low"

    input:
    tuple val(metric), path(pairwise_tsv)

    output:
    tuple val(metric), path("sample_distance_matrix.${metric}.tsv"), path("sample_ordination.${metric}.tsv"), path("sample_ordination_variance.${metric}.tsv"), emit: ordination

    script:
    """
    sample_ordination.py \
        --metric ${metric} \
        --pairwise ${pairwise_tsv} \
        --distance-output sample_distance_matrix.${metric}.tsv \
        --coordinates-output sample_ordination.${metric}.tsv \
        --variance-output sample_ordination_variance.${metric}.tsv
    """
}

process REPORT_POSSIBLE_SAMPLE_MIXUPS {

    tag "${metric}"
    label "low"

    input:
    tuple val(metric), path(pairwise_tsv)

    output:
    tuple val(metric), path("nearest_sample_neighbors.${metric}.tsv"), path("possible_sample_mixup_candidates.${metric}.tsv"), emit: reports

    script:
    """
    report_possible_sample_mixups.py \
        --metric ${metric} \
        --pairwise ${pairwise_tsv} \
        --nearest-neighbors-output nearest_sample_neighbors.${metric}.tsv \
        --candidate-output possible_sample_mixup_candidates.${metric}.tsv
    """
}

process SUMMARIZE_SAMPLE_SIMILARITY_CANDIDATE_EVIDENCE {

    tag "candidate_evidence"
    label "low"

    input:
    path(candidate_reports)

    output:
    path("sample_similarity_candidate_evidence.tsv"), emit: evidence

    script:
    def candidate_list = candidate_reports.collect { report -> "\"${report}\"" }.join(" ")
    """
    sample_similarity_candidate_evidence.py \
        --candidates ${candidate_list} \
        --output sample_similarity_candidate_evidence.tsv
    """
}

process PLOT_SAMPLE_ORDINATION {

    tag "${metric}"
    label "low"

    input:
    tuple val(metric), path(ordination_tsv), path(variance_tsv)

    output:
    tuple val(metric), path("sample_ordination.${metric}.html"), path("sample_ordination.${metric}.png"), emit: plots

    script:
    """
    plot_sample_ordination.py \
        --metric ${metric} \
        --coordinates ${ordination_tsv} \
        --variance ${variance_tsv} \
        --html-output sample_ordination.${metric}.html \
        --png-output sample_ordination.${metric}.png
    """
}
