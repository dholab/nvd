include { SOURMASH_COLLECT_QUERY_SKETCHES ; SOURMASH_COMPARE_QUERY_SKETCHES } from "../modules/sourmash"
include { COMPUTE_SAMPLE_SKETCH_DISTANCES ; COMPUTE_SAMPLE_ORDINATION ; REPORT_POSSIBLE_SAMPLE_MIXUPS ; PLOT_SAMPLE_ORDINATION } from "../modules/sample_similarity"

workflow SAMPLE_SIMILARITY_QC {
    take:
    ch_query_sketches  // tuple(sample_id, platform, read_structure, sourmash_signature)

    main:
    ch_query_sketch_files = ch_query_sketches
        .map { _sample_id, _platform, _read_structure, signature -> signature }
        .collect()

    SOURMASH_COLLECT_QUERY_SKETCHES(ch_query_sketch_files)

    ch_compare_inputs = Channel.of("abund", "noabund")
        .combine(SOURMASH_COLLECT_QUERY_SKETCHES.out.collection)

    SOURMASH_COMPARE_QUERY_SKETCHES(ch_compare_inputs)

    COMPUTE_SAMPLE_SKETCH_DISTANCES(SOURMASH_COMPARE_QUERY_SKETCHES.out.matrices)

    COMPUTE_SAMPLE_ORDINATION(COMPUTE_SAMPLE_SKETCH_DISTANCES.out.pairwise_distances)

    REPORT_POSSIBLE_SAMPLE_MIXUPS(COMPUTE_SAMPLE_SKETCH_DISTANCES.out.pairwise_distances)

    ch_plot_inputs = COMPUTE_SAMPLE_ORDINATION.out.ordination
        .map { metric, _distance_matrix, ordination, variance -> tuple(metric, ordination, variance) }

    PLOT_SAMPLE_ORDINATION(ch_plot_inputs)

    emit:
    query_sketch_collection = SOURMASH_COLLECT_QUERY_SKETCHES.out.collection
    similarity_matrices = SOURMASH_COMPARE_QUERY_SKETCHES.out.matrices
    pairwise_distances = COMPUTE_SAMPLE_SKETCH_DISTANCES.out.pairwise_distances
    ordination = COMPUTE_SAMPLE_ORDINATION.out.ordination
    possible_mixups = REPORT_POSSIBLE_SAMPLE_MIXUPS.out.reports
    plots = PLOT_SAMPLE_ORDINATION.out.plots
}
