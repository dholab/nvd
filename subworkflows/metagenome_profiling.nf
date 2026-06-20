include { SOURMASH_SKETCH_QUERY_METAGENOME } from "../modules/sourmash"

workflow METAGENOME_PROFILING {
    take:
    ch_preprocessed_reads  // tuple(sample_id, platform, read_structure, fastq)

    main:
    SOURMASH_SKETCH_QUERY_METAGENOME(ch_preprocessed_reads)

    emit:
    query_sketches = SOURMASH_SKETCH_QUERY_METAGENOME.out.query_sketches
}
