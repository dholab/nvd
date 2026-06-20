include { SOURMASH_SKETCH_QUERY_METAGENOME } from "../modules/sourmash"

workflow METAGENOME_PROFILING {
    take:
    ch_preprocessed_reads  // tuple(sample_id, platform, read_structure, fastq)

    main:
    // Avoid scheduling sourmash for read files that cannot produce a query
    // sketch. Empty gzip FASTQ outputs are not zero-byte files, so this is a
    // cheap scheduling guard rather than a FASTQ record-count predicate.
    ch_sketchable_reads = ch_preprocessed_reads.filter { _sample_id, _platform, _read_structure, reads ->
        def read_file = file(reads)
        read_file.name.endsWith(".gz") ? read_file.size() > 28 : read_file.size() > 0
    }

    SOURMASH_SKETCH_QUERY_METAGENOME(ch_sketchable_reads)

    emit:
    query_sketches = SOURMASH_SKETCH_QUERY_METAGENOME.out.query_sketches
}
