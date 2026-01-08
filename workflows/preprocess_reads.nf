include { DEDUP_WITH_CLUMPIFY ; TRIM_ADAPTERS ; FILTER_READS ; REPAIR_PAIRS } from "../modules/bbmap"
include { SCRUB_HOST_READS } from "../modules/hostile"

workflow PREPROCESS_READS {
    take:
    ch_fastq_tuple  // tuple(sample_id, platform, read_structure, reads)

    main:
    // Resolve optional step flags: explicit param wins, otherwise fall back to master switch
    def should_dedup = params.dedup ?: params.preprocess
    def should_trim = params.trim_adapters ?: params.preprocess
    def should_scrub = params.scrub_host_reads ?: params.preprocess
    def should_filter = params.filter_reads ?: params.preprocess

    // 1. Dedup
    ch_after_dedup = should_dedup
        ? DEDUP_WITH_CLUMPIFY(ch_fastq_tuple)
        : ch_fastq_tuple

    // 2. Adapter trim (Illumina only)
    ch_branched_for_trim = ch_after_dedup.branch { _id, platform, _read_structure, _reads ->
        illumina: platform == "illumina"
        other: true
    }

    ch_trimmed_illumina = should_trim
        ? TRIM_ADAPTERS(ch_branched_for_trim.illumina)
        : ch_branched_for_trim.illumina

    ch_after_trim = ch_trimmed_illumina.mix(ch_branched_for_trim.other)

    // 3. Host scrub
    ch_after_scrub = should_scrub
        ? SCRUB_HOST_READS(ch_after_trim)
        : ch_after_trim

    // 4. Quality/length filter (with platform-specific quality threshold)
    ch_with_quality_threshold = ch_after_scrub.map { sample_id, platform, read_structure, reads ->
        def min_qual = platform == "illumina"
            ? params.min_read_quality_illumina
            : params.min_read_quality_nanopore
        tuple(sample_id, platform, read_structure, reads, min_qual)
    }

    ch_after_filter = should_filter
        ? FILTER_READS(ch_with_quality_threshold)
        : ch_after_scrub

    // 5. Repair pairs (interleaved reads only) - fixes orphans from upstream steps
    ch_branched_for_repair = ch_after_filter.branch { _id, _platform, read_structure, _reads ->
        interleaved: read_structure == "interleaved"
        other: true
    }

    ch_repaired = REPAIR_PAIRS(ch_branched_for_repair.interleaved)
    ch_preprocessed = ch_repaired.mix(ch_branched_for_repair.other)

    emit:
    ch_preprocessed
}
