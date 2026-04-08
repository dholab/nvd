include { CLUMP_READS          } from "../modules/bbmap"
include { DEACON_FETCH_INDEX   } from "../modules/deacon"
include { DEACON_DEPLETE       } from "../modules/deacon"

workflow CLUMPIFY_WORKFLOW {
    take:
    ch_gathered_reads  // tuple(sample_id, platform, read_structure, fastq)
    ch_start_gate

    main:
    def gated_reads = ch_gathered_reads
        .combine(ch_start_gate)
        .map { it[0..-2] }

    // Reorder reads by kmer similarity for better compression.
    // Full 4-element tuple passes through so DEACON_DEPLETE can consume it directly.
    CLUMP_READS(gated_reads)

    // --- Host read depletion with deacon ---
    // Resolve deacon index: local path takes priority over URL download
    ch_local_index = params.deacon_index
        ? Channel.fromPath(params.deacon_index)
        : Channel.empty()

    ch_fetch_url = (!params.deacon_index && params.deacon_index_url)
        ? Channel.of(params.deacon_index_url)
        : Channel.empty()

    DEACON_FETCH_INDEX(ch_fetch_url)

    ch_index = ch_local_index.mix(DEACON_FETCH_INDEX.out.index)

    // DEACON_DEPLETE expects tuple(sample_id, platform, read_structure, reads, index)
    DEACON_DEPLETE(
        CLUMP_READS.out.combine(ch_index)
    )
}
