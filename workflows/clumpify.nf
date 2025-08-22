include { CLUMP_READS } from "../modules/bbmap"
include { SCRUB_HUMAN_READS } from "../modules/stat"

workflow CLUMPIFY_WORKFLOW {

    take:
    ch_gathered_reads
    ch_start_gate
        
    main:
    def gated_reads = ch_gathered_reads
        .combine(ch_start_gate)   // pairs (read_tuple, token)
        .map { read_tuple, _ -> read_tuple }

    // Add human read scrubbing step for public posting of clumpified data to SRA
    CLUMP_READS(
        gated_reads.map { id, _platform, reads -> tuple(id, file(reads)) }
    )

    // Add human read scrubbing if specified
    SCRUB_HUMAN_READS(
        CLUMP_READS.out
    )
}
