include { CLUMP_READS } from "../modules/bbmap"
include { SCRUB_HUMAN_READS } from "../modules/stat"

workflow CLUMPIFY_WORKFLOW {

    take:
    ch_gathered_reads
    _completion_tokens
        
    main:

    // Add human read scrubbing step for public posting of clumpified data to SRA

    SCRUB_HUMAN_READS(
        ch_gathered_reads.map { id, _platform, reads -> tuple(id, file(reads)) }
    )

    CLUMP_READS(
        SCRUB_HUMAN_READS.out
    )

}
