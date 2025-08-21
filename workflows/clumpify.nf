include { CLUMP_READS } from "../modules/bbmap"
include { SCRUB_HUMAN_READS } from "../modules/stat"

workflow CLUMPIFY_WORKFLOW {

    take:
    ch_gathered_reads
    _completion_tokens
        
    main:

    // Add human read scrubbing step for public posting of clumpified data to SRA

    CLUMP_READS(
        ch_gathered_reads.map { id, _platform, reads -> tuple(id, file(reads)) }
    )

    // Add human read scrubbing if specified
    if (params.human_read_scrub) {
            SCRUB_HUMAN_READS(
                CLUMP_READS.out
                )
            }
}