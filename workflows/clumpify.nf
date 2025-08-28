include { CLUMP_READS       } from "../modules/bbmap"
include { SCRUB_HUMAN_READS } from "../modules/stat"

workflow CLUMPIFY_WORKFLOW {
    take:
    ch_gathered_reads
    ch_start_gate

    main:
    def gated_reads = ch_gathered_reads
        .combine(ch_start_gate)
        .map { read_tuple, _rest -> read_tuple }

    // Add human read scrubbing step for public posting of clumpified data to SRA
    CLUMP_READS(
        gated_reads.map { id, _platform, reads -> tuple(id, file(reads)) }
    )

    // Check if human_read_scrub is a valid file path
    human_db_file = file(params.human_read_scrub)

    if (!human_db_file.exists()) {
        error("Error: Human database file does not exist: ${params.human_read_scrub}")
    }

    if (!human_db_file.isFile()) {
        error("Error: Human database path is not a file: ${params.human_read_scrub}")
    }

    SCRUB_HUMAN_READS(
        CLUMP_READS.out
    )
}
