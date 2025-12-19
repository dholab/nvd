include { CLUMP_READS       } from "../modules/bbmap"
include { SCRUB_HUMAN_READS } from "../modules/stat"

workflow CLUMPIFY_WORKFLOW {
    take:
    ch_gathered_reads
    ch_start_gate

    main:
    def gated_reads = ch_gathered_reads
        .combine(ch_start_gate)
        .map { it[0..-2] }

    // Add human read scrubbing step for public posting of clumpified data to SRA
    CLUMP_READS(
        gated_reads.map { id, _platform, reads -> tuple(id, file(reads)) }
    )

    // Check if sra_human_db is a valid file path
    if (params.sra_human_db != null && !file(params.sra_human_db).isFile()) {
        error("Error: Human database file does not exist: ${params.sra_human_db}")
    }

    ch_human_reads = params.sra_human_db ? Channel.fromPath( params.sra_human_db ) : Channel.empty()

    SCRUB_HUMAN_READS(
        CLUMP_READS.out.combine(ch_human_reads)
    )
}
