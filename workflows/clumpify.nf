include { CLUMP_READS } from "../modules/bbmap"

workflow CLUMPIFY_WORKFLOW {

    take:
    ch_gathered_reads
    _completion_tokens
        
    main:
    CLUMP_READS(
        ch_gathered_reads.map { id, _platform, reads -> tuple(id, file(reads)) }
    )

}
