include { DEACON_FILTER_CONTIGS } from "../modules/deacon"

workflow SCREEN_CONTIGS {
    take:
    ch_contigs          // tuple(sample_id, platform, read_structure, fasta, lookup) from PROCESS_CONTIGS
    ch_virus_index      // path: pre-built or freshly-built virus deacon index
    ch_depletion_index  // tuple(use_depletion, path): resolved host/contaminant depletion index or sentinel

    main:
    DEACON_FILTER_CONTIGS(
        ch_contigs
            .combine(ch_virus_index)
            .combine(ch_depletion_index)
    )

    emit:
    contigs = DEACON_FILTER_CONTIGS.out  // tuple(sample_id, platform, read_structure, screened_fasta, lookup)
}
