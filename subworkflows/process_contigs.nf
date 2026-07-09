include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"
include { COLLECT_CONTIGS } from "../modules/contigs"

workflow PROCESS_CONTIGS {
    take:
    ch_assembled_contigs  // tuple(sample_id, platform, read_structure, producer, fasta)

    main:

    COLLECT_CONTIGS(ch_assembled_contigs)

    MASK_LOW_COMPLEXITY(
        COLLECT_CONTIGS.out.contigs
    )

    FILTER_SHORT_CONTIGS(
        MASK_LOW_COMPLEXITY.out
    )

    emit:
    contigs = FILTER_SHORT_CONTIGS.out  // tuple(sample_id, platform, read_structure, fasta, lookup)
}
