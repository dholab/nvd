include { ASSEMBLE_WITH_SPADES } from "../modules/spades"
include { ASSEMBLE_WITH_MYLOASM ; ASSEMBLE_WITH_METAMDBG ; ASSEMBLE_WITH_METAFLYE } from "../modules/long_read_assembly"
include { NORMALIZE_CONTIGS ; FIND_CONTAINMENT_DUPLICATES ; EXTRACT_UNIQUE_CONTIGS } from "../modules/contig_union"

workflow LONG_READ_DENOVO_ENSEMBLY {
    take:
    ch_long_read_fastqs  // tuple(sample_id, platform, read_structure, fastq)

    main:

    if (params.experimental == true) {
        ASSEMBLE_WITH_MYLOASM(ch_long_read_fastqs)
        ASSEMBLE_WITH_METAMDBG(ch_long_read_fastqs)
        ASSEMBLE_WITH_METAFLYE(ch_long_read_fastqs)

        ch_assembler_contigs = ASSEMBLE_WITH_MYLOASM.out.contigs
            .mix(ASSEMBLE_WITH_METAMDBG.out.contigs)
            .mix(ASSEMBLE_WITH_METAFLYE.out.contigs)
            .groupTuple(by: [0, 1, 2], size: 3)

        NORMALIZE_CONTIGS(ch_assembler_contigs)
        FIND_CONTAINMENT_DUPLICATES(NORMALIZE_CONTIGS.out.prepared)
        EXTRACT_UNIQUE_CONTIGS(FIND_CONTAINMENT_DUPLICATES.out.candidates)
        ch_long_read_contigs = EXTRACT_UNIQUE_CONTIGS.out.contigs
    } else {
        ASSEMBLE_WITH_SPADES(ch_long_read_fastqs)
        ch_long_read_contigs = ASSEMBLE_WITH_SPADES.out
    }

    emit:
    contigs = ch_long_read_contigs  // tuple(sample_id, platform, read_structure, fasta)
}
