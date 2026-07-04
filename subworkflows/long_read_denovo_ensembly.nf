include { RUN_SPADES } from "../modules/spades"
include { RUN_MYLOASM ; RUN_METAMDBG ; RUN_METAFLYE } from "../modules/long_read_assembly"
include { NORMALIZE_CONTIGS ; FIND_CONTAINMENT_DUPLICATES ; EXTRACT_UNIQUE_CONTIGS } from "../modules/contig_union"

workflow LONG_READ_DENOVO_ENSEMBLY {
    take:
    ch_long_read_fastqs  // tuple(sample_id, platform, read_structure, fastq)

    main:

    if (params.experimental == true) {
        RUN_MYLOASM(ch_long_read_fastqs)
        RUN_METAMDBG(ch_long_read_fastqs)
        RUN_METAFLYE(ch_long_read_fastqs)

        ch_assembler_contigs = RUN_MYLOASM.out.contigs
            .mix(RUN_METAMDBG.out.contigs)
            .mix(RUN_METAFLYE.out.contigs)
            .groupTuple(by: [0, 1, 2], size: 3)

        NORMALIZE_CONTIGS(ch_assembler_contigs)
        FIND_CONTAINMENT_DUPLICATES(NORMALIZE_CONTIGS.out.prepared)
        EXTRACT_UNIQUE_CONTIGS(FIND_CONTAINMENT_DUPLICATES.out.candidates)
        ch_long_read_contigs = EXTRACT_UNIQUE_CONTIGS.out.contigs
    } else {
        RUN_SPADES(ch_long_read_fastqs)
        ch_long_read_contigs = RUN_SPADES.out
    }

    emit:
    contigs = ch_long_read_contigs  // tuple(sample_id, platform, read_structure, fasta)
}
