include { RUN_SPADES } from "../modules/spades"
include { RUN_MYLOASM ; RUN_METAMDBG ; RUN_METAFLYE ; CONCAT_LONG_READ_ASSEMBLIES } from "../modules/long_read_assembly"

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
            .map { sample_id, platform, read_structure, _assembler, fasta ->
                tuple(sample_id, platform, read_structure, fasta)
            }
            .groupTuple(by: [0, 1, 2])

        CONCAT_LONG_READ_ASSEMBLIES(ch_assembler_contigs)
        ch_long_read_contigs = CONCAT_LONG_READ_ASSEMBLIES.out.contigs
    } else {
        RUN_SPADES(ch_long_read_fastqs)
        ch_long_read_contigs = RUN_SPADES.out
    }

    emit:
    contigs = ch_long_read_contigs  // tuple(sample_id, platform, read_structure, fasta)
}
