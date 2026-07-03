include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"

workflow PREPROCESS_CONTIGS {
    take:
    ch_assembled_contigs  // tuple(sample_id, platform, read_structure, fasta)

    main:

    MASK_LOW_COMPLEXITY(
        ch_assembled_contigs
    )

    FILTER_SHORT_CONTIGS(
        MASK_LOW_COMPLEXITY.out
    )

    ch_filtered_contigs = FILTER_SHORT_CONTIGS.out
        .filter { _id, _platform, _read_structure, fasta ->
            def recordCount = fasta.countFasta()
            recordCount > 0
        }

    emit:
    contigs = ch_filtered_contigs  // tuple(sample_id, platform, read_structure, fasta)
}
