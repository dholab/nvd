include { ASSEMBLE_WITH_SPADES } from "../modules/spades"
include { ASSESS_LONG_READ_ASSEMBLY_ELIGIBILITY ; REPORT_LONG_READ_ASSEMBLY_ELIGIBILITY ; ASSEMBLE_WITH_MYLOASM ; ASSEMBLE_WITH_METAMDBG ; ASSEMBLE_WITH_METAFLYE } from "../modules/long_read_assembly"
include { NORMALIZE_CONTIGS ; FIND_CONTAINMENT_DUPLICATES ; EXTRACT_UNIQUE_CONTIGS } from "../modules/contig_union"

workflow LONG_READ_DENOVO_ENSEMBLY {
    take:
    ch_long_read_profiles  // tuple(meta, reads, profile_json, length_histogram)

    main:

    if (params.experimental == true) {
        ASSESS_LONG_READ_ASSEMBLY_ELIGIBILITY(
            ch_long_read_profiles.map { meta, reads, profile_json, _length_histogram ->
                tuple(meta, reads, profile_json)
            }
        )

        ASSEMBLE_WITH_MYLOASM(
            ASSESS_LONG_READ_ASSEMBLY_ELIGIBILITY.out.myloasm.map { meta, reads, _marker ->
                tuple(meta.id, meta.platform, meta.read_structure, reads)
            }
        )
        ASSEMBLE_WITH_METAMDBG(
            ASSESS_LONG_READ_ASSEMBLY_ELIGIBILITY.out.metamdbg.map { meta, reads, _marker ->
                tuple(meta.id, meta.platform, meta.read_structure, reads)
            }
        )
        ASSEMBLE_WITH_METAFLYE(
            ASSESS_LONG_READ_ASSEMBLY_ELIGIBILITY.out.metaflye.map { meta, reads, _marker ->
                tuple(meta.id, meta.platform, meta.read_structure, reads)
            }
        )

        REPORT_LONG_READ_ASSEMBLY_ELIGIBILITY(
            ASSESS_LONG_READ_ASSEMBLY_ELIGIBILITY.out.report
                .map { _sample_id, report -> report }
                .collect()
                .filter { reports -> !reports.isEmpty() }
        )

        ch_assembler_contigs = ASSEMBLE_WITH_MYLOASM.out.contigs
            .mix(ASSEMBLE_WITH_METAMDBG.out.contigs)
            .mix(ASSEMBLE_WITH_METAFLYE.out.contigs)
            .groupTuple(by: [0, 1, 2], size: 3, remainder: true)

        NORMALIZE_CONTIGS(ch_assembler_contigs)
        FIND_CONTAINMENT_DUPLICATES(NORMALIZE_CONTIGS.out.prepared)
        EXTRACT_UNIQUE_CONTIGS(FIND_CONTAINMENT_DUPLICATES.out.candidates)
        ch_long_read_contigs = EXTRACT_UNIQUE_CONTIGS.out.contigs
    } else {
        ASSEMBLE_WITH_SPADES(
            ch_long_read_profiles
                .filter { meta, _reads, _profile_json, _length_histogram -> meta.sequence_count >= 100 }
                .map { meta, reads, _profile_json, _length_histogram ->
                    tuple(meta.id, meta.platform, meta.read_structure, reads)
                }
        )
        ch_long_read_contigs = ASSEMBLE_WITH_SPADES.out
    }

    emit:
    contigs = ch_long_read_contigs  // tuple(sample_id, platform, read_structure, producer, fasta)
}
