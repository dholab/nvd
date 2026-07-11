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

        ch_available_assembler_contigs = ASSEMBLE_WITH_MYLOASM.out.contigs
            .mix(ASSEMBLE_WITH_METAMDBG.out.contigs)
            .mix(ASSEMBLE_WITH_METAFLYE.out.contigs)
            .groupTuple(by: [0, 1, 2], size: 3, remainder: true)
            .map { sample_id, platform, read_structure, assemblers, contigs ->
                tuple(sample_id, platform, read_structure, tuple(assemblers, contigs))
            }

        ch_post_qc_read_context = ch_long_read_profiles.map { meta, post_qc_reads, _profile_json, _length_histogram ->
            tuple(meta.id, meta.platform, meta.read_structure, post_qc_reads)
        }

        ch_long_read_query_candidates = ch_post_qc_read_context.join(
            ch_available_assembler_contigs,
            by: [0, 1, 2],
            remainder: true,
            failOnDuplicate: true,
        )

        ch_long_read_query_candidates_by_availability = ch_long_read_query_candidates.branch {
            _sample_id, _platform, _read_structure, _post_qc_reads, assembler_contigs ->
            with_contigs: assembler_contigs != null
            reads_only: true
        }

        ch_no_contigs = ch_long_read_query_candidates_by_availability.reads_only.map {
            sample_id, platform, _read_structure, _post_qc_reads, _assembler_contigs ->
            log.debug "nvd.contig_route sample_id=${sample_id} platform=${platform} outcome=no_contigs stage=long_read_assembly"
            tuple(sample_id, platform)
        }

        NORMALIZE_CONTIGS(
            ch_long_read_query_candidates_by_availability.with_contigs.map {
                sample_id, platform, read_structure, _post_qc_reads, assembler_contigs ->
                tuple(sample_id, platform, read_structure, assembler_contigs[0], assembler_contigs[1])
            }
        )
        FIND_CONTAINMENT_DUPLICATES(NORMALIZE_CONTIGS.out.prepared)
        EXTRACT_UNIQUE_CONTIGS(FIND_CONTAINMENT_DUPLICATES.out.candidates)
        ch_long_read_contigs = EXTRACT_UNIQUE_CONTIGS.out.contigs
    } else {
        ch_long_read_assembly_inputs = ch_long_read_profiles.branch { meta, _reads, _profile_json, _length_histogram ->
            eligible: meta.sequence_count >= 100
            ineligible: true
        }

        ch_no_contigs = ch_long_read_assembly_inputs.ineligible.map { meta, _reads, _profile_json, _length_histogram ->
            log.debug "nvd.contig_route sample_id=${meta.id} platform=${meta.platform} outcome=no_contigs stage=long_read_ineligible"
            tuple(meta.id, meta.platform)
        }

        ASSEMBLE_WITH_SPADES(
            ch_long_read_assembly_inputs.eligible
                .map { meta, reads, _profile_json, _length_histogram ->
                    tuple(meta.id, meta.platform, meta.read_structure, reads)
                }
        )
        ch_long_read_contigs = ASSEMBLE_WITH_SPADES.out
    }

    emit:
    contigs = ch_long_read_contigs  // tuple(sample_id, platform, read_structure, producer, fasta)
    no_contigs = ch_no_contigs      // tuple(sample_id, platform)
}
