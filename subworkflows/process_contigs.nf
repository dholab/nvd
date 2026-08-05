include { MASK_LOW_COMPLEXITY ; FILTER_SHORT_CONTIGS } from "../modules/bbmap"
include { COLLECT_CONTIGS } from "../modules/contigs"

workflow PROCESS_CONTIGS {
    take:
    ch_assembled_contigs  // tuple(sample_id, platform, read_structure, producer, fasta)

    main:

    COLLECT_CONTIGS(ch_assembled_contigs)

    ch_collected_by_content = COLLECT_CONTIGS.out.contigs.branch { _sample_id, _platform, _read_structure, fasta, _lookup ->
        nonempty: file(fasta).size() > 0
        empty: true
    }

    ch_empty_after_collection = ch_collected_by_content.empty.map { sample_id, platform, _read_structure, _fasta, _lookup ->
        log.debug "nvd.contig_route sample_id=${sample_id} platform=${platform} outcome=no_contigs stage=collect_contigs"
        tuple(sample_id, platform)
    }

    MASK_LOW_COMPLEXITY(
        ch_collected_by_content.nonempty
    )

    FILTER_SHORT_CONTIGS(
        MASK_LOW_COMPLEXITY.out
    )

    ch_length_filtered_by_content = FILTER_SHORT_CONTIGS.out.branch { _sample_id, _platform, _read_structure, fasta, _lookup ->
        nonempty: file(fasta).size() > 0
        empty: true
    }

    ch_empty_after_length_filter = ch_length_filtered_by_content.empty.map { sample_id, platform, _read_structure, _fasta, _lookup ->
        log.debug "nvd.contig_route sample_id=${sample_id} platform=${platform} outcome=no_contigs stage=length_filter"
        tuple(sample_id, platform)
    }

    emit:
    contigs = ch_length_filtered_by_content.nonempty  // tuple(sample_id, platform, read_structure, fasta, lookup)
    no_contigs = ch_empty_after_collection.mix(ch_empty_after_length_filter)  // tuple(sample_id, platform)
}
