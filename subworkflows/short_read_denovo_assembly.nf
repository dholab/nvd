include { ASSEMBLE_WITH_SPADES } from "../modules/spades"
include { CONCAT_READS_AS_FASTA } from "../modules/seqkit"
include { PROFILE_ASSEMBLY_FASTA_FOR_REPORT as PROFILE_SHORT_READ_ASSEMBLY_REPORT } from "../modules/fastx"

workflow SHORT_READ_DENOVO_ASSEMBLY {
    take:
    ch_profiled_batches_by_sample  // tuple(sample_meta, batches); each batch contains meta and reads

    main:

    ch_sample_reads = ch_profiled_batches_by_sample
        .map { meta, batches ->
            tuple(meta, batches.collect { batch -> batch.reads })
        }

    ch_assembly_inputs = ch_sample_reads.branch { meta, _reads ->
        eligible: meta.sequence_count >= 100
        ineligible: true
    }

    ch_no_contigs = ch_assembly_inputs.ineligible.map { meta, _reads ->
        log.debug "nvd.contig_route sample_id=${meta.id} platform=${meta.platform} outcome=no_contigs stage=short_read_ineligible"
        tuple(meta.id, meta.platform)
    }

    ch_reads_by_structure = ch_assembly_inputs.eligible.branch { meta, _reads ->
        interleaved: meta.read_structure == "interleaved"
        single: true
    }

    ch_interleaved_spades_inputs = ch_reads_by_structure.interleaved.map { meta, reads ->
        tuple(meta.id, meta.platform, meta.read_structure, reads[0])
    }

    ch_single_batches = ch_reads_by_structure.single.map { meta, reads ->
        tuple(meta.id, meta.platform, meta.read_structure, reads)
    }

    CONCAT_READS_AS_FASTA(ch_single_batches)

    ch_spades_inputs = ch_interleaved_spades_inputs.mix(CONCAT_READS_AS_FASTA.out)

    ASSEMBLE_WITH_SPADES(ch_spades_inputs)

    PROFILE_SHORT_READ_ASSEMBLY_REPORT(ASSEMBLE_WITH_SPADES.out)

    ch_assembly_eligibility_decisions = ch_assembly_inputs.ineligible.map { meta, _reads ->
        tuple(meta + [producer: 'spades'], 'skip', meta.sequence_count, 100, 'short_read_minimum_sequence_count')
    }
        .mix(ch_assembly_inputs.eligible.map { meta, _reads ->
            tuple(meta + [producer: 'spades'], 'run', meta.sequence_count, 100, 'short_read_minimum_sequence_count')
        })

    emit:
    contigs = ASSEMBLE_WITH_SPADES.out  // tuple(sample_id, platform, read_structure, producer, fasta)
    no_contigs = ch_no_contigs          // tuple(sample_id, platform)
    assembly_profiles = PROFILE_SHORT_READ_ASSEMBLY_REPORT.out.profiled
    eligibility_decisions = ch_assembly_eligibility_decisions
}
