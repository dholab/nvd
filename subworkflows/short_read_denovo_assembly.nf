include { ASSEMBLE_WITH_SPADES } from "../modules/spades"
include { CONCAT_READS_AS_FASTA } from "../modules/seqkit"

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

    emit:
    contigs = ASSEMBLE_WITH_SPADES.out  // tuple(sample_id, platform, read_structure, producer, fasta)
    no_contigs = ch_no_contigs          // tuple(sample_id, platform)
}
