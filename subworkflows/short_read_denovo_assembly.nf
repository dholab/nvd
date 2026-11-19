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

    ch_no_contigs_from_ineligibility = ch_assembly_inputs.ineligible.map { meta, _reads ->
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

    // ASSEMBLE_WITH_SPADES retries twice and then ignores the failure, and an
    // ignored task emits nothing. Joining the samples we *attempted* against the
    // contigs we actually got turns that silence into a null we can route on, so
    // a crashed assembler still reaches the read-level BLAST fallback. Keying on
    // the eligible samples rather than on ch_spades_inputs also covers an earlier
    // failure in CONCAT_READS_AS_FASTA.
    ch_attempted_assemblies = ch_assembly_inputs.eligible.map { meta, _reads ->
        tuple(meta.id, meta.platform)
    }

    ch_assemblies_by_outcome = ch_attempted_assemblies
        .join(
            ASSEMBLE_WITH_SPADES.out,
            by: [0, 1],
            remainder: true,
            failOnDuplicate: true,
        )
        // A matched sample arrives as the full assembly tuple; an unmatched one
        // arrives as the key plus a single null placeholder, so the two differ in
        // arity. Take the whole tuple as one list and test its last element rather
        // than destructuring, which would fail on the shorter tuple.
        .branch { assembly_outcome ->
            produced_output: assembly_outcome.last() != null
            produced_nothing: true
        }

    ch_no_contigs_from_assembly_failure = ch_assemblies_by_outcome.produced_nothing.map { assembly_outcome ->
        def sample_id = assembly_outcome[0]
        def platform = assembly_outcome[1]
        log.warn "Assembly produced no output for sample ${sample_id}; falling back to read-level BLAST queries."
        log.debug "nvd.contig_route sample_id=${sample_id} platform=${platform} outcome=no_contigs stage=assembly_failed"
        tuple(sample_id, platform)
    }

    ch_no_contigs = ch_no_contigs_from_ineligibility.mix(ch_no_contigs_from_assembly_failure)

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
