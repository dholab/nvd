include { ASSEMBLE_WITH_SPADES } from "../modules/spades"
include { CONCAT_READS_AS_FASTA } from "../modules/seqkit"

workflow SHORT_READ_DENOVO_ASSEMBLY {
    take:
    ch_short_read_shards  // tuple(meta, reads, profile_json, length_histogram); meta includes id, platform, read_structure, query_class, sequence_count

    main:

    ch_grouped_shards = ch_short_read_shards
        .map { meta, reads, _profile_json, _length_histogram ->
            tuple(meta.id, meta.platform, meta.read_structure, [meta: meta, reads: reads])
        }
        .groupTuple(by: [0, 1, 2])
        .map { sample_id, platform, read_structure, rows ->
            def count = rows.collect { row -> row.meta.sequence_count }.sum()
            def meta = [id: sample_id, platform: platform, read_structure: read_structure, sequence_count: count]
            def reads = rows.collect { row -> row.reads }
            tuple(meta, reads)
        }
        .filter { meta, _reads -> meta.sequence_count >= 100 }

    ch_shards_by_structure = ch_grouped_shards.branch { meta, _reads ->
        interleaved: meta.read_structure == "interleaved"
        single: true
    }

    ch_interleaved_spades_inputs = ch_shards_by_structure.interleaved.map { meta, reads ->
        tuple(meta.id, meta.platform, meta.read_structure, reads[0])
    }

    ch_single_shards = ch_shards_by_structure.single.map { meta, reads ->
        tuple(meta.id, meta.platform, meta.read_structure, reads)
    }

    CONCAT_READS_AS_FASTA(ch_single_shards)

    ch_spades_inputs = ch_interleaved_spades_inputs.mix(CONCAT_READS_AS_FASTA.out)

    ASSEMBLE_WITH_SPADES(ch_spades_inputs)

    emit:
    contigs = ASSEMBLE_WITH_SPADES.out  // tuple(sample_id, platform, read_structure, producer, fasta)
}
