include { ASSEMBLE_WITH_SPADES } from "../modules/spades"
include { CONCAT_READS_AS_FASTA } from "../modules/seqkit"

workflow SHORT_READ_DENOVO_ASSEMBLY {
    take:
    ch_short_read_shards  // tuple(sample_id, platform, read_structure, evidence_class, fastq)

    main:

    ch_grouped_shards = ch_short_read_shards
        .groupTuple(by: [0, 1, 2])
        .map { sample_id, platform, read_structure, _evidence_classes, reads ->
            def count = reads.collect { read_file -> file(read_file).countFastq() }.sum()
            tuple(sample_id, platform, read_structure, reads, count)
        }
        .filter { _sample_id, _platform, _read_structure, _reads, count -> count >= 100 }

    ch_shards_by_structure = ch_grouped_shards.branch { _sample_id, _platform, read_structure, _reads, _count ->
        interleaved: read_structure == "interleaved"
        single: true
    }

    ch_interleaved_spades_inputs = ch_shards_by_structure.interleaved.map { sample_id, platform, read_structure, reads, _count ->
        tuple(sample_id, platform, read_structure, reads[0])
    }

    ch_single_shards = ch_shards_by_structure.single.map { sample_id, platform, read_structure, reads, _count ->
        tuple(sample_id, platform, read_structure, reads)
    }

    CONCAT_READS_AS_FASTA(ch_single_shards)

    ch_spades_inputs = ch_interleaved_spades_inputs.mix(CONCAT_READS_AS_FASTA.out)

    ASSEMBLE_WITH_SPADES(ch_spades_inputs)

    emit:
    contigs = ASSEMBLE_WITH_SPADES.out  // tuple(sample_id, platform, read_structure, producer, fasta)
}
