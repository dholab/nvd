include {
    READ_COMPRESSION_PASSTHROUGH ;
    GOTTCHA2_PROFILE_NANOPORE ;
    GOTTCHA2_PROFILE_ILLUMINA ;
    GENERATE_FASTA
} from "../modules/gottcha2"

workflow GOTTCHA2_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel of sample IDs, platforms, and (interleaved) FASTQ files: tuple val(sample_id), val(platform), path(fastq)

    main:
    ch_gottcha2_db = Channel.fromPath("${params.gottcha2_db}*").collect(sort: true)

    ch_nanopore_fastqs = ch_sample_fastqs
        .filter { _sample_id, platform, _fastq -> platform == "nanopore" || platform == "ont" }
        .map { sample_id, _platform, fastq -> tuple(sample_id, file(fastq)) }

    ch_illumina_fastqs = ch_sample_fastqs
        .filter { _sample_id, platform, _fastq -> platform == "illumina" }
        .map { sample_id, _platform, fastq -> tuple(sample_id, file(fastq)) }

    READ_COMPRESSION_PASSTHROUGH(
        ch_nanopore_fastqs.combine(ch_gottcha2_db.collect(sort: true))
    )

    GOTTCHA2_PROFILE_NANOPORE(
        READ_COMPRESSION_PASSTHROUGH.out
    )

    GOTTCHA2_PROFILE_ILLUMINA(
        ch_illumina_fastqs.combine(ch_gottcha2_db.collect(sort: true))
    )

    GENERATE_FASTA(
        GOTTCHA2_PROFILE_NANOPORE.out.aligned.mix(GOTTCHA2_PROFILE_ILLUMINA.out.aligned)
    )

    ch_completion = GENERATE_FASTA.out.map{ _results -> "GOTTCHA2 complete!" }

    emit:
    completion = ch_completion
}


