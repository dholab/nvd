include {
    GOTTCHA2_PROFILE_NANOPORE ;
    GOTTCHA2_PROFILE_ILLUMINA ;
    GENERATE_FASTA
} from "../modules/gottcha2"

workflow GOTTCHA2_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel of sample IDs, platforms, and (interleaved) FASTQ files: tuple val(sample_id), val(platform), path(fastq)

    main:
    ch_gottcha2_db = Channel.fromPath("${params.gottcha2_db}*")

    ch_nanopore_fastqs = ch_sample_fastqs
        .filter { _sample_id, platform, _fastq -> platform == "nanopore" }
        .map { sample_id, _platform, fastq -> tuple(sample_id, file(fastq)) }

    ch_illumina_fastqs = ch_sample_fastqs
        .filter { _sample_id, platform, _fastq -> platform == "illumina" }
        .map { sample_id, _platform, fastq -> tuple(sample_id, file(fastq)) }

    GOTTCHA2_PROFILE_NANOPORE(
        ch_nanopore_fastqs.combine(ch_gottcha2_db.collect())
    )

    GOTTCHA2_PROFILE_ILLUMINA(
        ch_illumina_fastqs.combine(ch_gottcha2_db.collect())
    )

    GENERATE_FASTA(
        GOTTCHA2_PROFILE_NANOPORE.out.aligned.mix(GOTTCHA2_PROFILE_ILLUMINA.out.aligned)
    )
}
