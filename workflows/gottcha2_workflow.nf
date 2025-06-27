include {
    READ_COMPRESSION_PASSTHROUGH ;
    GOTTCHA2_PROFILE_NANOPORE ;
    GOTTCHA2_PROFILE_ILLUMINA ;
    GENERATE_FASTA
} from "../modules/gottcha2"
include { 
    VALIDATE_LK_GOTTCHA2
 } from '../subworkflows/validate_lk_gottcha2_lists.nf'
include {
    BUNDLE_GOTTCHA2_FOR_LABKEY
} from "../subworkflows/bundle_gottcha2_for_labkey"


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

    VALIDATE_LK_GOTTCHA2()

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

    BUNDLE_GOTTCHA2_FOR_LABKEY(
        GOTTCHA2_PROFILE_NANOPORE.out.full_tsv.mix(GOTTCHA2_PROFILE_ILLUMINA.out.full_tsv),
        GENERATE_FASTA.out
    )

    ch_completion = GENERATE_FASTA.out.map{ _results -> "GOTTCHA2 complete!" }

    emit:
    completion = ch_completion
}


