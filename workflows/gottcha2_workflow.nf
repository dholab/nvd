include {
    READ_COMPRESSION_PASSTHROUGH ;
    GOTTCHA2_PROFILE_NANOPORE ;
    GOTTCHA2_PROFILE_ILLUMINA ;
    GENERATE_FASTA
} from "../modules/gottcha2"
include { 
    VALIDATE_LK_GOTTCHA2
 } from '../subworkflows/validate_lk_gottcha2_lists'
include {
    BUNDLE_GOTTCHA2_FOR_LABKEY
} from "../subworkflows/bundle_gottcha2_for_labkey"
include { REMOVE_MULTIMAPS } from "../modules/bbmap"
include {
    VERIFY_WITH_BLAST ;
    STACK_VERIFIED_TABLES ;
    FILTER_DOWN_TO_SPECIES_STRAIN ;
    EXTRACT_WHITELISTED_FASTAS
} from "../modules/blast"


workflow GOTTCHA2_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel of sample IDs, platforms, and (interleaved) FASTQ files: tuple val(sample_id), val(platform), path(fastq)

    main:
    ch_gottcha2_db = Channel.fromPath("${params.gottcha2_db}*").collect(sort: true)

    if (params.labkey) {
        VALIDATE_LK_GOTTCHA2()
    }

    READ_COMPRESSION_PASSTHROUGH(ch_sample_fastqs)

    ch_nanopore_fastqs = READ_COMPRESSION_PASSTHROUGH.out
        .filter { _sample_id, platform, _fastq -> platform == "nanopore" || platform == "ont" }
        .map { sample_id, _platform, fastq -> tuple(sample_id, file(fastq)) }
        .filter { _id, fastq -> file(fastq).countFastq() >= params.min_gottcha_reads }

    ch_illumina_fastqs = READ_COMPRESSION_PASSTHROUGH.out
        .filter { _sample_id, platform, _fastq -> platform == "illumina" }
        .map { sample_id, _platform, fastq -> tuple(sample_id, file(fastq)) }
        .filter { _id, fastq -> file(fastq).countFastq() >= params.min_gottcha_reads }

    GOTTCHA2_PROFILE_NANOPORE(
        ch_nanopore_fastqs.combine(ch_gottcha2_db.collect(sort: true))
    )

    GOTTCHA2_PROFILE_ILLUMINA(
        ch_illumina_fastqs.combine(ch_gottcha2_db.collect(sort: true))
    )

    GENERATE_FASTA(
        GOTTCHA2_PROFILE_NANOPORE.out.aligned.mix(GOTTCHA2_PROFILE_ILLUMINA.out.aligned)
    )

    // Map the channel to match multimaps input
    ch_to_remove_multimaps = GENERATE_FASTA.out
        .map {id, fasta, full_tsv, _log, _lineages -> tuple(id, fasta, full_tsv)}

    REMOVE_MULTIMAPS(ch_to_remove_multimaps)

    if (params.labkey) {

        BUNDLE_GOTTCHA2_FOR_LABKEY(
            GOTTCHA2_PROFILE_NANOPORE.out.full_tsv.mix(GOTTCHA2_PROFILE_ILLUMINA.out.full_tsv),
            REMOVE_MULTIMAPS.out
        )
    }

    ch_completion = GENERATE_FASTA.out.map{ _results -> "GOTTCHA2 complete!" }

    // ch_completion = VERIFY_WITH_BLAST.out.map{ _results -> "GOTTCHA2 complete!" }

    emit:
    completion = ch_completion
    verification_results = STACK_VERIFIED_TABLES.out
}


