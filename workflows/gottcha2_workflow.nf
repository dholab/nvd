include {
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


workflow GOTTCHA2_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel: tuple val(sample_id), val(platform), val(read_structure), path(fastq)

    main:
    // Validate required params when GOTTCHA2 is selected
    if (params.tools && (params.tools.contains("gottcha") || params.tools.contains("all"))) {
        NvdUtils.validateLabkeyGottcha2(params)
    }

    // Guard channel declaration with ternary fallback to Channel.empty().
    // When --tools doesn't include 'gottcha', this param is null.
    // Empty channel maintains the DAG but results in no-op downstream processes.
    ch_gottcha2_db = params.gottcha2_db
        ? Channel
            .fromPath("${params.gottcha2_db}{.tax.tsv,.stats,.mmi}")
            .collect()
            .map { files ->
                def ref_mmi   = files.find { it.name.endsWith(".mmi") }
                def stats     = files.find { it.name.endsWith(".stats") }
                def tax_tsv   = files.find { it.name.endsWith(".tax.tsv") }
                tuple(ref_mmi, stats, tax_tsv)
            }
        : Channel.empty()

    // Gate sample channel: when gottcha2 is disabled, ch_gottcha2_enabled is empty,
    // so the combine produces nothing and no downstream operations execute.
    ch_gottcha2_enabled = params.gottcha2_db
        ? Channel.value(true)
        : Channel.empty()

    ch_gated_samples = ch_gottcha2_enabled
        .combine(ch_sample_fastqs)
        .map { _flag, sample_id, platform, read_structure, fastq ->
            tuple(sample_id, platform, read_structure, fastq)
        }

    // LabKey validation: no-input subworkflow, gated by if (matches stat_blast pattern)
    if (params.labkey) {
        VALIDATE_LK_GOTTCHA2()
    }

    // LabKey data gating: empty channel when disabled, value channel when enabled.
    // Combining data channels with ch_labkey_gate produces nothing when LabKey is off,
    // so BUNDLE_GOTTCHA2_FOR_LABKEY never executes without an if/else block.
    ch_labkey_gate = params.labkey
        ? Channel.value(true)
        : Channel.empty()

    ch_nanopore_fastqs = ch_gated_samples
        .filter { _sample_id, platform, _read_structure, _fastq -> platform == "nanopore" || platform == "ont" }
        .map { sample_id, _platform, _read_structure, fastq -> tuple(sample_id, file(fastq)) }
        .filter { _id, fastq -> file(fastq).countFastq() >= params.min_gottcha_reads }

    ch_illumina_fastqs = ch_gated_samples
        .filter { _sample_id, platform, _read_structure, _fastq -> platform == "illumina" }
        .map { sample_id, _platform, _read_structure, fastq -> tuple(sample_id, file(fastq)) }
        .filter { _id, fastq -> file(fastq).countFastq() >= params.min_gottcha_reads }

    GOTTCHA2_PROFILE_NANOPORE(
        ch_nanopore_fastqs.combine(ch_gottcha2_db)
    )

    GOTTCHA2_PROFILE_ILLUMINA(
        ch_illumina_fastqs.combine(ch_gottcha2_db)
    )

    GENERATE_FASTA(
        GOTTCHA2_PROFILE_NANOPORE.out.aligned.mix(GOTTCHA2_PROFILE_ILLUMINA.out.aligned)
    )

    ch_to_remove_multimaps = GENERATE_FASTA.out
        .map { id, fasta, full_tsv, _log, _lineages -> tuple(id, fasta, full_tsv) }

    REMOVE_MULTIMAPS(ch_to_remove_multimaps)

    // LabKey upload: gated by ch_labkey_gate. When LabKey is disabled, the combine
    // with ch_labkey_gate produces empty channels, so the bundle never executes.
    ch_gottcha2_results_for_labkey = GOTTCHA2_PROFILE_NANOPORE.out.full_tsv
        .mix(GOTTCHA2_PROFILE_ILLUMINA.out.full_tsv)
        .combine(ch_labkey_gate)
        .map { sample_id, full_tsv, ref_mmi, stats, tax_tsv, _flag ->
            tuple(sample_id, full_tsv, ref_mmi, stats, tax_tsv)
        }

    ch_extracted_for_labkey = REMOVE_MULTIMAPS.out
        .combine(ch_labkey_gate)
        .map { sample_id, fasta, full_tsv, _flag ->
            tuple(sample_id, fasta, full_tsv)
        }

    BUNDLE_GOTTCHA2_FOR_LABKEY(
        ch_gottcha2_results_for_labkey,
        ch_extracted_for_labkey
    )

    ch_completion = GENERATE_FASTA.out.map { _results -> "GOTTCHA2 complete!" }

    emit:
    completion = ch_completion
}
