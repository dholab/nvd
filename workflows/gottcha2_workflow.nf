include {
    GOTTCHA2_PROFILE_NANOPORE ;
    GOTTCHA2_PROFILE_ILLUMINA ;
    GENERATE_FASTA ;
    REGISTER_GOTTCHA2_HITS
} from "../modules/gottcha2"
include {
    VALIDATE_LK_GOTTCHA2
} from '../subworkflows/validate_lk_gottcha2_lists'
include {
    BUNDLE_GOTTCHA2_FOR_LABKEY
} from "../subworkflows/bundle_gottcha2_for_labkey"
include { REMOVE_MULTIMAPS } from "../modules/bbmap"
include { CHECK_RUN_STATE; COMPLETE_RUN } from "../modules/utils"


workflow GOTTCHA2_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel: tuple val(sample_id), val(platform), val(read_structure), path(fastq)

    main:
    // Validate required params when GOTTCHA2 is selected
    if (params.tools && (params.tools.contains("gottcha") || params.tools.contains("all"))) {
        NvdUtils.validateLabkeyGottcha2(params)
    }

    // Resolve directories for stateful vs stateless mode.
    // Returns absolute path strings (not Nextflow path objects) to avoid staging.
    dirs = NvdDirs.resolve(params, log)

    // Combine samplesheet, state_dir, and taxonomy_dir into a tuple for CHECK_RUN_STATE.
    ch_run_state_input = Channel.fromPath(params.samplesheet)
        .combine(Channel.value(dirs.state_dir))
        .combine(Channel.value(dirs.taxonomy_dir))

    // Check run state upfront (prevents duplicate processing of same sample set)
    CHECK_RUN_STATE(ch_run_state_input)

    // Convert run_context from queue channel to value channel so it can be
    // consumed by multiple processes (hit registration, LabKey uploads, etc.).
    ch_run_context = CHECK_RUN_STATE.out.run_context.first()

    // Extract state_dir from run_context for state-tracking processes
    ch_state_dir = ch_run_context.map { _sample_set_id, state_dir -> state_dir }

    // Create value channel for hits_dir
    ch_hits_dir = Channel.value(dirs.hits_dir)

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

    // Register GOTTCHA2 hits with idempotent keys in the state database.
    // Join extracted FASTA (post-multimap removal) with run context and hits_dir.
    ch_register_gottcha2_input = REMOVE_MULTIMAPS.out
        .combine(ch_run_context)
        .combine(ch_hits_dir)
        .map { sample_id, fasta, full_tsv, sample_set_id, state_dir, hits_dir ->
            tuple(sample_id, fasta, full_tsv, sample_set_id, state_dir, hits_dir)
        }

    REGISTER_GOTTCHA2_HITS(ch_register_gottcha2_input)

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

    if (params.labkey) {
        // Gate run completion on LabKey bundle completion
        ch_run_complete_gate = BUNDLE_GOTTCHA2_FOR_LABKEY.out.collect()
            .map { _logs -> true }
    } else {
        // Gate run completion on REGISTER_GOTTCHA2_HITS completion (all samples)
        ch_run_complete_gate = REGISTER_GOTTCHA2_HITS.out.collect()
            .map { _logs -> true }
    }

    // Mark the run as completed in the state database
    COMPLETE_RUN(
        ch_run_complete_gate,
        ch_state_dir,
        "completed"
    )

    ch_completion = GENERATE_FASTA.out.map { _results -> "GOTTCHA2 complete!" }

    emit:
    completion = ch_completion
}
