/*
 * STAT+BLAST Workflow
 * 
 * Human virus detection using NCBI STAT for initial classification
 * followed by two-phase BLAST verification (megablast + blastn).
 * 
 * Historical note: Previously called "NVD" (Novel Virus Detection) workflow.
 * Tool aliases: nvd, stat, blast, stat_blast, stast
 */

nextflow.enable.dsl=2

include { PREPROCESS_CONTIGS } from "../subworkflows/preprocess_contigs"
include { EXTRACT_HUMAN_VIRUSES } from "../subworkflows/extract_human_virus_contigs"
include { CLASSIFY_WITH_MEGABLAST } from "../subworkflows/classify_with_megablast"
include { CLASSIFY_WITH_BLASTN } from "../subworkflows/classify_with_blastn"
include { BUNDLE_BLAST_FOR_LABKEY } from "../subworkflows/bundle_blast_for_labkey"
include { COUNT_READS } from "../modules/count_reads"
include { CHECK_RUN_STATE } from "../modules/utils"
include { VALIDATE_LK_BLAST } from "../subworkflows/validate_lk_blast_lists.nf"
include { VALIDATE_LK_EXP_FRESH } from "../modules/validate_blast_labkey.nf"
include { REGISTER_LK_EXPERIMENT } from "../modules/validate_blast_labkey.nf"


workflow STAT_BLAST_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel of sample IDs, platforms, and (interleaved) FASTQ files: tuple val(sample_id), val(platform), path(fastq)

    main:
    // Check if STAT+BLAST workflow should run
    // Supports: nvd, stat, blast, stat_blast, stast (all equivalent for backward compatibility)
    if (params.tools && (params.tools.contains("nvd") || 
                         params.tools.contains("stat") || 
                         params.tools.contains("blast") || 
                         params.tools.contains("stat_blast") || 
                         params.tools.contains("stast") || 
                         params.tools.contains("all"))) {
        assert (
            params.blast_db              && file(params.blast_db).isDirectory()         &&
            params.stat_index            && file(params.stat_index).exists()            &&
            params.stat_dbss             && file(params.stat_dbss).exists()             &&
            params.stat_annotation       && file(params.stat_annotation).exists()       &&
            params.human_virus_taxlist   && file(params.human_virus_taxlist).exists()
        ) :
        """
        One or more required parameters are missing or point to non-existent files:

          blast_db            -> ${params.blast_db}
          stat_index          -> ${params.stat_index}
          stat_dbss           -> ${params.stat_dbss}
          stat_annotation     -> ${params.stat_annotation}
          human_virus_taxlist -> ${params.human_virus_taxlist}

        Please supply all of the above in your `-c nextflow.config` or via `-params-file`, and ensure each path exists.
        """

        // Validate LabKey params required for STAT+BLAST if LabKey is enabled
        NvdUtils.validateLabkeyBlast(params)
    }

    // Guard channel declarations with ternary fallback to Channel.empty().
    // When --tools doesn't include a BLAST alias, these params are null.
    // Empty channels maintain the DAG but result in no-op downstream processes.
    ch_blast_db_files = params.blast_db 
        ? Channel.fromPath(params.blast_db) 
        : Channel.empty()
    ch_stat_index = params.stat_index 
        ? Channel.fromPath(params.stat_index) 
        : Channel.empty()
    ch_stat_dbss = params.stat_dbss 
        ? Channel.fromPath(params.stat_dbss) 
        : Channel.empty()
    ch_stat_annotation = params.stat_annotation 
        ? Channel.fromPath(params.stat_annotation) 
        : Channel.empty()
    ch_human_virus_taxlist = params.human_virus_taxlist 
        ? Channel.fromPath(params.human_virus_taxlist) 
        : Channel.empty()

    // State directory for run tracking and upload deduplication.
    // Using Channel.fromPath with type:'dir' ensures Nextflow stages/mounts
    // this directory automatically in containerized environments.
    // The directory is created if it doesn't exist (matching db.py behavior).
    def state_dir_path = file(params.state_dir)
    if (!state_dir_path.exists()) {
        state_dir_path.mkdirs()
    }

    // Combine samplesheet and state_dir into a tuple for CHECK_RUN_STATE.
    // This ensures they travel together and we get visible cross-product failures
    // if either channel unexpectedly emits multiple values.
    ch_run_state_input = Channel.fromPath(params.samplesheet)
        .combine(Channel.fromPath(params.state_dir, type: 'dir', checkIfExists: true))

    // Check run state upfront (prevents duplicate processing of same sample set)
    // Reads sample IDs directly from samplesheet for immediate validation
    // Emits run_context tuple: [sample_set_id, state_dir] for downstream upload processes
    CHECK_RUN_STATE(ch_run_state_input)

    if (params.labkey) {
        VALIDATE_LK_BLAST()
    }

    // Count reads for each sample
    COUNT_READS(ch_sample_fastqs)

    PREPROCESS_CONTIGS(
        ch_sample_fastqs,
        ch_stat_dbss,
        ch_stat_annotation,
        ch_human_virus_taxlist
    )

    EXTRACT_HUMAN_VIRUSES(
        PREPROCESS_CONTIGS.out.contigs,
        PREPROCESS_CONTIGS.out.viral_reads,
        ch_stat_index,
        ch_stat_dbss,
        ch_stat_annotation
    )

    CLASSIFY_WITH_MEGABLAST(
        EXTRACT_HUMAN_VIRUSES.out.contigs,
        ch_blast_db_files
    )

    CLASSIFY_WITH_BLASTN(
        CLASSIFY_WITH_MEGABLAST.out.filtered_megablast,
        CLASSIFY_WITH_MEGABLAST.out.megablast_contigs,
        ch_blast_db_files
    )

    // Fork merged_results for multiple consumers.
    // Without this, .first() for the trigger would consume the first sample,
    // leaving only remaining samples for BUNDLE_BLAST_FOR_LABKEY.
    merged_results_for_labkey = CLASSIFY_WITH_BLASTN.out.merged_results
        .tap { merged_results_for_trigger }
        .tap { merged_results_for_completion }

    if (params.labkey) {
        // Check LabKey guard list - ensures experiment_id hasn't been uploaded before
        // Uses .first() on a forked channel to get a value channel trigger
        VALIDATE_LK_EXP_FRESH(
            merged_results_for_trigger.first()
        )

        // Bundle and upload - gates ensure both checks passed:
        // 1. CHECK_RUN_STATE.out.ready: local sample_set_id check passed
        // 2. VALIDATE_LK_EXP_FRESH.out.validated: LabKey guard list check passed
        // We combine both gates to ensure both validations complete before uploading
        ch_validation_gate = CHECK_RUN_STATE.out.ready
            .combine(VALIDATE_LK_EXP_FRESH.out.validated)
            .map { _ready, _validated -> true }

        // Convert run_context from queue channel to value channel so it can be
        // consumed by multiple processes (LABKEY_UPLOAD_BLAST and LABKEY_UPLOAD_FASTA).
        // This is safe because CHECK_RUN_STATE emits exactly once per pipeline run.
        ch_run_context = CHECK_RUN_STATE.out.run_context.first()

        BUNDLE_BLAST_FOR_LABKEY(
            merged_results_for_labkey,
            EXTRACT_HUMAN_VIRUSES.out.contigs,
            COUNT_READS.out.counts,
            params.experiment_id,
            workflow.runName,
            EXTRACT_HUMAN_VIRUSES.out.contig_read_counts,
            ch_validation_gate,
            ch_run_context  // value channel: [sample_set_id, state_dir]
        )

        // Register experiment ID in guard list after successful upload
        REGISTER_LK_EXPERIMENT(
            BUNDLE_BLAST_FOR_LABKEY.out.upload_log.collect()
        )

        labkey_log_ch = BUNDLE_BLAST_FOR_LABKEY.out.upload_log
    } else {
        // If no labkey upload then just pass an empty channel
        labkey_log_ch = Channel.empty()
    }

    ch_completion = merged_results_for_completion.map { _results -> "STAT+BLAST workflow complete!" }

    emit:
    completion = ch_completion
    labkey_log = labkey_log_ch
}
