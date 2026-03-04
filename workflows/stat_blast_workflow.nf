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
include { CHECK_RUN_STATE; REGISTER_HITS; COMPLETE_RUN; NOTIFY_SLACK } from "../modules/utils"
include { VALIDATE_LK_BLAST } from "../subworkflows/validate_lk_blast_lists.nf"
include { VALIDATE_LK_EXP_FRESH } from "../modules/validate_blast_labkey.nf"
include { REGISTER_LK_EXPERIMENT } from "../modules/validate_blast_labkey.nf"


workflow STAT_BLAST_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel: tuple val(sample_id), val(platform), val(read_structure), path(fastq)

    main:
    // Determine whether STAT+BLAST was selected by the user.
    def blast_selected = NvdUtils.isToolSelected(params, 'blast')

    // Validate required params when STAT+BLAST is selected
    if (blast_selected) {
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

    // Signal channel: emits true when STAT+BLAST is selected, empty otherwise.
    // This single signal gates both state-management processes (CHECK_RUN_STATE,
    // COMPLETE_RUN) and data-processing processes (COUNT_READS, PREPROCESS_CONTIGS).
    // When empty, combine() produces nothing, so downstream processes remain in
    // the DAG but never execute.
    ch_blast_enabled = blast_selected
        ? Channel.value(true)
        : Channel.empty()

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

    // Resolve directories for stateful vs stateless mode.
    // Returns absolute path strings (not Nextflow path objects) to avoid staging.
    // See lib/NvdDirs.groovy for validation logic and error messages.
    dirs = NvdDirs.resolve(params, log)

    // Gate CHECK_RUN_STATE: combine with the tool-selection signal so the input
    // tuple only emits when STAT+BLAST is selected. When empty, CHECK_RUN_STATE
    // remains in the DAG but never executes.
    ch_run_state_input = ch_blast_enabled
        .combine(Channel.fromPath(params.samplesheet))
        .combine(Channel.value(dirs.state_dir))
        .combine(Channel.value(dirs.taxonomy_dir))
        .map { _flag, samplesheet, state_dir, taxonomy_dir ->
            tuple(samplesheet, state_dir, taxonomy_dir)
        }

    // Check run state upfront (prevents duplicate processing of same sample set)
    CHECK_RUN_STATE(ch_run_state_input, "blast,blast_fasta")

    if (params.labkey) {
        VALIDATE_LK_BLAST()
    }

    // Gate sample channel: when STAT+BLAST is not selected, ch_blast_enabled is
    // empty, so the combine produces nothing and no downstream operations execute.
    ch_gated_samples = ch_blast_enabled
        .combine(ch_sample_fastqs)
        .map { _flag, sample_id, platform, read_structure, fastq ->
            tuple(sample_id, platform, read_structure, fastq)
        }

    // Count reads for each sample
    COUNT_READS(ch_gated_samples)

    PREPROCESS_CONTIGS(
        ch_gated_samples,
        ch_stat_dbss,
        ch_stat_annotation,
        ch_human_virus_taxlist
    )

    // Convert run_context from queue channel to value channel so it can be
    // consumed by multiple processes (hit registration, LabKey uploads, etc.).
    // This is safe because CHECK_RUN_STATE emits exactly once per pipeline run.
    ch_run_context = CHECK_RUN_STATE.out.run_context.first()

    // Extract state_dir from run_context for state-tracking processes
    // (run_context is a value channel: [sample_set_id, state_dir])
    ch_state_dir = ch_run_context.map { _sample_set_id, state_dir -> state_dir }
    
    // Create value channels for taxonomy_dir and hits_dir
    // These are resolved at workflow start and broadcast to all consumers
    ch_taxonomy_dir = Channel.value(dirs.taxonomy_dir)
    ch_hits_dir = Channel.value(dirs.hits_dir)

    EXTRACT_HUMAN_VIRUSES(
        PREPROCESS_CONTIGS.out.contigs,
        PREPROCESS_CONTIGS.out.viral_reads,
        ch_stat_index,
        ch_stat_dbss,
        ch_stat_annotation,
        ch_state_dir,
        ch_taxonomy_dir
    )

    CLASSIFY_WITH_MEGABLAST(
        EXTRACT_HUMAN_VIRUSES.out.contigs,
        ch_blast_db_files,
        ch_state_dir,
        ch_taxonomy_dir
    )

    CLASSIFY_WITH_BLASTN(
        CLASSIFY_WITH_MEGABLAST.out.filtered_megablast,
        CLASSIFY_WITH_MEGABLAST.out.megablast_contigs,
        ch_blast_db_files,
        ch_state_dir,
        ch_taxonomy_dir
    )

    // Register hits with idempotent keys in the state database.
    // Join contigs with merged BLAST results, then combine with run context and hits_dir.
    // DSL2 automatically forks EXTRACT_HUMAN_VIRUSES.out.contigs and 
    // CLASSIFY_WITH_BLASTN.out.merged_results for multiple consumers.
    ch_register_hits_input = EXTRACT_HUMAN_VIRUSES.out.contigs
        .join(CLASSIFY_WITH_BLASTN.out.merged_results, by: 0)
        .combine(ch_run_context)
        .combine(ch_hits_dir)
        .map { sample_id, contigs, blast_results, sample_set_id, state_dir, hits_dir ->
            tuple(sample_id, contigs, blast_results, sample_set_id, state_dir, hits_dir)
        }

    REGISTER_HITS(ch_register_hits_input)
    // Hits are registered in the state database for querying via `nvd hits export`.
    // Collect output to gate run completion when LabKey is disabled.

    if (params.labkey) {
        // Check LabKey guard list - ensures experiment_id hasn't been uploaded before
        // Uses .first() to get a value channel trigger
        VALIDATE_LK_EXP_FRESH(
            CLASSIFY_WITH_BLASTN.out.merged_results.first()
        )

        // Bundle and upload - gates ensure both checks passed:
        // 1. CHECK_RUN_STATE.out.ready: local sample_set_id check passed
        // 2. VALIDATE_LK_EXP_FRESH.out.validated: LabKey guard list check passed
        // We combine both gates to ensure both validations complete before uploading.
        // The .first() converts the result to a value channel so it broadcasts to all
        // consumers instead of being consumed by the first process that reads it.
        ch_validation_gate = CHECK_RUN_STATE.out.ready
            .combine(VALIDATE_LK_EXP_FRESH.out.validated)
            .map { _ready, _validated -> true }
            .first()

        BUNDLE_BLAST_FOR_LABKEY(
            CLASSIFY_WITH_BLASTN.out.merged_results,
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

        // Gate run completion on LabKey experiment registration
        // This ensures all uploads are complete before marking the run done
        ch_run_complete_gate = REGISTER_LK_EXPERIMENT.out.collect()
            .map { _logs -> true }
    } else {
        // If no labkey upload then just pass an empty channel
        labkey_log_ch = Channel.empty()

        // Gate run completion on REGISTER_HITS completion (all samples)
        ch_run_complete_gate = REGISTER_HITS.out.collect()
            .map { _logs -> true }
    }

    // Mark the run as completed in the state database
    COMPLETE_RUN(
        ch_run_complete_gate,
        ch_state_dir,
        "completed"
    )

    // Send Slack notification (if enabled)
    // Only runs when slack_enabled=true, slack_channel is set, and LabKey is enabled
    // (we need LabKey for the results URL in the notification)
    // Stats are queried from state database by the Python script using sample_set_id
    if (params.slack_enabled && params.slack_channel && params.labkey) {
        // Build LabKey URL to the hits list view
        // Format: https://{server}/{project}/list-grid.view?name={list_name}
        ch_labkey_url = Channel.value(
            "https://${params.labkey_server}/${params.labkey_project_name}/list-grid.view?name=${params.labkey_blast_meta_hits_list}"
        )

        NOTIFY_SLACK(
            COMPLETE_RUN.out.done,
            ch_run_context,  // tuple: (sample_set_id, state_dir)
            ch_labkey_url
        )
    }

    // Terminal channel: the last meaningful output in each execution path.
    ch_terminal = params.labkey
        ? REGISTER_LK_EXPERIMENT.out
        : REGISTER_HITS.out

    // Completion signal: always emits exactly one value regardless of whether
    // data flowed through the workflow.
    //
    // When STAT+BLAST is not selected: Channel.value() emits one token immediately.
    // When STAT+BLAST is selected: count() waits for ch_terminal to close, then
    // emits the number of items that flowed through (including 0 if all samples
    // were filtered out). Either way, exactly one emission.
    ch_completion = blast_selected
        ? ch_terminal.count().map { n -> "STAT+BLAST complete: ${n} samples processed" }
        : Channel.value("STAT+BLAST skipped: tool not selected")

    emit:
    completion = ch_completion
    labkey_log = labkey_log_ch
}
