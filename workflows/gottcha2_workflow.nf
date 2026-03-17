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
include { SANITIZE_EXTRACTED_FASTA } from "../modules/bbmap"
include { CHECK_RUN_STATE; COMPLETE_RUN } from "../modules/utils"
include { CHECK_MIN_READS } from "../modules/count_reads"


workflow GOTTCHA2_WORKFLOW {
    take:
    ch_sample_fastqs // Queue channel: tuple val(sample_id), val(platform), val(read_structure), path(fastq)

    main:
    // Determine whether GOTTCHA2 was selected by the user.
    def gottcha_selected = NvdUtils.isToolSelected(params, 'gottcha')

    // Validate required params when GOTTCHA2 is selected
    if (gottcha_selected) {
        assert (
            params.gottcha2_db &&
            file("${params.gottcha2_db}.mmi").exists() &&
            file("${params.gottcha2_db}.stats").exists() &&
            file("${params.gottcha2_db}.tax.tsv").exists()
        ) :
        """
        One or more required GOTTCHA2 database files are missing or gottcha2_db is not set:

          gottcha2_db          -> ${params.gottcha2_db}
          gottcha2_db.mmi      -> ${params.gottcha2_db ? params.gottcha2_db + '.mmi' : 'N/A'}
          gottcha2_db.stats    -> ${params.gottcha2_db ? params.gottcha2_db + '.stats' : 'N/A'}
          gottcha2_db.tax.tsv  -> ${params.gottcha2_db ? params.gottcha2_db + '.tax.tsv' : 'N/A'}

        Please supply --gottcha2_db in your config or command line, pointing to the
        GOTTCHA2 database prefix (without extension). The .mmi, .stats, and .tax.tsv
        companion files must exist at that prefix.
        """

        NvdUtils.validateLabkeyGottcha2(params)
    }

    // Signal channel: emits true when gottcha is selected, empty otherwise.
    // This single signal gates both state-management processes (CHECK_RUN_STATE,
    // COMPLETE_RUN) and data-processing processes (profiling, FASTA generation).
    // When empty, combine() produces nothing, so downstream processes remain in
    // the DAG but never execute.
    ch_gottcha_enabled = gottcha_selected
        ? Channel.value(true)
        : Channel.empty()

    // Resolve directories for stateful vs stateless mode.
    // Returns absolute path strings (not Nextflow path objects) to avoid staging.
    dirs = NvdDirs.resolve(params, log)

    // Gate CHECK_RUN_STATE: combine with the tool-selection signal so the input
    // tuple only emits when gottcha is selected. When empty, CHECK_RUN_STATE
    // remains in the DAG but never executes.
    ch_run_state_input = ch_gottcha_enabled
        .combine(Channel.fromPath(params.samplesheet))
        .combine(Channel.value(dirs.state_dir))
        .combine(Channel.value(dirs.taxonomy_dir))
        .map { _flag, samplesheet, state_dir, taxonomy_dir ->
            tuple(samplesheet, state_dir, taxonomy_dir)
        }

    // Check run state upfront (prevents duplicate processing of same sample set)
    CHECK_RUN_STATE(ch_run_state_input, "gottcha2,gottcha2_fasta")

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

    // Gate sample channel: when gottcha is not selected, ch_gottcha_enabled is
    // empty, so the combine produces nothing and no downstream operations execute.
    ch_gated_samples = ch_gottcha_enabled
        .combine(ch_sample_fastqs)
        .map { _flag, sample_id, platform, read_structure, fastq ->
            tuple(sample_id, platform, read_structure, fastq)
        }

    // LabKey validation: no-input subworkflow, gated by tool selection
    if (params.labkey && gottcha_selected) {
        VALIDATE_LK_GOTTCHA2()
    }

    // Check minimum read threshold on a compute node instead of using
    // countFastq() on the head node. countFastq() decompresses entire gzipped
    // FASTQs in the Nextflow process, blocking indefinitely for large files on
    // networked storage. CHECK_MIN_READS only reads up to min_reads records
    // before stopping — O(threshold) instead of O(file_size).
    CHECK_MIN_READS(ch_gated_samples)

    // Join threshold check back to samples and filter out those below minimum
    ch_passing_samples = ch_gated_samples
        .map { sample_id, platform, _read_structure, fastq -> tuple(sample_id, platform, fastq) }
        .join(CHECK_MIN_READS.out.counts)
        .filter { _sample_id, _platform, _fastq, read_count -> (read_count as Integer) >= params.min_gottcha_reads }

    ch_nanopore_fastqs = ch_passing_samples
        .filter { _sample_id, platform, _fastq, _read_count -> platform == "nanopore" || platform == "ont" }
        .map { sample_id, _platform, fastq, _read_count -> tuple(sample_id, file(fastq)) }

    ch_illumina_fastqs = ch_passing_samples
        .filter { _sample_id, platform, _fastq, _read_count -> platform == "illumina" }
        .map { sample_id, _platform, fastq, _read_count -> tuple(sample_id, file(fastq)) }

    GOTTCHA2_PROFILE_NANOPORE(
        ch_nanopore_fastqs.combine(ch_gottcha2_db)
    )

    GOTTCHA2_PROFILE_ILLUMINA(
        ch_illumina_fastqs.combine(ch_gottcha2_db)
    )

    GENERATE_FASTA(
        GOTTCHA2_PROFILE_NANOPORE.out.aligned.mix(GOTTCHA2_PROFILE_ILLUMINA.out.aligned)
    )

    ch_extracted_fasta = GENERATE_FASTA.out
        .map { id, fasta, full_tsv, _log, _lineages -> tuple(id, fasta, full_tsv) }

    SANITIZE_EXTRACTED_FASTA(ch_extracted_fasta)

    // Register GOTTCHA2 hits with idempotent keys in the state database.
    // Join sanitized FASTA with run context and hits_dir.
    ch_register_gottcha2_input = SANITIZE_EXTRACTED_FASTA.out
        .combine(ch_run_context)
        .combine(ch_hits_dir)
        .map { sample_id, fasta, full_tsv, sample_set_id, state_dir, hits_dir ->
            tuple(sample_id, fasta, full_tsv, sample_set_id, state_dir, hits_dir)
        }

    REGISTER_GOTTCHA2_HITS(ch_register_gottcha2_input)

    if (params.labkey && gottcha_selected) {
        // Build validation gate: both state check and LabKey validation must pass
        // before any uploads proceed. Matches stat_blast_workflow pattern.
        ch_validation_gate = CHECK_RUN_STATE.out.ready
            .combine(VALIDATE_LK_GOTTCHA2.out.validated)
            .map { _ready, _validated -> true }
            .first()

        // LabKey upload channels
        ch_gottcha2_results_for_labkey = GOTTCHA2_PROFILE_NANOPORE.out.full_tsv
            .mix(GOTTCHA2_PROFILE_ILLUMINA.out.full_tsv)

        ch_extracted_for_labkey = SANITIZE_EXTRACTED_FASTA.out

        BUNDLE_GOTTCHA2_FOR_LABKEY(
            ch_gottcha2_results_for_labkey,
            ch_extracted_for_labkey,
            params.experiment_id,
            workflow.runName,
            ch_validation_gate,
            ch_run_context
        )

        // Gate run completion on LabKey bundle completion
        ch_run_complete_gate = BUNDLE_GOTTCHA2_FOR_LABKEY.out.upload_log.collect()
            .map { _logs -> true }

        // Terminal channel for LabKey path
        ch_terminal = BUNDLE_GOTTCHA2_FOR_LABKEY.out.upload_log
    } else {
        // Gate run completion on REGISTER_GOTTCHA2_HITS completion (all samples)
        ch_run_complete_gate = REGISTER_GOTTCHA2_HITS.out.collect()
            .map { _logs -> true }

        // Terminal channel for non-LabKey path
        ch_terminal = REGISTER_GOTTCHA2_HITS.out
    }

    // Mark the run as completed in the state database
    COMPLETE_RUN(
        ch_run_complete_gate,
        ch_state_dir,
        "completed"
    )

    // Completion signal: always emits exactly one value regardless of whether
    // data flowed through the workflow.
    //
    // When gottcha is not selected: Channel.value() emits one token immediately.
    // When gottcha is selected: count() waits for ch_terminal to close, then
    // emits the number of items that flowed through (including 0 if all samples
    // were filtered out). Either way, exactly one emission.
    ch_completion = gottcha_selected
        ? ch_terminal.count().map { n -> "GOTTCHA2 complete: ${n} samples processed" }
        : Channel.value("GOTTCHA2 skipped: tool not selected")

    emit:
    completion = ch_completion
}
