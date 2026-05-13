/*
 * STAT+BLAST Workflow
 *
 * Human virus detection pipeline: deacon-based virus read extraction,
 * preprocessing, SPAdes assembly, and two-phase BLAST verification.
 *
 * Architecture: deacon virus extraction runs first on raw R1/R2 reads
 * (outputting interleaved), then preprocessing (dedup, trim, filter)
 * operates on the tiny virus subset, then SPAdes and BLAST.
 */

nextflow.enable.dsl=2

include { DEACON_BUILD_INDEX_FROM_STAT_K_MERS ; DEACON_FILTER_HUMAN_VIRUS_READS } from "../modules/deacon"
include { DEACON_FETCH_INDEX ; DEACON_DEPLETE } from "../modules/deacon"
include { DEDUP_WITH_CLUMPIFY ; TRIM_ADAPTERS ; FILTER_READS ; REPAIR_PAIRS } from "../modules/bbmap"
include { PREPROCESS_CONTIGS } from "../subworkflows/preprocess_contigs"
include { EXTRACT_HUMAN_VIRUSES } from "../subworkflows/extract_human_virus_contigs"
include { CLASSIFY_WITH_MEGABLAST } from "../subworkflows/classify_with_megablast"
include { CLASSIFY_WITH_BLASTN } from "../subworkflows/classify_with_blastn"
include { BUNDLE_BLAST_FOR_LABKEY } from "../subworkflows/bundle_blast_for_labkey"
include { CHECK_RUN_STATE; REGISTER_HITS; COMPLETE_RUN; NOTIFY_SLACK } from "../modules/utils"
include { VALIDATE_LK_BLAST } from "../subworkflows/validate_lk_blast_lists.nf"
include { VALIDATE_LK_EXP_FRESH } from "../modules/validate_blast_labkey.nf"
include { REGISTER_LK_EXPERIMENT } from "../modules/validate_blast_labkey.nf"


workflow STAT_BLAST_WORKFLOW {
    take:
    ch_sample_fastqs  // Queue channel: pre-interleave tuples from GATHER_READS
                      // Paired: tuple(sample_id, platform, R1, R2)
                      // Single: tuple(sample_id, platform, fastq)

    main:

    // Validate required params
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

    if (params.labkey) {
        NvdUtils.validateLabkeyBlast(params)
    }

    // Reference channels
    ch_blast_db_files = Channel.fromPath(params.blast_db)
    ch_stat_index     = Channel.fromPath(params.stat_index)
    ch_stat_dbss      = Channel.fromPath(params.stat_dbss)
    ch_stat_annotation     = Channel.fromPath(params.stat_annotation)
    ch_human_virus_taxlist = Channel.fromPath(params.human_virus_taxlist)

    // Resolve directories for stateful vs stateless mode.
    dirs = NvdDirs.resolve(params, log)

    // Check run state upfront (prevents duplicate processing of same sample set)
    ch_run_state_input = Channel.fromPath(params.samplesheet)
        .combine(Channel.value(dirs.state_dir))
        .combine(Channel.value(dirs.taxonomy_dir))

    CHECK_RUN_STATE(ch_run_state_input, "blast,blast_fasta")

    if (params.labkey) {
        VALIDATE_LK_BLAST()
    }

    // Normalize mixed-size tuples from GATHER_READS:
    //   Paired: (id, platform, R1, R2) -> (id, platform, R1, R2)
    //   Single: (id, platform, fastq)  -> (id, platform, fastq, NO_R2)
    // The sentinel file NO_R2 lets DEACON_FILTER_HUMAN_VIRUS_READS distinguish
    // paired from single-end input with a fixed-size tuple.
    ch_normalized_samples = ch_sample_fastqs
        .map { items ->
            if (items.size() == 4)
                tuple(items[0], items[1], file(items[2]), file(items[3]))
            else
                tuple(items[0], items[1], file(items[2]), file("NO_R2"))
        }

    // -------------------------------------------------------------------------
    // Step 1: Frontloaded deacon virus extraction
    // -------------------------------------------------------------------------
    // Build deacon virus index from STAT k-mers (runs once, not per-sample).
    DEACON_BUILD_INDEX_FROM_STAT_K_MERS(
        ch_stat_dbss,
        ch_stat_annotation,
        ch_human_virus_taxlist
    )

    // Extract virus reads — runs BEFORE any preprocessing. For paired reads,
    // deacon takes R1/R2 and outputs interleaved FASTQ in one step.
    DEACON_FILTER_HUMAN_VIRUS_READS(
        ch_normalized_samples.combine(DEACON_BUILD_INDEX_FROM_STAT_K_MERS.out.index)
    )

    // -------------------------------------------------------------------------
    // Step 2: Inlined preprocessing on virus-only reads
    // -------------------------------------------------------------------------
    ch_virus_reads = DEACON_FILTER_HUMAN_VIRUS_READS.out.reads

    // Extract total read counts from deacon summary JSON (replaces COUNT_READS).
    // The seqs_in field is the total input read count across R1+R2.
    ch_read_counts = DEACON_FILTER_HUMAN_VIRUS_READS.out.stats
        .map { sample_id, json_file ->
            def summary = new groovy.json.JsonSlurper().parse(json_file.toFile())
            tuple(sample_id, summary.seqs_in.toString())
        }

    // 2a. Dedup
    ch_after_dedup = params.dedup
        ? DEDUP_WITH_CLUMPIFY(ch_virus_reads)
        : ch_virus_reads

    // 2b. Adapter trim (Illumina only)
    ch_branched_for_trim = ch_after_dedup.branch { _id, platform, _rs, _reads ->
        illumina: platform == "illumina"
        other: true
    }
    ch_after_trim = params.trim_adapters
        ? TRIM_ADAPTERS(ch_branched_for_trim.illumina).mix(ch_branched_for_trim.other)
        : ch_after_dedup

    // 2c. Host scrub with deacon (optional)
    def has_deacon_host_config = params.deacon_index || params.deacon_index_url || params.deacon_contaminants_fasta
    if (params.scrub_host_reads && has_deacon_host_config) {
        ch_local_host_index = params.deacon_index
            ? Channel.fromPath(params.deacon_index)
            : Channel.empty()
        ch_host_fetch_url = (!params.deacon_index && params.deacon_index_url)
            ? Channel.of(params.deacon_index_url)
            : Channel.empty()
        DEACON_FETCH_INDEX(ch_host_fetch_url)
        ch_host_index = ch_local_host_index.mix(DEACON_FETCH_INDEX.out.index)
        ch_after_scrub = DEACON_DEPLETE(ch_after_trim.combine(ch_host_index)).reads
    } else {
        ch_after_scrub = ch_after_trim
    }

    // 2d. Quality/length filter
    ch_with_quality = ch_after_scrub.map { sample_id, platform, read_structure, reads ->
        def min_qual = platform == "illumina"
            ? params.min_read_quality_illumina
            : params.min_read_quality_nanopore
        tuple(sample_id, platform, read_structure, reads, min_qual)
    }
    ch_after_filter = params.filter_reads
        ? FILTER_READS(ch_with_quality)
        : ch_after_scrub

    // 2e. Repair pairs (interleaved only)
    ch_branched_for_repair = ch_after_filter.branch { _id, _p, read_structure, _r ->
        interleaved: read_structure == "interleaved"
        other: true
    }
    ch_repaired = params.repair_pairs
        ? REPAIR_PAIRS(ch_branched_for_repair.interleaved)
        : ch_branched_for_repair.interleaved
    ch_preprocessed = ch_repaired.mix(ch_branched_for_repair.other)

    // -------------------------------------------------------------------------
    // Step 3: Assembly and classification
    // -------------------------------------------------------------------------
    PREPROCESS_CONTIGS(ch_preprocessed)

    ch_run_context = CHECK_RUN_STATE.out.run_context.first()
    ch_state_dir = ch_run_context.map { _sample_set_id, state_dir -> state_dir }
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

    // Register hits
    ch_register_hits_input = EXTRACT_HUMAN_VIRUSES.out.contigs
        .join(CLASSIFY_WITH_BLASTN.out.merged_results, by: 0)
        .combine(ch_run_context)
        .combine(ch_hits_dir)
        .map { sample_id, contigs, blast_results, sample_set_id, state_dir, hits_dir ->
            tuple(sample_id, contigs, blast_results, sample_set_id, state_dir, hits_dir)
        }

    REGISTER_HITS(ch_register_hits_input)

    if (params.labkey) {
        VALIDATE_LK_EXP_FRESH(
            CLASSIFY_WITH_BLASTN.out.merged_results.first()
        )

        ch_validation_gate = CHECK_RUN_STATE.out.ready
            .combine(VALIDATE_LK_EXP_FRESH.out.validated)
            .map { _ready, _validated -> true }
            .first()

        BUNDLE_BLAST_FOR_LABKEY(
            CLASSIFY_WITH_BLASTN.out.merged_results,
            EXTRACT_HUMAN_VIRUSES.out.contigs,
            ch_read_counts,
            params.experiment_id,
            workflow.runName,
            EXTRACT_HUMAN_VIRUSES.out.contig_read_counts,
            ch_validation_gate,
            ch_run_context
        )

        REGISTER_LK_EXPERIMENT(
            BUNDLE_BLAST_FOR_LABKEY.out.upload_log.collect()
        )

        labkey_log_ch = BUNDLE_BLAST_FOR_LABKEY.out.upload_log

        ch_run_complete_gate = REGISTER_LK_EXPERIMENT.out.collect()
            .map { _logs -> true }

        ch_terminal = REGISTER_LK_EXPERIMENT.out
    } else {
        labkey_log_ch = Channel.empty()

        ch_run_complete_gate = REGISTER_HITS.out.collect()
            .map { _logs -> true }

        ch_terminal = REGISTER_HITS.out
    }

    COMPLETE_RUN(
        ch_run_complete_gate,
        ch_state_dir,
        "completed"
    )

    if (params.slack_enabled && params.slack_channel && params.labkey) {
        ch_labkey_url = Channel.value(
            "https://${params.labkey_server}/${params.labkey_project_name}/list-grid.view?name=${params.labkey_blast_meta_hits_list}"
        )

        NOTIFY_SLACK(
            COMPLETE_RUN.out.done,
            ch_run_context,
            ch_labkey_url
        )
    }

    emit:
    completion = ch_terminal.count().map { n -> "STAT+BLAST complete: ${n} samples processed" }
    labkey_log = labkey_log_ch
}
