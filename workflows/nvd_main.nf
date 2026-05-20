/*
 * NVD main workflow
 *
 * Human virus detection pipeline: deacon-based virus read extraction,
 * preprocessing, SPAdes assembly, and two-phase BLAST verification.
 *
 * Architecture: deacon virus extraction runs first on raw R1/R2 reads
 * (outputting interleaved), then preprocessing (dedup, trim, optional host depletion, filter)
 * operates on the tiny virus subset, then SPAdes and BLAST.
 */

nextflow.enable.dsl = 2

include { GATHER_READS            } from "./gather_reads"
include { PREPROCESS_READS        } from "../subworkflows/preprocess_reads"
include { PREPROCESS_CONTIGS      } from "../subworkflows/preprocess_contigs"
include { EXTRACT_HUMAN_VIRUSES   } from "../subworkflows/extract_human_virus_contigs"
include { CLASSIFY_WITH_MEGABLAST } from "../subworkflows/classify_with_megablast"
include { CLASSIFY_WITH_BLASTN    } from "../subworkflows/classify_with_blastn"
include { REPORTING               } from "../subworkflows/reporting"
include { BUNDLE_BLAST_FOR_LABKEY } from "../subworkflows/bundle_blast_for_labkey"
include { COMPUTE_RUN_CONTEXT ; ENSURE_TAXONOMY ; NOTIFY_SLACK } from "../modules/utils"
include { VALIDATE_LK_BLAST       } from "../subworkflows/validate_lk_blast_lists.nf"
include { VALIDATE_LK_EXP_FRESH   } from "../modules/validate_blast_labkey.nf"
include { REGISTER_LK_EXPERIMENT  } from "../modules/validate_blast_labkey.nf"


workflow NVD_MAIN {
  take:
  ch_samplesheet_row

  main:

  // Validate required params
  assert params.blast_db && file(params.blast_db).isDirectory() && (params.virus_index || params.virus_index_url || params.virus_reference_fasta) : """
    One or more required parameters are missing or point to non-existent files:

      blast_db                        -> ${params.blast_db}
      virus_index / virus_index_url / virus_reference_fasta: at least one must be set

    Please supply the above in your `-c nextflow.config` or via `-params-file`.
    """

  if (params.labkey) {
    NvdUtils.validateLabkeyBlast(params)
  }

  // Reference channels
  ch_blast_db_files = channel.fromPath(params.blast_db)

  // Resolve explicit taxonomy path; Python owns fallback taxonomy-cache rules.
  dirs = NvdDirs.resolve(params, log)

  // Compute run context upfront without writing workflow state.
  COMPUTE_RUN_CONTEXT(channel.fromPath(params.samplesheet))

  ENSURE_TAXONOMY(
    channel.value(dirs.taxonomy_dir)
  )

  if (params.labkey) {
    VALIDATE_LK_BLAST()
  }

  GATHER_READS(ch_samplesheet_row)

  PREPROCESS_READS(GATHER_READS.out)

  PREPROCESS_CONTIGS(PREPROCESS_READS.out.reads)

  ch_run_context = COMPUTE_RUN_CONTEXT.out.run_context.first()
  ch_taxonomy_dir = ENSURE_TAXONOMY.out.taxonomy_dir.first()
  EXTRACT_HUMAN_VIRUSES(
    PREPROCESS_CONTIGS.out.contigs,
    PREPROCESS_CONTIGS.out.viral_reads,
    PREPROCESS_READS.out.virus_index,
  )

  CLASSIFY_WITH_MEGABLAST(
    EXTRACT_HUMAN_VIRUSES.out.contigs,
    ch_blast_db_files,
    ch_taxonomy_dir,
  )

  CLASSIFY_WITH_BLASTN(
    CLASSIFY_WITH_MEGABLAST.out.filtered_megablast,
    CLASSIFY_WITH_MEGABLAST.out.megablast_contigs,
    ch_blast_db_files,
    ch_taxonomy_dir,
  )

  REPORTING(
    CLASSIFY_WITH_BLASTN.out.merged_results,
    PREPROCESS_READS.out.read_counts,
  )

  if (params.labkey) {
    VALIDATE_LK_EXP_FRESH(
      REPORTING.out.blast_results.first()
    )

    ch_validation_gate = COMPUTE_RUN_CONTEXT.out.ready
      .combine(VALIDATE_LK_EXP_FRESH.out.validated)
      .map { _ready, _validated -> true }
      .first()

    BUNDLE_BLAST_FOR_LABKEY(
      REPORTING.out.blast_results,
      EXTRACT_HUMAN_VIRUSES.out.contigs,
      PREPROCESS_READS.out.read_counts,
      params.experiment_id,
      workflow.runName,
      EXTRACT_HUMAN_VIRUSES.out.contig_read_counts,
      ch_validation_gate,
    )

    REGISTER_LK_EXPERIMENT(
      BUNDLE_BLAST_FOR_LABKEY.out.upload_log.collect()
    )

    labkey_log_ch = BUNDLE_BLAST_FOR_LABKEY.out.upload_log

    ch_run_complete_gate = REGISTER_LK_EXPERIMENT.out
      .collect()
      .map { _logs -> true }

    ch_terminal = REGISTER_LK_EXPERIMENT.out
  }
  else {
    labkey_log_ch = channel.empty()

    ch_run_complete_gate = REPORTING.out.blast_results
      .collect()
      .map { _logs -> true }

    ch_terminal = REPORTING.out.blast_results
  }

  if (params.slack_enabled && params.slack_channel && params.labkey) {
    ch_labkey_url = channel.value(
      "https://${params.labkey_server}/${params.labkey_project_name}/list-grid.view?name=${params.labkey_blast_meta_hits_list}"
    )

    NOTIFY_SLACK(
      ch_run_complete_gate,
      ch_run_context,
      ch_labkey_url,
    )
  }

  emit:
  completion = ch_terminal.count().map { n -> "NVD main workflow complete: ${n} samples processed" }
  labkey_log = labkey_log_ch
}
