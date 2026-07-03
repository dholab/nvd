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

include { GATHER_READS            } from "../subworkflows/gather_reads"
include { PREPROCESS_READS        } from "../subworkflows/preprocess_reads"
include { PREPROCESS_CONTIGS      } from "../subworkflows/preprocess_contigs"
include { EXTRACT_HUMAN_VIRUSES   } from "../subworkflows/extract_human_virus_contigs"
include { CLASSIFY_WITH_MEGABLAST } from "../subworkflows/classify_with_megablast"
include { CLASSIFY_WITH_BLASTN    } from "../subworkflows/classify_with_blastn"
include { METAGENOME_PROFILING    } from "../subworkflows/metagenome_profiling"
include { SAMPLE_SIMILARITY_QC    } from "../subworkflows/sample_similarity_qc"
include { REPORTING               } from "../subworkflows/reporting"
include { COMPUTE_RUN_CONTEXT ; ENSURE_TAXONOMY } from "../modules/utils"


workflow NVD_MAIN {
  take:
  ch_samplesheet

  main:

  // Validate required params. BLAST DB is only required when assembled contigs
  // can reach BLAST classification.
  def requires_blast_db = !(params.skip_assembly || params.skip_blast)
  def target_enrichment_enabled = NvdUtils.targetEnrichmentEnabled(params)
  def has_target_enrichment_index = NvdUtils.hasTargetEnrichmentIndex(params)
  assert (!requires_blast_db || (params.blast_db && file(params.blast_db).isDirectory())) && (!target_enrichment_enabled || has_target_enrichment_index) : """
    One or more required parameters are missing or point to non-existent files:

      blast_db                        -> ${requires_blast_db ? params.blast_db : 'not required when skip_assembly or skip_blast is enabled'}
      target enrichment index source  -> ${target_enrichment_enabled ? 'virus_index / virus_index_url / virus_reference_fasta: at least one must be set' : 'not required when target enrichment is disabled'}

    Please supply the above in your `-c nextflow.config` or via `-params-file`.
    """

  // Reference channels
  ch_blast_db_files = params.blast_db ? channel.fromPath(params.blast_db) : channel.empty()

  // Resolve explicit taxonomy path; Python owns fallback taxonomy-cache rules.
  dirs = NvdDirs.resolve(params, log)

  // Compute run context upfront without writing workflow state.
  COMPUTE_RUN_CONTEXT(channel.fromPath(params.samplesheet))

  ENSURE_TAXONOMY(
    channel.value(dirs.taxonomy_dir)
  )

  GATHER_READS(ch_samplesheet)

  PREPROCESS_READS(GATHER_READS.out.reads)

  ch_sourmash_tax_reports = channel.empty()

  if (params.experimental == true) {
    metagenome_profiling = METAGENOME_PROFILING(PREPROCESS_READS.out.reads)
    SAMPLE_SIMILARITY_QC(metagenome_profiling.query_sketches)
    ch_sourmash_tax_reports = metagenome_profiling.tax_reports
  }

  PREPROCESS_CONTIGS(PREPROCESS_READS.out.reads)

  ch_run_context = COMPUTE_RUN_CONTEXT.out.run_context
  ch_taxonomy_dir = ENSURE_TAXONOMY.out.taxonomy_dir
  EXTRACT_HUMAN_VIRUSES(
    PREPROCESS_CONTIGS.out.contigs,
    PREPROCESS_CONTIGS.out.viral_reads,
    PREPROCESS_READS.out.virus_index,
    PREPROCESS_READS.out.depletion_index,
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
    EXTRACT_HUMAN_VIRUSES.out.contigs,
    EXTRACT_HUMAN_VIRUSES.out.contig_read_counts,
    EXTRACT_HUMAN_VIRUSES.out.filtered_bam,
    PREPROCESS_READS.out.virus_enrichment_stats,
    ch_taxonomy_dir,
    COMPUTE_RUN_CONTEXT.out.ready,
    ch_run_context,
    ch_sourmash_tax_reports,
    workflow.runName,
  )

  emit:
  completion = REPORTING.out.blast_results.count().map { n -> "NVD main workflow complete: ${n} samples processed" }
  labkey_log = REPORTING.out.labkey_log
}
