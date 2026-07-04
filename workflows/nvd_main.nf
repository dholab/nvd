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
include { SHORT_READ_DENOVO_ASSEMBLY } from "../subworkflows/short_read_denovo_assembly"
include { LONG_READ_DENOVO_ENSEMBLY  } from "../subworkflows/long_read_denovo_ensembly"
include { PREPROCESS_CONTIGS      } from "../subworkflows/preprocess_contigs"
include { EXTRACT_HUMAN_VIRUSES   } from "../subworkflows/extract_human_virus_contigs"
include { CLASSIFY_WITH_MEGABLAST } from "../subworkflows/classify_with_megablast"
include { CLASSIFY_WITH_BLASTN    } from "../subworkflows/classify_with_blastn"
include { RAPID_SCREENING         } from "../subworkflows/rapid_screening"
include { SAMPLE_SIMILARITY_QC    } from "../subworkflows/sample_similarity_qc"
include { RAPID_SCREENING_EVAL    } from "../subworkflows/rapid_screening_eval"
include { REPORTING               } from "../subworkflows/reporting"
include { COLLECT_CONTIGS         } from "../modules/contigs"
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

  ch_sourmash_gather_csv = channel.empty()
  ch_sourmash_lineages = channel.empty()
  ch_sourmash_tax_reports = channel.empty()

  if (params.experimental == true) {
    rapid_screening = RAPID_SCREENING(PREPROCESS_READS.out.reads)
    SAMPLE_SIMILARITY_QC(rapid_screening.query_sketches)
    ch_sourmash_gather_csv = rapid_screening.gather_csv
    ch_sourmash_lineages = rapid_screening.lineages
    ch_sourmash_tax_reports = rapid_screening.tax_reports
  }

  // Filter out samples with fewer than 100 reads before assembly; they are
  // insufficient for de novo assembly and would otherwise fan out downstream.
  ch_reads_for_assembly = PREPROCESS_READS.out.reads
    .filter { _id, _platform, _read_structure, _fq -> !params.skip_assembly }
    .map { id, platform, read_structure, fq -> tuple(id, platform, read_structure, fq, file(fq).countFastq()) }
    .filter { _id, _platform, _read_structure, _fq, count -> count >= 100 }
    .map { id, platform, read_structure, fq, _count -> tuple(id, platform, read_structure, file(fq)) }

  ch_reads_by_platform = ch_reads_for_assembly.branch { _id, platform, _read_structure, _fq ->
    illumina: platform == "illumina"
    long_read: true
  }

  SHORT_READ_DENOVO_ASSEMBLY(ch_reads_by_platform.illumina)
  LONG_READ_DENOVO_ENSEMBLY(ch_reads_by_platform.long_read)

  ch_assembled_contigs = SHORT_READ_DENOVO_ASSEMBLY.out.contigs
    .mix(LONG_READ_DENOVO_ENSEMBLY.out.contigs)

  COLLECT_CONTIGS(ch_assembled_contigs)

  PREPROCESS_CONTIGS(COLLECT_CONTIGS.out.fasta)

  ch_run_context = COMPUTE_RUN_CONTEXT.out.run_context
  ch_taxonomy_dir = ENSURE_TAXONOMY.out.taxonomy_dir
  EXTRACT_HUMAN_VIRUSES(
    PREPROCESS_CONTIGS.out.contigs,
    PREPROCESS_READS.out.reads,
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
    COLLECT_CONTIGS.out.lookup,
    EXTRACT_HUMAN_VIRUSES.out.filtered_bam,
    PREPROCESS_READS.out.virus_enrichment_stats,
    ch_taxonomy_dir,
    COMPUTE_RUN_CONTEXT.out.ready,
    ch_run_context,
    ch_sourmash_tax_reports,
    workflow.runName,
  )

  if (params.experimental == true) {
    RAPID_SCREENING_EVAL(
      PREPROCESS_READS.out.read_counts,
      ch_sourmash_gather_csv,
      ch_sourmash_tax_reports,
      ch_sourmash_lineages,
      REPORTING.out.blast_results,
      REPORTING.out.crumbs_taxa,
      REPORTING.out.crumbs_contigs,
      workflow.runName,
    )
  }

  emit:
  completion = REPORTING.out.blast_results.count().map { n -> "NVD main workflow complete: ${n} samples processed" }
  labkey_log = REPORTING.out.labkey_log
}
