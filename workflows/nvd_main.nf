/*
 * NVD main workflow
 *
 * Human virus detection pipeline: deacon-based target read enrichment,
 * preprocessing, SPAdes assembly, and two-phase BLAST verification.
 *
 * Architecture: deacon target enrichment runs first on raw R1/R2 reads
 * (outputting interleaved), then preprocessing (dedup, trim, optional host depletion, filter)
 * operates on the enriched subset, then SPAdes and BLAST.
 */

nextflow.enable.dsl = 2

include { GATHER_READS            } from "../subworkflows/gather_reads"
include { PREPROCESS_READS        } from "../subworkflows/preprocess_reads"
include { SHORT_READ_DENOVO_ASSEMBLY } from "../subworkflows/short_read_denovo_assembly"
include { LONG_READ_DENOVO_ENSEMBLY  } from "../subworkflows/long_read_denovo_ensembly"
include { PROCESS_CONTIGS         } from "../subworkflows/process_contigs"
include { PREPARE_BLAST_QUERIES } from "../subworkflows/prepare_blast_queries"
include { CLASSIFY_WITH_MEGABLAST } from "../subworkflows/classify_with_megablast"
include { CLASSIFY_WITH_BLASTN    } from "../subworkflows/classify_with_blastn"
include { RAPID_SCREENING         } from "../subworkflows/rapid_screening"
include { SAMPLE_SIMILARITY_QC    } from "../subworkflows/sample_similarity_qc"
include { RAPID_SCREENING_EVAL    } from "../subworkflows/rapid_screening_eval"
include { MULTIQC_BUNDLING        } from "../subworkflows/multiqc_bundling"
include { REPORTING               } from "../subworkflows/reporting"
include { COMPUTE_RUN_CONTEXT ; ENSURE_TAXONOMY } from "../modules/utils"


workflow NVD_MAIN {
  take:
  ch_samplesheet
  ch_nvd_version_file

  main:

  // BLAST needs its database when contig queries may run or experimental
  // direct-read querying remains enabled while assembly is skipped.
  def requires_blast_db = !params.skip_blast && (!params.skip_assembly || params.experimental)
  def target_enrichment_enabled = NvdUtils.targetEnrichmentEnabled(params)
  def depletion_enabled = NvdUtils.depletionEnabled(params)
  def has_target_enrichment_index = NvdUtils.hasTargetEnrichmentIndex(params)
  assert (!requires_blast_db || (params.blast_db && file(params.blast_db).isDirectory())) && (!target_enrichment_enabled || has_target_enrichment_index) : """
    One or more required parameters are missing or point to non-existent files:

      blast_db                        -> ${requires_blast_db ? params.blast_db : 'not required by the enabled assembly/BLAST routes'}
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
  ch_risk_group_lookup = Channel.value(file("${projectDir}/assets/human_virus_risk_group_lookup.tsv"))

  if (params.experimental == true) {
    rapid_screening = RAPID_SCREENING(
      PREPROCESS_READS.out.profiled_batches_by_sample,
      ch_risk_group_lookup,
    )
    SAMPLE_SIMILARITY_QC(rapid_screening.query_sketches)
    ch_sourmash_gather_csv = rapid_screening.gather_csv
    ch_sourmash_lineages = rapid_screening.lineages
    ch_sourmash_tax_reports = rapid_screening.tax_reports
  }

  // Short reads retain their minimum-count gate; experimental long-read
  // assemblers use tool-specific length evidence from each cached profile.
  ch_short_read_batches_by_sample_for_assembly = PREPROCESS_READS.out.profiled_batches_by_sample
    .filter { meta, _batches -> !params.skip_assembly && meta.platform == "illumina" }

  ch_long_read_profiles_for_assembly = PREPROCESS_READS.out.profiled_read_batches
    .filter { meta, _reads, _profile_json, _length_histogram -> !params.skip_assembly && meta.platform != "illumina" }

  SHORT_READ_DENOVO_ASSEMBLY(ch_short_read_batches_by_sample_for_assembly)
  LONG_READ_DENOVO_ENSEMBLY(ch_long_read_profiles_for_assembly)

  ch_assembly_disabled = params.skip_assembly
    ? PREPROCESS_READS.out.profiled_batches_by_sample.map { meta, _batches ->
        log.debug "nvd.contig_route sample_id=${meta.id} platform=${meta.platform} outcome=no_contigs stage=assembly_disabled"
        tuple(meta.id, meta.platform)
      }
    : channel.empty()

  ch_assembled_contigs = SHORT_READ_DENOVO_ASSEMBLY.out.contigs
    .mix(LONG_READ_DENOVO_ENSEMBLY.out.contigs)

  PROCESS_CONTIGS(ch_assembled_contigs)

  ch_no_contig_samples = SHORT_READ_DENOVO_ASSEMBLY.out.no_contigs
    .mix(LONG_READ_DENOVO_ENSEMBLY.out.no_contigs)
    .mix(PROCESS_CONTIGS.out.no_contigs)
    .mix(ch_assembly_disabled)

  ch_run_context = COMPUTE_RUN_CONTEXT.out.run_context
  ch_taxonomy_dir = ENSURE_TAXONOMY.out.taxonomy_dir

  PREPARE_BLAST_QUERIES(
    PROCESS_CONTIGS.out.contigs,
    ch_no_contig_samples,
    PREPROCESS_READS.out.paired_reads_for_mapback,
    PREPROCESS_READS.out.single_reads_for_mapback,
    PREPROCESS_READS.out.target_index,
    PREPROCESS_READS.out.depletion_index,
  )

  CLASSIFY_WITH_MEGABLAST(
    PREPARE_BLAST_QUERIES.out.queries,
    ch_blast_db_files,
    ch_taxonomy_dir,
  )

  CLASSIFY_WITH_BLASTN(
    CLASSIFY_WITH_MEGABLAST.out.filtered_megablast,
    CLASSIFY_WITH_MEGABLAST.out.megablast_query_partition,
    ch_blast_db_files,
    ch_taxonomy_dir,
  )

  ch_sequence_flow_inputs = LONG_READ_DENOVO_ENSEMBLY.out.assembly_eligibility
    .map { _sample_id, report -> report }
    .mix(LONG_READ_DENOVO_ENSEMBLY.out.union_summaries.map { _sample_id, summary -> summary })
    .mix(PREPARE_BLAST_QUERIES.out.contig_filter_decisions.map { _sample_id, decision -> decision })
    .mix(PREPARE_BLAST_QUERIES.out.mapback_count_files.map { _sample_id, counts -> counts })
    .mix(PREPARE_BLAST_QUERIES.out.blast_query_summaries.map { _sample_id, summary -> summary })
    .mix(CLASSIFY_WITH_MEGABLAST.out.megablast_query_partition.map { _sample_id, _query_class, _accounted_ids, _blastn_candidates, summary -> summary })
    .mix(CLASSIFY_WITH_MEGABLAST.out.filter_decisions.map { _sample_id, _query_class, decision -> decision })
    .mix(CLASSIFY_WITH_BLASTN.out.filter_decisions.map { _sample_id, _query_class, decision -> decision })

  MULTIQC_BUNDLING(
    PREPROCESS_READS.out.raw_fastqc_packages,
    PREPROCESS_READS.out.raw_fastqc_zips,
    GATHER_READS.out.resolved_manifest,
    ch_nvd_version_file,
    channel.value(params.experimental == true),
    channel.value(target_enrichment_enabled),
    channel.value(depletion_enabled),
    channel.value(!params.skip_assembly),
    channel.value(!params.skip_blast),
    PREPROCESS_READS.out.target_enrichment_stats,
    PREPROCESS_READS.out.depletion_stats,
    PREPROCESS_READS.out.processed_read_profiles,
    PREPROCESS_READS.out.processed_read_quality_histograms,
    PREPARE_BLAST_QUERIES.out.filtered_contig_profiles,
    SHORT_READ_DENOVO_ASSEMBLY.out.assembly_profiles.mix(LONG_READ_DENOVO_ENSEMBLY.out.assembly_profiles),
    SHORT_READ_DENOVO_ASSEMBLY.out.eligibility_decisions.mix(LONG_READ_DENOVO_ENSEMBLY.out.eligibility_decisions),
    LONG_READ_DENOVO_ENSEMBLY.out.eligibility_summaries,
    LONG_READ_DENOVO_ENSEMBLY.out.union_summaries,
    PREPARE_BLAST_QUERIES.out.blast_query_summaries,
    CLASSIFY_WITH_MEGABLAST.out.megablast_query_partition.map { sample_id, query_class, _accounted_ids, _blastn_candidates, summary -> tuple(sample_id, query_class, summary) },
    channel.value(file("${projectDir}/assets/multiqc_config.yaml")),
  )

  REPORTING(
    CLASSIFY_WITH_BLASTN.out.merged_results,
    PREPROCESS_READS.out.read_counts,
    PREPARE_BLAST_QUERIES.out.contig_sequences,
    PREPARE_BLAST_QUERIES.out.query_lookups,
    PREPARE_BLAST_QUERIES.out.contig_read_counts,
    PREPARE_BLAST_QUERIES.out.filtered_bam,
    PREPARE_BLAST_QUERIES.out.no_contigs,
    PREPROCESS_READS.out.target_enrichment_stats,
    ch_taxonomy_dir,
    COMPUTE_RUN_CONTEXT.out.ready,
    ch_run_context,
    ch_sourmash_tax_reports,
    ch_risk_group_lookup,
    ch_sequence_flow_inputs,
    MULTIQC_BUNDLING.out.fastqc_zips,
    MULTIQC_BUNDLING.out.inputs,
    MULTIQC_BUNDLING.out.config,
    workflow.runName,
  )

  ch_scientific_completion = REPORTING.out.blast_results.count()

  if (params.experimental == true) {
    RAPID_SCREENING_EVAL(
      PREPROCESS_READS.out.read_counts,
      ch_sourmash_gather_csv,
      ch_sourmash_tax_reports,
      ch_sourmash_lineages,
      REPORTING.out.blast_results,
      REPORTING.out.crumbs_taxa,
      REPORTING.out.crumbs_queries,
      workflow.runName,
    )
  }

  emit:
  completion = ch_scientific_completion.map { n -> "NVD main workflow complete: ${n} samples processed" }
  labkey_log = REPORTING.out.labkey_log
}
