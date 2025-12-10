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
include { RETRIEVE_GETTAX } from "../modules/utils"
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
    }

    ch_blast_db_files = Channel.fromPath(params.blast_db)
    ch_stat_index = Channel.fromPath(params.stat_index)
    ch_stat_dbss = Channel.fromPath(params.stat_dbss)
    ch_stat_annotation = Channel.fromPath(params.stat_annotation)
    ch_human_virus_taxlist = Channel.fromPath(params.human_virus_taxlist)
    ch_gettax = RETRIEVE_GETTAX()


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
        ch_stat_annotation,
        ch_gettax
    )

    CLASSIFY_WITH_MEGABLAST(
        EXTRACT_HUMAN_VIRUSES.out.contigs,
        ch_blast_db_files,
        EXTRACT_HUMAN_VIRUSES.out.sqlite
    )

    CLASSIFY_WITH_BLASTN(
        CLASSIFY_WITH_MEGABLAST.out.filtered_megablast,
        CLASSIFY_WITH_MEGABLAST.out.megablast_contigs,
        ch_blast_db_files,
        EXTRACT_HUMAN_VIRUSES.out.sqlite
    )

    if (params.labkey) {
        // First: Check experiment ID is fresh (not already uploaded)
        // Use .first() to get a value channel trigger without consuming the queue channel
        VALIDATE_LK_EXP_FRESH(
            CLASSIFY_WITH_BLASTN.out.merged_results.first()
        )

        // Bundle and upload - depends on validation passing
        BUNDLE_BLAST_FOR_LABKEY(
            CLASSIFY_WITH_BLASTN.out.merged_results,
            EXTRACT_HUMAN_VIRUSES.out.contigs,
            COUNT_READS.out.counts,
            params.experiment_id,
            workflow.runName,
            EXTRACT_HUMAN_VIRUSES.out.contig_read_counts,
            VALIDATE_LK_EXP_FRESH.out.validated
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

    ch_completion = CLASSIFY_WITH_BLASTN.out.merged_results.map { _results -> "STAT+BLAST workflow complete!" }

    emit:
    completion = ch_completion
    labkey_log = labkey_log_ch
}
