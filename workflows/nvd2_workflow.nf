include { PREPROCESS_CONTIGS } from "../subworkflows/preprocess_contigs"
include { EXTRACT_HUMAN_VIRUSES } from "../subworkflows/extract_human_viruses"
include { CLASSIFY_WITH_MEGABLAST } from "../subworkflows/classify_with_megablast"
include { CLASSIFY_WITH_BLASTN } from "../subworkflows/classify_with_blastn"
// include { } from "../subworkflows/bundle_for_labkey"


workflow NVD2_WORKFLOW  {
    take:
    ch_sample_fastqs // Queue channel of sample IDs, platforms, and (interleaved) FASTQ files: tuple val(sample_id), val(platform), path(fastq)

    main:
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

    ch_blast_db_files = Channel.fromPath(params.blast_db)
    ch_stat_index = Channel.fromPath(params.stat_index)
    ch_stat_dbss = Channel.fromPath(params.stat_dbss)
    ch_stat_annotation = Channel.fromPath(params.stat_annotation)
    ch_human_virus_taxlist = Channel.fromPath(params.human_virus_taxlist)

    PREPROCESS_CONTIGS(
        ch_sample_fastqs,
        ch_stat_dbss,
        ch_stat_annotation,
        ch_human_virus_taxlist        
    )

    EXTRACT_HUMAN_VIRUSES(
        PREPROCESS_CONTIGS.out,
        ch_stat_index,
        ch_stat_dbss,
        ch_stat_annotation
    )

    CLASSIFY_WITH_MEGABLAST(
        EXTRACT_HUMAN_VIRUSES.out.contigs,
        ch_blast_db_files,
        EXTRACT_HUMAN_VIRUSES.out.sqlite
    )

    CLASSIFY_WITH_BLASTN(
        CLASSIFY_WITH_MEGABLAST.out.contigs,
        CLASSIFY_WITH_MEGABLAST.out.megablast,
        ch_blast_db_files,
        EXTRACT_HUMAN_VIRUSES.out.sqlite
    )

    ch_completion = CLASSIFY_WITH_BLASTN.out.merged_results.map { _results -> "NVD complete!" }

    emit:
    completion = ch_completion

}
