include {
    LABKEY_VALIDATE_BLAST_HITS_LIST ;
    LABKEY_VALIDATE_BLAST_FASTA_LIST ;
    LABKEY_PREPARE_BLAST ;
    LABKEY_PREPARE_FASTA ;
    LABKEY_CONCAT_ALL_SAMPLE_BLAST_RESULTS ;
    LABKEY_WEBDAV_UPLOAD_BLAST ;
    LABKEY_WEBDAV_UPLOAD_QUERY_FASTA ;
    LABKEY_WEBDAV_UPLOAD_CONCATENATED ;
    LABKEY_UPLOAD_BLAST ;
    LABKEY_UPLOAD_FASTA
} from "../modules/labkey"

workflow LIMS_INTEGRATION {
    take:
    blast_results         // queue channel: [ sample_id, query_class, batch_final_tsv ] — enriched with mapped_reads, total_reads, blast_db_version, nextflow_run_id; per (sample_id, query_class) batch
    sample_blast_results  // queue channel: [ sample_id, blast_tsv ] - per-sample stack of every query_class batch; one per sample
    query_fastas          // queue channel: [ sample_id, query_class, batch_fasta ] - per-read-type query FASTAs actually BLASTed; per (sample_id, query_class) batch
    experiment_id         // value channel: experiment ID (the one LabKey-specific field)
    run_id                // value channel: workflow run ID (needed for FASTA prep and uploads)
    run_ready             // value channel: gate ensuring upstream preflight passed

    main:
    ch_labkey_blast_results = params.labkey
        ? blast_results
        : channel.empty()
    ch_labkey_sample_blast_results = params.labkey
        ? sample_blast_results
        : channel.empty()
    ch_labkey_query_fastas = params.labkey
        ? query_fastas
        : channel.empty()

    ch_labkey_has_hits = ch_labkey_blast_results
        .first()
        .map { _first_result -> true }

    LABKEY_VALIDATE_BLAST_HITS_LIST(ch_labkey_has_hits)

    LABKEY_VALIDATE_BLAST_FASTA_LIST(ch_labkey_has_hits)

    ch_labkey_list_validation = LABKEY_VALIDATE_BLAST_HITS_LIST.out.validated
        .combine(LABKEY_VALIDATE_BLAST_FASTA_LIST.out.validated)
        .map { _hits, _fasta -> true }

    ch_validation_gate = run_ready
        .combine(ch_labkey_list_validation)
        .map { _ready, _list_valid -> true }
        .first()

    // BLAST row insertion does not require a contig FASTA. This keeps valid
    // read-only samples eligible for the reduced-schema LIMS table.
    // Runs per (sample_id, query_class) batch so each read type can be
    // uploaded eagerly downstream.
    ch_blast_labkey = ch_labkey_blast_results.map { sample_id, query_class, blast_tsv ->
        tuple(sample_id, query_class, blast_tsv)
    }

    // This upload now carries the stacked BLAST TSV alone. It used to send the
    // contig FASTA alongside it, which was redundant: the per-query-class
    // FASTAs uploaded below are built by partitioning those same contigs, with
    // no filtering, so for a single-platform Illumina sample the two files were
    // byte-identical and LIMS received the same contigs twice.
    //
    // Dropping the FASTA also drops the inner join this used to need, which is
    // a second fix: that join silently excluded read-only samples, since a
    // sample with no contigs had nothing to join against and so never had its
    // BLAST TSV uploaded at all. Every sample with results now gets one.
    ch_webdav_upload = ch_labkey_sample_blast_results

    // FASTA rows now follow the BLAST hits pattern: one batch per
    // (sample_id, query_class), covering every class that was actually
    // queried rather than contigs alone. A queue channel cannot be consumed
    // twice, so the WebDAV artifact upload and the LabKey row insert each take
    // their own branch of the same stream.
    ch_query_fasta_split = ch_labkey_query_fastas
        .multiMap { sample_id, query_class, fasta ->
            webdav: tuple(sample_id, query_class, fasta)
            labkey_rows: tuple(sample_id, query_class, fasta)
        }

    LABKEY_PREPARE_BLAST(
        ch_blast_labkey,
        experiment_id,
        ch_validation_gate,
    )

    LABKEY_WEBDAV_UPLOAD_BLAST(
        ch_webdav_upload,
        ch_validation_gate,
    )

    // Per-read-type query FASTAs (contig, merged, single) — the sequences
    // that were actually BLASTed — published as file artifacts alongside the
    // LabKey rows built from the same batches.
    LABKEY_WEBDAV_UPLOAD_QUERY_FASTA(
        ch_query_fasta_split.webdav,
        ch_validation_gate,
    )

    LABKEY_PREPARE_FASTA(
        ch_query_fasta_split.labkey_rows,
        experiment_id,
        run_id,
        ch_validation_gate,
    )

    ch_prepared_blast_csvs = LABKEY_PREPARE_BLAST.out.csv
        .map { _sample_id, _query_class, csv -> csv }
        .collect()
        .filter { files -> files.size() > 0 }

    LABKEY_CONCAT_ALL_SAMPLE_BLAST_RESULTS(
        ch_prepared_blast_csvs,
        experiment_id,
        ch_validation_gate,
    )

    LABKEY_WEBDAV_UPLOAD_CONCATENATED(
        LABKEY_CONCAT_ALL_SAMPLE_BLAST_RESULTS.out.concatenated_csv,
        ch_validation_gate,
    )

    LABKEY_UPLOAD_BLAST(
        LABKEY_PREPARE_BLAST.out.csv,
        experiment_id,
    )

    LABKEY_UPLOAD_FASTA(
        LABKEY_PREPARE_FASTA.out.csv,
        experiment_id,
    )

    // Upload completion replaces experiment registration as the Slack gate.
    // Collects completion of ALL per-batch uploads: the WebDAV raw-file
    // uploads as well as the row-level BLAST/FASTA inserts.
    ch_uploads_done = LABKEY_WEBDAV_UPLOAD_BLAST.out.done
        .mix(LABKEY_WEBDAV_UPLOAD_QUERY_FASTA.out.published)
        .mix(LABKEY_WEBDAV_UPLOAD_CONCATENATED.out.done)
        .mix(LABKEY_UPLOAD_BLAST.out.log)
        .mix(LABKEY_UPLOAD_FASTA.out.log)
        .collect()
        .map { _events -> true }

    ch_final_labkey_log = params.labkey
        ? LABKEY_UPLOAD_BLAST.out.log
            .mix(LABKEY_UPLOAD_FASTA.out.log)
            .collectFile(
                name: 'final_labkey_upload.log',
                storeDir: params.labkey_uploads + '/upload_logs',
            )
        : channel.empty()

    emit:
    upload_log = LABKEY_UPLOAD_BLAST.out.log.mix(LABKEY_UPLOAD_FASTA.out.log)
    final_labkey_log = ch_final_labkey_log
    uploads_done = ch_uploads_done
}
