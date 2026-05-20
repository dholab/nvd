include {
    VALIDATE_LABKEY_BLAST_HITS_LIST ;
    VALIDATE_LABKEY_BLAST_FASTA_LIST ;
    VALIDATE_LABKEY_EXPERIMENT_FRESH ;
    REGISTER_LABKEY_EXPERIMENT ;
    PREPARE_BLAST_LABKEY ;
    PREPARE_FASTA_LABKEY ;
    CONCAT_ALL_SAMPLE_BLAST_RESULTS ;
    WEBDAV_UPLOAD_BLAST ;
    WEBDAV_UPLOAD_CONCATENATED ;
    LABKEY_UPLOAD_BLAST ;
    LABKEY_UPLOAD_FASTA
} from "../modules/labkey"

workflow LABKEY_BLAST_REPORTING {
    take:
    blast_results        // queue channel: [ sample_id, csv ] - one per sample with retained BLAST hits
    contig_sequences     // queue channel: [ sample_id, fasta ] - one per sample
    read_counts          // queue channel: [ sample_id, total_reads ] - one per sample
    experiment_id        // value channel: experiment ID
    run_id               // value channel: workflow run ID
    contig_read_counts   // queue channel: [ sample_id, mapped_counts.tsv ] - one per sample
    run_ready            // value channel: gate ensuring upstream preflight passed

    main:
    ch_labkey_blast_results = params.labkey
        ? blast_results
        : channel.empty()
    ch_labkey_contigs = params.labkey
        ? contig_sequences
        : channel.empty()
    ch_labkey_read_counts = params.labkey
        ? read_counts
        : channel.empty()
    ch_labkey_contig_read_counts = params.labkey
        ? contig_read_counts
        : channel.empty()

    ch_labkey_has_hits = ch_labkey_blast_results
        .first()
        .map { _first_result -> true }

    VALIDATE_LABKEY_BLAST_HITS_LIST(ch_labkey_has_hits)

    VALIDATE_LABKEY_BLAST_FASTA_LIST(ch_labkey_has_hits)

    ch_labkey_list_validation = VALIDATE_LABKEY_BLAST_HITS_LIST.out.validated
        .combine(VALIDATE_LABKEY_BLAST_FASTA_LIST.out.validated)
        .map { _hits, _fasta -> true }

    VALIDATE_LABKEY_EXPERIMENT_FRESH(ch_labkey_has_hits)

    ch_validation_gate = run_ready
        .combine(ch_labkey_list_validation)
        .combine(VALIDATE_LABKEY_EXPERIMENT_FRESH.out.validated)
        .map { _ready_and_list_validation, _fresh -> true }
        .first()

    ch_all_sample_data = ch_labkey_blast_results
        .join(ch_labkey_contigs, by: 0)
        .join(ch_labkey_read_counts, by: 0)
        .join(ch_labkey_contig_read_counts, by: 0)

    ch_split = ch_all_sample_data
        .multiMap { sample_id, blast_csv, fasta, total_reads, mapped_counts ->
            def blast_output = "${sample_id}_blast_labkey.csv"
            def fasta_output = "${sample_id}_fasta_labkey.csv"
            blast_labkey: [sample_id, blast_csv, total_reads, mapped_counts, blast_output]
            webdav_upload: [sample_id, blast_csv, fasta]
            fasta_labkey: [sample_id, fasta, fasta_output]
        }

    PREPARE_BLAST_LABKEY(
        ch_split.blast_labkey,
        experiment_id,
        run_id,
        ch_validation_gate,
    )

    WEBDAV_UPLOAD_BLAST(
        ch_split.webdav_upload,
        ch_validation_gate,
    )

    PREPARE_FASTA_LABKEY(
        ch_split.fasta_labkey,
        experiment_id,
        run_id,
        ch_validation_gate,
    )

    ch_prepared_blast_csvs = PREPARE_BLAST_LABKEY.out.csv
        .map { _sample_id, csv -> csv }
        .collect()
        .filter { files -> files.size() > 0 }

    CONCAT_ALL_SAMPLE_BLAST_RESULTS(
        ch_prepared_blast_csvs,
        experiment_id,
        ch_validation_gate,
    )

    WEBDAV_UPLOAD_CONCATENATED(
        CONCAT_ALL_SAMPLE_BLAST_RESULTS.out.concatenated_csv,
        ch_validation_gate,
    )

    LABKEY_UPLOAD_BLAST(
        PREPARE_BLAST_LABKEY.out.csv,
        experiment_id,
        run_id,
    )

    LABKEY_UPLOAD_FASTA(
        PREPARE_FASTA_LABKEY.out.csv,
        experiment_id,
        run_id,
    )

    ch_upload_complete = WEBDAV_UPLOAD_BLAST.out.done
        .mix(WEBDAV_UPLOAD_CONCATENATED.out.done)
        .mix(LABKEY_UPLOAD_BLAST.out.log)
        .mix(LABKEY_UPLOAD_FASTA.out.log)
        .collect()
        .filter { events -> events.size() > 0 }
        .map { _events -> true }

    REGISTER_LABKEY_EXPERIMENT(ch_upload_complete)

    ch_final_labkey_log = LABKEY_UPLOAD_BLAST.out.log
        .mix(LABKEY_UPLOAD_FASTA.out.log)
        .collectFile(
            name: 'final_labkey_upload.log',
            storeDir: params.results,
        )

    emit:
    upload_log = LABKEY_UPLOAD_BLAST.out.log.mix(LABKEY_UPLOAD_FASTA.out.log)
    final_labkey_log = ch_final_labkey_log
    registered = REGISTER_LABKEY_EXPERIMENT.out.registered
}
