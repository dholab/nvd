include {
    LABKEY_VALIDATE_BLAST_HITS_LIST ;
    LABKEY_VALIDATE_BLAST_FASTA_LIST ;
    LABKEY_VALIDATE_EXPERIMENT_FRESH ;
    LABKEY_REGISTER_EXPERIMENT ;
    LABKEY_PREPARE_BLAST ;
    LABKEY_PREPARE_FASTA ;
    LABKEY_CONCAT_ALL_SAMPLE_BLAST_RESULTS ;
    LABKEY_WEBDAV_UPLOAD_BLAST ;
    LABKEY_WEBDAV_UPLOAD_CONCATENATED ;
    LABKEY_UPLOAD_BLAST ;
    LABKEY_UPLOAD_FASTA
} from "../modules/labkey"

workflow LIMS_INTEGRATION {
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

    LABKEY_VALIDATE_BLAST_HITS_LIST(ch_labkey_has_hits)

    LABKEY_VALIDATE_BLAST_FASTA_LIST(ch_labkey_has_hits)

    ch_labkey_list_validation = LABKEY_VALIDATE_BLAST_HITS_LIST.out.validated
        .combine(LABKEY_VALIDATE_BLAST_FASTA_LIST.out.validated)
        .map { _hits, _fasta -> true }

    LABKEY_VALIDATE_EXPERIMENT_FRESH(ch_labkey_has_hits)

    ch_validation_gate = run_ready
        .combine(ch_labkey_list_validation)
        .combine(LABKEY_VALIDATE_EXPERIMENT_FRESH.out.validated)
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

    LABKEY_PREPARE_BLAST(
        ch_split.blast_labkey,
        experiment_id,
        run_id,
        ch_validation_gate,
    )

    LABKEY_WEBDAV_UPLOAD_BLAST(
        ch_split.webdav_upload,
        ch_validation_gate,
    )

    LABKEY_PREPARE_FASTA(
        ch_split.fasta_labkey,
        experiment_id,
        run_id,
        ch_validation_gate,
    )

    ch_prepared_blast_csvs = LABKEY_PREPARE_BLAST.out.csv
        .map { _sample_id, csv -> csv }
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
        run_id,
    )

    LABKEY_UPLOAD_FASTA(
        LABKEY_PREPARE_FASTA.out.csv,
        experiment_id,
        run_id,
    )

    ch_upload_complete = LABKEY_WEBDAV_UPLOAD_BLAST.out.done
        .mix(LABKEY_WEBDAV_UPLOAD_CONCATENATED.out.done)
        .mix(LABKEY_UPLOAD_BLAST.out.log)
        .mix(LABKEY_UPLOAD_FASTA.out.log)
        .collect()
        .filter { events -> events.size() > 0 }
        .map { _events -> true }

    LABKEY_REGISTER_EXPERIMENT(ch_upload_complete)

    ch_final_labkey_log = LABKEY_UPLOAD_BLAST.out.log
        .mix(LABKEY_UPLOAD_FASTA.out.log)
        .collectFile(
            name: 'final_labkey_upload.log',
            storeDir: params.results,
        )

    emit:
    upload_log = LABKEY_UPLOAD_BLAST.out.log.mix(LABKEY_UPLOAD_FASTA.out.log)
    final_labkey_log = ch_final_labkey_log
    registered = LABKEY_REGISTER_EXPERIMENT.out.registered
}
