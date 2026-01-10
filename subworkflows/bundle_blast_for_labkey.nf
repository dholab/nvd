// meta is equivalent to sample_id
workflow BUNDLE_BLAST_FOR_LABKEY {
    take:
    blast_results        // queue channel: [ meta, csv ] - one per sample
    contig_sequences     // queue channel: [ meta, fasta ] - one per sample
    read_counts          // queue channel: [ meta, total_reads ] - one per sample
    experiment_id        // value channel: experiment ID (single, broadcast)
    run_id               // value channel: workflow run ID (single, broadcast)
    contig_read_counts   // queue channel: [ meta, mapped_counts.tsv ] - one per sample
    validation_complete  // value channel: gate ensuring validation passed (single)
    run_context          // value channel: [ sample_set_id, state_dir ] - bundled tuple from CHECK_RUN_STATE, converted via .first()

    main:

    // MEGA-JOIN: Combine all input queue channels into a single stream.
    // Queue channels passed to subworkflows via take: are consumed by the first
    // operator, so we join everything upfront and then split with multiMap.
    ch_all_sample_data = blast_results
        .join(contig_sequences, by: 0)
        .join(read_counts, by: 0)
        .join(contig_read_counts, by: 0)
    // Result: [meta, blast_csv, fasta, total_reads, mapped_counts]

    // Use multiMap to split the joined channel into three separate streams.
    // This consumes ch_all_sample_data exactly once and creates independent
    // output channels that can each be consumed by their respective processes.
    ch_all_sample_data
        .multiMap { meta, blast_csv, fasta, total_reads, mapped_counts ->
            def blast_output = "${meta}_blast_labkey.csv"
            def fasta_output = "${meta}_fasta_labkey.csv"
            blast_labkey: [meta, blast_csv, total_reads, mapped_counts, blast_output]
            webdav_upload: [meta, blast_csv, fasta]
            fasta_labkey: [meta, fasta, fasta_output]
        }
        .set { ch_split }

    PREPARE_BLAST_LABKEY(
        ch_split.blast_labkey,
        experiment_id,
        run_id,
        validation_complete
    )

    WEBDAV_UPLOAD_BLAST(
        ch_split.webdav_upload,
        validation_complete
    )

    PREPARE_FASTA_LABKEY(
        ch_split.fasta_labkey,
        experiment_id,
        run_id,
        validation_complete
    )

    // Concatenate all sample BLAST results into single file
    CONCAT_ALL_SAMPLE_BLAST_RESULTS(
        PREPARE_BLAST_LABKEY.out.csv.map { _meta, path -> path }.collect(),
        experiment_id,
        validation_complete
    )

    // Upload concatenated BLAST results to experiment root via WebDAV
    WEBDAV_UPLOAD_CONCATENATED(
        CONCAT_ALL_SAMPLE_BLAST_RESULTS.out.concatenated_csv,
        validation_complete
    )

    // Upload BLAST results to LabKey
    LABKEY_UPLOAD_BLAST(
        PREPARE_BLAST_LABKEY.out.csv.map { _meta, path -> path }.collect(),
        experiment_id,
        run_id,
        run_context
    )

    // Upload FASTA results to LabKey
    LABKEY_UPLOAD_FASTA(
        PREPARE_FASTA_LABKEY.out.csv.map { _meta, path -> path }.collect(),
        experiment_id,
        run_id,
        run_context
    )

    emit:
    blast_csv = PREPARE_BLAST_LABKEY.out.csv
    fasta_csv = PREPARE_FASTA_LABKEY.out.csv
    blast_concatenated_csv = CONCAT_ALL_SAMPLE_BLAST_RESULTS.out.concatenated_csv
    blast_concatenated_uploaded = WEBDAV_UPLOAD_CONCATENATED.out.done
    blast_upload_log = LABKEY_UPLOAD_BLAST.out.log
    fasta_upload_log = LABKEY_UPLOAD_FASTA.out.log
    upload_log = LABKEY_UPLOAD_BLAST.out.log.mix(LABKEY_UPLOAD_FASTA.out.log)
}

process WEBDAV_UPLOAD_BLAST {
    tag "${sample_id}"

    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(blast_csv), path(fasta)
    val validation_complete

    output:
    val(sample_id)

    script:
    """
    gzip -c ${fasta} > ${fasta}.gz
    gzip -c ${blast_csv} > ${blast_csv}.gz

    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${blast_csv}.gz ${params.experiment_id}/${sample_id}/nvd/${blast_csv}.gz

    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${fasta}.gz ${params.experiment_id}/${sample_id}/nvd/${fasta}.gz
    """
}

process WEBDAV_UPLOAD_CONCATENATED {
    label 'low'

    secret 'LABKEY_API_KEY'

    input:
    path concatenated_csv
    val _validation_complete  // Gate: ensures validation passed

    output:
    val true, emit: done

    script:
    """
    gzip -c ${concatenated_csv} > ${concatenated_csv}.gz

    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${concatenated_csv}.gz ${params.experiment_id}/${concatenated_csv}.gz
    """
}

process PREPARE_BLAST_LABKEY {
    tag "$meta"
    label 'low'

    input:
    tuple val(meta), path(blast_csv), val(total_reads), path(contig_mapped_read_counts), val(output_name)
    val experiment_id
    val run_id
    val validation_complete  // Ensures validation passed before preparing

    output:
    tuple val(meta), path(output_name), emit: csv

    script:
    """
        prepare_blast_labkey.py \\
        --blast-csv ${blast_csv} \\
        --contig-counts ${contig_mapped_read_counts} \\
        --output ${output_name} \\
        --meta '${meta}' \\
        --experiment-id ${experiment_id} \\
        --run-id '${run_id}' \\
        --total-reads ${total_reads} \\
        --blast-db-version '${params.blast_db_version}' \\
        --stat-db-version '${params.stat_db_version}'
    """
}

process CONCAT_ALL_SAMPLE_BLAST_RESULTS {
    label 'low'

    input:
    path "blast_results/*.csv"
    val experiment_id
    val _validation_complete  // Gate: ensures validation passed

    output:
    path "${experiment_id}_blast_concatenated.csv", emit: concatenated_csv

    script:
    """
    #!/usr/bin/env python3
    import polars as pl
    from pathlib import Path

    files = sorted(Path("blast_results").glob("*.csv"))

    if not files:
        # No samples to concatenate - create empty file
        Path("${experiment_id}_blast_concatenated.csv").touch()
    else:
        pl.concat([pl.scan_csv(f) for f in files]).sink_csv("${experiment_id}_blast_concatenated.csv")
    """
}

process PREPARE_FASTA_LABKEY {
    tag "$meta"
    label 'low'

    input:
    tuple val(meta), path(fasta), val(output_name)
    val experiment_id
    val run_id
    val validation_complete  // Ensures validation passed before preparing

    output:
    tuple val(meta), path("${output_name}"), emit: csv

    script:
    """
    #!/usr/bin/env python3

    import csv
    from Bio import SeqIO

    # Prepare FASTA data for LabKey
    fasta_data = []

    for record in SeqIO.parse('${fasta}', 'fasta'):
        labkey_row = {
            'experiment': ${experiment_id},
            'sample_id': '${meta}',
            'contig_id': record.id,
            'contig_sequence': str(record.seq),
            'notes': '',
            'snakemake_run_id': '${run_id}'
        }
        fasta_data.append(labkey_row)

    # Write formatted data
    if fasta_data:
        with open('${output_name}', 'w') as f:
            fieldnames = ['experiment', 'sample_id', 'contig_id',
                         'contig_sequence', 'notes', 'snakemake_run_id']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(fasta_data)
    else:
        # Create empty file if no data
        open('${output_name}', 'w').close()
    """
}

process LABKEY_UPLOAD_BLAST {
    label 'low'

    secret 'LABKEY_API_KEY'

    input:
    path csv_files        // collected list of all sample CSVs (value channel from .collect())
    val experiment_id
    val run_id
    tuple val(sample_set_id), val(state_dir)  // run_context: bundled from CHECK_RUN_STATE

    output:
    path "blast_labkey_upload.log", emit: log
    path csv_files        // preserve input files in work dir (no emit label - intentional)

    script:
    def sample_set_arg = sample_set_id ? "--sample-set-id '${sample_set_id}'" : ""
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    """
    labkey_upload_blast_results.py \
    --experiment-id '${experiment_id}' \
    --run-id '${run_id}' \
    ${sample_set_arg} \
    ${state_dir_arg} \
    --labkey-server '${params.labkey_server}' \
    --labkey-project-name '${params.labkey_project_name}' \
    --labkey-api-key \$LABKEY_API_KEY \
    --labkey-schema '${params.labkey_schema}' \
    --table-name '${params.labkey_blast_meta_hits_list}'
    """
}

// Uploading BLAST FASTA
process LABKEY_UPLOAD_FASTA {
    label 'low'

    secret 'LABKEY_API_KEY'

    input:
    path csv_files        // collected list of all sample CSVs (value channel from .collect())
    val experiment_id
    val run_id
    tuple val(sample_set_id), val(state_dir)  // run_context: bundled from CHECK_RUN_STATE

    output:
    path "fasta_labkey_upload.log", emit: log

    script:
    def sample_set_arg = sample_set_id ? "--sample-set-id '${sample_set_id}'" : ""
    def state_dir_arg = state_dir ? "--state-dir '${state_dir}'" : ""
    """
    labkey_upload_blast_fasta.py \
    --experiment-id '${experiment_id}' \
    --run-id '${run_id}' \
    ${sample_set_arg} \
    ${state_dir_arg} \
    --labkey-server '${params.labkey_server}' \
    --labkey-project-name '${params.labkey_project_name}' \
    --labkey-api-key \$LABKEY_API_KEY \
    --labkey-schema '${params.labkey_schema}' \
    --table-name '${params.labkey_blast_fasta_list}'
    """
}
