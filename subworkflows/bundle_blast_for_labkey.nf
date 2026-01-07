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

    // DSL2 automatically forks channels when used by multiple consumers.
    // No explicit tap() needed - just use the channels directly.

    // Prepare BLAST results for metagenomic_hits table
    ch_blast_labkey = blast_results
        .join(read_counts)
        .join(contig_read_counts)
        .map { meta, csv, total_reads, contig_mapped_reads  ->
            def output_file = "${meta}_blast_labkey.csv"
            [meta, csv, total_reads, contig_mapped_reads, output_file]
        }

    PREPARE_BLAST_LABKEY(
        ch_blast_labkey,
        experiment_id,
        run_id,
        validation_complete
    )

    // JOIN the channels to ensure same sample_id for WebDAV upload
    ch_webdav_upload = blast_results
        .join(contig_sequences)  // This ensures matching sample_ids
        .map { meta, blast_csv, fasta ->
            [meta, blast_csv, fasta]
        }

    // Upload BLAST CSV and fasta per sample to WebDAV
    WEBDAV_UPLOAD_BLAST(
        ch_webdav_upload,
        validation_complete
    )

    // Prepare FASTA sequences for fasta_hits table
    ch_fasta_labkey = contig_sequences
        .map { meta, fasta ->
            def output_file = "${meta}_fasta_labkey.csv"
            [meta, fasta, output_file]
        }

    PREPARE_FASTA_LABKEY(
        ch_fasta_labkey,
        experiment_id,
        run_id,
        validation_complete
    )

    // Upload BLAST results
    // run_context tuple [sample_set_id, state_dir] is passed intact and destructured in the process
    LABKEY_UPLOAD_BLAST(
        PREPARE_BLAST_LABKEY.out.csv.map { _meta, path -> path }.collect(),
        experiment_id,
        run_id,
        run_context
    )

    // Upload FASTA results
    // run_context tuple [sample_set_id, state_dir] is passed intact and destructured in the process
    LABKEY_UPLOAD_FASTA(
        PREPARE_FASTA_LABKEY.out.csv.map { _meta, path -> path }.collect(),
        experiment_id,
        run_id,
        run_context
    )

    emit:
    blast_csv = PREPARE_BLAST_LABKEY.out.csv
    fasta_csv = PREPARE_FASTA_LABKEY.out.csv
    blast_upload_log = LABKEY_UPLOAD_BLAST.out.log
    fasta_upload_log = LABKEY_UPLOAD_FASTA.out.log
    upload_log = LABKEY_UPLOAD_BLAST.out.log.mix(LABKEY_UPLOAD_FASTA.out.log)
}

process WEBDAV_UPLOAD_BLAST {

    tag "${sample_id}"

    secret 'LABKEY_API_KEY'

    input:
    tuple val(sample_id), path(blast_csv), path(fasta)
    val validation_complete  // Ensures validation passed before uploading

    output:
    val(sample_id)

    script:
    """
    gzip -c ${fasta} > ${fasta}.gz
    gzip -c ${blast_csv} > ${blast_csv}.gz

    # Upload tabular Gottcha2 report
    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${blast_csv}.gz ${params.experiment_id}/${sample_id}/nvd/${blast_csv}.gz

    # Upload fasta file
    webdav_CLIent.py \
        --password \$LABKEY_API_KEY \
        --server ${params.labkey_webdav} \
        upload ${fasta}.gz ${params.experiment_id}/${sample_id}/nvd/${fasta}.gz
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
