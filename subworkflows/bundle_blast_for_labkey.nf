// meta is equivalent to sample_id
workflow BUNDLE_BLAST_FOR_LABKEY {
    take:
    blast_results        // channel: [ meta, csv ]
    contig_sequences     // channel: [ meta, fasta ]
    read_counts          // channel: [ meta, total_reads ]
    experiment_id        // value: experiment ID
    run_id               // value: workflow run ID
    contig_read_counts   // channel: [ meta, mapped_counts.tsv ]

    main:

    // DEBUG: Add comprehensive channel debugging
    // blast_results.view { "DEBUG BLAST_RESULTS: $it" }
    // read_counts.view { "DEBUG READ_COUNTS: $it" }
    // contig_sequences.view { "DEBUG CONTIG_SEQUENCES: $it" }

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
        run_id
    )

    // JOIN the channels to ensure same sample_id for WebDAV upload
    ch_webdav_upload = blast_results
        .join(contig_sequences)  // This ensures matching sample_ids
        .map { meta, blast_csv, fasta ->
            [meta, blast_csv, fasta]
        }

    // Upload BLAST CSV and fasta per sample to WebDAV
    WEBDAV_UPLOAD_BLAST(
        ch_webdav_upload
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
        run_id
    )

    // Upload BLAST results
    LABKEY_UPLOAD_BLAST(
        PREPARE_BLAST_LABKEY.out.csv.map { _meta, path -> path }.collect(),
        experiment_id,
        run_id
    )

    // Upload FASTA results
    LABKEY_UPLOAD_FASTA(
        PREPARE_FASTA_LABKEY.out.csv.map { _meta, path -> path }.collect(),
        experiment_id,
        run_id
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

    secret 'nvd2'

    input:
    tuple val(sample_id), path(blast_csv), path(fasta)

    output:
    val(sample_id)

    script:
    """
    gzip -c ${fasta} > ${fasta}.gz
    gzip -c ${blast_csv} > ${blast_csv}.gz

    # Upload tabular Gottcha2 report
    webdav_CLIent.py \
        --password \$nvd2 \
        --server ${params.labkey_webdav} \
        upload ${blast_csv}.gz ${params.experiment_id}/${sample_id}/nvd/${blast_csv}.gz

    # Upload fasta file
    webdav_CLIent.py \
        --password \$nvd2 \
        --server ${params.labkey_webdav} \
        upload ${fasta}.gz ${params.experiment_id}/${sample_id}/nvd/${fasta}.gz
    """
}

process PREPARE_BLAST_LABKEY {
    tag "$meta"
    label 'process_low'

    input:
    tuple val(meta), path(blast_csv), val(total_reads), path(contig_mapped_read_counts), val(output_name)
    val experiment_id
    val run_id

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
    label 'process_low'

    input:
    tuple val(meta), path(fasta), val(output_name)
    val experiment_id
    val run_id

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
    label 'process_low'

    secret 'nvd2'

    input:
    path csv_files
    val experiment_id
    val run_id

    output:
    path "blast_labkey_upload.log", emit: log

    script:
    """
    labkey_upload_blast_results.py \
    --experiment-id '${experiment_id}' \
    --run-id '${run_id}' \
    --labkey-server '${params.labkey_server}' \
    --labkey-project-name '${params.labkey_project_name}' \
    --labkey-api-key \$nvd2 \
    --labkey-schema '${params.labkey_schema}' \
    --table-name '${params.labkey_blast_meta_hits_list}'
    """
}

// Uploading BLAST FASTA
process LABKEY_UPLOAD_FASTA {
    label 'process_low'

    secret 'nvd2'

    input:
    path csv_files
    val experiment_id
    val run_id

    output:
    path "fasta_labkey_upload.log", emit: log

    script:
    """
    labkey_upload_blast_fasta.py \
    --experiment-id '${experiment_id}' \
    --run-id '${run_id}' \
    --labkey-server '${params.labkey_server}' \
    --labkey-project-name '${params.labkey_project_name}' \
    --labkey-api-key \$nvd2 \
    --labkey-schema '${params.labkey_schema}' \
    --table-name '${params.labkey_blast_fasta_list}'
    """
}
