

workflow BUNDLE_FOR_LABKEY {
    take:
    blast_results     // channel: [ meta, csv ]
    contig_sequences  // channel: [ meta, fasta ]
    read_counts       // channel: [ meta, total_reads ]
    experiment_id   // value: experiment ID
    run_id           // value: workflow run ID

    main:

    // DEBUG: Add comprehensive channel debugging
    // blast_results.view { "DEBUG BLAST_RESULTS: $it" }
    // read_counts.view { "DEBUG READ_COUNTS: $it" }
    // contig_sequences.view { "DEBUG CONTIG_SEQUENCES: $it" }

    // Prepare BLAST results for metagenomic_hits table
    ch_blast_labkey = blast_results
        .join(read_counts)
        .map { meta, csv, total_reads ->
            def output_file = "${meta}_blast_labkey.csv"
            [meta, csv, total_reads, output_file]
        }

    PREPARE_BLAST_LABKEY(
        ch_blast_labkey,
        experiment_id,
        run_id
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
        PREPARE_BLAST_LABKEY.out.csv.map { meta, path -> path }.collect(),
        experiment_id,
        run_id
    )

    // Upload FASTA results
    LABKEY_UPLOAD_FASTA(
        PREPARE_FASTA_LABKEY.out.csv.map { meta, path -> path }.collect(),
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

process PREPARE_BLAST_LABKEY {
    tag "$meta"
    label 'process_low'

    input:
    tuple val(meta), path(blast_csv), val(total_reads), val(output_name)
    val experiment_id
    val run_id

    output:
    tuple val(meta), path("${output_name}"), emit: csv

    script:
    """
    #!/usr/bin/env python3

    import csv
    import os

    # Read BLAST results
    blast_data = []
    with open('${blast_csv}', 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('qseqid', '').lower() == 'qseqid':
                continue  # skip accidental header row
            # Format data for LabKey metagenomic_hits table
            labkey_row = {
                'experiment': ${experiment_id},
                'blast_task': row.get('task', 'megablast'),
                'sample_id': row.get('sample', '${meta}'),
                'qseqid': row.get('qseqid', ''),
                'qlen': row.get('qlen', ''),
                'sseqid': row.get('sseqid', ''),
                'stitle': row.get('stitle', ''),
                'tax_rank': row.get('rank', ''),
                'length': row.get('length', ''),
                'pident': row.get('pident', ''),
                'evalue': row.get('evalue', ''),
                'bitscore': row.get('bitscore', ''),
                'sscinames': row.get('sscinames', ''),
                'blast_db_version': 'unknown',
                'snakemake_run_id': '${run_id}',
                'mapped_reads': '0',
                'total_reads': '${total_reads}',
                'stat_db_version': 'unknown'
            }
            blast_data.append(labkey_row)

    # Write formatted data
    if blast_data:
        with open('${output_name}', 'w') as f:
            writer = csv.DictWriter(f, fieldnames=blast_data[0].keys())
            writer.writeheader()
            writer.writerows(blast_data)
    else:
        # Create empty file if no data
        open('${output_name}', 'w').close()
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
