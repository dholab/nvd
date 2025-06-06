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
    publishDir params.outdir, mode: params.publish_dir_mode

    input:
    path csv_files
    val experiment_id
    val run_id

    output:
    path "blast_labkey_upload.log", emit: log

    script:
    """
    #!/usr/bin/env python3

    import os
    import csv
    import sys
    from datetime import datetime

    # LabKey configuration from nextflow params
    LABKEY_SERVER = "${params.labkey_server}"
    PROJECT_NAME = "${params.labkey_project_name}"
    API_KEY = "${params.labkey_api_key}"
    SCHEMA = "${params.labkey_schema}"
    TABLE_NAME = "metagenomic_hits_test_nvd2"

    log_entries = []
    log_entries.append(f"LabKey BLAST Upload Log - {datetime.now()}")
    log_entries.append(f"Experiment ID: ${experiment_id}")
    log_entries.append(f"Run ID: ${run_id}")
    log_entries.append(f"Server: {LABKEY_SERVER}")
    log_entries.append(f"Project: {PROJECT_NAME}")
    log_entries.append(f"Target Table: {TABLE_NAME}")
    log_entries.append("=" * 80)

    # Check if we have LabKey credentials
    upload_enabled = bool(LABKEY_SERVER and PROJECT_NAME and API_KEY)

    if upload_enabled:
        try:
            from labkey.api_wrapper import APIWrapper
            import more_itertools
            log_entries.append("LabKey API wrapper found - attempting real upload")

            # Initialize API wrapper
            api = APIWrapper(LABKEY_SERVER, PROJECT_NAME, api_key=API_KEY)

        except ImportError as e:
            log_entries.append(f"ERROR: LabKey API wrapper not available - {str(e)}")
            log_entries.append("Falling back to simulation mode")
            upload_enabled = False
    else:
        log_entries.append("No LabKey credentials provided - running in simulation mode")

    # Get BLAST CSV files
    csv_files = [f for f in os.listdir('.') if f.endswith('.csv') and 'blast' in f.lower()]
    log_entries.append(f"Found {len(csv_files)} BLAST CSV files: {csv_files}")

    total_records_processed = 0
    total_records_uploaded = 0

    for csv_file in sorted(csv_files):
        log_entries.append(f"\\nProcessing BLAST file: {csv_file}")

        # Count records
        record_count = 0
        if os.path.getsize(csv_file) > 0:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                records = list(reader)
                record_count = len(records)
                total_records_processed += record_count

                if record_count > 0:
                    log_entries.append(f"  Records: {record_count}")
                    log_entries.append(f"  Sample fields: {', '.join(list(records[0].keys())[:5])}")

                    if upload_enabled:
                        # Upload using LabKey API wrapper in chunks of 1000
                        batch_size = 1000
                        record_chunks = list(more_itertools.chunked(records, batch_size))
                        success_count = 0

                        for i, chunk in enumerate(record_chunks, 1):
                            try:
                                api.query.insert_rows(
                                    schema_name=SCHEMA,
                                    query_name=TABLE_NAME,
                                    rows=chunk
                                )
                                log_entries.append(f"    Batch {i}/{len(record_chunks)}: SUCCESS ({len(chunk)} records)")
                                success_count += len(chunk)

                            except Exception as e:
                                log_entries.append(f"    Batch {i}/{len(record_chunks)}: ERROR - {str(e)}")

                        total_records_uploaded += success_count
                        log_entries.append(f"  Upload complete: {success_count}/{record_count} records successful")
                    else:
                        # Simulation mode
                        num_batches = (record_count + 999) // 1000  # Round up division
                        log_entries.append(f"  Would upload in {num_batches} batch(es) (SIMULATION)")
                        total_records_uploaded += record_count  # For simulation counting
                else:
                    log_entries.append(f"  No records found in file")
        else:
            log_entries.append(f"  Empty file - no records to upload")

    log_entries.append("\\n" + "=" * 80)
    log_entries.append("BLAST DATA UPLOAD SUMMARY")
    log_entries.append(f"Files processed: {len(csv_files)}")
    log_entries.append(f"Total records processed: {total_records_processed}")
    if upload_enabled:
        log_entries.append(f"Total records uploaded: {total_records_uploaded}")
        log_entries.append("BLAST UPLOAD COMPLETE")
    else:
        log_entries.append(f"Total records that would be uploaded: {total_records_uploaded}")
        log_entries.append("BLAST SIMULATION COMPLETE - No actual upload performed")

    # Write log file
    with open('blast_labkey_upload.log', 'w') as f:
        f.write('\\n'.join(log_entries))
    """
}

process LABKEY_UPLOAD_FASTA {
    label 'process_low'
    publishDir params.outdir, mode: params.publish_dir_mode

    input:
    path csv_files
    val experiment_id
    val run_id

    output:
    path "fasta_labkey_upload.log", emit: log

    script:
    """
    #!/usr/bin/env python3

    import os
    import csv
    import sys
    from datetime import datetime

    # LabKey configuration from nextflow params
    LABKEY_SERVER = "${params.labkey_server}"
    PROJECT_NAME = "${params.labkey_project_name}"
    API_KEY = "${params.labkey_api_key}"
    SCHEMA = "${params.labkey_schema}"
    TABLE_NAME = "fasta_hits_test_nvd2"

    log_entries = []
    log_entries.append(f"LabKey FASTA Upload Log - {datetime.now()}")
    log_entries.append(f"Experiment ID: ${experiment_id}")
    log_entries.append(f"Run ID: ${run_id}")
    log_entries.append(f"Server: {LABKEY_SERVER}")
    log_entries.append(f"Project: {PROJECT_NAME}")
    log_entries.append(f"Target Table: {TABLE_NAME}")
    log_entries.append("=" * 80)

    # Check if we have LabKey credentials
    upload_enabled = bool(LABKEY_SERVER and PROJECT_NAME and API_KEY)

    if upload_enabled:
        try:
            from labkey.api_wrapper import APIWrapper
            import more_itertools
            log_entries.append("LabKey API wrapper found - attempting real upload")

            # Initialize API wrapper
            api = APIWrapper(LABKEY_SERVER, PROJECT_NAME, api_key=API_KEY)

        except ImportError as e:
            log_entries.append(f"ERROR: LabKey API wrapper not available - {str(e)}")
            log_entries.append("Falling back to simulation mode")
            upload_enabled = False
    else:
        log_entries.append("No LabKey credentials provided - running in simulation mode")

    # Get FASTA CSV files
    csv_files = [f for f in os.listdir('.') if f.endswith('.csv') and 'fasta' in f.lower()]
    log_entries.append(f"Found {len(csv_files)} FASTA CSV files: {csv_files}")

    total_records_processed = 0
    total_records_uploaded = 0

    for csv_file in sorted(csv_files):
        log_entries.append(f"\\nProcessing FASTA file: {csv_file}")

        # Count records
        record_count = 0
        if os.path.getsize(csv_file) > 0:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                records = list(reader)
                record_count = len(records)
                total_records_processed += record_count

                if record_count > 0:
                    log_entries.append(f"  Records: {record_count}")
                    log_entries.append(f"  Sample fields: {', '.join(list(records[0].keys())[:5])}")

                    if upload_enabled:
                        # Upload using LabKey API wrapper in chunks of 1000
                        batch_size = 1000
                        record_chunks = list(more_itertools.chunked(records, batch_size))
                        success_count = 0

                        for i, chunk in enumerate(record_chunks, 1):
                            try:
                                api.query.insert_rows(
                                    schema_name=SCHEMA,
                                    query_name=TABLE_NAME,
                                    rows=chunk
                                )
                                log_entries.append(f"    Batch {i}/{len(record_chunks)}: SUCCESS ({len(chunk)} records)")
                                success_count += len(chunk)

                            except Exception as e:
                                log_entries.append(f"    Batch {i}/{len(record_chunks)}: ERROR - {str(e)}")

                        total_records_uploaded += success_count
                        log_entries.append(f"  Upload complete: {success_count}/{record_count} records successful")
                    else:
                        # Simulation mode
                        num_batches = (record_count + 999) // 1000  # Round up division
                        log_entries.append(f"  Would upload in {num_batches} batch(es) (SIMULATION)")
                        total_records_uploaded += record_count  # For simulation counting
                else:
                    log_entries.append(f"  No records found in file")
        else:
            log_entries.append(f"  Empty file - no records to upload")

    log_entries.append("\\n" + "=" * 80)
    log_entries.append("FASTA DATA UPLOAD SUMMARY")
    log_entries.append(f"Files processed: {len(csv_files)}")
    log_entries.append(f"Total records processed: {total_records_processed}")
    if upload_enabled:
        log_entries.append(f"Total records uploaded: {total_records_uploaded}")
        log_entries.append("FASTA UPLOAD COMPLETE")
    else:
        log_entries.append(f"Total records that would be uploaded: {total_records_uploaded}")
        log_entries.append("FASTA SIMULATION COMPLETE - No actual upload performed")

    # Write log file
    with open('fasta_labkey_upload.log', 'w') as f:
        f.write('\\n'.join(log_entries))
    """
}
