/*
Perform rapid sequence similarity search using MEGABLAST.

This rule aligns the assembled contigs against a comprehensive nucleotide database
to identify similar known sequences. It's crucial for initial, broad_scale
identification of viral sequences.
*/
process MEGABLAST {
    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 6

    input:
    tuple val(sample_id), path(human_virus_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.megablast.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def blast_task = "megablast"
    """
    export BLASTDB=${blast_db}
    blastn -task ${blast_task} \
    -db ${blast_db}/${params.blast_db_prefix} \
    -query ${human_virus_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    > ${sample_id}.megablast.txt
    """
}

//Create file with megablast results annotated with megablast task and full taxonomic rank.
process ANNOTATE_MEGABLAST_RESULTS {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(sample_id), path(blast_txt), path(sqlite_file)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_megablast.txt"), emit: hits

    script:
    """
    annotate_blast_results.py \
    --sqlite_cache ${sqlite_file} \
    --input_file ${blast_txt} \
    --output_file ${sample_id}.annotated_megablast.txt \
    --sample_name ${sample_id} \
    --task 'megablast'
    """
}

/*
Remove any contigs that do not have at least one hit corresponding to viruses.

This rule handles cases where the input file is empty and filters out
blast query read groups where none of the entries are from superkingdom: Viruses
or where the only viral hits are phages. This will leave in any query seq that has things that are 
non virus as long as one of the qseqid reads has a viral hit.
*/
process FILTER_NON_VIRUS_MEGABLAST_NODES {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.mb_virus_only.txt")

    script:
    """
    filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.mb_virus_only.txt
    """
}

/*
Create a list of contigs that are not classified by megablast.

This rule generates two outputs:
1. A list of classified contigs
2. A pruned FASTA file containing only unclassified contigs

The pruned file is used as input for more sensitive blastn searching.
*/
process REMOVE_MEGABLAST_MAPPED_CONTIGS {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(sample_id), path(filtered_megablast_nodes), path(human_virus_contigs)

    output:
    tuple val(sample_id), path("${sample_id}.classified.txt"), path("${sample_id}.pruned.fasta")

    script:
    """
    remove_megablast_mapped_contigs.py \
    --megablast_results ${filtered_megablast_nodes} \
    --contigs_fasta ${human_virus_contigs} \
    --pruned_contigs ${sample_id}.pruned.fasta \
    --classified_contigs ${sample_id}.classified.txt
    """
}

/*
Perform sequence similarity search using BLASTN.

This rule aligns the assembled contigs against a comprehensive nucleotide database
to identify similar known sequences with BLASTN
*/
process BLASTN_CLASSIFY {
    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 8

    input:
    tuple val(sample_id), path(classified), path(pruned_contigs), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}.blastn.txt")

    script:
    def outfmt = "6 qseqid qlen sseqid stitle length pident evalue bitscore sscinames staxids"
    def max_target_seqs = 5
    def blast_task = "blastn"
    """
    export BLASTDb=${blast_db}/
    blastn -task ${blast_task} \
    -db ${blast_db}/${params.blast_db_prefix} \
    -query ${pruned_contigs} \
    -num_threads ${task.cpus} \
    -outfmt "${outfmt}" \
    -max_target_seqs ${max_target_seqs} \
    > ${sample_id}.blastn.txt
    """
}

// Create file with blastn results annotated with blastn task and full taxonomic rank.
process ANNOTATE_BLASTN_RESULTS {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(sample_id), path(blastn_txt), path(taxa_sqlite)

    output:
    tuple val(sample_id), path("${sample_id}.annotated_blastn.txt")

    script:
    """
    annotate_blast_results.py \
    --sample_name ${sample_id} \
    --sqlite_cache ${taxa_sqlite} \
    --input_file ${blastn_txt} \
    --output_file ${sample_id}.annotated_blastn.txt \
    --task 'blastn'
    """
}

/*
Remove any contigs that do not have at least one hit corresponding to viruses.

This rule handles cases where the input file is empty and filters out
groups where none of the entries are from superkingdom: Viruses
or where the only viral hits are phages.
*/
process FILTER_NON_VIRUS_BLASTN_NODES {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(sample_id), path(annotated_blast)

    output:
    tuple val(sample_id), path("${sample_id}.nt_virus_only.txt")

    script:
    """
    filter_non_virus_blast_nodes.py ${annotated_blast} ${sample_id}.nt_virus_only.txt
    """
}

// Create file with merged megablast and blastn results.
process MERGE_FILTERED_BLAST_RESULTS {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(sample_id), path(megablast_viruses)
    tuple val(_sample_id), path(blastn_viruses)

    output:
    tuple val(sample_id), path("${sample_id}_blast.merged.tsv")

    script:
    """
    if [ -s ${megablast_viruses} ] || [ -s ${blastn_viruses} ]; then
        cat ${megablast_viruses} ${blastn_viruses} > ${sample_id}_blast.merged.tsv
    else
        touch ${sample_id}_blast.merged.tsv
        echo "Both input files are empty. Created an empty output file."
    fi
    """

}

// Extract unclassified contigs that do not have BLAST hits.
// These could be highly divergent viruses, but are most often junk.
// Not used currently, removed from logic of nvd-2 by WG
process EXTRACT_UNCLASSIFIED_CONTIGS {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(sample_id), path(classified_list), path(contigs), path("${sample_id}_megablast.txt"), path("${sample_id}_blastn.txt")

    output:
    tuple val(sample_id), path("${sample_id}.unclassified.fasta")

    script:
    """
    extract_unclassified_contigs.py \
    --contigs ${contigs} \
    --megablast_results ${sample_id}_megablast.txt \
    --blastn_results ${sample_id}_blastn.txt \
    --unclassified ${sample_id}.unclassified.fasta
    """
}

// Step to filter down the extracted fasta to have only species and strain level
// representative sequences.
process FILTER_DOWN_TO_SPECIES_STRAIN {
    
    tag "${sample_id}"

    cpus 4
    
    input:
    tuple val(sample_id), path(extracted_fasta), path(full_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.fasta"), path(full_tsv)

    script:
    """
    filter_gottcha_fasta.py --input ${extracted_fasta} --output ${sample_id}.filtered.fasta
    """
}

// Extract sequences for whitelisted taxids and create per-strain FASTAs
process EXTRACT_WHITELISTED_FASTAS {

    tag "${sample_id} - ${taxid}"

    input:
    tuple val(sample_id), path(extracted_fasta), path(full_tsv), val(taxid), path(blast_db)

    output:
    tuple val(sample_id), val(taxid), path(full_tsv), path("${sample_id}_${taxid}.pre_blast_check.fasta"), path("${taxid}_taxids.txt"), path("*taxon_specific.fasta")

    script:
    """
    extract_whitelisted_reads.py \
    --taxid ${taxid} \
    --nodes ${blast_db}/nodes.dmp \
    --fasta ${extracted_fasta} \
    --output ${sample_id}_${taxid}.pre_blast_check.fasta \
    --taxid-list ${taxid}_taxids.txt \
    --sample-id ${sample_id}

    # Ensure at least one taxon_specific file exists for Nextflow
    if ! ls *taxon_specific.fasta 1> /dev/null 2>&1; then
        touch ${sample_id}-empty.taxon_specific.fasta
        echo "No taxon-specific files created, created empty placeholder"
    fi
    """
}

// It gets a bit confusing here but for every whitelisted TAXID it could be a species, a genus, a family.
// We then take that info and expand the nodes.dmp taxon file to find all strains (called no rank) that match.
// Then we extract each strain into its own fasta from the whitelisted TAXID so we can perform a BLAST run
// by strain. So we have expanded the process to run from per sample to per sample x white_list_taxis x child_taxid_fasta
process VERIFY_WITH_BLAST {

    tag "${sample_id} - ${parent_taxid} - ${strain_taxid}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), val(parent_taxid), path(full_tsv), path(descendants_list), val(strain_taxid), path(strain_fasta)

    output:
    tuple val(sample_id), val(strain_taxid), path("${sample_id}_${strain_taxid}.verified_hits.tsv"), path("${sample_id}_${strain_taxid}.classification.tsv")

    when:
    strain_fasta.size() > 0

    script:
    """
    export BLASTDB=${params.blast_db}

    # Only proceed if FASTA file is not empty
    if [ -s ${strain_fasta} ]; then
        blastn -task megablast -query ${strain_fasta} \
            -db ${params.blast_db}/${params.blast_db_prefix} \
            -word_size 64 \
            -num_threads ${task.cpus} \
            -outfmt '6 qseqid sseqid stitle pident qcovhsp evalue bitscore staxids sscinames' \
            -out ${sample_id}_${strain_taxid}.verified_hits.tsv

    # Classification script
    python3 - <<EOF
    import csv
    import sys

    blast_file = "${sample_id}_${strain_taxid}.verified_hits.tsv"
    output_file = "${sample_id}_${strain_taxid}.classification.tsv"
    expected_taxid = "${strain_taxid}"

    total_hits = 0
    correct_hits = 0

    try:
        with open(blast_file) as bf:
            reader = csv.reader(bf, delimiter='\\t')
            for row in reader:
                if len(row) >= 8:  # Ensure we have enough columns
                    total_hits += 1
                    staxids = row[7].split(';') if row[7] else []
                    if expected_taxid in staxids:
                        correct_hits += 1
    except FileNotFoundError:
        print(f"BLAST output file not found: {blast_file}")

    with open(output_file, "w") as outf:
        outf.write("strain_taxid\\ttotal_hits\\tmatching_hits\\tconfidence\\n")
        confidence = correct_hits / total_hits if total_hits > 0 else 0
        outf.write(f"{expected_taxid}\\t{total_hits}\\t{correct_hits}\\t{confidence:.3f}\\n")
    EOF
        else
            # Create empty output files for empty input
            touch ${sample_id}_${strain_taxid}.verified_hits.tsv
            echo -e "strain_taxid\\ttotal_hits\\tmatching_hits\\tconfidence\\n${strain_taxid}\\t0\\t0\\t0.000" > ${sample_id}_${strain_taxid}.classification.tsv
        fi
    """
}

// We have the values returned from the megablast
// We also have the number of representative fasta for each strain.
// We can check through the rep strain and check that 

process STACK_VERIFIED_TABLES {

    tag "${sample_id}"

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 2

    input:
    tuple val(sample_id), path("hits/*"), path("classification/*"), path(full_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.verified_hits.tsv"), path("${sample_id}.verified_classification.tsv"), path(full_tsv)

    script:
    """
    #!/usr/bin/env python3

    from pathlib import Path
    import pandas as pd
    import sys

    # Collect all TSV files
    hit_tsvs = list(Path("hits").glob("*.tsv"))
    class_tsvs = list(Path("classification").glob("*.tsv"))

    print(f"Found {len(hit_tsvs)} hit files and {len(class_tsvs)} classification files")

    # Process hits files
    if hit_tsvs:
        hit_dfs = []
        for tsv in hit_tsvs:
            try:
                df = pd.read_csv(tsv, sep='\\t', header=None)
                if not df.empty:
                    hit_dfs.append(df)
            except Exception as e:
                print(f"Error reading {tsv}: {e}")
        
        if hit_dfs:
            hit_df = pd.concat(hit_dfs, ignore_index=True).drop_duplicates()
            hit_df.to_csv("${sample_id}.verified_hits.tsv", sep='\\t', index=False, header=False)
        else:
            # Create empty file with proper structure
            open("${sample_id}.verified_hits.tsv", 'w').close()
    else:
        open("${sample_id}.verified_hits.tsv", 'w').close()

    # Process classification files
    if class_tsvs:
        class_dfs = []
        for tsv in class_tsvs:
            try:
                df = pd.read_csv(tsv, sep='\\t')
                if not df.empty:
                    class_dfs.append(df)
            except Exception as e:
                print(f"Error reading {tsv}: {e}")
        
        if class_dfs:
            class_df = pd.concat(class_dfs, ignore_index=True).drop_duplicates()
            class_df.to_csv("${sample_id}.verified_classification.tsv", sep='\\t', index=False)
        else:
            # Create empty file with header
            pd.DataFrame(columns=['strain_taxid', 'total_hits', 'matching_hits', 'confidence']).to_csv(
                "${sample_id}.verified_classification.tsv", sep='\\t', index=False)
    else:
        pd.DataFrame(columns=['strain_taxid', 'total_hits', 'matching_hits', 'confidence']).to_csv(
            "${sample_id}.verified_classification.tsv", sep='\\t', index=False)
    """
}
