process COLLECT_CRUMBS_TAXIDS {

    label "low"

    input:
    tuple val(profile_id), path("blast_results/*")

    output:
    tuple val(profile_id), path("crumbs.taxids.txt"), emit: taxids

    script:
    '''#!/usr/bin/env python3
import csv
from pathlib import Path

output = Path("crumbs.taxids.txt")
newline = chr(10)

seen = set()
taxids = []
for blast_tsv in sorted(Path("blast_results").glob("*.tsv")):
    with blast_tsv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "adjusted_taxid" not in reader.fieldnames:
            raise SystemExit(f"{blast_tsv}: missing adjusted_taxid column")

        for row in reader:
            taxid = (row.get("adjusted_taxid") or "").strip()
            if not taxid or taxid in seen:
                continue
            seen.add(taxid)
            taxids.append(taxid)

output.write_text("".join(f"{taxid}{newline}" for taxid in taxids), encoding="utf-8")
'''
}

process PREPARE_NCBI_PROFILE_TAXONOMY {

    label "low"

    input:
    tuple val(profile_id), path(taxids)
    val taxonomy_dir

    output:
    tuple val(profile_id), path("crumbs.profile_taxonomy.tsv"), emit: profile_taxonomy

    script:
    def taxonomy_mode_arg = params.taxonomy_mode ? "--taxonomy-mode '${params.taxonomy_mode}'" : ""
    def taxonomy_max_age_arg = params.taxonomy_max_age_days ? "--taxonomy-max-age-days ${params.taxonomy_max_age_days}" : ""
    """
    prepare_ncbi_profile_taxonomy.py \
        --taxids ${taxids} \
        --taxonomy-dir '${taxonomy_dir}' \
        --output crumbs.profile_taxonomy.tsv \
        ${taxonomy_mode_arg} \
        ${taxonomy_max_age_arg}
    """
}

process ESTIMATE_CRUMBS_PROFILE {

    tag "${sample_id}"
    label "low"

    input:
    tuple val(sample_id), path(blast_tsv), path(coverage_tsv), path(profile_taxonomy_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.crumbs.contigs.tsv"), emit: contigs
    tuple val(sample_id), path("${sample_id}.crumbs.taxa.tsv"), emit: taxa
    tuple val(sample_id), path("${sample_id}.crumbs.bioboxes.profile.tsv"), emit: bioboxes_profile
    tuple val(sample_id), path("${sample_id}.crumbs.qc.json"), emit: qc

    script:
    """
    estimate_crumbs_profile.py \
        --sample-id ${sample_id} \
        --blast-tsv ${blast_tsv} \
        --coverage-tsv ${coverage_tsv} \
        --profile-taxonomy-tsv ${profile_taxonomy_tsv} \
        --output-dir .
    """
}

process EXPORT_CRUMBS_TAXONOMIC_REPORTS {

    tag "${sample_id}"
    label "low"

    input:
    tuple val(sample_id), path(taxa_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.crumbs.krona.tsv"), emit: krona
    tuple val(sample_id), path("${sample_id}.crumbs.kreport.txt"), emit: kreport

    script:
    """
    export_crumbs_taxonomic_reports.py \
        --sample-id ${sample_id} \
        --taxa-tsv ${taxa_tsv} \
        --output-dir .
    """
}

process RENDER_CRUMBS_TAXBURST {

    tag "${sample_id}"
    label "low"

    input:
    tuple val(sample_id), path(taxa_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.crumbs.taxburst.html"), path("${sample_id}.crumbs.taxburst.json"), emit: reports

    script:
    """
    render_multisample_taxburst.py \
        --input-format crumbs \
        --summary "${sample_id}=${taxa_tsv}" \
        --output "${sample_id}.crumbs.taxburst.html" \
        --output-json "${sample_id}.crumbs.taxburst.json"
    """
}
