include { BUILD_MULTIQC_INPUTS; PACKAGE_NVD_DEACON_REPORT as PACKAGE_TARGET_ENRICHMENT_REPORT; PACKAGE_NVD_DEACON_REPORT as PACKAGE_DEPLETION_REPORT; PACKAGE_NVD_PROFILE_REPORT as PACKAGE_FASTX_REPORT; PACKAGE_NVD_PROFILE_WITH_QUALITY_REPORT as PACKAGE_FASTX_QUALITY_REPORT; PACKAGE_NVD_PROFILE_REPORT as PACKAGE_ASSEMBLY_PROFILE_REPORT; PACKAGE_NVD_ELIGIBILITY_REPORT as PACKAGE_ASSEMBLY_ELIGIBILITY_REPORT; PACKAGE_NVD_LONG_READ_ELIGIBILITY_REPORT as PACKAGE_ASSEMBLY_LONG_READ_ELIGIBILITY_REPORT; PACKAGE_NVD_UNION_REPORT as PACKAGE_ASSEMBLY_UNION_REPORT; PACKAGE_NVD_PREPARED_QUERY_BATCHES_REPORT; PACKAGE_NVD_MEGABLAST_QUERY_PARTITION_REPORT; PACKAGE_NVD_TAXON_BIG_TABLE_REPORT } from "../modules/multiqc"

workflow MULTIQC_BUNDLING {
    take:
    ch_raw_fastqc_packages
    ch_raw_fastqc_zips
    ch_resolved_reads
    ch_nvd_version
    ch_experimental_enabled
    ch_target_enrichment_enabled
    ch_depletion_enabled
    ch_assembly_enabled
    ch_blast_enabled
    ch_target_enrichment_stats
    ch_depletion_stats
    ch_processed_read_profiles
    ch_processed_read_quality_histograms
    ch_filtered_contig_profiles
    ch_assembly_profiles
    ch_assembly_eligibility_decisions
    ch_long_read_eligibility_summaries
    ch_long_read_union_summaries
    ch_prepared_query_batch_summaries
    ch_megablast_query_partition_summaries
    ch_taxon_big_tables
    ch_multiqc_config

    main:
    ch_fastqc_packages = ch_raw_fastqc_packages.toList()
    ch_fastqc_zips = ch_raw_fastqc_zips.toList()

    PACKAGE_TARGET_ENRICHMENT_REPORT(
        ch_target_enrichment_stats.map { sample_id, deacon_json ->
            tuple('target_enrichment', sample_id, '', deacon_json)
        }
    )

    PACKAGE_DEPLETION_REPORT(
        ch_depletion_stats.map { meta, deacon_json ->
            tuple('depletion', meta.id, meta.query_class ?: '', deacon_json)
        }
    )

    ch_processed_read_quality_by_key = ch_processed_read_quality_histograms.map { meta, quality_histogram -> tuple(meta.profile_key, quality_histogram) }
    ch_processed_read_profiles_with_quality = ch_processed_read_profiles
        .map { meta, profile_json, length_histogram, _sequence_count -> tuple(meta.profile_key, meta, profile_json, length_histogram) }
        .join(ch_processed_read_quality_by_key, by: 0, remainder: true)
        .map { _profile_key, meta, profile_json, length_histogram, quality_histogram -> tuple(meta, profile_json, length_histogram, quality_histogram) }

    ch_processed_profiles_by_quality = ch_processed_read_profiles_with_quality.branch {
        with_quality: it[3] != null
        without_quality: true
    }

    PACKAGE_FASTX_REPORT(
        ch_processed_profiles_by_quality.without_quality
            .map { meta, profile_json, length_histogram, _quality_histogram ->
                tuple('fastx', meta.id, meta.profile_stage, meta.query_class ?: '', meta.producer ?: '', profile_json, length_histogram)
            }
            .mix(ch_filtered_contig_profiles.map { meta, profile_json, length_histogram, _sequence_count ->
                tuple('fastx', meta.id, meta.profile_stage, meta.query_class ?: '', meta.producer ?: '', profile_json, length_histogram)
            })
    )

    PACKAGE_FASTX_QUALITY_REPORT(
        ch_processed_profiles_by_quality.with_quality.map { meta, profile_json, length_histogram, quality_histogram ->
            tuple('fastx', meta.id, meta.profile_stage, meta.query_class ?: '', meta.producer ?: '', profile_json, length_histogram, quality_histogram)
        }
    )

    PACKAGE_ASSEMBLY_PROFILE_REPORT(
        ch_assembly_profiles.map { sample_id, _platform, _read_structure, producer, profile_json, length_histogram, _sequence_count ->
            tuple('assembly', sample_id, 'assembly', '', producer, profile_json, length_histogram)
        }
    )

    PACKAGE_ASSEMBLY_ELIGIBILITY_REPORT(
        ch_assembly_eligibility_decisions.map { meta, decision, sequence_count, threshold, reason ->
            tuple(meta.id, meta.producer, decision, sequence_count, threshold, reason)
        }
    )

    PACKAGE_ASSEMBLY_LONG_READ_ELIGIBILITY_REPORT(
        ch_long_read_eligibility_summaries.map { sample_id, summary_json -> tuple(sample_id, summary_json) }
    )

    PACKAGE_ASSEMBLY_UNION_REPORT(
        ch_long_read_union_summaries.map { sample_id, summary_json -> tuple(sample_id, summary_json) }
    )

    PACKAGE_NVD_PREPARED_QUERY_BATCHES_REPORT(
        ch_prepared_query_batch_summaries
    )

    PACKAGE_NVD_MEGABLAST_QUERY_PARTITION_REPORT(
        ch_megablast_query_partition_summaries
    )

    PACKAGE_NVD_TAXON_BIG_TABLE_REPORT(ch_taxon_big_tables)

    ch_target_enrichment_packages = PACKAGE_TARGET_ENRICHMENT_REPORT.out.packages.toList()
    ch_depletion_packages = PACKAGE_DEPLETION_REPORT.out.packages.toList()
    ch_fastx_packages = PACKAGE_FASTX_REPORT.out.packages.mix(PACKAGE_FASTX_QUALITY_REPORT.out.packages).toList()
    ch_assembly_packages = PACKAGE_ASSEMBLY_PROFILE_REPORT.out.packages
        .mix(PACKAGE_ASSEMBLY_ELIGIBILITY_REPORT.out.packages)
        .mix(PACKAGE_ASSEMBLY_LONG_READ_ELIGIBILITY_REPORT.out.packages)
        .mix(PACKAGE_ASSEMBLY_UNION_REPORT.out.packages)
        .toList()
    ch_query_preparation_packages = PACKAGE_NVD_PREPARED_QUERY_BATCHES_REPORT.out.packages.toList()
    ch_blast_packages = PACKAGE_NVD_MEGABLAST_QUERY_PARTITION_REPORT.out.packages.toList()
    ch_taxonomy_packages = PACKAGE_NVD_TAXON_BIG_TABLE_REPORT.out.packages.toList()

    BUILD_MULTIQC_INPUTS(
        ch_resolved_reads,
        ch_nvd_version,
        ch_fastqc_packages,
        ch_target_enrichment_packages,
        ch_depletion_packages,
        ch_fastx_packages,
        ch_assembly_packages,
        ch_query_preparation_packages,
        ch_blast_packages,
        ch_taxonomy_packages,
        ch_experimental_enabled,
        ch_target_enrichment_enabled,
        ch_depletion_enabled,
        ch_assembly_enabled,
        ch_blast_enabled,
    )

    emit:
    fastqc_zips = ch_fastqc_zips
    inputs = BUILD_MULTIQC_INPUTS.out.inputs
    config = ch_multiqc_config
}
