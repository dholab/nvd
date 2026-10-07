/**
 * Utility functions for the NVD Nextflow pipeline.
 *
 * This class provides static helper methods for LabKey configuration
 * validation. Methods are automatically available in .nf files due to
 * Nextflow's implicit import of classes in lib/.
 */
class NvdUtils {

    // -------------------------------------------------------------------------
    // LabKey parameter definitions
    // -------------------------------------------------------------------------

    private static final List<String> LABKEY_COMMON_PARAMS = [
        'labkey_server',
        'labkey_project_name',
        'labkey_webdav',
        'labkey_schema',
    ]

    private static final List<String> LABKEY_BLAST_PARAMS = [
        'labkey_blast_meta_hits_list',
        'labkey_blast_fasta_list',
    ]

    // -------------------------------------------------------------------------
    // Public methods
    // -------------------------------------------------------------------------

    /**
     * Validates LabKey parameters required for NVD BLAST reporting.
     * Checks common params plus BLAST-specific params.
     *
     * @param params The Nextflow params object
     * @throws IllegalStateException if labkey is enabled but required params are missing
     */
    public static void validateLabkeyBlast(params) {
        if (!params.labkey) {
            return
        }
        def requiredParams = LABKEY_COMMON_PARAMS + LABKEY_BLAST_PARAMS
        validateLabkeyParams(params, requiredParams, 'NVD BLAST reporting')
    }

    /**
     * Returns true when any target-enrichment index source has been configured.
     */
    public static boolean hasTargetEnrichmentIndex(params) {
        return params.virus_index || params.virus_index_url || params.virus_reference_fasta
    }

    /**
     * Resolve the effective target enrichment mode.
     *
     * Target enrichment is enabled when an index source is configured unless
     * no_enrichment explicitly disables it.
     */
    public static boolean targetEnrichmentEnabled(params) {
        def disabled = parseOptionalBool(params.no_enrichment) ?: false
        return hasTargetEnrichmentIndex(params) && !disabled
    }

    /**
     * Resolve whether host/contaminant depletion is configured.
     */
    public static boolean depletionEnabled(params) {
        return params.host_index || params.host_index_url || params.host_contaminants_fasta
    }

    /**
     * Returns true when a background depletion index is configured. Its
     * presence switches the first Deacon pass from enrichment to depletion.
     */
    public static boolean backgroundDepletionEnabled(params) {
        return params.background_index ? true : false
    }

    /**
     * Resolve the single filter that step 1 (and the contig filter) runs.
     *
     * Enrichment wins when enabled; otherwise a configured background index
     * selects depletion with the background thresholds; otherwise the empty
     * passthrough index is used and the thresholds are irrelevant.
     */
    public static Map stepOneFilterPolicy(params) {
        def enrichment = targetEnrichmentEnabled(params)
        def background = !enrichment && backgroundDepletionEnabled(params)
        return [
            target_enrichment_enabled: enrichment,
            background_depletion_enabled: background,
            abs_threshold: background ? params.background_abs_threshold : params.virus_abs_threshold,
            rel_threshold: background ? params.background_rel_threshold : params.virus_rel_threshold,
        ].asImmutable()
    }

    /**
     * Policy for DEACON_FILTER_CONTIGS: the same step-one filter the reads
     * saw, plus the optional host/contaminant depletion pass.
     */
    public static Map contigFilterPolicy(params, boolean use_depletion) {
        def step_one = stepOneFilterPolicy(params)
        return [
            target_enrichment_enabled: step_one.target_enrichment_enabled,
            target_abs_threshold: step_one.abs_threshold,
            target_rel_threshold: step_one.rel_threshold,
            depletion_enabled: use_depletion,
            depletion_abs_threshold: use_depletion ? params.host_abs_threshold : null,
            depletion_rel_threshold: use_depletion ? params.host_rel_threshold : null,
        ].asImmutable()
    }

    /**
     * Stop the run when step 1 cannot resolve to one filter.
     *
     * @throws IllegalStateException on a background index alongside enabled
     *         target enrichment, a bare --background_index flag, or a missing
     *         index file.
     */
    public static void validateStepOneFilter(params) {
        def background = params.background_index
        if (background == null || background == false) {
            return
        }
        if (!(background instanceof CharSequence)) {
            throw new IllegalStateException(
                "background_index must be a path to a prebuilt Deacon .idx file; received '${background}'. " +
                "Pass --background_index /path/to/background.idx."
            )
        }
        def path = new File(background.toString())
        if (!path.isFile()) {
            throw new IllegalStateException(
                "background_index points to a file that does not exist: ${background}"
            )
        }
        if (targetEnrichmentEnabled(params)) {
            def sources = ['virus_index', 'virus_index_url', 'virus_reference_fasta'].findAll { name -> params[name] }
            throw new IllegalStateException(
                "background_index cannot be combined with target enrichment (${sources.join(', ')} also set). " +
                "Step 1 runs one Deacon filter. Pass --no_enrichment true to deplete the background instead, " +
                "or drop background_index."
            )
        }
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Core validation logic shared by workflow-specific validators.
     *
     * @param params The Nextflow params object
     * @param requiredParams List of parameter names to validate
     * @param workflowName Name of the workflow for error messages
     * @throws IllegalStateException if any required params are missing
     */
    private static void validateLabkeyParams(params, List<String> requiredParams, String workflowName) {
        def paramValues = requiredParams.collectEntries { name ->
            [(name): params[name]]
        }

        def missing = paramValues.findAll { k, v -> v == null }.keySet()

        if (missing.isEmpty()) {
            return
        }

        def table = paramValues.collect { k, v ->
            def displayVal = v != null ? v.toString() : '(not set)'
            if (displayVal.length() > 42) {
                displayVal = displayVal.take(39) + '...'
            }
            "| ${k.padRight(40)} | ${displayVal.padRight(42)} |"
        }.join('\n            |')

        def message = """
            |
            |LabKey integration enabled (--labkey) but some required parameters are not set.
            |Workflow: ${workflowName}
            |
            |Required LabKey configuration for ${workflowName}:
            |+------------------------------------------+--------------------------------------------+
            || Parameter                                | Value                                      |
            |+------------------------------------------+--------------------------------------------+
            |${table}
            |+------------------------------------------+--------------------------------------------+
            |
            |Missing: ${missing.join(', ')}
            |
            |Set these in your user.config (~/.nvd/user.config), a preset, or via CLI flags.
            |See: nvd run --help
            |""".stripMargin()

        throw new IllegalStateException(message)
    }

    private static Boolean parseOptionalBool(value) {
        if (value == null) {
            return null
        }
        if (value instanceof Boolean) {
            return value
        }
        def normalized = value.toString().trim().toLowerCase()
        if (['true', '1', 'yes', 'y', 'on'].contains(normalized)) {
            return true
        }
        if (['false', '0', 'no', 'n', 'off'].contains(normalized)) {
            return false
        }
        throw new IllegalArgumentException("Expected boolean-like value for no_enrichment, got '${value}'")
    }
}
