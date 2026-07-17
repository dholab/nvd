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
        'labkey_exp_id_guard_list',
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
