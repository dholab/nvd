/**
 * Utility functions for the NVD2 Nextflow pipeline.
 *
 * This class provides static helper methods for configuration validation
 * and other common operations. Methods are automatically available in
 * .nf files due to Nextflow's implicit import of classes in lib/.
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

    private static final List<String> LABKEY_GOTTCHA2_PARAMS = [
        'labkey_gottcha_fasta_list',
        'labkey_gottcha_full_list',
        'labkey_gottcha_blast_verified_full_list',
    ]

    // -------------------------------------------------------------------------
    // Public validation methods
    // -------------------------------------------------------------------------

    /**
     * Validates LabKey parameters required for the STAT+BLAST workflow.
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
        validateLabkeyParams(params, requiredParams, 'STAT+BLAST')
    }

    /**
     * Validates LabKey parameters required for the GOTTCHA2 workflow.
     * Checks common params plus GOTTCHA2-specific params.
     *
     * @param params The Nextflow params object
     * @throws IllegalStateException if labkey is enabled but required params are missing
     */
    public static void validateLabkeyGottcha2(params) {
        if (!params.labkey) {
            return
        }
        def requiredParams = LABKEY_COMMON_PARAMS + LABKEY_GOTTCHA2_PARAMS
        validateLabkeyParams(params, requiredParams, 'GOTTCHA2')
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
}
