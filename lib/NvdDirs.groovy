/**
 * Directory resolution for NVD2 pipeline stateful and stateless modes.
 *
 * This class handles the logic for resolving state_dir, taxonomy_dir, and hits_dir
 * based on whether the pipeline is running in stateless mode or not. It validates
 * that required directories exist or can be created, failing early with clear
 * error messages if not.
 *
 * Usage in workflow:
 *     dirs = NvdDirs.resolve(params, log)
 *     // dirs.state_dir    - empty string "" in stateless mode, absolute path string otherwise
 *     // dirs.taxonomy_dir - absolute path string if explicit, empty string "" to derive from state_dir
 *     // dirs.hits_dir     - absolute path string for parquet output
 */
class NvdDirs {

    /** State directory path (empty string "" in stateless mode) */
    final String state_dir

    /** Taxonomy directory path (empty string "" if derived from state_dir) */
    final String taxonomy_dir

    /** Hits directory path for parquet output */
    final String hits_dir

    /** Whether running in stateless mode */
    final boolean stateless

    private NvdDirs(String state_dir, String taxonomy_dir, String hits_dir, boolean stateless) {
        this.state_dir = state_dir
        this.taxonomy_dir = taxonomy_dir
        this.hits_dir = hits_dir
        this.stateless = stateless
    }

    /**
     * Resolve directory paths based on pipeline parameters.
     *
     * In stateless mode:
     *   - state_dir is empty string "" (no state tracking)
     *   - taxonomy_dir is required and must exist
     *   - hits_dir is {results}/hits
     *
     * In stateful mode (default):
     *   - state_dir is created if needed
     *   - taxonomy_dir is optional (derived from state_dir if not set)
     *   - hits_dir is {state_dir}/hits
     *
     * @param nfParams The Nextflow params object
     * @param nfLog The Nextflow log object for info messages
     * @return NvdDirs instance with resolved paths
     * @throws AssertionError if required directories are missing or cannot be created
     */
    static NvdDirs resolve(nfParams, nfLog) {
        if (nfParams.stateless) {
            return resolveStateless(nfParams, nfLog)
        } else {
            return resolveStateful(nfParams, nfLog)
        }
    }

    // -------------------------------------------------------------------------
    // Private resolution methods
    // -------------------------------------------------------------------------

    private static NvdDirs resolveStateless(nfParams, nfLog) {
        // taxonomy_dir is required in stateless mode
        if (!nfParams.taxonomy_dir) {
            throw new AssertionError("""
                |ERROR: --taxonomy-dir is required when --stateless is enabled.
                |
                |In stateless mode, taxonomy data must be provided explicitly since
                |there is no state directory to derive it from.
                |
                |Example: nextflow run main.nf --stateless --taxonomy-dir /path/to/taxdump
                |""".stripMargin())
        }

        def taxonomy_dir_file = new File(nfParams.taxonomy_dir.toString())
        if (!taxonomy_dir_file.exists()) {
            throw new AssertionError("taxonomy_dir does not exist: ${nfParams.taxonomy_dir}")
        }
        def taxonomy_dir_str = taxonomy_dir_file.getAbsolutePath()

        // In stateless mode, hits go to results directory
        def results_file = new File(nfParams.results.toString())
        ensureDirectoryExists(results_file, 'results', nfParams.results.toString())
        def hits_dir_str = results_file.getAbsolutePath() + "/hits"

        nfLog.info "Running in STATELESS mode: no state tracking, taxonomy from ${taxonomy_dir_str}"

        // NOTE: We use "" instead of null because Channel.value(null) in Nextflow hangs forever.
        // Empty string is falsy in both Groovy and Python, so downstream conditionals work unchanged.
        return new NvdDirs(
            "",              // state_dir is empty string in stateless mode
            taxonomy_dir_str,
            hits_dir_str,
            true
        )
    }

    private static NvdDirs resolveStateful(nfParams, nfLog) {
        if (!nfParams.state_dir) {
            throw new AssertionError("params.state_dir must be set (check nextflow.config defaults)")
        }

        def state_dir_file = new File(nfParams.state_dir.toString())
        ensureDirectoryExists(state_dir_file, 'state', nfParams.state_dir.toString())
        def state_dir_str = state_dir_file.getAbsolutePath()

        // taxonomy_dir: use explicit param if set, otherwise empty string (scripts derive from state_dir)
        // NOTE: We use "" instead of null because Channel.value(null) in Nextflow hangs forever.
        // Empty string is falsy in both Groovy and Python, so downstream conditionals work unchanged.
        def taxonomy_dir_str = ""
        if (nfParams.taxonomy_dir) {
            def taxonomy_dir_file = new File(nfParams.taxonomy_dir.toString())
            if (!taxonomy_dir_file.exists()) {
                throw new AssertionError("taxonomy_dir does not exist: ${nfParams.taxonomy_dir}")
            }
            taxonomy_dir_str = taxonomy_dir_file.getAbsolutePath()
        }

        // hits_dir is inside state_dir
        def hits_dir_str = state_dir_str + "/hits"

        return new NvdDirs(
            state_dir_str,
            taxonomy_dir_str,
            hits_dir_str,
            false
        )
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Ensure a directory exists, creating it if necessary.
     * Fails with a clear error message if creation fails.
     */
    private static void ensureDirectoryExists(File dir, String dirType, String originalPath) {
        if (dir.exists()) {
            return
        }

        def created = dir.mkdirs()
        if (!created) {
            throw new AssertionError("""
                |ERROR: Cannot create ${dirType} directory: ${originalPath}
                |
                |Check that you have write permissions to the parent directory.
                |""".stripMargin())
        }
    }
}
