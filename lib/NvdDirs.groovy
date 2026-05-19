/**
 * Directory resolution for explicit NVD pipeline paths.
 *
 * The v3 pipeline no longer carries a generic state directory through the
 * Nextflow graph. Taxonomy cache fallback resolution belongs to the Python
 * taxonomy helpers; Nextflow only validates and forwards an explicit
 * taxonomy_dir when the user supplies one.
 *
 * Usage in workflow:
 *     dirs = NvdDirs.resolve(params, log)
 *     // dirs.taxonomy_dir - absolute path string if explicit, empty string "" otherwise
 */
class NvdDirs {

    /** Taxonomy directory path (empty string "" if Python should resolve fallback) */
    final String taxonomy_dir

    private NvdDirs(String taxonomy_dir) {
        this.taxonomy_dir = taxonomy_dir
    }

    /**
     * Resolve explicit taxonomy path from pipeline parameters.
     *
     * @param nfParams The Nextflow params object
     * @param nfLog The Nextflow log object for info messages
     * @return NvdDirs instance with resolved paths
     * @throws AssertionError if an explicit taxonomy directory is missing
     */
    static NvdDirs resolve(nfParams, nfLog) {
        if (!nfParams.taxonomy_dir) {
            nfLog.info "Taxonomy directory not provided; Python taxonomy preflight will resolve fallback location"
            return new NvdDirs("")
        }

        def taxonomy_dir_file = new File(nfParams.taxonomy_dir.toString())
        if (!taxonomy_dir_file.exists()) {
            throw new AssertionError("taxonomy_dir does not exist: ${nfParams.taxonomy_dir}")
        }
        def taxonomy_dir_str = taxonomy_dir_file.getAbsolutePath()

        nfLog.info "Using explicit taxonomy directory: ${taxonomy_dir_str}"
        return new NvdDirs(taxonomy_dir_str)
    }
}
