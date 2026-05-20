/**
 * Directory resolution for explicit NVD pipeline paths.
 *
 * The v3 pipeline no longer carries a generic state directory through the
 * Nextflow graph. Taxonomy is an explicit runtime reference because cluster
 * worker nodes need a path they can access directly.
 *
 * Usage in workflow:
 *     dirs = NvdDirs.resolve(params, log)
 *     // dirs.taxonomy_dir - absolute path string
 */
class NvdDirs {

    /** Taxonomy directory path */
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
     * @throws AssertionError if the taxonomy directory is missing or invalid
     */
    static NvdDirs resolve(nfParams, nfLog) {
        if (!nfParams.taxonomy_dir) {
            throw new AssertionError(
                "taxonomy_dir is required. Set params.taxonomy_dir in user.config " +
                "or pass --taxonomy-dir with a directory visible to worker nodes."
            )
        }

        def taxonomy_dir_file = new File(nfParams.taxonomy_dir.toString())
        if (!taxonomy_dir_file.exists()) {
            throw new AssertionError("taxonomy_dir does not exist: ${nfParams.taxonomy_dir}")
        }
        def taxonomy_dir_str = taxonomy_dir_file.getAbsolutePath()

        return new NvdDirs(taxonomy_dir_str)
    }
}
