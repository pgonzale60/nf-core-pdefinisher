class PdePublish {

    static String pathFor(Map params, Map meta, String defaultDir) {
        if (!params.per_species_outdir) {
            return defaultDir
        }
        if (!meta.assembly_name) {
            throw new IllegalArgumentException("per_species_outdir requires assembly_name in the samplesheet (sample: ${meta.id})")
        }
        def base = meta.species ?: meta.id
        return "${base}/analysis/${meta.assembly_name}/pdefinisher/${defaultDir}"
    }
}
