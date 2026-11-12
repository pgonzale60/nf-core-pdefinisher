class PdePublish {

    static String pathFor(Map params, Map meta, String defaultDir) {
        if (!params.per_species_outdir) {
            return defaultDir
        }
        def specimen = meta.tolid ?: (meta.specimen_id ?: meta.id)
        if (meta.species && specimen && meta.species != specimen) {
            return "${meta.species}/${specimen}/analysis/pdefinisher/${defaultDir}"
        } else if (meta.species) {
            return "${meta.species}/analysis/pdefinisher/${defaultDir}"
        } else {
            return "${meta.id}/analysis/pdefinisher/${defaultDir}"
        }
    }
}
