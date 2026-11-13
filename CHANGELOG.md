# nf-core/pdefinisher: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [unreleased<!-- TODO nf-core: replace with date on release -->]

### `Added`

- `DELFIES` publishes `breakpoint_foci__S2G.tsv`; `PDE_IDENTIFYPUTATIVEGRS` builds breakpoints from its per-strand support instead of the collapsed `breakpoint_locations.bed`, so reciprocal micro-GRS cuts are kept (#37).
- `recon_motif.nf`: entry running `PDE_RECONNAISSANCE` and `PDE_MOTIF_DISCOVERY` (GRS, DELFIES, MEME/FIMO) without assembly surgery, evaluation or validation (#27).
- `PDE_MOTIF_DISCOVERY`: core grid and auxiliary stage synced from `nf-core-cbmdiscovery` (SSOT-driven radii/widths, header-derived distal background, `REJECTED_TELOMERIC_ECHO` filter); aux telomere repeat taken from `meta.telomere_repeat` (#31).
- Per-species publishing now writes to `<species>/analysis/<assembly_name>/pdefinisher/<dir>`; new optional samplesheet column `assembly_name`, required when `per_species_outdir` is set (#35).
- `recon_motif.nf` is the standalone recon + motif entry (assembly + BAM in; mosdepth outputs, CBM motifs and DELFIES S2G out); the reconnaissance core now runs `MOSDEPTHCLASS` and no longer runs `SOFTCLIPENDS`, and mosdepth outputs are published (#70).

## v1.0.0 - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/pdefinisher, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
