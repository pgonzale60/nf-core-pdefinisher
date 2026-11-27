# nf-core/pdefinisher

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)

## Introduction

**nf-core/pdefinisher** is a bioinformatics pipeline for **automated T2T finishing, gap resolution, and programmed DNA elimination (PDE) analysis** in species undergoing chromatin diminution.

In species with PDE, mixed germline/somatic sequencing libraries produce complex assembly bifurcations at germline-soma elimination boundaries, internal scaffold gaps across eliminated satellites, uncapped neo-telomeric chromosome ends, and unplaced germline-restricted sequences (GRS). `nf-core/pdefinisher` automates the end-to-end finishing, curation, and validation of these assemblies from a draft genome FASTA and long-read sequencing data (PacBio HiFi / ONT).

### Pipeline Architecture

1. **`PDE_RECONNAISSANCE`** — Parallel evidence collection on the input assembly:
   - Long-read alignment (`minimap2`), read-depth profiling (`mosdepth`), structural variant calling (`Sniffles2`), neo-telomere breakpoint detection (`DELFIES`), ribosomal DNA profiling (`barrnap`), tandem repeat masking (`RED`), and unplaced-to-chromosome whole-genome alignment.
   - Putative Germline-Restricted Sequence (`IDENTIFYPUTATIVEGRS`) classification from multi-modal coverage, telomere, and somatic/germline read signatures.
2. **`CURATIONSYNTHESIS`** — Automated synthesis of evidence into an explicit curation plan (`curation_plan.json`) identifying scaffold gaps, uncapped chromosome terminals, rDNA array fractures, and trans-chromosomal germline links.
3. **`CURATION_RESOLUTION`** — Targeted local de novo assembly and graph walking:
   - `TARGETEDGAPASSEMBLY`: Local `hifiasm` assembly and graph bridging across internal scaffold gaps.
   - `TELOMEREWALKER`: Iterative long-read walking and capping of uncapped chromosome ends to telomere-to-telomere (T2T) or neo-telomeric termination.
   - `RDNACASSETTERESOLVER` & `TRANSCHROMOSOMALLINKRESOLVER`: Resolution of rDNA array boundaries and germline-bridged chromosome ends.
4. **`ASSEMBLY_SURGERY`** — Deterministic AGP and FASTA reconstruction (`CURATIONAGPGENERATOR`) governed by declarative policy thresholds (`assets/curation_policy.json`).
5. **`CURATION_EVALUATION` & `CURATION_VALIDATION`** — Post-surgery audit of gap closures, telomere extensions, read-support verification, and multi-evidence somatic/germline architecture classification (`CURATIONSUMMARY`).
6. **`PDE_CURATED_VALIDATION`** — Independent re-alignment (`minimap2`), coverage (`mosdepth`), neo-telomere (`DELFIES`), and GRS validation on the final curated assembly.

## Usage

### Samplesheet Input

Provide one row per sample in a CSV or TSV. For BAM-based curation, supply the assembly and reads aligned to that exact assembly:

```csv
sample,assembly,bam,bai,chromosomes,assembly_name
oscheius_onirici,/path/to/curated.fa.gz,/path/to/hifi.bam,/path/to/hifi.bam.bai,I;II;III;IV;V;X,nxOscOnir1.1
```

For graph-based curation, supply the whole-genome primary Hifiasm GFA, phased telomere-family archive, and nuclear chromosome names:

```csv
sample,assembly,gfa,end_archive,chromosomes
DF5120,/path/to/curated.fa.gz,/path/to/asm.p_utg.gfa.gz,/path/to/end_clusters,I;II;III;IV;V;X
```

### Running the Pipeline

```bash
nextflow run nf-core/pdefinisher \
   -profile singularity \
   --input samplesheet.csv \
   --telomere_repeat 'TTAGGC' \
   --outdir ./results
```

Optional parameters:

- `--telomere_repeat`: Canonical telomeric monomer repeat (default: `'TTAGGC'`).
- `--chromosomes`: Optional comma-separated list of primary chromosomes to restrict curation to (default: auto-detected from the assembly).
- `--curation_policy`: Path to custom JSON curation policy (default: `assets/curation_policy.json`).

## Recon + motif only

To run mosdepth coverage, DELFIES S2G, GRS segment definition and CBM motif discovery without assembly surgery, rDNA assembly, read fishing or softclip analysis:

```bash
nextflow run recon_motif.nf --input samplesheet.csv --outdir results
```

## Package a reviewed assembly

Run `python3 scripts/package_reviewed_assembly.py --help` to package an accepted
FASTA without editing it. The package contains the original FASTA, AGP, component
FASTA, an end inventory, unresolved join coordinates, and checksums. The script
verifies that the AGP and components reconstruct the input sequence. Components
are slices of the reviewed FASTA. An optional reviewed recipe can preserve source
graph and read evidence if its expected FASTA checksum matches.

## Graph-based curation

The main workflow runs `gfapdefinisher` for rows with `gfa` and `end_archive`. Run it with an environment that provides Python 3, `gfapdefinisher`, minimap2, and seqkit. Set `--gfapdefinisher` to select a specific installed executable. The archive is the output directory of the telomeric read clustering workflow. Optional samplesheet columns are `raw_gfa`, `tiara`, `barrnap`, `ont_reads`, `ont_paf`, `hifi_reads`, and `graph_recipe`. The ONT PAF must use the primary unitig graph as its reference. A reviewed recipe can draw sequence from the raw GFA, the curated FASTA, and original HiFi or ONT reads.

`gfapdefinisher` proposes paths from the graph. The workflow projects telomere-family consensuses onto both graphs and the curated assembly, records competing placements, and retains local haplotype patterns and original molecule IDs. It passes the [phased-end evidence bundle](docs/phased_end_evidence.md) to executables that support `--end-evidence`. With older executables the workflow associates the phased families with the graph candidates for review. A family consensus alignment by itself establishes a candidate placement.

Without `graph_recipe`, the workflow delivers the curated FASTA with an AGP and the candidate report. With a recipe, it verifies the reviewed nuclear assembly and retains curated records outside the recipe, such as the mitochondrial scaffold. The workflow round-trips the delivery FASTA against its component AGP. BAM-based surgery runs for rows with BAM input alone.

To collect both sources of proposals for the same sample, supply `assembly_name` and a BAM aligned to the input curated FASTA on the graph row. The BAM resolvers run before surgery and contribute their blueprint and manifests to `curation_edit_plan.json`. The plan records curated-reference coordinates, competing graph placements, and the decision for each proposal. Optional motif discovery uses the same input FASTA and BAM coordinates.

### Approve HiFi gap patches

Review `evidence/hifi_gap_patches.json` and the sequence payloads in `evidence/hifi_patch_sequences.fasta`. On the next run, set the samplesheet column `hifi_patch_approvals` to a JSON file containing the curated reference checksum and the IDs and sequence checksums you approve. Copy these values from the gap-patch report. Run with `-resume` to reuse the BAM resolver tasks.

```json
{
  "version": 1,
  "reference_sha256": "<curated reference file SHA-256>",
  "patches": [{"gap_id": "CARD-02", "patch_sha256": "<approved patch sequence SHA-256>"}]
}
```

The workflow applies approved patches after any reviewed graph replay. A patch must target an internal N-gap in the curated reference. At least one resolver-reported flank must uniquely anchor the locus in each assembly, and both flanks must delimit an intact gap in the accepted assembly in the expected order and orientation. The workflow records patches that are already present or withheld, including attempts to replace accepted graph sequence, anchored rDNA placeholders, or unresolved satellite ports. The final FASTA and AGP include the applied patches, and the assembly plan records their source and final coordinates.

The graph evidence bundle includes `curated_end_review.tsv`, `graph_end_proposals.tsv`, and `family_placements.tsv`. `sequence_gaps.tsv` records both linkage and missing sequence. For an anchored rDNA array the flanking connection is supported, while the sequence and copy order within the array remain unresolved. `unresolved_joins.tsv` records gaps with unresolved or unassessed linkage, including tandem repeats with competing distal placements.

## Credits

`nf-core/pdefinisher` was originally written by **Pablo M. Gonzalez de la Rosa**.

## Citations

If you use `nf-core/pdefinisher` for your analysis, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) and see [`CITATIONS.md`](CITATIONS.md) for the full list of tools used by the pipeline.
