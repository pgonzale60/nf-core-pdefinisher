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

Prepare a comma-separated `samplesheet.csv` with your sample ID, input draft assembly FASTA, and long reads (FASTA/FASTQ):

```csv
sample,assembly,reads
oscheius_onirici,/path/to/draft_assembly.fasta.gz,/path/to/hifi_reads.fasta.gz
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

## Graph finishing pilot

`graph_finish.nf` runs `gfapdefinisher` against an existing Hifiasm primary GFA and compares its proposed ends with the curated assembly. It also projects the telomere-family consensuses onto the primary graph, the optional raw graph, and the curated assembly. The evidence bundle keeps cluster membership, original molecule IDs, local haplotype patterns, and competing alignments. A consensus alignment identifies a candidate and cannot establish a physical germline closure.

```bash
nextflow run graph_finish.nf --curated curated.fa.gz --gfa asm.p_utg.gfa.gz \
    --end_archive path/to/end_clusters --chromosomes I,II,III,IV,V,X --outdir results
```

Use a worker environment with Python 3, `gfapdefinisher`, minimap2, and seqkit. Set `--gfapdefinisher` to select an installed checkout or executable. The end archive is the output directory of the telomeric read clustering workflow. Optional graph inputs are `--raw_gfa`, `--tiara`, `--barrnap`, `--ont_reads`, and `--ont_paf`. An ONT PAF must use the supplied primary unitig graph as its reference. The raw GFA and `--hifi_reads` can also supply sequence used by a reviewed recipe.

Discovery publishes the curated FASTA byte-for-byte with AGP and graph proposals for review. Adding `--recipe` replays a reviewed graph assembly. The replay checks its component source intervals against the spelled FASTA, publishes the source plan, and round-trips the delivered AGP and component FASTA. See `evidence/curated_end_review.tsv` for assignment and continuity and `evidence/graph_end_proposals.tsv` for the graph candidates. This pilot requires a prebuilt GFA and runs independently of the existing surgery and motif workflows.

## Credits

`nf-core/pdefinisher` was originally written by **Pablo M. Gonzalez de la Rosa**.

## Citations

If you use `nf-core/pdefinisher` for your analysis, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) and see [`CITATIONS.md`](CITATIONS.md) for the full list of tools used by the pipeline.
