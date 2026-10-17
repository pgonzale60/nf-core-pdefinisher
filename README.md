# nf-core/pdecuration

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)

## Introduction

**nf-core/pdecuration** is a bioinformatics pipeline for **automated genome assembly curation and validation in species undergoing Programmed DNA Elimination (PDE)** (chromatin diminution).

In species with PDE, mixed germline/somatic sequencing libraries produce complex assembly bifurcations at germline-soma elimination boundaries, internal scaffold gaps across eliminated satellites, uncapped neo-telomeric chromosome ends, and unplaced germline-restricted sequences (GRS). `nf-core/pdecuration` automates the end-to-end curation and validation of these assemblies from a draft genome FASTA and long-read sequencing data (PacBio HiFi / ONT).

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
nextflow run nf-core/pdecuration \
   -profile singularity \
   --input samplesheet.csv \
   --telomere_repeat 'TTAGGC' \
   --outdir ./results
```

Optional parameters:
- `--telomere_repeat`: Canonical telomeric monomer repeat (default: `'TTAGGC'`).
- `--chromosomes`: Optional comma-separated list of primary chromosomes to restrict curation to (default: auto-detected from the assembly).
- `--curation_policy`: Path to custom JSON curation policy (default: `assets/curation_policy.json`).

## Credits

`nf-core/pdecuration` was originally written by **Pablo M. Gonzalez de la Rosa** (Wellcome Sanger Institute).

## Citations

If you use `nf-core/pdecuration` for your analysis, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) and see [`CITATIONS.md`](CITATIONS.md) for the full list of tools used by the pipeline.
