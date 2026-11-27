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

## Audit chromosome ends on a curated assembly

`scripts/audit_chromosome_ends.py` takes a curated FASTA, a matching indexed long-read
BAM, and SatelliteProfiler arrays and family monomers. Use `--help` for the required
chromosome names, telomere motif, and optional DELFIES and graph inputs. The script
requires Python with pysam and reports in the input assembly's coordinates.

`ends.tsv` lists physical tip tracts, read support, gaps in the end window, and
places to review. `end_repeat_links.tsv` gives the replicated diagnostic sequence
footprint at each end. `repeat_families.tsv` counts ends with anchored molecules
for the families in the supplied panel. A footprint measures bases observed on a
read and provides a lower bound. `unplaced_ending_candidates.tsv` reports capped
assembly or graph fragments sharing a repeat family. Candidate fragments may be
compatible with several ends. The evidence tables retain read or fragment IDs so
the underlying sequence can be inspected.

The audit writes reports only. The assembly, alignments, and graph remain inputs.

## Credits

`nf-core/pdefinisher` was originally written by **Pablo M. Gonzalez de la Rosa**.

## Citations

If you use `nf-core/pdefinisher` for your analysis, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) and see [`CITATIONS.md`](CITATIONS.md) for the full list of tools used by the pipeline.
