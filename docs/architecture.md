# Architectural Upgrade Specification (`v2-tool-first`) — `nf-core/pdecuration`

## 1. Executive Summary & Branch Discipline

1. **Preserve `master` (`0129e8b`) Untouched**:
   - `master` is our verified, passing baseline. Before making any edits, create and switch to a new feature branch:
     ```bash
     git checkout -b v2-tool-first
     ```
2. **Goal of `v2-tool-first`**:
   - Replace custom Python sequence-manipulation wheels across `PDE_RECONNAISSANCE`, `CURATION_RESOLUTION`, and `ASSEMBLY_SURGERY` with official `nf-core/modules` and verified BioContainers (`bbduk.sh`, `vsearch`, `teloclip`, `abpoa`, `mindthegap`, `gfatools`, `moddotplot`, and `ragtag`), while keeping the PDE-specific biological policy/synthesis architecture (`CURATIONSYNTHESIS` blueprint, `curation_policy.json`, and Sanger AGP v2.1 audit) intact.

---

## 2. Verified `nf-core/modules` & Galaxy Singularity Depot Container Matrix

**Hard Rule**: Never guess or invent container tag hashes. Every container below has been verified (`HTTP 200`) against Seqera Community Wave or `https://depot.galaxyproject.org/singularity/`:

| Tool | Source | Verified Singularity Image | Verified Docker Image |
| :--- | :--- | :--- | :--- |
| **`bbmap/bbduk`** | Official `nf-core/modules` | `'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data'` | `'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9'` |
| **`samtools` + `bbmap` + `pigz`** | Seqera Wave | `'community.wave.seqera.io/library/samtools_bbmap_pigz:c8f695696758ff77'` | `'community.wave.seqera.io/library/samtools_bbmap_pigz:c8f695696758ff77'` |
| **`vsearch/cluster`** | Official `nf-core/modules` | `'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/363cb3f920602951d1dee4e0d248eb7171b4dfd2611db3a581a480500ac25996/data'` | `'community.wave.seqera.io/library/samtools_vsearch:36d56b7ebe9ea248'` |
| **`samtools` + `minimap2` + `vsearch` + `python`** | Seqera Wave | `'community.wave.seqera.io/library/samtools_minimap2_vsearch_python:1aa370c14225a820'` | `'community.wave.seqera.io/library/samtools_minimap2_vsearch_python:1aa370c14225a820'` |
| **`gfatools`** | Official `nf-core/modules` | `'https://depot.galaxyproject.org/singularity/gfatools:0.5--h577a1d6_5'` | `'quay.io/biocontainers/gfatools:0.5--h577a1d6_5'` |
| **`ragtag`** | Official `nf-core/modules` | `'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e783bbc5474a09f51a0f93142fb94a864715e77b40d9b33e11b53b9c6c4d08cf/data'` | `'community.wave.seqera.io/library/ragtag:2.1.0--6125126fe9d553ab'` |
| **`teloclip`** | BioContainers / Galaxy Depot | `'https://depot.galaxyproject.org/singularity/teloclip:0.4.0--pyhdfd78af_0'` | `'quay.io/biocontainers/teloclip:0.4.0--pyhdfd78af_0'` |
| **`teloclip` + `samtools`** | Seqera Wave | `'community.wave.seqera.io/library/teloclip_samtools:2e39bcbf9e597d36'` | `'community.wave.seqera.io/library/teloclip_samtools:2e39bcbf9e597d36'` |
| **`mindthegap`** | BioContainers / Galaxy Depot | `'https://depot.galaxyproject.org/singularity/mindthegap:2.3.0--hdcf5f25_3'` | `'quay.io/biocontainers/mindthegap:2.3.0--hdcf5f25_3'` |
| **`moddotplot`** | BioContainers / Galaxy Depot | `'https://depot.galaxyproject.org/singularity/moddotplot:0.9.9--pyhdfd78af_0'` | `'quay.io/biocontainers/moddotplot:0.9.9--pyhdfd78af_0'` |
| **`hifiasm` + `abpoa` + `minimap2` + `samtools` + `python`** | Existing `TARGETEDGAPASSEMBLY` | `'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_python:7a9feac98901a9ac'` | `'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_python:7a9feac98901a9ac'` |
| **`hifiasm` + `abpoa` + `minimap2` + `samtools` + `mindthegap` + `python`** | Seqera Wave | `'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_mindthegap_python:4e99aa7a05e907bf'` | `'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_mindthegap_python:4e99aa7a05e907bf'` |
| **`hifiasm` + `abpoa` + `minimap2` + `samtools` + `gfatools` + `python`** | Seqera Wave | `'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_gfatools_python:b88b84d24c308348'` | `'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_gfatools_python:b88b84d24c308348'` |

---

## 3. The 5 Architectural Upgrades to Implement on `v2-tool-first`

### Upgrade 1: Replace `seqkit grep` with `bbduk.sh` in `EXTRACTTELOMERICREADS`
- **Current bottleneck**: [`modules/local/extracttelomericreads/main.nf`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-pdecuration/modules/local/extracttelomericreads/main.nf) pipes `samtools fastq` into `seqkit grep -s -i -p`, performing slow linear regex scanning across the whole-genome BAM.
- **Required change**:
  - Pipe `samtools fastq` into `bbduk.sh` (`in=stdin.fq outm=${prefix}.telomeric_reads.fastq.gz literal="\$fwd_repeat" k=18 rcomp=t hdist=0 threads=${task.cpus}`) or add an upstream `BBMAP_BBDUK` step (note: when adapting `nf-core/bbmap/bbduk`, use `outm=` rather than `out=` so only matched reads are emitted).

### Upgrade 2: Replace Custom Python $O(N^2)$ Clustering in `TELOMERICJUNCTIONPROFILER` with `vsearch`
- **Current bottleneck**: In [`pde_telomeric_junction_profiler.py:159-212`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-pdecuration/modules/local/telomericjunctionprofiler/resources/usr/bin/pde_telomeric_junction_profiler.py#L159-L212), `cluster_junction_flanks()` uses a naive $O(N^2)$ Python Hamming loop on the first 60 bp (ignoring indels) and picks a single raw read (`flanks[len(flanks) // 2]`) to map back to the reference with `minimap2`.
- **Required change**:
  - Replace `cluster_junction_flanks()` with a call to **`vsearch`** (requiring `shutil.which("vsearch")` with no silent fallback):
    ```bash
    vsearch \
        --cluster_fast raw_flanks.fasta \
        --id 0.97 \
        --iddef 0 \
        --qmask none \
        --strand plus \
        --consout clusters_cons.fasta \
        --uc clusters.uc \
        --sizeout \
        --threads ${threads}
    ```
  - Align the error-corrected MSA consensus sequences (`clusters_cons.fasta`) to `ref_fasta` using `minimap2`.

### Upgrade 3: Replace Custom Python Overhang Walking in `TELOMEREWALKER` with `teloclip` + `abpoa` + `gfatools`
- **Current bottleneck**: [`pde_telomere_walker.py`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-pdecuration/modules/local/telomerewalker/resources/usr/bin/pde_telomere_walker.py) and `resources/usr/lib/telowalking/` (>1,500 lines of custom Python) implement a custom `UnmappedReadIndex`, iterative prefix-suffix read walking (`run_iterative_walker`), and single-read tip overhang selection (`maximal_tip_overhang`).
- **Required change**:
  1. **Switch `TELOMEREWALKER`'s container** to `'community.wave.seqera.io/library/hifiasm_abpoa_minimap2_samtools_python:7a9feac98901a9ac'` (which already includes `abpoa`, `hifiasm`, `minimap2`, `samtools`, and `python`), and add a dedicated `TELOCLIP_EXTRACT` module (`https://depot.galaxyproject.org/singularity/teloclip:0.4.0--pyhdfd78af_0`) in `CURATION_RESOLUTION` (or `PDE_RECONNAISSANCE`) that runs:
     ```bash
     teloclip filter --ref-idx ${ref}.fai --min-clip 100 --max-break 50 --motifs ${telomere_repeat} ${bam} \
       | samtools view -bS - | samtools sort -o ${prefix}.teloclip.bam
     teloclip extract --ref-idx ${ref}.fai --min-clip 100 --out-dir ${prefix}_teloclip_overhangs ${prefix}.teloclip.bam
     ```
  2. **Consensus Overhang Generation via `abpoa` (`telocapping.py`)**:
     - Instead of picking a single raw read in `maximal_tip_overhang`, run **`abpoa <overhangs.fa> -r 0`** on the extracted terminal soft-clipped overhangs to compute the SIMD-POA consensus extension (which naturally follows the dominant 85% somatic telomere branch at PDE forks and corrects raw read errors).
  3. **Subtelomeric Repeat Array Bridging via `hifiasm` Unitigs / `gfatools` (`telowalking.py`)**:
     - Replace the custom Python iterative read-to-read prefix-suffix walker (`UnmappedReadIndex` / `run_iterative_walker`) with local unitig bridging from `telo_unitigs` (`HIFIASM`) and `abpoa` consensus stitching, cutting `telowalking.py` down by >70%.

### Upgrade 4: Add Dentist's Improper-Alignment Guard ($I_{\max}$) & Tier-2 `MindTheGap fill` to `TARGETEDGAPASSEMBLY`
- **Why**:
  1. **Keep Tier 1 (`minimap2` + relative strand check + `abpoa`)** for gaps physically spanned by single 15 kb HiFi reads (because full-length spanning reads bridge $2\text{–}5\text{ kb}$ tandem repeats that collapse fixed-$k$ de Bruijn graphs).
  2. **Add Dentist's Improper-Alignment Guard ($I_{\max}$) before bridging**:
     - At each gap boundary, count reads aligned to the flank (`MAPQ >= 20`) that exhibit $>100\text{ bp}$ of unaligned soft-clipping diverging from the opposite flank ($I_{\text{improper}}$).
     - Compare $I_{\text{improper}}$ against Dentist's threshold:
       $$I_{\max} = \lfloor 0.5 \cdot x + \exp(0.1875 \cdot (8 - x)) \rfloor$$
       where $x$ is the local read depth at the flank. If $I_{\text{improper}} > I_{\max}$ on both flanks and no high-confidence single-copy anchor spans the repeat, flag the gap as `REPEAT_BOUNDARY_VETO` to prevent ectopic transposable-element chimeras.
  3. **Add `MindTheGap fill` (`mindthegap:2.3.0--hdcf5f25_3`) for Tier-2 Multi-Read Gaps (`n_spanning == 0`)**:
     - When no single HiFi read spans both flanks (`n_spanning == 0`), use `MindTheGap fill` (`mindthegap fill -in baited_reads.fq -bkpt gap_breakpoints.fa -kmer-size 63 -max-nodes 500`) alongside `hifiasm` unitig bridging to find bounded de Bruijn paths between the Left and Right flank anchors.

### Upgrade 5: Wire `HIFIASM_RDNA` GFA (`rd:i` Multiplicity) + `gfatools` + `ModDotPlot` into `RDNACASSETTERESOLVER` & Terminal `N`-Gap Stripping in `CURATIONAGPGENERATOR`
1. **Pass `.p_utg.gfa` into `RDNACASSETTERESOLVER`**:
   - In [`subworkflows/local/pde_reconnaissance/main.nf:342`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-pdecuration/subworkflows/local/pde_reconnaissance/main.nf#L342), `rdna_assembly_gfa` is already emitted, but `automated_pde_curation/main.nf` and `curation_resolution/main.nf` never pass it to `RDNACASSETTERESOLVER`!
   - Pass `rdna_assembly_gfa` into `RDNACASSETTERESOLVER` and parse `hifiasm`'s native `rd:i:<int>` unitig read-depth tags (`S <utg> <seq> LN:i:<len> rd:i:<depth>`) to compute the rDNA monomer copy-number multiplier:
     $$M_B = \left\lfloor \frac{D_{\text{monomer}}}{D_{\text{anchor}}} + 0.5 \right\rfloor$$
   - Add an optional `MODDOTPLOT` process (`https://depot.galaxyproject.org/singularity/moddotplot:0.9.9--pyhdfd78af_0`, `moddotplot static -f <cassette.fa> -o <outdir>`) to generate self-identity dotplots for reconstructed rDNA cassettes.
2. **Terminal `N`-Gap Stripping & `ragtag.py agpcheck` in `ASSEMBLY_SURGERY`**:
   - In `CURATIONAGPGENERATOR`, when a scaffold is cleaved at a PDE breakpoint (`SPLIT_G_TO_S`), ensure that if the cut site borders or lands inside a pre-existing scaffolding `N`-gap (`N` or `U` component), any exposed leading or trailing `N`-gap is stripped from both resulting fragments so no chromosome or GRS contig begins or ends with `N`s.
   - Add `ragtag.py agpcheck` validation (`community.wave.seqera.io/library/ragtag:2.1.0--6125126fe9d553ab`) in `ASMFORMATQC` / `ASSEMBLY_SURGERY` alongside Sanger's `agp-tpf-utils`.

### Upgrade 6: Fix `CURATIONVALIDATION` AGP Coordinate Projection, Retired-Unloc Handling & Subtelomeric Satellite Indel Classification
Forensic audit of `nxOscSper1_DF5120.curation_validation.md` and `nxOscDolc1_PS1017.curation_validation.md` revealed three bugs in [`modules/local/curationvalidation/resources/usr/lib/curationvalidation/`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-pdecuration/modules/local/curationvalidation/resources/usr/lib/curationvalidation/):
1. **Use `agp_entries` to Project Edited Coordinates (`loaders.py:52-100`)**:
   - **Bug**: `parse_curation_audit` checked `if edit_type == "TELOMERE_CAPPING_5P":`, whereas the audit TSV emits `"TELO_5P"`, `"TELO_3P"`, `"TELO_RDNA_CAPPED_3P"`, `"ALREADY_CAPPED"`, etc. Consequently, every `TELO_5P` was evaluated at `0–0` (`I:0-0`, `0.0x` depth `WARN`), and prepending a 5' extension of $+L_{5'}\text{ bp}$ shifted all downstream coordinates on that chromosome in the curated FASTA without updating `cur_start`/`cur_end` for internal gaps or `TELO_3P`!
   - **Required Fix**: In `parse_curation_audit`, use the parsed `agp_entries` (`*.curated.agp`) to resolve the **exact curated object coordinates (`start`, `end`)** of each surgical component (`*_5p_telo_ext`, `*_3p_telo_ext`, `*_patch`, `*_rdna_cassette`, or the stitch junction between adjacent `W` components after `EXCISE_GAP`). Specifically, for a 5' extension (`TELO_5P`) of length $\Delta$, the stitch boundary to validate for spanning reads and SVs is `[Δ - 250, Δ + 250]`; for a 3' extension (`TELO_3P`), it is `[W_3p.start - 250, W_3p.start + 250]`.
2. **Treat `RETIRE_ABSORBED` Like `PURGE_*` (`verifier.py:80-111`)**:
   - **Bug**: `verifier.py:80` only checked `"PURGE" in r.action.upper()`, causing absorbed unlocalized scaffolds (`RETIRE_ABSORBED`, e.g., `I_unloc_1` in `nxOscDolc1_PS1017` and `scaffold_MT` when `bam_path` was checked after line 85) to be flagged as `WARN: Low read coverage (0.0x < 5.0x floor)`.
   - **Required Fix**: Include `"RETIRE" in r.action.upper() or "RETIRE" in r.edit_type.upper() or "ABSORBED" in r.action.upper()` in `is_purged_action` so retired/absorbed scaffolds that are properly absent from the curated FASTA/BAM receive **`PASS: Purged/retired sequence verified absent from curated assembly`**.
3. **Distinguish True Structural Assembly Errors from Subtelomeric/GRS Tandem-Satellite Indels & Compensating Alignment Shifts (`verifier.py:158-225`)**:
   - **Problem Observed**: In `nxOscSper1_DF5120`, `find_flagrant_errors` flagged `X:3,841-3,899` (`DEL -58 bp`, `VAF=0.833`, `10/12` reads) as a `HIGH` homozygous assembly defect, even though 153 bp upstream at `X:3,688` Sniffles called a compensating `INS +59 bp` (`VAF=0.692`, `9/12` reads) of the exact same 152-bp tandem satellite monomer (net shift = $1\text{ bp}$). Similarly, `I:23,902` (`INS +117 bp`, `VAF=0.800`, `16/20` reads) is a small polymorphic indel inside the $20\times$ eliminated GRS tail ($6.1\text{ kb}$ from the original contig tip).
   - **Required Fix**:
     - Detect **compensating tandem-repeat indel pairs**: if a `DEL` of length $L$ and an `INS` of length $L \pm 5\text{ bp}$ occur within $500\text{ bp}$ of each other, downgrade/filter them as `TANDEM_REPEAT_ALIGNMENT_SHIFT` rather than `HOM_SV_DEFECT`.
     - Condition VAF on realized local depth (requiring the lower bound of the 95% binomial confidence interval of VAF to exceed $0.70$, or separating low-depth $\le 25\times$ polymorphic germline GRS tails from $\ge 40\times$ somatic euchromatin), and require `sv_len >= 500` (or `INV`/`BND`) for `CRITICAL`/`HIGH` euchromatic structural defects while reporting $50\text{–}499\text{ bp}$ indels separately as `MINOR_INDEL_POLYMORPHISM`.

### Upgrade 7: Replace `IDENTIFYPUTATIVEGRS` Segmentation with the Hybrid `mosdepth`-Quantized + `DELFIES`-Bounded GRS Definition & Rename Pipeline to `nf-core/pdefinisher`

#### 7.1 Hybrid `mosdepth` + `DELFIES`-Bounded GRS Definition (`pde_identify_grs.py` & `MOSDEPTH`)
- **Problem Observed**:
  1. Keep the existing somatic core baseline (`somatic_core_baseline`) and GRS threshold ($T_{\text{GRS}}$, `t_grs`) estimator in [`pde_identify_grs.py:513-590`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-pdecuration/modules/local/identifyputativegrs/resources/usr/bin/pde_identify_grs.py#L513-L590) (`fit_em_2comp` for bimodal genomes + DELFIES-interval geometric mean calibration for unimodal genomes), which accurately identifies $T_{\text{GRS}}$.
  2. However, the downstream interval segmentation in `pde_identify_grs.py` (`run_pushdown_automaton` + all-pairs `Fallback` + `build_chromosome_segments`) had two failure modes:
     - When `mosdepth` is run with `MAPQ >= 20` (`--mapq 20 --flag 1796`) to exclude satellite multimappers, repetitive blocks *inside* the retained somatic core appear to have `0x`–low coverage purely due to `MAPQ < 20` filtering; `build_chromosome_segments` and the all-pairs fallback misclassified these internal somatic repeats as `putative_GRS`.
     - Conversely, when an internal or subtelomeric GRS had a `DELFIES` breakpoint on only **one** side (an "orphan `+`" or "orphan `-`" whose opposite end terminated at an assembly gap or uncaptured repeat), `run_pushdown_automaton` discarded the orphan breakpoint and missed the GRS.
- **Required Fix (Port the Proven Logic from `PDE_origin/scripts/cbm/generate_coverage_grs_definitions.py`)**:
  1. **`MOSDEPTH` Filtering & Quantization**:
     - Run `MOSDEPTH` with `--mapq 20 --flag 1796` (excluding secondary/supplementary/QC-fail/duplicate and `MAPQ < 20` multimappers).
     - Inside `pde_identify_grs.py`, after computing `t_grs` ($T_{\text{GRS}}$), quantize the `mosdepth` coverage windows into low-coverage drop intervals (`depth <= t_grs`) and merge adjacent drop intervals separated by $\le 1,000\text{ bp}$ (`c_merged_drops`).
  2. **Chromosomal Scaffolds — Strict `DELFIES`-Bounded Rule (Repeats Excluded from Decisions)**:
     - On chromosomal scaffolds, **no interval may be classified as `putative_GRS` unless it is bounded on at least one side by a `DELFIES` neo-telomere breakpoint**:
       - **5' Distal GRS (`[0, first_neg]`)**: If the first breakpoint on the chromosome is `-` strand and $\ge 15\%$ of `[0, first_neg]` overlaps a quantized coverage drop (or `first_neg < 1,000 bp`).
       - **Paired Intercalary GRS (`[b_s(+), next_e(-)]`)**: Evaluated between consecutive `(`+`, `-`)` breakpoints: classified as `Micro_GRS_Paired` if `inter_span < 1,000 bp` (including $0\text{–}12\text{ bp}$ reciprocal cuts!), or `Macro_GRS_Paired` if $\ge 15\%$ of `[b_s, next_e]` overlaps a quantized coverage drop.
       - **Concentric Nested Variable GRS (`+1, +2, -1, -2` on Chr X)**: If two `+` breakpoints are followed by two `-` breakpoints and $\ge 15\%$ of the outer span `[+1, -2]` overlaps a coverage drop (or exhibits stepped somatic coverage reduction), emit the outer boundary `[+1, -2]` as `putative_GRS` (sub-typed `Nested_Variable_GRS_Outer` with inner `[+2, -1]` recorded in `grs_summary.tsv`).
       - **3' Distal GRS (`[last_pos(+), chr_len]`)**: If the last breakpoint on the chromosome is `+` strand and $\ge 15\%$ of `[last_pos, chr_len]` overlaps a quantized coverage drop (or `chr_len - last_pos < 1,000 bp`).
       - **Single-Anchored Intercalary (Orphan `+` or Orphan `-` Closed by Coverage Transition)**:
         - For an unmatched internal `+` breakpoint (`b_s`), if a merged coverage-drop interval `(ds, de)` in `c_merged_drops` starts within $\le 2,000\text{ bp}$ of `b_s` (or spans `b_s`) and `de - b_s >= 1,000 bp`, emit `[b_s, de]` (`Start_DELFIES_Anchored_End_Coverage_Closed`).
         - For an unmatched internal `-` breakpoint (`b_e`, with `i > 0`), if a merged coverage-drop interval `(ds, de)` ends within $\le 2,000\text{ bp}$ of `b_e` (or spans `b_e`) and `b_e - ds >= 1,000 bp`, emit `[ds, b_e]` (`Start_Coverage_Opened_End_DELFIES_Anchored`).
       - **All other chromosomal intervals** are classified as `putative_core` (preventing low-MAPQ somatic repeats from ever being called GRS without a `DELFIES` boundary).
  3. **Preserve Base-Resolution Reciprocal Micro-GRS (`0–12 bp`) in `DELFIES`**:
     - Ensure `DELFIES` clusters `+` and `-` neo-telomere addition reads **independently by strand** (never collapsing opposing `+` and `-` additions separated by $0\text{–}12\text{ bp}$ into a single-strand BED entry) so all reciprocal $0\text{–}12\text{ bp}$ Micro-GRS cuts are preserved in `*.breakpoint_locations.bed` and `*.pde_segments.bed`.
  4. **Unplaced / Non-Chromosomal Scaffolds — Repeat-Aware Uncertainty Bounds (`[min_elim_bp, max_elim_bp]`)**:
     - Using `--red-bed` (`RED.out.bed`) strictly on non-chromosomal scaffolds:
       - Compute `non_rep_bp = length - repeat_bp`.
       - If `non_rep_bp >= 1,000 bp`: classify deterministically (`drop_bp / length >= 0.50` $\rightarrow$ `putative_GRS`; else `putative_core`).
       - If `non_rep_bp < 1,000 bp` (fully repetitive unplaced contig): classify as `putative_GRS_uncertain`, and report `canon_grs_bp`, `unplaced_grs_bp`, `unplaced_uncertain_bp`, `min_elim_bp`, `max_elim_bp`, and `est_elim_bp` in `*.grs_summary.tsv` and `cohort.curation_summary.tsv`.

#### 7.2 Rename Pipeline from `nf-core/pdecuration` to `nf-core/pdefinisher`
- Because "curation" implies manual human judgment whereas this workflow performs automated T2T and PDE finishing, update the pipeline name and top-level metadata (`nextflow.config` manifest name `nf-core/pdefinisher`, `.nf-core.yml`, `README.md`, `main.nf` workflow entry points) from `pdecuration` to `pdefinisher`.

---

## 4. Required Verification & Parity Tests

1. **Container URL Verification Test**:
   - Write a test verifying that every container URL across `modules/**/*.nf` is a real, non-fabricated image tag.
2. **Unit & Self-Tests (`--self-test`)**:
   - Every modified module script must pass `--self-test` and unit tests covering:
     - `vsearch` consensus clustering in `TELOMERICJUNCTIONPROFILER`.
     - `abpoa` multi-read overhang consensus in `TELOMEREWALKER`.
     - Dentist's $I_{\max}$ improper-alignment veto in `TARGETEDGAPASSEMBLY`.
     - `rd:i` GFA multiplicity calculation in `RDNACASSETTERESOLVER`.
     - Terminal `N`-gap stripping upon breakpoint cleavage in `CURATIONAGPGENERATOR`.
3. **Pre-Commit Hook**:
   - Keep `.git/hooks/pre-commit` active so all `Script-SHA` stamps in `modules/local/*/main.nf` stay synchronized on every commit to `v2-tool-first`.
