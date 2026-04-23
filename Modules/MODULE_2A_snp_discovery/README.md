# MODULE_2A — SNP Discovery & Structure Foundation

Callable-site masking, ANGSD biSNP discovery, distance-thinned panels, BEAGLE GL generation, PCAngsd, NGSadmix, evalAdmix, ngsRelate relatedness pruning, and canonical best-seed-by-K selection. Produces the SNP panels, BEAGLE files, ancestry assignments, and pruned sample lists consumed by every downstream module.

## Why this module exists (for the inversion paper)

Every ancestry-flavored inversion signal is anchored to this module's biallelic SNP panel. ANGSD calls biSNPs from genotype likelihoods at cohort scale — the only statistically appropriate way to work with low-coverage data without hard genotype calls. The thinned panel produced here feeds NGSadmix for global ancestry Q, and Engine B / instant_q for per-window local Q. The downstream inversion-discovery pipeline (`inversion_modules/phase_2_discovery/`) computes against those Q estimates through its 4-layer evidence framework (Layer A dosage / local PCA, Layer C GHSL haplotype contrast, integrative scoring) — without this module's panel, there is nothing to compute against.

NAToRA first-degree pruning (A07) also produces the 81-unrelated subset used for Hardy-Weinberg-sensitive tests including per-group Hobs confirmation (phase_qc_shelf STEP_Q07b + Q07c, previously MODULE_5E) and inversion genotype-frequency tests.

Note: ANGSD biSNPs are **not** the same catalog as Clair3 SNPs (MODULE_4A). ANGSD works in GL space for population-level Q/Fst; Clair3 produces per-sample hard genotypes for per-sample analyses. Both are needed, for different downstream consumers.

## Pipeline

```
Filtered BAMs from MODULE_1 (226 samples, ~9× mean)
  │
  ├─ A. Callable mask & SFS ──────────────────────────────────────────
  │   A01  mask_regions_from_fasta.py → callable/non-callable BED
  │        RF/chromosome chunk lists for SLURM arrays
  │        per-chunk SAF (ANGSD -doSaf 1, folded)
  │        merge SAF → global folded SFS → pest prior
  │
  ├─ B. biSNP discovery ──────────────────────────────────────────────
  │   A02  ANGSD biSNP per chunk (-doMaf, -SNP_pval 1e-6, -pest)
  │
  ├─ C. Panel construction ───────────────────────────────────────────
  │   A03  merge per-chunk mafs → distance-thin (200/500/1k/5k/10k/25k bp)
  │        ANGSD sites index each panel
  │
  ├─ D. BEAGLE generation ────────────────────────────────────────────
  │   A04  per-RF BEAGLE GL (-doGlf 2, -doMajorMinor 3)
  │   A04  whole-genome BEAGLE GL
  │   A05  merge per-RF → whole-genome BEAGLEs
  │
  ├─ E. Structure (all 226 samples) ──────────────────────────────────
  │   A06  PCAngsd (eigen=6, per panel)
  │        NGSadmix K=2..12 × 3 seeds per panel
  │        evalAdmix per K × seed
  │        best-seed-by-K selection (loglik → mean|resid| → seed)
  │
  ├─ F. Relatedness & pruning ────────────────────────────────────────
  │   A07  ngsRelate (thin-500 panel)
  │        NAToRA multi-cutoff classification
  │        greedy first-degree pruning → pruned_samples.txt
  │        3-panel relatedness figure
  │
  ├─ G. Structure (pruned unrelated set) ─────────────────────────────
  │   A06  same as E but on pruned sample list
  │
  └─ H. Merge summaries ─────────────────────────────────────────────
      A08  combined_best_seed_by_K.tsv (sample_set = all | pruned)
  │
  └─→ SNP panels + BEAGLEs + ancestry Q tables + pruned sample list
        consumed by MODULE_2B, MODULE_3, MODULE_4A, MODULE_5A, MODULE_6
```

## Directory layout

```
MODULE_2A_SNP_Discovery/
  00_module2a_config.sh                ← central config (source this everywhere)
  README.md
  docs/
    MODULE_2A_methods.md               ← manuscript-ready methods prose
  steps/
    STEP_A01_masks_sfs.sh              ← callable mask → SAF → SFS → pest
    STEP_A02_bisnps.sh                 ← biSNP discovery per chunk
    STEP_A03_panels.sh                 ← merge + distance-thin + ANGSD-index
    STEP_A04_beagles.sh                ← BEAGLE GL per-RF + WG
    STEP_A05_merge_beagles.sh          ← merge per-RF BEAGLEs → WG
    STEP_A06_structure.sh              ← PCAngsd + NGSadmix + evalAdmix + best-seed
    STEP_A07_relatedness.sh            ← ngsRelate + NAToRA + pruning + plot
    STEP_A08_merge_structure_summaries.sh ← combined best-seed table
  slurm/
    SLURM_A01a_saf_chunk.sh            ← SLURM array: one SAF per chunk
    SLURM_A01b_merge_sfs.sh            ← SLURM: merge SAF + fold SFS
    SLURM_A02_bisnp_chunk.sh           ← SLURM array: one biSNP per chunk
    SLURM_A04a_beagle_rf.sh            ← SLURM array: one BEAGLE per RF
    SLURM_A04b_beagle_wg.sh            ← SLURM: whole-genome BEAGLE
    SLURM_A06a_pcangsd.sh              ← SLURM: PCAngsd per panel
    SLURM_A06b_ngsadmix.sh             ← SLURM array: NGSadmix per K×seed
    SLURM_A06c_ngsadmix_global.sh      ← SLURM: NGSadmix global (full node)
  launchers/
    LAUNCH_module2a.sh                 ← main wrapper — the ONLY script you call
  utils/
    mask_regions_from_fasta.py         ← BED extraction from reference FASTA
    select_best_seed_by_K.R            ← best-seed selection + palette + exports
    prune_first_degree_pairs.py        ← greedy first-degree pruning
    NAToRA_Public.py                   ← NAToRA multi-cutoff relatedness
    plot_relatedness_3panel.R          ← 3-panel relatedness figure
```

## Naming convention

Matches inversion codebase v8.5:

- Config: `00_module2a_config.sh` (like `00_inversion_config.sh`)
- Steps: `STEP_{PHASE}{NN}_{description}.sh` — Phase A = compute pipeline
- SLURM: `SLURM_{PHASE}{NN}{sub}_{description}.sh` — submitted by steps
- Launchers: `LAUNCH_module2a.sh`
- Utils: standalone tools (R, Python) called by steps

## Usage

```bash
# Validate config, tools, inputs
./launchers/LAUNCH_module2a.sh init

# Run individual steps
./launchers/LAUNCH_module2a.sh masks_sfs
./launchers/LAUNCH_module2a.sh call_bisnps
./launchers/LAUNCH_module2a.sh build_panels
./launchers/LAUNCH_module2a.sh make_beagles
./launchers/LAUNCH_module2a.sh structure_all
./launchers/LAUNCH_module2a.sh relatedness
./launchers/LAUNCH_module2a.sh structure_pruned
./launchers/LAUNCH_module2a.sh merge_summaries

# Run full pipeline
./launchers/LAUNCH_module2a.sh all
```

## CLI options

```
--config <file>     Config file. Default: ../00_module2a_config.sh
--task <name>       Alternative to positional subcommand
--dry-run           Print commands without executing
--force             Rebuild outputs even if existing
--resume            Skip completed outputs (default)
--thin <W>          Override default thinning for ancestry/relatedness
--kmin <int>        Override K_MIN. Default: 2
--kmax <int>        Override K_MAX. Default: 12
--seeds <list>      Override seeds (comma-separated). Default: 1,2,3
--threads <int>     Override thread count
--scope <name>      Analysis scope: global, chromosome, or thin label
```

## Canonical outputs (per sample set)

| File | Description | Used by |
|------|-------------|---------|
| `all_seed_metrics_by_K.tsv` | Every K × seed: loglik, evalAdmix residuals, file paths | Internal QC |
| `best_seed_by_K.tsv` | One row per K: selected best seed + metrics | MODULE_2B, MODULE_5A |
| `best_seed_copied_files.tsv` | Provenance: which files were copied as `*_best` | Internal |
| `cluster_palette_by_K.tsv` | Stable 12-color palette per K × cluster | All figure modules |
| `sample_order_reference.tsv` | Canonical sample ordering | All downstream |
| `sample_main_ancestry_by_K.tsv` | Per-sample: dominant cluster, full Q, color per K | MODULE_5A, MODULE_6 |
| `combined_best_seed_by_K.tsv` | Merged all + pruned with sample_set column | Manuscript |

## Key design decisions

**Folded SFS throughout.** No suitable outgroup species; hatchery broodstock with limited founders violates HWE assumptions required for unfolded SFS. All SFS computations use `-fold 1`.

**Pest prior for SNP discovery.** Per-site allele frequency prior estimated from folded SFS via realSFS, applied to biSNP calling with `-pest`. This avoids the common ANGSD pitfall of flat priors biasing low-MAF calls.

**Structure runs twice.** `STEP_A06_structure.sh` is ONE script that takes `--samples <file>`. The launcher calls it twice — once with all 226 samples, once with pruned samples after relatedness filtering. Same logic, same outputs, different sample list.

**Best-seed selection rule.** Within each K: (1) highest log-likelihood, (2) lowest mean |evalAdmix off-diagonal residual|, (3) lowest max |evalAdmix residual|, (4) lowest seed number as tie-break.

## ANGSD shared filter parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| -GL | 1 (SAMtools) | Standard for low-coverage WGS |
| -minQ | 25 | Base quality |
| -minMapQ | 25 | Mapping quality |
| -baq | 1 | BAQ recalibration |
| -C | 50 | Excessive mismatch downweight |
| -setMinDepthInd | 3 | Min depth per individual (from MODULE_1 QC) |
| -setMaxDepthInd | 57 | Max depth per individual (from MODULE_1 QC) |
| -minInd | 200 | Min individuals with data (200/226 = 88%) |
| -SNP_pval | 1e-6 | SNP significance threshold |
| -minMaf | 0.05 | Minimum minor allele frequency |

## Thinning distances

| Distance | Panel type | Primary use |
|----------|-----------|-------------|
| 200 bp | Fine | BEAGLE for local LD analyses |
| 500 bp | Fine | Default for structure + relatedness |
| 1,000 bp | Fine | Conservative LD-pruned structure |
| 5,000 bp | Broad | Genome-wide PCA |
| 10,000 bp | Broad | Low-SNP structure QC |
| 25,000 bp | Broad | Ultra-sparse sanity check |

## Dependencies

ANGSD (with realSFS, sites index), PCAngsd, NGSadmix, evalAdmix, ngsRelate, samtools, python3, R (data.table, ggplot2)

## Sidecar convention

Every step writes two audit files alongside its outputs:
- **`{step}.arg`** — all parameters, tool versions, paths, exact commands (written at start)
- **`{step}.results`** — all output file paths + descriptions (written at end)

## Methods

See [`docs/MODULE_2A_methods.md`](docs/MODULE_2A_methods.md) for manuscript-ready methods prose.
