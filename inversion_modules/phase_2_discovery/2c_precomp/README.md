# `2c_precomp/` — SV evidence prior, window-level precomputation, landscape detector

Third block of phase 2. Takes the MDS output from `2b_mds/` plus the SV
catalogs from `MODULE_4D` (DELLY INV), `MODULE_4F` (DELLY BND), and
`MODULE_4H` (Manta), and produces everything downstream needs:

- A structured per-chromosome SV evidence prior (WHO carries what SV)
- Per-window annotations: robust z-scores, similarity matrices, seed
  eligibility, inversion-likeness, dosage heterozygosity, Q-stamps,
  SV-overlap stamps
- Seeded regions grown from MDS z-outlier seeds — consumed directly by
  `phase_4/4a/STEP_C01d` as one of its scoring-dimension inputs via
  `--cores_dir`
- A chromosome-level landscape detector (blocks, classified boundaries,
  blue-cross diagnosis, block concordance) consumed by
  `phase_4/4a/STEP_C01g` which builds `boundary_catalog_unified.tsv.gz`

Codebase: `inversion_modules` v8.5. Scripts `STEP_C00`, `STEP_C01a`,
`STEP_C01b_1`, `PHASE_01C_block_detect`.

## Layout

```
2c_precomp/
├── STEP_C00_build_sv_prior.R            # SV catalogs → per-chr sv_prior RDS
├── STEP_C01a_precompute.R               # MDS + sv_prior → per-chr annotations
├── STEP_C01b_1_seeded_regions.R         # seeded region-growing from z-outliers
├── PHASE_01C_block_detect.R             # block + boundary + blue-cross landscape
├── diags/                               # C01a / C01b diagnostic helpers
├── patches/                             # legacy patch files (terminology stale)
├── RENAMING.md                          # terminology migration tracker
└── README.md
```

## Workflow

The scripts form a small DAG. C00 and C01a are the heavy steps; C01b_1
and PHASE_01C are independent readers of C01a's precompute cache, and
all four feed phase_4 directly (no merge step in phase_2 — see
"Architectural note" below):

```
DELLY INV / BND / Manta INV / DEL VCFs
        │
        ▼
STEP_C00_build_sv_prior.R
        │
        ▼  $SV_PRIOR_DIR/sv_prior_<chr>.rds
        │
        │   ┌─── $MDS_PREFIX.mds.rds (from phase_2/2b)
        │   │
        │   │   ┌─── $DOSAGE_DIR/<chr>.dosage.tsv.gz (optional, from 2a)
        │   │   │
        ▼   ▼   ▼
STEP_C01a_precompute.R
        │
        ▼  $PRECOMP_DIR/precomp/<chr>.precomp.rds
           $PRECOMP_DIR/window_inv_likeness.tsv.gz
        │
        ├──────────────────────────────┐
        ▼                              ▼
STEP_C01b_1_seeded_regions.R    PHASE_01C_block_detect.R
(seed selection + extension)    (blocks + boundary + blue-cross)
        │                              │
        ▼                              ▼
seeded_regions_<chr>.rds        landscape/
                                  block_registry_<chr>.tsv.gz
                                  boundary_catalog_<chr>.tsv.gz
                                  blue_cross_verdicts_<chr>.tsv.gz
                                  block_concordance_<chr>.tsv.gz
                                  01C_window_pa.tsv.gz

Downstream consumption (all at phase_4 catalog-birth time):
  seeded_regions_*.rds   ──►  phase_4/4a/STEP_C01d --cores_dir
                              (D2/D5 scoring dimensions consume them)
  landscape/*.tsv.gz     ──►  phase_4/4a/STEP_C01g
                              (merged into boundary_catalog_unified.tsv.gz)
                         ──►  phase_4/4a/STEP_C01d --boundary_dir
                              (D11 boundary concordance dimension)
  sv_prior_<chr>.rds     ──►  phase_4/4b/STEP_C01i_decompose
                              (seeds k-means for genotype classes)
                         ──►  phase_4/4a/STEP_C01g (SV-breakpoint layer)
```

### Architectural note — why no `2d_seeded_merge`

An earlier design had a merge step between C01b_1 and phase_4
(`STEP_C01b_2_merge.R`, fuzzy max-min composition of membership and
geometric relations). It was retired: on test chromosomes the 1D
tail-to-head fuzzy merge overmerged across real boundaries visible in
the 2D sim_mat (LG19 ~8–15 Mb, LG28 ~7–15 Mb). The replacement design
is two parallel boundary detectors — the staircase in
`2d_candidate_detection/` and PHASE_01C here — with the seeded regions
serving as independent internal evidence consumed directly at
catalog-birth by `phase_4/4a/STEP_C01d`. C01d already accepts
`--cores_dir`, and the cores contribute to its D2 (bands) / D5
(sub-regime) scoring dimensions.

The retired merge script is at
`_archive_superseded/fuzzy_merge_abandoned/` with a README documenting
why it was dropped.

## What each script does

### C00 — `STEP_C00_build_sv_prior.R`

Parses DELLY and Manta VCFs once, packages everything needed for
sample- and locus-level SV evidence lookup. Tests 01/02/03/08 from the
test catalogue are built into this script:

- **Test 01** — Per-sample INV genotypes across DELLY + Manta (anchors)
- **Test 02** — Het-DELs within ±50 kb of INV breakpoints (markers)
- **Test 03** — Het-DELs internal to the INV (stored for decomposition)
- **Test 08** — BND paired-breakpoint triangulation (rescues INVs
  misclassified as BND)

One RDS per chromosome plus a genome-wide summary. See the script
header for the exact contents of each RDS.

### C01a — `STEP_C01a_precompute.R`

The heavy precomputation step. Runs once (~1 hour), so that C01b's
seeded-region step can be re-run in 2–5 minutes while tuning parameters.

Per-chromosome annotations computed:

- Robust z-scores (median/MAD) on each MDS axis
- Similarity matrices (derived from lostruct distance matrices)
- Seed nearest-neighbour distances (drives seed eligibility)
- Inversion-likeness score (weighted combination of het, trimodality,
  band discreteness; replaces the earlier circular PVE1 score)
- Family-likeness (diagnostic, not a gate)
- Band Q-stamps (per-window per-band dominant Q component — test_05
  precompute)
- Dosage-based het rate (second-pass upgrade when `--dosage_dir`
  provided; dosage CV replaces het_contrast in inv_likeness for
  windows that have it)
- Beta(α, β) adaptive seed threshold (test_07 — per-chromosome fit
  that produces an adaptive_seed flag as an alternative to the fixed
  inv_likeness threshold)

Last step: stamps SV prior columns (sv_inv_overlap, sv_inv_confidence,
sv_inv_af, sv_het_del_count, sv_n_anchors) onto the per-window table
by reading C00's output. These are pure annotations — they do not
alter inv_likeness or any upstream computation.

### C01b_1 — `STEP_C01b_1_seeded_regions.R`

Seeded region-growing. **Seeds come from MDS z-score outliers.** Three
parameter scale tiers (1S small / 1M medium / 1L large) run
independently per chromosome to capture different-width candidates.

The algorithm:

1. **Seed selection.** A window becomes a seed candidate iff its
   `|max_abs_z|` exceeds `seed_z_min` OR its `inv_likeness` is very
   high (≥ 0.90) OR it passes the test_07 Beta adaptive gate. AND its
   `seed_nn_dist` is finite and small. Z-score is the primary
   criterion; inv_likeness and the Beta gate are OR-gates for cases
   where z alone misses a signal (e.g. compact inversions with modest
   per-axis z).
2. **Seed ordering.** Highest priority first:
   `priority = max_abs_z + 5 × inv_likeness + 2 × spiky_inv_score`
3. **Extension.** From each seed, extend bidirectionally accepting
   adjacent windows whose continuity score exceeds a per-tier accept
   threshold. Windows below the tolerate threshold accumulate damage
   (a scalar budget). Extension halts when damage exceeds `DMG_MAX`
   or when a sharp drop (> 2×sd below rolling mean) is detected.

Optional test_26 kin-pruned retention: after all regions are found,
test each one's sim_mat block contrast on a kin-pruned sample subset.
Regions whose signal collapses under pruning are flagged, not dropped
(collapse under pruning is informative in a hatchery-founder context).

### PHASE_01C — `PHASE_01C_block_detect.R`

A **landscape detector**, independent of the seed-based track. Reads
the same `precomp/<chr>.precomp.rds` produced by C01a and answers three
related questions at the chromosome scale:

1. **Where are the blocks?** Elevated-similarity regions on the
   sim_mat, detected via row-profile clustering — not from seeds,
   not from staircase votes. Written to `block_registry_<chr>.tsv.gz`.
2. **What kind of boundary sits at each transition?** Classified
   into `clear_hard` / `soft` / `inner_hard` / `inner_soft` /
   `diffuse` based on the shape of the local sim decay. Written to
   `boundary_catalog_<chr>.tsv.gz`.
3. **Are suspected boundaries real or assembly artefacts?**
   "Blue-cross" diagnosis: does the drop occur at an assembly gap or
   AGP-scaffold edge? Written to `blue_cross_verdicts_<chr>.tsv.gz`.

It also produces `block_concordance_<chr>.tsv.gz` (which blocks share
the same sample-partitioning pattern) and `01C_window_pa.tsv.gz`
(per-window block-membership in a P/A matrix).

Mode-dependent thresholds — `hatchery` mode uses slightly lower
`BLOCK_THRESH_ABOVE` and `MIN_BLOCK_SIM` than `wild`, reflecting the
weaker genome-wide background contrast in hybrid hatchery populations.

PHASE_01C's outputs are consumed at phase_4 catalog-birth time by:

- `phase_4/4a/STEP_C01g_boundary_catalog` — one of five boundary
  sources (PHASE_01C, staircase from 2d, seeded regions from C01b_1,
  SV breakpoints from C00, blue-cross inner boundaries from PHASE_01C
  itself) that get merged by proximity into
  `boundary_catalog_unified.tsv.gz`
- `phase_4/4a/STEP_C01d_candidate_scoring` — the unified boundary
  catalog feeds C01d's D11 dimension (boundary concordance), and the
  block_concordance layer can inform sub-regime scoring

Usage:

```bash
Rscript PHASE_01C_block_detect.R <precomp_dir> <landscape_outdir> \
  [--mode hatchery|wild] [--gaps gaps.bed] [--agp scaffold.agp] \
  [--chrom C_gar_LG28] [--local_range 80]
```

The `--gaps` / `--agp` files are optional but enable higher-quality
blue-cross diagnosis — without them, a "blue cross" (abrupt local
similarity drop) cannot be distinguished from a small assembly gap.

## Two scripts, same VCFs, different questions

This needs to be stated because it is confusing and it is not a bug.

Both `STEP_C00` (in this folder) and `MODULE_5A2_breakpoint_validation`
(a separate module at `phase_3_refine/`, corresponding to the repo's
`MODULE_5A2_breakpoint_validation/`) read the same DELLY and Manta
VCFs. They answer different questions:

|   | C00 — `build_sv_prior` | MODULE_5A2 — breakpoint validation |
|---|---|---|
| Reads | DELLY / Manta INV / DEL / BND VCFs | DELLY / Manta INV / BND VCFs + BAMs |
| Question | **WHO** carries what arrangement? | **WHERE** are the breakpoints at bp resolution? |
| Output | per-sample genotypes, anchor assignments, het-DEL carriers | refined boundary coordinates, concordance scores, breakpoint confidence |
| Consumed by | C01a precompute (stamps per-window), phase 4b decompose (seeds k-means) | phase 4a C01g boundary catalog, `classify_inversions` |
| Granularity | sample-level | locus-level |

They parse the same VCFs but extract different slices. C00 is about
samples, breakpoint_validation is about coordinates. No duplicated
analysis work — just two independent readers serving different
phases.

**Future cleanup (not blocking the manuscript):** a shared
`parse_sv_catalogs.R` helper that both C00 and breakpoint_validation
source. Not doing this now because the two outputs serve very
different downstream consumers and each parser already works.

## The 4-layer independence framework

The pipeline detects inversions from four mutually-independent data
sources. The SV prior built here is Layer B.

| Layer | Source | Algorithm | Produced by |
|---|---|---|---|
| A | dosage (genotype likelihoods) | lostruct local PCA → MDS | phase_2/2a + 2b |
| **B** | **SV caller output** | **DELLY2 + Manta catalogs** | **STEP_C00 (this folder)** |
| C | Clair3 phased genotypes | GHSL haplotype contrast | phase_2/2e_ghsl |
| D | genotype–breakpoint association | Fisher odds ratio linking A, B, C | phase_4 |

Layers are designed to fail independently. A candidate that converges
across layers is high-confidence (Tier 1); a candidate visible only to
one layer is low-confidence or a layer-specific artefact (Tier 3–4 or
SV-only).

## SV prior RDS structure

```r
sv_prior <- list(
  chrom                 = character,
  n_inv_calls           = integer,
  n_het_dels_breakpoint = integer,
  n_het_dels_internal   = integer,

  inv_calls             = data.table,   # test_01: INV calls
  gt_matrix             = matrix,       # sample × inv_id genotypes
  breakpoint_dels       = data.table,   # test_02: het-DELs at breakpoints
  test02_verification   = data.table,   # test_02 concordance check
  internal_dels         = data.table,   # test_03: internal het-DELs
  bnd_triangulated      = data.table,   # test_08: BND-triangulated inversions
  sample_inv_states     = data.table,   # combined per-sample × per-inv state

  params                = list,         # parameters used to build
  built                 = POSIXct       # build timestamp
)
```

Read via `readRDS(paste0(SV_PRIOR_DIR, "/sv_prior_", chr, ".rds"))`.

## Terminology note

Several names in this codebase are historical nicknames that have been
replaced with scientific vocabulary. The full map is in `RENAMING.md`.
Key replacements relevant here:

- `flashlight` → `sv_prior`
- `snake` / `snake1` → seeded region-growing (seeds from **MDS z-score
  outliers**)
- `core` / `cores` → seeded regions (algorithm output)
- `core family` → scale tier (1S / 1M / 1L parameter sets)
- `cheat N` → `test_NN` (biological evidence tests, 01 through 26)

Grep recipes to audit other scripts are in `RENAMING.md` section 8.

**Caveat for the `test_05` rename specifically (2026-04-17 cleanup):**
C01a writes columns named `test05_fst_pc1`, `test05_fst_q_best`,
`test05_fst_q_best_k`, `test05_family_fst_ratio`. The downstream reader
in `phase_4/4a/STEP_C01d` still names its output column
`cheat5_family_fst_ratio` (to avoid forcing a rename across 4b/4c/4e
registry keys and test suites), but accepts either name on read. Don't
rename the C01d output until a follow-up sweep aligns it across the
downstream consumers (`compute_candidate_status.R`,
`test_registry_sanity.py`, manuscript tables).

## CLI flag notes (2026-04-17 cleanup)

Several flags in `phase_4/4a/` scripts were parsed but never read. On
2026-04-17 they were silenced — still accepted for back-compat with
existing launchers, but no longer have any effect:

- `STEP_C01d_candidate_scoring.R --flashlight_dir` — D7's SV info
  already comes from the staircase scoring table's `sv_overlap_pct` /
  `n_sv_hits` columns (populated by `phase_2/2d/STEP_D06_sv_overlap.R`).
  The sv_prior RDS built by `STEP_C00` here is consumed by C01a (which
  stamps SV columns onto windows), not by C01d directly.
- `STEP_C01e_candidate_figures.R --repeats` and `--het_dir` — never
  referenced in the script body.
- `STEP_C01g_boundary_catalog.R --ref_fasta` — never referenced.
- `STEP_C01g_boundary_catalog.R --scores` — archived with the Cheat 17
  fossil-detection block (see
  `_archive_superseded/cheat17_fossil_detection/README.md`). Fossil
  detection required a two-pass run of C01g → C01d → C01g that was
  never wired; removing the flag simplifies the DAG.

## Paths and config

Every path comes from `00_inversion_config.sh`. The relevant variables
for this folder:

```
$DOSAGE_DIR                      # per-chr dosage (optional input to C01a)
$MDS_DIR, $MDS_PREFIX            # MDS output (input to C01a)
$SV_PRIOR_DIR                    # C00 output, C01a input
$PRECOMP_DIR                     # C01a output, C01b input
$SIM_MATS_DIR                    # similarity matrices (C01a output)
```

See the "config architecture" section of the repo-level
`CONFIG_ARCHITECTURE.md` for how this folder's config fits into the
multi-module hierarchy.

## Troubleshooting

**C01a fails with "Missing window_id in <chr>"** — you jumped over
STEP_A03 Stage 2 (the master registry merge in phase_2/2a). Run
`LAUNCH_A03_dense_registry_stage2.slurm` first.

**C01a runs but sv_prior columns are all NA** — `--sv_prior_dir` was
not passed, or C00 has not produced output for the chromosome yet.
Check `ls $SV_PRIOR_DIR/sv_prior_*.rds`.

**C01b produces zero seeded regions for a chromosome** — check the
`[seeded_regions]` log: the `z-score breakdown` and `inv-likeness
breakdown` lines show how many windows pass each seed criterion. If
`z ≥ 2` has very few windows, either the chromosome genuinely has
low regime structure or the MDS mode is wrong (chunked vs
chromosome — see phase_2/2b README).

**C01b fails with "Missing precomp file"** — check
`ls $PRECOMP_DIR/precomp/<chr>.precomp.rds`. If absent, re-run C01a
for that chromosome.

**test_26 reports "insufficient_pruned" for most regions** — your
`--pruned_samples` file does not overlap the PCA sample list, or the
file's sample names are not in the same format as the PC column names
in the precomp RDS. Check by comparing one sample name from the
pruned list against `colnames(obj$pca)` from any per-chr PCA RDS.
