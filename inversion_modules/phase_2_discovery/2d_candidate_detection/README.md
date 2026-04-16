# `2d_candidate_detection/` — matrix-based candidate detection

Fourth block of phase 2 discovery. Takes the similarity matrices from
`2c_precomp/` and produces candidate blocks via a boundary-voting
staircase detector, a suite of matrix transforms, and a consensus layer
that cross-validates blocks across transforms and NN-smoothing scales.

Two detection tracks run in parallel in phase 2:

| Track | Folder | Seed / signal |
|---|---|---|
| Seed-based (MDS z-outlier extension) | `2c_precomp/STEP_C01b_1_seeded_regions.R` | per-window `max_abs_z` and `inv_likeness` |
| Matrix-based (this folder) | `2d_candidate_detection/run_all.R` | sim_mat geometry (block edges in similarity space) |

Both feed the same downstream scoring (`STEP_C01d_candidate_scoring.R`).
The bridge script (`STEP_D12_bridge_to_C01d.R`) converts this folder's
output to the legacy `triangle_intervals.tsv.gz` format that C01d already
reads, so integrating this track does not require touching C01d or
anything after it.

**This folder replaces only the old `STEP_C01c_triangle_regimes.R`.**
Everything from C01d onwards is unchanged.

## What this folder is NOT

- Not a replacement for the seed-based track. The two tracks answer the
  same "where are candidate blocks" question from independent signals.
  A candidate that appears in both is higher confidence than one that
  appears only in one.
- Not a replacement for precompute. It consumes `$SIM_MATS_DIR` (built
  by `STEP_C01a_precompute.R`) and optionally the per-chromosome
  `precomp/<chr>.precomp.rds` for the peeling diagnostic.
- Not a geometric thing. The "triangles" in the old naming were a
  nickname for elevated square blocks on the similarity matrix as they
  look on a heatmap; they are not triangle-shaped in any mathematical
  sense. The new code uses **block** everywhere. The only genuine
  triangles in the code are the geometric upper/lower triangles of a
  symmetric matrix in `STEP_D17_plot_marginal_tracks.R`.

## Layout

```
2d_candidate_detection/
├── 00_config.R                              # module-specific thresholds
├── run_all.R                                # orchestrator (9 phases)
│
├── STEP_D01_staircase_boundaries.R          # boundary vote + block registry
├── STEP_D02_nn_persistence.R                # NN-scale survival
├── STEP_D03_flatness.R                      # decay profile
├── STEP_D04_interior_cv.R                   # interior homogeneity
├── STEP_D05_ghsl_stability.R                # sample partition stability
├── STEP_D06_sv_overlap.R                    # SV breakpoint matching
├── STEP_D07_matrix_transforms.R             # 6 image-processing transforms
├── STEP_D08_block_scoring.R                 # per-block metrics + shape class
├── STEP_D08b_local_pca.R                    # local PCA / ICA
├── STEP_D09_nn_sweep_tree.R                 # interval tree across NN scales
├── STEP_D09n_peeling_diagnostic.R           # sample-zeroing diagnostic
├── STEP_D10_variant_consensus.R             # cross-variant block consensus
├── STEP_D11a_ngsrelate_perchr.sh            # per-chr kinship (optional, one-shot)
├── STEP_D11b_perchr_pruning_table.R         # per-chr pruning table
├── STEP_D12_bridge_to_C01d.R                # format adapter → C01d
├── STEP_D13_plot_annotated_simmat.R         # chromosome-wide heatmap plot
├── STEP_D14_landscape_classifier.R          # block classification + confidence
├── STEP_D15_plot_zoomed_regions.R           # saturation-driven zoomed plots
├── STEP_D17_plot_marginal_tracks.R          # marginal-track chromosome plot
│
├── LAUNCH_run_all.slurm                     # 28-chr array, all phases
├── LAUNCH_test_LG01.slurm                   # quick smoke test
├── LAUNCH_peel.slurm                        # Phase 9 only (peeling)
│
├── tests/
│   └── test_staircase.R                     # synthetic-data tests
│
└── _archive_old_names/                      # original v9.3 docs kept for reference
    ├── HOW_TO_RUN_v9.3.md
    ├── INTEGRATION_AUDIT_v9.3.md
    └── AUDIT_LABELS_DECISIONS_v9.3.md       # ← READ THIS for full label semantics
```

## Workflow

```
$SIM_MATS_DIR/<chr>.sim_mat_nn{0,20,40,80}.rds     (from 2c_precomp)
$PRECOMP_DIR/<chr>.precomp.rds                     (optional, for peeling)
           │
           ▼
run_all.R  (all 9 phases, per chromosome — SLURM array)
           │
           ▼
$D_OUTDIR/
   blocks_<chr>.tsv                   Phase 1 — block registry
   artifacts_<chr>.tsv                Phase 1 — artifact registry
   staircase_<chr>.rds                Phase 1 — profiles + votes
   variants_<chr>.rds                 Phase 2 — 6 matrix transforms
   block_scores_<chr>.tsv             Phase 3 — per-variant scores
   nn_persistence_<chr>.tsv           Phase 4 — NN-scale survival
   nn_tree_<chr>.tsv                  Phase 5 — interval tree
   consensus_<chr>.tsv                Phase 6 — cross-variant consensus
   ev_flatness_<chr>.tsv              Phase 7a
   ev_cv_<chr>.tsv                    Phase 7b
   ev_ghsl_<chr>.tsv                  Phase 7c
   ev_sv_<chr>.tsv                    Phase 7d
   ev_pca_<chr>.tsv                   Phase 7e
   scoring_table_<chr>.tsv            Phase 8 — everything joined
   peel_diagnostic_<chr>.tsv          Phase 9 — peel before/after
   peel_groups_<chr>.tsv              Phase 9 — which samples zeroed
           │
           ▼
STEP_D12_bridge_to_C01d.R
           │
           ▼
triangle_intervals.tsv.gz            ← legacy filename for C01d compatibility
(plus empty sample_composition / bridges / offdiag_linkage tables)
           │
           ▼
STEP_C01d_candidate_scoring.R        (unchanged from here on)
```

## The two detection tracks

| Criterion | Seed-based (2c) | Matrix-based (2d, this folder) |
|---|---|---|
| Input | per-window annotations (z-scores, inv_likeness, sim_mat) | sim_mat at multiple NN smoothing scales |
| Seed / starting point | MDS z-outlier windows | step-down votes across row profiles |
| Block definition | region grown from seed under damage-budget extension | interval between high-vote boundaries |
| Multi-scale | three parameter tiers (1S/1M/1L) | NN-scale sweep (0, 20, 40, 80) + matrix variants |
| Output | seeded_regions_<chr>.rds | blocks_<chr>.tsv + scoring_table_<chr>.tsv |
| Merges | no (previous fuzzy merge retired — see `_archive_superseded/fuzzy_merge_abandoned/`); seeded regions now consumed directly by phase_4/4a/C01d via `--cores_dir` | no merge — consensus across variants instead |

They are designed to fail independently. A candidate visible to both is
higher confidence.

## Inputs and paths

Every path comes from `$00_inversion_config.sh`. The variables this
folder reads:

```
$SIM_MATS_DIR           # from config — per-chr sim_mat_nn*.rds
$PRECOMP_DIR            # from config — per-chr precomp.rds (for Phase 9)
$INVDIR                 # for output location
$DELLY_INV_VCF          # for Phase 7d SV overlap (optional)
```

This folder writes under:
`$INVDIR/06_mds_candidates/snake_regions_multiscale/candidate_detection_out/`

(The `snake_regions_multiscale/` path component is legacy — see
`RENAMING.md` section 5.1 for the planned directory rename.)

## Minimum viable run

```bash
# 1. Quick smoke test (one chr, phases 1-3, no transforms)
sbatch LAUNCH_test_LG01.slurm

# 2. Full run, 28 chromosomes, all 8 core phases
sbatch LAUNCH_run_all.slurm

# 3. Optional: Phase 9 peeling (after 1-8 completed)
#    Build per-chr kin table first
sbatch STEP_D11a_ngsrelate_perchr.sh
Rscript STEP_D11b_perchr_pruning_table.R \
  --ngsrelate_dir ngsrelate_perchr/ \
  --outfile $D_OUTDIR/per_chr_pruned.tsv

sbatch LAUNCH_peel.slurm

# 4. Bridge to C01d (single chromosome at a time, or in a loop)
for chr in C_gar_LG0{1..9} C_gar_LG{10..28}; do
  Rscript STEP_D12_bridge_to_C01d.R \
    --detector_dir $D_OUTDIR \
    --outdir       ${D_OUTDIR}/triangle_from_detector
done
```

## Key concepts

### Phase 1 — staircase boundary voting

For each row of the sim_mat, walk left and right and record the
positions where similarity drops by at least `STAIR_MIN_DROP` and stays
dropped for `STAIR_PERSIST_N` consecutive bins. A boundary is a
position where many rows agree. Blocks are the intervals between
high-vote boundaries with elevated interior similarity.

### Phase 2 — matrix transforms

Six image-processing transforms of the raw sim_mat. Each exposes
different block structure:

| Transform | What it does | Why |
|---|---|---|
| `distcorr` | subtract median sim at each off-diagonal distance | removes baseline decay, blocks pop out |
| `localnorm` | local z-score (median/MAD in ±50 bin window) | highlights local contrast, flattens global scale |
| `denoised` | iterated 2D median filter (fused-lasso proxy) | cleans speckle, preserves edges |
| `resid_bg` | subtract SVD rank-3 reconstruction | removes large-scale background |
| `edge` | Sobel gradient magnitude | boundaries light up |
| `support` | adaptive threshold at 0.75 quantile | binarises strong signal |

Blocks are re-detected on each variant; the consensus layer (Phase 6)
counts how many variants see each block. HIGH = ≥3 variants,
MODERATE = 2, LOW = 1.

### Phase 4 — NN-scale persistence

Runs the staircase independently on each of sim_mat_nn{0,20,40,80}.
A real inversion survives moderate NN smoothing; family LD does not.
`survives_nn40`, `survives_nn80` and `nn_topology` (stable / disappears /
splits / complex) are the key columns.

### Phase 9 — peeling diagnostic

For each block, recompute sim_mat after zeroing out sample subsets and
compare block contrast before / after:

- **L1**: zero ngsRelate-pruned samples (genome-wide or per-chr)
- **L1b**: zero chr-local PC1 relatives (|cor| > 0.7)
- **L2**: cluster samples by PC1 trajectory within the block, zero
  dominant cluster
- **L3**: zero high-leverage samples and C01i HOM_INV carriers

Effect class: `stable` / `weakened` / `disappeared` / `revealed_child` /
`ambiguous`. Most informative column: `l1b_effect`. "disappeared" after
L1b suggests the block is driven by a few closely-related samples
(family LD). "stable" suggests genuine structural signal.

## For full label semantics

Read `_archive_old_names/AUDIT_LABELS_DECISIONS_v9.3.md` — this is the
authoritative document describing what every column means, every
threshold, every classifier. It has not been rewritten because the
content is correct; only the filenames it refers to have changed. A
mapping from old to new names is at the bottom of this README.

## Terminology map vs the v9.3 drop

The content of this folder is functionally identical to the uploaded
`inv_detect_v9.3/` drop. Only names have changed.

| Old (v9.3 drop) | New (this folder) |
|---|---|
| `inv_detect_v9.3/` | `phase_2_discovery/2d_candidate_detection/` |
| `00_run_all.R` | `run_all.R` |
| `01_staircase_detector.R` | `STEP_D01_staircase_boundaries.R` |
| `02_evidence_nn_persistence.R` | `STEP_D02_nn_persistence.R` |
| `03_evidence_flatness.R` | `STEP_D03_flatness.R` |
| `04_evidence_interior_cv.R` | `STEP_D04_interior_cv.R` |
| `05_evidence_ghsl.R` | `STEP_D05_ghsl_stability.R` |
| `06_evidence_sv_overlap.R` | `STEP_D06_sv_overlap.R` |
| `07_matrix_cheats.R` | `STEP_D07_matrix_transforms.R` |
| `08_bloc_scoring.R` | `STEP_D08_block_scoring.R` |
| `08_evidence_local_pca.R` | `STEP_D08b_local_pca.R` |
| `09_nn_sweep_tree.R` | `STEP_D09_nn_sweep_tree.R` |
| `10_consensus.R` | `STEP_D10_variant_consensus.R` |
| `11a_run_ngsrelate_perchr.sh` | `STEP_D11a_ngsrelate_perchr.sh` |
| `11b_build_perchr_pruning_table.R` | `STEP_D11b_perchr_pruning_table.R` |
| `12_bridge_to_codebase.R` | `STEP_D12_bridge_to_C01d.R` |
| `13_plot_annotated_simmat.R` | `STEP_D13_plot_annotated_simmat.R` |
| `14_landscape_classifier.R` | `STEP_D14_landscape_classifier.R` |
| `15_plot_zoomed_systems.R` | `STEP_D15_plot_zoomed_regions.R` |
| `16_test_staircase.R` | `tests/test_staircase.R` |
| `17_plot_marginal_tracks.R` | `STEP_D17_plot_marginal_tracks.R` |
| `STEP_C01n_blockwise_peeling_diagnostic.R` | `STEP_D09n_peeling_diagnostic.R` |
| `cheat_*` functions | `mtx_*` functions |
| `CFG$CHEAT_*` config vars | `CFG$MTX_*` |
| `CFG$BLOC_*` config vars | `CFG$BLOCK_*` |
| `--skip-cheats` CLI flag | `--skip-transforms` |
| `bloc_scores_<chr>.tsv` | `block_scores_<chr>.tsv` |

Three things intentionally preserved:

1. **Output filename `triangle_intervals.tsv.gz`** emitted by
   `STEP_D12_bridge_to_C01d.R`. That filename is what C01d reads; it's
   a legacy interface, flagged as a future-cleanup item.
2. **Geometric "triangle"** in `STEP_D17_plot_marginal_tracks.R` where
   it refers to upper/lower triangles of a symmetric matrix. That
   usage is correct math.
3. **`STAIR_*`, `NN_*`, `CONSENSUS_*`, `LANDSCAPE_*`, `PEEL_*` config
   prefixes** in `00_config.R`. These are descriptive of what the
   parameters control, not historical nicknames.

## Troubleshooting

**No sim_mat files found.**
`run_all.R` looks for `<chr>.sim_mat_nn<N>.rds` under `$SIM_MATS_DIR`,
with fallback to `$PRECOMP_DIR/<chr>.precomp.rds` for NN=0. Check
`ls $SIM_MATS_DIR/` — you should see per-chr files with `nn0`, `nn20`,
`nn40`, `nn80` suffixes, produced by `STEP_C01a_precompute.R`.

**Phase 2 (transforms) is very slow.**
Expected on 9000×9000 matrices. Local contrast and denoise are
O(N²·window²). Use `--skip-transforms TRUE` to get the staircase
output fast, then re-run with transforms overnight.

**Phase 9 skipped with "no precomp RDS found".**
Pass `--precomp-dir $PRECOMP_DIR` explicitly. The default
(`$SIM_MATS_DIR`) often points one level too deep.

**C01d doesn't see any candidates from the bridge output.**
Check that `triangle_intervals.tsv.gz` is non-empty under
`$D_OUTDIR/triangle_from_detector/`. Empty means Phase 1 found no
blocks — either thresholds too strict (increase `STAIR_MIN_BLOCK_WIDTH`
or decrease `STAIR_MIN_DROP` in `00_config.R`) or the chromosome
genuinely has no detectable blocks.

**Matrix transforms produce fewer candidates than raw.**
Expected. A treated variant may suppress a real block that only
appeared in the raw matrix. The consensus layer (Phase 6) is designed
to handle this: a block in 1 variant is LOW confidence. HIGH confidence
requires it to appear in ≥3 variants.

## See also

- `phase_2_discovery/2c_precomp/README.md` — the precompute step whose
  output this folder consumes
- `phase_2_discovery/2c_precomp/RENAMING.md` — the cross-codebase
  terminology migration tracker (updated to include this folder's
  renames)
- `_archive_old_names/AUDIT_LABELS_DECISIONS_v9.3.md` — full label
  semantics, thresholds, classifier definitions
- `_archive_old_names/INTEGRATION_AUDIT_v9.3.md` — how this folder
  replaces old STEP_C01c and integrates with C01d
- `_archive_old_names/HOW_TO_RUN_v9.3.md` — the original how-to,
  superseded by this README's "Minimum viable run" but kept for its
  tuning-parameter discussion
