# RUNBOOK — producing GHSL and θπ phase-2 JSONs for the atlas

This runbook walks through producing the two new precomp JSON files —
`<chr>_phase2_ghsl.json` and `<chr>_phase2_theta.json` — that the
patched atlas (`pca_scrubber_v4_full_page3.html`) consumes for page 3
(GHSL) and page 12 (θπ).

You start with the dosage-precomp JSON already loaded in the atlas, then
drag-drop the GHSL and θπ enrichment JSONs on top. After this runbook
runs successfully on one chromosome (start with LG28), the atlas's
empty-state placeholders flip to real renderers.

---

## Prerequisites

- Active conda env `assembly` (the `Rscript` binary at
  `/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript`)
- LANTA cluster access, account `lt200308`, partition `compute`
- Existing upstream production outputs:
  - GHSL: `STEP_C04_snake3_ghsl_v6.R` already run, producing per-chromosome
    `*.ghsl_v6.annot.rds`, `*.ghsl_v6.per_sample.rds`, `*.ghsl_v6.karyotypes.rds`
  - GHSL: `STEP_C04b` PASS/WEAK/FAIL classifier already run on annot RDS
  - θπ: ANGSD per-sample `pestPG` files at the configured scale (see
    `00_theta_config.sh` PESTPG_DIR/PESTPG_SCALE)
- Working directory: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery/`

If any of those upstream outputs is missing, this runbook can't proceed
to the JSON exporter step on that chromosome — you'd need to run the
upstream stage first.

---

## Pipeline overview

```
   GHSL track (page 3, 8 layers)                 θπ track (page 12, 5 layers)
   ─────────────────────────────                 ────────────────────────────
   STEP_C04 v6                  (existing)       STEP_TR_A                  (existing)
       │                                              │
       ▼                                              ▼
   STEP_C04b PASS-runs          (existing)       STEP_TR_B v4               (existing)
       │                                              │      writes <chr>_phase2_theta.json
       ▼                                              │      (5 layers including |Z| envelopes)
   STEP_C04c local-PCA          (turn 2 — NEW)        ▼
       │      <chr>.ghsl_v6_localpca.rds          STEP_TR_C D17 wrapper      (turn 4 — NEW)
       │                                              │      writes <chr>_theta_d17{L1,L2}*.tsv
       ▼                                              ▼
   STEP_C04d D17 wrapper        (turn 4 — NEW)   STEP_TR_D augmenter        (turn 5 — NEW)
       │      <chr>_ghsl_d17{L1,L2}*.tsv               in-place augments JSON, adds
       │                                                theta_d17_envelopes layer
       ▼                                              │
   export_ghsl_to_json_v3.R     (turn 5 — NEW)        ▼
          writes <chr>_phase2_ghsl.json           <chr>_phase2_theta.json (5 layers, augmented)
          (8 layers consolidated)
```

Wall time per chromosome (LG28 estimates):

| Stage              | Size    | Time      | Memory |
|--------------------|---------|-----------|--------|
| STEP_C04c          | N=4,300 | ~5–10 s   | <2 GB  |
| STEP_C04d wrapper  | N=4,300 | ~10–30 s  | <2 GB  |
| export_ghsl v3     | —       | ~1 min    | <2 GB  |
| STEP_TR_C wrapper  | N=16,500| ~3–5 min  | ~2 GB  |
| STEP_TR_D          | —       | ~10 s     | <1 GB  |

---

## SECTION 1 — Dry-run on LG28 (single chromosome, interactive)

Pick LG28 first so you can eyeball every output before launching the
28-chromosome SLURM array.

### 1.1 — Set up environment

```bash
source ~/.bashrc
mamba activate assembly
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery

# Convenience var
RSCRIPT=/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript
CHROM=C_gar_LG28
```

### 1.2 — GHSL track: STEP_C04c (local PCA on existing v6 GHSL)

The new turn-2 script reads the v6 production matrices, runs heteroscedastic-
weighted local PCA per window (using `n_phased_het` as weights), sign-aligns
PC1/PC2 across windows, computes `sim_mat = |cor(pc1[i], pc1[j])|`, and
emits |Z|-threshold secondary envelopes. Output: one RDS file per chrom.

```bash
mkdir -p 2e_ghsl_discovery/03_localpca

GHSL_V6_DIR="2e_ghsl_discovery/02_ghsl_v6_outputs"   # whatever your v6 outputs dir is
GHSL_LOCALPCA_DIR="2e_ghsl_discovery/03_localpca"

$RSCRIPT 2e_ghsl_discovery/STEP_C04c_ghsl_local_pca.R \
  --chrom         $CHROM \
  --ghsl_matrices $GHSL_V6_DIR/${CHROM}.ghsl_v6_matrices.rds \
  --out_rds       $GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds \
  --pad           1 \
  --smoothing-scale none \
  --anchor        auto \
  --min-samples-per-window 50 \
  --z-threshold           2.5 \
  --min-l2-windows        5 \
  --merge-gap             3
```

If the matrices RDS isn't at `*.ghsl_v6_matrices.rds` — check what
STEP_C04 actually emitted with `ls $GHSL_V6_DIR/${CHROM}.*.rds` and
adjust the flag. The script reads the file via `readRDS()` and expects
`$div_mat` and `$n_phased_het_mat` slots.

Verify the RDS has reasonable shape:

```bash
$RSCRIPT -e 'lp <- readRDS("'$GHSL_LOCALPCA_DIR'/'$CHROM'.ghsl_v6_localpca.rds")
            cat("n_samples:", lp$n_samples, "\n")
            cat("n_windows:", lp$n_windows, "\n")
            cat("sim_mat dim:", dim(lp$sim_mat), "\n")
            cat("z_profile range:",
                range(lp$z_profile, na.rm=TRUE), "\n")
            cat("n secondary L2:",
                if(!is.null(lp$secondary_l2_envelopes)) nrow(lp$secondary_l2_envelopes) else 0, "\n")'
```

### 1.3 — GHSL track: STEP_C04d (D17 wrapper)

Runs the existing STEP_D17 cross-block boundary detector against the
sim_mat from step 1.2. Writes 4 TSVs: L1 envelopes, L2 envelopes, L1
boundaries, L2 boundaries — these are the D17 cross-block detection
outputs to be loaded by the exporter in step 1.4.

```bash
mkdir -p 2e_ghsl_discovery/04_d17

GHSL_D17_DIR="2e_ghsl_discovery/04_d17"
D17_L1_SCRIPT="/path/to/STEP_D17_multipass_L1_only_v7.R"        # adjust
D17_L2_SCRIPT="/path/to/STEP_D17_multipass_L2_v8.R"             # adjust

$RSCRIPT 2e_ghsl_discovery/STEP_C04d_ghsl_d17_wrapper.R \
  --chrom          $CHROM \
  --localpca       $GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds \
  --out_dir        $GHSL_D17_DIR \
  --d17_l1_script  $D17_L1_SCRIPT \
  --d17_l2_script  $D17_L2_SCRIPT
```

If D17 wall time at GHSL N=4,300 is over a minute, something's wrong —
expected ~5–10 seconds in R. Check the temp `precomp.rds` and `sim_mat.rds`
the wrapper writes inside its temp dir for sanity (the wrapper logs its
temp paths).

If you want to tune the grow validator (the dosage-tuned default may be
too aggressive for fine-grid GHSL), pass:
```bash
  --d17_l1_args "--boundary_grow_W_pct 0.005,0.02,0.05"
```

Verify the 4 TSVs:

```bash
ls -l $GHSL_D17_DIR/${CHROM}_ghsl_d17{L1,L2}_{envelopes,boundaries}.tsv
head -5 $GHSL_D17_DIR/${CHROM}_ghsl_d17L2_envelopes.tsv
```

### 1.4 — GHSL track: export_ghsl_to_json_v3 (consolidate)

Reads everything (STEP_C04 v6 RDSes + STEP_C04c localpca RDS + STEP_C04d
TSVs) and emits the consolidated `<chr>_phase2_ghsl.json` with all 8 layers.

```bash
JSON_OUT=2e_ghsl_discovery/05_json_out
mkdir -p $JSON_OUT

$RSCRIPT 2e_ghsl_discovery/export_ghsl_to_json_v3.R \
  --chrom          $CHROM \
  --annot_rds      $GHSL_V6_DIR/${CHROM}.ghsl_v6.annot.rds \
  --persamp_rds    $GHSL_V6_DIR/${CHROM}.ghsl_v6.per_sample.rds \
  --karyo_rds      $GHSL_V6_DIR/${CHROM}.ghsl_v6.karyotypes.rds \
  --localpca_rds   $GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds \
  --d17_l1_env     $GHSL_D17_DIR/${CHROM}_ghsl_d17L1_envelopes.tsv \
  --d17_l2_env     $GHSL_D17_DIR/${CHROM}_ghsl_d17L2_envelopes.tsv \
  --d17_l1_bnd     $GHSL_D17_DIR/${CHROM}_ghsl_d17L1_boundaries.tsv \
  --d17_l2_bnd     $GHSL_D17_DIR/${CHROM}_ghsl_d17L2_boundaries.tsv \
  --out_dir        $JSON_OUT \
  --primary_scale  s50 \
  --max_k          6
```

Verify all 8 layers emitted:

```bash
$RSCRIPT -e 'js <- jsonlite::fromJSON("'$JSON_OUT'/'$CHROM'/'$CHROM'_phase2_ghsl.json")
            cat("layers_present:\n")
            print(unlist(js$"_layers_present"))
            cat("\nlayer sizes:\n")
            for (ln in unlist(js$"_layers_present")) {
              v <- js[[ln]]
              cat(sprintf("  %-30s %s\n", ln,
                  if(is.list(v)) paste(length(v),"items") else paste("len",length(v))))
            }
            cat("\nfile size: ",
                file.info("'$JSON_OUT'/'$CHROM'/'$CHROM'_phase2_ghsl.json")$size,"bytes\n")'
```

Expect 8 layers: `tracks`, `ghsl_panel`, `ghsl_kstripes`,
`ghsl_karyotype_runs`, `ghsl_envelopes`, `ghsl_local_pca`,
`ghsl_secondary_envelopes`, `ghsl_d17_envelopes`. File size will be
~5–30 MB depending on the number of envelopes.

### 1.5 — θπ track: STEP_TR_A then STEP_TR_B (existing pipeline)

The existing θπ stream stays unchanged. STEP_TR_A reads ANGSD pestPG and
emits the per-sample TSV; STEP_TR_B reads it and emits the 4-layer JSON.

```bash
cd 2f_theta_discovery
source 00_theta_config.sh   # exports THETA_TSV_DIR, JSON_OUT_DIR, PESTPG_*

$RSCRIPT STEP_TR_A_compute_theta_matrices.R $CHROM
$RSCRIPT STEP_TR_B_classify_theta.R         $CHROM
```

Verify intermediate JSON has the 4 baseline layers:

```bash
$RSCRIPT -e 'js <- jsonlite::fromJSON("'$JSON_OUT_DIR'/'$CHROM'/'$CHROM'_phase2_theta.json")
            print(unlist(js$"_layers_present"))'
```

Expect: `tracks`, `theta_pi_per_window`, `theta_pi_local_pca`,
`theta_pi_envelopes`. (5th `theta_d17_envelopes` lands after the next
two steps.)

### 1.6 — θπ track: STEP_TR_C (D17 wrapper, NEW turn 4)

```bash
THETA_D17_DIR=$OUTROOT/05_d17
mkdir -p $THETA_D17_DIR

$RSCRIPT STEP_TR_C_theta_d17_wrapper.R \
  --chrom          $CHROM \
  --theta_json     $JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json \
  --out_dir        $THETA_D17_DIR \
  --d17_l1_script  $D17_L1_SCRIPT \
  --d17_l2_script  $D17_L2_SCRIPT
```

This is the slow step (~3–5 min at N=16,500). The wrapper computes
`sim_mat = abs(cor(pc1_mat))` cluster-side from STEP_TR_B's
`theta_pi_local_pca.pc1_loadings` field, then invokes D17.

If it's too slow, pass `--d17_l1_args "--boundary_grow_W_pct 0.005,0.02,0.05"`
for ~2.4× speedup. (See STEP_C04d header for explanation.)

### 1.7 — θπ track: STEP_TR_D (augment JSON, NEW turn 5)

```bash
$RSCRIPT STEP_TR_D_augment_theta_json.R \
  --theta_json   $JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json \
  --d17_l1_env   $THETA_D17_DIR/${CHROM}_theta_d17L1_envelopes.tsv \
  --d17_l2_env   $THETA_D17_DIR/${CHROM}_theta_d17L2_envelopes.tsv \
  --d17_l1_bnd   $THETA_D17_DIR/${CHROM}_theta_d17L1_boundaries.tsv \
  --d17_l2_bnd   $THETA_D17_DIR/${CHROM}_theta_d17L2_boundaries.tsv
```

(Default: writes back in place; pass `--out_json <other>` to write
elsewhere if you want to keep the pre-augmented version.)

Verify final 5 layers:

```bash
$RSCRIPT -e 'js <- jsonlite::fromJSON("'$JSON_OUT_DIR'/'$CHROM'/'$CHROM'_phase2_theta.json")
            print(unlist(js$"_layers_present"))'
```

Expect: `tracks`, `theta_pi_per_window`, `theta_pi_local_pca`,
`theta_pi_envelopes`, `theta_d17_envelopes`.

### 1.8 — Smoke test in atlas

Open `pca_scrubber_v4_full_page3.html` in a browser. Click `Browse...`,
load your existing **dosage** precomp JSON (the one you've been using
for page 1). Then drag-drop the two new JSONs as enrichment:

1. `<chr>_phase2_ghsl.json` — page-3 layer indicators flip ⚪→🟢, the
   GHSL option in the line-color-mode picker auto-enables, page 3 default
   mode renders the sim_mat heatmap with three envelope-layer overlays.

2. `<chr>_phase2_theta.json` — page-12 layer indicators flip ⚪→🟢.
   (Page 12's panel renderers are not yet wired by this session — only
   the empty-state status indicators.)

Watch the browser console (F12). If you see "ghsl_envelopes is
undefined" or similar, the merge-fix patch (turn 8b) wasn't applied.
Re-run the patcher chain.

---

## SECTION 2 — Full 28-chromosome SLURM array (after dry-run passes)

Once LG28 produces sane output and the atlas renders correctly, scale
to all 28 chromosomes via SLURM array. Submit two arrays in parallel
(GHSL and θπ are independent).

### 2.1 — GHSL SLURM array

Save as `2e_ghsl_discovery/run_ghsl_array.sh`:

```bash
#!/bin/bash
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-28
#SBATCH --output=logs/ghsl_%A_%a.out
#SBATCH --error=logs/ghsl_%A_%a.err

source ~/.bashrc
mamba activate assembly

CHROM_LIST=(C_gar_LG01 C_gar_LG02 C_gar_LG03 C_gar_LG04 C_gar_LG05
            C_gar_LG06 C_gar_LG07 C_gar_LG08 C_gar_LG09 C_gar_LG10
            C_gar_LG11 C_gar_LG12 C_gar_LG13 C_gar_LG14 C_gar_LG15
            C_gar_LG16 C_gar_LG17 C_gar_LG18 C_gar_LG19 C_gar_LG20
            C_gar_LG21 C_gar_LG22 C_gar_LG23 C_gar_LG24 C_gar_LG25
            C_gar_LG26 C_gar_LG27 C_gar_LG28)
CHROM=${CHROM_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

RSCRIPT=/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript
BASE=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery
cd $BASE

GHSL_V6_DIR=2e_ghsl_discovery/02_ghsl_v6_outputs
GHSL_LOCALPCA_DIR=2e_ghsl_discovery/03_localpca
GHSL_D17_DIR=2e_ghsl_discovery/04_d17
JSON_OUT=2e_ghsl_discovery/05_json_out
D17_L1_SCRIPT=/path/to/STEP_D17_multipass_L1_only_v7.R
D17_L2_SCRIPT=/path/to/STEP_D17_multipass_L2_v8.R

mkdir -p $GHSL_LOCALPCA_DIR $GHSL_D17_DIR $JSON_OUT

set -euo pipefail
echo "[$(date)] [$CHROM] STEP_C04c"
$RSCRIPT 2e_ghsl_discovery/STEP_C04c_ghsl_local_pca.R \
  --chrom $CHROM \
  --ghsl_matrices $GHSL_V6_DIR/${CHROM}.ghsl_v6_matrices.rds \
  --out_rds       $GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds

echo "[$(date)] [$CHROM] STEP_C04d"
$RSCRIPT 2e_ghsl_discovery/STEP_C04d_ghsl_d17_wrapper.R \
  --chrom $CHROM \
  --localpca $GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds \
  --out_dir $GHSL_D17_DIR \
  --d17_l1_script $D17_L1_SCRIPT --d17_l2_script $D17_L2_SCRIPT

echo "[$(date)] [$CHROM] export_ghsl_to_json_v3"
$RSCRIPT 2e_ghsl_discovery/export_ghsl_to_json_v3.R \
  --chrom $CHROM \
  --annot_rds    $GHSL_V6_DIR/${CHROM}.ghsl_v6.annot.rds \
  --persamp_rds  $GHSL_V6_DIR/${CHROM}.ghsl_v6.per_sample.rds \
  --karyo_rds    $GHSL_V6_DIR/${CHROM}.ghsl_v6.karyotypes.rds \
  --localpca_rds $GHSL_LOCALPCA_DIR/${CHROM}.ghsl_v6_localpca.rds \
  --d17_l1_env $GHSL_D17_DIR/${CHROM}_ghsl_d17L1_envelopes.tsv \
  --d17_l2_env $GHSL_D17_DIR/${CHROM}_ghsl_d17L2_envelopes.tsv \
  --d17_l1_bnd $GHSL_D17_DIR/${CHROM}_ghsl_d17L1_boundaries.tsv \
  --d17_l2_bnd $GHSL_D17_DIR/${CHROM}_ghsl_d17L2_boundaries.tsv \
  --out_dir $JSON_OUT

echo "[$(date)] [$CHROM] DONE"
```

Submit:

```bash
sbatch 2e_ghsl_discovery/run_ghsl_array.sh
squeue -u $USER
```

### 2.2 — θπ SLURM array

Save as `2f_theta_discovery/run_theta_array.sh`:

```bash
#!/bin/bash
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-28
#SBATCH --output=logs/theta_%A_%a.out
#SBATCH --error=logs/theta_%A_%a.err

source ~/.bashrc
mamba activate assembly
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery

source 00_theta_config.sh   # exports CHROM_LIST and all paths
CHROM=${CHROM_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

D17_L1_SCRIPT=/path/to/STEP_D17_multipass_L1_only_v7.R
D17_L2_SCRIPT=/path/to/STEP_D17_multipass_L2_v8.R
THETA_D17_DIR=$OUTROOT/05_d17
mkdir -p $THETA_D17_DIR

set -euo pipefail
echo "[$(date)] [$CHROM] STEP_TR_A"
$RSCRIPT STEP_TR_A_compute_theta_matrices.R $CHROM

echo "[$(date)] [$CHROM] STEP_TR_B"
$RSCRIPT STEP_TR_B_classify_theta.R $CHROM

echo "[$(date)] [$CHROM] STEP_TR_C"
$RSCRIPT STEP_TR_C_theta_d17_wrapper.R \
  --chrom $CHROM \
  --theta_json $JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json \
  --out_dir $THETA_D17_DIR \
  --d17_l1_script $D17_L1_SCRIPT --d17_l2_script $D17_L2_SCRIPT

echo "[$(date)] [$CHROM] STEP_TR_D"
$RSCRIPT STEP_TR_D_augment_theta_json.R \
  --theta_json $JSON_OUT_DIR/$CHROM/${CHROM}_phase2_theta.json \
  --d17_l1_env $THETA_D17_DIR/${CHROM}_theta_d17L1_envelopes.tsv \
  --d17_l2_env $THETA_D17_DIR/${CHROM}_theta_d17L2_envelopes.tsv \
  --d17_l1_bnd $THETA_D17_DIR/${CHROM}_theta_d17L1_boundaries.tsv \
  --d17_l2_bnd $THETA_D17_DIR/${CHROM}_theta_d17L2_boundaries.tsv

echo "[$(date)] [$CHROM] DONE"
```

Submit:

```bash
sbatch 2f_theta_discovery/run_theta_array.sh
squeue -u $USER
```

### 2.3 — Verify the array completed cleanly

```bash
# Count successful chromosomes
ls 2e_ghsl_discovery/05_json_out/*/C_gar_LG??_phase2_ghsl.json | wc -l
ls $JSON_OUT_DIR/*/C_gar_LG??_phase2_theta.json | wc -l
# Both should be 28.

# Check each JSON has all expected layers
for chr in C_gar_LG{01..28}; do
  $RSCRIPT -e 'js <- jsonlite::fromJSON("'2e_ghsl_discovery/05_json_out/$chr/${chr}_phase2_ghsl.json'")
              cat("'$chr': ", length(unlist(js$"_layers_present")), "layers\n")' 2>/dev/null
done
```

Should show `8 layers` for every chromosome on the GHSL side and
`5 layers` on the θπ side.

---

## SECTION 3 — Loading into the atlas

### 3.1 — Start with the patched atlas

Use `pca_scrubber_v4_full_page3.html` from this session's outputs.
This is the v4 base + 4 patches applied:

```
v4 base (md5 37eb28da..)
  + apply_page3_scaffold.py        (turn 6: page 15 tab + DOM + layer detection)
  + apply_ghsl_color_mode.py       (turn 7: COLORMAP_KSTRIPE + GHSL color resolver)
  + apply_page3_panels.py          (turn 8: 5 canvases + drawGhSim + 4 stubs + cursor sync)
  + apply_enrichment_merge_fix.py  (turn 8b: 7 new layer cases in merge switch)
= pca_scrubber_v4_full_page3.html (md5 c6c3ba48..)
```

If you ever need to re-apply on a fresh v4 base from scratch, the chain is:

```bash
python3 2e_ghsl_discovery/apply_page3_scaffold.py        --atlas pca_scrubber_v4.html      --out a1.html
python3 2e_ghsl_discovery/apply_ghsl_color_mode.py       --atlas a1.html --out a2.html --allow-base-mismatch
python3 2e_ghsl_discovery/apply_page3_panels.py          --atlas a2.html --out a3.html --allow-base-mismatch
python3 2e_ghsl_discovery/apply_enrichment_merge_fix.py  --atlas a3.html --out final.html --allow-base-mismatch
```

All four patches are idempotent.

### 3.2 — Drag-drop workflow

1. Open the patched HTML in a browser.
2. Sidebar `Browse...` → load your existing **dosage** precomp JSON
   (the one you've been using for page 1 inversion calling). This
   defines the dosage windows that the GHSL/θπ layers map onto.
3. Page 4 → "+ load enrichment" button → select `<chr>_phase2_ghsl.json`
   for the matching chrom. Page 3 (tab "15 local PCA GHSL") layer
   indicators flip ⚪→🟢. The line-color-mode picker's GHSL option
   becomes selectable. Page 1 lines and PCA can now be colored by
   GHSL kstripe.
4. Same enrichment button → select `<chr>_phase2_theta.json`. Page 12
   layer indicators flip ⚪→🟢. (Page 12 panel renderers are scaffold-
   only this session; only the empty-state indicators update.)

### 3.3 — Sanity checks in browser

Open the dev console (F12). Useful queries:

```javascript
// What layers are loaded?
Array.from(state.layersPresent).sort()
// → should include ghsl_local_pca, ghsl_envelopes, ghsl_secondary_envelopes,
//   ghsl_d17_envelopes, theta_pi_local_pca, theta_pi_envelopes, theta_d17_envelopes

// Inspect GHSL sim_mat shape
state.data.ghsl_local_pca.sim_mat_n
state.data.ghsl_local_pca.sim_mat.length
// → second number should equal n*(n+1)/2 for the first

// Check window-mapping built
state.ghslWindowMap?.length   // = state.data.windows.length (dosage N)
state.ghslWindowCount         // = ghsl panel n_windows (M)
```

---

## Troubleshooting

**`STEP_C04c says n_phased_het_mat slot missing`**: The v6 production
RDS may have a different field name. Check `names(readRDS(...))`. Map
the actual name to what STEP_C04c expects with the `--weight-field-name`
flag (if present in the script header), or rename in the RDS.

**`STEP_C04d D17 reports zero peaks`**: This is informative output, not
an error — it means the GHSL sim_mat has no detectable cross-block
boundaries. May be real (chromosome-wide noise) or a threshold issue.
Try `--d17_l1_args "--boundary_score_min 1.0"` for a more permissive
threshold.

**`exporter v3 says "annot RDS structure not recognized"`**: STEP_C04
v6 may emit annot data inside a list with a different slot name. The
exporter probes `annot$dt`, `annot$annot`, and the bare data.frame.
Check `str(readRDS("...annot.rds"))` to see what's actually there and
edit the exporter's `# Load STEP_C04 v6 annot RDS` block.

**`STEP_TR_C says "JSON missing theta_pi_local_pca.pc1_loadings"`**:
STEP_TR_B emits `theta_pi_local_pca` with `pc1_loadings` as a flat array
(per-window-flattened), not nested. If STEP_TR_C expects nested, you'll
see this error. Both shapes are valid — check current STEP_TR_B
emit format.

**`Atlas says no layers found after enrichment merge`**: The merge-fix
patch (turn 8b) wasn't applied. Verify with:
```bash
grep -c "turn 8b: GHSL page-3 layers" pca_scrubber_v4_full_page3.html
```
Should print `1`.

**`Page 3 panels stay empty even though layer indicators are 🟢`**:
The layer name was added to `state.layersPresent` but the corresponding
data wasn't copied to `state.data`. This is the merge-fix bug if the
patch is missing. If it IS present, check field names — particularly
that `state.data.ghsl_panel.start_bp` is an array (the cursor-sync
helper reads from there).

---

## Files in this bundle

```
2e_ghsl_discovery/
├── STEP_C04c_ghsl_local_pca.R              (turn 2)  local PCA on GHSL v6 matrices
├── STEP_C04d_ghsl_d17_wrapper.R            (turn 4)  D17 cross-block boundary wrapper
├── export_ghsl_to_json_v3.R                (turn 5)  consolidates 4 inputs → 8-layer JSON
├── apply_page3_scaffold.py                 (turn 6)  atlas tab + DOM + layer detection
├── apply_ghsl_color_mode.py                (turn 7)  COLORMAP_KSTRIPE + GHSL color resolver
├── apply_page3_panels.py                   (turn 8)  5 canvases + drawGhSim + 4 stubs
├── apply_enrichment_merge_fix.py           (turn 8b) merge-switch fix for new layers
├── test_*.py (5 files)                     50/50 structural tests
└── README.md

2f_theta_discovery/
├── STEP_TR_A_compute_theta_matrices.R       (existing, heavy precomp)
├── STEP_TR_B_classify_theta.R               (existing v4, 4-layer JSON emit)
├── STEP_TR_C_theta_d17_wrapper.R            (turn 4)  D17 wrapper for θπ
├── STEP_TR_D_augment_theta_json.R           (turn 5)  in-place augment with D17 layer
├── 00_theta_config.sh                       (existing, single source of truth)
└── README.md

[Atlas variants]
pca_scrubber_v4_full_page3.html             (v4 base + all 4 turn-6/7/8/8b patches)

[Documentation]
phase2_ghsl_arrangement_v1.md               architectural commit, v1 plan
schema_v2_addendum_ghsl_page3.md            schema addendum proposal §22b
RUNBOOK_produce_phase2_jsons.md             this file
```
