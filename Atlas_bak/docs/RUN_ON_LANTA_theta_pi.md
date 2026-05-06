# RUN_ON_LANTA_theta_pi.md — practical 28-chromosome θπ production

**One-page guide.** Step-by-step. Assumes you just dropped the Atlas folder
onto LANTA and nothing else. Read top-to-bottom; copy-paste each block in
order. Do not skip the sanity check.

This produces 28 × `<CHROM>_phase2_theta.json` files that the atlas
loads into page 12 (local PCA θπ).

Companion docs (already in this folder, deeper detail):
- `cluster_R/2f_theta_discovery/README_theta_pi_scaling.md` — the bug fix
  context (why we divide tP by nSites)
- `RUNBOOK_produce_phase2_jsons.md` — the full pipeline overview
- `cluster_R/2f_theta_discovery/00_theta_config.sh` — every parameter

---

## What this does

For each of 28 chromosomes (LG01..LG28), one SLURM array task runs:

```
STEP_TR_A   reads 226 pestPG files → one wide TSV at θπ-native grid
STEP_TR_B   local PCA + |Z| → writes <CHROM>_phase2_theta.json
STEP_TR_C   D17 multi-pass on the JSON's sim_mat (NEEDS D17 SCRIPTS — see step 4)
STEP_TR_D   in-place augments the JSON with theta_d17_envelopes
```

Output: `${OUTROOT}/json_out/<CHROM>/<CHROM>_phase2_theta.json` × 28.
Each file is ~5–10 MB, drag-droppable into the atlas.

---

## Step 0 — landing checks

After dragging the Atlas folder to LANTA, confirm where it sits and pick
a working location for the SLURM job. The recommended layout follows the
pattern used by sibling modules:

```bash
# 1. Where is the Atlas folder right now?
pwd
ls -la

# 2. Move (or symlink) the cluster_R contents into the canonical
#    discovery codebase tree, where the config expects it:
PROJ=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04
DEST=$PROJ/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery
mkdir -p "$DEST"

# Copy the four R scripts + sanity check + array launcher + config
cp -v Atlas/_scripts/cluster_R/2f_theta_discovery/*.R       "$DEST/"
cp -v Atlas/_scripts/cluster_R/2f_theta_discovery/*.sh      "$DEST/"
cp -v Atlas/_scripts/cluster_R/2f_theta_discovery/*.md      "$DEST/"

# Spot-check
ls -la "$DEST/"
```

You should see:
- `STEP_TR_A_compute_theta_matrices.R` (with the per-site scaling fix, --chrom CLI)
- `STEP_TR_B_classify_theta.R` (two-mode JSON emit, --chrom CLI)
- `STEP_TR_C_theta_d17_wrapper.R`
- `STEP_TR_D_augment_theta_json.R`
- `00_theta_config.sh`
- `00_sanity_check_pestPG_scaling.sh`
- `LAUNCH_TR_theta_pi.slurm` (B01-style launcher — primary path)
- `chrom.list` (companion list — header + 28 chromosomes)
- `run_theta_array.sh` (legacy full pipeline incl. D17)
- `run_theta_array_NO_D17.sh` (legacy slim — fallback if central config absent)
- `README_theta_pi_scaling.md` (the bug fix context)

---

## Step 1 — confirm pestPG data is reachable

The patched STEP_TR_A reads from `${PESTPG_DIR}/${SAMPLE}.${PESTPG_SCALE}.pestPG`.
Per the config, that's:

```
PESTPG_DIR   = ${BASE}/het_roh/02_heterozygosity/03_theta/multiscale
PESTPG_SCALE = win10000.step2000
```

Check that the files exist on LANTA:

```bash
cd $PROJ/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery
source 00_theta_config.sh

# Print resolved paths
echo "PESTPG_DIR  = $PESTPG_DIR"
echo "Looking for: $PESTPG_DIR/*.${PESTPG_SCALE}.pestPG"

# Count files (expect 226 at win10000.step2000)
ls "$PESTPG_DIR"/*."${PESTPG_SCALE}".pestPG 2>/dev/null | wc -l
```

**If the count is 0 or far below 226**, the data isn't on LANTA at the
expected path. Two options:

(a) **Upload from your local machine first.** From your laptop:
```bash
# Run from your laptop, NOT from LANTA
rsync -avh --progress \
    /mnt/e/01-catfish_assembly_manuscript_CGA/het_roh_02_heterozygosity_2026-04-20/02_heterozygosity/03_theta/multiscale/ \
    USER@LANTA-HOST:$PROJ/het_roh/02_heterozygosity/03_theta/multiscale/
```
Total = ~52 GB at win10000.step2000 + win5000.step1000 + win50000.step10000.
**At win10000.step2000 alone (which is what we need): ~16 GB. Upload only that scale to save time:**
```bash
rsync -avh --progress \
    /mnt/e/.../03_theta/multiscale/*.win10000.step2000.pestPG \
    USER@LANTA-HOST:$PROJ/het_roh/02_heterozygosity/03_theta/multiscale/
```

(b) **Adjust PESTPG_DIR in the config** if the data is already on LANTA at
a different path. Edit `00_theta_config.sh` line 76.

Re-run the count check until it shows 226.

---

## Step 2 — sanity check: confirm the bug exists in your real data

**Do not skip this.** Takes 30 seconds. Confirms that the per-site scaling
fix in STEP_TR_A is actually doing the right thing on YOUR pestPG files,
before we burn 28 SLURM jobs.

```bash
cd $PROJ/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery
chmod +x 00_sanity_check_pestPG_scaling.sh
bash 00_sanity_check_pestPG_scaling.sh
```

Read the output. The key lines:
```
mean(tP)            : <number>
mean(nSites)        : <number>
mean(tP / nSites)   : <number>
```

**Expected:**
- `mean(tP)` between 0.5 and 50 (sum-scale — the bug is real)
- `mean(tP / nSites)` between 1e-4 and 1e-2 (per-site, vertebrate-typical)
- Verdict line says: `tP looks like a SUM. Patch STEP_TR_A.`

If `mean(tP)` is already in the 1e-3 range, ANGSD changed behavior or
this isn't standard pestPG output — STOP and ask before running, the
patch would over-divide.

---

## Step 3 — pick a SLURM strategy

Two options:

### Option A: Skip D17 first (fastest path — recommended for first run)

D17 is the boundary refinement step that needs two scripts I don't have
canonical paths for (`STEP_D17_multipass_L1_only_v7.R`,
`STEP_D17_multipass_L2_v8.R`). The atlas works WITHOUT D17 — page 12
just won't have the `theta_d17_envelopes` 5th layer; the 4 core layers
(`theta_pi_per_window`, `theta_pi_local_pca`, `theta_pi_envelopes`,
plus tracks contribution) are what page 12 actually renders from.

**Run STEP_TR_A + STEP_TR_B only**. This produces a usable JSON with 4 of
5 layers loaded; you can drag it into the atlas right away.

### Option B: Full pipeline including D17

Requires you to find/place the D17 scripts and edit lines 41–42 of
`run_theta_array.sh` to point at them.

---

## Step 4 — submit the array (Option A path)

This θπ pipeline uses the same SLURM launcher idiom as `LAUNCH_B01_mds_stage1.slurm`
(the canonical 2b_mds inversion-codebase launcher): `sbatch --array=N-M LAUNCH_*.slurm chrom.list`.

The launcher (`LAUNCH_TR_theta_pi.slurm`) and its companion `chrom.list` are
already in this directory from Step 0. The list has a header line, then 28
chromosomes (one per row) — `SLURM_ARRAY_TASK_ID=N` reads line N+1, exactly
as B01 does. Submit all 28 chromosomes with:

```bash
sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list
```

You should get back: `Submitted batch job <JOB_ID>`.

Watch the queue:
```bash
squeue -u $USER
```

Each task should finish in 5–15 minutes. With 28 tasks running in parallel
(LANTA queue permitting) total wall time is ~10–20 minutes.

**To run a single chromosome** (e.g. just LG28 for a dry run, position 28
in chrom.list):

```bash
sbatch --array=28 LAUNCH_TR_theta_pi.slurm chrom.list
```

**Fallback if `00_inversion_config.sh` is missing on LANTA**: the older
`run_theta_array_NO_D17.sh` is still in the directory. It hard-codes paths
and doesn't depend on the central config. Submit with `sbatch run_theta_array_NO_D17.sh`.

---

## Step 5 — verify the output

While running:
```bash
# Tail one log to see it progress
tail -f logs/theta_*_1.out
```

After all tasks complete:
```bash
# Are all 28 JSONs there?
ls $JSON_OUT_DIR/*/C_gar_LG*_phase2_theta.json | wc -l    # expect 28

# Spot-check one with jq (structural dump):
jq 'keys, .schema_version, ._layers_present, .n_windows' \
   $JSON_OUT_DIR/C_gar_LG28/C_gar_LG28_phase2_theta.json
```

Expected JSON keys include `theta_pi_per_window`, `theta_pi_local_pca`,
`theta_pi_envelopes`, `available_modes` (`["per_site", "raw_sum"]`),
`default_mode` (`"per_site"`).

---

## Step 6 — drag into the atlas

On your laptop, sync the JSONs back:
```bash
# From laptop
rsync -avh USER@LANTA-HOST:${JSON_OUT_DIR}/*/C_gar_LG*_phase2_theta.json \
    ./theta_jsons/
```

Open `Atlas/Inversion_atlas.html` in a browser. Drag the LG28 dosage precomp JSON
in first to populate page 1. Then drag `C_gar_LG28_phase2_theta.json` on
top — it should fly into the green schema badge with the `+4 layers loaded`
animation, and the badge counter should jump from `4 layers` to `8 layers`
(or whatever the dosage precomp had + 4).

Switch to tab "2 local PCA θπ" — the empty-state should now light up the
loaded layers (⚪ → 🟢) and the panels should populate.

---

## Failure modes — what to look at if a task fails

If a task fails, check `logs/theta_<JOB_ID>_<TASK_ID>.err`:

| Error | Likely cause |
|---|---|
| `Could not load any sample's pestPG for <CHROM>` | PESTPG_DIR wrong, or files missing for that chrom |
| `Error in fread: ... not enough columns` | pestPG file truncated/malformed for one sample |
| `cannot allocate vector of size N MB` | Bump `--mem=64G` and resubmit |
| `STEP_TR_B: Error in svd: ...` | Local PCA crash — usually a chrom with too few callable windows |
| Job pending 2 hours | LANTA queue full; consider `--partition=memory` if available |

If a single chrom fails, re-run just that one:
```bash
sbatch --array=28 run_theta_array_NO_D17.sh   # just LG28
```

---

## What about D17? (optional, do later)

After you have the 4-layer JSONs working in the atlas and you want the
5th `theta_d17_envelopes` layer:

1. Locate `STEP_D17_multipass_L1_only_v7.R` and `STEP_D17_multipass_L2_v8.R`
   (these are from MODULE_2A/2B's local-PCA suite — Quentin knows where).
2. Edit `run_theta_array.sh` lines 41–42 with the absolute paths.
3. Run on already-produced JSONs:
   ```bash
   for chrom in C_gar_LG{01..28}; do
     $RSCRIPT STEP_TR_C_theta_d17_wrapper.R \
       --chrom "$chrom" \
       --theta_json "$JSON_OUT_DIR/$chrom/${chrom}_phase2_theta.json" \
       --out_dir "$OUTROOT/05_d17" \
       --d17_l1_script  "$D17_L1_SCRIPT" \
       --d17_l2_script  "$D17_L2_SCRIPT"
     $RSCRIPT STEP_TR_D_augment_theta_json.R \
       --theta_json "$JSON_OUT_DIR/$chrom/${chrom}_phase2_theta.json" \
       --d17_l1_env "$OUTROOT/05_d17/${chrom}_theta_d17L1_envelopes.tsv" \
       --d17_l2_env "$OUTROOT/05_d17/${chrom}_theta_d17L2_envelopes.tsv" \
       --d17_l1_bnd "$OUTROOT/05_d17/${chrom}_theta_d17L1_boundaries.tsv" \
       --d17_l2_bnd "$OUTROOT/05_d17/${chrom}_theta_d17L2_boundaries.tsv"
   done
   ```

The augmenter is in-place — re-loading the JSON in the atlas after this
should add the 5th layer to the schema badge count.

---

## tl;dr — the four commands

```bash
# 1. land scripts on LANTA
cp -v Atlas/_scripts/cluster_R/2f_theta_discovery/* \
   $PROJ/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery/

# 2. cd + source + sanity-check
cd $PROJ/inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery
source 00_theta_config.sh
bash 00_sanity_check_pestPG_scaling.sh   # confirm bug, ~30 sec

# 3. submit the 28-chrom array (B01 idiom: array N maps to chrom.list line N+1)
sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list

# 4. wait ~15 min, then verify
ls $JSON_OUT_DIR/*/C_gar_LG*_phase2_theta.json | wc -l   # expect 28
```

That's it. Drag the JSONs into the atlas, page 12 lights up.
