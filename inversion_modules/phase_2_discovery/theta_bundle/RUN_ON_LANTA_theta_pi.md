# Run θπ on LANTA — step-by-step

This bundle contains the patched θπ pipeline scripts and a preflight checker.
Both fixes (pestPG per-site scaling + atlas-canonical JSON shape) are baked in.

The bundle places everything under `inversion_modules/phase_2_discovery/2f_theta_discovery/`. Atlas stays viewer-only — no scripts under `Atlas/`.

---

## Path conventions

**LOCAL** (your laptop):
```
C:\Users\quent\Desktop\inversion-popgen-toolkit\inversion_modules\phase_2_discovery\2f_theta_discovery\
```

**LANTA** (mirror):
```
/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit/inversion_modules/phase_2_discovery/2f_theta_discovery/
```

The principle: **`Atlas/` is the viewer**, **`inversion_modules/phase_*/` is where pipelines live**, organized by phase. Pages 1, 2, 3, 12 are all rendered by the atlas, but the scripts that produce their JSONs all live under `inversion_modules/phase_2_discovery/`.

---

## Step 1 — Push this bundle to LANTA

From your laptop:

```bash
# tar up the bundle and ship it
tar -czf theta_bundle.tar.gz -C theta_bundle_v2 .
scp theta_bundle.tar.gz lanta:~/

# on LANTA — extract at the repo root
ssh lanta
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit
tar -xzf ~/theta_bundle.tar.gz

# Files land at:
#   inversion_modules/phase_2_discovery/2f_theta_discovery/
#       STEP_TR_A_compute_theta_matrices.R   (per-site fix)
#       STEP_TR_B_classify_theta.R           (turn-114 schema-v2 patch)
#       STEP_TR_C_theta_d17_wrapper.R        (optional D17 envelope detector)
#       STEP_TR_D_augment_theta_json.R       (optional D17 → JSON augmenter)
#       00_sanity_check_pestPG_scaling.sh
#       00_theta_config.sh                   (CODEBASE points at new layout, fallback to legacy)
#       LAUNCH_TR_theta_pi.slurm             (config resolution: new path first, legacy fallback)
#       chrom.list
#       README_theta_pi_scaling.md
#       verify_theta_setup.sh                (preflight check)
```

If you have the legacy `inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery/` tree on LANTA and want to keep the old TR_A / TR_B versions there:
```bash
cp /scratch/.../inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery/STEP_TR_A_compute_theta_matrices.R \
   /scratch/.../inversion_codebase_v8.5/phase_2_discovery/2f_theta_discovery/STEP_TR_A_compute_theta_matrices.R.bak.$(date +%F)
# (similar for TR_B)
```
Note: the **new** location overrides the old one — you don't need to copy patched scripts into the legacy tree. Just leave the legacy tree alone.

---

## Step 2 — Preflight check

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit/inversion_modules/phase_2_discovery/2f_theta_discovery
mamba activate assembly
bash verify_theta_setup.sh
```

Checks 7 things and fails loudly if any aren't ready:

1. You're in the right directory
2. All required files are present
3. `00_theta_config.sh` sources cleanly and exports the expected vars
4. `PESTPG_DIR` exists and contains ~226 `.pestPG` files at `PESTPG_SCALE`
5. `STEP_TR_A` has the per-site scaling fix (looks for `tP / nSites`)
6. `STEP_TR_B` has the turn-114 schema-v2 patch (looks for `schema_version = 2L`)
7. Runs `00_sanity_check_pestPG_scaling.sh` on a real pestPG to confirm `tP` is a SUM (`VERDICT: SUM`)

Expected end-of-output:
```
STATUS: READY ✓
```

The check submits nothing to SLURM — read-only.

---

## Step 3 — Dry run on LG28

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit/inversion_modules/phase_2_discovery/2f_theta_discovery
source 00_theta_config.sh
RSCRIPT="$(which Rscript)"

# build the per-sample × per-window TSV (3-5 min)
$RSCRIPT STEP_TR_A_compute_theta_matrices.R --chrom C_gar_LG28

# emit the 4-layer JSON (5-10 min)
$RSCRIPT STEP_TR_B_classify_theta.R --chrom C_gar_LG28
```

Verify:
```bash
JSON="${JSON_OUT_DIR}/C_gar_LG28/C_gar_LG28_phase2_theta.json"
ls -lh "$JSON"
jq '.schema_version, ._layers_present, .n_windows, (.theta_pi_per_window.values | length)' "$JSON"
```

Expected:
```
2
[ "theta_pi_per_window", "theta_pi_local_pca", "theta_pi_envelopes", "tracks" ]
16500
3729000
```

(`16500 × 226 = 3,729,000`, the flat row-major matrix length.)

If `schema_version` is `1`, the patched TR_B isn't being used. Re-check step 1.

---

## Step 4 — Full 28-chromosome SLURM array

```bash
sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list

# monitor
squeue -u $USER
tail -f logs/theta_pi.*.out
```

Each task runs `STEP_TR_A` then `STEP_TR_B` for one chromosome. Walltime per task: 5-15 min on LANTA.

`SLURM_ARRAY_TASK_ID=N` reads line `N+1` of `chrom.list` (line 1 is the header). `--array=1-28` covers `C_gar_LG01` through `C_gar_LG28`.

After all tasks finish, JSONs land at:
```
${JSON_OUT_DIR}/C_gar_LG01/C_gar_LG01_phase2_theta.json
...
${JSON_OUT_DIR}/C_gar_LG28/C_gar_LG28_phase2_theta.json
```

---

## Step 5 — (Optional) Add D17 envelopes

Skip on first pass. Only needed for the `theta_d17_envelopes` layer.

```bash
$RSCRIPT STEP_TR_C_theta_d17_wrapper.R --chrom C_gar_LG28

$RSCRIPT STEP_TR_D_augment_theta_json.R \
    --theta_json "${JSON_OUT_DIR}/C_gar_LG28/C_gar_LG28_phase2_theta.json" \
    --d17_l1_env "${OUTROOT}/d17/C_gar_LG28_theta_d17L1_envelopes.tsv" \
    --d17_l2_env "${OUTROOT}/d17/C_gar_LG28_theta_d17L2_envelopes.tsv" \
    --d17_l1_bnd "${OUTROOT}/d17/C_gar_LG28_theta_d17L1_boundaries.tsv" \
    --d17_l2_bnd "${OUTROOT}/d17/C_gar_LG28_theta_d17L2_boundaries.tsv"
```

(Filenames may vary — `ls "${OUTROOT}/d17/"` after TR_C runs to confirm.)

---

## Step 6 — Run STEP_T05 CUSUM on the JSON

Once the JSON exists, run the per-sample CUSUM observation step from earlier in this session:

```bash
$RSCRIPT /path/to/STEP_T05_theta_cusum.R \
    --json "${JSON_OUT_DIR}/C_gar_LG28/C_gar_LG28_phase2_theta.json" \
    --out-dir /tmp/cusum_LG28/ \
    --mode whole_chrom \
    --lib /path/to/lib_persample_cusum.R
```

`--mode whole_chrom` gives one cp_bp per sample, anywhere on the chromosome. Outputs:
- `theta_cusum_per_sample.tsv.gz`
- `theta_cusum_summary.tsv`

Look at these to see what the per-sample changepoint distribution actually looks like, before designing consensus / atlas rendering.

---

## Layout note — moving forward

The same principle should be applied to GHSL (page 3) when you're rested:
- Move scripts that currently live under `Atlas/_scripts/cluster_R/2e_ghsl_discovery/` into `inversion_modules/phase_2_discovery/2e_ghsl/` (your existing folder).
- Then `Atlas/_scripts/` can be deleted entirely. Atlas becomes pure viewer.

GHSL move is its own session — don't bundle it with this θπ work.

---

## What to do if something fails

| Symptom | Likely cause | Fix |
|---|---|---|
| `verify_theta_setup.sh` says "PESTPG_DIR does not exist" | Upstream MODULE_3 hasn't run | Run `02_run_heterozygosity.sh` upstream first |
| `verify_theta_setup.sh` says "pestPG files missing at scale" | `PESTPG_SCALE` mismatch | Check `00_theta_config.sh` — currently expects `win10000.step2000` |
| Launcher says "Missing central config" | `00_inversion_config.sh` not at any of the 3 candidate paths | Either copy it to the new location or `export BASE_DIR=` to override |
| TR_A errors "Could not load any sample's pestPG" | Wrong filename pattern | Check `${SAMPLE}.${PESTPG_SCALE}.pestPG` exists in `PESTPG_DIR` |
| `jq` says `schema_version: 1` | Pre-patch TR_B is in place | Re-extract bundle, verify with `grep "schema_version = 2L" STEP_TR_B_classify_theta.R` |
| `jq: error: theta_pi_per_window has no .values field` | Pre-turn-114 TR_B | Patched TR_B emits both `.values` and the legacy `.samples`-as-objects |

---

## Reference

- **Per-site scaling fix:** ANGSD GitHub issue #329; Korunes & Samuk 2021 (pixy)
- **Background:** `README_theta_pi_scaling.md` in this bundle
- **Atlas patches (turns 114-118):** `HANDOFF_turns_114_to_118_phase2_enrichment_plumbing.md` from earlier session
