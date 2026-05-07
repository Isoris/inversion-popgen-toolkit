# Phase-2 enrichment bundle — θπ (page 12) + GHSL (page 3)

This bundle delivers both phase-2 discovery streams that feed the atlas's
enrichment pages, organized under `inversion_modules/phase_2_discovery/`.

```
inversion_modules/phase_2_discovery/
├── 2e_ghsl/                        ← page 3 (GHSL haplotype contrast)
│   ├── README.md
│   ├── STEP_C04_snake3_ghsl_v6.R                  Heavy compute (~1 hr/chrom)
│   ├── STEP_C04b_snake3_ghsl_classify.R           Light classifier (~30 s)
│   ├── STEP_C04c_ghsl_local_pca.R                 Page-3 local PCA enrichment
│   ├── STEP_C04d_ghsl_d17_wrapper.R               Page-3 D17 boundary detector
│   ├── export_ghsl_to_json_v3.R                   Page-3 JSON consolidator
│   ├── LAUNCH_STEP_C04_ghsl_v6_compute.slurm      Stage 1 launcher
│   ├── LAUNCH_STEP_C04b_ghsl_v6_classify.slurm    Stage 2 launcher
│   └── LAUNCH_STEP_C04cd_ghsl_enrichment.slurm    Stage 3 launcher (NEW)
│
└── 2f_theta_discovery/             ← page 12 (θπ diversity)
    ├── README_theta_pi_scaling.md
    ├── STEP_TR_A_compute_theta_matrices.R         Heavy: pestPG → matrix
    ├── STEP_TR_B_classify_theta.R                 Classifier + JSON emitter
    ├── STEP_TR_C_theta_d17_wrapper.R              (optional) D17 envelope detector
    ├── STEP_TR_D_augment_theta_json.R             (optional) D17 → JSON augmenter
    ├── 00_sanity_check_pestPG_scaling.sh
    ├── 00_theta_config.sh
    ├── LAUNCH_TR_theta_pi.slurm                   28-chrom array launcher
    ├── chrom.list
    └── verify_theta_setup.sh                      7-check preflight
```

**No `Atlas/_scripts/`.** Atlas is the viewer. Pipelines live under `phase_2_discovery/`.

---

## Push to LANTA

```bash
# from your laptop
tar -czf bundle.tar.gz -C combined_bundle .
scp bundle.tar.gz lanta:~/

# on LANTA
ssh lanta
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit
tar -xzf ~/bundle.tar.gz

# Two folders land in inversion_modules/phase_2_discovery/
ls inversion_modules/phase_2_discovery/2e_ghsl/        # 9 files
ls inversion_modules/phase_2_discovery/2f_theta_discovery/   # 10 files
```

---

## Order of operations

### θπ (page 12) — easier, do this first

Single launcher does TR_A → TR_B per chrom:

```bash
cd inversion_modules/phase_2_discovery/2f_theta_discovery
mamba activate assembly
bash verify_theta_setup.sh                           # preflight, ~30 s
# If STATUS: READY ✓ →
sbatch --array=1-28 LAUNCH_TR_theta_pi.slurm chrom.list
```

Per-chrom walltime: 5–15 min. Outputs at `${JSON_OUT_DIR}/<chr>/<chr>_phase2_theta.json`.

For details: `2f_theta_discovery/README_theta_pi_scaling.md` and the standalone `RUN_ON_LANTA_theta_pi.md` from the prior bundle.

### GHSL (page 3) — three stages

```bash
cd inversion_modules/phase_2_discovery/2e_ghsl
mamba activate assembly

# Stage 1: heavy compute (~1 hr/chrom, all 28 in parallel = ~1 hr wall)
sbatch LAUNCH_STEP_C04_ghsl_v6_compute.slurm

# Stage 2: light classifier (~30 s/chrom, single job loops 28 chroms)
# Wait for Stage 1 first (or use --dependency=afterok)
sbatch LAUNCH_STEP_C04b_ghsl_v6_classify.slurm

# Stage 3: page-3 enrichment (C04c → C04d → export per chrom)
# Wait for Stage 2 first
sbatch LAUNCH_STEP_C04cd_ghsl_enrichment.slurm
```

Stage 3 outputs at `${GHSL_DIR}/json_out/<chr>/<chr>_phase2_ghsl.json`.

---

## What was renamed / refactored from the old bundle

The old bundle had two GHSL folders:
- `2e_ghsl/` — Family A: STEP_C04 + STEP_C04b + their LAUNCH_* scripts
- `2e_ghsl_discovery/` — Family B: STEP_C04c + STEP_C04d + export_v3 + `run_ghsl_array.sh`

This bundle merges them into one folder (`2e_ghsl/`) per your earlier choice. The only file with substantive changes:

| Old name | New name | What changed |
|---|---|---|
| `run_ghsl_array.sh` (Family B) | `LAUNCH_STEP_C04cd_ghsl_enrichment.slurm` | Now uses canonical SLURM-launcher idiom matching `LAUNCH_STEP_C04*` siblings. Sources `inversion_modules/00_inversion_config.sh` (with legacy fallback). All paths come from config (`GHSL_DIR`, `BASE`, `RSCRIPT_BIN`) — no hardcoded `inversion_codebase_v8.5/phase_2_discovery/2e_ghsl_discovery/` paths. Renamed for consistency with C04 / C04b SLURM scripts. |

All other R scripts (`STEP_C04*.R`, `export_ghsl_to_json_v3.R`) are byte-identical to what you uploaded — they were already fine.

The README in `2e_ghsl/` got a new section ("Page-3 atlas enrichment (C04c / C04d / export_v3)") explaining the three additional scripts and the run order.

---

## Config requirements

Both pipelines source `00_inversion_config.sh` (or the θπ one sources `00_theta_config.sh` after it). They expect these env vars to be exported:

| Var | Used by | Purpose |
|---|---|---|
| `BASE` | both | repo root (`/scratch/.../`) |
| `RSCRIPT_BIN` | both | path to Rscript in the assembly env |
| `MDS_DIR` | GHSL C04 | upstream precomp dir |
| `GHSL_DIR` | GHSL all | output root for GHSL artifacts |
| `INVDIR` | both | `${BASE}/inversion-popgen-toolkit/inversion_modules` |
| `PESTPG_DIR` | θπ | upstream pestPG files |
| `JSON_OUT_DIR` | θπ | output dir for `<chr>_phase2_theta.json` |

The θπ launcher tries the new path first (`inversion-popgen-toolkit/inversion_modules/00_inversion_config.sh`), falls back to the legacy locations. The GHSL launcher does the same.

If a launcher fails with "Missing config", set `CONFIG=...` to override.

---

## After both run — feeding the atlas + STEP_T05

Once you have all 28 θπ JSONs and 28 GHSL JSONs:

```bash
# Atlas (in browser)
# Drag-drop dosage precomp JSON first (page 1), then:
#   <chr>_phase2_theta.json   → page 12 indicators flip 🟢
#   <chr>_phase2_ghsl.json    → page 3 indicators flip 🟢

# STEP_T05 CUSUM on the θπ JSON (per-sample observations for boundary refinement)
$RSCRIPT /path/to/STEP_T05_theta_cusum.R \
    --json ${JSON_OUT_DIR}/C_gar_LG28/C_gar_LG28_phase2_theta.json \
    --out-dir /tmp/cusum_LG28/ \
    --mode whole_chrom \
    --lib /path/to/lib_persample_cusum.R
```

The matching `STEP_R02_ghsl_cusum.R` driver (mirror of T05 for GHSL matrices) is **not yet built** — it's deferred per your earlier note. When you're ready, see the handoff doc from earlier in this session for the architecture.

---

## Verification before submitting

For θπ: `bash verify_theta_setup.sh` (in the 2f folder).

For GHSL: no preflight script yet — but the LAUNCH scripts each do upstream-RDS existence checks and fail early with clear error messages. Stage 3 in particular checks for all 4 RDSes from C04 + the D17 scripts being available. If anything's missing, you'll get a clear error before any compute runs.
