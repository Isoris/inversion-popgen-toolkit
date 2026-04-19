# INSTALL — setting up MODULE_QC_ShelfDiagnosis from scratch

One-page walkthrough for a clean install on LANTA. Assumes `${BASE}` is the
project root at `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/`.

## 1. Place the module in your repo

If you have a local clone of the toolkit:

```bash
cd ${BASE}/inversion-popgen-toolkit            # your existing repo clone
# Extract the module tarball or rsync from your laptop:
tar xzf /path/to/MODULE_QC_ShelfDiagnosis.tar.gz -C Modules/
# Result: Modules/MODULE_QC_ShelfDiagnosis/

git add Modules/MODULE_QC_ShelfDiagnosis
git commit -m "Add MODULE_QC_ShelfDiagnosis"
git push
```

If you're cloning fresh on a new HPC node:

```bash
cd ${BASE}
git clone https://github.com/Isoris/inversion-popgen-toolkit.git
cd inversion-popgen-toolkit/Modules/MODULE_QC_ShelfDiagnosis
```

## 2. Run the installer

```bash
cd ${BASE}/inversion-popgen-toolkit/Modules/MODULE_QC_ShelfDiagnosis
bash install.sh
```

What it does:
- Verifies that `BEAGLE_DIR`, `PRECOMP_DIR`, `HET_DIR`, etc. actually exist.
- Compiles Engine B (`instant_q.cpp`) if source exists but the binary doesn't.
- Creates `config.local.sh` — a gitignored overlay where you can override any
  path without touching the committed `00_config.sh`.
- Makes `results/tracks`, `results/figures`, `results/logs` directories.

Flags:
- `bash install.sh --check-only` — validate only, don't build or write config
- `bash install.sh --compile-only` — recompile `instant_q` only

## 3. Edit config.local.sh if anything is off

```bash
vim config.local.sh
```

Uncomment and edit the variables that don't match your layout. Common ones:

```bash
export BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
export BEAGLE_DIR="${BASE}/inversion_localpca_v7/02_snps_beagle"
export PRECOMP_DIR="${BASE}/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/precomp"
export HET_DIR="${BASE}/het_roh"
export UNIFIED_ANCESTRY_DIR="${BASE}/unified_ancestry"
export LOCAL_Q_DIR="${UNIFIED_ANCESTRY_DIR}/local_Q"
```

## 4. Precompute Engine B (only if you don't have local_Q/ caches yet)

Engine B = `instant_q` C++ binary. After `install.sh` compiles it, build the
caches for all 28 chromosomes:

```bash
cd ${BASE}/unified_ancestry
sbatch launchers/LAUNCH_instant_q_precompute.slurm
```

This writes `${LOCAL_Q_DIR}/<CHR>.local_Q_summary.tsv.gz` and
`<CHR>.local_Q_samples.tsv.gz` for every chromosome. Takes ~10 minutes on a
LANTA compute node.

### Multi-scale precompute (1x, 5x, 10x window sizes)

`STEP_Q06_multiscale.sh` expects caches at three scales:

```
${LOCAL_Q_DIR}/scale_1x/<CHR>.local_Q_summary.tsv.gz
${LOCAL_Q_DIR}/scale_5x/<CHR>.local_Q_summary.tsv.gz
${LOCAL_Q_DIR}/scale_10x/<CHR>.local_Q_summary.tsv.gz
```

To produce them, rerun the precompute with different `WIN_SIZE` values and
move the outputs into the scale subdirectories. Example:

```bash
# After 1x precompute:
mkdir -p ${LOCAL_Q_DIR}/scale_1x
mv ${LOCAL_Q_DIR}/*.local_Q_summary.tsv.gz ${LOCAL_Q_DIR}/scale_1x/
mv ${LOCAL_Q_DIR}/*.local_Q_samples.tsv.gz ${LOCAL_Q_DIR}/scale_1x/

# 5x: resubmit with 5× the WIN_SIZE
WIN_SIZE=500 sbatch launchers/LAUNCH_instant_q_precompute.slurm
mkdir -p ${LOCAL_Q_DIR}/scale_5x
mv ${LOCAL_Q_DIR}/*.local_Q_summary.tsv.gz ${LOCAL_Q_DIR}/scale_5x/
# ... and so on for 10x
```

Backward compatibility: if you only have the flat 1x cache, `STEP_Q06_multiscale.sh`
will use it as scale_1x automatically — no need to move files.

## 5. Quick test on LG28

```bash
cd ${BASE}/inversion-popgen-toolkit/Modules/MODULE_QC_ShelfDiagnosis
SHELF_START_MB=15 SHELF_END_MB=18 bash run_all.sh C_gar_LG28
```

Expected wall time: 2–5 minutes. Outputs:
- `results/tracks/*.C_gar_LG28.*` (TSV per Q-step)
- `results/figures/diagnostic.C_gar_LG28.pdf` (8-track composite)
- Console summary printed at the end with shelf-vs-ref ratios

## 6. Full 28-chrom run via SLURM

```bash
sbatch slurm/array_28chrom.sh
# If you want the shelf highlighted on all 28 (unusual):
sbatch --export=ALL,SHELF_START_MB=15,SHELF_END_MB=18 slurm/array_28chrom.sh
```

## 7. Build the scrubber JSON for interactive use

After `run_all.sh` finishes on a chromosome:

```bash
Rscript export_precomp_to_json.R \
  --precomp  ${PRECOMP_DIR}/C_gar_LG28.precomp.rds \
  --theta    results/tracks/theta.C_gar_LG28.win50000.step10000.tsv.gz \
  --ancestry results/tracks/ancestry_window.C_gar_LG28.tsv \
  --ancestry_samples results/tracks/ancestry_sample.C_gar_LG28.maxQ_wide.tsv.gz \
  --out      LG28_scrubber.json
```

Download the JSON + the `pca_scrubber.html` file (kept at the repo root, not
inside this module — since it's reusable across modules). Open the HTML
locally, load the JSON.

## 8. Troubleshooting

**"No pos file" for a chromosome**: Q01 needs `main_qcpass.<CHR>.pos.fixed`.
Check `BEAGLE_DIR`. The naming convention can differ — inspect:
`ls ${BEAGLE_DIR} | head`

**Engine B cache missing**: Q06 is skipped with a warning. Fix by running the
precompute in step 4. Not fatal — the rest of the pipeline still runs.

**mosdepth output not found**: Q03 drops the coverage track. Run mosdepth
yourself or set `RUN_MOSDEPTH=1` in env to have Q03 invoke it.

**"pestPG header mismatch"**: Q05 expects 14-column ANGSD thetaStat output.
If you regenerated with a different ANGSD version, the column count may
differ. Check:  `zcat file.pestPG | head -1 | awk '{print NF}'`
