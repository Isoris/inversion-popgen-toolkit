# Scrubber JSON pipeline — refresh checkpoint (April 28, 2026)

This is a checkpoint of every R script and JSON-producing pipeline file in the
inversion scrubber stack. It tells you **what version each file is at**, **what
JSON layer it produces**, and **the exact command to regenerate that layer
from a clean precomp .rds**.

If your scrubber shows missing data (PC variance not displaying, ancestry
panel says "data not loaded", relatedness icons absent, etc.), trace it to
the table below — it tells you which exporter you'd need to re-run.

---

## File index — what's in /mnt/user-data/outputs/

| File | Version | Produces | Size |
|------|---------|----------|------|
| `pca_scrubber_v3.html` | **v3.41** | the scrubber itself (single-file webapp) | ~509 KB |
| `STEP_A02_local_pca_windows_by_chr.R` | v9.4 (npc=4) | `<chrom>_local_pca_windows.rds` (raw lostruct windows) | 5.7 KB |
| `STEP_A03_dense_registry_stage1.R` | v9.4 (npc=4) | `<chrom>_dense_registry_stage1.rds` | 9.9 KB |
| `STEP_C01a_precompute.R` | v9.4 (multi-PC) | `<chrom>_precomp.rds` (the precomp `dt` data.table) | 97 KB |
| `export_precomp_to_json_v3.R` | v3 | `<chrom>.json` (the **base** scrubber JSON) | 32 KB |
| `export_ghsl_to_json_v2.R` | v2 (per-K stripes + karyotype runs + heatmap) | `<chrom>_phase4_ghsl.json` | 22 KB |
| `export_ancestry_to_json_v1.R` | v1 (pos_bp contract) | `<chrom>_phase4_ancestry.json` | 21 KB |
| `export_relatedness_to_json_v1.R` | v1 | `<chrom>_relatedness.json` | 15 KB |
| `inspect_ghsl_v6.R` | v6 | (diagnostic tool — prints GHSL summary; doesn't write JSON) | 8 KB |
| `precomp_n_snps_patch.txt` | patch | adds real `n_snps` per window to STEP_C01a + exporter | 6 KB |
| `SCHEMA_V2.md` | v2 | documents the JSON schema for every layer | 25 KB |
| `example_manuscript_bundle.md` | — | example methods-section text bundling all the layers | 9 KB |
| `05_pca_scrubber_contingency.tar.gz` | — | bundled R sources for L3 contingency precomp | 189 KB |

---

## Why "PC variance is missing" — most likely causes

The scrubber's PCA panel shows "PC1 (λ=…)" / "PC2 (λ=…)" axis labels using
`window.lam1` / `window.lam2` from the JSON. If those values are NA or zero
in your JSON, the labels degrade to "PC1" / "PC2" with no λ.

Things to check, in order of likelihood:

1. **Your precomp .rds is from STEP_C01a v9.3 or older** — earlier versions
   removed `lam_1`/`lam_2` computation as "circular" and didn't always restore
   passthrough. **Fix:** rerun STEP_C01a v9.4 (the version in `outputs/`),
   confirming `out_dt` contains `lam_1` and `lam_2` columns before the
   exporter runs.

2. **Your JSON was produced before lam1/lam2 were added to the exporter
   schema.** The current exporter (`export_precomp_to_json_v3.R`) writes
   `lam1`/`lam2` per window. If your old JSON predates this, regenerate.

3. **You're looking at a window in the middle of a coverage gap or N-stretch
   where the precomp couldn't resolve eigenvalues.** This is real and
   correct. Check `n_snps` for that window — if `n_snps < 20`, lam1/lam2
   aren't trustworthy. With the n_snps patch (next bullet), the scrubber
   shows n_snps and you'll see the gap visually.

4. **You haven't applied `precomp_n_snps_patch.txt`.** This patch makes
   STEP_C01a write a real `n_snps` per window to the precomp, and makes the
   exporter pass it through. Without it, n_snps is defaulted to 100 and you
   don't see the gaps that explain missing variance.

---

## Pipeline order — fresh build of one chromosome

The general flow on LANTA at
`/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04`:

```
raw BCF ──► [your existing variant pipeline] ──► chr-split BCFs
                                                       │
                                                       ▼
                                       STEP_A02_local_pca_windows_by_chr.R
                                       (lostruct windows, npc=4)
                                                       │
                                                       ▼
                                              <chrom>_local_pca_windows.rds
                                                       │
                                                       ▼
                                       STEP_A03_dense_registry_stage1.R
                                       (dense overlapping registry)
                                                       │
                                                       ▼
                                              <chrom>_dense_registry_stage1.rds
                                                       │
                                                       ▼
                                       STEP_C01a_precompute.R
                                       (the master precomp; writes dt)
                                                       │
                                                       ▼
                                              <chrom>_precomp.rds
                                                       │
                            ┌──────────────────────────┼─────────────────────────────┐
                            ▼                          ▼                              ▼
            export_precomp_to_json_v3.R    export_ghsl_to_json_v2.R    export_ancestry_to_json_v1.R
            (base scrubber JSON)           (GHSL layer)                  (ancestry layer)
                            │                          │                              │
                            ▼                          ▼                              ▼
                    <chrom>.json          <chrom>_phase4_ghsl.json      <chrom>_phase4_ancestry.json
                                                                                       │
                                                                       export_relatedness_to_json_v1.R
                                                                       (relatedness layer; cohort-level)
                                                                                       │
                                                                                       ▼
                                                                       <chrom>_relatedness.json (single
                                                                       file shared across chroms)
```

You always need at least the base JSON (`<chrom>.json`). The other three
layers are optional — load them via the scrubber's drag-and-drop after the
base JSON has loaded.

---

## Step-by-step rebuild commands

All commands assume the working directory is your project folder on LANTA
and the R environment has lostruct, data.table, jsonlite installed.

### Step 1 — STEP_A02: lostruct windows

Required inputs:
- chromosome BCF (e.g. `data/bcf/LG28.bcf`)
- sample whitelist (TSV with one CGA name per line)

```bash
Rscript STEP_A02_local_pca_windows_by_chr.R \
    --bcf data/bcf/LG28.bcf \
    --samples data/samples_226.tsv \
    --chrom LG28 \
    --window_snps 100 \
    --npc 4 \
    --out work/LG28_local_pca_windows.rds
```

**What you'll see:** progress logs like `[STEP_A02] window 1000 / 4302 …`.
At completion, the .rds contains a list with `$pcs` (windows × samples ×
npc) and `$snps` (windows × snp_count).

### Step 2 — STEP_A03: dense registry

```bash
Rscript STEP_A03_dense_registry_stage1.R \
    --windows work/LG28_local_pca_windows.rds \
    --chrom LG28 \
    --npc 4 \
    --out work/LG28_dense_registry_stage1.rds
```

The dense registry is the index lookup that lets STEP_C01a access any
window quickly during precompute.

### Step 3 — STEP_C01a: precompute

This is the heavy step. On the 226-sample cohort with npc=4 it takes
30-60 minutes per chromosome on LANTA's normal queue.

```bash
Rscript STEP_C01a_precompute.R \
    --registry work/LG28_dense_registry_stage1.rds \
    --windows work/LG28_local_pca_windows.rds \
    --chrom LG28 \
    --npc 4 \
    --out work/LG28_precomp.rds
```

The output is a list `prec` with:
- `prec$dt` — the per-window data.table with one row per window. Critical
  columns: `window_id`, `start_bp`, `end_bp`, `center_mb`, `z_pc1` (or
  `robust_z_pc1`), `lam_1`, `lam_2`, `pc1_mat` (or per-sample columns),
  `pc2_vec`, optionally `n_snps`.
- `prec$samples` — the sample table aligned to the columns of `pc1_mat`.
- `prec$meta` — chromosome metadata (chrom, length_bp, n_windows, etc.)

**Sanity check before proceeding:**

```r
prec <- readRDS("work/LG28_precomp.rds")
str(prec$dt[1:3])         # one row per window
table(is.na(prec$dt$lam_1))   # should be 0
range(prec$dt$lam_1)          # should be > 0
range(prec$dt$z_pc1)          # absolute values up to ~5-7 typical
```

If `lam_1` is all NA, your STEP_C01a is too old — pull the v9.4 in
`outputs/` and rerun. If `z_pc1` doesn't exist, check for `robust_z_pc1`
or `lam_1_z` (the exporter accepts any of these names).

### Step 4 — base scrubber JSON

Required inputs:
- the precomp from step 3
- the L1/L2 envelope tables from your envelope-discovery pipeline

```bash
Rscript export_precomp_to_json_v3.R \
    --precomp work/LG28_precomp.rds \
    --l1_envelopes data/envelopes/LG28_d17L1_envelopes.tsv \
    --l2_envelopes data/envelopes/LG28_d17L2_envelopes.tsv \
    --l2_boundaries data/envelopes/LG28_d17L2_boundaries.tsv \
    --out json/LG28.json
```

**What gets written:**
- per-window: `start_bp`, `end_bp`, `center_mb`, `z`, **`lam1`, `lam2`**
  (← these drive the variance display you noticed missing), `pc1[]`, `pc2[]`,
  optional `theta[]`
- per-envelope: candidate id, window range, mean similarity
- per-sample: `cga`, `ind`, `family_id`, `ancestry`
- multi-scale sim_mat thumbnails at each smoothing scale

**Sanity check the JSON:**

```bash
python3 -c "
import json
with open('json/LG28.json') as f: d = json.load(f)
print('chrom:', d['chrom'])
print('n_windows:', d['n_windows'])
print('window 0:', list(d['windows'][0].keys()))
print('lam1 range:', min(w['lam1'] for w in d['windows']),
                    max(w['lam1'] for w in d['windows']))
print('n_snps present:', 'n_snps' in d['windows'][0])
"
```

Expected output:
```
chrom: LG28
n_windows: 4302
window 0: ['start_bp', 'end_bp', 'center_mb', 'z', 'lam1', 'lam2', 'pc1', 'pc2']
lam1 range: 0.0xxx 0.0xxx
n_snps present: False    (← True only if you applied the n_snps patch)
```

### Step 5 (optional) — apply n_snps patch

If `n_snps present: False` and you want real per-window SNP counts in the
windows summary table, apply `precomp_n_snps_patch.txt`:

```bash
# Read the patch file. It contains two diff blocks:
#   1. a small change to STEP_C01a_precompute.R to record n_snps per window
#   2. a small change to export_precomp_to_json_v3.R to pass n_snps through
# Apply manually (the patches are short and clearly commented).
cat precomp_n_snps_patch.txt
# After applying, rerun STEP_C01a (heavy) and the exporter (fast).
# OR if you only want to update the existing precomp's exporter without
# rerunning C01a, just apply patch block 2 and rerun the exporter — it
# will default n_snps=100 per window with a warning, which the scrubber
# annotates as "(N×100 default)".
```

### Step 6 (optional) — GHSL layer

GHSL (Group-level Haplotype Similarity by Locus) provides per-K stripes and
karyotype-call-by-window data. The export script reads from your existing
GHSL outputs in the `inversion_codebase_v8.5/snake_*` results.

```bash
Rscript export_ghsl_to_json_v2.R \
    --ghsl_dir results/snake3_ghsl_LG28/ \
    --chrom LG28 \
    --precomp work/LG28_precomp.rds \
    --out json/LG28_phase4_ghsl.json
```

The v2 exporter writes:
- `panel`: top-level GHSL panel (per-window scalar)
- `per_K_stripes`: separate panels for K=3, K=6
- `per_sample_karyotype_runs`: contiguous karyotype call segments per sample
- `heatmap`: full GHSL similarity heatmap data

Drop the resulting `<chrom>_phase4_ghsl.json` onto the scrubber after the
base JSON loads — additional tracks appear in the lines panel.

### Step 7 (optional) — ancestry layer

If you have NGSadmix or PCAngsd runs:

```bash
Rscript export_ancestry_to_json_v1.R \
    --ngsadmix_qopt results/ngsadmix/LG28_K6.qopt \
    --ngsadmix_log  results/ngsadmix/LG28_K6.log \
    --samples data/samples_226.tsv \
    --chrom LG28 \
    --window_size_bp 50000 \
    --out json/LG28_phase4_ancestry.json
```

This populates page 7 (ancestry) with admixture-along-chromosome heatmap +
per-sample dot panel + per-window Q-mean tracks.

### Step 8 (optional) — relatedness layer

This is **cohort-level** — one file shared across all chromosomes. Run once.

```bash
Rscript export_relatedness_to_json_v1.R \
    --ngsrelate results/ngsrelate/catfish_226_relatedness.res \
    --samples data/samples_226.tsv \
    --thresh_1st 0.177 \
    --thresh_2nd 0.0884 \
    --thresh_3rd 0.0442 \
    --out json/catfish_226_relatedness.json
```

The exporter expects ngsrelate output with theta in **column 18** (the
historical bug in older ngsrelate versions placed theta in column 5 — the
exporter has a fallback that detects this and warns).

Output JSON has `hub_id_1st`, `hub_id_2nd`, `hub_id_3rd` per sample (which
relatedness hub they belong to at each kinship threshold).

Drop the relatedness.json onto the scrubber once and it persists across
chromosome switches in the same session.

---

## Common failure modes

### "PC1 (λ=)" with empty λ value
- Your JSON's `windows[i].lam1` is null/missing. Regenerate base JSON (Step 4)
  with the current exporter against a precomp that has `lam_1` columns.

### Z panel shows nothing but blank
- Your JSON's `windows[i].z` is null. Either your precomp has no z column,
  or the exporter picked the wrong column name. The exporter accepts:
  `z_pc1`, `robust_z_pc1`, `lam_1_z`, `z_lam1`. Check which one your precomp
  actually uses and confirm the exporter found it (look at the exporter's
  log line that says "[exporter] Using z column: <name>").

### Lines panel shows nothing
- Your JSON's `windows[i].pc1` arrays are empty or have wrong length.
  Confirm `prec$dt$pc1_mat` has shape (n_windows × n_samples) and the
  exporter is iterating correctly.

### Ancestry tab says "Load a precomp JSON, then drop in <chrom>_phase4_ancestry.json"
- You haven't loaded the ancestry layer. This is normal — the scrubber's
  base JSON doesn't contain ancestry. Run Step 7 to produce
  `<chrom>_phase4_ancestry.json` and drop it onto the scrubber.

### Relatedness icons all grey
- Either you haven't loaded `catfish_226_relatedness.json` (run Step 8), or
  your sample CGA names in the relatedness data don't match the sample
  CGAs in the precomp. The exporter has a name-canonicalization step but
  if your sample TSV uses different name conventions, you may need to
  patch the matcher.

---

## Quick-restart cheat sheet

If you just want a single chromosome scrubbing session up and running, the
minimum pipeline is:

```bash
# Assumes you have STEP_A02 outputs in work/ already
Rscript STEP_C01a_precompute.R --chrom LG28 --npc 4 \
        --registry work/LG28_dense_registry_stage1.rds \
        --windows  work/LG28_local_pca_windows.rds \
        --out      work/LG28_precomp.rds

Rscript export_precomp_to_json_v3.R \
        --precomp work/LG28_precomp.rds \
        --l1_envelopes data/envelopes/LG28_d17L1_envelopes.tsv \
        --l2_envelopes data/envelopes/LG28_d17L2_envelopes.tsv \
        --l2_boundaries data/envelopes/LG28_d17L2_boundaries.tsv \
        --out json/LG28.json

# Open pca_scrubber_v3.html in a browser, drop json/LG28.json onto it.
```

For a full multi-layer session add Steps 6, 7, 8 once each.
