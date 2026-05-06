# PATCH — add `lam_top_k` to per-window precomp JSON (schema 2.16)

**Target file**: `export_precomp_to_json_v3.R` on LANTA at
`/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/`
(or wherever your current copy lives; this is the script the empty-state
help text on the atlas references).

**Atlas-side counterpart**: turn 120, `_renderScreeInsetHTML` reads
`w.lam_top_k` (or falls back to `[w.lam1, w.lam2]` with a hint when only
the legacy fields are loaded). So you can ship this patch on whatever
schedule fits — atlas already works with both old and new precomp.

---

## What changes

Per window, the JSON gains one new field:

```json
"lam_top_k": [0.82, 0.31, 0.18, 0.12, 0.08]
```

Sorted descending, length 5 by default (configurable). `lam_top_k[0]`
must equal `lam1` and `lam_top_k[1]` must equal `lam2` (atlas asserts
this on load).

`meta$schema_version` bumps to `"2.16"`.

---

## The change (unified diff style)

The exact line numbers depend on your current copy. Search for the
per-window emit loop (where `lam1` and `lam2` are written today). The
loop body looks roughly like:

```r
windows[[i]] <- list(
  idx       = i - 1L,
  start     = win_starts[i],
  end       = win_ends[i],
  center_mb = (win_starts[i] + win_ends[i]) / 2 / 1e6,
  z         = z_scores[i],
  lam1      = local_pca[[i]]$eigenvalues[1],
  lam2      = local_pca[[i]]$eigenvalues[2],
  ...
)
```

Add ONE line:

```r
windows[[i]] <- list(
  idx       = i - 1L,
  start     = win_starts[i],
  end       = win_ends[i],
  center_mb = (win_starts[i] + win_ends[i]) / 2 / 1e6,
  z         = z_scores[i],
  lam1      = local_pca[[i]]$eigenvalues[1],
  lam2      = local_pca[[i]]$eigenvalues[2],
+ lam_top_k = head(local_pca[[i]]$eigenvalues, 5),
  ...
)
```

Conditions:

1. **`local_pca[[i]]$eigenvalues` must already contain ≥ 5 eigenvalues.**
   If your local-PCA worker uses `RSpectra::eigs_sym(..., k = 2)` for
   speed (only top 2), you need to bump `k = 5` first. Check the worker
   call site.

2. **If `local_pca[[i]]$eigenvalues` doesn't exist** (older pipelines
   stored only `lam1`/`lam2` as scalars), wrap the new line in a guard:

   ```r
   lam_top_k = if (length(local_pca[[i]]$eigenvalues) >= 5)
                 head(local_pca[[i]]$eigenvalues, 5) else NULL
   ```

   `jsonlite::toJSON(NULL, null = "null")` writes JSON `null`, which the
   atlas treats as "fall back to [lam1, lam2]" — same behavior as the
   field being absent entirely.

3. **Bump the schema version**:

   ```r
   meta$schema_version <- "2.16"
   ```

   Search for the existing `schema_version <- "2.15"` (or whatever
   you're at today) and bump.

---

## If the pipeline only computes top-2 eigenvalues today

This is the harder case. Your local-PCA worker likely lives in a
function called something like `compute_local_pca_window()` or
`window_pca()` — search for `RSpectra` or `eigs_sym` or `eigen` in your
codebase. The fix:

```r
# OLD:
ev <- RSpectra::eigs_sym(window_cov, k = 2, which = "LM")
list(eigenvalues = ev$values, pc1 = ev$vectors[, 1], pc2 = ev$vectors[, 2])

# NEW:
ev <- RSpectra::eigs_sym(window_cov, k = 5, which = "LM")
list(eigenvalues = ev$values, pc1 = ev$vectors[, 1], pc2 = ev$vectors[, 2])
```

Per-window slowdown: ~2× (going from 2 to 5 eigenvalues on
n_samples × n_samples matrices is cheap). Total precomp runtime impact:
likely under a minute even on full-genome runs.

If your matrix is small (n_samples ≤ 30, like the F1 hybrid cohort) you
can use the full `eigen()` instead and not worry about k:

```r
ev <- eigen(window_cov, symmetric = TRUE, only.values = FALSE)
list(eigenvalues = ev$values[1:5],
     pc1 = ev$vectors[, 1],
     pc2 = ev$vectors[, 2])
```

For the 226-sample cohort, stick with `RSpectra::eigs_sym(k = 5)` — full
`eigen()` on a 226×226 matrix per window × thousands of windows is
slower than necessary.

---

## Validation — re-run on one chrom, verify the inset

After patching, regenerate ONE chromosome's precomp:

```bash
Rscript export_precomp_to_json_v3.R \
  --precomp /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/precomp/LG28.rds \
  --l1_envelopes ..._d17L1_envelopes.tsv \
  --l2_envelopes ..._d17L2_envelopes.tsv \
  --l2_boundaries ..._d17L2_boundaries.tsv \
  --out LG28.json
```

Then load `LG28.json` in the atlas, enable the **scree** checkbox in
the tracked-samples panel (right side of the PCA scatter), and scrub
through windows. You should see 5 bars (no fallback hint). If you still
see 2 bars + the "k≥3 needs precomp ≥2.16" hint, the JSON didn't pick
up the new field — check that `meta$schema_version == "2.16"` and that
`windows[[1]]$lam_top_k` is a 5-element array (`jq '.windows[0].lam_top_k' LG28.json`).

---

## Output size impact

Per window: 5 doubles × ~8 bytes = ~40 extra bytes (JSON encoding adds
some). Per chrom (LG28 ~3000 windows): ~120 KB extra on a ~40 MB
precomp. **0.3% size increase. Negligible.**

If you go to k = 7 instead of 5 for richer scree: ~170 KB extra
(0.4%). Still negligible.

---

## Backward compatibility

The atlas already handles missing `lam_top_k` gracefully (legacy 2-bar
inset + hint). So:

- Old precomps continue to work.
- Mixed sessions (some chroms loaded with 2.16, some with 2.15) work
  per-window: each window's inset uses whatever fields that window has.
- No need to regenerate everything at once. Patch, regenerate one chrom
  for testing, regenerate the rest at your convenience.

---

## SCHEMA_V2.md doc snippet

Add to `SCHEMA_V2.md` under §windows or wherever per-window field docs
live:

> ### `lam_top_k` *(schema 2.16, optional)*
>
> **Type**: array of doubles, length 5 (default) or up to 8.
> **Sort**: descending.
> **Constraint**: when present, `lam_top_k[0] === lam1` and
> `lam_top_k[1] === lam2`.
>
> Top-k eigenvalues from the per-window local PCA covariance
> decomposition. Used by the atlas's optional scree-plot inset on the
> tracked-sample plane (page 1) to surface the elbow shape — "is this
> window's variance dominated by PC1, or smeared across many PCs?"
>
> Backward compat: omitting the field is fine; the atlas falls back to
> a 2-bar inset built from `lam1` and `lam2` with a hint that schema
> ≥2.16 is needed for k ≥ 3.
