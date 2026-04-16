# Cheat 17 fossil detection — archived 2026-04-17

Removed from `phase_4_postprocessing/4a_existence_layers/STEP_C01g_boundary_catalog_*.R`
on 2026-04-17 after finding that the feature was structurally unreachable
in single-pass runs.

## Why archived

**Circular dependency.** The removed block reads C01d's
`candidate_scores.tsv.gz` via `--scores` to find boundaries lying
*outside* any candidate (the "orphaned" boundaries that become
fossil candidates). But C01g runs *before* C01d in the canonical
flow because C01d reads C01g's `boundary_catalog_unified.tsv.gz` via
`--boundary_dir` for D11 scoring. The two scripts therefore require
a two-pass orchestration:

    pass 1: C01g (no --scores) → boundary_catalog_unified_pass1.tsv.gz
    pass 2: C01d --boundary_dir=pass1 → candidate_scores.tsv.gz
    pass 3: C01g --scores=pass2      → boundary_catalog_unified_pass2.tsv.gz
    pass 4: (optional) C01d again with the updated boundary catalog

The orchestrator never wired this up. In single-pass runs the
`--scores` flag was always an empty path, the entire fossil block
fell through to the `else` branch at former line 1263-1264
("No --scores file — Cheat 17 fossil detection skipped"), and
`boundary_activity` defaulted to `"active"` for every boundary
regardless.

**Not in the manuscript main path.** The fossil-vs-active biological
distinction was a diagnostic annotation ("here are structurally
plausible boundaries with no currently-segregating inversion behind
them"), not a main-path deliverable. Removing it is simpler than
orchestrating the two passes.

## What was removed

- `--scores` CLI flag (now accepted silently for back-compat)
- The entire Cheat 17 block that populated `cheat17_class` and
  `cheat17_inv_likeness`
- The conditional logic inside `boundary_activity` that mapped
  `cheat17_class` values to activity states

## What was kept

- The two columns `cheat17_class` and `cheat17_inv_likeness` are
  still emitted in `boundary_catalog_unified.tsv.gz`, as all-NA, so
  downstream readers expecting them don't break
- `boundary_activity` column still emitted, now unconditionally
  `"active"` on every row

## Restoring the feature

If the fossil distinction becomes needed again, the cleanest way is
a **post-hoc annotation script** that reads the already-written
`boundary_catalog_unified.tsv.gz` + `candidate_scores.tsv.gz` and
writes a separate `boundary_activity.tsv.gz`. No need to re-enter
the C01g main flow.

The original code block is below, preserved verbatim from the
2026-04-17 working copy.

---

## Original code (verbatim)

```r
# =============================================================================
# CHEAT 17: FOSSIL BREAKPOINT CLASSIFICATION
# =============================================================================
# Boundaries that look structural in sim_mat but have NO active inversion:
#   - no trimodal PCA (no 3-class grouping)
#   - no SV calls at the position
#   - inv_likeness may be elevated (residual signal from old event)
# These are relics of inversions that broke apart. Diagnostic annotation.

message("\n[C01g] Running Cheat 17 (fossil breakpoints)...")

# Load candidate regions if available (to identify orphaned boundaries)
cand_regions <- NULL
if (!is.null(scores_file) && file.exists(scores_file)) {
  cand_regions <- tryCatch(fread(scores_file), error = function(e) NULL)
  if (!is.null(cand_regions)) {
    if ("start_mb" %in% names(cand_regions) && !"start_bp" %in% names(cand_regions))
      cand_regions[, start_bp := start_mb * 1e6]
    if ("end_mb" %in% names(cand_regions) && !"end_bp" %in% names(cand_regions))
      cand_regions[, end_bp := end_mb * 1e6]
    if (!"chr" %in% names(cand_regions) && "chrom" %in% names(cand_regions))
      cand_regions[, chr := chrom]
    message("[C01g] Loaded ", nrow(cand_regions), " candidate regions for fossil check")
  }
}

deduped[, cheat17_class := NA_character_]
deduped[, cheat17_inv_likeness := NA_real_]

if (!is.null(cand_regions) && nrow(cand_regions) > 0) {
  for (chr in unique(deduped$chrom)) {
    pc <- precomp_list[[chr]]
    if (is.null(pc)) next
    dt_chr <- pc$dt
    chr_idx <- which(deduped$chrom == chr)
    chr_cands <- cand_regions[chr == (..chr) | chrom == (..chr)]
    if (length(chr_idx) == 0) next

    for (bi in chr_idx) {
      bp <- deduped$boundary_bp[bi]
      inside_candidate <- FALSE
      if (nrow(chr_cands) > 0) {
        inside_candidate <- any(chr_cands$start_bp <= bp & chr_cands$end_bp >= bp)
      }
      if (inside_candidate) next  # not orphaned

      mid_bps <- (dt_chr$start_bp + dt_chr$end_bp) / 2
      near_w <- which(abs(mid_bps - bp) < 100000)
      if (length(near_w) < 3) next

      near_dt <- dt_chr[near_w]

      il <- if ("inv_likeness" %in% names(near_dt))
        mean(near_dt$inv_likeness, na.rm = TRUE) else NA_real_
      il_elevated <- !is.na(il) && il > 0.3

      trimodal <- FALSE
      if ("inv_dip_p" %in% names(near_dt)) {
        trimodal <- any(near_dt$inv_dip_p < 0.05, na.rm = TRUE)
      }

      sv_overlap <- FALSE
      if ("sv_inv_overlap" %in% names(near_dt)) {
        sv_overlap <- any(near_dt$sv_inv_overlap > 0, na.rm = TRUE)
      }

      classification <- if (il_elevated && trimodal && sv_overlap) "MISSED_INVERSION"
        else if (il_elevated && !trimodal && !sv_overlap) "FOSSIL_CANDIDATE"
        else if (!il_elevated) "FALSE_BOUNDARY"
        else "AMBIGUOUS"

      deduped[bi, cheat17_class := classification]
      deduped[bi, cheat17_inv_likeness := round(il, 4)]
    }
  }
}
```
