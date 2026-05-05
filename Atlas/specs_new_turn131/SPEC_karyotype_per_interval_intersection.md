# SPEC — Karyotype × interval intersection (split GHSL calls per inversion system)

**Status**: drafted turn 130 final session. Quentin's "30-line
post-processing script" from chat `5b793a68`. Not yet implemented
atlas-side. May exist as LANTA-side post-processing — confirm before
re-implementing.

**Trigger** (Quentin, chat `5b793a68`):
> *"The 51 INV/INV calls are from ALL inversions mixed together. Sample
> CGA045 might be INV/INV at the 0.4-6.9 Mb inversion AND INV/nonINV at
> a different inversion at 20 Mb. The current karyotype calling just
> looks at rank position without knowing WHICH inversion system the
> window belongs to."*
>
> *"The fix is simple: split karyotype calls BY triangle interval."*

GHSL produces per-window karyotype calls (`INV_INV` / `INV_nonINV` /
`unclear`) without knowing which inversion system any given window
belongs to. This spec defines the post-processing step that intersects
GHSL calls with triangle/L2 intervals to produce a clean per-sample ×
per-inversion karyotype matrix.

---

## 1. The problem

Current GHSL output:

```
sample_id  start_window  end_window  call         mean_rank  n_windows
CGA045     1             50          INV_INV      0.05       45
CGA045     50            100         INV_nonINV   0.82       38
CGA045     100           150         INV_INV      0.08       42
```

This is **per run** of consecutive windows with consistent karyotype.
Each run might span several different inversion systems mixed together,
or sit entirely inside one system, or partially overlap multiple
systems. The consumer (decomposition, breeding analysis, karyotype
catalogue) needs to know **which inversion system this karyotype call
belongs to**.

## 2. The fix — intersection with triangle/L2 intervals

```python
for sample, run in ghsl_runs:
  for interval in triangle_intervals_per_chrom[run.chrom]:
    overlap = compute_window_overlap(run, interval)
    if overlap >= MIN_OVERLAP_WINDOWS:
      emit_per_sample_per_interval(sample, interval.id, run.call,
                                    n_windows=overlap,
                                    mean_rank=run.mean_rank)
```

Output: per-sample × per-interval karyotype:

```
sample_id  interval_id  interval_type       call         mean_rank  n_windows
CGA045     I1           strong_triangle     INV_INV      0.05       45
CGA045     I2           moderate_triangle   INV_nonINV   0.82       38
CGA045     I3           strong_triangle     INV_INV      0.08       42
CGA102     I1           strong_triangle     INV_nonINV   0.78       45
CGA102     I2           moderate_triangle   INV_INV      0.08       38
```

This is the **per-sample × per-inversion genotype matrix** that:
- Decomposition (`STEP_C01i`) consumes for marker co-segregation
- Breeding analysis consumes for pairing recommendations
- The atlas catalogue uses for per-candidate karyotype counts
- The manuscript bundle uses for Supplementary Table S2

## 3. Conflict-resolution rules

What if a single GHSL run spans **multiple** intervals with different
karyotype implications? Three cases:

### Case A — Run entirely inside one interval
Trivial: assign the run's call to that interval. No conflict.

### Case B — Run spans two adjacent intervals with same call
Most common case. Split the run at the interval boundary, emit two
rows (one per interval), preserving the call.

### Case C — Run spans intervals with different karyotype implications
Should be rare with proper L1/L2 boundaries. If it happens, emit the
run with both interval IDs and a `boundary_crossing: true` flag.
Downstream decides whether to:
- Trust the GHSL call (run consistency suggests no real boundary)
- Trust the L1/L2 boundary (interval definition suggests a structural
  boundary that GHSL missed)

## 4. Edge cases

- **Run too short for an interval**: if the overlap is below `MIN_OVERLAP_WINDOWS`
  (default 5), don't emit a call for that interval. The sample's
  karyotype at that interval is unknown.
- **Sample never appears in any run for an interval**: emit `unknown`.
- **Two GHSL runs from the same sample partially overlap with the same
  interval**: pick the run with higher confidence (mean_rank closer to
  the karyotype centroid).

## 5. Atlas integration

The intersection step runs **upstream** of the atlas (LANTA-side R
script). Atlas consumes the resulting per-sample × per-interval
karyotype TSV directly. Existing infrastructure to surface this:

- The candidate registry (`state.candidateList`) already has per-candidate
  `locked_labels` (one per sample). The intersection's per-interval
  call should populate this.
- The page-3 catalogue's per-candidate karyotype counts are derived
  from `locked_labels`.
- The breeding-readiness card (Atlas 5 Part B) consumes the per-interval
  karyotype.

## 6. Implementation

### LANTA-side: 30-line R script
```r
# split_ghsl_by_interval.R
ghsl_runs        <- read_tsv("ghsl_v6_karyotype_runs.tsv")
triangle_intervals <- read_tsv("triangle_intervals.tsv")

result <- list()
for (i in 1:nrow(ghsl_runs)) {
  run <- ghsl_runs[i, ]
  intervals_in_chrom <- triangle_intervals %>% filter(chrom == run$chrom)
  for (j in 1:nrow(intervals_in_chrom)) {
    interval <- intervals_in_chrom[j, ]
    overlap_windows <- compute_window_overlap(run, interval)
    if (overlap_windows >= MIN_OVERLAP_WINDOWS) {
      result[[length(result) + 1]] <- data.frame(
        sample_id = run$sample_id,
        interval_id = interval$id,
        interval_type = interval$type,
        call = run$call,
        mean_rank = run$mean_rank,
        n_windows = overlap_windows,
        boundary_crossing = !run_entirely_inside(run, interval)
      )
    }
  }
}
write_tsv(bind_rows(result), "per_sample_per_interval_karyotype.tsv")
```

This is **the entire script** Quentin promised in chat `5b793a68`. ~30
lines including I/O.

### Atlas-side: layer detection
Add `per_sample_per_interval_karyotype_v1` to `detectSchemaAndLayers`.
When present, atlas can:
- Cross-reference candidate karyotype calls with this matrix for QC
- Surface "GHSL agrees with PCA bands" as a tier-validation flag
- Color the dosage heatmap rows by GHSL-derived karyotype (an
  alternative to PCA-derived karyotype)

## 7. Open questions

1. **Where does this run?** LANTA-side after GHSL emits its TSV.
   Atlas just consumes the output. No atlas compute.
2. **Should atlas auto-detect agreement/disagreement** between GHSL
   and PCA-derived karyotype, and flag candidates where they diverge?
   Probably yes — high-divergence candidates need manual review.
3. **What about sub-bands (K=6 substructure)?** This spec is K=3 only.
   K=6 substructure is a separate dimension; the per-interval intersection
   doesn't handle it.

## 8. Tests

- Synthetic GHSL runs across two adjacent intervals with same karyotype
  → split into two rows correctly.
- GHSL run wholly inside one interval → single row.
- GHSL run spanning interval boundary with different karyotype on
  either side → boundary_crossing=true, two rows.
- Sub-`MIN_OVERLAP_WINDOWS` overlap → no row emitted.

## 9. Dependencies

- LANTA-side GHSL output (existing per `STEP_C04_snake3_ghsl_v3.R`).
- LANTA-side triangle intervals (existing per `STEP_C01c_triangle_regimes.R`).
- Atlas-side: nothing new; just consumes the resulting TSV via the
  existing layer-detection pattern.

## 10. What this is NOT

- Not a karyotype-calling algorithm. It post-processes already-called
  karyotypes.
- Not a refactor of GHSL. GHSL stays as it is; this script intersects
  its output with another data source.
- Not a sub-cluster detector. It works at K=3 karyotype resolution
  (REF/HET/INV).

## 11. Cross-references

- `SPEC_marker_panel_design_atlas.md` — consumes per-interval karyotype
  for control sample selection.
- `SPEC_per_candidate_breeding_readiness_card.md` — the per-arrangement
  burden table needs per-interval karyotype to compute.
- Chat `5b793a68` — original design conversation.
