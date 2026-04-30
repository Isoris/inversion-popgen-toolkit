# Handoff — Regime change / band-switch event detector

**Status**: ⚪ specified, not implemented
**Triggers to start**: this is a self-contained piece of work — can be done
   in any session. Doesn't conflict with the parallel Diversity Atlas
   (touches a separate track in MODULE_5A2 and a separate panel in the
   Inversion Atlas).
**Estimated scope**: 1 day for the detector script + JSON contract;
   half a day for the atlas-side panel; another half day for plot/figure
   polish for manuscript.
**Manuscript framing**: this is a flagship novelty result for the
   inversion paper — at 5-kb per-window per-sample resolution in a
   226-sample 9×-coverage hatchery cohort, no prior animal study has
   resolved gene-conversion-like vs double-crossover-like events
   inside polymorphic inversions. Both are RARE because inversions
   suppress recombination — detecting them at all is the point.

---

## What's being detected

Within an inversion candidate's L2 envelope, each sample has a
window-by-window band assignment (from `kmeans-K3` or whichever K2/K3
clustering is active). Most samples stay in one band across the entire
candidate. A small minority briefly **leave their band and return** —
visible in the per-sample-lines panel as a sample whose PC1 trajectory
suddenly dips/rises, traverses through another band's territory, then
comes back.

The biological interpretation depends on TRACT LENGTH:

| Tract length | Likely mechanism | Notes |
|---|---|---|
| `< 10 kb` (1-2 windows at 5-kb step) | Gene-conversion-like | Non-crossover homologous repair tract; classic GC tracts in animals are 50 bp–10 kb |
| `10-100 kb` | Ambiguous (long GC tract OR short double-CO) | Flag as such; suggest manual inspection |
| `> 100 kb` | Double-crossover-like | Two reciprocal recombinations flanking the swapped segment |
| `> 500 kb` | Very strong double-CO signal | Reportable as "near full inversion-region replacement" |

**Why this matters**: inversions suppress recombination to first order.
Detecting any switching event at all means recombination suppression is
leaky in specific places — direct evidence supporting the boundary
erosion / recombination-suppression-relaxation framework you've already
developed in the founder null model. Per-event localization tells you
WHERE in the inversion suppression breaks down.

---

## Spec — the detector

### Inputs

The detector operates on the per-window per-sample band assignments
that the L2/L3 layer already computes. Specifically:

- A 2D matrix `band[s, w]` of integer band IDs (1, 2, 3 for K=3) for
  each sample s and window w, restricted to the windows in the
  candidate's L2 envelope (or in a manually-specified region)
- Window coordinates (start_bp, end_bp, center_mb)
- Candidate metadata: id, scan_start_bp, scan_end_bp, K
- Sample metadata: per-sample known karyotype assignment (the "majority
  band" within this candidate)

### Algorithm

For each sample s:

1. Take the band sequence `band[s, :]` across all windows in the candidate
2. Find the sample's "home band" = the modal band across all windows
3. Run-length-encode the sequence: produce a list of `(band, w_start,
   w_end, length_bp)` runs
4. **Mark events**: any run where `band ≠ home_band` AND
   `length_windows >= MIN_EVENT_WINDOWS` (default 2)
5. For each event, record:
   - `sample_id`
   - `home_band`
   - `event_band` (the band the sample switched into)
   - `w_start`, `w_end` (window indices)
   - `bp_start`, `bp_end`
   - `length_windows`
   - `length_bp`
   - `flanking_runs`: did the sample return to home_band on both sides?
     (`true` = double-CO-like; `false` = boundary effect or
     end-of-candidate)
   - `event_class`: derived from `length_bp`:
     - `gene_conversion_like` if `length_bp < GC_MAX_BP` (default 10000)
     - `ambiguous` if `length_bp < DC_MIN_BP` (default 100000)
     - `double_crossover_like` otherwise
   - `flanking_polarity`: if `flanking_runs == true`, did the sample
     come back to the SAME band on both sides (`return`) or settle in a
     DIFFERENT band on one side (`partial`)? Pure double-CO requires
     `return` on both sides.

### Edge case handling

- **CRITICAL: Cramér's V (or NMI/AMI/ARI) gating per window**. The K-means
  band assignment at a window is ONLY meaningful when the window's local
  PCA is genuinely structured into K bands aligned with the focal
  candidate. Per window, compute a contingency-strength metric (Cramér's
  V is the conventional choice; NMI / AMI / ARI are alternatives — the
  multi-metric chip on the per-sample-lines panel computes all four)
  comparing the window's local K-means clustering against the focal
  candidate's reference clustering. The V heat-strip already visible in
  the per-sample-lines panel uses this exact metric (blue = low V =
  local PCA at this window doesn't reflect candidate structure;
  yellow/red = high V = local clustering aligns with candidate karyotype
  assignment).

  The detector MUST mask out low-V windows. Define `MIN_V` (default 0.5,
  needs empirical calibration). Windows with `V < MIN_V` are NOT counted
  as evidence of any band — they're treated as missing data in the run-
  length encoding.

  This means an "event" is only flagged when:
  ```
  high-V window(s) in home_band  →  high-V window(s) in event_band  →  high-V window(s) in home_band
  ```
  ALL three runs must be composed entirely of high-V windows. If a low-V
  window falls in the middle, it's a gap, not part of the event.

  Equivalently: a flagged event requires that at every window in the
  event tract, the band assignment is statistically meaningful (the
  local PCA at that window genuinely separates samples into K=3 bands
  the way the focal candidate does). Otherwise the "switch" is noise
  from a window that doesn't have karyotype-resolving power.

  This is rigour, not just conservatism — it ensures every flagged
  event is one where the detector can defend "the sample's local
  karyotype was X, then was Y, then was X" with statistical backing
  at each point.

- **Single-window flips**: with 9× coverage, occasional band misassignment
  per window happens (~1-3% of windows in noisy regions). The
  `MIN_EVENT_WINDOWS = 2` threshold filters these. For ambiguous cases,
  also report `single_window_flip_count` per sample so the user can
  decide whether 1-window flips are signal vs noise.
- **Candidate boundaries**: at the inversion flanks, samples drift
  between bands as the inversion's local-PCA effect fades. The detector
  should restrict to the candidate's INTERIOR (default: skip outer 10%
  of windows on each end). Events in the flank zones go into a
  `flank_events` separate list, NOT counted as evidence of internal
  recombination suppression breakdown.
- **K=2 vs K=3**: at K=3, "switching from het to homozygous-A" is a
  different biological event than "switching from homo-A to homo-B".
  Record `event_band` AND `home_band` so both can be inspected. For the
  manuscript, focus on **homozygous-to-homozygous switches** (the most
  strongly diagnostic of recombination event).
- **Karyotype dropout**: some samples have ambiguous PC1 throughout
  (low confidence in any band). Skip these. Use
  `confidence[s, w] >= MIN_CONFIDENCE` (default 0.6) — confidence comes
  from the existing K-means soft-assignment / silhouette score.

### Outputs

One JSON per chromosome per candidate:

```json
{
  "tool": "regime_change_detector_v1",
  "schema_version": 1,
  "generated_at": "2026-04-30T...",
  "cohort": "C_gariepinus_226_hatchery",
  "chrom": "C_gar_LG28",
  "candidate_id": "cand_028_001",
  "candidate_scan_bp": [7900000, 8400000],
  "candidate_K": 3,
  "n_windows_in_candidate": 26,
  "n_windows_interior": 22,
  "n_windows_high_v": 19,
  "min_event_windows": 2,
  "v_metric": "cramers_v",
  "v_threshold": 0.5,
  "gc_max_bp": 10000,
  "dc_min_bp": 100000,
  "events": [
    {
      "sample_id": "CGA191",
      "home_band": 1,
      "event_band": 3,
      "w_start": 1118,
      "w_end": 1119,
      "bp_start": 7910000,
      "bp_end": 7920000,
      "length_windows": 2,
      "length_bp": 10000,
      "flanking_runs": true,
      "flanking_polarity": "return",
      "event_class": "gene_conversion_like",
      "confidence_min": 0.78,
      "v_min_in_event": 0.62,
      "v_mean_in_event": 0.71,
      "v_min_in_flanks": 0.58,
      "v_threshold_used": 0.5
    },
    ...
  ],
  "flank_events": [...],
  "summary": {
    "n_samples_with_any_event": 8,
    "n_events_total": 11,
    "events_by_class": {
      "gene_conversion_like": 6,
      "ambiguous": 3,
      "double_crossover_like": 2
    },
    "samples_by_class": { ... }
  }
}
```

---

## Implementation location

The detector is a post-processor on the existing L2/L3 band assignments.
Two natural places:

1. **As a new step inside MODULE_5A2** (where breakpoint validation
   already lives — `MODULE_5A2_breakpoint_validator_standalone.py` is
   the precedent). Add `MODULE_5A2_regime_change_detector.py` (or `.R`).
2. **As an Atlas-side computation** triggered when a candidate is
   selected, reading band assignments from the chrom JSON. Faster
   iteration but doesn't get persisted across sessions.

Recommend option 1 — it persists across atlas reloads, can be batch-run
across all candidates, and the JSON output integrates with the existing
MODULE_5A2 / Atlas data flow.

The script reads the chrom JSON's L2/L3 layer (which already has
`per_sample_per_window_band` matrices), identifies all candidates in
`state.candidates`, runs the detector per candidate, emits one JSON per
candidate (or one JSON per chrom containing all candidates).

---

## Atlas-side panel

When an event JSON is loaded for the active candidate, the boundaries
page (or a new page) should show:

### Panel A — event summary table

Same idiom as the class summary auto-scan panel (which already exists).
Sortable table:

```
flag | sample      | event_class           | length_bp | from→into→back | confidence
●    | CGA191      | gene_conversion_like  | 10 kb     | 1→3→1          | 0.78
●    | CGA112      | double_crossover_like | 350 kb    | 2→1→2          | 0.85
○    | CGA045      | ambiguous             | 45 kb     | 1→3→1          | 0.62
```

Click a row → highlights that sample in the per-sample-lines panel
above (sets it as a tracked sample with a special "regime change"
marker glyph).

### Panel B — per-sample-lines integration

When event-detector data is loaded:

- Samples with events get a small triangle/diamond marker overlaid on
  their PC1 trajectory at the event window position
- Marker color: gene-conversion-like = green, ambiguous = amber,
  double-CO-like = red
- Hover a marker → tooltip with sample id, event class, length, plus
  per-event V values (`v_min_in_event`, `v_mean_in_event`)
- The existing V heat-strip below the per-sample-lines panel is
  retained — events should ONLY appear in windows where the heat-strip
  is yellow/red (high V), confirming visually that the V-gating is
  working. If a marker appears over a blue (low-V) zone, the gating
  failed and the event is suspect.

### Panel C — manuscript figure mode

A standalone view (or export option) that plots:

- Top: candidate position on chromosome (small ideogram)
- Middle: per-sample-lines panel with event markers
- Bottom: event summary table

Export to SVG with academic theme styling. Goal: drop directly into
the manuscript as Figure X.

---

## What this means for the manuscript

This is a STORY worth a dedicated figure + supplementary section in the
inversion manuscript. Suggested structure:

### Main text

A subsection in Results titled something like:

> "Sub-megabase-resolution detection of recombination events within
> polymorphic inversions"

Two-three paragraphs:
1. Methodology recap (per-window per-sample band assignment at 5 kb,
   detector logic, GC vs DC classification by tract length)
2. Catalogue of detected events across all candidates (numbers, length
   distribution, per-sample distribution)
3. Implications for inversion theory (boundary erosion already invoked;
   this is the per-window evidence of suppression breakdown)

### Figure

Three panels:
- (a) per-sample-lines for one striking candidate showing 1-2 events
- (b) histogram of event tract lengths colored by class
- (c) cumulative count of events by candidate, ordered by chromosome

### Supplementary

- Per-candidate event tables (one row per event)
- Per-candidate event summary stats
- Notes on noise floor, window-flip filtering, K=3 vs K=2 differences

---

## Rigor / things to watch

To avoid overclaiming:

1. **The V-gating is the cornerstone of detector rigour**. Every flagged
   event must be defensible as "the band assignment was statistically
   meaningful at every window in the event tract". Without V-gating, a
   "switch" in a window with collapsed local PCA (no karyotype-
   resolving power) is just noise. WITH V-gating, the manuscript can
   defend each event as "at every window in this 50-kb tract, the
   local genotype data clustered samples into the same K=3 groups as
   the focal candidate, and this sample's group assignment changed
   coherently". That's a strong claim. Without V-gating, it's
   indistinguishable from a single-sample band-misassignment cascade.

   For the manuscript: report the V threshold used, the per-event V
   distribution, and demonstrate that lowering the V threshold
   inflates the event count (showing the detector is correctly noise-
   filtering rather than missing signal).

2. **The detector must distinguish noise from signal**. Run the detector
   ALSO on shuffled control data (randomize window order or sample
   labels). If shuffled data gives ≥10% as many events as real data,
   the noise floor is too high to trust the signal. Report this control.

3. **The "double crossover" claim is mechanistic**. The atlas-detected
   pattern is "sample's band changed for X windows then changed back"
   — that's CONSISTENT with a double crossover, but other mechanisms
   could produce it (e.g. assembly error in the reference, mismapping
   in a TE-rich region, contamination). For each event reported as
   double-CO-like in the manuscript:
   - Spot-check on IGV: is read coverage uniform? Discordant pairs?
     Split reads at the event flanks?
   - Check if event flanks overlap repeat-rich regions (use the new
     TEfull layer!)
   - Confirm the sample passes overall QC (genome-wide H, F_ROH, etc.)

4. **Sample size justifies the conclusion**. With 226 samples and
   typical inversion candidates having 10-100 windows interior, you have
   ~2200-22000 sample-window observations per candidate. Even a 0.1%
   noise floor would produce 2-22 spurious events per candidate. The
   detector's output must be filtered against this floor.

5. **Cross-validate against the existing breakpoint validator**.
   MODULE_5A2 has Fisher's/chi-square/Cochran-Armitage tests for
   karyotype-genotype association at candidate flanks. Events near
   either flank should agree with what the breakpoint validator
   reports. Disagreements are diagnostic of either detector or
   validator bugs.

6. **The "no one has ever seen this" claim** is true at this scale and
   resolution combination, BUT relevant prior work exists in:
   - Drosophila inversion genetics (gene conversion within inversions:
     Korunes & Noor 2017 Genetics)
   - Sunflower inversions (Todesco et al. 2020 Nature: "Massive haplotypes
     underlie ecotypic differentiation in sunflowers" — has some
     within-inversion recombination evidence)
   - Avian inversion polymorphisms (ruff, white-throated sparrow)
   Cite them in the manuscript with a clear "this work extends to
   higher resolution and larger sample size in animals" framing rather
   than "no one has ever seen this".

---

## Open questions for the next session that takes this on

1. Implement detector in MODULE_5A2 (Python) or as standalone Atlas-
   side JS? Recommend MODULE_5A2.
2. **What V threshold to use?** The default `MIN_V = 0.5` needs
   empirical calibration. Plot the per-candidate V distribution and
   pick a threshold that captures the "yellow/red zone" — likely 0.4-
   0.6. Could also use a chromosome-specific or candidate-specific
   threshold (e.g. top 75% of windows) rather than a fixed value.
3. **Which V metric?** Cramér's V is the conventional choice and it's
   what the existing per-sample-lines panel V heat-strip uses. NMI /
   AMI / ARI are alternatives that the multi-metric chip on that panel
   computes. Recommend Cramér's V for consistency with the existing
   visualization, but if any of the alternatives have better
   discrimination at this scale, defer to that.
4. What thresholds for `MIN_EVENT_WINDOWS`, `GC_MAX_BP`, `DC_MIN_BP`?
   Defaults proposed above (2 / 10 kb / 100 kb) need empirical
   calibration against shuffled controls.
5. Should the detector run on all candidates in batch, or per-candidate
   on demand? Recommend batch — small per-candidate compute means
   "all candidates" is still fast.
6. What confidence metric to use for the band assignment? Soft
   K-means assignment probability is the obvious choice but silhouette
   width or per-window L2 z-score might be better.
7. Cross-link to TEfull layer: do events disproportionately occur in
   TE-rich regions? Quick analysis once both layers exist.

Bring these to the implementation session.
