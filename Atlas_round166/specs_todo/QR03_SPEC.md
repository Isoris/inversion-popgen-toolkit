# STEP_QR03 — Q regime detector spec

**Status:** specification recovered from past-chat archaeology, April 30 2026.
**Sources:** the "Nested subgroup extraction" chat (April 25, 2026,
`https://claude.ai/chat/92fef806-5a81-45f7-bbdc-eaf0e1df6791`) is the
canonical reference; supplementary detail from the
"Insulation scoring" chat (April 10) and the LG28 Q-figure
analysis chat (April 25).

**Purpose:** turn STEP_QR02's per-window Q metrics into candidate
intervals on a chromosome. This is the discovery step that promotes
ancestry-Q from "tracks-only" to "candidate detector," at which point
the Q page can render envelopes (per ADR-13).

---

## 1. Why the current dom.Q + Δ12 approach fails on LG28

From the LG28 figure analysis: the canonical inversion at 15–18 Mb
shows **no shift in dom.Q** and a Δ12 pattern that's "distinct" but
not clearly elevated relative to chromosomal background. Two possible
explanations:

**A. Inversion shared across breeding lines.** If the inversion
predates breeding-line divergence, all carriers across all Q
components carry it at similar frequencies; cohort-mean dominant
ancestry stays flat. Signal is in **per-sample Q variance**
(carriers have different Q profiles than non-carriers, even if
cohort-mean dominant doesn't shift).

**B. K=8 doesn't separate breeding lines finely enough.** If the
carrier line is grouped with non-carrier lines into a single
K-component, the cohort-mean Q for that component is diluted by
non-carriers and doesn't shift at the inversion.

**Either way: dom.Q is the wrong metric. The right metrics are
dispersion / bimodality / per-sample deviation, not central-tendency.**

---

## 2. STEP_QR02 outputs (the inputs to QR03)

QR02 must emit, per window:

### 2.1 Cohort-level central-tendency tracks (existing)
- `dom_Q` — dominant ancestry component per window (the existing track)
- `Δ12` — gap between top-1 and top-2 ancestries
- `Δ_ij` for top 3–4 components — `Δ12, Δ13, Δ14, Δ23, Δ24` (subset of C(K,2) pairs; the rest are noise at K=8)
- `cohort_mean_Q[K]` — per-component cohort-mean (K values per window)
- `exp_H` — Shannon-based effective ancestry diversity scalar (on cohort-mean Q)
- `ENA` — Simpson-based effective ancestry diversity scalar (on cohort-mean Q)

### 2.2 Cohort-level dispersion tracks (NEW — the high-priority additions)
- `cohort_var_Q[K]` — per-component cohort-variance across samples (K values per window)
- `L1_dist_baseline` — L1 distance of window's cohort-mean Q vector from chromosomal-background mean Q vector

### 2.3 Cohort-level bimodality tracks (NEW — the most informative single track)
For each component k = 1..K:
- `dip_stat_k` — Hartigan's dip statistic on the cohort vector of sample-level Q_k values at this window
- `dip_p_k` — p-value from `diptest::dip.test`
- `is_bimodal_k` — boolean from `dip_p_k < 0.05`

Aggregate: `any_bimodal` = TRUE if any k has `dip_p_k < 0.05`.

### 2.4 Per-sample tracks (preserve sample identity)
- `sample_Q_deviation[N_samples]` — per-window: for each sample i, deviation_i = ‖q_iw − q_i_baseline‖_1 where q_i_baseline is sample i's genome-wide-mean Q vector
- `mean_sample_deviation` — cohort mean of the above
- `max_sample_deviation` — cohort max
- `frac_samples_high_deviation` — fraction with deviation > threshold (default threshold: chromosomal 90th percentile)

---

## 3. STEP_QR03 detector — multi-track Z-score + combined score

### 3.1 The detection model

Walk the chromosome window-by-window. For each window w, compute a
**combined outlier score** from the QR02 tracks. Detect intervals where
the combined score is sustained above a threshold for ≥ N consecutive
windows.

This is **the same architectural pattern as the dosage |Z| detector
and the GHSL secondary detector** (per ADR-2). Same kind of detector,
different input track family.

### 3.2 Per-track Z-scoring

For each scalar track T (dom_Q, Δ12, exp_H, ENA, L1_dist_baseline,
mean_sample_deviation, etc.):

```
T_z[w] = (T[w] - median(T over chromosome)) / MAD(T over chromosome)
```

Use median + MAD, not mean + sd, to be robust against the inversion
itself shifting the chromosomal moments. This matches the dosage
detector's robust-Z convention.

For per-component tracks (cohort_var_Q[k], dip_stat_k):
- Compute per-component robust Z within the chromosome
- Aggregate via `max_k(|T_z_k|)` — i.e., take the most-extreme
  component's z-score. This catches inversions enriched in any single
  K-component without requiring all components to flag the same window.

### 3.3 Combined score

The combined score is a **weighted sum of |z|-scores across track
families**. Initial weights (to be tuned on LG28):

```
combined_z[w] =
    0.30 * |z(cohort_var_Q_max_k)|       # dispersion family
  + 0.25 * |z(dip_stat_max_k)|           # bimodality family
  + 0.20 * |z(mean_sample_deviation)|    # per-sample deviation family
  + 0.10 * |z(L1_dist_baseline)|         # whole-Q-vector shift
  + 0.10 * |z(Δ12)|                      # central-tendency (existing track)
  + 0.05 * |z(exp_H)|                    # diversity
```

The weighting prioritizes **dispersion + bimodality + per-sample
deviation** because those are the hypothesis-driven additions that
would catch the LG28 case the existing dom.Q + Δ12 misses.

### 3.4 Threshold selection

Per past-chat discussion (the v5.2 vs v5.3 GHSL scoring problem in
the LG28 PASS=5/9192 chat): **don't pick a fixed Z threshold from
nothing**. Use one of:

- **Quantile-based:** flag windows where `combined_z[w]` exceeds the
  chromosomal Nth percentile (e.g., top 10–15%, like the GHSL fix
  proposed)
- **Adaptive:** use the same adaptive thresholding D17 uses
  (chromosome-specific quantiles of observed scores)

This is the same lesson learned from the GHSL z-score scoring failure:
fixed thresholds break when the baseline distribution is unusual.
Quantile-based is more robust to chromosomes that legitimately have
elevated baseline Q activity.

### 3.5 Run-length aggregation

A single high-z window is noise. Real candidates show sustained
elevation. Apply the same envelope-detection logic as dosage / GHSL:

- **L2 envelopes:** contiguous runs of `combined_z[w] > threshold`,
  minimum 5 windows (default), merge gap of 3 windows
- **L1 envelopes:** merge L2's within 9 windows of each other into
  L1 groups (matches the GHSL convention)

Output records per L2:
- `q_l2_id`, `win_start`, `win_end`, `start_bp`, `end_bp`, `n_windows`
- `peak_combined_z`, `mean_combined_z`
- `dominant_signal` — which track family contributed most to the peak
  (dispersion / bimodality / sample_deviation / central_tendency)
- `bimodal_components` — list of K-components flagged bimodal in the window range

---

## 4. Honest scope

The "Nested subgroup extraction" chat (April 25) commits to:

> "STEP_QR03 turns those metrics into candidate intervals. The
> current dom.Q + Δ12 approach is one input to QR03, not the whole
> thing. ... QR03 is about 100 lines of detection logic."

That's the right scope. **QR03 is not algorithmically novel.** It's
the same envelope-detection pattern as dosage / GHSL / θπ
(robust Z + run-length + multi-scale aggregation), applied to a
different family of input tracks. The novelty is in QR02
(computing the right Q metrics), not QR03 (detecting envelopes from
them).

This means QR03 is **probably a 100–200 line R script**, not a
2000-line endeavor. Most of the work is in:
1. Loading QR02 outputs
2. Computing robust Z per track family
3. Combining via weighted sum
4. Run-length envelope detection (mirror dosage/GHSL convention)
5. Writing out L1 + L2 envelope TSVs (mirror dosage/GHSL output schema)

The output schema must match the existing envelope TSV conventions
(start_bp, end_bp, n_windows, peak_score, status fields) so that
`export_q_to_json.R` can pack them into the Q page's JSON layer
identically to how GHSL/θπ envelopes are packed.

---

## 5. Validation against LG28 (the killer test)

Per the April 25 chat: there's a one-afternoon test that validates
the whole approach before committing to it. Take the existing local
Q matrix for LG28, restrict to:

- Candidate region: 15–18 Mb
- Flanking control regions: 13–15 Mb and 18–20 Mb

Compute per-region:
- Cohort-mean dom_Q (existing — should be flat across all 3 regions)
- Cohort-variance per component (NEW — should jump in the candidate)
- Dip test per component (NEW — should be bimodal in the candidate
  for whichever component carries the inversion)
- Per-sample Q deviation from baseline (NEW — should jump for carriers
  in the candidate, stay flat for non-carriers)

**Predicted result if the spec is right:**
- variance and bimodality elevated in the candidate, flat in flanks
- mean and exp(H) flat across all three regions
- per-sample deviation bimodal across samples in the candidate

If that prediction holds → ship QR02 + QR03 with the dispersion/
bimodality emphasis. If it fails → reconsider before building.

This validation is **a one-script test, runs in seconds**, and
should be the gating step before QR03 goes into production.

---

## 6. Deferred / open questions

- **Multi-scale Δ12.** The canonical figure shows Δ12 single-scale and
  multi-scale (presumably the rolling-smoothed version at multiple
  window sizes, like GHSL's s10/s50/s100). Whether QR02 emits
  multi-scale or just single-scale matters for QR03's detection — if
  multi-scale is in, QR03 should pick the most-informative scale per
  region, similar to how the GHSL detector picks primary_scale.
  **Decision:** TBD; depends on QR02 implementation choice.

- **K-selection.** The K=8 choice is from MODULE_2B's NGSadmix run.
  If K=8 dilutes the carrier line (option B in §1), QR03 might
  benefit from re-running QR01 at higher K or per-LG K. **Out of
  scope for QR03 itself**; this is a MODULE_2B decision.

- **Phase 4 dual-clustering integration.** Per ADR-13, the Ancestry
  page renders the 11-track shelf QC view that includes a per-sample
  ancestry_karyotype label. That label comes from clustering per-
  sample × K local Q vectors *inside the candidate* — k-means or
  k-medoids on samples in K-dimensional Q space, with k=2..5
  silhouette selection, Hungarian-matched against dosage_karyotype
  to align labels. **This is phase-4 work, not QR03.** QR03 only
  needs to produce candidate envelopes; the per-candidate
  per-sample Q clustering happens in `STEP_AR01_local_q_per_candidate.R`
  in the phase-4 sibling folder.

- **Confounder-aware classification.** The "ancestry as resolution
  layer" framing means QR01's per-sample × per-window matrix feeds
  the family/Q confounder test in phase 5. QR03 is not the consumer
  of that test; it's the discovery side. The same QR01 matrix
  underpins both.

---

## 7. Build order when implementing

1. **Validation script first** (§5) — one afternoon. Test predictions
   on LG28's local Q matrix. Decision gate.
2. **STEP_QR02 implementation** — produces all the new tracks (§2).
   Most of the work, ~2 days. Math is straightforward; the
   per-component dip test is the only non-trivial piece (depends on
   `diptest` R package being available in `assembly` env).
3. **STEP_QR03 implementation** — ~100–200 lines (§3). Matches the
   existing envelope-detection conventions of dosage/GHSL. Wall time
   should be seconds per chromosome.
4. **export_q_to_json.R** — packs QR02 tracks + QR03 envelopes into
   `LG##_q.json`. Mirrors `export_ghsl_to_json_v3.R` structure.
5. **Atlas Q page** — clones page 3/12 panel layout with
   `q*` prefixed canvases. Renders dom_Q track, Δ12 track,
   variance/bimodality tracks, envelope overlays.

Total estimate: **2–3 days for QR02, 1 day for QR03, 1 day for
exporter, 2 days for atlas page**. About a week if all data is
available and the validation gates pass cleanly.

---

## 8. What this spec does NOT commit to

- Exact field names in the output JSON (those follow `JSON_CONTRACT.md`
  conventions when implemented)
- Exact threshold values (these come from the LG28 validation, not
  from this spec)
- Whether to use raw Q or Q deviation as the primary discovery track
  (the spec recommends deviation; LG28 validation confirms or refutes)
- Multi-scale handling (deferred per §6)

These get resolved when the validation runs against real data, not
in advance.

---

## 9. tl;dr

**QR03 is a thin envelope detector** layered on top of QR02's metric
tracks. Same algorithmic pattern as the dosage and GHSL secondary
detectors (robust Z + run-length + L1/L2 aggregation). The novelty
is in QR02's metric choices: **dispersion + bimodality + per-sample
deviation**, not central-tendency.

The LG28 validation in §5 is the gate. Run it before committing to
the rest of the build.
