# STEP_QR02 + QR03 — Q regime metrics and detector spec (v2)

**Status:** v2, April 30 2026. Supersedes v1 by expanding the QR02 metric
section with more detail on what each metric catches and why.
**Sources:**

- v1: "Nested subgroup extraction" chat (April 25,
  `https://claude.ai/chat/92fef806-5a81-45f7-bbdc-eaf0e1df6791`)
- v2 additions: the screenshots-batch on detailed QR02 metrics
  ("Component-by-component tracks per sample", "Δ-track of the Q vector",
  "Why your figure doesn't show the regime shift", "What we actually
  decided about NGSadmix")

**Purpose:** turn STEP_QR01's per-sample × per-window local Q matrix into
candidate intervals on a chromosome. This is the discovery step that
promotes ancestry-Q from "tracks-only" to "candidate detector,"
at which point the Q page can render envelopes (per ADR-13).

## Changes from v1

- **§2 expanded:** the QR02 metric set is now detailed with explicit
  rationale for each track, including the worked example showing why
  Δ_ij beyond Δ12 matters
- **§2 reweighted:** `exp(H)` flagged as the cleanest single track
  (regime shift independent of component swap)
- **§2 new sub-section:** Component-by-component per-sample tracks —
  K independent matrices, one per Q component, asking whether the
  spread of that component's loading across samples shifts at the
  window. **This is what catches inversions originating in one breeding
  line and now segregating** (the LG28 case)
- **§2 new sub-section:** Δ-track of Q vector formalized as the L1/L2
  distance from chromosomal-background cohort-mean Q vector — the
  analog of GHSL Part A's contrast metric but in Q space
- **NEW §6:** NGSadmix global/baseline rationale (why genome-wide K=8
  Q is appropriate as baseline, why inversions don't bias it)
- **NEW §7:** hatchery-cohort marketable angle for the manuscript
- **NEW §8:** 4-vs-5 stream resolution after rare-MAF demotion

---

## 1. Why the current dom.Q + Δ12 approach fails on LG28

From the LG28 figure analysis: the canonical inversion at 15–18 Mb
shows **no shift in dom.Q** and a Δ12 pattern that's "distinct" but
not clearly elevated relative to chromosomal background. Two
possible explanations:

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
That said, central-tendency metrics are still computed and shipped
because they're cheap (ADR-12) and they're informative when they DO
shift.

---

## 2. STEP_QR02 outputs (the inputs to QR03)

QR02 must emit, per window, the following metric families:

### 2.1 Cohort-level central-tendency tracks (existing)
- `dom_Q` — dominant ancestry component per window (the existing track)
- `Δ12` — gap between top-1 and top-2 ancestries
- `cohort_mean_Q[K]` — per-component cohort-mean (K values per window)

### 2.2 Cohort-level diversity tracks (the cleanest single tracks)

**`exp_H` — per-window expected heterozygosity of the Q distribution.**
Treat the cohort-mean Q vector as a probability distribution and
compute Shannon entropy `H = -Σ q_k log q_k`, then `exp(H) ∈ [1, K]`.

- Low `exp(H)` → one component dominates (assignment is concentrated)
- High `exp(H)` → ancestry is mixed across multiple components
- A regime change shows up as a **step in exp(H), independent of
  which component is dominant**

**This is probably the cleanest single track for catching ancestry-
regime shifts** because it doesn't depend on a component swap — only
on the sharpness of the assignment changing.

**`ENA` — effective number of ancestries.**
`ENA = 1 / Σ q_k²` (inverse Simpson index of the Q vector). Bounded
between 1 (one component carries all weight) and K (uniform across
components). Reports the same biology as exp(H) but with different
sensitivity to small components — Simpson is less sensitive to rare
components than Shannon. **Worth reporting both since they can
disagree informatively.**

### 2.3 Pairwise gap tracks

**`Δ_ij` for the full top-set.** Pairwise differences across the full
Q vector, not just top-vs-second. With K=8 there are C(8,2) = 28
pairwise gaps. Most aren't useful, but a small subset is. Specifically:

- `Δ12` (1st minus 2nd) — existing, retained
- `Δ13` (1st minus 3rd) — tells you whether the sub-dominant signal
  is a clean second component or a tie between several components
- `Δ14`, `Δ23`, `Δ24` — generally useful

**Worked example showing why Δ13 matters:**
- Window A: Q1=0.5, Q2=0.25, Q3=0.25 → Δ12=0.25, Δ13=0.25
- Window B: Q1=0.5, Q2=0.45, Q3=0.05 → Δ12=0.05, Δ13=0.45

Same dominant component (K1 at 0.5 in both), totally different
ancestry structure. Δ12 alone misses this. The full Δ_ij set
catches it.

**Default emission:** `Δ12, Δ13, Δ14, Δ23, Δ24`. The other 23 pairwise
gaps are noise at K=8 — emit only if a future case justifies.

### 2.4 Cohort-level dispersion tracks (high-priority additions)

**`cohort_var_Q[K]` — per-component cohort-variance.** For each
component k, the variance of sample-level Q_k values across the
cohort at this window. **Jumps at inversions where the carrier
set has a different Q profile than non-carriers**, even when the
cohort-mean is unchanged.

### 2.5 Cohort-level bimodality tracks

For each component k = 1..K:
- `dip_stat_k` — Hartigan's dip statistic on the cohort vector of
  sample-level Q_k values at this window
- `dip_p_k` — p-value from `diptest::dip.test`
- `is_bimodal_k` — boolean from `dip_p_k < 0.05`

Aggregate: `any_bimodal` = TRUE if any k has `dip_p_k < 0.05`.

**Per-component bimodality is the most informative single track when
the inversion has a discrete carrier set in one breeding line** (the
LG28 case). Carriers cluster around Q_k ≈ 0.6, non-carriers around
Q_k ≈ 0.1; the per-sample distribution of Q_k is bimodal at the
inversion windows and unimodal in flanks.

### 2.6 Component-by-component per-sample tracks (NEW v2)

Beyond cohort-mean Q, compute a **per-sample × per-window matrix
per Q component** — K=8 such matrices, one per component. For each
component independently, ask:

> Does the spread of this component's loading across samples shift
> at this window?

Operationally, for each component k and each window w:
- Compute `mean_k_w` = cohort mean of Q_k at window w
- Compute `var_k_w` = cohort variance of Q_k at window w
- Compute `frac_present_k_w` = fraction of samples with Q_k > threshold
  (e.g. > 0.05 or > 0.10)

**Why this catches LG28-type inversions:**
Component K3 might be present in 40% of samples in flanks but in 100%
of samples inside an inversion — a ~60-point increase in cohort-mean
K3 with reduced variance. **This is what catches inversions that
originated in one breeding line and are now segregating** in the cohort:
K3 carriers are the inversion carriers, the inversion's discovery signal
is the per-sample variance-reduction in K3 loading inside the candidate.

This is a different metric family from §2.4 (cohort-variance summed
across all components) — it's K independent variance tracks, one per
component, so the discovery step can identify which component is
involved.

### 2.7 Δ-track of the Q vector (NEW v2 — the GHSL-contrast analog)

For each window, compute the L1 or L2 distance between the window's
cohort-mean Q vector and the chromosomal-background cohort-mean Q
vector:

```
L1_dist[w] = Σ_k | mean_Q_k[w] − mean_Q_k_background |
L2_dist[w] = sqrt( Σ_k ( mean_Q_k[w] − mean_Q_k_background )² )
```

This collapses the K-dimensional vector to a single scalar — distance
from chromosome background — and lets you walk the chromosome looking
for regions where ancestry composition deviates from background. Big
spike = candidate.

**This is the analog of what GHSL Part A's contrast metric does, but
in Q space.** The architectural symmetry matters: GHSL detects local
divergence in haplotype space, θπ detects local divergence in
diversity space, this Δ-track detects local divergence in ancestry-Q
space. Three streams, three contrast metrics, same envelope-detection
pattern in QR03 / GHSL detector / θπ detector.

### 2.8 Per-sample tracks (preserve sample identity)

- `sample_Q_deviation[N_samples]` — per-window: for each sample i,
  deviation_i = ‖q_iw − q_i_baseline‖_1 where q_i_baseline is sample
  i's genome-wide-mean Q vector
- `mean_sample_deviation` — cohort mean of the above
- `max_sample_deviation` — cohort max
- `frac_samples_high_deviation` — fraction with deviation > threshold
  (default threshold: chromosomal 90th percentile)

These are the **per-sample resolution layer** — they preserve sample
identity for downstream phase-4 work. Note: in the per-sample QR01
matrix, this is the input layer for CUSUM-style per-carrier
boundary detection (analogous to R02 / T05 for GHSL / θπ). **Whether
to add a fourth CUSUM (Q-regime per-sample CUSUM) is deferred** —
CUSUM_SPEC.md §11 explicitly rejected this for now because Q-regime
detection wants dispersion/bimodality not level shifts. Re-evaluate
once QR02 produces real data.

---

## 3. STEP_QR03 detector — multi-track Z-score + combined score

### 3.1 The detection model

Walk the chromosome window-by-window. For each window w, compute a
**combined outlier score** from the QR02 tracks. Detect intervals
where the combined score is sustained above a threshold for ≥ N
consecutive windows.

This is **the same architectural pattern as the dosage |Z| detector
and the GHSL secondary detector** (per ADR-2). Same kind of detector,
different input track family.

### 3.2 Per-track Z-scoring

For each scalar track T (`dom_Q`, `Δ12`, `Δ13`, `exp_H`, `ENA`,
`L1_dist`, `mean_sample_deviation`, etc.):

```
T_z[w] = (T[w] - median(T over chromosome)) / MAD(T over chromosome)
```

Use median + MAD, not mean + sd, to be robust against the inversion
itself shifting the chromosomal moments. This matches the dosage
detector's robust-Z convention.

For per-component tracks (`cohort_var_Q[k]`, `dip_stat_k`, K
component-by-component matrices from §2.6):
- Compute per-component robust Z within the chromosome
- Aggregate via `max_k(|T_z_k|)` — i.e., take the most-extreme
  component's z-score. This catches inversions enriched in any single
  K-component without requiring all components to flag the same window

### 3.3 Combined score (revised v2 weights)

The combined score is a **weighted sum of |z|-scores across track
families**. Initial weights (to be tuned on LG28):

```
combined_z[w] =
    0.25 * |z(exp_H)|                    # diversity (cleanest single track)
  + 0.20 * |z(component_var_max_k)|      # component-by-component dispersion
  + 0.15 * |z(dip_stat_max_k)|           # per-component bimodality
  + 0.15 * |z(L1_dist_baseline)|         # Q-vector contrast (GHSL analog)
  + 0.10 * |z(mean_sample_deviation)|    # per-sample deviation
  + 0.05 * |z(Δ12)| + 0.05 * |z(Δ13)|    # central-tendency (existing)
  + 0.05 * |z(ENA)|                      # alternative diversity
```

Reweighting rationale vs v1:
- `exp(H)` promoted to top weight because it's the cleanest single
  track and doesn't depend on component swap
- Component-by-component variance promoted because it's specifically
  what catches breeding-line-segregating inversions (LG28 case)
- L1_dist promoted to mid-tier because it's the GHSL-contrast analog
  and architecturally important
- Δ12 demoted to low weight (still emit, but it's known to fail
  on LG28-type cases)

**These are starting weights, to be tuned on the LG28 validation
test (§5).**

### 3.4 Threshold selection

Same as v1: quantile-based or adaptive thresholding, not fixed Z.
See v1 §3.4 for full rationale (the GHSL z-score lesson).

### 3.5 Run-length aggregation

Same as v1: L1 and L2 envelopes matching dosage / GHSL / θπ
conventions. See v1 §3.5.

Output records per L2:
- `q_l2_id`, `win_start`, `win_end`, `start_bp`, `end_bp`, `n_windows`
- `peak_combined_z`, `mean_combined_z`
- `dominant_signal` — which track family contributed most to the peak
  (diversity / component_dispersion / bimodality / contrast /
  sample_deviation / central_tendency)
- `bimodal_components` — list of K-components flagged bimodal in the
  window range
- `enriched_components` — list of K-components flagged via component-
  by-component variance reduction (NEW v2 field — links the
  candidate to a specific breeding line via its dominant K-component)

The `enriched_components` field is what makes phase-4 dual-clustering
straightforward: the candidate's primary K-component identifier
points directly at the carrier-line ancestry, which is the input to
the family/Q confounder test.

---

## 4. Honest scope

The "Nested subgroup extraction" chat commits to:

> "STEP_QR03 turns those metrics into candidate intervals. The
> current dom.Q + Δ12 approach is one input to QR03, not the whole
> thing. ... QR03 is about 100 lines of detection logic."

That scope is unchanged in v2. **QR03 itself is still ~100–200 lines
of R**, same envelope-detection pattern as dosage / GHSL / θπ.

**QR02 grows in v2.** With per-component matrices added (§2.6),
QR02 needs to handle K=8 separate per-sample × per-window matrices
plus the cohort-mean tracks. That's still tractable — the output is
a single TSV with K-fold-replicated columns for component-specific
metrics.

Updated scope estimate:
- **STEP_QR02:** ~3 days now (was 2). Component-by-component matrix
  computation adds the most overhead.
- **STEP_QR03:** ~1 day (unchanged from v1)
- **export_q_to_json.R:** ~1 day (unchanged)
- **Atlas Q page:** ~2 days (unchanged)
- **diptest R package availability check** in `assembly` env

Total ~7–8 days if validation passes. Up from ~6 in v1, mostly
because §2.6 component-by-component matrices add real work.

---

## 5. Validation against LG28 (the killer test)

Same as v1 — restrict the local Q matrix to:
- Candidate region: 15–18 Mb
- Flanking control regions: 13–15 Mb and 18–20 Mb

Compute per-region:
- Cohort-mean dom_Q (existing — should be flat across all 3 regions)
- Cohort-variance per component (should jump in the candidate)
- **Component-by-component cohort-mean and variance per K**
  (NEW in v2 — should show one specific component shifting up
  with reduced variance, identifying the carrier line)
- Dip test per component (should be bimodal in the candidate for
  the carrier-line component)
- Per-sample Q deviation from baseline (should jump for carriers)
- **L1 distance of cohort-mean Q vector from genome-wide baseline**
  (NEW in v2 — should spike in the candidate)
- exp(H) per region (NEW emphasis — should show a clear step
  at the inversion boundary)

**Predicted result if the spec is right (v2 update):**
- Cohort-mean dom_Q: flat across all 3 regions (matches the figure)
- Mean exp(H): **steps at the candidate boundary** — this is the new
  prediction; v1 said exp(H) was "probably flat too," but v2 says
  exp(H) should actually shift cleanly because the assignment
  sharpness changes
- Component-by-component K3 (or whichever): mean jumps, variance
  drops inside the candidate
- L1 distance from baseline: spikes inside the candidate
- Dip test per component: bimodal for the carrier-line component
- Per-sample deviation: bimodal across samples in the candidate

If exp(H) and L1_dist BOTH show the predicted shift → ship.
If only one shows it → still ship but reweight the combined score.
If neither shows it → reconsider the entire framework.

This validation is **a one-script test, runs in seconds**, and is
the gating step before QR03 goes into production.

---

## 6. NGSadmix baseline rationale (NEW v2)

Why the genome-wide Q (K=8 NGSadmix on NAToRA-pruned 81 unrelated
individuals) is the appropriate baseline:

**The genome-wide Q is a clean ancestry assignment that isn't
contaminated by the structure inside any specific inversion candidate.**
This is the **family-bias-resistance argument**: inversions are a
tiny fraction of the genome (a handful of regions across 28
chromosomes), so they contribute negligibly to the genome-wide Q
decomposition. The genome-wide signal averages over all loci, and
any single inversion contributes negligibly.

Operational consequences:
- **Q baseline is genome-wide K=8**, not per-chromosome
- **`q_i_baseline`** in §2.8 (per-sample deviation from baseline) is
  sample i's genome-wide-mean Q vector, not per-chromosome
- **`mean_Q_k_background`** in §2.7 (Δ-track L1/L2 distance) is
  the chromosomal-background cohort-mean Q vector — i.e., the
  cohort mean across all windows on this chromosome that aren't
  flagged as candidates by the dosage / GHSL / θπ streams
  (use existing candidate masks from those streams to define
  "background")

**Why use chromosomal-background rather than genome-wide for the
Δ-track baseline:** chromosomal-background controls for chromosome-
level structure (e.g., chromosomes vary in their mean K-component
profile). Genome-wide baseline would conflate chromosomal-level
ancestry signal with within-chromosome regime shifts.

Compromise: use **genome-wide baseline for per-sample deviation
(§2.8)** — preserves sample-identity-meaningful baseline — and
**chromosomal-background for the Δ-track (§2.7)** — controls for
chromosome-level structure when looking for windows that deviate
from the rest of the chromosome.

**Don't re-run NGSadmix per-chromosome.** That's a temptation that
arises from "K=8 doesn't separate breeding lines finely enough"
(§1 explanation B). The right response is not to re-run NGSadmix at
higher K or per-chromosome (which would be expensive and would
contaminate the baseline with the very inversions we're trying to
discover). The right response is to use **per-component variance
reduction** (§2.6) which catches the LG28-type case even when the
carrier line is grouped with non-carrier lines into a single K=8
component.

---

## 7. Hatchery-cohort marketable angle (NEW v2)

For the manuscript: ancestry-Q regime as a primary discovery signal
makes this pipeline the **first published inversion pipeline to
properly handle hatchery / breeding-program structure**. Most
published inversion methods are built for wild populations where
ancestry is unknown or unstructured. A pedigree-aware pipeline that
treats Q-regime as a first-class signal is a **methodological
contribution and a marketable angle**.

This positioning supports the manuscript's framing: "we developed
this pipeline specifically for aquaculture genomics in cohorts with
known pedigree structure." Reviewers like methods that match the
data they're applied to.

The 11-track shelf-QC view (per ADR-13's Ancestry page) reinforces
this: eleven aligned tracks per chromosome × candidate, generated in
seconds, is the kind of figure that makes reviewers trust the
methodology because they can see all the evidence at once and read it
themselves.

---

## 8. 4-vs-5 stream resolution after rare-MAF demotion (NEW v2)

The other-Claude conversation said "discovery-stream count to five,
and the architecture absorbs it without modification" when ancestry-Q
was added. **That count was 5: dosage + GHSL + θπ + rare-MAF +
ancestry-Q.**

This session demoted rare-MAF to the diversity atlas (per the
session-opening rare-MAF triage). **Final committed count: 4 streams:
dosage + GHSL + θπ + ancestry-Q.**

The "absorbs without modification" claim still holds — the
architecture is a per-stream pattern (each stream has its own
discovery page mirroring the others). Whether the count is 4 or 5
doesn't change the pattern. Demoting rare-MAF strengthens the
story (no noisy 5th stream diluting the convergent-evidence framing)
without breaking the architecture.

**Walk-back protection:** if a future session proposes "add rare-MAF
back as a 5th stream," reject. Reasons: rare_inversion_score has
family-structure noise issues (best-of-top-5 outliers rationalizes
noise, run-length boost amplifies sibling clusters), and
rare_sfs_pairwise is genome-wide / not per-window so it doesn't fit
the per-stream discovery page pattern. See the rare-MAF triage
discussion at the top of this session's transcript.

---

## 9. Deferred / open questions (was §6 in v1)

- **Multi-scale Δ12.** Whether QR02 emits multi-scale or just
  single-scale matters for QR03's detection. Decision deferred to
  QR02 implementation.

- **K-selection.** K=8 from MODULE_2B's NGSadmix is the baseline.
  Don't re-run NGSadmix per-LG or at higher K to "fix" LG28; use
  component-by-component variance reduction instead (§2.6 + §6).

- **Phase 4 dual-clustering integration.** ADR-13 covers this; QR03
  is the discovery side, `STEP_AR01_local_q_per_candidate.R` in
  the phase-4 sibling folder handles per-candidate per-sample
  ancestry karyotype clustering.

- **Whether to add a 4th CUSUM (Q-regime per-sample CUSUM).**
  CUSUM_SPEC.md §11 rejected this. Reconsider if QR02 produces
  per-sample matrices that show clean CUSUM-detectable level shifts.
  Default: no.

- **`enriched_components` field semantics.** A candidate with K3
  enriched suggests the K3 breeding line carries the inversion,
  but this needs the Hungarian-matched mapping from K-component
  to breeding-line label that's done in MODULE_2B. The field is
  emitted as raw K-component IDs; the breeding-line label
  translation happens at atlas render time.

---

## 10. Build order when implementing (was §7 in v1)

Same order as v1, with QR02 estimate revised upward:

1. **Validation script first** (§5) — one afternoon. Test predictions
   on LG28's local Q matrix. Decision gate.
2. **STEP_QR02 implementation** — ~3 days (was 2). Math is
   straightforward; the per-component dip test and the K
   component-by-component matrices are the non-trivial parts.
3. **STEP_QR03 implementation** — ~1 day (~150–200 lines).
4. **export_q_to_json.R** — packs QR02 tracks + QR03 envelopes
   into `LG##_q.json`. ~1 day.
5. **Atlas Q page** — ~2 days. Clones page 3/12 panel layout.

Total estimate: **~7–8 days** if all data is available and
validation gates pass cleanly.

---

## 11. JSON layer expansion (NEW v2)

`LG##_q.json` schema gains these fields per window:

```json
{
  "windows": [
    {
      "win_idx": 1234,
      "start_bp": 14000000,
      "end_bp": 14005000,
      "dom_Q": 3,
      "delta_12": 0.32, "delta_13": 0.45, "delta_14": 0.55,
      "delta_23": 0.13, "delta_24": 0.23,
      "exp_H": 2.34,
      "ENA": 2.18,
      "cohort_mean_Q": [0.12, 0.08, 0.55, 0.05, 0.10, 0.05, 0.03, 0.02],
      "cohort_var_Q":  [0.01, 0.01, 0.04, 0.02, 0.01, 0.01, 0.01, 0.01],
      "L1_dist_chr_bg": 0.42,
      "L2_dist_chr_bg": 0.21,
      "mean_sample_dev": 0.18,
      "max_sample_dev": 0.55,
      "frac_samples_high_dev": 0.24,
      "dip_stat":   [0.02, 0.03, 0.08, 0.02, 0.03, 0.02, 0.01, 0.02],
      "dip_p":      [0.45, 0.31, 0.001, 0.62, 0.28, 0.51, 0.78, 0.66],
      "is_bimodal": [false, false, true, false, false, false, false, false]
    }
  ],
  "candidate_envelopes": [
    {
      "q_l2_id": "LG28_Q_001",
      "win_start": 1234, "win_end": 1567,
      "start_bp": 14000000, "end_bp": 18500000,
      "n_windows": 333,
      "peak_combined_z": 4.7,
      "mean_combined_z": 2.8,
      "dominant_signal": "diversity",
      "bimodal_components": [3],
      "enriched_components": [3]
    }
  ]
}
```

Per-component arrays (cohort_mean_Q, cohort_var_Q, dip_stat, dip_p,
is_bimodal) have length K (=8 for the current MODULE_2B). The atlas
Q page renders selected components based on user choice (default:
the dominant or most-bimodal component for the candidate).

The component-by-component per-sample matrices from §2.6 are NOT
shipped in the JSON — they're too big (8 × 226 × 16,500 × Float32
≈ 120 MB per chromosome). They remain on the cluster, available for
re-computation if the user changes sample subset (same cluster-side
on-demand pattern as CUSUM in CUSUM_SPEC.md §8).

---

## 12. Walk-back protection (revised v2)

- **"Re-run NGSadmix per-chromosome to fix LG28"** — reject. Use
  per-component variance reduction (§2.6); see §6.
- **"Drop exp(H), use only Δ_ij"** — reject. exp(H) is the cleanest
  single track per the v2 reweighting; Δ_ij is supplementary.
- **"Component-by-component matrices are overkill"** — reject. They're
  what catches LG28-type inversions; without them the spec is back
  to v1's failure mode.
- **"Just ship dom.Q and Δ12, the manuscript only needs central
  tendency"** — reject. The manuscript's marketable angle (§7) is
  "first pipeline to handle hatchery structure properly" — that
  requires the dispersion/component-level metrics that don't fail
  on LG28.
- **"Add rare-MAF back as a 5th stream"** — reject. See §8.
- **"Add a 4th CUSUM for Q-regime per-sample"** — defer, default no.
  See §9.

---

## 13. tl;dr (revised v2)

- **QR02 emits 8 metric families:** central-tendency (dom_Q, Δ_ij),
  diversity (exp(H), ENA), cohort dispersion (cohort_var_Q),
  bimodality (dip_stat per K), component-by-component per-sample
  (mean_k, var_k, frac_present_k for each K), Δ-track (L1/L2 from
  chromosomal background), per-sample deviation (from genome-wide Q).
- **`exp(H)` is the cleanest single track** — regime shift independent
  of which component is dominant.
- **Component-by-component variance reduction is the LG28-killer** —
  catches inversions that originated in one breeding line and now
  segregate, where dom.Q is flat.
- **The Δ-track from chromosomal background is the GHSL-contrast
  analog in Q space** — architectural symmetry with the other streams.
- **NGSadmix stays genome-wide K=8 on NAToRA-pruned individuals** —
  family-bias-resistance argument; don't re-run per-chromosome.
- **QR03 is still ~150-200 lines** — same envelope detection pattern
  as dosage / GHSL / θπ.
- **4 discovery streams, not 5** — rare-MAF demoted to diversity atlas.
- **Manuscript angle:** first inversion pipeline for hatchery cohorts
  with known pedigree structure.
- **Build order unchanged:** validation → QR02 (~3 days) → QR03
  (~1 day) → exporter → atlas page. ~7-8 days total.
- **Don't build until LG28 dry-run validates the GHSL / θπ atlas.**
  Same gating as everything else.
