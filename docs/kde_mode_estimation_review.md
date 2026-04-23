# On KDE mode estimation for ancestral-fragment boundaries

**Your memory is correct.** Silverman-bandwidth KDE mode estimation is
the wrong tool for the fragment-boundary distribution you're working
with. Specifically:

- Silverman's rule minimises MISE under the assumption the underlying
  density is Gaussian. Standard reference: Silverman 1986. The R
  default `bw.nrd0` implements it.
- For non-Gaussian densities, and **especially multimodal or skewed
  densities**, Silverman oversmooths. This is well documented:
  "If f₀ is multimodal, Silverman's rule of thumb is known to
  oversmooth the density estimate" (Nielsen, CSwR).
- For a distribution with a **sharp peak at the true breakpoint plus
  a right-skewed recombinant tail** — which is exactly what your
  per-carrier fragment boundaries look like on the left side of an
  inversion — Silverman will:
  (a) choose a bandwidth dominated by the tail's spread
  (b) smear the sharp mode peak across several kb
  (c) shift the KDE mode interior to the true breakpoint (away from
      the sharp peak, toward the tail mean)
  (d) widen the bootstrap CI because each resample's KDE is over-
      smoothed the same way

The existing `METHODOLOGY.md §3.6` acknowledges this partially — it
restricts bandwidth computation to the 5–95 percentile range ("avoiding
the recombinant tail biasing the bandwidth"). That is a half-fix.
The bandwidth is still derived as if the core mass were Gaussian, which
it isn't — it's the concentrated peak of carriers whose ancestors
never recombined, and the peak is **sharper** than Gaussian on the
scale of SNP spacing.

## Three alternatives that fit your use case

### Option A — Half-sample mode (HSM, Bickel & Frühwirth 2006)

The half-sample mode is bandwidth-free, has breakdown point 0.5, is
low-bias on asymmetric distributions, and is specifically designed for
"the sharp-peak-with-outlier-tail" case you have.

**Algorithm** (iterative):
1. Sort the boundary positions.
2. Find the shortest interval containing ⌈n/2⌉ points.
3. Recurse on that interval (now with n' ≈ n/2 points).
4. Stop when down to 1–2 points; return the midpoint.

**R implementation:** `modeest::hsm()` (CRAN). Also available as
`genefilter::half.range.mode`.

**Properties on your distribution:**
- If 90% of carriers have fragments extending to the true breakpoint
  (the sharp peak) and 10% have shorter fragments (the recombinant
  tail), HSM will converge to the peak because that's where points
  are densest.
- No bandwidth, so no oversmoothing.
- Robust to the number of recombinant-tail carriers (breakdown 0.5).
- Fast: O(n log n) once sorted.

**Downside:** HSM converges on the densest half of the sample
by construction. If the actual mode is very narrow (say, 10 carriers
out of 166 all vote exactly at bp 15,115,243 while the rest are
scattered), HSM's interval still contains ~83 points and its centre
is pulled toward the bulk. This matches the biology (breakpoint
estimate is where dense support exists, not where a handful of
carriers are identical), so it's usually a feature, not a bug. But it
is not the same estimator as KDE-mode.

### Option B — Shorth midpoint (Grübel 1988, or LMS location)

The shorth is the midpoint of the shortest interval containing half
the sample. Single-pass, O(n log n). Not iterative like HSM.

For asymmetric distributions, the shorth is more biased than HSM
(Bickel 2002a §3.4), so **prefer HSM over shorth** for your case.
Mentioning only for completeness.

### Option C — Grenander's mode estimator (Grenander 1965)

Uses a bandwidth-like parameter k (number of closest pairs to average)
but the mode is the weighted average of interval midpoints:

  M̂ = Σ ((x_{i+k} + x_i)/2) × (x_{i+k} - x_i)^{-p}  /  normaliser

with p, k chosen such that k > 2p. Lower bias than HSM for clean
unimodal distributions, but less robust to outliers. **Not recommended
for your case** because the recombinant tail IS an outlier pattern
in Grenander's sense.

### Option D — Keep KDE but fix the bandwidth

If you want to keep the KDE framework (mode of the estimated density
function), switch bandwidth selector:

- **Sheather-Jones plug-in** (`bw.SJ` in R). Plug-in methods estimate
  the curvature of the true density from the data rather than assuming
  Gaussian. Published as the best-performing rule for multimodal and
  skewed distributions.
- **Least-squares cross-validation** (`bw.ucv`). Data-driven, no
  distribution assumption. Slower than SJ but equally principled.

SJ on your use case:
- Correctly picks a smaller bandwidth near the sharp peak (lower local
  curvature assumed by Silverman → correction downward).
- Still has the general KDE disadvantage that the mode can shift
  based on kernel shape, and that the estimator is continuous (good
  for smoothing, bad when the underlying ECDF has a near-step at the
  true breakpoint).

## Recommendation

**Primary estimator: switch to HSM.**

It matches the biology directly:
- Ancestral fragments whose boundaries cluster at the true breakpoint
  → densest half of the data → HSM converges to them.
- Recombinant-tail fragments → by definition not in the densest half
  once n > 20 carriers are in the tail-free region → HSM is not
  shifted by them.
- No bandwidth to choose or defend to reviewers.
- Published in Computational Statistics & Data Analysis (peer-reviewed,
  ~800 citations). Will not be questioned.

**Report the KDE mode as a secondary estimator** for sensitivity check.
If HSM and KDE-mode disagree by >CI width, flag it. In practice they
usually agree closely when the distribution is clean, and disagree
(HSM being correct) when the distribution has substantial recombinant
mass.

**CI: bootstrap-on-HSM** instead of bootstrap-on-KDE-mode. The
computation is identical (resample carriers, recompute HSM, take 2.5–
97.5 percentile of bootstrap HSM distribution). HSM itself is
O(n log n), same as or faster than KDE-mode, so 1000 bootstraps cost
the same.

## The concrete patch to `02_ancestral_fragments.R`

Two changes to `summarize_fragment_side()`:

1. Replace Silverman-bandwidth KDE with `modeest::hsm()`.
2. Keep KDE as a secondary estimator for reporting.

See the patched version in `scripts/02_ancestral_fragments_hsm_patch.R`.

The patch:
- Adds `modeest` as a dependency (install.packages("modeest") —
  already on CRAN, no compilation).
- Reports both `mode_hsm_bp` (primary) and `mode_kde_bp` (secondary).
- Reports `mode_agreement_kb` = |HSM − KDE| / 1000.
- Bootstrap uses HSM.
- All existing fields remain for backwards compatibility.

## Why this is an upgrade, not a deletion

Reading the methodology doc honestly:

> §3.6 … The modal position of this distribution estimates the true
> breakpoint. The spread estimates uncertainty. The tail quantifies
> recombinant history.

The *principle* is right. The *implementation choice* (Silverman-bandwidth
KDE) is suboptimal for the distribution shape your biology produces.
Switching to HSM sharpens the mode estimate while preserving the spread
and tail diagnostics (those come from the raw boundary positions, not
from the KDE, so they're independent of this change).

The KDE was not wrong to the point of breaking LG28 — §5 of METHODOLOGY.md
says "final_left_bp within 30 kb of 15.115 Mb." That's a defensible
precision. But the bootstrap CI will be narrower and the mode more
stable on complex cases (like your LG28 interior with its banding) if
you switch to HSM.

## Other small issues in the current implementation

While reading the code I noticed two minor issues worth mentioning:

### Issue 1 — `bw.nrd0` on core mass is double-truncation

Lines 154–156:
```r
q_lo <- quantile(boundary_positions, 0.05)
q_hi <- quantile(boundary_positions, 0.95)
core_mass <- boundary_positions[boundary_positions >= q_lo & boundary_positions <= q_hi]
bw <- bw.nrd0(core_mass)
```

You trim to 5–95%, then call `bw.nrd0` which internally uses
`min(IQR/1.34, sd)`. The IQR step already robustifies against tails.
Trimming first then taking IQR-based bandwidth slightly over-shrinks
the bandwidth (you've already removed the tails). Not fatal, but
inconsistent. If you keep KDE as secondary, either trim OR use
`bw.nrd0`, not both.

### Issue 2 — bootstrap mode estimator uses the same bandwidth

Lines 175–179:
```r
bw <- <computed above from full sample's core mass>
for (b in 1:bootstrap_reps) {
  bs <- sample(boundary_positions, n, replace = TRUE)
  boot_modes[b] <- kde_fit(bs)  # uses `bw` from outer scope
}
```

The bandwidth is computed once from the original sample, then re-used
across all bootstrap resamples. This is a legitimate choice (keeps
bandwidth constant as a "feature of the data"), but it means the
bootstrap CI reflects only sampling variability of the mode under a
fixed bandwidth — it does NOT include bandwidth uncertainty. The CI
is narrower than a "full" bootstrap that recomputes bw on each resample.

For HSM this issue disappears entirely (no bandwidth).

For KDE secondary reporting, either document this choice explicitly
in the output as "fixed-bandwidth bootstrap CI" or recompute bw per
bootstrap (slower but more honest CI).

## A note on LG28 specifically

From your METHODOLOGY.md §5 expected outcomes:

> 1. `final_left_bp` within 30 kb of 15.115 Mb
> 2. `final_right_bp` within 30 kb of 18.005 Mb
> 3. Left CI width < 100 kb; right CI width < 100 kb

30 kb is a defensible mode-estimation accuracy but it's larger than
the ~5 kb you should achieve with SNP-resolution fragment boundaries
(mean SNP spacing ~230 bp; ~20 informative markers near the boundary
means a natural precision on the order of a few kb if the peak is
sharp).

Switching to HSM is one of the factors that should tighten that to
~5–15 kb mode accuracy with correspondingly narrower bootstrap CI.
I'd predict LG28 HSM mode at ~15,115,000 ± 5 kb on the left and
~18,005,000 ± 10 kb on the right, based on the shelf boundaries
you already have confirmed.

When HPC is back, run the patched version once and compare to the
existing KDE output to validate. If HSM mode agrees with the
KDE mode to within kilobases, great — both methods confirm each
other. If they disagree by 20+ kb, the KDE was oversmoothed and
HSM is the correct answer.

## One caveat I want to flag

HSM has a subtle property: it converges to the densest HALF of the
distribution. If you have an inversion where 40% of carriers are
recombinant and 60% are clean (so the clean peak is barely above half
the sample), HSM's first iteration will cover both the peak and part
of the recombinant tail, then narrow to the densest sub-region. This
is fine for clean peaks but can bias slightly interior when the
recombinant fraction approaches 50%.

In those cases the v7 `complex_rearrangement_out_of_scope` gate
already fires (threshold 70% recombinant) OR `deep_interior_recombinants`
(40–70% threshold). So the rare scenario where HSM might misbehave is
already flagged.

## Reference list (for the manuscript)

- Silverman BW (1986). *Density Estimation for Statistics and Data
  Analysis*. Chapman & Hall.
- Bickel DR (2002). Robust estimators of the mode and skewness of
  continuous data. *Computational Statistics & Data Analysis* 39:
  153–163.
- Bickel DR, Frühwirth R (2006). On a fast, robust estimator of the
  mode: comparisons to other robust estimators with applications.
  *Computational Statistics & Data Analysis* 50: 3500–3530.
- Sheather SJ, Jones MC (1991). A reliable data-based bandwidth
  selection method for kernel density estimation. *JRSSB* 53: 683–690.

The Bickel & Frühwirth paper is the one to cite in your Methods if
you switch to HSM.
