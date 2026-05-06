# U-shape evolution — class rules (plain language)

This module classifies inversion candidates by the **shape of their
between-arrangement divergence profile**, not by their causal mechanism.

The phrasing matters. "neutral-like" ≠ "neutral". "Locally-adapted-like" ≠
"adapted". The shape is the verdict; the cause is a separate question that
needs additional data we usually don't have.

## What the shape data is

For each candidate, the module computes per-window dXY, FST, π_HOMO1 and
π_HOMO2 across the inversion plus flanking regions, then summarises the
profile with a small set of scores:

- **dXY inside/flank ratio** — how much higher between-arrangement
  divergence is *inside* the inversion vs the surrounding genome.
- **U-score** — edge dXY / center dXY. High when both edges are taller
  than the center (the "suspension bridge" pattern of Berdan et al. /
  Guerrero et al.).
- **internal peak score** — max dXY in the center / mean dXY at edges.
  High when the center has a sharp peak that dominates the edge signal.
- **flatness score** — sd(inside dXY) / mean(inside dXY). Low = flat.
- **asymmetry** — left-edge dXY vs right-edge dXY.
- **oldness score** — composite of inside/flank ratios + private SNP
  burden. Acts as the **age gate**: shape statistics are uninformative
  on young inversions (Berdan/Guerrero figure 2, top row), so the
  classifier refuses to assign a shape verdict below threshold.

## The rules

The classifier is a **two-step gate-then-shape**:

1. If the candidate fails the age gate (oldness too low, or dXY inside
   too low, or fewer than ~8 inside windows), it gets
   `young_weak_divergence` or `insufficient_data` and stops there.
2. Otherwise the candidate enters one of the shape classes below.

### 1. neutral_like_U_shape

**When:** dXY higher near both edges than in the center.
Main score: `dxy_u_score = mean(edge dXY) / mean(center dXY)`.

**Reading:** consistent with long-suppressed recombination producing a
U-shaped coalescence-time profile. **Not proof of neutrality** — the same
pattern can arise from epistatic balancing selection or any long-lived
inversion (Berdan et al. caption explicitly warns about this).

### 2. locally_adapted_internal_peak_like

**When:** dXY or FST has an internal peak that dominates edge signal.
Main scores: `dxy_internal_peak_score`, `fst_internal_peak_score`.

**Reading:** consistent with selected alleles or a differentiated
functional region inside the inversion (Berdan/Guerrero figure 2c).
**Not proof of adaptation** — peaks can be drift-driven in small
populations, hatchery family structure can mimic them, and the test
doesn't tell us *which* alleles matter.

### 3. locally_adapted_breakpoint_like

**When:** strong edge signal in dXY/FST without a corresponding
center peak. Main score: `breakpoint_adaptation_score`.

**Reading:** consistent with breakpoint-proximal selection or
breakpoint-linked divergence. **Not proof of breakpoint mechanism** —
this is a profile signature, not a molecular finding.

### 4. flat_high_deep_structural_haplotype

**When:** dXY is broadly elevated across the inversion with low
flatness (small variance). Looks "flat-deep" rather than U-shaped.

**Reading:** consistent with long-term suppressed recombination over a
deep/old structural haplotype that has accumulated divergence
everywhere. Reviewer-safe statement; does not commit to a mechanism.

### 5. young_weak_divergence

**When:** local-PCA grouping exists but dXY/FST is weak everywhere
and/or oldness is below threshold.

**Reading:** young inversion, weakly diverged arrangement, or the
candidate may not be a real structural polymorphism. The shape
classifier deliberately stops here — there isn't enough divergence to
read a profile from.

### 6. asymmetric_edge

**When:** one edge is much hotter than the other
(`abs_log2_asymmetry ≥ ~1.0`).

**Reading:** one boundary may be sharper, older or more differentiated;
or the candidate may be misbounded (a wide call with one true edge and
one drifted edge); or one breakpoint sits on different architecture from
the other. **Worth manual inspection** — the call is honest about the
asymmetry rather than pretending it's a clean profile.

### 7. complex_mixed

**When:** multiple shape rules fire simultaneously above threshold.

**Reading:** nested or overlapping inversions, recombination/gene
conversion smearing the profile, or genuinely mixed history. Don't
over-interpret. The atlas exposes which signals fired in `flags`.

### 8. insufficient_data

**When:** too few inside windows, ratios are NaN, or HOMO_1/HOMO_2
sample sizes are below the minimum.

**Reading:** the module declines to call. Get more data or revisit the
karyotype assignment.

## What this layer is NOT

- It is **not** an adaptation test. It cannot answer "is this inversion
  adaptive?". Use cargo annotation + deleterious load (the breeding-relevant
  core) for the main paper claims.
- It does **not** include any TE / repeat / SD / mappability content.
  Breakpoint architecture lives in the comparative-fragility module
  (`phase_8_comparative_breakpoint_fragility/`); this is intentional
  — they answer different questions and shouldn't be conflated.
- It does **not** confirm karyotype groups. The module trusts the input
  groups (HOMO_1 / HOMO_2 / HET) and classifies the resulting shape.
  Wrong groups → wrong shape; the upstream PCA is the source of truth.

## Safe phrasings for the manuscript

Use these, in order of preference:

> "Some candidates show a U-shaped between-arrangement divergence
> profile consistent with long-suppressed recombination."

> "Several candidates display an internal peak in differentiation
> consistent with locally-adapted alleles, although our hatchery design
> cannot test adaptation directly."

> "The inversion atlas catalogues structural haplotypes and their
> arrangement-specific deleterious load; evolutionary-shape
> classification is offered as a supplementary descriptor, not as
> evidence of selection."

Avoid:
- "This inversion is adaptive."
- "This inversion is neutral."
- "Recombination is suppressed inside this inversion." (true *if* it's
  a real heterokaryotype-segregating inversion, but that's a separate
  proof.)
