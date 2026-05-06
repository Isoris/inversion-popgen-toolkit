# S3 — Karyotype-stratified Bayesian breakpoint evidence scoring

**Status:** spec only. Future module — explicitly do not implement
this turn (your instruction).
**Used by:** P5.1 Section 8 (read-evidence clusters) when this
module ships.

---

## What this module is

A read-evidence scoring framework for breakpoint candidates in the
inversion atlas. Operates on **markdup BAMs**, not VCFs. For each
candidate boundary zone, asks:

> Which H_class / karyotype group carries this breakpoint-like
> signal — and how strongly?

Output: a posterior probability per (cluster × H_class) plus a label
from a controlled vocabulary.

## Why Bayesian, why beta-binomial

Read counts vary wildly with coverage — a 30× sample contributes
more raw evidence than a 10× sample. Naive raw-count comparison
between H_classes inflates whichever class happens to contain the
deeper-coverage samples.

Beta-binomial framing solves this:
- `y_g` = breakpoint-supporting reads in group g
- `c_g` = local coverage denominator in group g (capped per fish)
- `p_g ~ Beta(α₀, β₀)` prior
- posterior: `p_g | data ~ Beta(α₀ + y_g, β₀ + c_g − y_g)`

Reports a posterior **rate** per group, not a raw count, and gives
credible intervals so downstream comparisons are honest.

## Inputs

### Per candidate boundary zone

For each candidate `cand` and each boundary side `side ∈ {left, right}`:

- A genomic window around the proposed boundary, e.g.
  `[boundary_bp − Δ, boundary_bp + Δ]` with Δ ≈ 1-5 kb.
- The candidate's `H_class` assignment per sample (from
  `candidate_group_membership.tsv` — see spec S1).

### Per sample, per cluster

From a markdup BAM aligned to the reference, for each sample and
each evidence cluster within the boundary zone:

| field | description |
|---|---|
| `n_split_support` | reads with SA (supplementary alignment) tag landing across the boundary |
| `n_pair_support` | discordant mate pairs with one mate inside, one outside, in unexpected orientation |
| `n_clip_support` | soft-clipped reads with clip position at the boundary |
| `n_mapq0_support` | reads with MAPQ=0 in the cluster (ambiguity proxy) |
| `cov_local` | local coverage (mean depth within the cluster window) — the **denominator** |

`cov_local` is the per-sample denominator. We cap it at a reasonable
max (e.g. 50× truncation) so one ultra-deep sample doesn't dominate
the H_class average.

**Cluster** = a set of read evidence linked by being within ε bp of
each other and / or by sharing SA / MATEID tags. The clustering
step merges co-located evidence into one event so we don't
double-count one breakpoint as multiple "supports".

## Model

For each cluster `k` and H_class `g`, with samples
`{1, …, n_g}`:

- Sufficient statistics:
  - `Y_{k,g} = Σ_{i in g} y_{k,g,i}` (total breakpoint-like reads)
  - `C_{k,g} = Σ_{i in g} cov_local_i` (total local coverage,
    capped per sample)

- Posterior on the rate `p_{k,g}`:
  - `p_{k,g} | data ~ Beta(α₀ + Y_{k,g}, β₀ + C_{k,g} − Y_{k,g})`

- Prior: weakly informative `Beta(α₀=1, β₀=99)` — assumes default
  rate of breakpoint reads is low (~1%). Can be tuned per study.

## Reported summaries

Per (cluster × H_class):

- `posterior_mean` = `(α + Y) / (α + β + C)`
- `cred_interval_95` = `[Beta.ppf(0.025), Beta.ppf(0.975)]`

Per cluster (across H_classes):

- `P(top_group > bottom_group)` — Monte Carlo from the two
  posteriors
- `enrichment_ratio` = `posterior_mean(top) / posterior_mean(bottom)`
- For 3-class (H1/H1, H1/H2, H2/H2):
  - `P(intermediate)` = `P(p_HET between p_H1H1 and p_H2H2)`
  - This tests the canonical "het carriers are intermediate"
    expectation. If it's high, the cluster behaves like a true
    breakpoint marker.

## Pattern labels

The cluster gets one of these labels:

| label | when |
|---|---|
| `breakpoint_resolved_candidate` | One H_class strongly carries the cluster (`P(top > bottom) > 0.99`), het is intermediate (`P(intermediate) > 0.95`), and split-read or assembly evidence is present (Tier-1 SV from spec S1) |
| `breakpoint_supported_interval` | Same as above but only Tier-2/3 evidence (no split-read confirmation yet) |
| `boundary_support_signal` | Cluster has H_class enrichment but no single-class dominance — the boundary is real but the cluster spans it ambiguously |
| `linked_internal_marker` | Cluster has H_class enrichment but is INSIDE the candidate body, not at the boundary — it's cargo, not a breakpoint |
| `mapq0_repeat_ambiguity` | Cluster's signal is dominated by MAPQ0 reads — likely repeat-mediated, not a real breakpoint |
| `family_artifact` | Cluster correlates with hatchery family rather than H_class — exclude |
| `unresolved_noise` | No H_class enrichment, no sharp signal — default fallback |

The labels match those in P5.1 Section 10 (final classification)
and feed into the candidate's overall tier.

## Pipeline shape (future)

A standalone script `STEP_R45_breakpoint_evidence.py` (or R
equivalent) that:

1. For each candidate in `candidate_catalogue.tsv`:
2. For each boundary side (left, right):
3. Pull markdup BAM slices for all 226 samples in
   `[boundary_bp − Δ, boundary_bp + Δ]`.
4. Cluster reads by SA tag / MATEID / position-proximity.
5. Per cluster, count `n_split_support`, `n_pair_support`,
   `n_clip_support`, `n_mapq0_support` per sample.
6. Compute `cov_local` per sample (capped).
7. Join with `candidate_group_membership.tsv` to get H_class.
8. Aggregate to (cluster × H_class) sufficient statistics.
9. Compute Beta posteriors and summaries.
10. Apply pattern label rules.
11. Emit `breakpoint_evidence_per_cluster.tsv` and
    `breakpoint_evidence_per_candidate.tsv` (rolled-up to
    candidate level — best label across clusters).

## Atlas-side rendering (when this lands)

In P5.1 Section 8 (read-evidence clusters), the renderer expects
the per-cluster TSV. For each cluster:

- Box plot of posterior mean per H_class with credible intervals
- Pattern label as a chip
- Drill-down: per-sample raw counts (audit trail)

## Why this is spec-only

You explicitly said:
> Future read-evidence algorithm can be a spec or placeholder
> unless we explicitly ask to code it now.

It needs:
1. Markdup BAMs locally (currently only on LANTA — not local).
2. The clustering algorithm tuned per evidence type.
3. The β prior calibrated against a known truth set.
4. Validation on at least one resolved candidate (e.g. LG28 shelf)
   before claiming to score new candidates.

None of these are blockers for the immediate UI / data-bridge work.
The placeholder section in P5.1 Section 8 is enough scaffolding for
when this ships.

## Risk notes

- **Coverage cap matters a lot.** Without it, a single 60× outlier
  sample makes its H_class look like it has more breakpoint
  evidence. Cap at 30× or 50× per sample before summing into `C_g`.
- **Family confounding** (same as S2) — apply family correction or
  stratify by family before claiming H_class enrichment.
- **MAPQ0 is context, not signal.** `n_mapq0_support` enters the
  model as a separate cluster type or as a confounder, not as raw
  breakpoint evidence. The pattern-label `mapq0_repeat_ambiguity`
  fires when MAPQ0 dominates the cluster.
- **Don't over-claim.** A `breakpoint_supported_interval` label
  means "the data is consistent with a breakpoint here at the
  resolution of split reads" — not "the breakpoint is at coordinate
  X". The actual base-pair coordinate still needs the assembly
  junction or long-read confirmation.
