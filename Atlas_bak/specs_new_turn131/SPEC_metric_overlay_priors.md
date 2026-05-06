# SPEC — Metric overlay priors (secondary boundary detectors)

**Status**: drafted turn 130 final session. Most metrics already
computed by upstream LANTA modules. This spec defines how to **surface
them in the atlas** as overlay tracks aligned with the boundary-consensus
aggregator (`SPEC_boundary_consensus_aggregator.md`).

**Trigger** (Quentin, chat `5b793a68`):
> *"From the flashlight system and our discussions, the metrics for
> Engine B boundary detection are: θ_pi / θ_w (Module 3), heterozygosity
> per sample, Fst between PC1 bands (Module 2B/2C), Admixture Q shift,
> Dosage variance, DELLY/Manta breakpoint positions, BND triangulation,
> read depth dips/spikes, GHSL within-sample haplotype contrast, sim_mat
> insulation drop, squareness ratio, eigenvalue λ1/λ2 transitions, PC1
> band composition switch, inv_likeness step function, Tajima's D sign
> change."*
>
> *"The iteration: Pass 1 finds boundaries from genotype patterns. Pass 2
> confirms/adds from metrics. Pass 3 re-resolves with updated boundary
> set. Repeat until stable. Each pass peels one layer of nesting."*

Each metric is **independently noisy** but **convergent at real
structural boundaries**. The atlas overlays them visually so the user
can see (or an aggregator can compute) where signals agree.

---

## 1. The metric inventory

### 1.1 Primary boundary signals (already feeds aggregator)

| Metric | Source script | What changes at a boundary |
|---|---|---|
| θπ-CUSUM | `STEP_T05_theta_cusum.R` | per-carrier diversity changepoint |
| Dosage changepoint | `STEP_C01h_recombinant_scanner.R` | per-carrier dosage step |
| L3 contingency-table | (planned L3 step) | adjacent-L2 karyotype concordance |

### 1.2 Secondary (independent) boundary signals

| Metric | Source | What changes at a boundary |
|---|---|---|
| **Fst step** | Engine B `region_popstats.c` | Fst spikes between groups across a real boundary |
| **dXY step** | Engine B | absolute divergence steps |
| **Tajima's D sign change** | Engine B `theta_pi` + thetaStat | swept ↔ neutral region transitions |
| **Heterozygosity drop** | Module 3 | het rate steps at inversion edges |
| **Q shift** | Engine A NGSadmix | local ancestry composition changes |
| **GHSL run boundary** | snake 3 | within-sample haplotype contrast switch |
| **sim_mat insulation drop** | C01a precomp | ±20/40/80/160/320 multi-scale boundary |
| **Squareness ratio** | C01c triangles | far/near similarity geometry |
| **λ₂ transition** | C01a | eigenvalue 2 changes between clean and complex |
| **inv_likeness step** | C01a | drops at boundaries |
| **DELLY/Manta breakpoints** | Module 4 | hard breakpoint coordinates with CIPOS/CIEND |
| **Read depth dips** | Module 4 / cheat 10 | depth signature at breakpoints |
| **Clipped-read pileup** | Module 4 / cheat 11 | soft-clip enrichment |
| **Mutual information (MI)** | `region_popstats.c` (added per chat `260860f9`) | pairwise informativeness step |
| **Cumulative indel slope** | Cheat 9 (Clair3) | het_slope/hom_slope changes at boundaries |

## 2. Atlas surface — overlay tracks

Each metric becomes an optional track in the **boundary-consensus view**
(see `SPEC_boundary_consensus_aggregator.md`):

```
[Stacked KDE of primary streams] ────────────────────────
   θπ-CUSUM observations ▮▮▮▮  ▮▮▮▮▮▮▮▮▮▮▮▮▮▮ (single-mode at 5')
   Dosage observations    ▮▮▮▮▮▮▮▮▮  ▮▮▮▮▮▮  (matches)
   Contingency calls      ▮       ▮            (sparse)
   ───────────────────────────────────────────────────────
   Mode peak: 15.115 Mb   (from KDE)

[Secondary metric overlay strip — one row per metric, colored where
 signal exceeds threshold]
   Fst (band1↔band3)     ─────────█████─────████──────  steps at 15.1 + 17.99
   Tajima's D             ─────────────████─────────────  sign change at ~16
   Heterozygosity         ─────────███─────────████─────  drops at 15 + 18
   Q shift (K=8)          ─────────█████──────────────  shift at 15.1
   GHSL boundary          ─────────█████─────████──────  matches Fst
   sim_mat ±80            ─────────█████─────████──────  geometric
   ...
   ───────────────────────────────────────────────────────
   Convergent signal: 15.115 Mb (8/12 metrics) → HIGH CONFIDENCE
   Convergent signal: 17.985 Mb (7/12 metrics) → HIGH CONFIDENCE
```

Click any metric → opens that metric's per-window plot beneath.

## 3. Signal aggregation into "convergent boundary call"

For each candidate boundary (already detected by the primary aggregator),
score how many secondary metrics **support** it:

```
For each metric M in secondary list:
  if M has a discontinuity within ±20 kb of the candidate boundary:
    convergence_count += 1
    contributing_metrics.append(M)
```

Output:
- **HIGH CONFIDENCE** boundary: ≥6 of 12 secondary metrics support
- **MODERATE**: 3-5
- **LOW**: 1-2 — flag for review (single-stream may not be real)
- **DISPUTED**: 0 secondary metrics — boundary may be artifact

This is **not** a hard threshold; it's a confidence label.

## 4. Iterative refinement (the Pass 1/2/3 logic)

Per chat `5b793a68`, the pass logic:

**Pass 1**: Primary streams (θπ-CUSUM + dosage + contingency) call
boundaries via consensus aggregator.

**Pass 2**: Secondary metrics confirm or add. If a secondary metric
shows a step at a position the primary streams missed, add it as a
candidate boundary for Pass 3.

**Pass 3**: Re-run consensus aggregator including any new candidate
boundaries. Iterate until stable (no new boundaries added in a pass).

Atlas-side: this could be implemented as a single-pass overlay (compute
once, render). Iterative refinement is a LANTA-side concern; the atlas
just renders the final set of boundary calls + their convergence
scores.

## 5. Implementation slices

### Slice 1 — overlay strip rendering (~0.5 turn)
- Per-metric horizontal strip in the boundary view
- Heat-coding for signal strength (z-score thresholded)
- Click handler to open per-metric detail

### Slice 2 — convergence count (~0.3 turn)
- Per-boundary convergence aggregator
- Confidence label (HIGH/MOD/LOW/DISPUTED)
- Tooltip listing contributing metrics

### Slice 3 — per-metric detail panel (~0.5 turn)
- When user clicks a metric, render its per-window plot beneath the
  consensus view (e.g., Fst track, Tajima's D track)

## 6. Open questions

1. **Threshold per metric**: each metric has different scale (Fst 0–1,
   λ₂ unbounded, depth 0–60). Per-metric z-score threshold (e.g., |z|≥2)
   is the simplest standardization.
2. **Smoothing window**: each metric needs different smoothing.
   Heterozygosity is per-SNP and very noisy; θπ is per-window and
   already smooth. Defer to per-metric defaults from upstream.
3. **What's a "real" discontinuity vs noise?** This is the calibration
   challenge. Chat `5b793a68` discussed using boundary-relative null
   distributions. Atlas-side, just use z-score thresholds; LANTA-side
   producers can pre-flag confident-step positions.

## 7. Tests

- Synthetic candidate with planted Fst step at 15.1 Mb → overlay
  shows step at correct position.
- Convergence count: 8 metrics step at the same position → HIGH
  CONFIDENCE label.
- Disputed boundary: only 1 metric supports → LOW label with review
  note.

## 8. Dependencies

- Engine B `region_stats_dispatcher.R` for Fst, dXY, θπ, Tajima's D,
  Q (existing per chat `260860f9`, `53db009e`).
- Module 3 het track (existing).
- Module 4 SV calls (existing).
- GHSL panel (existing).
- C01a precomp for sim_mat / λ / inv_likeness (existing).
- `SPEC_boundary_consensus_aggregator.md` Stage 2 output (the
  candidate boundaries to confirm).

## 9. What this is NOT

- Not a primary boundary detector. It supplements the aggregator.
- Not a recombinant-detection method. That's
  `SPEC_recombinant_dosage_changepoint_detector.md`.
- Not a method paper in itself. It's an aggregation strategy that uses
  existing metrics; the novelty is in the multi-stream consensus, not
  in any individual metric.

## 10. Cross-references

- `SPEC_boundary_consensus_aggregator.md` — primary aggregator that
  this spec extends with secondary signals.
- Chat `5b793a68` — original metric inventory.
- Chat `41c41d1f` — the iterative refinement design.
- Chat `260860f9` — Engine B (region_popstats.c) producing many of
  these metrics.
