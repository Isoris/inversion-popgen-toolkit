# SPEC — Boundary consensus aggregator

**Status**: drafted turn 130 final session. Not yet implemented.
Estimated ~1–1.5 turns once producer streams (CUSUM + contingency +
dosage) are stable. **This is the spec that fixes the "stack-then-aggregate"
problem Quentin flagged in chat `487c7f04`.**

**Trigger** (Quentin, chat `487c7f04`):
> *"Combine evidence streams — currently each stream pre-aggregates then
> merges; the right design is stack-then-aggregate, but that requires
> the redesign we deferred."*
>
> *"It will use all 3 pages CUSUMS + contingency tables + dosage. And
> identify every inversion for every chromosomes all automatic when we
> load all JSONs at once."*

This spec defines how to combine three independent boundary-detection
streams (θπ-CUSUM, contingency-table merge, dosage changepoint) into
a single boundary call per candidate, without the flaw of the existing
phase 4c `03_consensus_merge.R` (which assumes one-inversion-equals-one-start-and-end).

---

## 1. The problem

The atlas + LANTA pipeline has three boundary signals:

| Stream | What it sees | Where it lives | Strength |
|---|---|---|---|
| **θπ-CUSUM** | per-carrier diversity changepoint | LANTA `STEP_T05_theta_cusum.R` | sharp at inversion edges, noisy in low-diversity regions |
| **L3 contingency-table merge** | adjacent-L2 karyotype concordance | LANTA L3 step (planned) | recovers over-splits; can't recover under-splits |
| **Dosage changepoint** | per-carrier dosage step | LANTA `STEP_C01h_recombinant_scanner.R` | direct signal at recombinant breakpoints |

Plus several **secondary** boundary detectors that should feed this
aggregator (per `SPEC_metric_overlay_priors.md`): Fst step, GHSL
karyotype-run boundaries, sim_mat insulation drops, λ₂ transitions,
Tajima's D sign change, Q shift.

The current phase 4c logic pre-aggregates each stream's per-carrier
boundaries into a single mean+CI per inversion, then merges across
streams. That fails on hatchery LD multimodal cases where one inversion
has two distinct boundary clouds (genuinely bimodal carrier observations
because of recombinant tracts).

The fix: **stack first, aggregate after**. Pile all per-carrier
observations into one position-density distribution per inversion,
then do mode detection on the stacked distribution.

## 2. The aggregator algorithm

### 2.1 Stage 1 — stack

For inversion candidate `cid`:

```
all_obs[cid] = []
for stream in [theta_cusum, contingency, dosage_changepoint]:
  for carrier_call in stream.calls_for(cid):
    all_obs[cid].append({
      stream,
      carrier_id,
      bp,
      weight,           # stream-specific reliability score
      side,             # 5'-end or 3'-end
    })
```

The stack is a flat list of observations, not pre-aggregated.

### 2.2 Stage 2 — mode detection on stacked distribution

For each side (5' and 3' separately):

```
positions = [obs.bp for obs in all_obs[cid] if obs.side == side]
weights   = [obs.weight for obs in all_obs[cid] if obs.side == side]

kde      = gaussian_kde(positions, weights=weights)
modes    = find_local_maxima(kde)
```

If `n_modes == 1` → standard unimodal boundary call:
`bp = mode_position`, `CI = bootstrap KDE`.

If `n_modes ≥ 2` → emit each mode as a sub-candidate boundary:
- Mode 1 at bp_a → candidate A spans (bp_a_5, bp_a_3)
- Mode 2 at bp_b → candidate B spans (bp_b_5, bp_b_3)

Mode separation gate: modes within 50 kb collapse into one. Wider
separation → genuinely separate boundaries that may indicate
**two overlapping inversions** rather than one inversion with bimodal
boundaries.

### 2.3 Stage 3 — per-stream attribution

Track which streams contributed to each mode:

```
mode_attribution[mode_a] = {
  theta_cusum: 12 carriers,
  contingency: 1 boundary call,
  dosage: 8 carriers,
}
```

A mode supported by ≥2 streams → high confidence. Single-stream mode →
flag as `single_stream_only`.

## 3. Decision: collapse vs split sub-candidates

When ≥2 modes are detected, the aggregator does NOT auto-decide whether
to split into two candidates. It emits both modes with their attribution
and a **flag for human review**: `boundary_multimodal: true`.

The user (or a downstream rule) decides:
- Collapse if the modes are within typical CI (~20 kb at 9× coverage)
- Split if the modes are >100 kb apart and supported by ≥2 streams each

This is the key learning from chat `487c7f04`: don't bake in the
mode-decision before the user has seen real data.

## 4. Stream weights (calibration)

Default weights (subject to calibration on real LG28):

```
theta_cusum:       1.0  (clean signal at inversion edges)
contingency:       0.7  (works on adjacent L2s only)
dosage:            1.5  (direct, per-carrier, sharpest signal)
fst_step:          0.6  (secondary; noisy outside inversions)
ghsl_run_boundary: 0.8  (only when run length is sufficient)
sim_mat_drop:      0.5  (geometric, multi-scale; can be confounded)
lambda2_trans:     0.4  (eigenvalue, indirect)
tajimas_d_sign:    0.3  (independent confirmation, sometimes flat)
q_shift:           0.4  (admixture; depends on K resolution)
```

Calibration: re-tune after Quentin's diagnostic protocol on real LG28.

## 5. Output schema

Per candidate, the aggregator emits:

```json
{
  "cid": "LG28_inv5",
  "boundary_5": {
    "best_bp": 15115000,
    "ci_lo": 15095000,
    "ci_hi": 15135000,
    "n_modes": 1,
    "modes": [
      { "bp": 15115000, "weight": 18.4, "streams": ["theta_cusum", "dosage"], "n_carriers": 24 }
    ],
    "multimodal_flag": false
  },
  "boundary_3": {
    "best_bp": 17985000,
    "ci_lo": 17920000,
    "ci_hi": 18050000,
    "n_modes": 2,
    "modes": [
      { "bp": 17920000, "weight": 12.1, "streams": ["theta_cusum"], "n_carriers": 14 },
      { "bp": 18050000, "weight": 14.5, "streams": ["dosage", "ghsl_run_boundary"], "n_carriers": 18 }
    ],
    "multimodal_flag": true,
    "review_note": "3'-end has two modes 130 kb apart; possible recombinant cloud or two overlapping inversions."
  }
}
```

## 6. Implementation slices

### Slice 1 — stack + KDE (~0.5 turn, atlas or LANTA)
- Implement stacking pipeline reading per-stream JSON layers
- KDE-based mode detection
- Tests against synthetic per-carrier observations

### Slice 2 — multi-mode flagging (~0.3 turn)
- Mode-separation gate
- Stream-attribution computation
- Multimodal flag with review note

### Slice 3 — atlas surface (~0.5 turn)
- Per-candidate boundary view: render the stacked KDE with mode peaks
  highlighted
- Multi-mode candidates show both modes side-by-side
- Review action: collapse / split / accept

## 7. Open questions

1. **Where does this run — atlas or LANTA?** The compute is
   straightforward (KDE on a few hundred observations per side). Atlas
   side is fine if all per-stream observations are bundled into the
   chromosome JSON. LANTA-side is needed if observations are too large
   for the JSON.
2. **What's a "mode" exactly?** Local maximum of the KDE above a
   minimum threshold. Avoid spurious low-density bumps from single
   observations. Threshold = 5% of max density.
3. **Bandwidth selection**: Silverman's rule gives a default; per-stream
   bandwidth could be tuned (CUSUM is sharp, dosage even sharper, Fst
   coarser).
4. **What about 4+ modes?** Theoretically possible with multiple
   nested or overlapping inversions in the same region. Keep all modes,
   flag for human review.

## 8. Tests

- Synthetic unimodal: 50 observations centered at 15 Mb → KDE recovers
  one mode at 15 Mb.
- Synthetic bimodal: 25 obs at 15 Mb + 25 obs at 16 Mb → two modes.
- Synthetic single-stream: only theta_cusum has observations →
  `single_stream_only` flag.
- Calibration: re-run on real LG28 once data is available; verify
  modes match Quentin's eyeballed boundaries.

## 9. Dependencies

- LANTA-side `STEP_T05_theta_cusum.R` (already exists per chat
  `487c7f04`).
- LANTA-side `STEP_C01h_recombinant_scanner.R` (already exists per
  chat `b7bc2608`).
- L3 contingency-table step (planned per
  `SPEC_l2_sweep_inheritance.md` and chat `92caef41`).
- Per-stream JSON layer emission (consumer side; defined by each
  producer's JSON contract).

## 10. What this is NOT

- Not a candidate detector. It refines boundaries on already-detected
  candidates.
- Not a final boundary call. Multi-mode candidates need human review
  before the boundary becomes the manuscript number.
- Not a replacement for biological judgment. The aggregator says
  "two modes 130 kb apart"; the user decides if that means "two
  inversions" or "one inversion with recombinant cloud."

## 11. Cross-references

- `SPEC_metric_overlay_priors.md` — secondary boundary detectors that
  feed this aggregator.
- `SPEC_recombinant_dosage_changepoint_detector.md` — one of the three
  primary streams.
- `SPEC_karyotype_per_interval_intersection.md` — provides per-interval
  karyotype calls that the contingency-table stream consumes.
- Chat `487c7f04` — the design conversation that produced this spec.
- Chat `5b793a68` — the metric-overlay list (Engine B).
