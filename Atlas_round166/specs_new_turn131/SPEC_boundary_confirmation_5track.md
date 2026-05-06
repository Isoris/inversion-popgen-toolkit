# SPEC — 5-track boundary confirmation composite (the convicting figure)

**Status**: drafted turn 130 final session. Referenced in chat
`79f2b4d7` as `SPEC_boundary_confirmation_5track.md` — that spec was
not in the spec folder, so this is the consolidated version. Figure F4
panel A of the manuscript depends on it.

**Trigger** (Quentin, chat `79f2b4d7`):
> *"Panel A: 5-track boundary-confirmation composite for prototype LG28
> candidate (local PCA |z|, Hudson Fst, dXY, θπ by karyotype, repeat
> composite). The figure that 'convicts' LG28 as a real inversion."*

A per-candidate composite figure that overlays five independent tracks
aligned at the candidate's boundaries. When all five tracks show
coherent signal at the same Mb position, the inversion is "convicted."

---

## 1. The 5 tracks

### Track 1 — Local PCA |z|
The original z-score envelope from the staircase detector. The signal
that identified the candidate. Shows L1/L2 boundary positions.

### Track 2 — Hudson Fst (HOM_REF vs HOM_INV)
Computed by Engine B `region_popstats.c` with the karyotype-derived
HOM_REF and HOM_INV groups as inputs. Should spike inside the
inversion, drop outside.

### Track 3 — dXY (HOM_REF vs HOM_INV)
Absolute divergence between arrangements. Different from Fst — Fst is
relative differentiation, dXY is absolute. Together they distinguish
"groups are similar but well-separated" from "groups are very different."

### Track 4 — θπ stratified by karyotype
Three sub-tracks overlaid:
- θπ in HOM_REF samples
- θπ in HET samples
- θπ in HOM_INV samples

Expected pattern: HET θπ is elevated inside the inversion (carrying both
arrangements simultaneously). HOM_REF and HOM_INV θπ are similar to
each other and lower than HET inside the inversion.

### Track 5 — Repeat composite
A density track showing repeat element density (TE classes, satellite
DNA) at the breakpoints. Inversion breakpoints frequently sit on
repeat elements (NAHR substrates) — this track tests the mechanism.

## 2. Visual layout

```
  [Mb axis: 14 ─ 15 ─ 16 ─ 17 ─ 18 ─ 19]
  ─────────────────────────────────────────
  Track 1: Local PCA |z|
    │ ▮▮▮▮▮▮▮▮ (envelope shown as filled area)
  ─────────────────────────────────────────
  Track 2: Hudson Fst (REF vs INV)
    │       ────────  (high inside, low outside)
  ─────────────────────────────────────────
  Track 3: dXY (REF vs INV)
    │       ────────
  ─────────────────────────────────────────
  Track 4: θπ by karyotype
    │         ─── HET (red, elevated)
    │       ─── REF (blue)
    │         ─── INV (green)
  ─────────────────────────────────────────
  Track 5: Repeat composite
    │ ▮     ▮ TE density spikes at breakpoints
  ─────────────────────────────────────────
  Boundary alignment: ↓ ↓ (5'-end at 15.115 Mb, 3'-end at 17.985 Mb)
```

## 3. Atlas-side surface

Atlas page-2 candidate-focus has a "boundary confirmation" sub-tab.
When all 5 tracks have data, the composite renders. Missing tracks
show a placeholder.

Each track is a thin canvas (~40 px high). Total composite ~250 px tall.
Highly compressed for screen reading; expands to full-page SVG on export.

## 4. JSON layer schema

The atlas reads several existing layers + a new one:

| Track | Source layer | Status |
|---|---|---|
| Local PCA |z| | `local_pca_z_v1` (precomp) | shipped |
| Hudson Fst | `popstats_per_window_v1` (Engine B output) | producer needed |
| dXY | `popstats_per_window_v1` | same as above |
| θπ by karyotype | `theta_pi_per_karyotype_v1` (new) | producer needed |
| Repeat composite | `te_density_v1` (atlas TEfull layer) | shipped |

## 5. Producer (LANTA-side)

```r
# emit_popstats_per_window.R
# Run region_popstats.c per window for each candidate
# Emit Hudson Fst + dXY + per-group θπ in one JSON

for (cand in candidate_list) {
  groups <- list(
    HOM_REF = read_karyotype_samples(cand, "HOM_REF"),
    HOM_INV = read_karyotype_samples(cand, "HOM_INV"),
    HET     = read_karyotype_samples(cand, "HET")
  )
  for (window in cand$windows) {
    stats <- region_popstats(beagle, window, groups)
    emit_to_layer(stats)
  }
}
```

Per-candidate JSON layer:

```json
{
  "schema": "popstats_per_window_v1",
  "per_candidate": {
    "LG28_inv5": {
      "windows": [
        { "start_bp": 14000000, "end_bp": 14010000, "fst": 0.04, "dxy": 0.012, "theta_pi": { "HOM_REF": 0.0034, "HET": 0.0089, "HOM_INV": 0.0028 } },
        ...
      ]
    }
  }
}
```

## 6. Convicting score

A summary score per candidate:

```
n_tracks_supporting = count(
  fst_step_at_5p_boundary && fst_step_at_3p_boundary,
  dxy_step_at_both_boundaries,
  het_theta_pi_elevated_inside_inversion,
  repeat_density_spike_at_boundaries,
  local_pca_z_envelope_consistent
)

verdict = (n_tracks_supporting >= 4) ? "CONVICTED" : "PARTIAL_SUPPORT" : "WEAK"
```

This score appears on the verdict panel.

## 7. Implementation slices

### Slice 1 — JSON layer detection + composite rendering (~0.7 turn)
- Detect `popstats_per_window_v1` and `theta_pi_per_karyotype_v1`
- Composite canvas with all 5 tracks
- Aligned Mb axis
- Boundary indicators

### Slice 2 — convicting score (~0.3 turn)
- Per-track signal-at-boundary check
- 4-of-5 majority verdict
- Surface in verdict panel + catalogue

### Slice 3 — manuscript SVG export (~0.5 turn)
- High-resolution SVG with manuscript typography
- Configurable color scheme
- Caption template

## 8. Open questions

1. **What's a "step" at boundary?** A z-score >2 jump in track value
   between adjacent windows (or windows within ±20 kb of the boundary).
   Calibrate on real LG28.
2. **HET θπ elevation threshold**: should be detectably above HOM
   tracks. Effect-size threshold like Cohen's d > 0.5.
3. **What if no HOM_INV samples exist (rare inversions, n=1-3)?**
   Skip Tracks 2, 3 (Fst and dXY between arrangements not computable)
   and just show Tracks 1, 4, 5. Verdict adjusted accordingly.

## 9. Tests

- Synthetic candidate with planted 5-track agreement → CONVICTED.
- Synthetic with only 2-track agreement → PARTIAL_SUPPORT.
- Real LG28 should be CONVICTED (validation target).

## 10. Cross-references

- `SPEC_boundary_consensus_aggregator.md` — uses the same per-track
  evidence but for boundary refinement, not visual composite.
- `SPEC_metric_overlay_priors.md` — broader metric inventory; this
  spec is the focused 5-track version for the manuscript figure.
- `SPEC_per_candidate_breeding_readiness_card.md` — the convicting
  score appears in the per-candidate card.
- `SPEC_manuscript_bundle_export.md` — exports the SVG for F4 panel A.
- Chat `79f2b4d7` for the original design.
- Chat `260860f9` for `region_popstats.c` (the Engine B C binary
  producing Tracks 2, 3, 4).

## 11. What this is NOT

- Not a candidate detector. It validates already-detected candidates.
- Not a replacement for the hypothesis test framework
  (`SPEC_hypothesis_test_framework_atlas.md`). The hypothesis tests
  are statistical; this composite is visual evidence.
- Not the only confirming figure. Atlas 5 (breeding card) is for
  practitioners; this 5-track composite is for the manuscript reader.
