# SPEC — Age/origin atlas panel (cheat30 GDS surface on candidate focus)

**Status**: drafted 2026-05-03, scope = page 3 (DOM `page2`) candidate focus.
**Sibling specs**: `SPEC_inversion_age_atlas_surface.md` (between-inversion
ranking via dXY → time, manuscript-grade Task A + Task B). This spec is
distinct: cheat30 answers a *different question* about each candidate —
"is the inversion real, single-origin or recurrent, and roughly how
diverged are the arrangements?" — using genotype-stratified pairwise GDS
(Porubsky et al. 2022 Cell, Fig 3C/D).

## 1. What cheat30 produces

Source: `cheat30_gds_by_genotype.R` (Quentin's HPC pipeline, runs on the
BEAGLE dosage matrix per candidate). Per-candidate output structure
(directly from `run_cheat30()` return value):

```
{
  candidate_id: "LG28_cand_01",
  chrom: "LG28",
  start_bp: 15_010_000,
  end_bp: 18_240_000,
  n_ref:   60,                      # HOM_REF sample count
  n_het:  106,                      # HET sample count
  n_inv:   60,                      # HOM_INV sample count
  separation_p:        4.8e-12,     # Wilcoxon: same-genotype GDS > diff-genotype GDS
  separation_effect:   0.082,       # mean(GDS_same) - mean(GDS_diff)
  age_proxy:           0.067,       # mean(GDS_REF/REF) - mean(GDS_REF/INV); larger = older
  dip_stat:            0.041,       # Hartigan's dip on GDS_INV/INV distribution
  dip_p:               0.012,       # dip test p-value
  is_bimodal:          true,        # I/I bimodal → recurrent
  origin_class:        "recurrent", # "single_origin" / "recurrent" / "weak_signal" / "inconclusive"
  mean_ibs_same:       0.847,       # mean GDS within REF/REF + within INV/INV
  mean_ibs_diff:       0.765,       # mean GDS across REF/INV
  pair_table:          [ ... ]      # optional per-pair rows for diagnostic plots
}
```

The `pair_table` is the data needed for the Porubsky-style violin /
ridgeline plot (one density curve per genotype-pair type:
`HOM_REF_HOM_REF`, `HOM_REF_HOM_INV`, `HOM_INV_HOM_INV`, optionally
`HET_HET`, `HOM_REF_HET`, `HET_HOM_INV`).

## 2. Wire format — `cheat30_gds_results.json`

The HPC writes one JSON file per chromosome (matching the ncRNA / repeat
density convention). Atlas accepts it via the existing drop-zone or via a
new "+ load enrichment" button (whichever is easier).

```json
{
  "schema_version": "cheat30_v1",
  "generated_at": "2026-05-15T03:14:22Z",
  "species": "Clarias gariepinus (Gar subgenome)",
  "chrom": "LG28",
  "candidates": {
    "LG28_cand_01": {
      "n_ref": 60, "n_het": 106, "n_inv": 60,
      "separation_p": 4.8e-12,
      "separation_effect": 0.082,
      "age_proxy": 0.067,
      "dip_stat": 0.041, "dip_p": 0.012, "is_bimodal": true,
      "origin_class": "recurrent",
      "mean_ibs_same": 0.847,
      "mean_ibs_diff": 0.765,
      "pair_summaries": {
        "HOM_REF_HOM_REF": { "n": 1770, "mean": 0.851, "sd": 0.022,
                             "quantiles": [0.78, 0.81, 0.85, 0.87, 0.91] },
        "HOM_INV_HOM_INV": { "n": 1770, "mean": 0.842, "sd": 0.041,
                             "quantiles": [0.72, 0.79, 0.84, 0.87, 0.93] },
        "HOM_REF_HOM_INV": { "n": 3600, "mean": 0.765, "sd": 0.018,
                             "quantiles": [0.71, 0.74, 0.77, 0.79, 0.82] }
      },
      "pair_density": {
        "HOM_REF_HOM_REF": { "x": [0.70, 0.75, ...], "y": [0.01, 0.04, ...] },
        "HOM_INV_HOM_INV": { "x": [0.65, 0.70, ...], "y": [0.02, 0.03, ...] },
        "HOM_REF_HOM_INV": { "x": [0.65, 0.70, ...], "y": [0.05, 0.18, ...] }
      }
    },
    "LG28_cand_02": { ... }
  }
}
```

**Density**: pre-computed on HPC as a 50–100 point KDE (so the atlas
doesn't redo the KDE in JS — it just plots line segments). Pair_summaries
quantiles are the 5-number summary (min, Q1, median, Q3, max) for an
optional violin/box overlay.

## 3. Atlas-side state slot

```
state.cheat30Results = {
  LG28: {
    schema_version: "cheat30_v1",
    generated_at: "2026-05-15T03:14:22Z",
    candidates: {
      "LG28_cand_01": { ... },
      ...
    },
  },
  ...
};
```

Same shape pattern as `state.repeatDensity` and `state.ncRNADensity`. Per-chrom
keying so a multi-chrom atlas session can carry results for all chromosomes
loaded.

Per-candidate access:
```
function _cheat30ForCandidate(c) {
  if (!c || !state.cheat30Results) return null;
  const block = state.cheat30Results[c.chrom];
  if (!block || !block.candidates) return null;
  const cid = (c.id != null) ? c.id : c.candidate_id;
  return block.candidates[cid] || null;
}
```

## 4. Panel layout (page 3 / DOM page2)

Inserted into the candidate-focus stack between
`candidateRegimeRowHtml(c)` and `candidateNotesHtml(c)`. New builder:
`candidateAgeOriginHtml(c)`.

Layout (vertical stack, ~300px tall when populated, ~60px empty state):

```
┌─────────────────────────────────────────────────────────────────────┐
│  Age & origin   (cheat30 / GDS by genotype)            [Q5 ⓘ chip]  │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  ┌─────── ridgeline plot ────────┐  ┌──── verdict & numbers ────┐   │
│  │  HOM_REF / HOM_REF  ▂▄█▆▂      │  │ origin: RECURRENT          │   │
│  │  HOM_INV / HOM_INV  ▁▂█▄█▃     │  │  ⚠ I/I distribution bimodal │   │
│  │  HOM_REF / HOM_INV    ▄█▆▁     │  │     (dip P = 0.012)         │   │
│  │                                │  │                             │   │
│  │  GDS axis 0.6 ... 0.95         │  │ separation P:  4.8 × 10⁻¹²  │   │
│  │  (vertical line at mean_diff)  │  │ separation effect: 0.082    │   │
│  │                                │  │ age proxy: 0.067            │   │
│  │                                │  │   (D/D vs D/I gap; larger   │   │
│  │                                │  │    = older arrangements)    │   │
│  │                                │  │                             │   │
│  │                                │  │ n: REF=60 HET=106 INV=60    │   │
│  └────────────────────────────────┘  └─────────────────────────────┘ │
│                                                                     │
│  Caveat: cheat30 measures genotype-pair GDS — this is "shared       │
│  haplotype background per arrangement", NOT a clock. Age proxy is   │
│  ordinal (compare candidates), not absolute (no Mya).               │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

Empty state (no `cheat30_gds_results.json` loaded for this chrom):

```
┌─────────────────────────────────────────────────────────────────────┐
│  Age & origin   (cheat30 / GDS by genotype)             [layer ?]   │
├─────────────────────────────────────────────────────────────────────┤
│  No cheat30 results loaded for LG28.                                │
│  Run STEP_cheat30_gds_by_genotype.R on this chromosome and load     │
│  cheat30_gds_results_LG28.json via the schema badge drop-zone.      │
└─────────────────────────────────────────────────────────────────────┘
```

Empty state when results are loaded but this candidate isn't in them:

```
│  No cheat30 result for this candidate.                              │
│  cheat30 ran on LG28 but did not return results for LG28_cand_01.   │
│  Likely cause: too few samples per class (need ≥5 in REF and INV).  │
│  Counts: REF=3 HET=4 INV=2  (below MIN_SAMPLES_PER_CLASS)           │
```

## 5. Ridgeline plot details

SVG, 300×180px, three ridges stacked vertically:

- **HOM_REF / HOM_REF** ridge in the cohort REF colour (canonical green
  from karyotype palette)
- **HOM_INV / HOM_INV** ridge in the cohort INV colour (canonical red)
- **HOM_REF / HOM_INV** ridge in a neutral grey

Each ridge plotted from `pair_density.{type}.{x, y}` arrays. X-axis is
GDS (0.5 ... 1.0 typical range, autoscale). A dashed vertical line at
each ridge's mean. Hover any ridge → tooltip with n, mean, sd, quartiles.

If `is_bimodal: true` for the I/I ridge, show a small `⚠` icon next to
the ridge label and color the ridge with a striped pattern fill (visual
cue for "the recurrent inversion signal").

## 6. Verdict-and-numbers panel

A simple grid of (label, value) rows with the origin-class pill at the
top. Class colour palette:

| origin_class       | colour   | label text                          |
|--------------------|----------|-------------------------------------|
| `single_origin`    | green    | SINGLE ORIGIN (one founding event)  |
| `recurrent`        | red      | RECURRENT (multiple origins)        |
| `weak_signal`      | amber    | WEAK SIGNAL (genotype-GDS link present but soft)|
| `inconclusive`     | grey     | INCONCLUSIVE (insufficient power)   |

Below the pill: a small explanation paragraph derived from the inputs:

- `single_origin` →  "I/I distribution is unimodal (dip P = X). All
  HOM_INV carriers share haplotype background — consistent with one
  founding inversion event."
- `recurrent` →  "I/I distribution is bimodal (dip P = X). HOM_INV
  carriers cluster on multiple haplotype backgrounds — consistent with
  ≥2 independent inversion events at this locus."
- `weak_signal` → "Same-genotype pairs are more similar than
  different-genotype pairs (Wilcoxon P = X) but the effect is small
  (Δ = Y). Consistent with either young polymorphism, segregating
  introgression, or technical artifact."
- `inconclusive` → "Could not classify: separation test n.s. or class
  size below threshold."

Numeric block below: separation_p, separation_effect, age_proxy, n_ref,
n_het, n_inv, mean_ibs_same, mean_ibs_diff. Each with a `?` tooltip
explaining what the number is.

## 7. Manuscript-bundle integration

New helper `_bundleAgeOriginBlock(c)` joins the existing 6 enrichment
helpers from turn 118. Returns null when no cheat30 result for this
candidate. Renders:

```markdown
**Age & origin** (cheat30 GDS by genotype):

- **origin_class**: `RECURRENT`
- **separation P**: 4.8e-12 (Wilcoxon, same-genotype > different-genotype GDS)
- **separation effect**: 0.082 (mean GDS gap)
- **age proxy**: 0.067 (mean(REF/REF) – mean(REF/INV); larger = older)
- **I/I bimodality** (dip test): dip = 0.041, P = 0.012, **bimodal**
- **n samples**: REF=60, HET=106, INV=60
- **mean GDS same-genotype**: 0.847; **different-genotype**: 0.765

*Interpretation*: HOM_INV samples carry multiple haplotype backgrounds
(I/I distribution is bimodal). Consistent with ≥2 independent inversion
events at this locus across the hatchery cohort, NOT a single ancient
polymorphism.
```

TSV scalar columns added to `_candidateTSV`:

- `cheat30_origin_class`
- `cheat30_separation_p`
- `cheat30_separation_effect`
- `cheat30_age_proxy`
- `cheat30_dip_p`
- `cheat30_is_bimodal`  (boolean → 0/1)
- `cheat30_n_ref`, `cheat30_n_het`, `cheat30_n_inv`
- `cheat30_mean_gds_same`, `cheat30_mean_gds_diff`

All nullable; empty string when cheat30 not loaded.

## 8. Q5 chip (deferred to a follow-up spec)

The "Q5 ⓘ chip" in the panel header is a *placeholder* for the
question-completion chip system. The chip itself is a 24×24px square
button labelled "Q5" that hovers up a 39-row tooltip listing each Q5
key (resolved ✓, missing ✗, aspirational ◯) with overall % progress.

That's a cross-cutting feature (every per-question panel could have
one) and needs a registry-key TSV loader. **Spec'd separately in
`SPEC_question_completion_chips.md`** — not implemented in this
turn.

## 9. Card-flip interaction (deferred)

User mentioned "maybe the candidate page is too crowded — flip the
card on click to see the details". Defer until the panel ships and we
see whether crowding is actually a problem. CSS card-flip is a known
pattern (perspective + transform) but adds complexity for a problem we
might not have. Empirically check first.

If we do add it later: the front of the card stays the existing
candidate-focus stack; the back shows ALL the per-question detail
(tier, age/origin, completion chips for Q1–Q7) in a denser grid layout.
Click anywhere on the card (or a corner ⤺ flip button) to toggle.

## 10. Implementation checklist for next turn

- [ ] `_isCheat30JSON(parsed)` — schema validator (~20 lines)
- [ ] `_storeCheat30Results(parsed)` — chrom-keyed store (~30 lines)
- [ ] Drop-zone integration — extend the existing schema-badge file
  router so dropping a `cheat30_gds_results_*.json` lands in the right
  store (~10 lines)
- [ ] `_cheat30ForCandidate(c)` — accessor (~10 lines)
- [ ] CSS for `.ageorig-*` classes (~40 lines)
- [ ] `candidateAgeOriginHtml(c)` — main panel HTML builder (~80 lines)
- [ ] `_drawCheat30Ridgeline(svg, density, opts)` — SVG renderer
  (~80 lines: 3 ridges, mean lines, axis, legend)
- [ ] Wire into `renderCandidateMetadata()` between
  `candidateRegimeRowHtml(c)` and `candidateNotesHtml(c)` (~2 lines)
- [ ] `_bundleAgeOriginBlock(c)` — manuscript bundle markdown
  (~50 lines)
- [ ] TSV scalar columns in `_bundleEnrichmentScalarsForCandidate` +
  `_candidateTSV` headers (~20 lines)
- [ ] Caveats addition (1 line)
- [ ] Source-summary status row (1 line)
- [ ] Test: `test_age_origin_panel.js` with synthetic cheat30 results,
  verify panel mounts, verdict pill colour matches origin_class,
  ridgeline SVG renders without throwing, bundle block emits expected
  markdown (~150 lines, ~20 assertions)

**Estimated atlas line count**: +400 lines (single focused turn).
