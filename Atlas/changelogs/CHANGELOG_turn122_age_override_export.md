# CHANGELOG — turn 122: age model + override + export + TE fragility (page16b)

**Atlas Δ:** +666 lines (60,139 → 60,805)
**Tests:** 79 new (407/407 cumulative green) · No regressions on turns 117/118/119/120/121.

---

## What shipped — all four deferred items from turn 121

### 1. Age-model auto-suggest chip (5 classes)

Computes age model from the comparative-evidence signals (lineage_distribution + dXY fold-elevation), with rules ordered by specificity:

- **MULTI-AGE-HOTSPOT** — mixed rearrangement types across ≥3 species. Strongest mechanism class (Quentin's words: *"the strongest Nature-style mechanism class"*).
- **LINEAGE-KARYO** — fission/fusion event_type AND ≥1 sister species shows fission/fusion. Karyotype event polarized to a lineage.
- **OLD-BP-YOUNG-INV** — boundary shared in ≥3 species AND dXY fold not strongly elevated (<1.5×). The "fragile region is old, current inversion is young" case.
- **OLD-POLY** — dXY fold ≥1.5× AND ≥2 species share the boundary. Old polymorphism; standard and inverted haplotypes have been diverging long.
- **YOUNG-POP** — dXY fold <1.2× AND boundary present in ≤1 species. Recent within-Cgar polymorphism.

Returns `{age_model, confidence, rationale, signals}`. Rationale is the manuscript-ready prose — same pattern as architecture chip.

### 2. Manual classification override

Editor at the top of the center column on page 16b:

| Field | Type | Behaviour |
|---|---|---|
| Architecture | dropdown (— auto — / A / B / C / D / E / F) | overrides architecture chip |
| Age model | dropdown (— auto — / 5 age tags) | overrides age chip |
| Confidence | dropdown (— default — / low / medium / high / manual) | overrides confidence label |
| Notes | textarea | free-form manuscript-style note |
| `save override` | button | persists to `state.classifications[bp_id]` |
| `clear override` | button | removes the override (auto-suggest takes back over) |

Empty fields don't shadow auto-suggest — only the fields you set are sticky. Persisted across sessions in `localStorage` under `inversion_atlas.classifications.v1`. The displayed chip says `(manual)` when the value came from the override.

### 3. Manuscript-export TSV

Bottom of the classification block — `⬇ export classification TSV — all breakpoints`. One row per `cs_breakpoint`. Columns:

```
bp_id  gar_chr  gar_pos_mb  event_type
architecture_class  architecture_source  architecture_confidence
age_model           age_source            age_confidence
n_species_total  n_species_present  n_species_fission_fusion
fold_elevation_dxy  dxy_within
classification_notes
lineage_<species_id>...      ← one column per species in the manifest
```

`*_source` is `manual | auto | ""`. Manual override always wins; falls back to auto-suggest. Filename: `multispecies_classification_<timestamp>.tsv`.

### 4. TE-fragility focal-vs-bg panel

Right column on page 16b, below the dotplot + boundary KV. Reuses `popgenFocalVsBg.makeFocalVsBgPanel({ page_id: 'page16b', anchor_type: 'cs_breakpoint', ... })` — same widget that drives page 16's TE fragility panel. Shares `state._focalVsBg.csRadiusBp` so radius preferences travel between pages.

The widget surfaces |Z| / theta_pi / FST at the breakpoint focal vs chromosome background. When a TE density layer is loaded (`repeat_density_v2`), `all_TE` joins the metric pool — that's the comparative fragility signal at the homologous region.

### 5. dxy_per_inversion_v1 layer (NEW)

Per-inversion dXY between standard and inverted haplotypes inside the inversion + flanking regions. Quantitative input to the age-model auto-suggest. Schema:

```json
{
  "tool": "dxy_per_inversion_v1",
  "schema_version": 1,
  "params": { "method": "Reich_dXY_from_dosage_matrix", "flank_size_bp": 1000000 },
  "per_inversion": [
    {
      "candidate_id": "cand_lg28_15_18",
      "chrom": "C_gar_LG28",
      "start_bp": 15115000, "end_bp": 18217000,
      "n_hom_ref": 60, "n_het": 106, "n_hom_inv": 60,
      "dxy_within_inversion_ref_vs_inv": 0.00420,
      "dxy_flank_left_ref_vs_inv":       0.00195,
      "dxy_flank_right_ref_vs_inv":      0.00210,
      "fold_elevation_inside_vs_flank":  2.10,
      "interpretation_default": "elevated_consistent_with_old_inversion"
    }
  ]
}
```

Detector + store + persist + restore + clear quintet, full dispatcher integration. localStorage key `inversion_atlas.dxyPerInversion.v1`. `_msGetDxyForBreakpoint` matches by `candidate_id` first, position-window fallback (chrom + within span).

---

## State slots added

| Slot | Purpose | Persisted |
|---|---|---|
| `state.dxyPerInversion` | parsed `dxy_per_inversion_v1.json` | localStorage `inversion_atlas.dxyPerInversion.v1` |
| `state.classifications` | per-breakpoint manual overrides | localStorage `inversion_atlas.classifications.v1` |

---

## Demo files

- `dxy_per_inversion_v1.demo.json` — 2 inversions: `cand_lg28_15_18` (fold 2.10×, → OLD-POLY when paired with the synteny demo) and `cand_lg12_8_11` (fold 1.07×, → YOUNG-POP).

---

## Test infrastructure fix

The vm-sandbox `pullFunction` helper used in tests had a bug — it didn't skip `//` line comments before processing string literals. A function comment containing `// We don't need...` (with an apostrophe) would start a runaway string scan.

Fixed in both `test_turn121_multispecies.js` and `test_turn122_age_override_export.js`:
- Skip `//` comments to end of line
- Skip `/* ... */` comments to closing `*/`
- Bail out of single/double quoted strings on unescaped newlines (defensive — prevents runaway when source contains weird quote patterns)

---

## What's still NOT shipped

These would be incremental polish, not core scope:

- **busco_anchors_v1 layer** — BUSCO ticks on the dotplot, especially valuable for the deep-divergence (Trichomycterus, Silurus) backbone.
- **comparative_te_breakpoint_fragility_v1 layer** — pre-computed per-species TE density at homologous breakpoints (atlas would consume this directly into the focal-vs-bg metric pool).
- **Tree internal-node click → MRCA comparison** — currently leaves only.
- **Hypothesis registry page integration** — the framework's hooks (`predictions[].layer === "dxy_per_inversion"`) for marking hypotheses as supported/refuted aren't wired yet.

---

## Files

```
Inversion_atlas.html               +666 lines (60,139 → 60,805)
test_turn122_age_override_export.js  NEW (79 tests)
test_turn121_multispecies.js         FIXED (pullFunction comment-skip)
dxy_per_inversion_v1.demo.json     NEW
CHANGELOG_turn122_age_override_export.md  THIS FILE
```

---

## Operating notes

The age-model auto-suggest is **deliberately conservative** when evidence is missing. Without a `dxy_per_inversion_v1` layer loaded, it can still reach LINEAGE-KARYO (which doesn't need dXY) and MULTI-AGE-HOTSPOT (also doesn't, because mixed rearrangement types across many species is itself the signal). All other classes require dXY. So the page populates progressively as you load layers — there's no "garbage in, garbage out" failure mode where it confidently emits YOUNG-POP without seeing any dXY.

The override editor explicitly does NOT save empty fields, so a user can override JUST architecture without touching age model — and the auto-suggest stays in charge of the other label. This is what you want for the manuscript pipeline: the auto-suggest takes care of 80% of cases, and you spend manual time only on the contested calls.

The TSV exports both classification + the auto-suggest's evidence-count signals (`n_species_total`, `n_species_present`, `n_species_fission_fusion`, `fold_elevation_dxy`, `dxy_within`). That's what reviewers will need to audit a class call. Plus the per-species lineage status columns at the end for full reproducibility.
