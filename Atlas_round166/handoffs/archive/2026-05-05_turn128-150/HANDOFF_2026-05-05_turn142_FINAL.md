# HANDOFF — turn 142 — Breeding-readiness card, Turn A (data plumbing)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (68,679 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC.
**Supersedes**: turn 141 candidate-bands handoff.

This is **Turn A** of a 4-turn build of
`SPEC_per_candidate_breeding_readiness_card.md` (Atlas 5 Part B from
chat `c03fc41e`, "the deliverable that converts the paper from a
population genomics study into a hatchery management resource"). Turn A
ships the foundation: the atlas-side data layer for cohort F_ROH / H /
K=8 ancestry, joined to the per-chrom sample roster by CGA id.

> Quentin (chat `c03fc41e`):
> *"Atlas 5 — Per-individual integration atlas / Per-candidate
> breeding-readiness atlas (Supplementary). [...] Part B —
> Per-candidate breeding-readiness card (one page per marker-ready
> candidate)."*

---

## 0. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent surname.

---

## 1. Why this is a multi-turn build

The breeding-readiness card spec has many fields with different
prerequisite levels:

| Field | Prereq | Available? |
|---|---|---|
| chrom / coordinates / span / tier | `state.candidate` | ✓ |
| PCA panel | existing `drawPCA` | ✓ |
| Karyotype counts (REF/HET/INV) | turn 139 h_classification + PC1 sort | ✓ |
| **Carrier × K=8 ancestry table** | F_ROH layer + K=8 from Diversity atlas | ✓ once Turn A lands |
| **Per-arrangement F_ROH + Wilcoxon** | F_ROH per-sample × candidate karyotype | ✓ once Turn A lands |
| Per-arrangement damaging-load | MODULE_CONSERVATION (deleterious calls) | ✗ (R-side not shipped) |
| Recombinant fraction | dosage changepoint detector | ✗ (separate spec) |
| Marker panel | marker panel TSV | ✗ (separate spec) |
| Cross-species orthology | `state.crossSpecies` | ✓ |
| Atlas reviewer link | `cand.id` | ✓ |
| Auto-pairing advice | computed from above | partial (F_ROH-only viable now) |

The four-turn split:

- **Turn A (THIS)** — Data plumbing: cohort_diversity_v1 loader, byCGA
  index, sample-idx → diversity row resolver, fixture file, file-picker
  dispatch, localStorage round-trip.
- **Turn B (next)** — Computation: per-arrangement F_ROH stats,
  Wilcoxon rank-sum (using existing `normalCDF` at line ~27315),
  carrier × K=8 cross-tab, candidate-karyotype derivation from
  h_classification + PC1 ordering, breeding-card data assembly,
  pairing-advice rule engine.
- **Turn C** — Render: print-ready HTML one-pager (matches SPEC §1
  visual contract), inline SVG PCA reuse, A4 print stylesheet, "Print"
  button on candidate focus (page 3).
- **Turn D** — Catalogue integration: "Generate breeding cards" button
  on page 5, per-tier filter, per-card download as standalone HTML,
  bulk export.

---

## 2. What turn A ships

### 2.1 JSON loader (~290 LOC, after `_clearCrossSpecies` ~line 20034)

```
const COHORT_DIVERSITY_TOOL    = 'cohort_diversity_v1';
const COHORT_DIVERSITY_LS_KEY  = 'inversion_atlas.cohort_diversity';
const _COHORT_DIVERSITY_NUMERIC_COLS = [
  'h', 'f_hom', 'f_roh',
  'roh_total_bp', 'roh_n', 'roh_longest_bp', 'roh_mean_bp',
  'th_in', 'th_out', 'th_ratio', 'callable_bp',
];
function _isCohortDiversityJSON(data)         { ... }   // dual-shape detector
function _normalizeCohortDiversityRow(row)    { ... }   // NaN/wrong-type → null
function _storeCohortDiversity(parsed)        { ... }   // build state.cohortDiversity + byCGA Map
function _persistCohortDiversity()            { ... }   // localStorage write (drops byCGA)
function _restoreCohortDiversity()            { ... }   // round-trip rebuild
function _clearCohortDiversity()              { ... }
function _diversityForSampleIdx(si)           { ... }   // chrom-sample → diversity row
function _diversityForCGA(cga)                { ... }   // direct lookup, case-insensitive
function _cohortDiversityCoverageOnCurrentChrom() { ... }   // diagnostic
```

All 9 helpers exposed on `window.*` for testability.

**Two accepted JSON shapes**:

1. **Wrapped form** (preferred for new producers):
```jsonc
{
  "tool": "cohort_diversity_v1",
  "schema_version": 1,
  "generated_at": "2026-05-05T...",
  "cohort": { "n_samples": 226, "species": "Clarias gariepinus", ... },
  "samples": [
    { "sample_id": "CGA009", "k8": "K1", "pruned81": true,
      "h": 0.00469888, "f_hom": -0.0323, "f_roh": 0.25413,
      "roh_total_bp": 146609035, "roh_n": 3190, ...
      "callable_bp": 576895264 },
    ...
  ]
}
```

2. **Raw-array form** (Diversity-atlas dt_S1 paste convenience):
```jsonc
[ { "sample": "CGA009", "k8": "K1", "h": 0.0046, "f_roh": 0.254, ... }, ... ]
```

The detector auto-discriminates. Raw-array form lets the user open
`Diversity_atlas.html`, view source, find `<script id="dt_S1">`, save
the inner blob as `cohort_diversity_v1.json`, drop into the file
picker — zero pre-processing.

### 2.2 File-picker dispatch

Two integration points in `loadJSON` / `_classifyJSONKind`:

1. `_classifyJSONKind` returns `'cohort_diversity'` early (above
   `'chromosome'`), so the recent-files history reflects the right kind.
2. `loadJSON` adds a branch right after the cross-species branch
   (~line 48226) that calls `_storeCohortDiversity` + `_persistCohortDiversity`.

### 2.3 Startup restore

`try { _restoreCohortDiversity(); } catch (_) {}` added to the startup
restore block (~line 62270), right after `_restoreCrossSpecies`.

### 2.4 Real fixture: `json/cohort_diversity_v1.json`

65 KB, 226 samples derived from Diversity-atlas `dt_S1`, packed into
the canonical wrapped form. Useful as:
- A one-shot drop-in for users who want F_ROH context now
- The reference fixture for Turn B/C/D tests
- The template producer pipelines should target

```jsonc
{
  "tool": "cohort_diversity_v1",
  "schema_version": 1,
  "generated_at": "<iso>",
  "cohort": {
    "n_samples": 226,
    "species": "Clarias gariepinus",
    "reference_haplotype": "fClaHyb_Gar_LG.fa",
    "module": "MODULE_3",
    "source": "Diversity_atlas dt_S1 paste",
    "cohort_label": "226-sample pure C. gariepinus hatchery"
  },
  "samples": [ <226 rows>, ... ]
}
```

---

## 3. Tests

`tests/test_turn142_cohort_diversity_loader.js` — **133 / 0** across
14 sections:

1. Source-level definitions present (15)
2. Window exports (5)
3. File-picker dispatch wired (5)
4. Startup restore wired (1)
5. Detection behavior (14)
6. Row normalization (16)
7. Store + index (16)
8. byCGA lookup (case-insensitive, defensive) (7)
9. `_diversityForSampleIdx` (chrom-sample → row, 11)
10. Coverage diagnostic (8)
11. localStorage round-trip (10)
12. Clear (4)
13. Real fixture end-to-end (11)
14. Raw-array form (Diversity dt_S1 paste) (6)

Behavioral coverage includes:
- Both shapes (wrapped / raw-array) detected; non-detection inputs
  rejected (null, empty, chromosome JSON, cs_breakpoints JSON, missing
  f_roh, missing schema_version)
- Row normalization: dt_S1 `sample` key → canonical `sample_id`; k8 int
  → string; pruned81 truthy/falsy/garbage handling; NaN → null;
  non-numeric → null; whitespace trimming
- Indexed lookup: case-insensitive (CGA009 / cga009 / Cga009 all hit);
  null / non-string / "" args → null without throwing
- localStorage: byCGA Map dropped before persist (Maps don't serialize),
  rebuilt on restore; corrupt localStorage doesn't throw
- Defensive: missing state.data, missing state.cohortDiversity, sample
  index out of range (-1, ∞, non-int, string), no cga on a sample
- Real fixture: 226 / 0 drops, 0 duplicates, CGA009 f_roh ≈ 0.254
  matches Diversity atlas; CGA322 outlier (very low f_roh) preserved
- dt_S1 paste from Diversity atlas: detected, stored, indexed verbatim

Adjacent suites unchanged: turn 141 candidate-bands (62/0), turn 140
H-label chip (45/0), turn 139 H-label classifier (143/0), turn 132
lines renderer (38/0). Full turn-numbered suite **1620 / 0** (up from
1487 / 0 at turn 141).

---

## 4. Atlas state

| | LOC | Tests | Files |
|---|---|---|---|
| Pre-turn (turn 141 baseline) | 68,343 | 1487 | 40 |
| Post-turn (this) | 68,679 | 1620 | 41 |
| Δ | +336 | +133 | +1 |

Plus `json/cohort_diversity_v1.json` (65 KB, 226 samples) added to
the bundle.

---

## 5. The breeding-readiness build map

For Turn B and beyond, here's what's needed and what already exists:

### 5.1 Per-candidate karyotype derivation (Turn B Slice 1)

The h-label classifier (turn 139, `_classifyHLabelBands`) produces
per-band classification (`HOM` / `HET` / `AMBIGUOUS` / `NO_DOSAGE`)
and per-band `median_pc1`. To get **per-sample** REF / HET / INV:

```
karyotypeBySample[si] = (() => {
  const lab = candidate.locked_labels[si];
  const band = h_classification.bands[lab];
  if (band.classification === 'HET') return 'HET';
  if (band.classification === 'HOM') {
    // Sort all HOM bands by median_pc1 ascending; lower = REF, higher = INV
    return /* leftmost-PC1 HOM = REF, rightmost = INV */;
  }
  return null;  // AMBIGUOUS / NO_DOSAGE
})();
```

For K=3 with 2 HOM + 1 HET this is unambiguous. For K=4+ multi-haplotype
cases, defer to ordinal labels (Slice 2 of the H-label spec is supposed
to resolve this).

### 5.2 Per-arrangement F_ROH (Turn B Slice 2)

```
function _perArrangementFROH(cand) {
  const karyo = _candidateKaryotypePerSample(cand);   // §5.1
  const groups = { HOM_REF: [], HET: [], HOM_INV: [], all_carriers: [] };
  for (let si = 0; si < karyo.length; si++) {
    const div = _diversityForSampleIdx(si);
    if (!div || !Number.isFinite(div.f_roh)) continue;
    if (karyo[si] === 'HOM_REF') groups.HOM_REF.push(div.f_roh);
    else if (karyo[si] === 'HET') {
      groups.HET.push(div.f_roh);
      groups.all_carriers.push(div.f_roh);
    } else if (karyo[si] === 'HOM_INV') {
      groups.HOM_INV.push(div.f_roh);
      groups.all_carriers.push(div.f_roh);
    }
  }
  // Compute n / mean / median / sd / IQR per group, plus Wilcoxon
  // rank-sum REF vs INV (using existing normalCDF at ~line 27315)
  ...
}
```

### 5.3 Wilcoxon rank-sum (Turn B Slice 3)

Atlas already has `normalCDF` at line ~27315 (Abramowitz & Stegun
7.1.26). Implementation sketch:

```
function _wilcoxonRankSumP(a, b) {
  const n_a = a.length, n_b = b.length;
  if (n_a < 1 || n_b < 1) return { U: NaN, z: NaN, p_two_sided: NaN };
  // Pool, rank with average ties handling
  const pooled = a.map(v => ({ v, g: 'a' })).concat(b.map(v => ({ v, g: 'b' })));
  pooled.sort((x, y) => x.v - y.v);
  // ... standard rank-sum with tie correction ...
  const U = ...;
  const mu = n_a * n_b / 2;
  const sigma2 = n_a * n_b * (n_a + n_b + 1) / 12;
  // tie correction: subtract Σ(t³−t) / (12 (n+1))
  const z = (U - mu) / Math.sqrt(sigma2);
  const p = 2 * (1 - normalCDF(Math.abs(z)));
  return { U, z, p_two_sided: p };
}
```

### 5.4 Carrier × K=8 cross-tab (Turn B Slice 4)

Trivial once §5.1 + cohort data are in place: walk samples, bucket by
(`karyo[si]`, `_diversityForSampleIdx(si).k8`), build a K8 × {REF, HET,
INV} contingency table.

### 5.5 Card data assembly + pairing advice (Turn B Slice 5)

Pure data builder per SPEC §3:
```
function _buildBreedingCard(cand) {
  return {
    chrom, coordinates, tier, confidence,
    karyotype_counts: { REF, HET, INV },
    ancestry_table: <K8 × karyotype>,
    burden: { REF, HET, INV, wilcoxon_p },
    cross_species: state.crossSpecies ? <look up> : null,
    auto_advice: _generatePairingAdvice(card),
    atlas_url,
    coverage: _cohortDiversityCoverageOnCurrentChrom(),
    // Turn C (render) consumes this object.
  };
}
```

Pairing advice: SPEC §2 table, but only the F_ROH-asymmetric rule fires
in Turn B (the others need data we don't have yet).

---

## 6. What this is NOT

- **Not the card render.** Turn C ships the HTML one-pager. Turn A is
  data only.
- **Not the catalogue export.** Turn D ships bulk generation.
- **Not a producer.** This loader expects a `cohort_diversity_v1.json`
  to exist — the Diversity atlas already has the data, the loader just
  brings it across via paste-and-drop. A formal R-side producer
  (`STEP_R34_emit_cohort_diversity.R` or similar) is a future spec.
- **Not damaging-load.** MODULE_CONSERVATION outputs aren't in the
  cohort_diversity_v1 schema (intentionally — that's a separate per-
  sample/per-candidate axis). When MODULE_CONSERVATION lands, it gets
  its own `damaging_load_v1.json` layer with the same loader pattern.
- **Not the recombinant detector.** SPEC `SPEC_recombinant_dosage_changepoint_detector.md`
  is its own thing.
- **Not the marker panel.** SPEC `SPEC_marker_panel_design_atlas.md` is
  its own thing.

---

## 7. Backups

```
Inversion_atlas.html.bak_pre_breeding_readiness_turnA   (pre-Turn A baseline)
```

`.bak_*` files NOT in bundle. Re-derivable.

---

## 8. Where to start the next chat

### Option 8a — Turn B (RECOMMENDED, continues the build)

Build the computation layer per §5. Slices 1–5 fit comfortably in one
turn if scoped tight: candidate-karyotype derivation +
`_perArrangementFROH` + `_wilcoxonRankSumP` + `_carrierByK8Table` +
`_buildBreedingCard` data builder + `_generatePairingAdvice` rules.
Tests use the real fixture from Turn A.

Estimate: ~250–400 LOC + ~80 tests. The Wilcoxon implementation is the
trickiest part (tie correction, continuity correction); reuse the
existing `normalCDF`.

### Option 8b — Pivot to a different spec

The user can decide other priorities matter more. Tier 1 specs still in
queue (per `SPECS_TIER_INDEX.md`):

- **`SPEC_inversion_age_origin_atlas.md`** — Porubsky 2022, Hartigan,
  Corbett-Detig. GDS classification + dip test + age proxy. ~400 LOC.
- **`SPEC_boundary_consensus_aggregator.md`** — KDE multimodal mode
  detection. Handles unimodal / bimodal carrier-distribution.
- **`SPEC_recombinant_dosage_changepoint_detector.md`** — changepoint
  literature, per-carrier dosage-step detection.
- **`SPEC_hypothesis_test_framework_atlas.md`** — T1–T11 with BH
  correction.

### Recommendation

**8a** — finish what we started. Turn A's foundation is concretely
useful only after Turn B's computation lands; shipping just Turn A
leaves the F_ROH data resident in state but no surface that uses it.

---

## 9. Honest framing

**What turn 142 actually delivered:**
- A clean, dual-shape JSON loader that respects the existing atlas
  conventions (cs_breakpoints / dotplot / repeat_density patterns)
- A foundation that lets Turn B compute "INV-arrangement carriers have
  F_ROH 0.04 above REF mean (Wilcoxon p=0.003)" — the SPEC §2 first
  pairing-advice rule, which is the manuscript-critical claim
- A real, working fixture (`json/cohort_diversity_v1.json`, 226
  samples) so Turn B has something to test against immediately
- 133 tests covering detection, normalization, indexing, persistence,
  and end-to-end fixture round-trip

**What it deliberately didn't deliver:**
- Any computation. Turn A is pure plumbing on purpose: the right unit
  to ship + test in a single turn, leaving Turn B to focus on the
  Wilcoxon math without also debugging the loader.
- Any surface change. No new pages, no new buttons, no card render.
  The user notices Turn A's shipping only via the file picker now
  accepting a new JSON kind. That's intentional — the foundation
  should be invisible until Turn C wires it to a visible surface.
- A formal R-side producer spec. The Diversity-atlas-paste path covers
  the immediate need; a producer spec is a follow-up if Quentin wants
  the cohort_diversity layer regenerated reproducibly from MODULE_3
  outputs.

**Manuscript impact (when Turns B + C + D land):**
- Atlas 5 Part B becomes a real deliverable: per-candidate one-pager
  PDFs with PCA, karyotype counts, broodline ancestry breakdown,
  per-arrangement F_ROH + Wilcoxon, pairing advice
- The chat-`c03fc41e` framing — *"converts the paper from a population
  genomics study into a hatchery management resource"* — gets a
  concrete artifact
- Hatchery managers / aquaculture-genomics readers who don't read the
  full manuscript can scan one PDF per inversion and decide whether to
  use it for breeding decisions

Walk the map carefully, respect cohort discipline, don't break the
test suite. Turn A is a clean foundation; Turn B is the obvious next
step.
