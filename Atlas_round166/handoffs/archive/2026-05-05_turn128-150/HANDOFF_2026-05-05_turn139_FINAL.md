# HANDOFF — turn 139 — H-label classifier Slice 1 (no labels yet, calibration mode)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (68,001 lines, 1380/0 tests passing)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC (account `lt200308`).
**Supersedes**: `HANDOFF_2026-05-05_turn136_FINAL.md` (G-panel karyotype Slice 2).

This turn ships **Slice 1** of `SPEC_observable_allele_h_label_system.md`
(the spec promoted from `specs_new_turn131/` to `specs_todo/` at the
start of this turn). Slice 1 ships **classification + regime inference
+ calibration TSV export — NO H-LABELS YET**, by deliberate design (the
spec's empirical-first ordering: observe before labelling).

---

## 0. Cohort discipline (NEVER conflate)

Same as before. Three separate cohorts:
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent
surname.

---

## 1. What this turn shipped

### 1.1 Per-band classifier (`_classifyHLabelBands`)

For a focused candidate, classifies each K-band as one of:

| Class | Trigger |
|---|---|
| **HOM** | dosage-HET fraction < 0.30 across band members, n_dosage_non_na ≥ 5 |
| **HET** | dosage-HET fraction ≥ 0.70 across band members, n_dosage_non_na ≥ 5 |
| **AMBIGUOUS** | 0.30 ≤ fraction < 0.70, n ≥ 5 |
| **NO_DOSAGE** | n_dosage_non_na < 5 OR no dosage data at all |

Reads dosage per-sample state via `_computeHetRateForL2(ref_l2)` (the
existing turn-128d / 129 pipeline), thresholding rates ∈ [0.30, 0.70]
as HET, others as HOM, NaN as NA.

Caches result on `candidate.h_classification` (matches SPEC §5 JSON
shape minus the `label` field — labels arrive in Slice 2).

### 1.2 Regime inference

From `(n_hom, n_het)` band counts, derives:

- `implied_H` from het count via SPEC §3.1 lookup (0:1, 1:2, 2:3, 3:3, ≥4:4)
- `n_inferred_systems = max(0, implied_H - 1)`
- `consistency` ∈ {clean, partial_obs_minor, partial_obs_major, suspect, suspect_extra_hom}
  per SPEC §3.2

Per-sample verdicts (AGREE / MISMATCH_FALSE_HOM / MISMATCH_FALSE_HET /
NO_CALL) only generated for samples in HOM/HET-classified bands.
Samples in AMBIGUOUS or NO_DOSAGE bands skip — verdict requires both
axes resolved.

### 1.3 G-panel karyotype tab additions

**Regime summary line** (above per-band rows): shows count breakdown,
implied H, system count, consistency outcome (color-coded — green for
`clean`, orange for partial/ambiguous, red for `suspect*`), and overall
dosage-agreement fraction.

**Per-band classification chip** (new column in band rows, between
count and CGA preview): shows `HOM`/`HET`/`AMBIGUOUS`/`NO_DOSAGE` with
color coding and tooltip explaining the dosage-HET fraction.

**Calibration export button** in tab footer: `📊 calibration TSV
(regime histogram)`. Iterates `state.candidateList`, classifies any
candidate not yet classified, emits a TSV with one row per candidate:

```
candidate_id  K  n_hom  n_het  n_ambiguous  n_no_dosage  implied_H  n_systems  consistency  agree_fraction  n_called
```

Quentin runs this on real LG28+LG12+a few more, eyeballs the
`(n_hom, n_het)` histogram, and reports back. **Slice 2 ships H-label
assignment using the empirical regime distribution as the constraint
set.**

### 1.4 Cache invalidation hook

`_invalidateHetRateCache()` (the existing dosage-rotation hook) now
also calls `_invalidateHLabelClassificationCache()`, dropping the
classifier cache on every candidate when the dosage layer reloads.
Prevents stale classifications when the user reloads dosage data.

---

## 2. What this turn did NOT ship

- **H-labels** (`H1/H1`, `H1/H2`, etc.) — explicitly deferred to Slice 2
  per the spec's empirical-first ordering. Bands keep their existing
  ordinal/legacy/detailed labels for now (via the unchanged
  `getKaryotypeLabel(bandIdx, K, mode)`).
- **PCA mismatch highlighting** — SPEC §6 Slice 3, deferred.
- **Inheritance tab H-label propagation** — Slice 4, deferred.
- **Missing-class identification** — `missing_hom_classes` /
  `missing_het_classes` arrays in the JSON output are present but
  empty. Populated in Slice 2 once we know what regimes are real.
- **R-side STEP_HL01 producer** — separate spec
  (`SPEC_STEP_HL01_h_label_classifier.md`), not yet drafted. Atlas
  consumes the dosage_chunks layer directly via `_computeHetRateForL2`,
  no R-side gating.

---

## 3. Files changed

```
Inversion_atlas.html
  ├─ +block at ~line 37561 (after _gpKaryoColor window-export):
  │    "H-LABEL CLASSIFIER — Slice 1 (turn 139)" — adds:
  │      - 5 constants (_HLABEL_HET_FRACTION_HIGH/LOW, MIN_DOSAGE_N,
  │        SAMPLE_HET_RATE_LO/HI)
  │      - _hlabelDosageStateFromRate(rate)
  │      - _hlabelImpliedHFromHetCount(nHet)
  │      - _hlabelRegimeConsistency(nHom, impliedH)
  │      - _classifyHLabelBands(candidate, opts?)
  │      - _invalidateHLabelClassificationCache()
  │      - _hlabelBuildRegimeHistogramTSV(opts?)
  │      - _hlabelExportRegimeHistogramTSV()
  │      - 10 window.X exports
  │
  ├─ ~modified _gPanelRenderTabKaryotype (added ~70 LOC):
  │    - Regime summary line above band rows (count breakdown,
  │      implied H, n_systems, consistency, agree fraction)
  │    - Empty-state hint when classifier returns null
  │    - Per-band classification chip column (grid changed from
  │      5 to 6 columns: 18px 80px 60px 90px 1fr auto)
  │    - Calibration export button #gpHlabelCalibrationBtn in footer
  │
  ├─ +calibration button wiring (~5 LOC) in karyotype-tab interaction
  │    branch of _renderGPanelModal
  │
  └─ +cache invalidation hook (~3 LOC) in _invalidateHetRateCache

specs_new_turn131/SPEC_observable_allele_h_label_system.md
  → moved to specs_todo/ (the spec promoted at the start of this turn)

tests/test_turn139_h_label_classifier_slice1.js  [NEW, 143/0]
  Sections:
    1. Source-level: 7 functions + 5 constants + 4 window exports
    2. Helper extraction utility
    3. _hlabelDosageStateFromRate: NaN/null/Infinity/boundary cases
       (rate=0.30 → HET inclusive, rate=0.70 → HET inclusive,
        rate=0.85 → HOM)
    4. _hlabelImpliedHFromHetCount: 0..6 + defensive (-1, null, NaN)
    5. _hlabelRegimeConsistency: clean / partial_obs_minor /
       partial_obs_major / suspect / suspect_extra_hom + defensive
    6. _classifyHLabelBands sandbox (10 sub-scenarios):
       a. K=3 clean (HOM/HET/HOM, all bands ≥6 samples)
       b. K=4 partial_obs_minor (HOM/HET/HOM/HET → 2 hom + 2 het → H=3)
       c. K=4 suspect_extra_hom (3 hom + 1 het)
       d. AMBIGUOUS band (het_fraction = 0.5)
       e. NO_DOSAGE (n<min_n=5)
       f. NO_DOSAGE (all NaN dosage rates)
       g. Per-sample verdicts (MISMATCH_FALSE_HOM / MISMATCH_FALSE_HET)
       i. Caching (cache hit, force=true bypass)
       j. Null/missing inputs → null
       k. No dosage layer (graceful all-NO_DOSAGE fallback)
    7. _hlabelBuildRegimeHistogramTSV: header + per-candidate rows,
       runMissing flag classifies on the fly
    8. _hlabelExportRegimeHistogramTSV: sandbox returns false, real
       Blob mock returns true with correct filename pattern
    9. _invalidateHLabelClassificationCache: focused + list candidates
    10. Karyotype tab additions (source-level): chip column, regime
        summary, calibration button, all 4 chip styles
    11. Cache invalidation hook plumbing
```

LOC delta: 67,474 → 68,001 = **+527**. Of that, ~440 is the new
classifier block, ~70 is the karyotype tab modifications, ~10 is
the wiring + cache hook, ~5 is window exports.

---

## 4. Tests

```
tests/test_turn139_h_label_classifier_slice1.js — 143 / 0
```

Coverage map highlights:

| Sub-scenario | What it locks down |
|---|---|
| 3 — dosage state | rate=0.30 → HET (boundary inclusive); rate=0.85 → HOM (above HI) |
| 4 — implied_H | n_het=2 → H=3 (3 alleles, one missing cross); n_het=3 → H=3; n_het≥4 → H=4 capped |
| 5 — consistency | (n_hom=2, H=3) → partial_obs_minor; (n_hom=3, H=2) → suspect_extra_hom |
| 6a — K=3 clean | HOM/HET/HOM correctly classified, agreement = 100% |
| 6b — K=4 partial | 2 hom + 2 het → implies H=3, partial_obs_minor |
| 6c — K=4 extra hom | 3 hom + 1 het → implies H=2, suspect_extra_hom |
| 6d — AMBIGUOUS | 50% het_fraction → AMBIGUOUS classification |
| 6e — NO_DOSAGE n<min | n=4 with all dosage HET → NO_DOSAGE despite het_fraction=1 |
| 6f — NO_DOSAGE NaN | All-NaN rates → NO_DOSAGE, het_fraction=NaN |
| 6g — verdicts | sample 5 (HOM band, dosage HET) → MISMATCH_FALSE_HOM |
| 6i — cache | First call caches; mutating cache reflects on second call; force:true bypasses |
| 7 — regime TSV | Header line + 1 row per candidate; runMissing classifies on the fly |
| 8 — TSV export | Sandbox returns false; real Blob mock returns true with `h_label_regime_histogram.<chrom>.tsv` filename |

| | Tests | Files |
|---|---|---|
| Turn 132 baseline | 879 | 33 |
| Turn 135 G-panel S1 | 1136 | 36 |
| Turn 136 G-panel S2 | 1237 | 37 |
| Turn 139 H-label S1 | 1380 | 38 |
| Δ this turn | +143 | +1 |

Zero regressions. Turn 136's 100/0 still passes despite the
karyotype-tab grid changing from 5 to 6 columns (the test was looking
for source-level patterns that survived the change).

---

## 5. Sandbox vs. real-data note

What the sandbox does NOT exercise:

- **Real `_computeHetRateForL2` integration**. Tests stub it via
  `rateMap`. Browser smoke test: pin a candidate on real data, open
  G-panel karyotype tab, verify regime summary appears with correct
  numbers (load dosage data first via the heatmap/stripe-quality
  toolbar to populate the dosage chunk LRU cache).
- **Real Blob download for the calibration TSV**. Tests verify the
  call path returns true with a Blob mock; actual byte content of
  the TSV isn't inspected.
- **Real cache invalidation chain**. Tests verify the source-level
  hook in `_invalidateHetRateCache` calls
  `_invalidateHLabelClassificationCache`. Browser smoke test: load
  dosage data, classify a candidate (open karyotype tab), reload
  dosage data via the heatmap reset action, verify the regime
  summary recomputes (rather than showing stale numbers).

Browser smoke test sequence (~2 minutes):

1. Pin a candidate on page 1 (lock + promote on a real LG28 region)
2. Load dosage data via the page-1 heatmap toolbar or stripe-quality
   button
3. Press `g` to open G-panel
4. Click `karyotype` tab → see:
   - Regime summary line: `counts: X hom · Y het · Z ambiguous /
     no-dosage` + `implied: H=N · M systems` + `consistency: clean`
     (or whatever)
   - Per-band rows with classification chips (HOM in green, HET in
     gold/yellow, AMBIGUOUS in orange, NO_DOSAGE in dim grey)
   - Footer with the new `📊 calibration TSV (regime histogram)`
     button
5. Click the calibration button → file
   `h_label_regime_histogram.<chrom>.tsv` downloads with one row per
   saved candidate
6. Open the TSV in a viewer or `head` it on the terminal — the
   `consistency` column histogram tells Quentin which regimes appear

If the regime summary shows all NO_DOSAGE: dosage data isn't loaded,
or the candidate's `ref_l2` doesn't have a chunk in the LRU cache.
Load via the heatmap to populate.

---

## 6. State / window slots claimed this turn

```
candidate.h_classification    object — cache for classifier output
                               (cleared by _invalidateHLabelClassificationCache)
```

```
window-exposed (this turn):
_classifyHLabelBands
_invalidateHLabelClassificationCache
_hlabelBuildRegimeHistogramTSV
_hlabelExportRegimeHistogramTSV
_hlabelDosageStateFromRate
_hlabelImpliedHFromHetCount
_hlabelRegimeConsistency
_HLABEL_HET_FRACTION_HIGH       // 0.70
_HLABEL_HET_FRACTION_LOW        // 0.30
_HLABEL_MIN_DOSAGE_N            // 5
```

```
DOM ids claimed (this turn):
#gpHlabelCalibrationBtn       — calibration TSV export button
```

```
localStorage keys claimed: (none — Slice 1 has no persistent state.
Classifier results live on candidate objects which already persist
via candidateList localStorage path)
```

---

## 7. Backups present

```
Inversion_atlas.html.bak_pre_l2_sweep_slice1     (turn 133 baseline)
Inversion_atlas.html.bak_post_l2_sweep_slice1    (turn 133 final)
Inversion_atlas.html.bak_post_l2_sweep_inspector (turn 134 final)
Inversion_atlas.html.bak_pre_g_panel_slice1      (turn 135 baseline)
Inversion_atlas.html.bak_post_g_panel_slice1     (turn 135 final)
Inversion_atlas.html.bak_pre_g_panel_slice2      (turn 136 baseline)
Inversion_atlas.html.bak_post_g_panel_slice2     (turn 136 final)
Inversion_atlas.html.bak_pre_h_label_slice1      (this turn baseline)
Inversion_atlas.html.bak_post_h_label_slice1     (this turn final)
```

`.bak_*` files NOT in bundle.

---

## 8. Where to start the next chat

The natural next step is **Slice 2** of the H-label spec, but it's
**gated on Quentin running the calibration TSV on real data**. Until
that's done, Slice 2 doesn't have the empirical regime list it needs.

### Option 8a — Wait on calibration, ship something else first

Recommended. Quentin runs Slice 1 on LG28+LG12+a few more next time
he's at LANTA, exports the regime histogram TSV, and reports back.
Meanwhile, parallel options:

- **G-panel Slice 3 (inheritance tab content)** — the third tab of
  the unified groups popup. ~1 turn. Slice 1 placeholder still in
  place. Doesn't block on H-label calibration.
- **L2-sweep calibration on real LANTA data** — also waiting on LANTA
  runs. Same deferral.
- **Trajectory matrix viewer Slice 3** — separate priority axis,
  standalone.

### Option 8b — H-label Slice 3 (PCA mismatch highlighting)

Doesn't depend on Slice 2. ~0.5 turn. Adds:
- `state.hlMismatchHighlight: bool` toggle
- A new G-panel button `⚠ flag mismatches` that sets it
- Repaint hook in `drawPCA` that paints MISMATCH samples with a 2px
  red ring (no fill change)
- Tests

### Option 8c — H-label Slice 2 directly (assume biological priors)

If Quentin doesn't want to wait on the empirical histogram, Slice 2
can ship with a **biological-prior regime list** assumed: just `[(2 hom +
1 het), (3 hom + 3 het)]` (the two clean cases) plus `partial_obs_minor`
fallback. Less honest than the empirical-first design but unblocks
labelling now. Recommend NOT doing this unless Quentin explicitly
asks — the calibration histogram is the whole point.

### Recommendation

**8a**, with 8b as a small bonus while waiting on calibration.
Inheritance tab (G-panel Slice 3) is probably more useful than the
PCA mismatch highlighting if you have to pick one.

---

## 9. Honest framing

**What turn 139 actually delivered:**

- A pure-JS classifier that derives band-level regime classifications
  + per-sample dosage agreement verdicts from the existing
  dosage_chunks layer, without H-labels.
- Visible surfacing of the regime inference in the G-panel karyotype
  tab — the user sees `2 hom · 1 het · implies H=2 · 1 system · 94.3%
  agreement` immediately on opening the tab.
- A calibration TSV export so Quentin can build the empirical regime
  distribution that Slice 2 needs.
- 143 tests covering all classification axes, regime derivation,
  per-sample verdicts, caching, TSV construction, and the cache
  invalidation hook.

**What it deliberately doesn't deliver:**

- H-labels (`H1/H1`, `H1/H2`, etc.). The placeholder ordinal labels
  (`band 1 (lo)` in legacy mode, `H1/H1` in detailed mode but per
  the existing fixed table — NOT data-driven) are unchanged. Slice 2
  ships data-driven H-label assignment after the calibration histogram
  is in.
- Missing-class identification. The `implied_regime.missing_hom_classes`
  array is present in the JSON but empty until Slice 2.
- Any modification to the `_KARYO_DETAILED_LABELS` legacy fallback
  path. That stays as-is for sessions without dosage data.

**What this means for the manuscript:**

Methods section can now claim: *"Per-candidate band classifications
were computed by cross-referencing K-means cluster membership against
per-sample dosage HET/HOM calls; bands with dosage-HET fraction ≥
0.70 were classified as heterozygous classes, bands < 0.30 as
homozygous classes, and intermediate or low-coverage bands as
low-confidence."*

Once the calibration TSV is in hand, the Methods section gains a
quantitative claim about modal regimes and multi-system locus
prevalence (per spec §10).

---

## 10. Bundle contents

```
Inversion_atlas.html              (current, 68,001 lines)
tests/                            (all *.js)
  ├─ test_turn139_h_label_classifier_slice1.js  [NEW, 143/0]
  ├─ test_turn136_g_panel_karyotype_slice2.js   (unchanged, 100/0)
  ├─ test_turn135_g_panel_slice1.js             (unchanged, 73/0)
  ├─ test_turn134_l2_sweep_inspector.js         (unchanged, 104/0)
  ├─ test_turn133_l2_sweep_auto_promote.js      (unchanged, 81/0)
  └─ test_turn132_*.js                          (unchanged, 321/0)
specs_todo/                       (active build queue)
  └─ SPEC_observable_allele_h_label_system.md   (promoted from review queue this turn)
specs_new_turn131/                (pending review queue, minus the H-label spec)
HANDOFF_2026-05-05_turn139_FINAL.md  (this file)
all previous handoffs             (kept for history)
OBSERVATIONS_TO_FIX.txt
```

`.bak_*` files NOT in bundle.

---

Walk the map carefully, respect cohort discipline, don't break the
test suite. Slice 1 is in a clean stopping state — Quentin's next
move is the calibration TSV run on real LANTA data. While that's
pending, Slice 3 (PCA mismatch highlighting) or G-panel Slice 3
(inheritance tab) are the recommended parallel work.
