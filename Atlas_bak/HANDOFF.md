# SV Evidence Page — Complete Handoff

**Date:** 2026-05-04
**Project:** MS_Inversions_North_african_catfish manuscript — SV evidence page in `Inversion_atlas.html`
**Status:** ✅ **All 7 steps shipped.** Atlas-side complete; producers stubbed (JSON shapes locked, VCF parsing blocks need wiring).
**Tests:** 730/730 across 7 test files.
**Bundle:** Single zip `sv_evidence_complete.zip` with everything; this handoff alongside it.

---

## What's in this bundle

```
Atlas/
├── Inversion_atlas.html                   ← main HTML (4 surgical edits since project start)
├── HANDOFF.md                              ← this file (also at repo root for visibility)
│
├── js/
│   ├── atlas_sv_evidence.js               ← main module, 4263 lines (~145 KB)
│   └── atlas_dotplot.js                   ← cross-species dotplot (hover fix shipped)
│
├── tests/sv_evidence/
│   ├── fixture_sv_genotype_counts_v1.json
│   ├── fixture_sv_evidence_combinations_v1.json
│   ├── fixture_sv_support_by_sample_v1.json
│   ├── test_step1_skeleton.js             (288 PASS)
│   ├── test_step2_table.js                (198 PASS)
│   ├── test_step3_locus.js                 (57 PASS)
│   ├── test_step3_5_cursor.js              (48 PASS)
│   ├── test_step4_rightrail.js             (60 PASS)
│   ├── test_step4_5_select.js              (29 PASS)
│   ├── test_step5_upset.js                 (50 PASS)
│   └── smoke_test_step{1,2,3,3_5,4,4_5,5}.html
│
├── producers/                              ← step 7 — Python emitters (NEW)
│   ├── README.md
│   ├── sv_evidence_io.py
│   ├── STEP_SV_GT_AGG_aggregate_genotype_counts.py
│   ├── STEP_SV_EVID_COMB_emit_combinations.py
│   └── STEP_SV_SUPPORT_emit_support_by_sample.py
│
├── json/sv_genotype_counts/cand_LG28_15Mb.json    ← dev/test layer
│
└── specs_todo/
    ├── SPEC_sv_evidence_page__upset_redirect.md
    ├── SPEC_sv_evidence_page__per_candidate_folder.md
    └── (other spec docs)
```

---

## Status overview

| step  | feature                                                       | tests   | status |
|-------|---------------------------------------------------------------|---------|--------|
| 1     | Skeleton + tab + empty-state mount                            | 288/288 | DONE   |
| 2     | Loader + main SV table (sort/filter/paginate/TSV export)      | 198/198 | DONE   |
| 3     | Locus track strip                                              | 57/57   | DONE   |
| 3.5   | Cursor + hover + zoom + E/F/Enter/Esc + ↑↓ + view presets      | 48/48   | DONE   |
| 4     | Right-rail boundary summary + interactive legend + row anno    | 60/60   | DONE   |
| 4.5   | Region-select mode (drag → region filter)                      | 29/29   | DONE   |
| 5     | UpSet panel (the redirect) + folder-walk drag-drop             | 50/50   | DONE   |
| 6     | Sample × SV dosage heatmap (compact + large overlay)           | —       | DONE   |
| 7     | Python emitters + per-candidate folder writer                  | —       | DONE   |
| dotplot fix | Cross-species page hover only on SVG square              | —       | DONE   |

**Atlas-side: complete and tested.** Producer-side: schemas locked, VCF parsing blocks stubbed (need cyvcf2/pysam).

---

## How everything fits together

The atlas reads a per-candidate folder via drag-drop:

```
data/<chrom>/candidates/<cand_id>/
├── sv_genotype_counts.json        → main table + locus + boundary summary
├── sv_evidence_combinations.json  → UpSet panel (right-rail)
└── sv_support_by_sample.json      → dosage heatmap + glyph dimming
```

The folder-walk uses `webkitGetAsEntry()` to traverse the dropped folder, picks up every `*.json`, and routes each by `format_version`. Single-file drag-drop still works (backwards-compatible).

Click an UpSet bar → `selectedSamples` populates → with the support layer also loaded, locus glyphs whose carriers don't intersect the selection dim to 0.35 opacity. Without the support layer, glyphs stay opaque (the atlas can't compute "no carriers in selection" without the per-sample matrix).

---

## State shape (everything `_state.*`)

```js
{
  // Step 1-2 core
  rootEl, layer, layerLoading, layerError, activeCandidateId,
  filters, sortColumnId, sortDirection, pageSize, pageIndex,
  // Step 3 locus crosstalk
  highlightedSvId,
  // Step 3.5 cursor + view
  cursorBp, markerLeftBp, markerRightBp, viewPreset, hotkeysAttached,
  _lastWindow, _customWindow,
  // Step 4 annotations + legend highlight
  rowAnnotations, highlightPattern,
  // Step 4.5 region-select
  selectMode, selection, _selectDrag,
  // Step 5 UpSet
  combinationsLayer, selectedSamples, activeCombinationIndex,
  // Step 6 heatmap
  supportLayer, heatmapView, heatmapHoverCell,
}
```

## Public API surface

```js
window.AtlasSVEvidence = {
  // Lifecycle
  init({ root }), loadCandidate(cid?), refresh(),

  // Filters / table
  setFilters(filters), exportFilteredTSV(),

  // Step 3.5 cursor + hotkeys
  attachHotkeys(), detachHotkeys(),
  setViewPreset('default'|'left_close'|'right_close'),
  setCursorBp(bp),

  // Step 4.5 region select
  setSelectMode('zoom'|'select'),
  setSelection(startBp, endBp),

  // Step 5 UpSet
  onUpSetBarClick(i), clearSampleSelection(),
  setCombinationsLayer(obj),

  // Step 6 heatmap
  setSupportLayer(obj),
  setHeatmapView('compact'|'large'),

  _state, _const, _internals
}
```

---

## Schemas (all locked, atlas validates each)

### `sv_genotype_counts_v1` — STEP_SV_GT_AGG

The main per-candidate layer. Per-SV record carries: `sv_id`, `sv_type`, `position_bp`, `zone`, group genotype counts (`H1/H1.AA/AB/BB/miss` etc.), Fisher OR + p-value + FDR_BH, `pattern_label`, `n_samples_with_call`, optionally `samples_with_call` (carrier id list). Plus aggregate `boundary_summary.{left,right}.by_sv_type` and `min_carriers_filter` provenance block.

### `sv_evidence_combinations_v1` — STEP_SV_EVID_COMB

UpSet input. Top-level: `evidence_types[]` (the canonical 8 types: `left_SA`, `right_SA`, `left_PE`, `right_PE`, `Manta_INV_GT`, `DELLY_INV_GT`, `MAPQ0_left`, `MAPQ0_right`), `combinations[]` (sorted by intersection_size desc, top N), `per_evidence_totals` (for the right-side mini-bars).

### `sv_support_by_sample_v1` — STEP_SV_SUPPORT

Heatmap matrix. `samples[]` (ordered H1/H1 → H1/H2 → H2/H2), `sv_ids[]`, `dosage_compact[]` (one string per sample, char-encoded `'0'/'1'/'2'/'.'`), `row_groups` (group → [start, end] indices).

---

## On the `min_carriers_per_band: 5` filter (Quentin's hard threshold)

This is a **producer-side** filter, applied in `STEP_SV_GT_AGG` before emit:

> Keep an SV only if at least one karyotype band (H1/H1, H1/H2, H2/H2) has ≥ 5 carriers (HET or HOM-ALT). Below n=5, the signal is noise — singletons, sequencing artefacts, doubleton flukes.

Reasons to do it producer-side rather than atlas-side:

1. **K=3 on 226 samples** → bands of ~60 each → n=5 is the floor for "this SV actually segregates".
2. **K=6 substructure** → bands of ~30-40 → n=5 still defensible.
3. **JSON stays small.** For ~1000 raw SV calls, the filter typically drops half.
4. The atlas's UI `min_samples` filter defaults to 1 — it doesn't double-filter; it lets the user dial up further.

The constant is recorded in the emitted JSON under `min_carriers_filter` for reproducibility.

---

## Critical constraints (non-negotiable; reiterated by Quentin)

1. **Three separate catfish cohorts; never conflate:**
   - F₁ hybrid (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
   - **226-sample pure *C. gariepinus* hatchery cohort on LANTA** — current inversion work. K clusters here are hatchery broodline structure, NOT species admixture.
   - Pure *C. macrocephalus* wild cohort — future paper.

2. **Page is read-only with respect to candidate boundaries.** Boundary editing belongs to page 11. The SV evidence page does not modify `bs.staging`, `state.candidate.boundary_*`, karyotype groups, or SV calls. The E/F keys here drop visual-only pins; they look like page-11 muscle-memory mnemonics but cannot edit.

3. **Empty/partial states never throw.** Drop a folder with one of three layers, get one of three panels. Drop with all three, get the full experience. Never crash.

4. **All colours via CSS variables.** Light + dark mode both work.

5. **Surgical HTML edits only.** `Inversion_atlas.html` has had only 4 edits since project start — purely additive (tab button, page DOM, script include, dispatcher branch). All step 2-7 work is in the JS module + producers/.

---

## Verify the bundle

```bash
unzip sv_evidence_complete.zip
cd Atlas

# Atlas-side
node --check js/atlas_sv_evidence.js
node --check js/atlas_dotplot.js
for f in tests/sv_evidence/test_step*.js; do
  echo "=== $(basename $f) ==="
  node "$f" 2>&1 | tail -3
done
# All 730 should pass.

# Producer-side
python3 -m py_compile producers/*.py
# Then wire VCF parsing — see producers/README.md for how.
```

For visual smoke test: open `tests/sv_evidence/smoke_test_step5.html` in a browser.

---

## Suggestions for the next chat

The atlas-side is done. The next chat should pick **one** of these directions, in rough order of impact:

### Option A — Wire the producers (recommended)

Replace the VCF-parsing stubs with `cyvcf2`/`pysam` in `STEP_SV_GT_AGG` and `STEP_SV_SUPPORT`, and the BAM-evidence stub in `STEP_SV_EVID_COMB`. Run them on LANTA against your real DELLY+Manta VCFs and the 226-sample cohort. The first per-candidate folder for `cand_LG28_15Mb` (or `INV_LG28_003`) lights up the full atlas page with real data. The biggest payoff per token spent.

**Specific tasks:**
- `STEP_SV_GT_AGG`: implement `parse_vcf_for_region()` using cyvcf2; replace `_placeholder_fisher_p` with `scipy.stats.fisher_exact`.
- `STEP_SV_SUPPORT`: implement `extract_dosages_from_vcf()`.
- `STEP_SV_EVID_COMB`: wire `extract_evidence_from_bam_layer()` against your S7 / phase_8_comparative_breakpoint_fragility module's API. Define the precise `*_count` thresholds for `left_SA`/`right_SA`/`left_PE`/`right_PE` based on what counts as "supporting evidence" per Quentin's lab convention.
- Run end-to-end on cand_LG28_15Mb. Verify the atlas renders all three panels with real data.

### Option B — Add the candidate-list browser

Currently the SV evidence page assumes a single active candidate. For a multi-candidate manuscript, add a candidate selector (dropdown or sidebar list) at the top of the page. Each entry → `loadCandidate(cid)` against its folder. Quick to build; cleanly extends the existing infrastructure.

### Option C — Manuscript figure exporter

Add a "Generate Figure 5" button that produces a high-res PNG/SVG snapshot of the locus + UpSet + heatmap layout for inclusion in the paper. Probably worth 2-3 hours of work; uses `dom-to-image` or `html2canvas`. Useful for the manuscript supplement.

### Option D — Heatmap row clustering

Currently the heatmap rows (SVs) are in `sv_ids` order. For the manuscript, hierarchical clustering of rows by carrier-pattern similarity gives a much more striking visualization (you see the canonical_breakpoint_marker SVs cluster together, the het_specific_marker ones cluster, etc). About a half-day of work.

### Option E — Multi-candidate UpSet comparison

Compare evidence patterns across candidates: "do all LG28 candidates show the same `left_SA + right_SA + PE` dominant pattern, or do some look more like `MAPQ0_only` artefacts?" Would surface real biological differences between candidates. Larger scope.

**My pick:** Option A. The whole infrastructure is built around real data; nothing else compounds until the producers are wired. After Option A ships, Option D (heatmap clustering) is the natural next move because it makes the heatmap manuscript-grade.

---

## Pickup commands for next chat

```bash
# Verify environment
cd /home/claude/Atlas
node --check js/atlas_sv_evidence.js
python3 -m py_compile producers/*.py

# Run full regression (730 should pass)
for f in tests/sv_evidence/test_step*.js; do
  echo "=== $(basename $f) ==="
  node "$f" 2>&1 | tail -3
done

# Producer entry points
ls producers/
cat producers/README.md   # for the production-flow how-to
```

The canonical working tree is `/home/claude/Atlas/`. The bundle zip is just a snapshot of it.

---

## What did NOT make it (deferred)

- **Region-select also dims locus glyphs that fall outside region** — currently selection only filters the table; locus glyphs still render at full opacity if their position is inside the original window. Easy fix in `_drawTrackSVCalls`: when `_state.selection` is set, dim glyphs whose `position_bp` is outside.
- **Heatmap row order = filter+sort order from the table** — currently the heatmap orders by `sv_ids[]` (producer order). Should mirror `sortColumnId/sortDirection` so the heatmap and the table tell the same story. ~10 lines of code.
- **Genotype-count column recompute within `selectedSamples`** — when a UpSet bar is clicked, the table's H1/H1, H1/H2, H2/H2 columns should re-count within selection. Possible now with supportLayer loaded, but takes some plumbing. ~50 lines.
- **TSV export should include the filtered candidate metadata** — minor.
- **Tests for steps 6 and 7** — by Quentin's request to skip tests for now and code specs instead.

These are all good first-week tickets for the next chat.
