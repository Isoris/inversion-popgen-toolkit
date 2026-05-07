# Migration inventory — splitting the legacy Atlas

**Source**: `Inversion_atlas.html` (the legacy single-file, ~76k LOC,
turn 165 final + turn 166 spec extension; the spec extension was
spec-only, no code change, so the binary is the turn 165 close).

**Target**: four sub-atlases per the canonical naming:
- `inversion_discovery.html` + `inversion_discovery/`
- `inversion_review.html`    + `inversion_review/`
- `inversion_catalogue.html` + `inversion_catalogue/`
- `inversion_comparative.html` + `inversion_comparative/`

Plus shared modules under `shared/` for primitives used by ≥2 sub-atlases.

**Rule**: nothing in the legacy file gets dropped. Every page, every
shipped turn, every active spec must land somewhere on this inventory.
If you can't find an item on the list below, the inventory is
incomplete — flag it before extraction proceeds.

---

## 1. Pages — every `data-page=...` tab in the legacy atlas

22 pages total. The split follows the existing `data-stage` attribute
on each tab.

| # | data-page | Tab label | Stage (legacy) | Lands in | Notes |
|---|---|---|---|---|---|
| 1 | `page1` | local PCA \|z\| | discovery | `inversion_discovery` | The main scrubber. Sim_mat heatmap, robust \|Z\|, per-sample lines, K-means PCA, L3 contingency. |
| 2 | `page12` | local PCA θπ | discovery | `inversion_discovery` | Mirror of page1 but on θπ. Empty-state until R-side `STEP_R39/R40/R41` ship. |
| 3 | `page15` | local PCA GHSL | discovery | `inversion_discovery` | Mirror of page1 but on GHSL. Empty-state until R-side `STEP_C04c/d` ship. |
| 4 | `page2` | candidate focus | discovery | `inversion_discovery` | Deep-dive on one promoted candidate (sim_mini, karyogram, dosage heatmap, ancestry strip). |
| 5 | `page8` | windows | refinement | `inversion_discovery` | Per-window summary table. **Note**: legacy stage is "refinement" but it's a discovery diagnostic — moving to discovery. |
| 6 | `page19` | negative regions | discovery | `inversion_discovery` | Inversion-negative inventory (complement of the catalogue). |
| 7 | `page11` | boundaries | refinement | `inversion_review` | Refine each candidate's [start_bp, end_bp] interval. |
| 8 | `page_sv_evidence` | SV evidence (5b) | refinement | `inversion_review` | SV calls clustered around boundaries; karyotype-group genotype counts. |
| 9 | `page4` | karyotype / tier | refinement | `inversion_review` | Two views: per-candidate sample-level regime + 14-axis tier classification. |
| 10 | `page7` | ancestry | classification | `inversion_review` | Per-window ancestry view (K-cluster + delta12 heatmap). |
| 11 | `page6` | popstats | classification | `inversion_review` | Stack of pop-genetic tracks (Z, SNP density, BEAGLE, coverage, θπ, FST, Hobs/Hexp). |
| 12 | **NEW** | inversion review | refinement | `inversion_review` | **NEW** — the SPEC BLOCK 2 page (auto-promote, bulk actions, sample-concordance candidate proposals). Adds, doesn't replace existing pages. |
| 13 | `page3` | catalogue | discovery | `inversion_catalogue` | Sortable/filterable catalogue of all L2 envelopes / L1-merged inversions. **Legacy stage is "discovery" but functionally it's synthesis** — moving to catalogue. |
| 14 | `page9` | confirmed | synthesis | `inversion_catalogue` | Carousel walk-through of confirmed candidates. |
| 15 | `page21` | annotation cockpit | synthesis | `inversion_catalogue` | Enlarged per-sample-lines panel with cursor-driven candidate/band selection. |
| 16 | `page10` | markers | synthesis | `inversion_catalogue` | Diagnostic PCR marker panels per candidate. Activates with `marker_panel_summary`. |
| 17 | `page17` | stats profile | synthesis | `inversion_catalogue` | Statistical profile of inversion-associated genomic features. |
| 18 | `page18` | marker readiness | synthesis | `inversion_catalogue` | Marker readiness panel (private-indel architecture, tier 1-4). |
| 19 | `page_overview` | overview | synthesis | `inversion_catalogue` | One-page summary of the inversion atlas. |
| 20 | `page16` | cross-species breakpoints | compare | `inversion_comparative` | Cross-species breakpoints between Cgar and Cmac. |
| 21 | `page16b` | multi-species classification | compare | `inversion_comparative` | Multi-species classification cockpit. |
| 22 | `page5` | help | help | `inversion_comparative` | Vocabulary/hotkeys/pipeline. **Lives anywhere; parking in comparative.** |

### Decisions / disagreements with legacy stage labels

- **`page8` (windows)**: legacy stage `refinement`; moving to `inversion_discovery`. Justification: it's a per-window diagnostic table the user opens DURING scrubbing to investigate a window's robust |Z| and eigenvalue ratios. It's not a refinement step.
- **`page3` (catalogue)**: legacy stage `discovery`; moving to `inversion_catalogue`. Justification: it's the manuscript-grade output table; users build it AFTER discovery, not during.
- **`page6` (popstats)** and **`page7` (ancestry)**: legacy stage `classification`; moving to `inversion_review`. Justification: they're per-candidate evidence views the user consults while reviewing whether to confirm. Could also live in `inversion_catalogue` for cross-candidate stats — flagging for re-decision once we touch them.

If any of these decisions look wrong, flag now. They're cheap to swap pre-extraction.

---

## 2. Shipped turn-by-turn features — what's in the binary today

Distilled from the handoff archive (`HANDOFF_2026-05-05_turn164_*.md`,
`turn165_*.md`, etc.) and turn-marker grep on the binary. Listed by
turn so we don't lose anything.

### Pre-turn-100 era (foundation)
The bulk of the scrubber, PCA, L3 contingency, candidate strip,
manuscript-grade vocabulary discipline, Hungarian projection,
candidate registry, and per-page UIs were established turns 1–100.

| Turn | Feature | Lands in |
|---|---|---|
| early | `_buildContingency`, `_ssContingencyTableHtml`, K×K table renderer | `shared/contingency.js` |
| early | `_hungarianChainProjection`, `_concordanceMatrix` | `shared/hungarian.js` |
| early | K-means impl (`getL2Cluster`, `getL2ClusterAt`, `getL2ClusterByMode`) | `shared/kmeans.js` (factor out) |
| early | `state` global object — single source of truth for current candidate, K, layers, scrubber pos | each sub-atlas re-establishes its own; cross-atlas via JSON snapshots |
| early | DOM scaffolding for tab bar, side bar, sim_mat, |Z| panel, lines panel, PCA panel, L3 panel | split per-page across sub-atlases |
| early | candidate strip + auto-promote machinery | mostly `inversion_discovery` |
| early | `state.candidateList`, `state.candidate`, candidate registry | `shared/candidate_registry.js` |
| early | localStorage keys (`pca_scrubber_v3.*`) | each sub-atlas namespaces its own |

### Turn 101–129 (popgen depth)

| Turn | Feature | Lands in |
|---|---|---|
| 101 | structural-haplotype transition graph | `inversion_review/transition_graph.js` |
| 103 | band-reach + regime-breadth per L2 envelope | `inversion_discovery/band_reach.js` |
| 104 | regime-breadth strip on PC1 sub-panel | `inversion_discovery/pc1_strip.js` |
| 114 | cross-species breakpoint overlay vertical lines | `inversion_comparative/breakpoint_overlay.js` |
| 115–127 | (many small UI/popgen turns — to be itemised when we touch them) | various |
| 128 | per-sample θπ lines mode plumbing | `inversion_discovery` (and `shared/het_rate.js`) |
| 129 | `_HET_RAMP`, `_hetRateColor`, `_computeHetRateForL2` | `shared/het_rate.js` |

### Turn 130 (review surfaces foundation)

| Turn | Feature | Lands in |
|---|---|---|
| 130 | `_isAutoCandidate`, lineage compute foundation, dashed-outline auto candidates | `shared/candidate_registry.js` (the predicate); UI in `inversion_review` |
| 130 | `_gatherActiveCandidatesForInheritance` | `inversion_review` |
| 130 | `state.lineageResult` foundation (Slice 3 of `SPEC_review_surfaces_auto_and_lineages.md` deferred) | `inversion_review` |
| 130 follow-up | filter-modifiers shipped (the auto-candidate-aware filter) | `inversion_review` |

### Turn 133–134 (L2-sweep auto-promote)

| Turn | Feature | Lands in |
|---|---|---|
| 133 | L2-sweep auto-promote pipeline (`source = 'auto_l2_sweep'`, `confirmed = false`) | `shared/candidate_registry.js` (compute) + `inversion_discovery` (UI) |
| 134 | `_AUTO_PROMOTE_MIN_BAND_SIZE` and related thresholds | `shared/candidate_registry.js` |

### Turn 135–141 (band-trace + lines panel)

| Turn | Feature | Lands in |
|---|---|---|
| 138–140 | various band-trace plumbing (TBD: extract specifics on touch) | `inversion_review/band_trace.js` |
| 141 | `_paintCandidateBands` (lines-panel candidate alpha intervals) | `inversion_review/lines_panel.js` (also reused by `inversion_discovery`) |

### Turn 142–157 (G-panel + V-shape modal + cross-candidate matrix)

| Turn | Feature | Lands in |
|---|---|---|
| 152 | G-panel inheritance tab (cross-candidate matrix) | `inversion_review/g_panel.js` |
| 156 | V-shape modal pattern | `shared/modal.js` (reused by `inversion_review` chip popovers per SPEC BLOCK 2 §7.3) |

### Turn 159–164 (band-trace pipeline complete)

| Turn | Feature | Lands in |
|---|---|---|
| 160 | `_bandTraceForFishSet`, `_bandTraceRegimeRuns` (single-cohort tracker, L2-envelope granularity) | `shared/band_trace.js` |
| 161 | `setBandTraceFishSet`, band-trace UI (Slice 4 UI half), `_drawBandTraceStrip` | `shared/band_trace.js` (compute) + `inversion_review/band_trace_panel.js` (UI) |
| 162 | per-L2 hit rectangles for hover tooltip | `shared/band_trace.js` |
| 163 | chain-break tick color | `shared/band_trace.js` |
| 164 | lasso-linkage cache invalidation | `shared/candidate_registry.js` |

### Turn 165 (G-panel `auto` review tab)

| Turn | Feature | Lands in |
|---|---|---|
| 165 | `_GPANEL_TABS` extended with `auto` entry, `_gPanelHasAutoCandidates` predicate | `inversion_review/g_panel.js` |
| 165 | per-row Confirm/Dismiss/Inspect actions for auto candidates, bulk actions | `inversion_review/g_panel.js` |
| 165 follow-up | (Slice 1 of `SPEC_review_surfaces_auto_and_lineages.md` discovered already-shipped) | (no new code, just verified) |

### Turn 166 (spec-only, no code change)

Spec extension to `SPEC_band_track_extraction_and_l3_single_band_rows.md`
adding §1.4–1.5, §4A (het-anchored mode), §6 TSVs E and F, §7.4–7.5,
§8 slices update, §10 tests update. No binary change.

### Stub and scaffold work that LOOKS shipped but is actually stubbed

These need to be migrated as stubs (not as completed features) so the
sub-atlas inheriting them knows the body still needs writing.

| Item | What's there | What's stubbed | Lands in |
|---|---|---|---|
| `_resolveSampleColorByMode` (v3.99 turn 14e+) | UI scaffold (mode picker dropdown), state slot, validator, greyed-out tooltips | the per-mode renderer bodies (stub returns null) | `shared/lines_panel_color.js` (the resolver hook) + `inversion_discovery/lines_panel.js` (the picker UI) |
| `_resolveSampleScopeColor` (v3.99 turn 14e+) | same as above | stub returns null | same as above |
| `theta_pi_per_window` layer detection | `detectSchemaAndLayers` knows the name; mode picker option auto-enables when layer arrives | layer never arrives until R-side `STEP_R36` ships | `shared/state_io.js` (already done) |
| `_bundleTierBlock` (truncated in our copy) | reads `state.data.final_classification`; renders 14-axis table | partial (truncated mid-function in the tarball; recover from full file when uploaded) | `inversion_review/tier_block.js` |
| Page12 (θπ scrubber) | full empty-state UI with per-layer status indicator | empty until R-side `STEP_R39/R40/R41` ships `theta_pi_*` layers | `inversion_discovery/page12_theta_pi.js` |
| Page15 (GHSL scrubber) | full empty-state UI | empty until `STEP_C04c/d`/`export_ghsl_to_json_v3` ships | `inversion_discovery/page15_ghsl.js` |
| Page10 (markers) | UI hidden until `marker_panel_summary` loaded | activates on layer load | `inversion_catalogue/page10_markers.js` |

---

## 3. Active specs — what's still pending

Some specs are shipped, some half-shipped, some not started. These
need to land somewhere so the work isn't forgotten.

| Spec | Status | Owner sub-atlas |
|---|---|---|
| `SPEC_band_track_extraction_and_l3_single_band_rows.md` (1769 LOC) | unshipped (3 slices to go) | `inversion_review` (slice 1: per-band rows in L3 contingency); `inversion_review` + `shared/band_tracks.js` (slice 2: track extraction + skeleton + lines coloring); `inversion_review` (slice 3: TSV exports + classification chip) |
| `SPEC_distant_band_concordance_fish_trajectory.md` | partly shipped (slices 1+2 in turn 130) | `inversion_review/fish_trajectory.js` (remaining slices) |
| `SPEC_l2_sweep_inheritance.md` | shipped (turn 133–134) | `inversion_discovery` (compute) + `shared/candidate_registry.js` |
| `SPEC_review_surfaces_auto_and_lineages.md` | slice 0+1+2 shipped (turn 130/165); slice 3 (lineages tab) deferred | `inversion_review/g_panel.js` |
| `SPEC_lines_panel_candidate_bands.md` | slice 1 shipped (turn 141); slice 2 alpha intervals also shipped | `inversion_review/lines_panel.js` |
| `SPEC_g_panel_unified_groups.md` (karyotype tab) | shipped (turn 152) | `inversion_review/g_panel.js` |

### NEW spec from this chat — SPEC BLOCK 1 (the R-module)

Quentin's `SPEC BLOCK 1` (the band-continuity / het-backbone candidate
proposer) is an HPC R-module, NOT atlas-side code. It produces the 8
output files that land in `data/precomp/<chrom>/`:

```
LG28.band_nodes.tsv
LG28.band_edges.tsv
LG28.transition_events.tsv
LG28.band_trajectories.tsv
LG28.het_band_backbones.tsv
LG28.candidate_track_proposals.tsv
LG28.candidate_tracks.json
LG28.manual_review_queue.tsv
```

Lives in `inversion_codebase_v8.5/MODULE_BAND_TRACKS/` on LANTA. Not
part of the atlas-side split, but referenced here because
`inversion_review` consumes its outputs.

### NEW spec from this chat — SPEC BLOCK 2 (the Inversion Review page)

Quentin's `SPEC BLOCK 2` (the new Inversion Review page) — 5 panels,
auto-promote, bulk actions, sample-concordance proposals. Lands in
`inversion_review`. Built on top of (= imports from) the existing
shared primitives.

---

## 4. Shared primitives (extracted to `shared/`)

Anything used by ≥2 sub-atlases moves to `shared/`. First-pass
inventory; finalised when we touch the legacy file:

```
shared/
├── state_io.js            ✓ shipped (turn N — this chat)
├── contingency.js         _buildContingency, _ssContingencyTableHtml, partition primitives
├── hungarian.js           _hungarianChainProjection, _concordanceMatrix
├── kmeans.js              getL2Cluster, getL2ClusterAt, getL2ClusterByMode, K-means impl
├── pca.js                 PCA solver used by the three scrubbers
├── het_rate.js            _HET_RAMP, _hetRateColor, _computeHetRateForL2 (turn 129)
├── band_trace.js          _bandTraceForFishSet, _bandTraceRegimeRuns (turns 160-164)
├── band_tracks.js         (NEW) _buildBandTracksForCandidate, _extractHetSkeleton (SPEC BLOCK 2)
├── candidate_registry.js  state.candidateList, _isAutoCandidate, lineage compute, L2-sweep auto-promote
├── modal.js               V-shape modal pattern (turn 156); chip popovers
├── lines_panel_color.js   _resolveSampleColorByMode, _resolveSampleScopeColor (with their stubs)
├── color_ramps.js         all the color scales (warm/cold, K-band palette, etc.)
└── dom_utils.js           shared DOM helpers (createElement helper, throttle, etc.)
```

## 5. JSON layers — who reads/writes what

Cross-atlas JSON-roundtrip contract (so cross-atlas state plumbing has
explicit producer/consumer agreements):

| File | Workflow | Producer | Consumers |
|---|---|---|---|
| `data/precomp/<chr>/<chr>.json` | core | upstream R pipeline | all four sub-atlases |
| `data/precomp/<chr>/<chr>.repeat_density.scrubber_windows.json` | core | upstream R pipeline | all four |
| `data/precomp/<chr>/<chr>.{band_nodes,band_edges,...}` (8 files) | inversion | SPEC BLOCK 1 R-module | `inversion_review` |
| `data/cohort/relatedness.json` | core | upstream R pipeline | all four |
| `data/cohort/cohort_diversity_v1.json` | core | upstream R pipeline | all four |
| `data/cohort/sample_froh.json` | core | upstream R pipeline | all four |
| `data/candidates/<cid>/sv_genotype_counts.json` | inversion | upstream R pipeline | `inversion_review`, `inversion_catalogue` |
| `data/candidates/<cid>/boundaries_refined.json` | inversion | upstream R pipeline OR `inversion_review` (when user manually refines) | `inversion_review`, `inversion_catalogue` |
| `data/candidates/<cid>/final_classification.json` | inversion | upstream R pipeline | `inversion_catalogue` |
| `data/comparative/cs_breakpoints_v1.json` | inversion | upstream R pipeline | `inversion_comparative` |
| `data/review/inversion/manual_overrides.json` | inversion | `inversion_discovery`, `inversion_review` (any user override) | `inversion_review` (read on next session) |
| `data/review/inversion/candidate_review_decisions.json` | inversion | `inversion_review` (accept/reject/split/merge) | `inversion_catalogue` (only confirmed shows in catalogue) |
| `data/review/inversion/locked_karyotype_groups.json` | inversion | `inversion_review` (karyotype/tier tab) | `inversion_catalogue`, `inversion_comparative` |
| `data/review/inversion/confirmed_candidates.json` | inversion | `inversion_catalogue` (re-derives from decisions) | manuscript export |

Cross-atlas state contract: state is **only** shared via files under
`data/`. No global JS object. Sub-atlases that need to know about
cross-atlas decisions reload the relevant file on (a) page open, and
(b) when the user clicks a "refresh" button.

---

## 6. Extraction order (proposed)

Surgical risk-managed sequence. Each step ships an atlas that runs
end-to-end before moving on. **One step per turn** is realistic.

| Step | Turn | Deliverable | Risk |
|---|---|---|---|
| **0** | done | foundation (shared/state_io.js, build/, data/ skeleton, smoke tests) | done |
| **1** | next | extract `shared/` primitives (contingency, hungarian, kmeans, pca, het_rate, color_ramps) — JS only, no UI yet. Each module exports + has unit tests. | medium (need to identify what's used by ≥2 callers) |
| **2** | +1 | `inversion_discovery.html` + `inversion_discovery/` — page1 + page12 + page15 + page2 + page8 + page19. Behavioural parity with legacy. | high (the scrubber is the bulk of the legacy code) |
| **3** | +2 | `inversion_review.html` + `inversion_review/` — page11 + page_sv_evidence + page4 + page7 + page6 + the existing band-trace UI + G-panel auto tab. Migration only, no new SPEC BLOCK 2 features yet. | high (cross-page state plumbing kicks in) |
| **4** | +3 | `inversion_catalogue.html` + `inversion_catalogue/` — page3 + page9 + page21 + page17 + page18 + page10 + page_overview. | medium |
| **5** | +4 | `inversion_comparative.html` + `inversion_comparative/` — page16 + page16b + page5. | low |
| **6** | +5 | parity verification — load every JSON layer the legacy atlas loads, click every button, confirm same behaviour | low (testing) |
| **7** | +6 | NEW work begins: SPEC BLOCK 2 features (auto-promote / bulk actions / sample-concordance proposals) added to `inversion_review`, on top of the migrated baseline | self-contained |
| **8** | +7 | NEW work: SPEC BLOCK 1 R-module on LANTA, producing `data/precomp/<chr>/*.{band_nodes,...}` | self-contained |

Steps 1–6 are pure migration: nothing new, just sorted into the right
homes. Step 7+ is where the auto-promote / bulk-action functionality
finally lands.

If any step blows up the budget (e.g. step 2 turns out to be 3 turns),
we re-plan rather than skipping things.

---

## 7. What's NOT on this inventory and shouldn't be lost either

If you spot any of these missing from the lists above, flag it:

- All `state.*` slots (`state.k`, `state.activePage`, `state.candidate`,
  `state.candidateList`, `state.layersPresent`, `state.linesColorMode`,
  `state.bandTraceFishSet`, `state.bandTraceOn`, `state.lineageResult`,
  `state.gPanelTab`, `state.l3MoreExpanded`, `state.atlasMode`, ...)
  → these need a `shared/state.js` module
- All localStorage keys (`pca_scrubber_v3.*`, `inversion_atlas.*`)
  → preserved verbatim per sub-atlas
- All keyboard shortcuts (E/F for boundaries on page11, ←/→ for cursor
  on page21, etc.) → migrated with their pages
- All inline event handlers (`onclick="..."`) → must be migrated to
  `addEventListener` since modules don't expose top-level names to
  inline handlers
- All `window.foo = foo` debug exposures (`window._isAutoCandidate`,
  `window._computeHetRateForL2`, `window.setBandTraceFishSet`, ...)
  → preserved at module entry, exposed for console debugging
- `_BTRACE_*` constants (turn 161–163)
- All schema-version migration code in `detectSchemaAndLayers`
- The 14-axis tier classification block (`_bundleTierBlock`)
- All breeding-card export logic (`_breedingCardPrintHTML`)
- `_resolveSampleColorByMode` stub (when filled in) per turn 14e+ note
- `<style>` blocks (4,706 LOC of CSS) — split per sub-atlas; shared
  CSS goes to `shared/styles.css` linked from each HTML

---

## 8. Open questions / decisions needed before extraction starts

1. **CSS sharing strategy.** All 4 sub-atlases will reuse a lot of CSS
   (the panel layout, the candidate strip, the chip styles). Three
   options: (a) `shared/styles.css` linked from each HTML; (b) inline
   per-HTML; (c) per-component CSS in modules. **I recommend (a)** —
   one shared file plus a per-page additions file.
2. **`state` reconstruction.** The legacy `state` global has dozens of
   slots. After the split, each sub-atlas has its own `state` (different
   in-memory object). For cross-atlas slots (`state.candidate`,
   `state.candidateList`), they round-trip via JSON. **I recommend**:
   `shared/state.js` defines the canonical slot list + serializer +
   deserializer. Each sub-atlas instantiates one and loads it from
   `data/review/inversion/last_session.json` on startup.
3. **Page-tab navigation.** Today, clicking a tab swaps a panel within
   one HTML. After the split, clicking "discovery" from the review
   atlas should… open `inversion_discovery.html`? Or be a top-level
   nav we add? **I recommend**: a top-level `index.html` that lists
   the four sub-atlases as cards. Within each sub-atlas, the tab bar
   only shows that sub-atlas's pages.
4. **Tarball strategy after extraction.** Today the legacy atlas
   tarball is too big for the upload channel (LG28.json is 17 MB). The
   extracted version, with `data/precomp/LG28/LG28.json` excluded
   from source-tree tarballs, drops to ~2-3 MB total. **I recommend**:
   the source repo never includes `data/precomp/<chrom>/<chrom>.json`
   in tarballs; users fetch them from LANTA via a one-off script.

---

## 9. Sign-off

This inventory is the contract for extraction. Quentin reviews,
flags missing pieces, signs off. After sign-off, extraction proceeds
step-by-step per §6.

**Status**: DRAFT — awaiting review.

**Reviewer**: Quentin Andres
**Date drafted**: 2026-05-05
**Source**: turn 165 close binary (= what we have intact in the
truncated tarball, modulo the last ~800 LOC), plus turn 166 spec
extension (no binary change).
