# HANDOFF — full Atlas bundle — session 2026-05-05

**Bundled at**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (73,773 lines)
**Working dir at bundle time**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

This is a **full-tree bundle**, not a delta. The archive contains the
entire `/home/claude/Atlas/Atlas/` working directory at the close of
this session: `Inversion_atlas.html`, all 90 test files, all 41
handoff markdown files, all subdirectories (`tests/`, `js/`, `json/`,
`docs/`, `producers/`, `_scripts/`, `dosage_viewer/`, `engine_fast_ld/`,
`server_turn1/`, `specs_new_turn131/`, `specs_todo/`, `demo/`,
`changelogs/`, `handoffs_to_implement/`, `patches/`).

---

## 0. Full-suite status at bundle time

| metric                       | value                |
|---                           |---                   |
| `Inversion_atlas.html` LOC   | 73,773               |
| Test files (`test_turn*.js`) | 60                   |
| Total assertions             | **2991 / 0**         |
| JS syntax check              | clean                |
| HTML parser                  | 0 errors             |
| Bundle size                  | ~1 MB compressed     |

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly
   paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse,
pragmatic. Tarball is the standard handoff format.

---

## 2. Session arc — turns 151 through 158

This session ran across multiple parallel chats with Quentin running
simultaneous tracks. I (this chat) shipped odd-numbered turns; a
parallel track shipped some even-numbered turns and one of the two
turn-157 deliverables. Both tracks committed into the same working
tree. The summary below covers everything that landed.

### Turn 151 — Slab het-coloring parity (mine)
- Wired heterozygosity-rate coloring into `drawSlabMiniPCA` so slab
  L3 panes match the L2 panes shipped in turn 129.
- Refactored `_computeHetRateForL2` into a thin wrapper around new
  `_computeHetRateForRange(start_bp, end_bp, cacheKey)`.
- Added `_computeHetRateForSlab(start_w, end_w)` translating window
  indices to bp.
- +109 LOC. **2557 / 0** at end.

### Turn 152 — G-panel Slice 3 inheritance tab (mine)
- Filled the previously-stub Slice 3 (inheritance tab) with the real
  cross-candidate Jaccard cluster view.
- New state slot: `state.gPanelInheritanceThreshold` (default 0.15,
  persisted to localStorage).
- Eight new helpers: `_gpInhEnsureThreshold`, `_gpInhSetThreshold`,
  `_gpInhSiToCga`, `_gpInhMembersForGroup`, `_gpInhMakeManualGroup`,
  `_gpInhItemLabel`, `_gpInhRenderEmpty`, `_gpInhRenderResult`.
- Modal wiring: 220ms-debounced threshold slider, recompute button,
  matrix-view shortcut (opens existing `openInheritanceMatrix`),
  per-group `+ make manual` buttons.
- Cleaned stale "(slice N)" copy from header subtitle + tab strip.
- Reuses `runInheritanceCompute`, `inheritanceGroupClustering`,
  `_gatherActiveCandidatesForInheritance`, `addToManualGroup`,
  `openInheritanceMatrix` — all from turn 115 / 122. Pure wiring layer.
- +393 LOC. **2630 / 0** at end.

### Turn 153 — Auto-register inheritance on candidate-list mutations (mine)
- Hooked `_autoRegisterInheritanceOnCandidateChange()` into
  `persistCandidateList()`. Compares cache key vs expected; on
  mismatch, invalidates `state.inheritanceResult` +
  `state.inheritanceCacheKey`, then debounces a notify pass via
  `queueMicrotask` (with `setTimeout` fallback).
- New `_notifyInheritanceConsumers()` re-renders all open inheritance
  consumers: `requestRepaint` (lines-strip pills), G-panel modal if
  on inheritance tab, matrix popup if visible. Each guarded with
  try/catch — partial notify is OK.
- Coverage: every candidate-mutation path in the codebase reaches
  `persistCandidateList` (`addCandidateToList`, `removeCandidateFromList`,
  `removeCandidateFully`, `commitL3Draft`, JSON import, L2-sweep
  auto-promote, regime updates). One hook covers all.
- Idempotent: no-op when cache key matches (avoids spurious churn on
  favorites toggles in Save-Session payloads).
- +109 LOC. **2685 / 0** at end.

### Turn 154 — G-panel inheritance compute UX hardening (parallel track)
- New state: `state.inheritanceLastStatus` carrying
  `{ ok, reason, n_items, n_groups, threshold, error, message, at }`.
- Five status reasons captured by `runInheritanceCompute`:
  `computed` / `cached` / `insufficient_items` / `null_result` /
  `exception`. `inheritanceGroupClustering` call wrapped in try/catch
  so exceptions become status updates instead of bubbling.
- Tab body auto-computes on mount (cache-stale check, no `force: true`).
- Status pill in tab UI: green/●, grey/○, amber/⚠ depending on outcome.
- Distinct messages for 0-groups vs null-result vs never-computed.
- Recompute button click feedback — `computing…` flash + disable +
  re-render.
- +159 LOC. **2746 / 0** at end.

### Turn 155 — Threshold folded into inheritance cache key (mine)
- The honesty fix flagged across turns 152/154 §6: threshold was NOT
  in the cache key, meaning two computes at different thresholds
  shared the same cached entry.
- Extended `_inheritanceCacheKey(items, mode)` →
  `_inheritanceCacheKey(items, mode, threshold)`. Threshold optional;
  defaults through `state.gPanelInheritanceThreshold` then
  `_IGC_DEFAULT_COSINE_DIST_THRESHOLD`. Backwards-compatible — all 6
  existing callers keep working without code change.
- `runInheritanceCompute` resolves threshold *before* the cache key
  and forwards it explicitly into `inheritanceGroupClustering` via
  `Object.assign`.
- Default threshold in `runInheritanceCompute` now reads from
  `state.gPanelInheritanceThreshold` (the slider value), so
  lines-strip auto-fire and slider agree by construction. Single
  source of truth.
- 4-decimal rounding before key build to coalesce slider FP jitter.
- `_gpInhSetThreshold`'s `invalidateInheritanceCache` retained as
  defense-in-depth.
- +56 LOC. **2776 / 0** at end.

### Turn 156 — V-shape (u, agreement_fraction) per-candidate diagnostic (mine)
- Port of the FIG_C30 V-shape signature from STEP29's coherence
  pipeline. Per-sample scatter of `agreement_fraction` (from
  `computeStripeQuality`) against `u` (rotated PC1 from
  `_getOrComputeUVRotation`).
- Five new functions: `_buildVShapeData(candidate, chunk, sqRows)`,
  `_vShapeColor(group)`, `_drawVShapePlot(canvas, data, opts)`,
  `openVShapePlot(candidate)`, `closeVShapePlot()`. All on `window`.
- DPR-aware canvas renderer with axis labels, gridlines at agreement
  = 0/0.25/0.5/0.75/1, dashed reference lines at `het_u` / `hom1_u`
  / `hom2_u`, baseline at `agreement = 0.5`. Draw order
  `HOMO_1 → HOMO_2 → HET` so the dip is visible.
- Modal popup with Escape + click-outside, lazy DOM creation.
  Graceful states for missing chunk / missing stripe-quality /
  missing candidate.
- `🔍 V-shape diagnostic` button in the G-panel karyotype tab footer.
- Pure read-only diagnostic. No state mutation.
- +518 LOC. **2840 / 0** at end.

### Turn 157 — TWO parallel deliverables (different tracks)

**Turn 157A — V-shape tooltip on hover (mine)**
- Wired the `hits[]` array from turn 156 into a hover tooltip.
- Six new helpers: `_vsTooltipEnsureEl`, `_vsTooltipBuildHtml`,
  `_vsTooltipShow`, `_vsTooltipHide`, `_vsHitTest`,
  `_wireVShapeTooltip`. Module-level `_vsLastRender` snapshot
  populated inside `openVShapePlot`.
- Tooltip content: CGA + group color chip, u (4 dec), agreement (4
  dec), coherence_class (semantic color), stripe_quality tier
  (semantic color).
- **CSS-pixel coord space** for hit-testing (intentional, documented):
  the V-shape canvas uses `setTransform(dpr,…)` *before* drawing, so
  `hits[]` are CSS-pixel coords. Hover handler compares raw
  `ev.clientX - rect.left` deltas without re-applying DPR scaling
  (different from inheritance pill pattern).
- Closest-point hit test (handles overlapping HET-on-HOMO points in
  the dip region — first-wins would resolve to wrong group).
- `closeVShapePlot` also calls `_vsTooltipHide` so no orphan tooltip.
- +232 LOC. **2926 / 0** after.

**Turn 157B — dosage shim load fix + page6/page7 try/catch (parallel track)**
- Fixed dosage.js shim load (Quentin's reported issue #1).
- Added try/catch guards on page6/page7 (likely root cause for
  ancestry/popstats tabs blocking other tabs from clicking — issues
  #2 / #3).
- Triaged 8 issues from Quentin's bulk message; this turn ships #1 +
  #2/#3 partial.
- +58 LOC. Items #4–#8 deferred to subsequent turns with concrete scoping.

### Turn 158 — Popstats track categorization: QC vs popstats split + off-by-default (parallel track)
- Closes Quentin's items #5 and #6 from the turn-157B handoff.
- New `category` field on each `POPSTATS_TRACKS` entry: `'always'` |
  `'qc'` | `'popstats'`.
- Helper `_popstatsCategoryOf(track)` defaults to `'other'` for
  auto-discovered tracks (legacy behaviour preserved).
- `loadPopstatsView` now returns `{ hidden, shown }` instead of a
  bare hidden Set. Shown set stores explicit user opt-ins for
  qc/popstats categories.
- `savePopstatsView` accepts the same shape (or a bare Set for
  back-compat with any caller that hadn't migrated).
- `renderPopstatsPage`'s chip strip groups by category with section
  labels ("QC:", "popstats:"). Visibility per-category: `always` =
  show if data; `qc`/`popstats` = show iff explicitly shown; `other`
  = legacy.
- Click handler dispatches based on category — `qc`/`popstats` toggle
  the shown set; `always`/`other` toggle the legacy hidden set.
- **2991 / 0** at bundle time.

---

## 3. Cumulative session deltas

|  | LOC start | LOC end | Δ LOC | Tests start | Tests end | Δ tests |
|---|---|---|---|---|---|---|
| **Session total** | 72,090 | 73,773 | **+1,683** | 2504 | **2991** | **+487** |

JS sanity: clean throughout. HTML parser: 0 errors throughout. No
existing tests left in a broken state.

---

## 4. What's in the bundle

The full archive contains:

```
Atlas/
├── Inversion_atlas.html             ← main file, 73,773 lines
├── Diversity_atlas.html
├── Genome_atlas.html
├── HANDOFF*.md                      ← 41 handoffs (incl. all turn-15x)
├── DEPLOY_TO_LANTA.md
├── README.md (if present)
├── tests/
│   ├── test_turn*.js                ← 60 turn-numbered test suites
│   └── (additional helper subdirs)
├── js/                              ← stand-alone JS modules
├── json/                            ← analysis-side JSON layers
│   ├── comparative_breakpoint_fragility/
│   ├── ushape_evolution/
│   └── sv_genotype_counts/
├── _scripts/                        ← R-side / cluster-side scripts
├── docs/
├── changelogs/
├── specs_new_turn131/               ← active specs
├── specs_todo/                      ← deferred specs
├── handoffs_to_implement/
├── producers/                       ← producer-side helpers
├── dosage_viewer/                   ← dosage viewer assets
├── engine_fast_ld/                  ← LD engine
├── server_turn1/                    ← server-side Python
├── demo/
└── patches/
```

---

## 5. What Quentin should exercise this session

### From turns 151–157A (the inheritance / V-shape stack)

1. **Slab het coloring (turn 151)** — open a slab L3 pane on a region
   where you already trust a candidate's broodline structure; the
   het swatch should look the same as it does on the L2 panes.

2. **G-panel inheritance tab (turn 152)** — `g` to open, switch to
   inheritance. With ≥ 2 locked candidates, drag the threshold
   slider; per-group cards re-render. Click `+ make manual` on a
   group to promote it to the manual list.

3. **Auto-register on candidate change (turn 153)** — leave the
   G-panel open on inheritance tab. Promote a new candidate from page
   1. The tab should re-render itself with the new candidate in the
   group cards. Same with the matrix popup.

4. **Compute UX (turn 154)** — clicking ↻ recompute now flashes
   `computing…` + status pill below the slider. The pill text
   describes what compute did or didn't do (computed/cached/insufficient/
   null_result/exception).

5. **Threshold cache fix (turn 155)** — drag slider to a non-default
   value (e.g. 0.20). Click ↻ recompute → status: `computed`. Click
   again at same value → status: `cache hit`. Drag back to 0.15 and
   click again → status: `cache hit` (because keys match again).
   Pre-turn-155 this would have been silent stale-result.

6. **V-shape diagnostic (turn 156)** — page 1 promote → page 2
   dosage heatmap → "Compute stripe quality" → G-panel karyotype tab
   → 🔍 V-shape diagnostic. Look for the V dip at u≈het_u with high
   homo plateaus. The status footer reports per-group n + mean
   agreement.

7. **V-shape tooltip (turn 157A)** — hover any point in the V-shape
   plot. Tooltip shows CGA, group, u, agreement, coherence class,
   stripe-quality tier with semantic colors. Cursor switches to
   pointer over hits.

### From turns 157B–158 (the popstats / page-fixes track)

8. **Dosage shim load (turn 157B)** — should now load on page mount
   without 404.

9. **Tab navigation (turn 157B)** — clicking ancestry or popstats
   should no longer block other tabs.

10. **QC vs popstats split (turn 158)** — open the popstats page.
    Right-side selector should show two sections ("QC:", "popstats:")
    instead of one. Both off by default. Toggle on the tracks you
    actually want.

---

## 6. What's NOT done — queued options

From the standing queue at the end of this session:

- **A** — UV refactor on locked groups. ~200 LOC. Locked-labels
  caveat applies (burned us before).
- **C** — Per-sample-lines het coloring. ~300–500 LOC + perf
  concerns at 226 × 30k.
- **F** — Per-group "show fish" expand toggle in G-panel inheritance
  card. ~50 LOC. Useful for spot-checking before manual promotion.
- **H** — Cross-candidate V-shape gallery / sparkline matrix.
  Triages many candidates at once.
- **D** — Pivot based on what het / inheritance / V-shape have
  revealed on real data.

From the parallel-track turn-157B handoff (Quentin's items 4–8 not
yet shipped):
- **#4** — G-bar / I·g / candidate karyotypes consumed directly by
  ancestry & popstats panels. Spec deferred.
- **#7** — LD track double-heatmap (common vs rare karyotype groups).
- **#8** — Plot export (PNG / SVG) from the various canvases.

---

## 7. Test discipline log

This session followed the convention codified in turn 150's handoff:
**pre-existing tests that codified flagged-but-not-fixed behaviour
must be inverted when the underlying behaviour is corrected.**

Inversions made this session:
- Turn 152 inverted assertions in
  `test_turn135_g_panel_slice1.js` (placeholder text expectations).
- Turn 155 inverted the `_inheritanceCacheKey` 2-arg signature
  assertion in `test_turn153_inheritance_auto_register.js`.

Regex-tightening (not inversion) on:
- Turn 155: tightened the try/catch around-the-call regex in
  `test_turn154_compute_ux_hardening.js` after a comment was
  inserted between `try {` and the call.

No tests deleted. No tests left failing.

---

## 8. Files in the bundle (top level)

```
HANDOFF_FULL_BUNDLE_2026-05-05.md           ← this file
HANDOFF_2026-05-05_turn157_vshape_tooltip.md
HANDOFF_2026-05-05_turn157_dosage_bridge_load.md
HANDOFF_2026-05-05_turn156_vshape_diagnostic.md
HANDOFF_2026-05-05_turn155_threshold_in_cache_key.md
HANDOFF_2026-05-05_turn154_compute_ux_hardening.md
HANDOFF_2026-05-05_turn153_inheritance_auto_register.md
HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md
HANDOFF_2026-05-05_turn151_slab_het_coloring.md
HANDOFF_2026-05-05_turn151_session_handoff.md
... (32 more older handoffs preserved for reference)
```

Note: there is NO standalone turn 158 handoff in the tree — the
parallel track shipped turn 158 with code + tests but did not write a
dedicated `.md`. Coverage is inside `tests/test_turn158_popstats_qc_split.js`'s
top docstring and in this consolidated handoff §2.

---

## 9. Honest framing

**What's solid:**
- 2991 / 0 sweep. Clean JS + HTML throughout the session.
- The inheritance stack (turns 152–155) is now coherent: tab works,
  candidate changes auto-propagate, threshold is honestly part of the
  cache key, status pill surfaces what compute did.
- The V-shape stack (turns 156–157A) is a complete read-only
  diagnostic with hover detail, no state mutation, built on existing
  primitives (UV rotation, computeStripeQuality, dosage chunk).
- Popstats UI (turns 157B–158) addresses concrete user-reported
  blockers.

**What's risky:**
- **Two parallel tracks both used "turn 157" as a label.** Distinct
  features, distinct test files (`test_turn157_vshape_tooltip.js` vs
  `test_turn157_dosage_bridge_load.js`), distinct handoffs. No code
  collision detected — they touch different parts of the file —
  but future archaeology will be easier if Quentin's mental model
  treats them as 157A / 157B from now on.
- Items #4 / #7 / #8 from the turn-157B triage are still open.

**What's queued (full list):**
- **A, C, F, H, D** from my standing queue (see §6).
- **#4, #7, #8** from the parallel track (G-bar consumed by ancestry/
  popstats; LD double heatmap; plot export).

End of consolidated handoff.
