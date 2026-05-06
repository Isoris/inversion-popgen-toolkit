# HANDOFF — turn 131 final spec-writing pass + bundle (v2)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (64,243 lines, 1900/0 tests passing)
**Working dir**: `/home/claude/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC (account `lt200308`).
**This handoff supersedes** all prior `HANDOFF_*.md` files in this folder.

---

## 0. Critical reminders before starting

### Cohort discipline (NEVER conflate)
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — genome
   assembly paper only. NOT the data here.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

### User
**Quentin Andres** (Kasetsart University Bangkok). Never invent surname.
Communication style: terse, direct, pushes back precisely. Prefers
honest "this won't work" over speculative yes.

### Test runner pattern (preserved from turn 130)
- Tests live in `/home/claude/Atlas/tests/`
- Run with `node tests/<test_file>.js`
- Test path: `path.resolve(__dirname, '..', 'Inversion_atlas.html')`
- Output format MUST be `PASS: N   FAIL: N` for runner grep
- vm.createContext typed-array gotcha: use `.constructor.name === 'Float32Array'`

### File hygiene
- Backups: `.bak_pre_<feature>` before, `.bak_post_<feature>` after
- JS-string footgun: literal `</script>` → use `'<' + '/script>'`

---

## 1. Folder layout for specs (NEW after Quentin's feedback)

Two separate folders to keep the build queue clean:

```
specs_todo/             # Active build queue — specs already reviewed
                        # by Quentin, queued for slicing.
                        # 27 .md files. Includes 4 of the 4+1 inheritance
                        # methods (fish-trajectory partly shipped, L2-sweep,
                        # sliding-window, unification decision-gated).

specs_new_turn131/      # Pending review queue — 19 new specs from turn 131
                        # written in one pass from past-chat mining.
                        # Includes SPECS_TIER_INDEX.md and a README.
                        # When Quentin approves a spec for build, MOVE it
                        # to specs_todo/ (not copy).
```

**Why separate**: many existing `specs_todo/` specs overlap with shipped
atlas features. Mixing 19 new pending-review specs into that folder
would pollute the active build queue.

**Read `specs_new_turn131/SPECS_TIER_INDEX.md` first** — it sorts ALL
specs (across BOTH folders) by tier and gives recommended build
priority order.

---

## 2. The 4+1 inheritance methods — all spec'd, status confirmed

You asked: *"Our 4+1 methods in this chat for inheritance is it spec
I know that we coded 2 of 5 in slice."*

Confirmed: all 4+1 methods are spec'd. 2/5 shipped, 3/5 spec-only,
+1 decision-gated.

| # | Method | Spec location | Status |
|---|---|---|---|
| 1 | **Candidate-only** (cross-Cramér's V on `state.candidates`) | in atlas already | SHIPPED earlier |
| 2 | **Fish-trajectory** | `specs_todo/SPEC_distant_band_concordance_fish_trajectory.md` | **Slices 1+2 SHIPPED turn 130** (compute + UI). Slices 3+4 pending. |
| 3 | **L2-sweep** | `specs_todo/SPEC_l2_sweep_inheritance.md` | spec only |
| 4 | **Sliding-window** | `specs_todo/SPEC_sliding_window_inheritance.md` | spec only |
| 5 | **Cross-chromosome lineages** | `specs_new_turn131/SPEC_cross_chromosome_lineages.md` | spec only (new this turn) |
| +1 | **Unification dispatcher** | `specs_todo/SPEC_inheritance_unification.md` | **decision-gated** (don't build until 3/4 producers ship) |

So when you said "2 of 5 in slice" — that matches. Methods 1 and 2 are
the ones that have shipped code; methods 3, 4, 5 are spec-only; the
unification +1 is deferred until you've seen the producers.

---

## 3. What shipped this session (turn 130)

### Slice 1: Fish-trajectory lineage compute
Hungarian-aligned per-chromosome lineage clustering. State + helpers:
`state.lineageResult`, `state.lineageCacheKey`, `state._lineageComputeScheduled`,
`_hungarianChainProjection`, `_concordanceMatrix`, `_lineageClustering`,
`_lineageCacheKey`, `runLineageCompute`, `invalidateLineageCache`,
`_LINEAGE_DEFAULT_THRESHOLD`, `_LINEAGE_CHAIN_BREAK_AGREEMENT`. Spec:
`specs_todo/SPEC_distant_band_concordance_fish_trajectory.md`. Tests:
`tests/test_turn130_lineage_compute.js`.

### Slice 2: Lineage UI (color mode + strip)
`_lineageColor`, `_drawLineageStrip`, `setLinesLineageStripOn`,
`_LINES_LINEAGE_STRIP_KEY`, color-mode integration. Same spec.

### Slice 0: Auto-review dashed-row infrastructure
`_isAutoCandidate(cand)` predicate. CSS `.cand-list-item.is-auto`
(dashed border, opacity 0.85). 🤖 prefix in candidate list. Sort
modifier puts auto last. `_gatherActiveCandidatesForInheritance` skips
auto candidates. Spec:
`specs_todo/SPEC_review_surfaces_auto_and_lineages.md` Slice 0. Tests:
`tests/test_turn130_auto_review_infra.js` (32 tests).

---

## 4. What this session (turn 131) accomplished

**Spec-writing pass — no code shipped.** User instruction: *"stop coding,
finish writing all specs."*

19 new specs in `specs_new_turn131/`:

### Tier 1 (literature-anchored)
1. **`SPEC_MASTER_full_automatic_pipeline.md`** — umbrella vision.
   ~6-9 turns total to ship the click-button → all-inversions pipeline.
2. **`SPEC_recombinant_dosage_changepoint_detector.md`** — per-carrier
   dosage step. LANTA `STEP_C01h_recombinant_scanner.R` exists; this
   spec covers atlas-side viz + review.
3. **`SPEC_boundary_consensus_aggregator.md`** — KDE on stacked
   per-carrier observations from θπ-CUSUM + L3 contingency + dosage.
   Multi-mode flag when 2+ peaks separated >100 kb.
4. **`SPEC_metric_overlay_priors.md`** — 12+ secondary boundary
   detectors (Fst, dXY, Tajima's D, Q shift, GHSL, sim_mat insulation,
   λ₂, inv_likeness, DELLY/Manta breakpoints, depth, MI, indel slope).
5. **`SPEC_karyotype_per_interval_intersection.md`** — Quentin's 30-line
   post-processing. Splits GHSL karyotype runs by triangle/L2 intervals.
6. **`SPEC_hypothesis_test_framework_atlas.md`** — T1-T11 with BH.
   H1-H4 verdict synthesis. Per-candidate verdict panel + decision
   trail.
7. **`SPEC_multi_evidence_regime_framework.md`** — doctrine: K-means
   provisional, biological regimes need multi-layer agreement. K=6
   NOT biological by default.
8. **`SPEC_per_band_evidence_layer.md`** — per-band H, θπ, F_ROH,
   family. F_ROH-confounder, family-confounder, θπ-fail flags.
9. **`SPEC_inversion_age_origin_atlas.md`** — Porubsky 2022 GDS,
   Hartigan dip test, age via REF-INV gap.
10. **`SPEC_boundary_confirmation_5track.md`** — F4 panel A "convicting"
    composite (Hudson Fst + dXY + θπ-by-karyotype + repeat + |z|).
    4-of-5 majority verdict.
11. **`SPEC_marker_panel_design_atlas.md`** — private indel tier system,
    F1 hybrid controls.
12. **`SPEC_per_candidate_breeding_readiness_card.md`** — Atlas 5 Part B.
    Per-candidate one-pager with auto-pairing advice.
13. **`SPEC_gene_annotation_overlay.md`** — GFF3 integration.

### Tier 2 (UI / orchestration)
14. **`SPEC_cross_chromosome_lineages.md`** — Hungarian-align across
    chroms (also serves as Method 5 of the 4+1 inheritance system).
15. **`SPEC_multichrom_load_orchestrator.md`** — bulk-load 28 JSONs lazy.
16. **`SPEC_genome_wide_ideogram.md`** — 28-chrom F1 ideogram page.
17. **`SPEC_manuscript_bundle_export.md`** — .zip with README +
    catalogue + boilerplate + figures.
18. **`SPEC_lasso_inheritance_backgrounds.md`** — alpha intervals on
    per-sample lines (chat `f74cf5d4`).
19. **`SPEC_interval_collector.md`** — passive accumulator from
    tracked samples × contingency.

### Index
**`SPECS_TIER_INDEX.md`** — sorts all 47 specs (across both folders)
by tier with recommended build priority order.

---

## 5. Where to start the next chat

### Priority sequence (from `SPECS_TIER_INDEX.md`)

1. **L2-sweep auto-promote Slice 1** (Tier 1 producer; unlocks the
   dashed-row infra shipped in turn 130 follow-up).
   Spec: `specs_todo/SPEC_l2_sweep_inheritance.md` Slice 1.

2. **G-panel scaffold Slice 1** (Tier 2 review surface; enables
   auto/lineages tabs to render).
   **Pre-Slice fix**: `_rebuildCandidateRegistries()` (~0.3 turn) MUST
   land before Slice 1 — see `OBSERVATIONS_TO_FIX.txt`.
   Spec: `specs_todo/SPEC_g_panel_unified_groups.md` Slice 1.

3. **Trajectory matrix viewer Slice 3** (Tier 1 surface for shipped
   lineage compute).
   Spec: `specs_todo/SPEC_distant_band_concordance_fish_trajectory.md`
   Slice 3.

4. **Cross-chromosome lineages Slice 1** (Tier 2 producer, depends on (1)).
   Spec: `specs_new_turn131/SPEC_cross_chromosome_lineages.md`. **When
   you start this**, MOVE the spec from `specs_new_turn131/` to
   `specs_todo/`.

After this critical-path quartet, the rest of the master pipeline
(genome-wide ideogram, manuscript bundle export, breeding-readiness
cards) becomes possible.

### Honest scope
The "click-button → all-inversions" vision needs ~6-9 turns total.
None of the pieces are show-stoppers. Quentin reviews per-candidate
boilerplate prose, not a fully-written manuscript.

---

## 6. Open observations (`OBSERVATIONS_TO_FIX.txt`)

### Highest priority
- **Bug 1**: L3 arrow up/down promote/demote sub-panel highlighting +
  merge broken.
- **Pre-Slice fix for G-panel**: `_rebuildCandidateRegistries()` auto-sync.

### Lower priority
- Bug 3: page-1 PCA legend overlap (blocked on screenshot).
- Boundaries Mb mismatch (~14 vs 15.75).
- Top-bar hover gap, L3 detailed click no-op, diversity atlas plot
  placement, ST2 plot bottom margin, ROH heatmap S8b last-row collapse,
  Loess fit at zero coverage.

---

## 7. State / window slots through turn 130

```
state.lineageResult, state.lineageCacheKey, state._lineageComputeScheduled
state.l3HetColoring, state.linesLineageStripOn, state.__hetRateCache
state.dismissedAutoIds (planned)
state.chromSummary, state.activeChrom, state.bulkLoadProgress (planned)
```

```
window: _rebuildCandidateRegistries (planned), _hetRateColor,
_computeHetRateForL2, _invalidateHetRateCache, _getHetRateCache,
_refreshL3HetToggleAvailability, _hungarianChainProjection,
_concordanceMatrix, _lineageClustering, _lineageCacheKey,
runLineageCompute, invalidateLineageCache, _LINEAGE_DEFAULT_THRESHOLD,
_LINEAGE_CHAIN_BREAK_AGREEMENT, _lineageColor, _drawLineageStrip,
setLinesLineageStripOn, _LINES_LINEAGE_STRIP_KEY, _isAutoCandidate
```

---

## 8. Backups present

- `Inversion_atlas.html.bak_pre_lineage_compute` / `.bak_post_lineage_compute`
- `Inversion_atlas.html.bak_pre_lineage_ui` / `.bak_post_lineage_ui`
- `Inversion_atlas.html.bak_pre_auto_review_infra` / `.bak_post_auto_review_infra`
- Earlier: `.bak_post_het_coloring`, `.bak_pre/post_registry_bridge`,
  `.bak_post_turn129`

---

## 9. Bundle contents

`/mnt/user-data/outputs/Atlas_full_bundle_2026-05-05_turn131.tar.gz`

Contains:
- `Inversion_atlas.html` (current, 64,243 lines)
- `tests/` (all *.js, *.py)
- `specs_todo/` (27 .md — active queue)
- `specs_new_turn131/` (20 .md including tier index + README — pending review)
- `OBSERVATIONS_TO_FIX.txt`
- This handoff
- All other previous handoffs (kept for history)

`.bak_*` files NOT in bundle to keep size down.

---

## 10. Test count

**1900/0 PASS/FAIL** across all suites. No regressions in turn 131
(no code changed).

---

## 11. Honest framing for the next chat

The vision is real, the path is mapped, most producers exist on LANTA.
What's needed atlas-side:

1. **L2-sweep auto-promote** is the single biggest unlock. Without it
   the dashed-row review infrastructure has nothing to display.
2. **G-panel scaffold** is the review surface that turns auto-promoted
   candidates into reviewable items.
3. **Genome-wide ideogram + manuscript bundle export** are the
   user-facing deliverables that complete the click-button → manuscript
   promise.

Interpretive paragraphs of the manuscript — what inversions mean
biologically, what selection might be doing, what recurrent breakpoints
suggest — those are **always Quentin's job**. No pipeline writes those.
The atlas's role is to produce numbers Quentin trusts and present them
in a form that drops cleanly into a manuscript skeleton.

Walk the map carefully, respect cohort discipline, don't break the
test suite.
