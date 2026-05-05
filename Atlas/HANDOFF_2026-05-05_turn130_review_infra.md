# HANDOFF — turn 130 (auto-review infra follow-up)

**Date**: 2026-05-05
**Atlas line count**: 64,243 (was 64,160 → +83)
**Test status**: **1900 PASS / 0 FAIL** across all suites (was 1868 → +32)

---

## What shipped (Slice 0 of review-surfaces)

Quentin asked for review surfaces for auto-promoted candidates + lineage
blocs — either a G-panel tab, or last rows in the catalogue with dashed
outlines. After mapping the request to the codebase, the catalogue
(page 3) is L2-row-based, not candidate-based, so it doesn't naturally
host candidate-list rows. The right surface is the **page-2 candidate
list** (`#candListContainer`) — which IS candidate-list-based.

Slice 0 ships the **plumbing infrastructure** so when L2-sweep lands its
auto-promote mechanic next turn, every review surface lights up
automatically. Nothing is user-visible YET (no producer exists), but
the next L2-sweep producer hooks straight in.

### Code

- **`_isAutoCandidate(cand)`** — central predicate. Returns true when
  `cand.source` starts with `'auto_'` AND `!cand.confirmed`. Source tags
  reserved: `'auto_l2_sweep'`, `'auto_lineage'`. Confirmed candidates
  (regardless of original source) return false.

- **`_gatherActiveCandidatesForInheritance` skips auto candidates** —
  unreviewed auto-promoted candidates don't pollute the inheritance
  Jaccard compute. I·g pills + cross-candidate matrix only see
  user-confirmed structure. Auto candidates show up only after the user
  confirms (`cand.confirmed = true`), at which point they re-enter the
  pipeline.

- **Defensive fallback** — when `_isAutoCandidate` isn't in scope (e.g.
  the turn-129 registry-bridge test sandboxes only the gather function),
  the inline source check still skips auto candidates. The older test
  passes unchanged.

- **CSS treatment** — `.cand-list-item.is-auto` gets dashed borders +
  opacity 0.85; `.cli-auto-prefix` styles the 🤖 emoji.

- **Sort modifier in `refreshCandidateListUI`** — auto candidates land
  at the bottom; non-auto sorted by `created_at` desc as before. Within
  each group, ordering preserved.

- **Render modifier** — auto rows get the `is-auto` class +
  `<span class="cli-auto-prefix">🤖&nbsp;</span>` before the label.

### Tests

`tests/test_turn130_auto_review_infra.js` — 32 tests:
- 11 source-level (predicate definition + window export, gather skip,
  defensive fallback, sort modifier, class emission, prefix emission,
  CSS rules)
- 21 behavioural sandboxed (predicate edge cases for null/empty/manual/
  auto/confirmed; gather skips unconfirmed auto but keeps confirmed
  auto; defensive fallback works without `_isAutoCandidate` in scope;
  sort comparator semantics — non-auto first by created_at desc, auto
  last by created_at desc)

### Spec status

`SPEC_review_surfaces_auto_and_lineages.md` — Slice 0 marked SHIPPED.
Slices 1–3 wait on L2-sweep auto-promote producer.

---

## What this means in practice

When `SPEC_l2_sweep_inheritance.md` ships its Slice 1 (auto-promote from
sweep) next turn, it just needs to:

1. Call something like `state.candidateList.push({ ..., source: 'auto_l2_sweep', confirmed: false })`
2. Call `persistCandidateList()` + `refreshCandidateListUI()`

Then immediately:
- The auto candidate appears as a dashed row at the bottom of the
  page-2 saved-candidate list with a 🤖 prefix.
- I·g pills + cross-candidate matrix continue to show only user-confirmed
  structure.
- When the user clicks "confirm" (any path that sets
  `cand.confirmed = true`), the dashed treatment falls away and the
  candidate re-enters the inheritance compute.

No additional UI work needed for Slice 1 of L2-sweep to be useful.

---

## What still doesn't exist

- **Candidate strip on page 1** — that's canvas-rendered via
  `drawCandidateBar`. Auto candidates render in solid form there for now;
  per-track dashed canvas treatment is a separate ~0.5 turn of canvas
  rendering work. Low priority — page-2 list is the primary review
  surface.
- **G-panel `auto` tab** — the popup review queue. Slice 2 of this spec.
  Bulk confirm/dismiss + Inspect shortcut. Estimated 0.5 turn after
  G-panel scaffold ships (`SPEC_g_panel_unified_groups.md` Slice 1).
- **G-panel `lineages` tab** — same idea for lineage blocs. Slice 3.
- **Cross-candidate matrix filter** — currently the matrix builds from
  whatever `inheritanceGroupClustering` outputs. Since gather already
  skips auto, the matrix is implicitly auto-free; explicit filter check
  could be added defensively but isn't strictly needed.
- **Lines-panel candidate-bands** — that spec is still blocked on
  Quentin's 1/3-vertical-vs-full-height design question.

---

## File deltas

```
Inversion_atlas.html              64160 → 64243   (+83 lines)

NEW tests:
  tests/test_turn130_auto_review_infra.js                32

UPDATED specs:
  specs_todo/SPEC_review_surfaces_auto_and_lineages.md
    Slice 0 marked SHIPPED in status header + slice list

Backups (drop before bundling):
  Inversion_atlas.html.bak_pre_auto_review_infra
  Inversion_atlas.html.bak_post_auto_review_infra
```

---

## Recommended next-turn priorities

1. **L2-sweep Slice 1** — actually produce auto-promoted candidates so
   the infrastructure shipped here has data to work on. Spec
   `SPEC_l2_sweep_inheritance.md`. After this lands, the dashed-outline
   review surface on page 2 lights up immediately.

2. **G-panel scaffold (Slice 1 of `SPEC_g_panel_unified_groups.md`)** —
   prerequisite for the `auto` and `lineages` review tabs. Currently the
   bridge from turn 129 unblocks the popup scaffold.

3. **Slice 3 of fish-trajectory** — co-membership matrix viewer modal.
   Standalone visual investigation tool for the 226×226 concordance
   matrix.

4. **Diagnostic protocol** — Quentin opens the atlas with the lineage
   strip and explores LG28; reports which thresholds need calibration.
