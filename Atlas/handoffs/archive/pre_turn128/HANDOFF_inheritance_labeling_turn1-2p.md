# Atlas inheritance/labeling work — handoff

**Status as of end of session.** 16 turns shipped (turn 1, 2a–2p).
54,202 lines of atlas HTML. 305 tests across 16 suites, all green.

This document captures what's done, what's pending, and the technical
specs for each pending item so future sessions (or another developer)
can pick up cleanly.

---

## What's shipped (16 turns)

| Turn | Feature | Tests |
|------|---------|-------|
| 1    | Band-reach bimodality detector | 12 |
| 2a   | Cross-candidate Cramér's V matrix (compute) | 33 |
| 2b   | Inheritance-group clustering (Jaccard + average linkage) | 19 |
| 2c   | UI labels strip ("I1·3g" pills above PC1 panel) | 20 |
| 2d   | Per-band haplotype annotations panel (page 2) | 23 |
| 2e+2g| Atlas export contract v1 → v2 (full atlas context dump) | 32 |
| 2f   | Auto-classifier + vocab system (binary/standard/multi2/multi3/free) | 18 |
| 2h   | Cramér's V matrix UI (modal with cell drill-down) | 16 |
| 2i   | Tracked-linkage projection (genome-wide band shading) | 18 |
| 2j   | Page-2 click-to-shade + family-confound warning | 11 |
| 2k   | Annotation cockpit page (page 21, full-width cursor canvas) | 24 |
| 2l   | Cockpit labeling tools (digit keys + linkage + hap panel) | 15 |
| 2m   | Atlas-tools toolbar (export/matrix/auto-fill/labels-toggle) | 12 |
| 2n   | Cache invalidation hardening (locked_labels hash + bp boundaries) | 16 |
| 2o   | Two-mode architecture (atlasServer adapter + reference Python server) | 23 |
| 2p   | Inheritance pill tooltip (hover for full group composition) | 13 |

---

## Still on the deferred list (NOT shipped)

### D1 · Cross-chromosome cohort aggregation (~3–5 turns, large)

**Problem.** All cross-candidate analyses today (Cramér's V matrix,
inheritance clustering, linkage projection) operate within a single
chromosome's loaded candidates. A real cohort spans 28 chromosomes;
nothing in atlas currently lets you ask "is candidate I3 on LG28
co-inherited with candidate I7 on LG14?"

**Why deferred.** This is genuinely a different scope.
Single-chromosome state lives in `state.candidates`, which is replaced
on every chromosome load. Cross-chrom would need either:

- (a) Persistent cohort-wide state (`state.cohort_candidates`,
  populated by walking all chromosomes' JSONs at startup), or
- (b) A server-mode-only feature that asks the backend for the
  aggregate.

**Recommended path: (b), gated on server mode.**

Reasons:

1. The browser already strains under a single-chromosome cohort with
   ~40 candidates. 28 × 40 = 1120 candidates, plus the cross-candidate
   matrix is O(N²) = ~1.25M cells. Browser-side computation is risky.
2. Server mode (turn 2o) is built for exactly this kind of heavy
   compute. Add `/compute/cross_chrom_matrix` to `atlas_server.py`,
   call it from atlas, render the result.
3. The atlas-side change is small: a new tab "11c cohort" with a
   matrix viewer that calls `atlasServer.compute('cross_chrom_matrix')`
   and renders the JSON response. Empty state when not in server mode.

**Spec:**

- New page21b (or page22): "cohort matrix"
- Empty state shown unless `atlasServer.mode === 'server'`
- Calls `atlasServer.compute('cross_chrom_matrix', { metric: 'cramer_v' })`
- Server returns `{n_items, ids, chroms, matrix: [n×n flat array]}`
- Render with same logic as turn 2h's `renderInheritanceMatrix` but
  with chromosome-grouped axis labels (e.g. "LG01 ▸ LG28")
- Click cell → drill into per-pair contingency (server-side compute)

**Server-side implementation hints (Python):**

```python
def _compute_cross_chrom_matrix(args, project_root):
    metric = args.get('metric', 'cramer_v')
    # Walk project_root for all candidates_*.json files
    # Build per-candidate label arrays (sample ID alignment is critical)
    # Compute pairwise Cramér's V using vectorized numpy
    # Return {ids, chroms, matrix: numpy_arr.flatten().tolist()}
```

**Critical requirement for the cohort feature:** sample identity
alignment. Different chromosomes' JSONs may list samples in different
orders. The server compute MUST align by sample ID before computing
contingency. Atlas's existing per-chrom code doesn't worry about this
because all candidates on one chrom share the same sample ordering.

---

### D2 · `chrom_len_bp` robustness (~1 turn, small)

**Problem.** The cockpit's chromosome extent (`_annoCockpitChromExtent`)
falls back to candidate ranges when `state.data.chrom_len_bp` is
missing. This means if no candidates have been promoted yet, the
cockpit shows "no candidates" even when other parts of atlas (windows,
tracks) know perfectly well that the chromosome is N Mb long.

**Spec for fix:**

1. Read order in `_annoCockpitChromExtent`:
   - First: `state.data.chrom_len_bp` (most authoritative)
   - Second: `state.data.windows[last].end_bp` (already implemented)
   - Third: `state.data.tracks[*].end_bp` max (currently missing —
     atlas often loads tracks even before candidates exist)
   - Fourth: candidate range (current fallback)
2. Always show the chromosome axis (Mb tick marks) once any extent is
   known, even with 0 candidates.
3. Update the page21 empty-state copy to clarify: "promote candidates
   on page 1 to begin annotating; the chromosome extent is already
   known."

**Test plan:** add 3 cases to `test_annotation_cockpit.js`:
- chrom_len_bp present + 0 candidates → extent returns chrom range
- tracks present + 0 candidates → extent returns tracks range
- nothing loaded → returns null (current behavior)

---

### D3 · Cockpit live-sync to other pages (~1 turn, small)

**Problem.** Currently, if you edit a candidate on page 2 (change K,
re-cluster, edit haplotype labels) and then switch back to page 21
(annotation cockpit), the cockpit re-renders correctly via the tab
activation hook. But if both pages are visible in different tabs (or
via a hypothetical split-pane future feature), the cockpit doesn't
update until tab activation fires.

**Spec for fix:**

1. Wire a state-change observer pattern. Pick a single source of truth
   event: `_dispatchAtlasEvent('candidatesChanged')`.
2. Call it from every mutation site:
   - K-means re-cluster (page 1 lock-colors)
   - K change (page 1 K-picker)
   - Haplotype label edit (turn 2d wire function)
   - Boundary refinement (page 4 catalogue)
   - Add/remove candidate (page 1 promote, page 4 promote)
3. Cockpit subscribes to this event and calls
   `refreshAnnotationCockpit()` if its page is currently active.
4. Same pattern for the matrix modal — close on `candidatesChanged`
   to avoid showing stale matrices (or recompute if visible).

**Why I didn't ship this.** Wiring state-change events across all
mutation sites is high-risk for regressions in code I haven't audited.
The current tab-activation refresh pattern works for the
single-screen-at-a-time UI. This becomes important if/when atlas grows
multi-pane support.

---

### D4 · Real callers for `atlasServer` (~2–4 turns, medium)

**Problem.** Turn 2o shipped the wire but didn't reroute any existing
callers. Today, even when atlasServer is available:
- `downloadAtlasExport()` still triggers a browser download instead of
  writing to the project folder
- The Cramér's V matrix still computes in JS (slow on large cohorts)
- Haplotype labels still persist to localStorage instead of the
  project folder (lost when localStorage clears)

**Spec for rerouting (incremental):**

Phase A — Export rewiring:
- `downloadAtlasExport(opts)` checks `atlasServer.mode === 'server'`
- If yes: `atlasServer.write('exports/atlas_<chrom>_<timestamp>.json', payload)`
  and show a "saved to project folder" toast instead of triggering download
- If no: existing browser-download behavior

Phase B — Haplotype labels:
- `saveHaplotypeLabels(c)` writes to localStorage AND (if server mode)
  to `metadata/haplotype_labels.json` on disk
- `loadHaplotypeLabels(c)` checks the file first if available, falls
  back to localStorage
- Handles offline editing → online sync gracefully (last-write-wins
  on the file path is fine for single-user)

Phase C — Heavy compute:
- Cramér's V matrix: if server mode, call
  `atlasServer.compute('contingency_matrix', {candidates: [ids]})`
- Server returns {n, matrix: flat_array} — atlas renders unchanged
- Speedup: numpy vectorized contingency is ~50–100× faster than the
  JS implementation for 200×200 matrices

**Why I didn't ship this.** Real-callers rerouting is a follow-up
that benefits from having the actual server running with project data
loaded. It's a much better experience to write+test the full cycle
once you have the Python server doing the real compute, instead of
shipping wire-only stubs that nobody can verify against.

---

### D5 · Multi-chromosome state preservation (small to medium)

**Problem.** Switching chromosomes resets `state.candidates`. This
means if you label candidates on LG28, switch to LG14, then come back
to LG28, you have to reload the JSON. Annotations persist (localStorage
keyed by chrom+id), but the candidate roster does not.

**Spec for fix:**

- Add `state.candidatesByChrom: { [chrom]: { [id]: candidate } }`
- On chromosome load, populate this map
- Switch chromosome → save current `state.candidates` into
  `state.candidatesByChrom[currentChrom]`, restore from
  `state.candidatesByChrom[newChrom]` if cached
- Force-reload button to clear and re-fetch from JSON
- Memory budget: ~226 fish × ~40 candidates × 1 byte per locked label
  + metadata = ~10 KB per chromosome × 28 chroms = 280 KB. Negligible.

**Why I didn't ship this.** Tightly coupled to D4 (server-mode
architecture is the right home for full-cohort persistence — let the
filesystem be the source of truth, not the browser). Worth doing
together.

---

### D6 · Mobile / narrow-viewport polish (small)

**Problem.** The toolbar buttons added in turn 2m + the server badge
in turn 2o + the existing 16 tab buttons may overflow on narrow
viewports. Atlas already has aggressive flex-wrap rules for tabs but
the tools group is `margin-left:auto` and may push tabs into ellipsis.

**Spec for fix:**

- At viewport <1100px: collapse tools group into a single overflow
  menu button (⋮) that opens a dropdown
- Server badge stays visible in collapsed mode (it's a status
  indicator, not an action)

---

## Architectural notes for next dev

### Source of truth for state

Atlas has TWO state objects in some environments:
- Lexical `const state = ...` at the top of the script
- `window.state` set by `_nrEnsureState` to `{}` if undefined

In production browsers these can be made the same object via explicit
`window.state = state` assignment, but atlas does NOT do this — so
in JSDOM (and possibly some browsers) they differ. Functions like
`setLinesInheritanceLabelsOn` use `window.state` if truthy else
lexical, which is an ambiguous source of bugs.

**Recommendation:** Add `window.state = state;` near the top of the
script (after lexical state is fully initialized). One line, removes
all the two-state hedging code I had to write in turn 2m.

### Test infrastructure

Tests use JSDOM. `/tmp/atlas.html` is the working copy; `cp` from
`/home/claude/work/atlas.html` before running tests. The 16 test files
are `/tmp/test_*.js` — same names as the ones shipped in `/mnt/user-data/outputs/`.

To run all suites:
```bash
cd /tmp
for t in test_*.js; do NODE_PATH=/tmp/node_modules node $t; done
```

JSDOM canvas has no backend — tests use `makeStubCtx()` helpers to
avoid `getContext` returning null. Pattern is in every test file.

### Atlas size

Started at 49,558 lines (turn 1).
End of session: 54,202 lines (+4,644 lines, +9.4%).

This is a lot of growth for one session but the per-feature size is
reasonable (~290 lines/turn). The bulk of new code is in the
inheritance subsystem (turns 2a–2c), the cockpit (2k–2l), the export
contract (2g), and the server adapter (2o).

### Skill files / conventions

- All compute helpers prefer `window.state` if truthy, else lexical
- `_gatherActiveCandidatesForInheritance()` is the canonical "items
  list" projection — never iterate `state.candidates` directly in
  inheritance code
- Cache key for inheritance: `_inheritanceCacheKey(items, mode)` which
  hashes locked_labels + bp boundaries. This is the only correct way
  to detect "did the candidates change".
- Functions intended for testing get `window.X` exposes near their
  definition. Don't put them at the end of the script — JSDOM aborts
  on canvas init failures and never reaches the bottom.

---

## Files in /mnt/user-data/outputs/

| File | Purpose |
|------|---------|
| `Inversion_atlas.html` | Main atlas (54,202 lines) |
| `atlas_server.py` | Reference Python server (259 lines, stdlib only) |
| `HANDOFF.md` | This file |
| `test_band_reach_bimodality.js` | Turn 1 |
| `test_cross_candidate.js` | Turn 2a |
| `test_inheritance_clustering.js` | Turn 2b |
| `test_inheritance_ui.js` | Turn 2c |
| `test_haplotype_annotations.js` | Turn 2d |
| `test_atlas_export.js` | Turn 2e+2g |
| `test_auto_classifier.js` | Turn 2f |
| `test_inheritance_matrix.js` | Turn 2h |
| `test_tracked_linkage.js` | Turn 2i |
| `test_band_click_linkage.js` | Turn 2j |
| `test_annotation_cockpit.js` | Turn 2k |
| `test_cockpit_labeling.js` | Turn 2l |
| `test_atlas_tools.js` | Turn 2m |
| `test_cache_invalidation.js` | Turn 2n |
| `test_atlas_server.js` | Turn 2o |
| `test_pill_tooltip.js` | Turn 2p |

To run server-mode end-to-end:

```bash
cd /path/to/your/project_folder
python3 atlas_server.py .
# In another terminal: open Inversion_atlas.html in browser
# Toolbar should show "● server" in green
```

---

## Workflow as it stands today

The intended labeling workflow is now end-to-end functional:

1. Open atlas → load chromosome JSON
2. Promote candidates on page 1 (lock colors → promote)
3. Click tab "11b annotation" → cockpit opens
4. Use ←/→ to walk cursor across chromosome; pills above the lines
   panel show inheritance group counts; hover any pill for full
   composition tooltip (turn 2p)
5. When cursor enters a candidate, the right footer pane shows that
   candidate's haplotype-label panel with auto-classifier suggestions
6. Press 0/1/2/... to project that band's fish across the chromosome
   → genome-wide linkage shading appears on the canvas, linkage table
   appears in the left footer pane
7. Pick haplotype labels via dropdown (auto-fill button does cohort-
   wide if you want)
8. Click 📥 export in the toolbar → full atlas context as JSON
9. Click 🔗 matrix in the toolbar → Cramér's V matrix opens
10. Click 🪄 auto-fill → cohort-wide classifier run

Caching invalidates correctly across all paths (turn 2n). The
toolbar buttons handle missing functions gracefully. Server mode
detects automatically and shows status badge.

This is what was asked for at the start of the session. The deferred
list (D1–D6) is genuine future work, not regrets — each item has its
own cost/benefit and several are gated on data flows that don't yet
exist (server compute payloads, multi-chrom workflow).

Good luck with the Theta pi testing this evening.
