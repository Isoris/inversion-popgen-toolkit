# HANDOFF — turn 156 — V-shape (u, agreement_fraction) per-candidate diagnostic

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (73,435 lines, +518 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.

**Picked up from**: post-turn-155 state (2776 / 0). Turn 154 from a
parallel session was already in place; turn 155 (mine) did the
threshold-in-cache-key honesty fix. Turn 156 is the next bounded
scientific deliverable — option B from the queue.

---

## 0. What this turn ships

The V-shape diagnostic plot — port of FIG_C30 from STEP29's coherence
pipeline (`current_followup/STEP29_candidate_coherence_and_polarity.R`).

For one focal candidate, plots per-sample `(u, agreement_fraction)`:
- **u** = rotated PC1 from `_getOrComputeUVRotation(cand.ref_l2)`
- **agreement_fraction** = from `computeStripeQuality(...)` rows. The
  fraction of informative dosage markers where the sample is closer to
  its assigned group's mean than to the average of the other groups'
  means.

A clean candidate produces a V-shape:
- Two high plateaus at u ≈ hom1_u and u ≈ hom2_u (agreement ~0.85–1.0)
- A dip at u ≈ het_u (agreement ~0.5–0.7) — HET samples sit between
  the two homo means by construction, so the "closer to own than to
  other" test is closer to a coin flip
- A noisy, smeared, or no-dip plot ⇒ the assignment is wrong or
  there's no real biological inversion structure

Pure read-only diagnostic. No mutation of state. Reuses every existing
helper (UV rotation, stripe-quality, dosage chunk).

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse, pragmatic.
Tarball is the standard handoff format.

---

## 2. What shipped this turn

### 2.1 Five new functions

```
_buildVShapeData(candidate, chunk, sqRows)
  → { ok, points: [{u, agreement, group, sample_idx, cga, coherence_class, stripe_quality}],
      het_u, hom1_u, hom2_u, n_total, n_with_data, reason? }
```

Joins UV rotation output (cohort-indexed `us` array) with stripe-quality
output (chunk-indexed rows). Builds a per-chunk-sample-id → cohort-idx
map, looks up `us[cohortIdx]` per row, drops samples with non-finite u
or agreement. Reasons enum: `MISSING_REF_L2`, `UV_FAIL`,
`MISSING_CHUNK`, `MISSING_SQ`, `NO_POINTS`.

```
_vShapeColor(group) → '#hex'
```

Maps `HOMO_1` / `HET` / `HOMO_2` to the existing `_BR_KARYO_COLORS`
palette. Falls back to literal hex strings when palette isn't loaded
(sandbox tests). Unknown groups render grey `#888888`.

```
_drawVShapePlot(canvas, data, opts) → { ok, n_drawn, bounds, hits, reason? }
```

Canvas renderer with axis labels, gridlines at agreement=0/0.25/0.5/0.75/1,
reference vertical lines at `het_u` (dashed, green), `hom1_u` (dashed,
blue), `hom2_u` (dashed, orange), reference horizontal at
`agreement=0.5`. Points drawn in `HOMO_1, HOMO_2, HET` order so HET
sits on top. Legend in top-right corner.

DPR-aware. Returns `hits[]` (each `{x, y, r, point}`) for future
tooltip wiring; v0 doesn't use them.

```
openVShapePlot(candidate)
closeVShapePlot()
```

Modal popup. Lazy-creates `#vShapePlotModal` on first call. Escape key
+ click-outside both close. Reads dosage state via
`_ensureDosageHmState()`; if no chunk loaded or `last_sq.candidate_id`
doesn't match the focal candidate, renders a helpful inline message
instead of crashing. Status footer reports per-group `n` and
`mean_agreement` plus `het_u`.

### 2.2 G-panel karyotype tab integration

Added `🔍 V-shape diagnostic` button next to the existing Color
PCA / Export TSV / Calibration buttons in the karyotype tab footer.
Click opens the popup at the current focal candidate.

Tooltip explains the prerequisite: "Requires that you ran 'Compute
stripe quality' on the dosage heatmap first."

### 2.3 No state mutation

Diagnostic only. No new state slots, no localStorage persistence, no
side effects on candidate data.

---

## 3. What did NOT change

- **`computeStripeQuality`** — untouched. Just consumed by the new builder.
- **`_getOrComputeUVRotation`** — untouched.
- **Dosage chunk pipeline** — untouched.
- **Inheritance system** — untouched.
- **G-panel state slots** — untouched.

---

## 4. Test status

|                          | LOC     | Tests           | Files |
|---                       |---      |---              |---    |
| Pre-session (turn 155)   | 72,917  | 2776            | 59    |
| Post-turn-156 (current)  | 73,435  | 2840            | 60    |
| **Δ session**            | +518    | +64             | +1    |

Full sweep: **2840 / 0**. JS syntax: clean. HTML parser: 0 errors.

`tests/test_turn156_vshape_diagnostic.js` (64 tests) covers:

1. All 5 functions declared and exposed on `window`
2. `_buildVShapeData` source structure (early-return reasons, point
   shape, output keys)
3. `_drawVShapePlot` source structure (DPR, gridlines, reference
   lines, draw order, legend, hits return)
4. `openVShapePlot` source structure (lazy DOM, Escape, click-outside,
   `_ensureDosageHmState` lookup, candidate-id match, graceful empty
   states, status line)
5. G-panel karyotype tab integration (button + click handler)
6. Sandboxed `_buildVShapeData` behaviour:
   - All 5 early-return reasons
   - Happy path with 6 samples → 6 points, V-shape signal verified
     (mean HOMO > mean HET in both directions)
   - Unknown sample id silently dropped
   - NaN u or NaN agreement silently dropped
7. Sandboxed `_vShapeColor` palette mapping
8. Existing flow preserved

No existing test needed inverting.

---

## 5. What Quentin can do next session

Workflow:

1. **Page 1:** lock + promote a candidate (existing workflow).
2. **Page 2:** open the dosage heatmap on that candidate.
3. **Heatmap toolbar:** click "Compute stripe quality". This populates
   `_ensureDosageHmState().last_sq = { candidate_id, rows }`.
4. **G-panel** (`g`): switch to the **karyotype** tab.
5. **Footer:** click `🔍 V-shape diagnostic`. Popup opens with the
   scatter.

What to look for:
- **Strong V-shape:** dip at `het_u` (green dashed line), high homo
  plateaus. Indicates clean candidate.
- **No dip / smeared:** noisy candidate, possibly wrong K, possibly
  no real inversion structure.
- **Mean agreement per group** in the status footer:
  - HOMO_1 mean ≥ 0.85, HOMO_2 mean ≥ 0.85, HET mean ≈ 0.55–0.70 →
    healthy
  - Any group below ~0.6 → suspect

The plot is also a useful "calibration" view — eyeballing several
candidates at once gives a sense of which are genuinely structural.

---

## 6. Things I almost broke and fixed

- **Initially tried to use `_pc1OfSample` from `_buildSampleLookups`.**
  That returns `fc.u` from `candidate.fish_calls`, which is the
  cached-rotation u from when the candidate was promoted. For the
  V-shape I want the *current* L2-level rotation — `cand.ref_l2`'s
  rotation specifically. Fixed by calling
  `_getOrComputeUVRotation(cand.ref_l2)` directly inside
  `_buildVShapeData` and indexing `us[cohortIdx]`. This means the
  V-shape uses the freshest rotation, not whatever was baked into
  `fish_calls` at promotion time.
- **Cohort vs chunk indexing.** `us` is cohort-indexed; `sqRows` are
  chunk-indexed (each row carries a sample id string from
  `chunk.samples[si]`). Fixed by building an in-function id→cohort map
  by walking `state.data.samples`, then looking up cohort idx per row
  via the row's sample id.
- **HET draw order.** Initially drew points in iteration order, which
  meant HOMO_2 (drawn last) covered HET dots in the dip region. Fixed
  by enforcing `HOMO_1, HOMO_2, HET` draw order so HET sits on top.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 156 patched, 73,435 lines.
- `tests/test_turn156_vshape_diagnostic.js` — NEW, 64 tests.
- `HANDOFF_2026-05-05_turn156_vshape_diagnostic.md` — this file.

Plus prior handoffs (carried for reference):
- `HANDOFF_2026-05-05_turn155_threshold_in_cache_key.md`
- `HANDOFF_2026-05-05_turn154_compute_ux_hardening.md`
- `HANDOFF_2026-05-05_turn153_inheritance_auto_register.md`
- `HANDOFF_2026-05-05_turn152_g_panel_inheritance_slice3.md`

---

## 8. Honest framing

**What's solid:**
- Pure read-only diagnostic. Zero risk of breaking anything else.
- Reuses three existing primitives — UV rotation cache,
  computeStripeQuality, dosage chunk state — without touching any of
  them.
- Graceful empty states at every prerequisite step (no candidate, no
  chunk, no stripe-quality, build fail) with helpful inline messages.
- Test discipline maintained: source-pattern + sandboxed unit tests,
  no DOM emulation.

**What's NOT done (and why that's right):**
- **No tooltip / hover detail.** `_drawVShapePlot` returns a `hits[]`
  array specifically so a future turn can wire mousemove → "this
  sample is CGA xxx, group HET, agreement 0.42, stripe_quality
  peripheral". Not in v0 because per-sample drill-in is most useful
  paired with the karyotype/heatmap surfaces, which already have it.
- **No batch / cross-candidate view.** A future turn could iterate
  saved candidates, render V-shapes side-by-side or as a sparkline
  matrix. Useful for triage. Not in v0 — single-candidate first to
  validate the signal looks right.
- **No polarity overlay.** STEP29 also computes L1/L2 polarity blocks.
  That's a separate diagnostic and would clutter v0.
- **No export.** Could ship a "Download PNG" button. Not necessary
  yet — Quentin can screenshot or use the canvas's built-in
  `toDataURL` if he needs.

**What's queued:**
- **A** — UV refactor on locked groups. Locked-labels caveat still
  applies.
- **C** — Per-sample-lines het coloring. Perf concerns at 226 × 30k.
- **D** — Pivot based on what het / inheritance / V-shape reveals.
- **F** — Per-group "show fish" expand toggle in inheritance card.
- **G** — V-shape tooltip / hover (uses `hits[]` from this turn).
- **H** — Cross-candidate V-shape gallery / sparkline matrix.

End of handoff.
