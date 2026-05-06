# HANDOFF — turn 158 — popstats QC/popstats split + off-by-default

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (73,773 lines, +273 LOC)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.

**Closes Quentin's items #5 + #6** from the turn 157 handoff §8:

> "We separate the tracks on the right hand side selection for QC and
> for popstats. And turn off popstats track by default but QC tracks
> are off by default."

(Read literally: BOTH QC and popstats are off by default. Confirmed
by the design intent.)

**Picked up from**: post-turn-157 working tree, which now also contains
a parallel-session turn 157 ("turn 157b": V-shape tooltip — option G
from my prior queue, wired the `hits[]` array I returned in turn 156).
That track is clean and orthogonal to this one. **Treating my
dosage-bridge work as "turn 157a" and the V-shape tooltip as "turn
157b"** in the file history; turn 158 stacks on top of both.

---

## 0. What this turn ships

Three coordinated changes to the popstats page (page 8 in the UI):

1. **Track categorization.** Each entry in `POPSTATS_TRACKS` now carries
   a `category` field with one of three values:
   - `always` — structural / discovery anchors (ideogram, sim_collapse,
     z). Always visible when data is present.
   - `qc` — data-quality diagnostics (snp_density, beagle_unc,
     coverage, low_cov_count). Off by default.
   - `popstats` — population-genetic statistics (theta_invgt,
     fst_hom1_hom2, hobs_hexp, delta12, delta12_multi). Off by
     default.

   Auto-discovered tracks (from `state.data.tracks` and from
   `popgenGallery.discoverTracks`) keep their legacy "default-on
   when data present" behaviour — they're treated as `category: 'other'`.

2. **Storage schema extended.** `loadPopstatsView` previously returned
   a bare `Set` of explicitly-hidden track ids. Now returns
   `{ hidden: Set, shown: Set }`. The `shown` set is new — it stores
   explicit user opt-ins for QC/popstats categories. Storage JSON now
   carries both `hiddenChips` and `shownChips`. Back-compat:

   - Legacy storage with only `hiddenChips` loads cleanly (shown stays
     empty — meaning QC/popstats stay off, which is the intended
     default).
   - `savePopstatsView` accepts either the new `{hidden, shown}`
     object or a bare `Set` (treated as the hidden set with an empty
     shown set).
   - Corrupted JSON returns empty sets without throwing.

3. **Chip strip grouped + visibility per-category.** The right-hand
   selector now sorts chips by category (always → qc → popstats →
   other) and inserts small section labels (`QC:` and `popstats:`)
   between groups. Each chip carries a `data-category` attribute so
   future styling can colorize the sections. The visibility check is
   now:

   ```
   alwaysOn        → visible
   category=always → visible iff hasData
   category=qc     → visible iff in view.shown
   category=popstats → visible iff in view.shown
   other (legacy)  → visible iff hasData (default-on preserved)
   ```

   Click handler dispatches based on category: qc/popstats toggle
   `cur.shown`; always/other toggle `cur.hidden`. An explicit show
   also clears any prior hide of the same id, so the storage stays
   coherent across rapid toggling.

   Tooltips for qc/popstats chips include "(QC track — off by default;
   click to show.)" / "(popstats track — off by default; click to
   show.)" so users understand why the chip is dim on first visit.

---

## 1. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
3. **C. macrocephalus wild** — future paper.

Quentin Andres (Kasetsart University Bangkok). Direct, terse,
pragmatic. Tarball is the standard handoff format.

---

## 2. Specific category assignments

| Track id          | Category   | Rationale                              |
|---                |---         |---                                     |
| ideogram          | `always`   | Structural; chromosome cartoon ref.    |
| sim_collapse      | `always`   | Local-PCA similarity matrix; primary discovery. |
| z                 | `always`   | Robust \|Z\|; primary discovery signal. |
| snp_density       | `qc`       | Variant density — QC for coverage gaps. |
| beagle_unc        | `qc`       | Imputation confidence — QC artifact.   |
| coverage          | `qc`       | Read depth — QC for callable regions.  |
| low_cov_count     | `qc`       | Per-window low-cov sample count.       |
| theta_invgt       | `popstats` | π by inversion genotype.               |
| fst_hom1_hom2     | `popstats` | Differentiation between homo bands.    |
| hobs_hexp         | `popstats` | Heterozygosity ratio per group.        |
| delta12           | `popstats` | Top-1 vs top-2 ancestry confidence.    |
| delta12_multi     | `popstats` | Multi-scale ancestry Δ12.              |

**Note on `z`**: it's technically a population-genetic statistic
(robust outlier signal from local PCA), but it's the primary discovery
signal that everyone needs to see when opening page 8. Categorized as
`always` deliberately. If Quentin wants `z` off by default too, flip
its category to `popstats` in one line.

---

## 3. What did NOT change

- **`collectPopstatsTracks`** — untouched. Still returns the union of
  built-in tracks + auto-discovered + popstatsLive-discovered, with
  the gallery filter applied. Auto-discovered tracks keep their legacy
  default-on behaviour because they don't have a `category` field
  (and `_popstatsCategoryOf` returns `'other'` for them, which falls
  through to the legacy visibility branch).
- **The actual track renderers** (`drawPopstatsTracks`, etc.) —
  untouched. They iterate whatever `tracks` array has, with the
  visibility filter applied externally.
- **The `popgenPage6.wrapAllTrackDefs` wiring** at the end of
  `collectPopstatsTracks` — untouched. The category info is preserved
  through the wrap because the wrapper does an object spread.
- **All non-popstats pages** — untouched.
- **All turn 152–157 work** (G-panel, inheritance, V-shape, tooltip,
  dosage bridge) — preserved.

---

## 4. Test status

|                                | LOC     | Tests           | Files |
|---                             |---      |---              |---    |
| Pre-session (turn 157, post-b) | 73,500  | 2926 (incl. 157b) | 62  |
| Pre-turn-158 (this turn start) | 73,500  | 2926            | 62    |
| Post-turn-158 (current)        | 73,773  | 2991            | 63    |
| **Δ this turn**                | +273    | +65             | +1    |

Full sweep: **2991 / 0**. JS syntax: clean. HTML parser: 0 errors.

`tests/test_turn158_popstats_qc_split.js` (65 tests) covers:

1. Track categorization: every entry has a category; specific
   id→category mappings for all 12 entries; total entry count
   preserved (12).
2. `_popstatsCategoryOf` helper: declared, exposed, falls back to
   `'other'` for missing/non-string category.
3. Storage shape: `loadPopstatsView` returns `{ hidden, shown }`,
   reads both arrays from JSON; `savePopstatsView` accepts both
   shapes, persists both arrays.
4. `renderPopstatsPage` source structure: visibility branches,
   category sort order, section labels, chip click dispatch by
   category, off-by-default tooltip hint.
5. Sandboxed storage round-trip: empty, full, back-compat (bare Set,
   legacy hiddenChips-only), corrupted JSON.
6. Sandboxed `_popstatsCategoryOf` behaviour: all 6 cases.
7. Existing flow preserved: render function, exports,
   `collectPopstatsTracks`, alwaysOn flags, V-shape from turn 156,
   dosage bridge auto-install from turn 157a.

No existing test needed inverting — the schema change is additive
(legacy storage loads cleanly, legacy `Set` arg to `savePopstatsView`
still accepted).

---

## 5. What Quentin can verify after this turn

1. **Reload Inversion_atlas.html.** Open page 8 (popstats).
2. **Expected first-visit state**: chip strip shows only ideogram +
   sim_mat + Z as active (or whichever subset has data). The QC and
   popstats sections show small "QC:" / "popstats:" labels followed
   by dim chips. Hover any dim chip — tooltip ends with "(QC track —
   off by default; click to show.)" or the popstats variant.
3. **Click a QC chip** (e.g. SNP density) → it activates, the track
   appears below. **Reload the page** → state persists; the chip
   stays active.
4. **Click an active chip a second time** → it deactivates, track
   disappears, state persists.
5. **Existing user with prior storage** (where they had explicitly
   hidden some tracks): their hidden tracks remain hidden. Their
   previously-default-on tracks (QC + popstats) **become hidden** on
   first reload after this update — they need to re-opt-in via clicks.
   This is the intentional behaviour change.

If anything in step 2 is wrong (e.g. a popstats chip shows as active
on first visit despite no `shown` storage), check the browser console
and paste back. The visibility check has clear branching, so any bug
should be reproducible.

---

## 6. Things I deliberately did NOT do

- **Did not add a "show all QC" / "hide all QC" macro button.** Quentin
  can ask if it'd help. With 4 QC chips and 5 popstats chips, individual
  toggling is fast.
- **Did not touch the popstats live-server path.** The `popstatsLive`
  state (auto-discovered tracks from the request layer) keeps its
  legacy default-on behaviour because those tracks don't have a
  `category` field. If Quentin wants those off by default too, the
  fix is to assign `category: 'popstats'` inside the
  `popgenGallery.discoverTracks` adapter or in the wrap step. That's
  a turn for the live-server work (#3 from the queue).
- **Did not change `z`'s category.** It's primary discovery; defaults
  to `always`. Flip it to `popstats` if Quentin disagrees.
- **Did not move the chip strip's CSS layout.** It still uses
  `flex-wrap: wrap`. The section label uses a `.ps-chip-section`
  class — Quentin can style it from the existing CSS file.

---

## 7. Files in the bundle

- `Inversion_atlas.html` — turn 158 patched, 73,773 lines.
- `tests/test_turn158_popstats_qc_split.js` — NEW, 65 tests.
- `HANDOFF_2026-05-05_turn158_popstats_qc_split.md` — this file.

Plus prior handoffs (carried for reference):
- 157b (V-shape tooltip — parallel session)
- 157a (dosage bridge load — me)
- 156 (V-shape diagnostic)
- 155, 154, 153, 152

---

## 8. Honest framing

**What's solid:**
- Bounded, additive, fully back-compat. Legacy storage loads cleanly;
  legacy callers of `savePopstatsView(hiddenSet)` still work.
- Visibility logic is one function with clear category branches. Tests
  cover all 5 paths.
- Chip strip groups visually with section labels.

**What's NOT done (and why that's right):**
- **#3 from queue (server-connected tracks not computed)** —
  partially addressed by turn 157a (dosage bridge load), but the
  "tracks not computed when popstats live server is connected"
  symptom wasn't reproducible without a live server run. Quentin
  should verify after this turn lands. If the issue persists, console
  output from the now-installed bridge will pin the cause.
- **#4 (G-bar groups → popstats live recompute)** — biggest item from
  the queue, ~200-400 LOC, separate turn.
- **#7 (LD common-vs-rare double heatmap)** — separate turn.
- **#8 (plot export)** — separate turn.
- **A "z is also off by default" toggle** — held back. Trivial change
  if Quentin asks.

**What's queued:**
- **Verify defaults work as expected** (Quentin's hands-on next session).
- **#3 / #4** — popstats live server consistency. After Quentin
  confirms whether #1 from turn 157a fixed the symptom.
- **#7** — LD double-heatmap.
- **#8** — plot export.

End of handoff.
