# Handoff — TEfull layer atlas integration

**Status**: ⚪ data shipped · atlas changes pending
**Files on disk**:
- `Cgar_TE_density/C_gar_LG{NN}_repeat_density_TEfull.json` — 28 JSONs, ~6-8 MB each
- `Cmac_TE_density/C_mac_LG{NN}_repeat_density_TEfull.json` — 27 JSONs
- `RUN_REPORT_TE_FULL.txt`
**Triggers to start**: when the parallel Diversity Atlas chat finishes and merges, or
   when you have a free session and want to polish the boundaries-page UX.

---

## What works TODAY without any atlas changes

You can validate and use the TEfull layer right now by drag-drop. The
existing v2 loader accepts the new JSONs because they share the schema
(version=2, binning_source=scrubber_windows, by_class={...}). The atlas
discovers all classes from the JSON's `by_class` keys.

Specifically, all of these work without modification:

- Drag-drop one or many TEfull JSONs onto the boundaries page
- The class dropdown lists all ~85 classes alphabetically
- The `__young` / `__old` suffixes keep each class's three variants
  adjacent in the sort order (e.g. `Gypsy_LTR_retrotransposon`,
  `Gypsy_LTR_retrotransposon__old`, `Gypsy_LTR_retrotransposon__young`)
- Class summary auto-scan computes q99 and Δ-vs-chrom for every class,
  flags enriched ones, surfaces them at the top of the boundaries panel
- Tight LOESS bands (0.05 / 0.10) work on every class
- ↑/↓ keys cycle through classes in panel
- `default_class: "all_TE"` in the JSON makes the atlas land on the
  aggregate track first, which is the most useful starting view

This means **the layer is usable for manuscript-quality investigation
before any atlas code changes**.

---

## What has friction (the actual handoff work)

Five UX issues, in priority order. Each is a separable atlas change.

### 1. Class dropdown overflow (HIGH priority)

Going from 26 to 85 classes in a single dropdown is a usability cliff.
80 alphabetically-sorted entries scroll past the viewport on most screens.

**Suggested fix**: split into 3-tier dropdown:

- Top tier: aggregates (`all_TE`, `young_TE_all`, `old_TE_all`,
  `insertion_count`, `intact_element_count`)
- Middle tier: TSD layers (`target_site_duplication`,
  `target_site_duplication_all`)
- Lower tier: per-class with sub-strat (collapsed by default; click to
  expand a class to see its three variants `<class>` / `<class>__young`
  / `<class>__old`)

OR (simpler): replace the `<select>` with a search-as-you-type filter
input. Most classes are TE family names; typing 3-4 chars instantly
narrows. Same UI as the existing chromosome picker on page 1.

**Estimate**: 2-3 hours to implement either approach. The latter is
simpler.

### 2. Class summary panel: surface young/old separately (MEDIUM priority)

The class summary auto-scan currently treats `Gypsy_LTR_retrotransposon`
and `Gypsy_LTR_retrotransposon__young` as two separate classes. They both
get q99/Δ scores. That's correct but visually noisy — you see Gypsy listed
three times (all/young/old) in the flagged list.

**Suggested fix**: group the three variants into one row with sub-cells:

```
class                 q99(all) | q99(young) | q99(old) | Δ vs chrom
Gypsy_LTR_retrotrans   0.45    |   0.32     |   0.13   |  +0.18
```

Then the sort by q99 is unambiguous. Click a sub-cell to set THAT
variant as active in the panel below.

**Estimate**: 4-5 hours. Touches the `_classSummaryForCandidate` and
`_renderClassSummaryPanel` functions, plus a small CSS tweak.

### 3. Side-by-side panel (HIGH priority for cross-strat analysis)

You'll often want to compare `<class>` vs `<class>__young` on the same
window range. Currently you have to flip back and forth.

**Suggested fix**: extend the repeat density panel from one track to up
to three tracks, stacked. The class picker becomes a multi-select. Each
track gets its own y-scale option.

**Estimate**: 6-8 hours. Visual polish at the candidate flank is the
hard part — the existing tight-LOESS rendering gets called per-track.

### 4. IndexedDB cache management (MEDIUM priority)

55 JSONs × ~7 MB ≈ 385 MB total. Cache-restore on page load currently
does all-or-nothing — for the TEfull layer it'll add ~5 sec to each
page load. Plus IndexedDB has per-origin quota limits that some
browsers cap at 200 MB on certain devices.

**Suggested fix**: per-chromosome lazy load. The cache stores all 55
JSONs, but only loads the chrom JSON corresponding to the currently
viewed chromosome. Switching chroms triggers a fresh load (10-50 ms).

**Estimate**: 3-4 hours. Touches `_idbRestoreAll` and `_idbPersistChrom`.

### 5. Cross-species TEfull comparison (LOW priority — research feature)

The most ambitious extension: load both Gar and Mac TEfull JSONs for
syntenic regions and show difference tracks. Requires synteny BEDs
(from the F1 hybrid genome paper) wired in.

**Estimate**: half-day of work. Should defer until at least items 1-3
are done.

---

## Recommended implementation order

If you have one half-day session:

- **Item 1** (class dropdown filter) — 2-3 hours, by far the highest UX
  payoff. Ships immediately.

If you have a full day:

- **Item 1 + Item 3** (side-by-side panel) — the side-by-side panel makes
  the young/old split actually useful. Without it, the strat is mostly
  decorative. With it, you can directly see "young Gypsy spikes flank
  this candidate while old Gypsy is genome-wide" — that's a manuscript-
  quality observation.

If you have two days:

- **Item 1 + Item 2 + Item 3** — the trifecta. Class summary becomes
  young/old aware, side-by-side comparison works, picker is searchable.
  After this, the TEfull layer is fully usable to its potential.

Item 4 (IDB lazy load) only matters once all 55 JSONs are loaded
simultaneously — defer until you actually hit the quota issue.
Item 5 (cross-species) is a research-grade extension; defer until the
cross-species story is the bottleneck.

---

## What to do BEFORE starting any of these

Validate the layer is rendering correctly. Drag a single Gar JSON for
a chromosome with a known polymorphic inversion — the manuscript notes
LG06 (~17 kb inversion, EDTA TEs flanking) and LG22 (CACTA_TIR at
breakpoint). For each, check:

1. Does the candidate's class summary auto-scan flag the expected
   classes? (LG22 should flag CACTA_TIR_transposon and likely its
   `__young` variant.)
2. Does `young_TE_all` show spikes at the candidate flanks?
3. Does the per-class q99 of the flagged class look meaningful (≥ 0.10)
   or near-zero (signal too sparse to render)?

If the answer to any is "no", the underlying data has a problem and the
atlas-side polish work is premature. If "yes", proceed with item 1 of
the implementation order above.

---

## Files this handoff will touch

When the work happens, the changes will be in `Inversion_atlas.html`
sections:

- `_renderRepeatDensityPanel()` (~line 13700-14600)
- `_renderClassSummaryPanel()` (~line 14707-14814)
- `_classSummaryForCandidate()` (search the script for this name)
- `_idbRestoreAll`, `_idbPersistChrom` for item 4

No schema changes needed — the v2 schema accommodates everything in the
TEfull JSONs. No new tabs needed — this is all boundaries-page work.

---

## Open questions for next session

1. Does drag-drop validation show the layer rendering correctly? (gate
   for any further atlas work)
2. Search-as-you-type vs 3-tier dropdown for item 1?
3. Is item 3 (side-by-side panel) worth the 6-8 hours, or is item 1 +
   item 2 enough? Depends on how often you'd actually compare strats.
4. For item 5 (cross-species), does the F1 hybrid synteny BED already
   live somewhere on LANTA, or does it need to be re-emitted from the
   Cactus HAL?

Bring these to whichever session takes the integration work.
