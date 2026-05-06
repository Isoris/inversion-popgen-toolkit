# HANDOFF — turn 136 — G-panel karyotype tab content (Slice 2)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (67,474 lines, 1237/0 tests passing)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC (account `lt200308`).
**Supersedes**: `HANDOFF_2026-05-05_turn135_FINAL.md` (G-panel scaffold).

This turn ships **Slice 2** of `SPEC_g_panel_unified_groups.md` — the
karyotype tab content. Per turn 135's recommendation §9a.

---

## 0. Cohort discipline (NEVER conflate)

Same as turn 135. Three separate cohorts:
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent
surname.

---

## 1. What this turn shipped

### 1.1 Karyotype tab body (real, not placeholder)

Replaced `_gPanelRenderTabKaryotype()`'s Slice 1 placeholder with full
SPEC §4 content:

- **Header breadcrumb** showing candidate id, bp range (Mb), K, and
  active label vocab (`legacy` or `detailed`)
- **Per-band rows** (one per K-band of `state.candidate.locked_labels`):
  - K-cluster swatch (6-color palette matching `_ACK_BAND_PALETTE`)
  - Label from `getKaryotypeLabel(bandIdx, K)` honoring
    `state.labelVocab` (legacy: `'band 1 (lo)'`, etc.; detailed:
    `'H1/H1'`, `'H1/H2'`, etc.)
  - Sample count badge (`n=N`)
  - Collapsed sample-list preview (first 6 CGAs + count remaining,
    full list as title attribute)
  - Per-row `+ manual` button (disabled for empty bands)
- **Detailed-vocab caveat** — when `state.labelVocab === 'detailed'`,
  the header shows `⚠ operational labels` with the full caveat from
  `getKaryotypeLabelCaveat()` as title, since H-labels are PC1-ordered
  not biologically confirmed
- **Ungrouped sample count** — surfaces when any sample has lab=-1 or
  out-of-K-range
- **Footer actions**:
  - `🔒 color PCA by this candidate` — sets `state.lockedLabels =
    Int8Array(candidate.locked_labels)` + `state.lockedRefL2 =
    candidate.ref_l2`, mirrors the sidebar `#lockColorsBtn` handler,
    closes popup so PCA repaint is immediately visible
  - `⬇ export TSV` — downloads `karyotype.<candId>.<chrom>.tsv` with
    `cga\tband_label` rows

### 1.2 Five new helper functions

All exported on `window.X`:

- `_gpKaryoExtractBands(candidate, samples, K)` — pure derivation;
  returns `{bands: [{bandIdx, label, sampleIdxList, members, color}],
  nUngrouped, K}` or `null`
- `_gpKaryoMakeManualGroup(candidate, bandIdx)` — wraps
  `addToManualGroup(name, sampleIdxList, opts)`; default name
  `<candId>_<bandLabel>`
- `_gpKaryoColorByCandidate(candidate)` — mirrors `#lockColorsBtn`;
  triggers `drawPCA + renderL3Panel + drawLinesPanel + refreshLockBtn
  + refreshCandidateUI`
- `_gpKaryoExportTSV(candidate)` — fail-safe (returns false in
  non-DOM env, no throw)
- `_gpKaryoColor(bandIdx)` — palette indexer

### 1.3 Interaction wiring

`_renderGPanelModal()` now has a `karyotype` branch (sibling to
`manual`) that binds:
- Per-row `_gpKaryoMakeBtn` clicks → `_gpKaryoMakeManualGroup` +
  visual feedback (button text flashes "✓ added", then re-renders
  the modal so the band reflects "promoted to manual" state)
- `#gpKaryoColorByBtn` click → `_gpKaryoColorByCandidate` + popup
  close (so user immediately sees the colored PCA)
- `#gpKaryoExportBtn` click → `_gpKaryoExportTSV`

---

## 2. Vocabulary correction caught during audit

Turn 135's Slice 1 placeholder text said *"as named groups (HOM_REF /
HET / HOM_INV at K=3, six bands at K=6)"*. That was wrong — the atlas
does NOT use HOM_REF/HET/HOM_INV as the K-band vocabulary. The actual
vocab comes from `getKaryotypeLabel(bandIdx, K)`:

- **Legacy (default)**: `'band 1 (lo)'`, `'band 2 (mid)'`, `'band 3 (hi)'`
- **Detailed**: `'H1/H1'`, `'H1/H2'`, `'H2/H2'` (K=3) → up to
  `'H3/H3'` (K=6)

The HOM_REF/HET/HOM_INV vocabulary IS used elsewhere in the atlas
(the lines panel coloring, page 4 karyotype display) but as a
**color-binding convention** (band 0 → blue / HOM_REF, band 1 →
orange / HET, band 2 → purple / HOM_INV), NOT as the displayed band
label. Slice 2 keeps that color binding via the 6-color palette but
displays the vocab-correct label.

Importantly, the atlas already has a documented caveat that the
detailed `'H1/H1'`-style labels are *operational ordering by median
PC1, NOT a confirmed biological assignment* (line 32898). Slice 2
surfaces that caveat in the tab header when detailed-vocab is active.

---

## 3. Files changed

```
Inversion_atlas.html
  ├─ +block at ~line 36930 (after _GPANEL_TABS):
  │    "turn 136 — Slice 2: karyotype tab content" — adds:
  │      - _GPANEL_KARYO_PALETTE constant (6 hex colors)
  │      - _gpKaryoColor(bandIdx)
  │      - _gpKaryoExtractBands(candidate, samples, K)
  │      - _gpKaryoMakeManualGroup(candidate, bandIdx)
  │      - _gpKaryoColorByCandidate(candidate)
  │      - _gpKaryoExportTSV(candidate)
  │
  ├─ ~rewrote _gPanelRenderTabKaryotype (replaces ~38 LOC placeholder
  │    with ~110 LOC full Slice 2 body — header breadcrumb, per-band
  │    rows, ungrouped count, footer actions)
  │
  ├─ +karyotype interaction wiring in _renderGPanelModal (~line 37395):
  │    `if (activeTab === 'karyotype')` branch with per-row +manual
  │    button bindings, gpKaryoColorByBtn, gpKaryoExportBtn
  │
  └─ +5 window.X exports for the new helpers

tests/test_turn135_g_panel_slice1.js
  - Updated 3 test assertions that asserted on the Slice 1 placeholder
    text ("Slice 2 pending", "Currently focused candidate" header).
    These now assert the contracts that survive Slice 2 (no-candidate
    empty-state still says "No candidate currently focused"; Slice 2
    has shipped → no "Slice 2 pending" placeholder text). Net +1 test
    in this file (72 → 73).

tests/test_turn136_g_panel_karyotype_slice2.js  [NEW]
  100 tests across 9 sections + sub-scenarios:
    1. Source-level: 5 helpers + palette constant + window exports
    2. _gpKaryoColor palette indexing (band wraparound, null fallback)
    3. _gpKaryoExtractBands: K=3 happy path, K=6 + ungrouped, out-of-K
       counted as ungrouped, length mismatch, .ind fallback, null
       inputs, empty arrays
    4. _gpKaryoMakeManualGroup: wraps addToManualGroup with default
       name + color; empty band → null + no call; out-of-range → null;
       special-char sanitization in candidate id
    5. _gpKaryoColorByCandidate: snapshot semantics (mutating original
       doesn't affect state.lockedLabels), all 5 repaint hooks called,
       null/missing inputs → false
    6. _gpKaryoExportTSV: returns false in sandbox, real-Blob mock
       confirms one download with correct filename pattern
    7. _gPanelRenderTabKaryotype: no-candidate empty-state, full K=3
       rendering, detailed-vocab caveat, no-locked_labels empty-state,
       empty-band disabled buttons, ungrouped count display
    8. Interaction wiring source-level: karyotype branch present,
       all 3 bindings, color-by closes popup
```

LOC delta: 67,129 → 67,474 = **+345**. Of that, ~310 is the new
helpers + Slice 2 body, and ~35 is the karyotype branch in
`_renderGPanelModal`.

---

## 4. Tests

```
tests/test_turn136_g_panel_karyotype_slice2.js — 100 / 0
tests/test_turn135_g_panel_slice1.js           — 73 / 0  (was 72)
```

Coverage map highlights:

| Sub-scenario | What it locks down |
|---|---|
| 3a — extractBands K=3 | Happy path: 3 bands, 2 members each, sampleIdxList aligned with members |
| 3b — extractBands K=6 + ungrouped | -1 labels excluded, nUngrouped counted, empty bands preserved |
| 3c — out-of-K | Labels > K-1 also count as ungrouped (defensive) |
| 3d — length mismatch | `min(labels.length, samples.length)` — no OOB indexing |
| 3e — .ind fallback | Sample ID resolved via `.cga \|\| .ind` |
| 4d — candidate-id sanitization | Special chars stripped from group name |
| 5a — colorByCandidate | Snapshot is decoupled from original (Int8Array copy); all 5 repaint hooks fire |
| 5d — ref_l2 missing | `state.lockedRefL2 = null` (not `undefined`) |
| 6a — exportTSV sandbox | Returns false in non-DOM env (no throw) |
| 6d — exportTSV real Blob | Filename pattern `karyotype.<candId>.<chrom>.tsv` |
| 7c — detailed-vocab caveat | `⚠ operational labels` indicator + H-labels rendered |
| 7e — empty bands | All K bands rendered; empty bands' `+ manual` buttons disabled; `(no samples)` preview |
| 7f — ungrouped surfacing | "2 samples have label -1" message when present |

| | Tests | Files |
|---|---|---|
| Turn 132 baseline | 879 | 33 |
| Turn 135 G-panel S1 | 1136 | 36 |
| Turn 136 G-panel S2 | 1237 | 37 |
| Δ this turn | +101 | +1 |

Zero regressions in turns 132-134 (321 + 81 + 104 = 506 tests still
green).

---

## 5. Sandbox vs. real-data note

What the sandbox does NOT exercise:

- **Real `addToManualGroup` round-trip** — the test stubs it. Browser
  smoke test: open the popup with a focused candidate, click `+ manual`
  on a non-empty band, verify the new group appears in the manual tab
  AND in the sidebar, with members matching the band's CGAs.
- **Real PCA repaint** — `_gpKaryoColorByCandidate` calls `drawPCA()`,
  `renderL3Panel()`, `drawLinesPanel()` but tests use stubs. Browser
  smoke test: click `🔒 color PCA by this candidate`, verify popup
  closes and PCA dots colorize per the candidate's locked_labels.
- **Real Blob download** — test 6d uses a Blob mock that confirms the
  click handler fires with correct filename, but the actual byte
  content of the TSV isn't inspected (the mock retains `_parts` which
  could be inspected; left as future work if Quentin reports a bug).

Browser smoke test sequence (~1 minute):

1. Pin a candidate on page 1 (lock + promote)
2. Press `g` to open the popup
3. Click the `karyotype` tab → see N rows (K bands), each with swatch,
   label, count, CGA preview, `+ manual` button
4. Click `+ manual` on band 0 → button flashes "✓ added", then becomes
   disabled (because that band is now in the manual list)
5. Click the `manual` tab → see the new group `<candId>_<label>` with
   correct member count
6. Click back to `karyotype` → click `🔒 color PCA by this candidate`
   → popup closes, PCA dots take the candidate's band colors
7. Re-open popup, click `karyotype` → click `⬇ export TSV` → file
   `karyotype.<candId>.<chrom>.tsv` downloads with cga + label rows

If any step fails, document and feed back; most likely failure point
is step 4 (the `+ manual` button → re-render → disabled flow), which
relies on `addToManualGroup`'s side effect updating
`state.candidateList`-equivalent state.

---

## 6. State / window slots added this turn

```
state.lockedLabels       Int8Array — set by "color PCA by this candidate";
                                     SAME slot the existing #lockColorsBtn uses
state.lockedRefL2        number    — same convention as #lockColorsBtn
                                     (already documented in turn 47s+)
```

These are NOT new state slots — they're existing slots that Slice 2
now writes to, just like the sidebar's lock-colors button does. Same
contract.

```
window-exposed functions (added this turn):
_gpKaryoExtractBands       — pure derivation
_gpKaryoMakeManualGroup    — wraps addToManualGroup
_gpKaryoColorByCandidate   — mirrors #lockColorsBtn handler
_gpKaryoExportTSV          — TSV download
_gpKaryoColor              — palette indexer
```

```
DOM ids claimed (this turn):
#gpKaryoColorByBtn         — footer "color PCA by this" button
#gpKaryoExportBtn          — footer "export TSV" button
._gpKaryoMakeBtn (class)   — per-band "+ manual" buttons
._gpKaryoBandRow (class)   — per-band row container
[data-gpkbi="N"]           — band index attribute on per-row buttons
```

```
localStorage keys claimed: (none — Slice 2 only writes to existing
state.lockedLabels / state.manualGroups, both of which persist via
their existing paths)
```

---

## 7. Backups present

```
Inversion_atlas.html.bak_pre_l2_sweep_slice1     (turn 133 baseline)
Inversion_atlas.html.bak_post_l2_sweep_slice1    (turn 133 final)
Inversion_atlas.html.bak_post_l2_sweep_inspector (turn 134 final)
Inversion_atlas.html.bak_pre_g_panel_slice1      (turn 135 baseline)
Inversion_atlas.html.bak_post_g_panel_slice1     (turn 135 final)
Inversion_atlas.html.bak_pre_g_panel_slice2      (this turn baseline)
Inversion_atlas.html.bak_post_g_panel_slice2     (this turn final)
```

`.bak_*` files NOT in bundle.

---

## 8. What this is NOT

- **Not editable band names**. SPEC §4 mentioned editable group names
  in the karyotype tab. Slice 2 uses the canonical
  `getKaryotypeLabel(bandIdx, K)` output. Editable names are deferred
  to Slice 4 polish — the renaming would have to be a per-candidate
  override (since the canonical labels come from a global vocab),
  which is more state plumbing than fits this slice.
- **Not multi-track**. SPEC §4 mentioned `track 1 / 2` toggle for
  multi-track candidates (rare). Slice 2 only reads
  `candidate.locked_labels` (the focal track). Multi-track support
  needs `candidate.tracks[trackIdx].active_bands` plumbing; deferred.
- **Not biological-vs-operational disambiguation**. The detailed-vocab
  caveat surfaces the warning but doesn't offer a way to assert
  ground-truth biological labels. That's a manuscript-level call,
  not an atlas feature.
- **Not the inheritance tab**. Slice 3 still pending. Today's
  inheritance tab still shows the Slice 1 placeholder body
  (confirmed-count + "Slice 3 pending").
- **Not "Slice 4 polish"** — no badges on tab labels showing band
  counts, no transition animations, no per-tab keyboard shortcuts.

---

## 9. Where to start the next chat

### Option 9a — G-panel Slice 3 (inheritance tab content)

~1 turn. Builds on the same scaffold. Needs:
- Threshold slider for `_IGC_DEFAULT_COSINE_DIST_THRESHOLD` with
  explicit `[ Compute ]` button (no auto-recompute on every drag)
- Per-group rendering (group ID, name, swatch, per-candidate band
  breakdown, sample count badge)
- Same action buttons as karyotype tab (Make manual, Color by, Export)
- The placeholder body already counts confirmed candidates — Slice 3
  just needs the compute + render.

The Slice 1 placeholder for inheritance still says "Slice 3 pending".
This turn's pattern (replace placeholder, write helper functions,
wire interaction branch in `_renderGPanelModal`, add tests) directly
applies.

### Option 9b — Trajectory matrix viewer Slice 3

Different priority axis. Standalone, doesn't depend on G-panel.

### Option 9c — Cross-chromosome lineages Slice 1

Bigger scope.

### Option 9d — Calibrate L2-sweep on real data (LANTA)

Out of scope for this chat.

### Recommendation

**9a — G-panel Slice 3 (inheritance tab)**. The G-panel was estimated
at ~3 turns total across 4 slices in SPEC §13; Slice 3 is mid-stack.
Finishing it makes the popup feature-complete (modulo Slice 4 polish).
The inheritance compute is the only "expensive" operation in the
spec, but `_IGC_DEFAULT_COSINE_DIST_THRESHOLD` already exists and the
Compute button is the cheap-defer pattern from SPEC §5.

After Slice 3, the recommended pivot is back to L2-sweep calibration
on real LANTA data (turn 134 §8 / 135 §9e — the toolchain has been
ready since turn 134).

---

## 10. Honest framing

**What turn 136 actually delivered:**

- The karyotype tab body is now FULL — from the candidate's
  locked_labels through to per-band rows, action buttons, ungrouped
  surfacing, vocab caveat handling.
- Five small pure-derivation / wrapper functions, each independently
  testable, each window-exported for future surface re-use.
- Vocab discipline caught at audit time — Slice 1's placeholder
  used HOM_REF/HET/HOM_INV which is the COLOR convention not the
  LABEL vocab. Slice 2 uses the canonical `getKaryotypeLabel`
  output and surfaces the operational-label caveat when
  detailed-vocab is active.
- 100 new tests covering all five helpers' edge cases, full renderer
  output, and the wiring branch's source-level presence.

**What it deliberately didn't deliver:**

- Editable band names, multi-track support, transition animations,
  badge counts on tab labels — Slice 4 polish.
- The inheritance tab — Slice 3.
- A "page-1 navigate to this band's representative" action. Could
  be added if Quentin wants click-row-to-jump-to-PCA, but the
  current design treats the tab as a triage surface, not a
  navigation surface.

**What this means for the manuscript:**

Same story as turn 135 — no direct manuscript impact, but the
karyotype tab makes the per-candidate band review faster. Once
Slice 3 lands and L2-sweep is calibrated on real data, the Methods
section can mention "candidate karyotypes and cross-candidate
inheritance groups were inspected via the unified G-panel surface,
allowing per-band promote-to-manual snapshots and PCA color-locking
for visual concordance checks."

---

## 11. Bundle contents

```
Inversion_atlas.html              (current, 67,474 lines)
tests/                            (all *.js)
  ├─ test_turn136_g_panel_karyotype_slice2.js  [NEW, 100/0]
  ├─ test_turn135_g_panel_slice1.js            (73/0, +1 from turn 135)
  ├─ test_turn134_l2_sweep_inspector.js        (unchanged, 104/0)
  ├─ test_turn133_l2_sweep_auto_promote.js     (unchanged, 81/0)
  └─ test_turn132_*.js                         (unchanged, 321/0)
specs_todo/                       (active build queue)
specs_new_turn131/                (pending review queue, unchanged)
HANDOFF_2026-05-05_turn136_FINAL.md  (this file)
all previous handoffs             (kept for history)
OBSERVATIONS_TO_FIX.txt
```

`.bak_*` files NOT in bundle.

---

Walk the map carefully, respect cohort discipline, don't break the
test suite. G-panel Slice 2 (karyotype) is in a clean stopping state
— Slice 3 (inheritance) is the natural next step, then Slice 4
polish, then pivot to LANTA L2-sweep calibration when bandwidth
allows.
