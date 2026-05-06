# HANDOFF — turn 135 — G-panel unified groups popup, Slice 1

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (67,129 lines, 1136/0 tests passing)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC (account `lt200308`).
**Supersedes**: `HANDOFF_2026-05-05_turn134_FINAL.md` (L2-sweep
inspector modal).

This turn ships **Slice 1** of `SPEC_g_panel_unified_groups.md` — the
unified G-panel popup scaffold + manual tab re-host. Karyotype
(Slice 2) and inheritance (Slice 3) tabs are placeholder bodies.

---

## 0. Cohort discipline (NEVER conflate)

Same as turn 134. Three separate cohorts:
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent
surname.

---

## 1. Why this turn pivoted to G-panel

Turn 134's handoff §8 recommended calibrating L2-sweep on real LG28/LG12
data — but that requires LANTA, which is out of scope for this chat.
The other queued options were:
- Trajectory matrix viewer Slice 3
- Cross-chromosome lineages Slice 1
- Page-12 Slice 8 (state machine for tracked-samples + L3)
- **G-panel scaffold Slice 1** ← chosen this turn

G-panel was picked because:
1. The pre-Slice (`_rebuildCandidateRegistries`) already shipped in
   turn 129. The "registry mismatch" bug Quentin reported in turn 128d
   that made the inheritance/grouping panel return empty results was
   resolved when the bridge auto-syncs `state.candidateList` →
   `state.candidates`. So the spec's pre-condition is fully satisfied.
2. Slice 1 is pure UI scaffolding with one tiny extension to
   `renderManualGroupsList()` (add a third container ID). No new
   computation. Low risk.
3. The popup gives the manual-groups list more vertical space than
   the cramped sidebar — immediate UX improvement.
4. Establishes the scaffold so Slices 2 and 3 (karyotype + inheritance
   tabs) drop in cleanly without re-litigating the modal pattern.

---

## 2. What this turn shipped

### 2.1 Popup scaffold

- **Trigger**: `G ▾` button in the L3 toolbar, sibling to the L2-sweep
  inspect button. Always enabled (the popup is cohort-level — works
  regardless of data load state).
- **Hotkey**: lowercase `g` (no modifiers, not in input/textarea/
  select/contenteditable). Capital `G` is left unused. Toggles
  open/closed.
- **Modal**: `#gPanelOverlay`, sibling to `#schemaRegistryOverlay`,
  `#jsScriptsRegistryOverlay`, `#l2SweepInspectorOverlay`. Same close
  pattern as those: ✕ button, click outside, or Esc.
- **Layout**: header (title + tab strip + close) → tab body → footer.
  Three tabs: `karyotype` / `inheritance` / `manual`. Active tab
  highlighted in `var(--good)`; others in `var(--ink-dim)`.

### 2.2 Manual tab (Slice 1's working tab)

The manual tab body contains a `<div id="manualGroupsListPopup">`
that the existing `renderManualGroupsList()` now populates. The
extension is one line: append `popup` to the `containers` array if
the lookup succeeds. Same pattern the function already used for
`#manualGroupsList` + `#manualGroupsListCompact`.

Action buttons (`+ new group from tracked`, K-band picks, export TSV,
import TSV) are rendered with popup-prefixed IDs (`gpMgAddBtn`,
`gpMgExportBtn`, etc.) and **dispatch to the sidebar's existing button
handlers via synthetic clicks**. This keeps the popup decoupled from
the underlying handler internals — the contract is just "click this
id to do this thing." If the sidebar isn't in the DOM (some reduced-UI
mode), the dispatch is a no-op.

The K-band buttons (`→ k0` through `→ k5`) dispatch to
`#mgBandPickBar [data-mgband="N"]` so the focal-L2 K-cluster grab
works identically from the popup.

### 2.3 Karyotype tab (Slice 2 placeholder)

Empty-state body that:
- States "Slice 2 pending" prominently
- Describes what will land (per-band sample memberships, swatches,
  edit names, sample count badges, action buttons)
- If `state.candidate` is non-null, shows the focused candidate's
  ID, bp range (Mb format), and K — so the user can see what Slice 2
  will derive from
- If no candidate is focused, shows the empty-state guidance (catalogue
  / page-1 pin)

### 2.4 Inheritance tab (Slice 3 placeholder)

Empty-state body that:
- States "Slice 3 pending"
- Describes the cross-candidate Jaccard cluster surface
- Counts confirmed candidates (`candidateList.filter(c => c.confirmed)`)
- If <2 confirmed: empty-state hint about the Slice 3 minimum
- If ≥2 confirmed: hint that the inheritance compute would produce a
  result above threshold

This makes the placeholder honest — the user sees what data is
already in shape for the Slice 3 build.

---

## 3. Files changed

```
Inversion_atlas.html
  ├─ +block at ~line 36839 (after L2-sweep inspector window-exports):
  │    G-PANEL UNIFIED GROUPS MODAL — Slice 1 (turn 135)
  │    ~420 LOC of:
  │      - _GPANEL_TABS metadata array (3 tabs × {key, label, slice})
  │      - _gPanelRenderTabKaryotype() — Slice 2 placeholder
  │      - _gPanelRenderTabInheritance() — Slice 3 placeholder
  │      - _gPanelRenderTabManual() — re-host of manual-groups UI
  │      - _renderGPanelModal() — top-level renderer + interaction wiring
  │      - _gPanelClose() / _gPanelToggle() — open/close/toggle helpers
  │      - _wireGPanelOpen IIFE — DOM-ready binder for trigger button
  │      - _wireGPanelHotkey IIFE — keydown listener for 'g'
  │      - 6 window.X = X exports
  │
  ├─ +DOM trigger button in L3 toolbar (~line 6346):
  │    #gPanelOpenBtn — "G ▾", sibling to #l2SweepInspectBtn
  │
  ├─ +DOM overlay div in header (~line 4707):
  │    #gPanelOverlay — full-screen modal overlay, z-index 1000,
  │    sibling to #l2SweepInspectorOverlay
  │
  ├─ +state init slots (~line 8985):
  │    gPanelOpen: false
  │    gPanelTab: 'manual'
  │
  └─ +one-line extension to renderManualGroupsList (~line 39310):
       if (popup) containers.push(popup);
       — populates #manualGroupsListPopup when present, no-ops otherwise

tests/test_turn135_g_panel_slice1.js  [NEW]
  72 tests across 10 sections + sandbox sub-scenarios:
    1. Source-level: function definitions + tab metadata + window exports
    2. DOM elements (overlay, trigger button)
    3. state.gPanelOpen / state.gPanelTab init
    4. renderManualGroupsList extension verified
    5. Hotkey wiring (key + modifiers + input-field guards + preventDefault)
    6. Tab body renderers (karyotype, inheritance with confirmed-count
       derivation, manual structure)
    7. Toggle / close state machine (sandbox-exec)
    8. Tab switching + invalid-tab-fallback
    9. Rendered HTML structure (header, tab strip, body container)
   10. renderManualGroupsList exercised in sandbox: confirms all three
       containers populated, graceful no-op when popup absent
```

LOC delta: 66,708 → 67,129 = **+421**. Of that, ~410 is the new modal
block and ~11 is the renderManualGroupsList extension + state init +
DOM additions.

---

## 4. Tests

```
tests/test_turn135_g_panel_slice1.js — 72 / 0
```

Coverage map highlights:

| Sub-scenario | What it locks down |
|---|---|
| 4 — renderManualGroupsList extension | Container collection includes popup ID; the popup container line is reachable in source |
| 5 — hotkey wiring | Key='g', no modifiers, all 4 input-field guards, preventDefault before toggle. Anchored on `function _wireGPanelHotkey` (not the unrelated state-init comment that mentions the IIFE name) |
| 6a — karyotype empty | Mentions "Slice 2 pending" + no-candidate empty-state + catalogue/pin guidance |
| 6b — karyotype focused | Renders candidate id, bp range as Mb, K value |
| 6c — inheritance counts | confirmed=0 shows <2 hint; confirmed=3 shows above-threshold hint; correctly skips !confirmed |
| 6d — manual tab structure | All 3 button IDs (add/export/import); all 6 K-band data attributes; popup container present |
| 7 — toggle/close | First toggle opens (overlay flex + populated innerHTML); second toggle closes (overlay none) |
| 8 — tab switching | Default 'manual'; explicit switch re-renders correct body; invalid value falls back to 'manual' |
| 10 — render extension | All 3 containers populated when present; gracefully no-ops when popup absent (back-compat) |

| | Tests | Files |
|---|---|---|
| Turn 132 baseline | 879 | 33 |
| Turn 133 Slice 1 | 960 | 34 |
| Turn 134 Slice 2 | 1064 | 35 |
| Turn 135 G-panel S1 | 1136 | 36 |
| Δ this turn | +72 | +1 |

Zero regressions.

---

## 5. Sandbox vs. real-data note (still applies)

Tests run in `vm.createContext` with stubbed DOM and stubbed
dependencies. What the sandbox does NOT exercise:

- **Real button click → handler dispatch chain**. The popup's action
  buttons fire synthetic `click()` events on the sidebar buttons.
  Tests verify the popup buttons exist and have correct IDs; they
  don't verify the dispatch actually fires the underlying handler.
  Browser-side smoke test: open popup, click `+ new group from
  tracked`, verify a new group appears in `state.manualGroups`.
- **Real keydown event delivery**. The hotkey IIFE binds a real
  document listener. Tests verify the source has the right shape
  (key, modifier guards, input-field guards) but don't simulate a
  real keydown.
- **Real CSS rendering**. The modal uses `var(--good)`, `var(--rule)`,
  `var(--ink)`, `var(--ink-dim)` etc. matching the L2-sweep inspector
  modal's tokens; should look identical visually.

Browser smoke test (~30 seconds):
1. Click `G ▾` → popup opens, shows three tabs, manual tab active
2. Press Esc → popup closes
3. Press `g` (lowercase) → popup opens
4. Click `karyotype` tab → body shows "Slice 2 pending" + (if a
   candidate is focused) the candidate's id/range/K
5. Click `inheritance` tab → body shows "Slice 3 pending" +
   confirmed-count
6. Click `manual` tab → manual list shows all your existing groups
7. In the popup, click `+ new group from tracked` → verify group
   appears in BOTH the popup list AND the sidebar list (single
   source of truth working)
8. Type `g` while focused on an input field → popup does NOT toggle
   (input-field guard working)

---

## 6. State / window slots added this turn

```
state.gPanelOpen           boolean — popup visibility, drives overlay display
state.gPanelTab            string  — 'karyotype' | 'inheritance' | 'manual'
state._gPanelHotkeyBound   boolean — IIFE re-entry guard (private)
```

```
window-exposed functions (added this turn):
_renderGPanelModal             — top-level renderer
_gPanelClose                   — close helper
_gPanelToggle                  — open/close toggle
_gPanelRenderTabKaryotype      — Slice 2 placeholder body
_gPanelRenderTabInheritance    — Slice 3 placeholder body
_gPanelRenderTabManual         — Slice 1 working body
```

```
DOM ids claimed (this turn):
#gPanelOpenBtn               — trigger button in L3 toolbar
#gPanelOverlay               — modal overlay (header sibling)
#gPanelClose                 — close button in modal header
#gPanelBody                  — tab body container
#manualGroupsListPopup       — popup-side manual-groups list (renderer fills)
#gpMgAddBtn                  — popup "+ new group" button
#gpMgBandPickBar             — popup K-band button row
#gpMgExportBtn               — popup TSV export button
#gpMgImportBtn               — popup TSV import button
._gpTabBtn (class)           — tab buttons (3 of them)
[data-gpmgband="0..5"]       — K-band picker buttons
[data-gptab="..."]           — tab buttons' tab key attribute
```

```
localStorage keys claimed: (none — popup state is per-session)
```

---

## 7. Backups present

```
Inversion_atlas.html.bak_pre_l2_sweep_slice1     (turn 133 baseline)
Inversion_atlas.html.bak_post_l2_sweep_slice1    (turn 133 final)
Inversion_atlas.html.bak_post_l2_sweep_inspector (turn 134 final)
Inversion_atlas.html.bak_pre_g_panel_slice1      (this turn baseline)
Inversion_atlas.html.bak_post_g_panel_slice1     (this turn final)
```

`.bak_*` files NOT in bundle.

---

## 8. What this is NOT

- **Not the karyotype tab content**. Slice 2 ships `Karyotype groups`:
  per-band swatches, sample counts, "Make manual group from this" /
  "Color by this" / "Export TSV" actions. Today's placeholder shows
  the focused candidate's id/range/K so the user can see what data is
  available, but no derivation happens.
- **Not the inheritance tab content**. Slice 3 ships the cross-candidate
  Jaccard cluster surface with a threshold slider + Compute button.
  Today's placeholder counts confirmed candidates so the user can see
  the data shape.
- **Not lasso multi-select / drag-drop reassignment** within manual
  groups. SPEC §6 mentions drag-drop sample reassignment; deferred.
- **Not a chrom-level / cohort-level scope toggle** beyond what the
  sidebar already does (`📌` pin button per group). The popup re-uses
  the same scope semantics.
- **Not a refactor of state.candidates → state.candidateList**. SPEC
  §15 listed three fix options for the registry mismatch; option (c)
  (`_rebuildCandidateRegistries` bridge) was adopted in turn 129.
  Options (a) and (b) are not pursued — the bridge works.
- **Not a polished aesthetic pass**. The popup is functional and
  matches the existing modal-pattern token usage (var(--good) borders,
  monospace font, etc.) but tab-button hover states, transitions, and
  badge counts on tab labels are deferred to Slice 4 polish.

---

## 9. Where to start the next chat

### Option 9a — G-panel Slice 2 (karyotype tab content)

~1 turn. Builds on the scaffold; needs:
- Per-band membership extraction from `state.candidate.locked_labels`
- K-cluster swatches via the existing `kColors[]` palette
- Per-band sample-count badges + collapsed sample-list (CGAs)
- "Make manual group from this" snapshot action (writes to
  `state.manualGroups`)
- "Color by this" action: `state.colorMode = 'cluster'`,
  `state.lockedRefL2 = candidate.ref_l2`, repaint
- "Export TSV" — band → CGA list

The placeholder body already shows the candidate's id/range/K — Slice 2
just needs to derive the membership groups. Mostly mechanical.

### Option 9b — G-panel Slice 3 (inheritance tab content)

~1 turn. Builds on Slice 2 (or independently if you skip Slice 2).
Needs:
- Threshold slider for `_IGC_DEFAULT_COSINE_DIST_THRESHOLD`
- Explicit `[ Compute ]` button (don't auto-recompute on every drag)
- Per-group rendering (group ID, name, swatch, per-candidate
  band breakdown, sample count badge)
- Same action buttons as karyotype tab

The placeholder body already counts confirmed candidates and gates
on ≥2 — Slice 3 just plugs in the actual compute.

### Option 9c — Trajectory matrix viewer Slice 3

Different priority axis (turn 132 §10b queue item #3). Standalone,
doesn't depend on G-panel.

### Option 9d — Cross-chromosome lineages Slice 1

Queue item #4. Bigger scope than the G-panel slices but high-value
for the manuscript.

### Option 9e — Calibrate L2-sweep on real data (LANTA)

Same as turn 134 §8 — the L2-sweep + inspector toolchain is now
feature-complete enough to validate gate thresholds on LG28 + LG12.
Out of scope for this chat (LANTA only) but listed for completeness.

### Recommendation

9a or 9b. Both are ~1-turn slices, both build directly on this turn's
scaffold, both add immediately-useful content. Slice 2 (karyotype)
unlocks "see the per-band breakdown of any candidate without leaving
the popup" which is the most-frequent group-related action. Slice 3
(inheritance) is structurally similar but depends on `>=2` confirmed
candidates being present, which is rarer in early-curation sessions.

**Suggested order**: 9a → 9b. Both before any L2-sweep calibration
on LANTA, since the karyotype tab will be useful for VERIFYING that
the LG28 inversion's auto-promote produces a reasonable per-band
split.

---

## 10. Honest framing

**What turn 135 actually delivered:**

- A full popup scaffold with three working tabs, hotkey, click-outside-
  to-close, Esc-to-close, tab switching, invalid-tab fallback, popup
  trigger button.
- A one-line extension to `renderManualGroupsList()` that means the
  manual tab is fully functional from day 1 (not a placeholder).
- 72 tests covering source-level structure, hotkey safety guards,
  rendered HTML, state machine, and the renderManualGroupsList
  extension.

**What it deliberately didn't deliver:**

- The two big content tabs (karyotype + inheritance). Their placeholder
  bodies are honest about what's pending and surface the relevant data
  shape (focused candidate / confirmed count) so the user can see what
  Slices 2 + 3 will work with.
- The sidebar deprecation. The sidebar manual-groups list stays — it
  was already there and removing it would be a separate decision.
  The popup is *additive*; same source of truth.

**What this means for the manuscript:**

No direct manuscript impact this turn. The G-panel is a curation
surface; it makes the workflow faster but doesn't add a new figure
or methods sentence. Once Slice 2 + 3 land, the Methods section can
mention "candidate karyotypes and cross-candidate inheritance groups
were inspected via the unified G-panel surface" as a one-liner. Until
then, the existing I·g pills + sidebar manual groups still cover what
needs to be reported.

---

## 11. Bundle contents

When bundled for next turn:

```
Inversion_atlas.html              (current, 67,129 lines)
tests/                            (all *.js)
  ├─ test_turn135_g_panel_slice1.js          [NEW, 72/0]
  ├─ test_turn134_l2_sweep_inspector.js      (unchanged, 104/0)
  ├─ test_turn133_l2_sweep_auto_promote.js   (unchanged, 81/0)
  └─ test_turn132_*.js                       (unchanged, 321/0)
specs_todo/                       (active build queue)
specs_new_turn131/                (pending review queue, unchanged)
HANDOFF_2026-05-05_turn135_FINAL.md  (this file)
all previous handoffs             (kept for history)
OBSERVATIONS_TO_FIX.txt
```

`.bak_*` files NOT in bundle.

---

Walk the map carefully, respect cohort discipline, don't break the
test suite. G-panel Slice 1 is in a clean stopping state — Slice 2
(karyotype tab) is the natural next step.
