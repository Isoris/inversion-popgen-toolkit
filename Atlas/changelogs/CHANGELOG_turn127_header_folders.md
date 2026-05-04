# CHANGELOG — turn 127: header reorganization with colour-coded folders

**Atlas Δ:** +99 lines (61,924 → 62,023)
**Tests:** 67 new (762/762 cumulative green across 12 test suites)
**No regressions** — all prior turns 117–126 still pass.

---

## What you asked for

From the screenshots you shared:

1. Group the top-row buttons into colour-coded folders (blue / green / yellow)
2. Move the right-hand-side row of buttons (export / matrix / labels / active /
   view) up to the top header row, same size as the others
3. Drop the redundant info from the metadata strip:
   - `active: 226/226` was duplicated (already shown via `226 samples`)
   - `L1=10 · L2=48 · L2 peaks=134` was unnecessary noise

All three shipped.

---

## What changed

### Top-row layout — single row now

Before (two rows):
```
┌──────────────────────────────────────────────────────────────────────────┐
│ Inversion Atlas  C_gar_LG28 · 4302 W · 226 samples · L1=10 · L2=48 · L2  │
│   peaks=134  schema v2 · 3 layers  candidate mode  save  load  free      │
│   reset layout  Inversion Atlas ▾  academic                              │
├──────────────────────────────────────────────────────────────────────────┤
│ 1 local  2 local  2b local  3 candi  ...  16 ho  | export matrix 0       │
│                                                  | filled I·g labels:on  │
│                                                  | active:226/226 ○ view │
└──────────────────────────────────────────────────────────────────────────┘
```

After (one row):
```
┌──────────────────────────────────────────────────────────────────────────┐
│ Inversion Atlas  C_gar_LG28 · 4302 W · 226 samples  schema v2 · 3 layers │
│       <flex spacer>   ● session ▾   ● mode ▾   ● data ▾   Inv. Atlas ▾   │
│                                                            academic      │
├──────────────────────────────────────────────────────────────────────────┤
│ 1 local  2 local  2b local  3 candi  ...  16 ho                          │
└──────────────────────────────────────────────────────────────────────────┘
```

### Three colour-coded folders

Each folder is a button that on **hover** (or `:focus-within`) reveals a
popover panel listing the grouped controls. Same idiom as the existing
`#atlasModeIndicator` four-atlas dropdown — consistent UX across the header.

🔵 **session** (blue) folder contains:
- `💾 save` — save session JSON
- `📂 load` — load session JSON
- `📐 fixed` / free — layout mode toggle
- `↺ reset layout`

🟢 **mode** (green) folder contains:
- `🏗 candidate mode` — toggle interval drafting
- `🪄 auto-fill` — fill haplotype labels for all candidates
- `I·g labels: on` — toggle inheritance-group labels strip

🟡 **data** (yellow) folder contains:
- `📥 export` — full atlas context export
- `🔗 matrix` — Cramér's V matrix modal
- `active: —` — active samples picker
- (server status badge appended dynamically — `○ view` / `● server` / `◌ checking...`)

The yellow folder's inner panel keeps the ID `atlasToolsGroup` so the
dynamic server-badge append code (line ~48140) finds the right target with
zero changes.

### Visual cues

- **Folder button colour** matches the dropdown's left-border tint:
  - blue ≈ `rgba(64, 144, 220, *)`
  - green ≈ `rgba(72, 200, 120, *)`
  - yellow ≈ `rgba(232, 196, 76, *)`
- **Hover** brightens the folder button background and reveals the panel.
- **Active state** — when the green folder contains an active mode button
  (e.g. `candidateModeBtn[data-active="1"]`), the folder button gets a
  subtle gold inset ring so you can see candidate mode is on without
  opening the folder. (Not yet wired by JS — CSS rule exists, the JS
  hookup is a small follow-up if you want it.)

### Metadata strip cleanup

Before: `C_gar_LG28 · 4302 W · 226 samples · L1=10 · L2=48 · L2 peaks=134`
After: `C_gar_LG28 · 4302 W · 226 samples`

The `L1=10 · L2=48 · L2 peaks=134` data is still present in the sidebar
data-status panel (line 41624 area) — just not duplicated in the header.

The `active: 226/226` button now lives in the yellow folder. When 226/226
samples are active it's redundant with `226 samples` in the metadata strip,
but when a subset is active (e.g. `active: 217/226`) the folder button
flips amber so you can see at a glance you're in subset mode. So it's not
exactly the same info, just usually-the-same info.

### Old `#atlasToolsGroup` removed from `<nav>`

The old div sitting at the end of the tab bar (with `margin-left:auto;
display:flex;height:22px;...`) is gone. All the children moved into the
yellow folder with their IDs preserved. No JS event listeners needed
re-wiring — they all attach by ID.

The old legacy CSS rules for `#themeToggleBtn { margin-left: auto }` and
`#resetLayoutBtn { margin-left: auto }` are removed. A single explicit
`<div style="flex: 1;"></div>` spacer in the markup now handles
right-alignment for everything past it (folders + atlas selector + theme
toggle), giving consistent gaps and predictable layout.

---

## Files changed

```
Inversion_atlas.html                      +99 lines (61,924 → 62,023)
test_turn127_header_folders.js            NEW (67 tests)
CHANGELOG_turn127_header_folders.md       THIS FILE
```

No other files touched. No new JSON layers, no schema changes. Pure UI.

---

## What's left for future turns

- **Active-folder-glow JS hookup**: the CSS rule
  `header .header-folder[data-has-active="1"][data-color="green"]`
  exists but no JS sets `data-has-active="1"` on the folder when
  `candidateModeBtn[data-active="1"]` flips on. Small follow-up. The folder
  works fine without it; you just need to open it to see candidate mode is
  active.
- **TSV export columns** for karyo verdict + wfmash refinement state
  (mentioned in turn 126 todo). Not part of this turn.
- **BUSCO integration** is the bigger next thing per your earlier note.
