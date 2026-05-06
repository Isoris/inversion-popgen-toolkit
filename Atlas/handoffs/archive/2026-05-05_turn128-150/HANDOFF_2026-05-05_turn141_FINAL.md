# HANDOFF — turn 141 — Per-sample-lines candidate-band highlights (Slice 1)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (68,343 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC (account `lt200308`).
**Supersedes**: turn 140 H-label chip handoff.

This turn ships **Slice 1 of `SPEC_lines_panel_candidate_bands.md`**
(spec lives at `specs_todo/`) — faint full-height vertical bands behind
the per-sample-lines panel, one band per confirmed candidate on the
current chrom, painted in distinguishable colors so the user can tell
multiple co-occurring inversion systems apart at a glance without
constantly looking up at the candidate strip.

> Quentin (turn 129):
> *"If we have 2 or 3 inversion systems in the per sample lines we must
> draw the interval of the candidate taking 1/3 vertical space for their
> highlights alpha background you could use like yellow same as now one
> green one a bit blue?"*

---

## 0. Cohort discipline (NEVER conflate)

Same as turn 132 handoff:
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok). Never invent surname.

---

## 1. SPEC reading resolution (§2.1)

The SPEC notes two readings of "1/3 vertical":

- **Reading 1**: stacked thirds (each band gets 1/3 of vertical space)
- **Reading 2**: full-height at low alpha, with "1/3" describing the
  three colors

This turn ships **Reading 2** per the SPEC's own resolved approach:
geometrically clean on overlap, scales past 3 candidates via
golden-angle palette cycle, doesn't bisect sample lines, matches the
existing yellow active-candidate convention. If Quentin meant Reading 1
the spec needs revision and this turn is wrong; that's the only open
SPEC question this turn doesn't resolve.

---

## 2. What's new

### 2.1 State + persistence (~3 LOC)

State init slot at line ~9132:

```js
linesPanelCandidateBands: true,   // turn 141 Slice 1
```

localStorage key: `pca_scrubber_v3.linesPanelCandidateBands`. Default-ON
contract identical to the lineage strip (turn 130 Slice 2): explicit
`'0'` disables, anything else (including missing key) keeps the
default. Setter `setLinesPanelCandidateBands(b)` writes the key and
fires `drawLinesPanel()` redraw.

### 2.2 Helpers (~150 LOC near `setLinesLineageStripOn`)

- `_LINES_CAND_BAND_PALETTE = ['#F5C518', '#4CAF50', '#3B82F6']`
  (yellow / green / blue, hex per SPEC §2.2 verbatim).
- `_candidateBandColor(idx, alpha=0.10)` — returns CSS rgba/hsla
  string. First three indices use the canonical palette; idx ≥ 3
  cycles via golden-angle HSL rotation (137.5° step from baseHue 217°)
  with saturation 60% / lightness 55% to match the canonical pastels.
  Defensive: negative / non-numeric idx normalises to 0; alpha clamped
  to [0, 1].
- `_paintCandidateBands(ctx, opts)` — paints full-height fillRects for
  every confirmed on-chrom candidate. Walks `state.candidateList` in
  array order (matches the candidate-strip order Quentin already
  reads). Returns the count of painted bands. Defensive against null
  ctx, missing toX, invalid mb range, malformed bp values, non-positive
  plot dims, and missing chrom (when chrom is null, cross-chrom filter
  is disabled — useful for tests).

All three exposed on `window.*` for testability.

### 2.3 Toolbar checkbox (~10 LOC HTML)

`<input id="linesCandBandsToggle" checked>` next to the lineage toggle
in `#linesYsourceBar`. Tooltip explains the colors and that drafts get
the existing yellow active-overlay (this layer is confirmed-only).

### 2.4 Handler wireup (~17 LOC)

Inside the existing IIFE that wires `linesLineageStripToggle`. Same
default-ON localStorage restore pattern.

### 2.5 `drawLinesPanel` integration (~17 LOC)

Slot the paint call **before** the frame stroke in the per-sub-panel
loop, just after `state.__linesGeom[source]` stash. Gated on
`state.linesPanelCandidateBands !== false`. Wrapped in try/catch so a
malformed candidate can't take out the rest of the panel — band
drawing is purely additive.

Draw order (SPEC §2.4) is now:

```
ctx.clearRect → bands → frame stroke → y-ticks → lines → tracking → cursor
```

Bands paint per sub-canvas (one canvas per source: PC1 / PC2 / GHSL /
het). Each candidate's bp range maps through that sub's local `toX`,
so vertical alignment with the candidate strip stays exact in every
sub-panel and across zoom.

### 2.6 Stable palette index

Off-screen candidates (clipped out of `[mbMin, mbMax]`) are not
painted but **still bump the palette index counter**, so palette
assignment stays stable when the user pans/zooms — yellow doesn't
silently jump to the second candidate when the first scrolls out of
view.

---

## 3. Tests

`tests/test_turn141_lines_panel_candidate_bands.js` — **62 / 0** across
6 sections:

1. State init + localStorage (3)
2. Setter + helper definitions (12)
3. Toolbar checkbox + change-handler wireup (5)
4. drawLinesPanel integration including draw-order check (6)
5. Behavioral execution via vm sandbox (33)
6. Spec discipline (3)

Behavioral coverage includes: palette content; canonical 3-color hex;
golden-angle hsla cycle for idx ≥ 3; alpha clamping; negative /
non-numeric idx normalisation; correct count for 3 confirmed on-chrom;
draft filtering; cross-chrom filtering; empty list; null ctx; missing
toX; invalid mb range; non-positive plotW; null chrom (cross-chrom
filter disabled); off-screen-candidate-bumps-index palette stability;
custom alpha=0.05 propagation to fillStyle; malformed bp graceful skip.

Adjacent test suites (turn 132 lines renderer, turn 130 lineage UI,
turn 139 H-label classifier, turn 140 H-label chip) still **0 fails**.

Full turn-numbered suite: **1487 PASS / 0 FAIL** across 40 files (up
from 1425 / 0 at turn 140; 1487 = 1425 + 62 new).

---

## 4. Atlas state

| | LOC | Tests | Files |
|---|---|---|---|
| Pre-turn (turn 140 baseline) | 68,120 | 1425 | 39 |
| Post-turn (this) | 68,343 | 1487 | 40 |
| Δ | +223 | +62 | +1 |

LOC budget came in close to the SPEC §5 estimate (helpers ~150 +
integration ~17 + state/setter ~30 + handler ~17 + HTML ~10 ≈ 220).

---

## 5. What this is NOT

- **Not a per-sub-panel toggle.** All sub-panels (PC1 / PC2 / GHSL /
  het) get the same band background. SPEC §6.2 noted "drafts vs
  confirmed" as the per-sub axis to think about; per-sub band toggling
  isn't requested.
- **Not a stripe overlay on top of lines.** Bands are pure background
  (z-order: paint first, before frame stroke). They never obscure data.
- **Not on the L3 mini-PCAs.** SPEC §7 explicit out-of-scope.
- **Not a Slice 2.** SPEC §9 lists potential Slice 2 work
  (cross-chromosome bands, per-band labels inside the band, hover
  tooltips, cursor-following labels). None of those ship here.
- **Not on page 12 (`#thLinesPanel`).** Page 12's `_drawThLinesPanel`
  is a separate render path (turn 132 Slice 5). If Quentin wants the
  same bands there, it's a small follow-up — call `_paintCandidateBands`
  with page-12's local `toX` / `mbMin` / `mbMax`. Estimate: ~0.2 turns.
- **Not Reading 1 of "1/3 vertical".** If Quentin meant stacked thirds
  (one third per candidate, vertically), this turn is wrong and SPEC §2
  needs to be revised. SPEC §6.1 flagged this as the open question.

---

## 6. SPEC §6 open questions — resolved status

| # | Question | Resolution |
|---|---|---|
| 1 | "1/3 vertical" reading | **Reading 2** (full-height low-alpha) per SPEC §2.1 default. Confirm with Quentin. |
| 2 | Drafts get a band? | **No.** Drafts use existing yellow active-overlay. SPEC §3 default. |
| 3 | Palette ordering | **Array order** of `state.candidateList` (matches candidate-strip ordering). SPEC §6.3 default. |
| 4 | >3 candidates | **Cycle palette via golden-angle HSL rotation** (137.5° step from baseHue 217°). SPEC §6.4 default. |
| 5 | Active-candidate yellow vs band | **Bands paint first, active overlay on top at higher alpha.** No special-casing — natural alpha-stacking handles it. SPEC §6.5 default. |

Items 1 and 4 are the ones most likely to need adjustment after
Quentin sees the result on real data.

---

## 7. Backups

```
Inversion_atlas.html.bak_pre_lines_cand_bands  (pre-turn baseline)
```

`.bak_*` files NOT in bundle to keep size down. Re-derivable from git.

---

## 8. Files in this turn

- `Inversion_atlas.html` — modified (5 sites, see §2)
- `tests/test_turn141_lines_panel_candidate_bands.js` — new (62/0)
- This handoff
- (No SPEC moved out of `specs_todo/` — Slice 1 ships but Slice 2
  scope still lives there; move when Quentin signs off on Reading 2.)

---

## 9. What's next

### Option 9a — Polish this slice based on browser smoke test

Open the atlas with LG28 + ≥2 confirmed candidates. Verify:
- Faint yellow / green / blue bands appear behind lines
- Bands track candidate bp ranges exactly
- Toggling `[x] cand bands` removes them and the choice persists across F5
- All sub-panels (PC1/PC2/GHSL/het if loaded) get the same bands
- Active-candidate yellow overlay still pops on top of any underlying band

If Reading 1 (stacked thirds) was actually wanted, revise SPEC §2 and
re-spec the geometry (the helpers stay; only `_paintCandidateBands`
needs new vertical math).

### Option 9b — Mirror to page 12 `#thLinesPanel`

Small (~0.2 turns) follow-up: call `_paintCandidateBands` from
`_drawThLinesPanel` (turn 132 Slice 5). Same helper, same palette,
different `toX` / mb range. Gives page 12 the same multi-candidate
visual on the θπ side.

### Option 9c — Pick up the next item from the build priority queue

Per `specs_new_turn131/SPECS_TIER_INDEX.md` recommended order, the
unbuilt items in priority order are:

1. **G-panel scaffold Slice 3+** (turns 135–136 shipped Slices 1–2)
2. **Trajectory matrix viewer Slice 3**
3. **Cross-chromosome lineages Slice 1**
4. **Multi-chrom load orchestrator**
5. **Genome-wide ideogram**
6. **Marker panel design atlas**
7. **Per-candidate breeding card**
8. **Manuscript bundle export**

L2-sweep auto-promote (was #1) shipped turns 133–134; G-panel was the
item built turns 135–136.

### Option 9d — Tier 1 producer specs (R-side / LANTA)

Several R-side specs are ready to build but blocked on this sandbox
because R isn't installed:

- `SPEC_STEP_T05_theta_cusum_emitter.md` — feeds the page-12 hero
  (renderer shipped turn 132 Slice 3)
- `SPEC_STEP_R39_theta_pi_per_window_emitter.md` — TODO companion
- `SPEC_inversion_age_atlas_surface.md` — needs producer + atlas-side
  build (the H-label chip's age slot is reserved for this)

These are best built on LANTA where R + the actual producer scripts
already exist.

### Recommendation

**9a** — confirm Reading 2 is right with Quentin via browser smoke
test. The helpers are stable; only the geometry choice is in question.
Once confirmed, **9b** (mirror to page 12) is a fast follow-up that
makes both pages consistent.

---

## 10. Honest framing

**What turn 141 actually delivered:**
- A clean, narrow visual layer that solves the stated diagnostic
  problem (telling co-occurring inversions apart in the lines panel
  at a glance) using the existing draw infrastructure
- Tested helpers that another atlas surface can call (page 12's
  `#thLinesPanel`, future cross-chromosome panels)
- A toggle the user can disable when the bands feel like noise

**What it deliberately didn't deliver:**
- The page-12 mirror (~0.2 turns; deferred to a follow-up so this
  turn's slice stayed atomic)
- Per-band labels, hover tooltips, sliding cursor labels — all SPEC §9
  Slice 2 candidates that need user feedback before they're worth
  building
- Resolution of SPEC §6.1 (Reading 1 vs Reading 2). Picked Reading 2
  per the spec's own default; if wrong, geometry needs a re-spec.

**Manuscript impact:**
- Multi-inversion chromosomes (LG12, LG28) become much easier to read
  in figures pulled from the atlas — the bp-to-feature mapping is now
  visible inside the panel itself, not just in the strip above. Useful
  for "Figure 3 — per-sample PC1 trajectories on LG28 with
  candidate-1, candidate-2, candidate-3 highlighted" -style figures.
- Reading 2 (full-height) is also closer to how IGV / UCSC genome
  browsers display multi-region highlights, so figure conventions stay
  familiar to readers from the genome-browser community.

Walk the map carefully, respect cohort discipline, don't break the
test suite. Page 1 lines panel is in a clean stopping state; pick up
wherever makes most sense.
