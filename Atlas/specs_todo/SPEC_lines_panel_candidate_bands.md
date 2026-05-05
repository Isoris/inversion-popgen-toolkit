# SPEC — Per-sample-lines candidate-interval vertical band highlights

**Status**: drafted turn 129. Not yet implemented. Estimated ~1–1.5 turns.

**Trigger**: Quentin's request (turn 129):
> "If we have 2 or 3 inversion systems in the per sample lines we must draw
> the interval of the candidate taking 1/3 vertical space for their
> highlights alpha background you could use like yellow same as now one
> green one a bit blue?
>
> Yellow / Green / Blue
>
> So its like 3 intervals but at the same time it's not corresponding to
> the bands below. So idk. But maybe we can understand because the
> candidates are on top in Cand L2 L1 Win."

---

## 1. The diagnostic question

When the user has 2–3 confirmed candidates on the same chromosome, the
per-sample-lines panel currently draws ALL samples × the current PC1 (or
metric of choice) across the whole chromosome — but the candidates'
genomic intervals are only marked in the top-row tracks (Cand / L2 / L1 /
Win). The dense per-sample-lines middle pane has no visual indication of
which bp range belongs to which candidate. When inversion #1 spans
(2–4 Mb) and inversion #2 spans (8–11 Mb), reading the per-sample
patterns at the right place requires constantly looking up to the
candidate strip, then back down to the lines.

The fix: paint a faint vertical band behind the lines, one band per
confirmed candidate, in distinguishable colors, so the user immediately
sees "the swirly stuff between bp 2M and 4M is candidate 1; the swirly
stuff between 8M and 11M is candidate 2".

## 2. Visual design

The vertical bands paint **across the full vertical extent of the
per-sample-lines plotting area** (not 1/3 — see §2.1 for the resolution
of Quentin's "1/3 vertical space" phrasing). Each candidate's band
covers the bp range `[candidate.start_bp, candidate.end_bp]` along the x
axis.

### 2.1 Resolution of the "1/3 vertical" phrasing

Quentin wrote: *"we must draw the interval of the candidate taking 1/3
vertical space for their highlights alpha background you could use like
yellow same as now one green one a bit blue"*.

Two readings:

1. **Stacked**: each band gets ONE-THIRD of the vertical space, stacked
   top-to-bottom (band1 = top third, band2 = middle third, band3 =
   bottom third).
2. **Full-height with low alpha**: each band paints the FULL vertical
   space at low alpha; "1/3" was Quentin describing one of three
   candidates and the colors he wants.

Reading (1) is geometrically awkward — when candidates overlap in bp,
their thirds don't combine cleanly, and the user loses the "this
sample's line is up here" cue because the line is now bisected by
band-stripe boundaries. It also breaks if there are 4+ candidates.

Reading (2) matches the existing yellow-highlight convention (the
candidate strip already paints the active candidate yellow at full
height) and scales to N candidates by cycling the palette. It's also
how genome browsers (IGV, UCSC) render multi-region highlights.

**Resolved approach**: Reading (2) — each candidate paints the FULL
vertical extent of the per-sample-lines plot at low alpha (~0.10). The
"yellow / green / blue" palette Quentin specified gives 3 distinguishable
colors. For 4+ candidates we cycle the palette with a tiny lightness/
alpha shift to prevent collision (or warn the user — see §6 open
questions).

If Quentin meant Reading (1), this spec needs revision — the open
question §6.1 makes it explicit so the implementation doesn't ship the
wrong interpretation.

### 2.2 Color palette

Three distinguishable, low-saturation pastels at α=0.10 background:

| # | Name | Hex | RGB | Rationale |
|---|---|---|---|---|
| 1 | Yellow | `#F5C518` | `(245, 197, 24)` | Matches the existing active-candidate highlight; visually warm |
| 2 | Green  | `#4CAF50` | `(76, 175, 80)` | Universal "second" hue distinct from yellow |
| 3 | Blue   | `#3B82F6` | `(59, 130, 246)` | "Third" hue, matches the cool-anchor in the dosage ramp |

At α=0.10 they look like washed pastels: pale yellow, mint, sky-blue.
Strong enough to read which region is which; weak enough that the
per-sample lines render legibly on top.

For the 4th, 5th, ... candidates: cycle with `hsl()` rotation by 137.5°
(golden angle) so consecutive candidates stay distinguishable, OR show
only the first 3 and a warning chip. **Default: cycle palette** for
graceful degradation.

### 2.3 Stripe behavior on overlap

If two candidates overlap on bp: each paints its full-height band at
α=0.10; the overlapping region naturally renders at the additive blend
(~α=0.19). That's a feature — overlap regions get visibly darker, so
the user sees "these two candidates fight over this stretch". No
special-casing needed.

### 2.4 Z-ordering

The bands paint **first**, before:
- any per-sample lines
- any sample-tracking highlights
- the cursor crosshair
- the inheritance-labels strip

So the bands are pure background. They MUST NOT obscure data lines.

## 3. Which candidates qualify for a band?

Default rule: every entry in `state.candidateList` that has
`confirmed === true` AND `start_bp + end_bp` defined for the current
chromosome. Drafts (`confirmed === false`) do NOT get a band — they're
already shown via the existing yellow draft highlight, and showing a
fourth color for drafts would add noise.

If `state.candidateList` is empty: no bands, lines render exactly as
today. (No-op, no behavior change, no risk.)

## 4. State + persistence

- New slot: `state.linesPanelCandidateBands` (boolean, default `true`).
  When `false`, behaves like today (no bands).
- localStorage key: `pca_scrubber_v3.linesPanelCandidateBands`.
- Toggle in the per-sample-lines panel header: `[x] cand bands`. Lives
  in the existing toolbar near the PC1/PC2 + colorMode controls.

## 5. Implementation checklist (Slice 1, ~1 turn)

- [ ] Add `state.linesPanelCandidateBands: true` to the state init block.
- [ ] Add `linesPanelCandidateBands` localStorage restore at startup.
- [ ] Add `[x] cand bands` checkbox to the per-sample-lines header.
- [ ] In `drawLinesPanel` (or its modern equivalent — verify the actual
      function name first), paint the candidate bands as the very first
      step after frame setup, before any line / dot / track render.
- [ ] Helper `_candidateBandColor(idx, alpha=0.10)` returning the cycled
      palette CSS color.
- [ ] Helper `_paintCandidateBands(ctx, plotRect, bpToX)` that iterates
      `state.candidateList.filter(c => c.confirmed && c.chrom === state.data.chrom)`,
      maps each `start_bp/end_bp` through the panel's bp-to-pixel
      function, and fills a rectangle.
- [ ] Toggle wireup with re-render on change.
- [ ] Tests: source-level (slot + checkbox + helper presence + draw-order)
      + behavioural (palette cycling, alpha format, bp clamp at panel
      edges).

## 6. Open design questions

1. **Reading of "1/3 vertical"**: confirm with Quentin that full-height
   low-alpha (Reading 2) is what's wanted. If he meant stacked thirds
   (Reading 1), revise §2 and re-spec the geometry.
2. **Drafts**: do drafts get a band too? Default answer: no (drafts get
   the existing yellow active-candidate highlight). Confirm.
3. **Palette ordering**: should bands always paint in `start_bp` order
   (so "the leftmost candidate gets yellow")? Or in `candidateList`
   array order (so the order matches the page-2 candidate strip)? Best
   default: array order, since that mirrors the strip the user is
   already reading.
4. **>3 candidates**: cycle palette with golden-angle rotation, OR
   render only first 3 + warning chip "+N more"? Default: cycle with
   rotation; let the user complain if it's too noisy.
5. **Interaction with the existing "active candidate" yellow highlight**:
   if candidate #2 (green band) is the active candidate, should green
   "win" over yellow, or both render? Default: candidate-band paints
   first; the active-candidate yellow highlight paints on top at higher
   alpha so the active one always pops.

## 7. What this is NOT

- **Not a replacement for the candidate strip** (Cand / L2 / L1 / Win
  rows in the top tracks). The strip stays — the bands are an
  additional, redundant cue at the lines-panel scope.
- **Not a track**. Bands are pure background paint, not a height-
  consuming track. The lines-panel vertical real estate is unchanged.
- **Not per-band coloring**. The vertical band corresponds to the
  candidate's **bp interval**, not to a sample's K-cluster band
  assignment. Quentin's note: *"its not corresponding to the bands
  below"* — exactly. The colors are candidate-identity colors, not
  band-identity colors.
- **Not for the L3 mini-PCAs**. Bands only show in the per-sample-lines
  middle pane.

## 8. Tests

Source-level + behavioural for Slice 1:

- `state.linesPanelCandidateBands` slot exists, default `true`.
- localStorage round-trip: set, re-restore.
- `_candidateBandColor(0)` → yellow; `(1)` → green; `(2)` → blue;
  `(3)` → palette cycles (golden-angle rotation).
- `_candidateBandColor(idx, 0.05)` → returns same hue at α=0.05.
- `_paintCandidateBands` filters out unconfirmed candidates.
- `_paintCandidateBands` filters out cross-chromosome candidates.
- `drawLinesPanel` references `_paintCandidateBands` BEFORE any line-
  render call (source-level draw-order check).

## 9. Out of scope (Slice 2 candidates, deferred)

- Cross-chromosome candidate display (currently the lines panel is
  per-chrom; bands are per-chrom for now).
- Per-band labels ("Cand1", "Cand2") inside the band itself — could be
  added as a top-edge text tag if Quentin wants attribution without
  reading the strip.
- Hover tooltips on the bands (which candidate, bp span, confirm
  status).
- Sliding labels that follow the cursor.
