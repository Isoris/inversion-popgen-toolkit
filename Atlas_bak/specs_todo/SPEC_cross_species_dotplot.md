# SPEC — Cross-species dot plot widget (mashmap multi-resolution + wfmash synteny)

**Status**: forward-looking spec. Not implemented. Created same session as
`SPEC_focal_vs_background_widget.md`, `SPEC_xpehh_track.md`,
`SPEC_ld_decay_overlay.md`, `SPEC_ncrna_density_layer.md`.

**Reading order**: this spec → existing `STEP_CS01_extract_breakpoints.py`
(the cross-species pipeline already emits `synteny_blocks` at schema_version 2
— that's the manuscript-grade primary data source) → atlas
`_renderCrossSpeciesSynteny` (line ~19888, the existing macro-synteny
permutation panel this spec sits next to, not replaces) →
`SPEC_OVERVIEW_multispecies_architecture.md` (multi-species framing) → the
Perl reference script Quentin pasted (we port only the layout algorithm —
LIS-based fattest-alignment diagonal layout — to JS; everything else of that
script's gnuplot/X11/clipboard architecture is dropped).

**Reference idiom**: mummerplot-style dot plot. Reference chromosome on x,
query chromosome on y, alignments rendered as colored line segments (forward
strand on the main diagonal, reverse strand on the anti-diagonal of each
chrom-pair box). Color encodes percent identity on a continuous ramp.

---

## 0. Settled infrastructure — REUSE

### 0.1 wfmash synteny_blocks already in `cs_breakpoints_v1.json` (v2)

`STEP_CS01_extract_breakpoints.py` (the script Quentin pasted) emits, at
schema_version 2:

```jsonc
{
  "synteny_blocks": [
    { "gar_chr": "LG28", "gar_start": 14100000, "gar_end": 16500000,
      "mac_chr": "Mac_27", "mac_start": 8200000, "mac_end": 10600000,
      "strand": "+", "block_size_bp": 2400000, "mapping_quality": 60 },
    ...
  ],
  "chrom_lengths_query":  { "LG28": 38500000, ... },
  "chrom_lengths_target": { "Mac_27": 41200000, ... }
}
```

This is exactly the shape a dot plot needs. **Day-1 mode reads this
directly — no new compute, no new server endpoint.** Each block is one
diagonal line segment in the rendered plot.

### 0.2 Page 16 hosts the existing macro-synteny panel

`_renderCrossSpeciesSynteny()` (atlas line ~19888) already renders a
synteny-block panel + permutation test on page 16 below the focus area.
The dot plot is a **complementary visual** (geometric dot plot vs the
existing block-list/perm-test view), not a replacement. Both render
when v2 data is loaded; user navigates between them via tabs inside the
synteny section.

### 0.3 Atlas page-16 styling, hover infrastructure

The existing `csUI` namespace handles catalogue rendering, focus panel,
chip styling. The dot-plot widget reuses the same CSS variables
(`--accent`, `--ink`, `--ink-dim`, `--rule`, `--panel-2`, `--mono`),
the same chip styles, and the same drag-drop loader pattern.

### 0.4 Hover/popup pattern is new — but small

No existing atlas widget has the hover-enlarge-with-click-to-pin
behavior. This spec defines it once (§5.4); it's small (~50 lines of JS)
and reusable for future widgets.

---

## 1. Purpose

A dot plot panel on page 16 cross-species, showing genome-wide synteny
between *C. gariepinus* and *C. macrocephalus* haplotypes. Two display
modes:

1. **Mini panel** (~280 × 280 px) inline below the cs_breakpoints
   focus area, showing the full Gar × Mac dot plot at low detail.
2. **Enlarged popup** (~720 × 720 px) on hover, with full detail and
   tooltips. Click-to-pin: clicking the popup keeps it open until
   either (a) clicking outside the popup, or (b) clicking the × close
   button in the upper-right corner of the popup.

Each chromosome pair appears as a sub-rectangle within the global axes;
within each sub-rectangle, alignments are line segments (start_bp,
start_bp) → (end_bp, end_bp) with slope ±1 depending on strand. Color
encodes percent identity (continuous viridis-style ramp; the original
Perl script's 7-stop palette is genuinely unpleasant to look at and we
do not port it).

---

## 2. Why this matters

The existing macro-synteny panel (block list + permutation test) shows
*statistical* synteny structure but not *geometric* synteny. A dot plot
makes visible at a glance:

- **Inversions**: anti-diagonal segments inside an otherwise diagonal
  chrom pair
- **Translocations**: off-diagonal blocks (alignments to a non-homologous
  chrom)
- **Fissions / fusions**: a Mac chrom that maps to two Gar chroms (or
  vice versa) — visible as parallel diagonals across two chrom-pair
  boxes
- **Tandem duplications**: short parallel diagonals adjacent to the main
  one

These are exactly what the manuscript needs to display. The block list
is the data; the dot plot is the figure.

---

## 3. What this widget DOES / does NOT claim

### DOES

- Geometric visualization of synteny / inversions / translocations across
  the two haplotypes
- At-a-glance identification of chrom-pair relationships (1-to-1, 1-to-2
  fission, 2-to-1 fusion)
- Per-block tooltip with strand, % identity, block size, mapping quality

### Does NOT

- Phylogenetic polarity ("which is ancestral") — that's a separate spec
  (`SPEC_phylogenetic_tree_integration.md`)
- Resolution beyond the alignment tool's capability — wfmash with
  default block size (~5 kb) won't show small inversions; the multi-
  resolution mashmap mode (§5.2) is the answer
- Functional annotation overlay — gene density / TE density / centromere
  positions are separate tracks; could be future overlays but not v1

---

## 4. Vocabulary discipline

| Allowed | Discouraged |
|---|---|
| "Anti-diagonal segment indicates an inversion in this chromosome pair" | "Confirms the inversion" |
| "Block sizes ≥ 50 kb at PI ≥ 90% (mashmap-50k-pi90)" | "high-quality alignment" without parameters |
| "Mac_27 maps to LG28 with one large inversion (anti-diagonal at 14.1–16.5 Mb)" | "the inversion is between …" |

Always state alignment parameters in figure captions: tool, segment-length
threshold, percent-identity threshold, one-to-one filter status.

---

## 5. Architecture

### 5.1 Data sources — hybrid (option 3 from Quentin's design call)

Two parallel data sources, atlas chooses the best available:

#### Source A — wfmash synteny_blocks (manuscript-grade primary)

- Already exists in `cs_breakpoints_v1.json` v2
- Produced by wfmash with whatever parameters the user chose for
  breakpoint detection (typically segment ≥ 50 kb, PI ≥ 95%)
- One run, one parameter set, deterministic
- This is the **default day-1 source** — zero new compute

#### Source B — mashmap multi-resolution (`.out` files)

A new pre-compute pipeline runs mashmap **once** in `--one-to-one` mode
across the full Gar × Mac genome × 3 PI levels × 4 segment lengths,
producing 12 `.out` files. mashmap is fast (10s of seconds at this
scale) so this is one-shot pre-compute, not on-demand.

Recommended parameter grid:

| Run | `-s` (segment) | `--pi` | Filename |
|---|---|---|---|
| 1 | 100 kb | 0.85 | `mashmap_s100k_pi85.out` |
| 2 | 100 kb | 0.95 | `mashmap_s100k_pi95.out` |
| 3 | 500 kb | 0.85 | `mashmap_s500k_pi85.out` |
| 4 | 500 kb | 0.95 | `mashmap_s500k_pi95.out` |
| 5 | 1 Mb   | 0.85 | `mashmap_s1M_pi85.out`   |
| 6 | 1 Mb   | 0.95 | `mashmap_s1M_pi95.out`   |
| 7 | 2 Mb   | 0.85 | `mashmap_s2M_pi85.out`   |
| 8 | 2 Mb   | 0.95 | `mashmap_s2M_pi95.out`   |

(8 runs, not 12 — the third PI level is absorbed by the 0.85 floor since
mashmap's default lower bound is already there. If a third PI tier is
desired, add `--pi 0.99` for the strict tier across all four segment
sizes.)

All runs use `--one-to-one` to suppress paralogy clutter. mashmap output
is the standard whitespace-separated format (the Perl script's regex
documents it); a small Python script `STEP_DP01_aggregate_mashmap.py`
reads all 8 `.out` files and emits one combined JSON:

```jsonc
// File: dotplot_mashmap_v1.json
{
  "schema_version": 1,
  "tool": "mashmap v3.x",
  "generated_at": "2026-MM-DD",
  "species_query":  { "name": "C. gariepinus",   "haplotype": "fClaHyb_Gar_LG" },
  "species_target": { "name": "C. macrocephalus", "haplotype": "fClaHyb_Mac_LG" },
  "chrom_lengths_query":  { "LG28": 38500000, ... },
  "chrom_lengths_target": { "Mac_27": 41200000, ... },
  "resolutions": [
    {
      "label": "s100k_pi85", "segment_length_bp": 100000, "pi_min": 0.85,
      "n_alignments": 4823,
      "alignments": [
        { "q_chr": "LG28", "q_start": 14100000, "q_end": 14250000,
          "t_chr": "Mac_27", "t_start": 8200000, "t_end": 8350000,
          "strand": "+", "pi": 0.97 },
        ...
      ]
    },
    { "label": "s500k_pi85", ... },
    ...
  ]
}
```

- One JSON per haplotype pair, ~2-5 MB at typical scales
- Loaded via the same drag-drop pattern as other layers
- Stored on `state.dotplotMashmap`

#### Atlas mode resolution

The widget toolbar has a **resolution dropdown** with entries:
- `wfmash synteny (precise)` — Source A, default when cs_breakpoints v2
  is loaded
- `mashmap s100k pi85` / `s100k pi95` / `s500k pi85` / ... (8 entries) —
  Source B, populated from `dotplot_mashmap_v1.json`

User picks the resolution. Switching is instant (no fetch — both layers
are in state).

### 5.2 Layout algorithm — port `LayoutIDs` / `SpanXwY` from the Perl script

The single useful idea in the Perl reference is the **LIS-based fattest-
alignment diagonal layout**: chromosomes are reordered on both axes so
the heaviest synteny falls on the main diagonal. The algorithm (from
`LayoutIDs` + `SpanXwY` in the script):

1. For each query chrom Q, find its dominant target chrom T(Q) — the
   target chrom carrying the most aligned bp from Q
2. For each target chrom T, find its dominant query chrom Q(T) — symmetric
3. Greedily span: place query chrom with the largest size first; for
   each placed Q, place its dominant T(Q) on the y-axis at the same
   ordinal; recurse on T's dominant queries
4. Flip target chroms whose dominant alignment is reverse-strand (the
   `slope == -1` clause in `SpanXwY`)
5. Append unplaced chroms (no significant alignment) at the end

This produces a clean diagonal layout. **Port this to JS** (~80 lines)
in a new module helper `_dotplotLayoutChroms(blocks, chrom_lens_q,
chrom_lens_t)`. Tested against a synthetic fixture with 5 query × 5
target chroms in a known fission pattern.

### 5.3 Renderer — SVG (mini) and Canvas (enlarged)

#### Mini panel (~280 × 280 px) — SVG

- Light enough to render inline without performance issues
- ~5000 alignment segments at typical scales — fine for SVG
- One `<svg>` with grid lines at chrom boundaries, alignment segments
  as `<line>` elements colored by PI

#### Enlarged popup (~720 × 720 px) — Canvas

- More alignments potentially visible (zoom + smaller segments)
- Canvas faster than SVG above ~10k segments
- Tooltips computed on hover via the same coordinate transforms

Color ramp: continuous viridis-equivalent. Map PI ∈ [0.85, 1.0] →
colormap-position. Use the `var(--accent)` family for the high-PI end so
the plot reads on-theme. Suggested 5-stop ramp (low PI to high PI):

```
#1f4e79  →  #2c7a39  →  #cc8a00  →  #b54b1e  →  #c01a1a
```

Forward strand: solid line, opacity 0.85. Reverse strand: same color,
opacity 0.55 + 1px dotted. (The Perl script uses different colors for
forward vs reverse — color is overloaded; use opacity/dash instead.)

### 5.4 Hover / click behavior

State machine:

```
       ┌─────────┐                    ┌──────────┐                    ┌─────────┐
       │  MINI   │ ──── hover-in ───▶ │ ENLARGED │ ─── hover-out ───▶ │  MINI   │
       │ (inline)│                    │ (popup)  │                    │ (inline)│
       └─────────┘                    └──────────┘                    └─────────┘
                                           │
                                          click
                                           │
                                           ▼
                                      ┌──────────┐
                                      │  PINNED  │
                                      │ (popup)  │
                                      └──────────┘
                                           │
                              ┌────────────┴────────────┐
                              │                         │
                       click outside              click × button
                              │                         │
                              ▼                         ▼
                         ┌─────────┐               ┌─────────┐
                         │  MINI   │               │  MINI   │
                         └─────────┘               └─────────┘
```

Specifics:
- Hover-in: enlarge **fast** (no delay — Quentin's choice). CSS
  transition 80ms or instant (start with instant, add transition only
  if it feels jarring).
- Hover-out: shrink back to mini, immediate.
- Click on enlarged popup: pins it. State transitions to PINNED. Hover
  no longer controls — popup stays open.
- Pinned popup gets a × close button in the upper-right corner.
- Click outside the popup (anywhere on page background): closes popup,
  back to mini.
- × button click: closes popup, back to mini.

Implementation:
- `state.dotplotPin = { active: false, anchorEl: null }` — global flag
- `mousemove` hover handler reads `dotplotPin.active`; if pinned, hover
  events are no-ops
- Document-level `click` listener (added on pin, removed on unpin) that
  closes popup if `event.target` is outside the popup element

### 5.5 Page 16 integration

Add below `#csSynteny` (which already renders the existing macro-synteny
permutation panel):

```html
<div id="csDotplot" style="display: none; margin-top: 14px;"></div>
```

Renderer:

```javascript
function _renderCrossSpeciesDotplot() {
  const slot = document.getElementById('csDotplot');
  if (!slot) return;
  // Show only when at least one source is loaded
  const blocks   = _csGetSyntenyBlocks();              // wfmash source
  const mashmap  = state.dotplotMashmap;               // mashmap source
  if (!blocks && !mashmap) {
    slot.style.display = 'none';
    return;
  }
  slot.style.display = 'block';
  slot.innerHTML = _csBuildDotplotMini(blocks, mashmap);
  _wireCsDotplotHover(slot);                           // §5.4 state machine
}
```

The mini panel has its own header with the resolution dropdown:

```
┌────────────────────────────────────────────────────────────┐
│ Dot plot · gar × mac · resolution: [wfmash synteny ▾]      │   ← header
├────────────────────────────────────────────────────────────┤
│ ┌────────────────────────────────────────────────────────┐ │
│ │                                                        │ │
│ │       (mini dot plot, ~280×280 px)                     │ │
│ │                                                        │ │
│ └────────────────────────────────────────────────────────┘ │
│ hover to enlarge · click to pin · ×  close pinned          │
└────────────────────────────────────────────────────────────┘
```

Enlarged popup is positioned `position: fixed` at viewport center with
backdrop dim (`background: rgba(10, 12, 14, 0.5)` overlay). z-index
above the page-16 sticky header.

### 5.6 Per-breakpoint highlight (optional v2)

When the user has a breakpoint focused in the catalogue, the
corresponding region(s) on the dot plot can be highlighted:
- Yellow rectangle outlining the chrom-pair box that the breakpoint
  belongs to
- Red marker at the breakpoint's position within that box

Day-1 implementation: skip this — keeps the dot plot independent of
catalogue selection. Day-2 nice-to-have once the basic widget is solid.

### 5.7 Performance budget

- Mini render (5000 segments, SVG): ~30 ms
- Enlarged render (5000 segments, Canvas): ~15 ms
- Layout computation (one-time per data load): ~50 ms for 28 × 27 chroms
- Hover transition: instant (no animation cost)
- Pinned popup mouse-tracking for tooltip: throttled to 16 ms (60 fps)

All within budget for a snappy UI.

---

## 6. Module structure

```
phase_X_dotplot/
  README.md
  STEP_DP01_run_mashmap_grid.sh           # NEW — 8 mashmap calls
  STEP_DP02_aggregate_mashmap_to_json.py  # NEW — .out files → dotplot_mashmap_v1.json

atlas/
  Inversion_atlas.html                    # MODIFIED — csDotplot mount,
                                          # _renderCrossSpeciesDotplot,
                                          # _csBuildDotplotMini,
                                          # _csBuildDotplotEnlarged,
                                          # _wireCsDotplotHover,
                                          # _dotplotLayoutChroms
  tests/
    test_dotplot.html                     # NEW — layout + render + hover

specs_todo/
  SPEC_cross_species_dotplot.md           # this file
```

---

## 7. Tests

- `layout_places_dominant_chroms_on_diagonal`: synthetic 5×5 chrom
  fixture with known fission pattern, layout puts dominant target chroms
  at matching ordinals to their query chroms
- `layout_flips_reverse_dominant_pairs`: chrom pair where dominant
  alignment is reverse-strand → target placed flipped
- `mini_renders_wfmash_blocks`: 50 blocks → 50 line segments visible
- `mini_renders_mashmap_with_resolution_dropdown`: 8 resolutions in
  state, dropdown shows 9 entries (1 wfmash + 8 mashmap)
- `hover_in_enlarges_immediately`: mouseover triggers enlarge, no delay
- `hover_out_returns_to_mini`: mouseleave triggers shrink
- `click_pins_popup`: click while enlarged → state.dotplotPin.active = true
- `click_outside_closes_pinned`: pinned + click on backdrop → close
- `close_button_closes_pinned`: pinned + click × → close
- `pi_color_ramp_continuous`: PI 0.85 → low end of ramp, PI 1.0 → high
  end, monotonic
- `forward_strand_solid_reverse_strand_dotted`: same color but different
  visual treatment per strand

Target: 11 tests.

---

## 8. Open questions / explicitly deferred

- **Multi-species generalization** — currently Gar × Mac only. When
  species 3 lands (per `SPEC_OVERVIEW_multispecies_architecture.md`),
  the dot plot widget can pivot via a species-pair dropdown in the
  header. Forward-compat: the JSON shape already has
  `species_query` / `species_target` fields, so adding multiple JSONs
  (one per pair) is cheap. Not specced here, just flagged.
- **Per-breakpoint highlight** (§5.6) — deferred to v2.
- **Gene/TE/centromere overlays** — deferred. The atlas already has these
  as separate layers; combining them on the dot plot is a separate
  spec.
- **Zoom + pan inside enlarged popup** — v1 popup is fixed-scale. Pan/zoom
  could be added (Canvas-friendly), but most use cases for "look at one
  region in detail" are better served by the page 11 boundaries view at
  the relevant chromosome. Deferred.
- **Cross-link from the dot plot to the breakpoint catalogue** — clicking
  an alignment segment could jump to the nearest breakpoint in the
  catalogue. Day-2 nice-to-have.
- **`--one-to-one` parameter sensitivity** — mashmap's `--one-to-one`
  filter is aggressive and may drop legitimate paralogous alignments
  that are biologically interesting (segmental duplications, lineage
  duplications). For day 1, accept the aggressive filter (cleaner plot);
  flag in the README that users wanting to see paralogous structure
  should run mashmap separately without `--one-to-one`.

---

## 9. Summary

Dot plot widget on page 16 with two data sources (wfmash synteny_blocks
already in `cs_breakpoints_v1.json` v2; mashmap multi-resolution
`dotplot_mashmap_v1.json` from a new pre-compute step). Mini panel
inline (~280 px), enlarges fast on hover, click-to-pin with × close +
outside-click close. Color ramp on-theme (5-stop accent ramp, NOT the
Perl script's 7-stop palette). Forward/reverse strand encoded via
opacity+dash, not color. Layout uses the LIS-based fattest-alignment
diagonal algorithm ported from the Perl `LayoutIDs`/`SpanXwY` routines.

End of spec.
