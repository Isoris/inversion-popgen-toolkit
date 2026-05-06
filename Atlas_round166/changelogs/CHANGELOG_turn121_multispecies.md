# CHANGELOG — turn 121: multi-species classification cockpit (page16b)

**Atlas Δ:** +1,068 lines  ·  **Tests:** 110 new (328/328 cumulative green)
**No regressions** — all 218 prior tests (turn117/117b/118/119/120) still pass.

---

## What shipped

A new synthesis-stage page **`page16b — multi-species`** that sits next
to page 16 cross-species. Its job: **place each Cgar↔Cmac breakpoint on
the catfish phylogeny, and auto-classify it** using the two-layer
schema Quentin codified in turn 113-tail (recorded on page 5 help).

The motivation, in one sentence: *the comparative species are not
inversion-detection targets — they are the evidence layer that
activates Quentin's age-model and architecture-class classification
chips for each cs_breakpoint*.

### Three-column layout

| Column | What's there | Reusable? |
|---|---|---|
| **Left (~26%)** | Clickable phylo tree (vertical cladogram, custom SVG, no library) | NEW |
| **Center (~42%)** | Architecture-class chip + auto-suggest rationale, lineage distribution table | NEW |
| **Right (~32%)** | Per-species detail: Cgar↔active species dotplot + boundary status | reuses `popgenDotplot.makeDotplotPanel` |

### Architecture-class auto-suggest

Reads `lineage_distribution` from the new `synteny_multispecies_v1`
layer and routes to one of A–F:

- **Class A** (simple inversion) — boundary in focal species only
- **Class B** (synteny-boundary) — exactly 2 species share the boundary
- **Class C** (fusion/fission-associated) — sister species shows fission/fusion
- **Class D** (terminal translocation) — derived from event_type when no lineage data
- **Class E** (recurrent rearrangement hotspot) — ≥3 species share — strongest mechanism class
- **Class F** (ambiguous) — fallback

When `lineage_distribution` is absent (no `synteny_multispecies_v1` JSON
loaded), the auto-suggest falls back to event-type heuristics with
`confidence: low` and a hint about what to load.

### Three layers consumed

1. **`synteny_multispecies_v1.json`** (NEW) — lineage_distribution +
   per-species synteny_blocks. Tool name accepted as either
   `synteny_multispecies_v1` or `catfish_synteny_toolkit_v1` (matches
   the toolkit's emit name).
2. **`phylo_tree_v1.json`** (NEW) — Newick + species_set + provenance.
3. **State from page 16** — active breakpoint via `state._crossSpeciesUI.active_id`.

### Graceful degradation hierarchy

1. **No active breakpoint** → empty state with "switch to page 13" link.
2. **Active breakpoint, no multi-species data** → fallback 9-species
   reference tree from `_MS_DEFAULT_SPECIES`, classification chip with
   low-confidence heuristic, lineage table shows "load JSON" hint.
3. **Active breakpoint + synteny_multispecies loaded** → full
   classification, lineage table populated.
4. **+ phylo_tree_v1 loaded** → tree on the left uses the published
   Newick instead of the reference flat tree.

### State slots added

| Slot | Purpose | Persisted |
|---|---|---|
| `state.syntenyMultispecies` | parsed `synteny_multispecies_v1.json` | localStorage `inversion_atlas.syntenyMultispecies.v1` |
| `state.phyloTree` | parsed `phylo_tree_v1.json` | localStorage `inversion_atlas.phyloTree.v1` |
| `state._multiSpeciesUI.active_species` | clicked species in tree | localStorage `inversion_atlas.multiSpeciesUI.v1` |

### JSON dispatcher

Two new detector + store + persist + restore + clear quintets, wired into:
- `_classifyJSONKind` (recognizes new tools)
- `loadMultipleJSONs` dispatch chain (stores + re-renders page16b on load)
- Restore-on-page-load chain (alongside `_restoreCrossSpecies`,
  `_restoreDotplotMashmap`, `_restoreNcRNADensity`)

### Demo files

- `synteny_multispecies_v1.demo.json` — 4 synteny blocks (Cgar↔Cmac,
  Cgar↔Cfus, Cgar↔Capus, Cgar↔Phyp) + 1 lineage distribution covering
  the LG28 subtelomeric breakpoint (10 species, hits Class E with 3
  species `boundary_present` + 2 `fission_in_target`). Drag into the
  atlas to see page 16b populate end-to-end.
- `phylo_tree_v1.demo.json` — Newick from the F1 manuscript
  Supplementary Data S4 lineage with 10 species, calibrated against
  TimeTree5.

### Cross-page navigation

- Page 16 → page 16b: the existing cs-breakpoint selection on page 16
  is honored automatically (the active_id state is shared).
- Page 16b active-header has a "view on page 13" link that switches
  back to page 16.
- "ⓘ classification framework" link in page 16b subtitle jumps to
  page 5 help anchor `#bp-classification-section` (same pattern as
  page 16).

---

## What was deferred to a future turn

- **Age-model auto-suggest** (YOUNG-POP / OLD-POLY / OLD-BP-YOUNG-INV /
  LINEAGE-KARYO / MULTI-AGE-HOTSPOT) — needs the `dxy_per_inversion_v1`
  layer + multi-species shared-breakpoint counter. The classification
  framework already shows a chip slot for this; the helper is the
  natural turn-122 follow-up.
- **TE fragility focal-vs-bg panel** in the right column — would reuse
  `popgenFocalVsBg.makeFocalVsBgPanel` once the
  `comparative_te_breakpoint_fragility_v1` layer is loaded.
- **Manuscript export** — TSV with one row per cs_breakpoint and both
  classification labels (`architecture_class`, `age_model`) + supporting
  evidence counts. Trivial helper once both auto-suggests are in.
- **Manual classification override** — editor at the top of the center
  column to set architecture and age model manually + save to
  `state.classifications[]`. ~150 lines.
- **Internal-node click on tree** — currently leaves-only. Internal
  nodes could show "comparison vs MRCA" view; deferred until needed.

---

## Files changed

```
Inversion_atlas.html          +1,068 lines (59,033 → 60,101)
test_turn121_multispecies.js  NEW (110 tests)
synteny_multispecies_v1.demo.json  NEW (demo)
phylo_tree_v1.demo.json            NEW (demo)
CHANGELOG_turn121_multispecies.md  THIS FILE
```

---

## Operating notes

The default `_MS_DEFAULT_SPECIES` list mirrors the 9-species
recommendation Quentin's prior chat work converged on (the 17-species
TE/BUSCO dataset filtered down to the species that are
phylogenetically informative for the Cgar↔Cmac breakpoint atlas, not
duplicates and not too-distant-to-bridge). It is **only** a fallback
for navigation when no real `phylo_tree_v1.json` is loaded — the real
tree from the F1 manuscript Supplementary Data S4 should drop in via
the file picker and replace it.

The `_msGetLineageDistribution` matcher tries (in order): exact `bp_id`,
then `gar_chr` + position-within-100kb. Loose enough to handle
breakpoint coordinate refinement between pipeline runs, strict enough
to never match the wrong breakpoint on the same chromosome.
