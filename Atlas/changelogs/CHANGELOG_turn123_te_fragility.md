# CHANGELOG — turn 123: comparative TE breakpoint fragility (page16b)

**Atlas Δ:** +218 lines (60,799 → 61,017)
**Tests:** 73 new (480/480 cumulative green across 8 test suites)
**No regressions** — turns 117/118/119/120/121/122 all still pass.

---

## What shipped

A new JSON layer **`comparative_te_breakpoint_fragility_v1`** plus a
per-species TE density strip at the bottom of page 16b's center column.
This implements **hypotheses 4 and 5** from the April analysis chat
(1b4d8e12):

> **H4** — Some Cgar polymorphic inversions reuse ancient synteny boundaries.
> **H5** — Some homologous old-species regions are predicted polymorphic hotspots.

Both hypotheses test **fragility**, not polymorphism. The vocabulary
discipline established in the SPEC for this layer (and now baked into
the demo JSON's `vocabulary_discipline` block) is preserved on screen
in the legend disclaimer:

> *Fragility proxy only — not evidence of polymorphism in non-focal species.*

This matters because reviewer 2 will flag overclaiming. The atlas
shows only what the data supports: TE density at the homologous region
in each comparative species, in fold-enrichment vs that species'
chromosome-wide median.

### What the strip looks like

Below the existing **lineage distribution** table, a new
**comparative TE fragility — per species** block. One row per species
in the effective species list, columns:

| Species | fold vs chrom bg | %ile | density |
|---|---|---|---|
| ★ Clarias gariepinus | 1.72× *(elevated)* | 91%ile | 58.7% |
| ★ Clarias macrocephalus | 1.92× *(elevated)* | 96%ile | 61.2% |
| Clarias fuscus | 2.31× *(elevated)* | 98%ile | 68.1% |
| Channallabes apus | 1.11× *(neutral)* | 62%ile | 42.2% |
| Pangasianodon hypophthalmus | 3.13× *(strong)* | 99%ile | 75.5% |
| Hemibagrus wyckioides | 0.94× *(depleted)* | 38%ile | 28.8% |
| ... | no data | | |

Color bands on the fold-enrichment chip:

| Band | Range | Meaning |
|---|---|---|
| **strong** (red) | ≥ 2.5× | Candidate hotspot |
| **elevated** (orange) | 1.5–2.5× | Fragility-suggestive |
| **neutral** (grey) | 1.0–1.5× | At baseline |
| **depleted** (light) | < 1.0× | Below chromosome median |

### Why this matters for the manuscript

The **architecture-class auto-suggest** chip (Class A–F, turn 121) is
already on the page. When the lineage distribution shows ≥3 species
sharing the boundary, it suggests **Class E (recurrent rearrangement
hotspot)**. The TE fragility strip is the *next layer of evidence* for
that classification: if Class E is suggested AND multiple species show
elevated TE density at the homologous region, the **mechanism story
gets stronger** — repeat-mediated chromosome evolution at a reused
fragile substrate.

This is the multi-evidence convergence the Nature Comms reviewers will
look for. The atlas now displays it on a single page per breakpoint.

### Layer schema

```jsonc
{
  "tool": "comparative_te_breakpoint_fragility_v1",
  "schema_version": 1,
  "params": { "focal_radius_bp": 1100000, "te_class": "all_TE", ... },
  "per_breakpoint_per_species": [
    {
      "bp_id": "csbp_lg28_15_1Mb",
      "species": "Cmac",
      "gar_chr": "C_gar_LG28",
      "gar_pos_bp": 15115000,
      "homologous_chrom": "C_mac_LG01",
      "focal_lo_bp": 9100000, "focal_hi_bp": 11300000,
      "focal_te_density": 0.612,
      "bg_te_density_chrom": 0.318,
      "fold_enrichment": 1.92,
      "percentile": 96
    },
    ...
  ]
}
```

Tool name accepts both `comparative_te_breakpoint_fragility_v1` (full)
and `te_fragility_v1` (legacy alias) so pipeline scripts can use
either.

### Lookup matcher

`_msGetTEFragilityForBreakpoint(bp, speciesId)` matches in this order:
1. Exact `bp_id` match
2. `gar_chr` + `gar_pos_bp` within 100 kb of `bp.gar_pos_mb`

Loose enough to survive breakpoint coordinate refinement between
pipeline runs, strict enough to never match the wrong breakpoint on
the same chromosome. Same matcher pattern as `_msGetLineageDistribution`
and `_msGetDxyForBreakpoint`.

### Wiring

- `_classifyJSONKind` recognizes `comp_te_fragility`
- `loadMultipleJSONs` dispatches to `_storeCompTEFragility` and
  re-renders page16b if active
- `_restoreCompTEFragility` runs on page load alongside the other
  multi-species layers

### State slot

| Slot | Purpose | Persisted |
|---|---|---|
| `state.teFragility` | parsed `comparative_te_breakpoint_fragility_v1.json` | localStorage `inversion_atlas.teFragility.v1` |

### Demo file

`comparative_te_breakpoint_fragility_v1.demo.json` — six species
(Cgar, Cmac, Cfus, Capus, Phyp, Hwyc) covering the LG28 subtelomeric
breakpoint with values that exercise all four fold bands:
- Cgar 1.72×, Cmac 1.92×, Cfus 2.31× → elevated trio (focal + sister
  context — the "boundary substrate is shared" pattern)
- Phyp 3.13× → strong (deep mid-outgroup hotspot — supports recurrent
  reuse claim)
- Capus 1.11× → neutral (boundary present in synteny but not TE-rich)
- Hwyc 0.94× → depleted (homologous region is *not* fragile in this
  bagrid)

The full demo set is now: `synteny_multispecies_v1.demo.json` +
`phylo_tree_v1.demo.json` + `dxy_per_inversion_v1.demo.json` +
`comparative_te_breakpoint_fragility_v1.demo.json`. Drop all four
into the atlas → page 16b shows the full classification cockpit
(architecture chip + age model chip + lineage table + dXY note +
TE fragility strip) end-to-end for the LG28 subtelomeric breakpoint.

---

## What's left for next turn

- **TSV export columns.** `_msBuildClassificationTSV` already writes
  one row per cs_breakpoint with architecture + age + lineage columns.
  Could add per-species TE fold columns next.
- **busco_anchors_v1 layer.** BUSCO single-copy ortholog ticks on the
  page-16/16b dotplot for chromosome-homology backbone. Useful for
  deep-divergence species (Phyp, Tfulv, Smer, Tros) where wfmash
  doesn't bridge but BUSCO does.
- **Tree internal-node click → MRCA comparison.** Currently leaves-only
  on the phylo tree.
- **Hypothesis registry page wiring.** The big spec-todo from turn 120.
- **Fold/percentile chip on the right column dotplot caption** so the
  active-species detail view shows the TE fragility number inline.

---

## Files changed

```
Inversion_atlas.html                          +218 lines
test_turn123_te_fragility.js                  NEW (73 tests)
comparative_te_breakpoint_fragility_v1.demo.json  NEW (demo)
CHANGELOG_turn123_te_fragility.md             THIS FILE
```
