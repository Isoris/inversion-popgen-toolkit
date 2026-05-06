# SPEC patch — UpSet panel redirect (step 5)

**Status:** SHIPPED — landed alongside the SV evidence page (see
`tests/test_turn145_server_unify.js` and the SV-evidence track handoffs in
`handoffs/archive/2026-05-05_turn128-150/`).
**Date:** 2026-05-04
**Patches:** `specs_todo/SPEC_sv_evidence_page.md` §3.6 "UpSet" subsection and §2.2 schema (adds new layer)
**Origin:** Quentin's step-4 chat turn — the v0.1 spec had UpSet showing co-occurrence of which SVs share carriers; revised vision is more useful.

## What changed

The original spec §3.6 said the UpSet panel rows were `sv_id`s and the bars showed how many samples carried each SV-combination — a "which SVs travel together" view.

The redirected design uses UpSet as an **evidence-pattern selector** for ONE candidate:

- **Rows = per-SV-evidence types** (not SV ids).
- **Top bars = how many samples (fish) carry that exact evidence combination.**
- **Click a bar** → the bar's set of samples becomes the "selected_samples" — the locus track and main table both highlight only those samples' SVs (translucent intervals + dim non-selected glyphs).

## Worked example (from Quentin)

For candidate `INV_LG28_003`:

```
Rows (8 evidence types):
  left_SA          split-read clipping at left boundary
  right_SA         split-read clipping at right boundary
  left_PE          paired-end discordant orientation, left
  right_PE         paired-end discordant orientation, right
  Manta_INV_GT     Manta INV genotype call passes
  DELLY_INV_GT     DELLY2 INV genotype call passes
  MAPQ0_left       low-MAPQ region at left boundary
  MAPQ0_right     low-MAPQ region at right boundary

Top bars (per evidence-combination):
  42 fish carry [left_SA + right_SA + PE]            ← strongest support
  21 fish carry [MAPQ0 only]                          ← suspicious — likely repeat artefact
  17 fish carry [Manta_INV_GT + DELLY_INV_GT]         ← caller-agreement only
   8 fish carry [left breakpoint only]                ← single-sided evidence
```

Each top-bar maps to a column in the dot-matrix below the bars; the column's filled dots tell you which evidence rows are part of that combination. Standard UpSet plot semantics, just with our specific row vocabulary.

## Why this is a better design

The original "which SVs share carriers" framing was already partly answered by the `boundary_summary` table. The new framing answers a clinically important question: **for the same candidate, which evidence layers agree, and which fish does each agreement pattern correspond to?** That's the question that decides whether a candidate is real and whether its breakpoints are well-localised.

The 21-fish "MAPQ0 only" bar in the example is the perfect illustration: those samples have **no genuine evidence**, only a low-MAPQ region near the boundary. They're the samples that should NOT be counted as carriers. Selecting that bar lets the user inspect them and decide whether they're false positives.

## Schema extension

The current `sv_genotype_counts_v1` JSON does not carry per-sample-per-evidence-type information. The spec needs a new sibling layer:

### New layer: `sv_evidence_combinations_v1`

```
json/sv_evidence_combinations/<cid>.json
```

```json
{
  "format_version": "sv_evidence_combinations_v1",
  "candidate_id": "cand_LG28_15Mb",
  "n_samples_total": 226,
  "evidence_types": [
    { "id": "left_SA",       "label": "Left split-read",  "side": "left",  "kind": "SA" },
    { "id": "right_SA",      "label": "Right split-read", "side": "right", "kind": "SA" },
    { "id": "left_PE",       "label": "Left PE",          "side": "left",  "kind": "PE" },
    { "id": "right_PE",      "label": "Right PE",         "side": "right", "kind": "PE" },
    { "id": "Manta_INV_GT",  "label": "Manta INV",        "side": null,    "kind": "caller" },
    { "id": "DELLY_INV_GT",  "label": "DELLY INV",        "side": null,    "kind": "caller" },
    { "id": "MAPQ0_left",    "label": "Left MAPQ0",       "side": "left",  "kind": "mapq0" },
    { "id": "MAPQ0_right",   "label": "Right MAPQ0",      "side": "right", "kind": "mapq0" }
  ],
  "combinations": [
    {
      "members":            ["left_SA", "right_SA", "left_PE", "right_PE"],
      "intersection_size":  42,
      "samples":            ["FL01_001", "FL01_002", "..." ]
    },
    {
      "members":            ["MAPQ0_left", "MAPQ0_right"],
      "intersection_size":  21,
      "samples":            ["FL01_120", "..." ]
    },
    {
      "members":            ["Manta_INV_GT", "DELLY_INV_GT"],
      "intersection_size":  17,
      "samples":            ["FL02_044", "..." ]
    },
    {
      "members":            ["left_SA"],
      "intersection_size":  8,
      "samples":            ["FL03_005", "..." ]
    }
  ],
  "per_evidence_totals": {
    "left_SA":       { "n_samples": 67 },
    "right_SA":      { "n_samples": 65 },
    "left_PE":       { "n_samples": 49 },
    "right_PE":      { "n_samples": 51 },
    "Manta_INV_GT":  { "n_samples": 73 },
    "DELLY_INV_GT":  { "n_samples": 69 },
    "MAPQ0_left":    { "n_samples": 28 },
    "MAPQ0_right":   { "n_samples": 31 }
  }
}
```

### Notes on the schema

- `evidence_types` is the **canonical row order** for the UpSet — display tier-by-tier (split-read first, then PE, then caller calls, then MAPQ0). Producer-side decision.
- `combinations` is **sorted by intersection_size descending**, top N (≤ 30) for the panel. The producer pre-bins; atlas just reads.
- `samples[]` per combination is the list the click-handler uses to set `selected_samples`. These ids must match the karyotype-locked sample ids in `state.candidateState[cid].locked_labels`.
- `per_evidence_totals` powers the right-hand "Set size" mini-bars beside each row (standard UpSet UI element).

### Producer side

A new `STEP_SV_EVID_COMB` script (in MODULE_5A2 or MODULE_4) builds this layer from:

- DELLY2 + Manta VCFs → `Manta_INV_GT`, `DELLY_INV_GT` per-sample
- BAM-evidence track JSON (S7 / `phase_8_comparative_breakpoint_fragility`) → `*_SA`, `*_PE`, `MAPQ0_*` per-sample at the candidate's boundaries
- Outer join on sample_id; one row per sample with 0/1 per evidence type; group-by-pattern; write the JSON.

This producer is **separate** from `STEP_SV_GT_AGG` (which builds `sv_genotype_counts_v1`). The two layers can ship independently.

## Atlas-side: handshake to track + table

When the user clicks a bar in the UpSet:

1. `_state.selectedSamples` ← `combinations[i].samples` (Set)
2. **Locus**: SV glyphs for samples NOT in selected_samples render at 0.18 opacity (same dim treatment as filtered-out SVs). Glyphs FOR samples in selected_samples get a translucent band background (alpha-25%) spanning the SV's interval, in the per-SV-type colour. Spec quote: *"its interval is alpha background highlighted"*.
3. **Table**: the SV-table is filtered to rows that have ≥ 1 sample in selected_samples carrying the SV. The genotype-count columns are recomputed within selected_samples (so `H1/H1 (n=42) AA AB BB miss` becomes the within-selection view).
4. **Right-rail boundary summary**: optionally re-counted within selection (probably yes — clarify in step 5 build chat).
5. **A "Clear selection" pill** appears at the top of the page when selected_samples is non-empty.

## Acceptance for step 5

- [ ] UpSet panel renders with rows = `evidence_types`, top bars = `combinations[].intersection_size`, dot-matrix below the bars.
- [ ] Per-evidence-type set-size mini-bars on the right of the row labels (from `per_evidence_totals`).
- [ ] Click a top bar → `selected_samples` populated; locus + table re-render.
- [ ] Click bar again → deselect (clears selection).
- [ ] "Clear selection" pill appears when active; one click clears everything.
- [ ] Panel handles missing layer gracefully ("Drop sv_evidence_combinations/<cid>.json to populate") — fits the local-first workflow.
- [ ] Schema-violating JSON rejected with a clear error.

## Region-select mode (step 4.5, between step 4 and step 5)

Quentin's other request — "for the selector maybe when we select the track its not a zoom but actually yes its a zoom but we can also have a mode that selects and dynamically updates the table and upsetR".

Implementation: a toggle button in the toolbar `Mode: zoom | select`. In `select` mode, click-drag on the locus draws a translucent rectangle. On mouseup, the rectangle's [bp_start, bp_end] becomes a sticky region filter; the table + locus + UpSet (when shipped) re-filter to SVs whose `position_bp` falls inside the region. Esc clears the region. Same pattern as page 11's manual-override pen-tool.

State addition:
```
_state.selectMode = 'zoom' | 'select'
_state.selection  = { startBp, endBp } | null
```

This is a small step — defer until after step 5 ships.
