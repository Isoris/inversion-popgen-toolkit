# SPEC — SV Evidence & Boundaries page

**Status:** spec only. Code-ready for the next chat.
**Authoritative inputs:** three mockups (`1000100491.png`, `1000100494.png`, `1000100496.png`) showing the intended page in light + dark mode, plus an alternate layout variant.
**Module name:** `atlas_sv_evidence`
**Page id:** `page_sv_evidence` (mounts via `<div data-page="page_sv_evidence">`)
**Tab label:** `SV Evidence` (or `SVs` in the compact variant from mockup 2)
**Author hand-off:** Quentin Andres · MS_Inversions_North_african_catfish

---

## 1. What this page is for

Surface, **per candidate**, every SV call (DELLY2 + Manta) within ±500 kb of
the candidate's left and right boundary, scored against karyotype-group
genotype counts (H1/H1, H1/H2, H2/H2), classified into a small set of
pattern labels (canonical breakpoint marker / dominant / het-specific /
sub-haplotype / internal linked / uninformative).

The page answers four questions in one screen:

1. *Where do SV calls cluster relative to the candidate's edges?* → locus track view
2. *How many calls are there, by SV type, per boundary?* → BOUNDARY SUMMARY (right rail)
3. *Which specific SVs distinguish karyotype groups?* → main SV table with per-group counts + Fisher OR + FDR
4. *Do those SVs co-occur in coherent combinations?* → UpSet panel + sample × SV heatmap

The page is **read-only**. It does not edit boundaries (that's page 11),
does not edit karyotype groups (that's the scrubber), does not edit SV
calls (those come from upstream MODULE_5A2). It joins three things the
atlas already has access to and makes them legible.

## 2. Inputs

### 2.1 Existing layers (already produced by the cluster pipeline)

| layer name | path convention | producer | required |
|---|---|---|---|
| `sv_evidence`            | `json/sv_evidence/<candidate_id>.json` | `MODULE_5A2 → export_sv_to_json` | **yes** |
| candidate registry block | embedded in the global candidate JSON, key `q3_*` | STEP_03B bridge   | **yes** (for boundary positions) |
| karyotype groups         | already in `state.candidateState[cid].locked_labels` (HOMO_1 / HET / HOMO_2) | scrubber locks | **yes** |
| dosage matrix            | `json/<chrom>.json` (existing layer) | dosage emitter | **yes** (for the dosage heatmap track) |
| local PCA traces         | already in `state.candidatePCA[cid]` | `STEP_R32`/D17  | **yes** (for the PCA-axis-1 track) |
| genes (GFF)              | already in `state.geneTrack` | GFF loader | optional |
| repeat density           | already in `state.repeatDensity` | `STEP_RD_TE_density_full.R` | optional |

### 2.2 New layer this page requires

The existing `sv_evidence` JSON has **per-call rows** (sv_id, chrom, pos,
type, qual, callers). It does **not** have **per-call × per-karyotype-group
genotype counts**, which is what the main table renders.

We need a new sibling layer `sv_genotype_counts_v1`:

```
json/sv_genotype_counts/<candidate_id>.json
```

Schema:

```json
{
  "format_version": "sv_genotype_counts_v1",
  "candidate_id": "cand_LG28_15Mb",
  "chrom": "C_gar_LG28",
  "boundary_left_bp": 15142000,
  "boundary_right_bp": 18124000,
  "zone_definitions_bp": {
    "left_flank":     [14140000, 14640000],
    "left_boundary":  [14640000, 15640000],
    "inversion_body": [15640000, 17620000],
    "right_boundary": [17620000, 18620000],
    "right_flank":    [18620000, 19120000]
  },
  "groups_used": {
    "H1/H1":  { "label": "H1/H1", "n": 61,  "members": ["..."] },
    "H1/H2":  { "label": "H1/H2", "n": 103, "members": ["..."] },
    "H2/H2":  { "label": "H2/H2", "n": 62,  "members": ["..."] }
  },
  "sv_calls": [
    {
      "sv_id":       "SV001",
      "sv_type":     "BND",
      "chrom":       "C_gar_LG28",
      "position_bp": 15142380,
      "end_bp":      18118200,
      "zone":        "left_boundary",
      "distance_to_edge_bp": 380,
      "n_samples_with_call": 226,
      "quality":     "PASS",
      "callers":     ["delly2", "manta"],
      "genotype_counts": {
        "H1/H1": { "AA": 55, "AB": 5,  "BB": 1,  "miss": 0 },
        "H1/H2": { "AA": 10, "AB": 81, "BB": 8,  "miss": 4 },
        "H2/H2": { "AA": 0,  "AB": 6,  "BB": 55, "miss": 1 }
      },
      "fisher": {
        "comparison":   "H1/H1_vs_H2/H2",
        "odds_ratio":   42.1,
        "p_value":      1.2e-12,
        "fdr_bh":       2.1e-9
      },
      "pattern_label": "canonical_breakpoint_marker",
      "notes":         "BND pair → 18.12 Mb"
    }
  ],
  "boundary_summary": {
    "left": {
      "interval_bp":  [14640000, 15640000],
      "by_sv_type":   {
        "BND":   { "n_total": 24, "n_associated_fdr_lt_0_05": 9 },
        "INV":   { "n_total": 11, "n_associated_fdr_lt_0_05": 5 },
        "DEL":   { "n_total": 38, "n_associated_fdr_lt_0_05": 7 },
        "DUP":   { "n_total": 22, "n_associated_fdr_lt_0_05": 4 },
        "Other": { "n_total": 7,  "n_associated_fdr_lt_0_05": 0 }
      }
    },
    "right": {
      "interval_bp":  [17620000, 18620000],
      "by_sv_type": {
        "BND":   { "n_total": 27, "n_associated_fdr_lt_0_05": 11 },
        "INV":   { "n_total": 12, "n_associated_fdr_lt_0_05": 6 },
        "DEL":   { "n_total": 41, "n_associated_fdr_lt_0_05": 8 },
        "DUP":   { "n_total": 24, "n_associated_fdr_lt_0_05": 5 },
        "Other": { "n_total": 5,  "n_associated_fdr_lt_0_05": 0 }
      }
    }
  },
  "upset_top_combinations": [
    { "members": ["SV001"],          "intersection_size": 42 },
    { "members": ["SV001","SV002"],  "intersection_size": 28 },
    { "members": ["SV001","SV004"],  "intersection_size": 16 }
  ]
}
```

The new producer for this layer is **STEP_SV_GT_AGG** in MODULE_5A2 — it
joins the per-call SV TSV with the candidate karyotype-group membership
and runs Fisher + BH-FDR per call across the three pairwise comparisons,
keeping the maximum-magnitude one as `fisher.comparison`. The script does
not exist yet; spec it in the same chat that builds the page or stub it
with a small Python emitter that the next data-side person can replace.

### 2.3 Pattern-label vocabulary

These six labels appear in the right-rail legend (mockup 1, ASSOCIATION
panel) and as colored text in the table's `Pattern (label)` column.
The mapping rule is deterministic; render it in `docs/sv_pattern_labels.md`:

| label                          | when                                                                                            | colour |
|---|---|---|
| `canonical_breakpoint_marker`  | sv near boundary AND `OR ≥ 25` AND `FDR < 1e-6` AND H1/H1+H2/H2 fully separated by AA vs BB     | violet `#7d4cdb` |
| `dominant_presence_marker`     | sv near boundary AND `OR ≥ 5` AND `FDR < 0.01` AND non-zero `AB` row                              | amber  `#d08770` |
| `het_specific_marker`          | sv where H1/H2 row > H1/H1+H2/H2 rows for `AB`                                                   | teal   `#5e81ac` |
| `sub_haplotype_marker`         | sv associated only inside one homozygote row (e.g., H1/H1 only)                                  | blue   `#88c0d0` |
| `internal_linked_marker`       | sv inside `inversion_body` zone, `OR ≥ 2`                                                        | yellow `#ebcb8b` |
| `uninformative`                | `FDR ≥ 0.05`                                                                                     | grey   `#888888` |

Tie-breakers + edge cases live in the docs file. The atlas trusts the
producer's label and **does not recompute** — it only colour-codes.

## 3. Page anatomy

### 3.1 Layout (matches mockup 1, the canonical one)

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│ Catfish Inversion Atlas — tabs: Overview | Karyotype/Tiers | Boundaries | Dosage │
│                          | PCA | [SV Evidence ✓] | Coherence | Polarity | Sim Mat │
├────────────────┬───────────────────────────────────────────────┬─────────────────┤
│ LEFT RAIL      │ MAIN: locus track view                         │ RIGHT RAIL      │
│  (~210px)      │                                                 │  (~290px)       │
│                │  ┌─────────────────────────────────────────┐    │                 │
│ Candidate      │  │ Zone bar:                                │    │ BOUNDARY        │
│ ─ chrom        │  │ [L flank][Left bound][Inv body][R bound][R flank] │ │ SUMMARY         │
│ ─ region       │  ├─────────────────────────────────────────┤    │ (left)          │
│ ─ span         │  │ Coordinate axis: 14.5 .. 19.0 Mb         │    │ ┌─SV/N/Assoc─┐ │
│                │  ├─────────────────────────────────────────┤    │ │ BND 24 9   │ │
│ Karyotype      │  │ Dosage heatmap (cap 200 samples) ─ ─ ─ ─│    │ │ INV 11 5   │ │
│  H1/H1  n=61   │  ├─────────────────────────────────────────┤    │ │ DEL 38 7   │ │
│  H1/H2  n=103  │  │ PCA local: PC1/PC2/PC3 line traces       │    │ │ DUP 22 4   │ │
│  H2/H2  n=62   │  ├─────────────────────────────────────────┤    │ │ Other 7 0  │ │
│                │  │ SV calls (BND/INV/DEL/DUP/Other) icons   │    │ │ Total 102 25│ │
│ Boundaries     │  ├─────────────────────────────────────────┤    │ └────────────┘ │
│  Left  …       │  │ Genes (GFF strand arrows)                │    │                 │
│  Right …       │  ├─────────────────────────────────────────┤    │ BOUNDARY (right)│
│                │  │ Repeat density (mini bars)               │    │ ┌────────────┐ │
│ Zones          │  └─────────────────────────────────────────┘    │ │ BND 27 11  │ │
│  L flank …     │                                                  │ │ ...        │ │
│  L boundary …  │                                                  │ └────────────┘ │
│  Body …        │  ┌── SVs near boundaries (±500 kb) ──────────┐    │                 │
│  R boundary …  │  │ tabs: Summary | Genotype counts by group   │    │ LEGEND          │
│  R flank …     │  │      | UpSet | Sample × SV heatmap         │    │ ─ FDR colours   │
│                │  ├────────────────────────────────────────────┤    │ ─ Pattern labels│
│ FILTERS        │  │ sv_id type zone pos d_edge n AA AB BB ...  │    │                 │
│  SV type [All] │  │ SV001 BND  L     15142k +380 226 …         │    │ UPSET (top 12)  │
│  Quality PASS  │  │ SV002 INV  L     15138k -3100 …            │    │ ─ vertical bars │
│  Zone  [Bnd…]  │  │ ...                                         │    │ ─ dot matrix    │
│  Min samples 5 │  │ 25 rows, paginated, 148 total              │    │                 │
│  ☑ Show only   │  └────────────────────────────────────────────┘    │ SAMPLE × SV     │
│    associated  │                                                  │ HEATMAP         │
│  FDR thr 0.05  │                                                  │ (top assoc SVs) │
│                │                                                  │                 │
│  [Apply][Reset]│                                                  │ Sort by [H sys] │
└────────────────┴───────────────────────────────────────────────┴─────────────────┘
```

### 3.2 Top toolbar (above the locus view)

- `Locus: <chrom>: <left.left_flank>–<right.right_flank>  (<span_mb> Mb window)`
- `View presets: [Default ▾]` — dropdown with three named viewport zooms:
  - `Default` (left flank → right flank, ~5 Mb)
  - `Left boundary close-up` (boundary_left_bp ± 250 kb)
  - `Right boundary close-up` (boundary_right_bp ± 250 kb)
- Zoom in / zoom out / fit / reset buttons (icon row)
- A small `(i)` info button → opens a compact help popover quoting the
  spec's TL;DR

### 3.3 Left rail — candidate context + filters

Both sections render from existing atlas state (`state.candidateState[cid]`,
`state.candidateList`) — no new data needed for the rail itself.

**Candidate context** (read-only, no inputs):

- `Candidate:` bold, with `TIER 1/2/3/4` chip (from `final_classification`
  layer if present, else hidden)
- chromosome, region, span (Mb)
- karyotype group buttons (3 total, colour-tinted teal/amber/violet),
  each clickable to **filter the heatmap and table to one group**
- inferred boundaries (`Left` and `Right` in bp)
- zone bp ranges (5 rows, colour swatches matching the locus zone bar)

**Filter SVs** — applies live to the table and the locus SV-icon track:

- `SV type` dropdown: All | BND | INV | DEL | DUP | Other
- `Quality` dropdown: PASS | All
- `Zone` dropdown: All | Boundary ±500 kb | Inversion body | Flanks only
- `Min samples` slider 1..226, default 5
- Checkbox: `Show only associated` (FDR < threshold), default ON
- `FDR threshold` dropdown: 1e-6 | 1e-4 | 0.01 | **0.05** | 0.10
- Buttons: `Apply filters` / `Reset`

State stored in `state.svEvidenceFilters[candidate_id]`, persisted to
localStorage (mirror the page-11 pattern).

### 3.4 Main locus tracks (top half)

A single shared x-axis (bp coordinates) with these vertically stacked tracks:

1. **Zone bar** — five colour-coded segments labelled `Left flank`,
   `Left boundary`, `Inversion body`, `Right boundary`, `Right flank`.
   Boundary labels appear above the corresponding zone, in the zone's
   accent colour.
2. **Coordinate axis** — bp ticks every 0.5 Mb, dashed verticals at
   `boundary_left_bp` (blue) and `boundary_right_bp` (red).
3. **Dosage heatmap** — sample × window matrix, capped to 200 displayed
   samples, palette `RdBu_r` clipped at ±2. Reuse the existing dosage
   renderer from page 11; just remount in this page's container.
4. **PCA local** — three line traces (`PC1`, `PC2`, `PC3`) per window,
   y-axis ~[-0.15, +0.15]. Same data as page 11 PCA panel; remount.
5. **SV calls (all types)** — one row per SV-type lane, glyph per call:
   `▼ BND` (red), `▼ INV` (violet), `▼ DEL` (orange), `▼ DUP` (teal),
   `▲ Other` (grey). Click a glyph → highlights the corresponding row in
   the bottom table + scrolls it into view. Hover → mini tooltip with
   `sv_id`, position, distance to nearest edge, FDR.
6. **Genes (GFF)** — strand arrows. Existing renderer.
7. **Repeat density** — mini bar chart. Existing renderer.

The track order and the legend strip beneath the SV-calls track
(`SV types: ▼ BND  ▼ INV  ▼ DEL  ▼ DUP  ▲ Other`) match mockup 1 exactly.

### 3.5 SVs-near-boundaries table (bottom half, main column)

Tabs across the top of the table:

- **Summary table** (default) — one row per SV, full column set
- **Genotype counts by group** — wider view, three sub-blocks (one per
  group) with AA/AB/BB/miss columns
- **UpSet (combinations)** — same data as the right-rail UpSet, but
  larger and interactive
- **Sample × SV heatmap** — same as right-rail heatmap, larger

**Default columns** (Summary tab):

```
sv_id | SV type | Zone | Position (bp) | Dist. to edge (bp) | Samples (n)
| H1/H1 (n=61): AA AB BB miss
| H1/H2 (n=103): AA AB BB miss
| H2/H2 (n=62): AA AB BB miss
| H1/H1 vs H2/H2: OR (Fisher) | p value | FDR
| Pattern (label) | Notes
```

Behaviour:

- `sv_id` rendered as a small star icon when starred (★), bold when the
  SV's icon in the locus view is highlighted.
- Sortable by every column. Default sort: `FDR` ascending.
- Pagination at the bottom: `Rows per page [25 ▾]`, page indicator
  `1 .. n`, prev/next.
- Footer counters: `AA = 0/0 (ref/ref)   AB = 0/1 (het)   BB = 1/1 (alt/alt)
  OR = odds ratio` — verbatim from mockup 1.
- `Export table (TSV)` button (top-right of the table region) — emits the
  currently-filtered, currently-sorted rows.

### 3.6 Right rail — summary, legend, UpSet, heatmap

**Boundary summary** — two stacked tables (left and right), columns
`SV type | # SVs | # Associated (FDR < 0.05)`. Click a row to filter the
main table to that SV type for that boundary.

**Legend** — two sub-sections side by side:

- ASSOCIATION (FDR colours): six dots
  ```
  • < 1e-6     red
  • 1e-6 – 1e-4 orange
  • 1e-4 – 0.01 yellow
  • 0.01 – 0.05 green
  • > 0.05     grey
  ```
  These dots colour the FDR cell background in the main table.
- PATTERN LABELS — six coloured words matching the §2.3 vocabulary.

**UpSet (top 12 associated SVs)** — vertical bars (intersection sizes) +
dot-matrix below. Only renders SVs with FDR < threshold. Click a bar →
highlights those SVs in the main table.

**Sample × SV heatmap (top associated SVs)** — top 6 SVs by lowest FDR,
samples columns coloured by group bar at top (H1/H1 teal | H1/H2 amber |
H2/H2 violet). Cell colours: red high alt-dosage, blue low. Sort dropdown:
`Sort samples by [H system | family | dosage at SV001]`. Checkbox: `Group
by family` (default ON).

## 4. Wiring contract

### 4.1 Tab registration

```html
<!-- in #tabBar near the existing tab buttons -->
<button data-page="page_sv_evidence" data-stage="refinement"
        title="SV calls clustered around the candidate's boundaries, scored
               against karyotype groups (H1/H1, H1/H2, H2/H2). Per-call
               genotype counts, Fisher OR, FDR, pattern labels, UpSet
               combinations, sample × SV heatmap. Read-only.">
  <span class="num">5</span> SV evidence
</button>
```

Pick the position in the tab bar between page11 (boundaries) and page3
(catalogue) — that's the visual order in mockup 1's top tab strip.

Renumber the existing tab `<span class="num">N</span>` labels to keep the
sequence contiguous. The mockup numbering shows `4 boundaries / 5 SV
evidence / 6 catalogue`. Treat the renumbering as a separate trivial
patch (one search-replace in the tab bar markup) and don't touch tab
*ids* — only the visible numeric label.

### 4.2 Page mount

```html
<!-- in the page container area, after the existing pages -->
<div class="page" data-page="page_sv_evidence">
  <!-- atlas_sv_evidence.js renders into this container -->
  <div id="sv_evidence_root"></div>
</div>
```

### 4.3 JS module

New file: `js/atlas_sv_evidence.js`

Exposes:

```js
window.AtlasSVEvidence = {
  init(opts),               // called once at atlas startup
  loadCandidate(cid),       // fetches sv_genotype_counts/<cid>.json, renders
  refresh(),                // re-renders the active candidate using state.* updates
  setFilters(filters),      // programmatic filter override (URL-state)
  exportFilteredTSV(),
  _state,                   // for debugging
};
```

Lifecycle:

1. On `pageActivate` event for `page_sv_evidence`, the atlas calls
   `AtlasSVEvidence.loadCandidate(activeCandidateId)`.
2. If no candidate is active, the page renders an empty-state
   ("Promote a candidate on page 1 or page 3, then come back here.").
3. When the user changes karyotype locks on the scrubber, the dispatched
   `karyotype_locks_changed` event triggers `AtlasSVEvidence.refresh()`
   (which re-aligns groups by reading `state.candidateState[cid].locked_labels`).
4. The locus tracks reuse the existing dosage / PCA / GFF / repeat-density
   renderers — pass them the same `region_bp` they receive on page 11.

### 4.4 Empty / partial states

| condition | behaviour |
|---|---|
| no `sv_evidence` JSON loaded                | left-rail context renders; main panel shows `Load <code>sv_evidence/&lt;cid&gt;.json</code> to populate this page.` |
| `sv_evidence` loaded but no `sv_genotype_counts` | locus tracks render with SV glyphs; main table shows columns for `sv_id, type, zone, position, samples` only; the genotype/Fisher/FDR/pattern columns show `—` and a single info pill `genotype counts layer not loaded` |
| candidate has no karyotype locks            | "lock karyotype groups on the scrubber first" empty state in the right rail |
| 0 SV calls within ±500 kb of either boundary | locus shows tracks with no SV row; table shows "no SV calls within filter window"; right-rail boundary summary still renders zeros |

Never throw. Render partials.

## 5. Data flow details

### 5.1 Joining sv_evidence with karyotype groups (atlas-side)

The atlas already has `locked_labels: Array<{sample_id, label}>` per
candidate from the scrubber's K=3 PCA cluster.

```js
// pseudo
const labels = state.candidateState[cid].locked_labels;  // [{sid, "HOMO_1"}, ...]
const remap = { HOMO_1: "H1/H1", HET: "H1/H2", HOMO_2: "H2/H2" };
const groupOfSample = new Map(labels.map(r => [r.sample_id, remap[r.label]]));
```

If `sv_genotype_counts_v1` is loaded the page **uses its precomputed
counts and Fisher results verbatim** — it does not recompute. If
`sv_genotype_counts_v1` is absent, the page may fall back to a JS-side
recomputation from `sv_evidence + dosage matrix + groupOfSample` for the
top 25 SVs by quality only (limit hardcoded to keep the browser
responsive). The fallback is opt-in behind a checkbox `Recompute in
browser (slow)` to surface that this is not the canonical numbers.

### 5.2 FDR colouring rule

```js
function fdrColour(fdr) {
  if (fdr < 1e-6)   return '#bf616a';  // red
  if (fdr < 1e-4)   return '#d08770';  // orange
  if (fdr < 0.01)   return '#ebcb8b';  // yellow
  if (fdr < 0.05)   return '#a3be8c';  // green
  return '#888888';                    // grey
}
```

Apply to the FDR cell background only. Other rows stay unstyled.

### 5.3 Pattern label colouring

Apply to the `Pattern (label)` cell text colour, using the §2.3 colours.
Same in the right-rail legend.

## 6. Cross-references

- This page does **not** render boundary refinement — that's `page11`.
- This page does **not** call SV variants — those come from upstream
  MODULE_5A2 / DELLY2 / Manta.
- This page **uses** the karyotype locks from the scrubber — but doesn't
  edit them.
- The S7 spec (`specs_todo/from_turn129/S7_karyotype_breakpoint_internal_evidence.md`)
  defines the read-evidence layer that the table's `Notes` column may
  cite ("BND pair → 18.12 Mb"). When S7's atlas JSON lands, this page
  promotes the Notes column from free text into a hover-pop with the
  cluster-level evidence breakdown. Until then, Notes is whatever the
  upstream emitter writes.

## 7. Implementation order (next chat)

1. Create `js/atlas_sv_evidence.js` with the lifecycle + empty-state
   skeleton + filter rail. Mount the empty page. **Verify** the tab
   appears, the candidate context populates from existing state, and
   `loadCandidate(cid)` no-ops gracefully when the JSON is absent.
2. Implement `sv_genotype_counts_v1` loader + the **main SV table** with
   sort, filter, pagination, FDR colouring, pattern-label colouring,
   TSV export. **Verify** with a tiny stub JSON.
3. Implement the **locus track strip** by reusing existing renderers
   (dosage, PCA, GFF, repeat density). Add the new SV-icons track from
   scratch. **Verify** click-to-highlight handshake between SV icons
   and the table.
4. Implement the **right-rail boundary summary tables** + **legend** —
   pure presentational, no new data path.
5. Implement the **UpSet** panel (right rail compact + main-area large
   variant). Reuse `d3.upset` if already available; else a tiny custom
   bar+dot SVG.
6. Implement the **sample × SV heatmap** (right rail compact + main-area
   large). Reuse the dosage-heatmap renderer parameterised by
   row-source.
7. Add the **STEP_SV_GT_AGG** stub Python emitter under
   `_scripts/svgt/STEP_SV_GT_AGG_aggregate_genotype_counts.py` so the
   data side can replace it later. Stub fills the schema with random
   values from the loaded `sv_evidence` + a karyotype TSV. Document the
   real producer's TODO inside the file.

Each step ends with a screenshot match against the corresponding mockup
region and a tiny test fixture committed under `tests/sv_evidence/`.

## 8. Out of scope (explicitly)

- BAM-side read-evidence rendering (clipping, split reads, MAPQ0) —
  that's S7. This page may **link** into S7 once it ships, but does not
  re-render its content here.
- LD heatmap of SVs — already on a different page.
- Editing SV calls / merging callers — that's MODULE_5A2.
- Cross-species SV comparison — page 16.
- Breakpoint-architecture classification (A–G classes) — that's
  `phase_8_comparative_breakpoint_fragility`.
- Karyotype editing — scrubber.

## 9. Acceptance test

The page is "done" when, given a real candidate's
`sv_evidence/<cid>.json` + `sv_genotype_counts/<cid>.json` + the existing
loaded atlas state:

- The locus track strip matches mockup 1 visually within 5% pixel
  variance (track order, colour, glyph types, dashed boundary verticals).
- The main SV table renders ≥ 100 rows, sorts and paginates correctly,
  exports the filtered set as TSV in a single click.
- Clicking an SV glyph in the locus view highlights its table row and
  scrolls it into view. Clicking a row's `sv_id` highlights the glyph.
- The boundary summary right-rail tables match the candidate's actual
  per-zone counts (verified against the JSON).
- Toggling `Show only associated` removes all FDR ≥ 0.05 rows and
  recomputes the right-rail UpSet + heatmap.
- Switching candidate via the candidate selector triggers a re-render
  with the new JSON, no stale rows from the previous candidate.
- All four table tabs (Summary / Genotype counts by group / UpSet /
  Sample × SV heatmap) render without browser console errors.
- The page works in light AND dark mode (mockup 2 is the dark variant).

## 10. Unresolved questions for Quentin

These are decisions deferred to the build chat — flag them at the top of
the next session:

1. **Tab numbering renumber:** the mockup shows `4 boundaries / 5 SV
   evidence / 6 catalogue`. Confirm you want the renumber, or keep the
   tab as `5b SV evidence` to avoid touching the existing numeric labels.
2. **Pattern-label authority:** are the labels emitted by the data side
   (preferred) or computed atlas-side from FDR + zone + group counts
   (faster to ship, harder to keep consistent across sessions)?
3. **`sv_genotype_counts` producer:** is there an existing R / Python
   script in MODULE_5A2 that already emits per-call × group counts (in
   any format), or does this need a new STEP_SV_GT_AGG?
4. **Notes column source:** for now it's a free text field on the SV row.
   Once S7 lands, do you want Notes to be the S7 cluster summary, or stay
   as free text and surface S7 elsewhere?
5. **Recompute-in-browser fallback:** keep it (transparent to users that
   they're seeing approximate numbers) or remove it (force the data-side
   layer)?

---

**One-line summary for the build chat:** "Build a new tab `SV evidence`
mounting at `<div data-page="page_sv_evidence">`, JS module
`js/atlas_sv_evidence.js`, consuming the existing `sv_evidence` JSON
plus a new `sv_genotype_counts/<cid>.json` layer, rendering the locus
track strip + main SV table + right-rail boundary summaries / UpSet /
sample×SV heatmap exactly as in mockup `1000100494.png`."
