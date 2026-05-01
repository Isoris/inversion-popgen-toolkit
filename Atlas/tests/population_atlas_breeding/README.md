# Population Atlas — Breeding page (page 8) + Hatchery health page (page 9)

Two new tabs added to `Population_atlas.html`. Companion to the existing
6 pages (samples / QC / families / het / diversity / inbreeding) and the
about page.

After this addition the visible tab order is:
**1 samples · 2 QC · 3 families & clusters · 4 heterozygosity · 5 diversity · 6 inbreeding · 7 hatchery health · 8 breeding · 9 about**

The hatchery health page sits between Inbreeding and Breeding by design
so the visual flow reads: how-inbred (page 6) → how-is-the-stock-doing-overall (page 7) → which-fish-to-act-on (page 8).

## What it is

A short, committee-friendly list of **flagged broodstock samples**. Each
row is one fish that stands out for one specific reason (high inbreeding,
member of a large family cluster, rare inversion haplotype carrier, marker
validation control, etc.). The page deliberately keeps interpretation
minimal — it is a **visual highlight list**, not a decision table.

Three traffic-light statuses:
- 🟢 **green** — useful / preserve / good control
- 🟡 **yellow** — review / balance / moderate concern
- 🔴 **red** — caution / high-risk flag

## Page layout

- **Tab label**: "7 breeding" (between "6 inbreeding" and "8 about")
- **DOM ID**: `page8`
- **Default view**: 4-column simple table (Sample · Highlight · Value · Source)
- **Detailed view** (toggle): adds Status text + Notes/origin columns (6 total)
- **Filters**: by status (green/yellow/red), by source (5 categories)
- **Buttons**: load `breeding_highlights.json`/.tsv, export CSV, reset
- **Source counter cards**: live counts per category + total
- **Empty state**: friendly message before any data is loaded

## Two paths to populate rows

### 1. Auto-derived from atlas state (defensive; activates when data layers load)

| Source | Trigger | Status | Reads from |
|--------|---------|--------|-----------|
| Diversity | top-1% F_ROH | red | `state.perSampleStats[i].F_ROH` |
| Diversity | top-1% θπ | green | `state.perSampleStats[i].theta_pi` |
| **Diversity** | **focal cohort hatchery health verdict (page 9)** | **green / yellow / red** (verdict-mapped) | `state._hatcheryHealth.cohorts` |
| Relatedness | members of largest family cluster | yellow | `state.familyClusters.clusters` |
| Burden | top-1% deleterious burden | red | `state.perSampleStats[i].deleterious_burden` |
| Inversions | rare-haplotype carrier (cohort freq ≤ 5%) | green | `state.inversionCarriers` |
| Markers | positive/HET/negative controls | green | `state.markerControls` |

The hatchery-health auto-row is **cohort-level**: the focal cohort name
appears in the Sample column and the colored verdict text (green / amber /
red, e.g. "Caution — historical erosion + recent inbreeding") appears in the
Highlight column via the row schema's optional `highlight_html` field. The
Value column shows `F_ROH 0.277 · H/H_ref 0.29 (H_ref 1.56e-2)`. The
ref_url links to `#page9`.

These hooks are no-ops when their state slot is missing (current scaffold
state). When future modules emit JSON layers (MODULE_3 per-sample F_ROH,
MODULE_2B family clusters, MODULE_CONSERVATION burden, Inversion Atlas
karyotype assignments + marker panel), the corresponding rows start
appearing automatically.

### 2. User upload via `breeding_highlights.json` (or `.tsv`)

```json
{
  "metadata": { "panel_name": "...", "date": "..." },
  "rows": [
    { "sample": "CGA001",
      "highlight": "High inbreeding on Chr17",
      "value": "top 1% F_ROH",
      "source": "Diversity",
      "status": "red",
      "notes": "F_ROH = 0.387",
      "ref_url": "Diversity_atlas.html#page5" }
  ]
}
```

User rows override auto-derived rows on `(sample, source)` match (the
overlay is annotated with `_origin: 'overlaid'`). User-only rows
(no auto match) append with `_origin: 'user'`.

TSV variant: header row with columns `sample · highlight · value · source ·
status · notes · ref_url`, tab-separated.

A complete demo is at `breeding_highlights.demo.json` in this directory
(8 rows: 2 red, 2 yellow, 4 green, across all 5 source categories).

## Tests

```bash
# Breeding page rendering + auto-derivation (24 checks)
NODE_PATH=/path/to/node_modules node scripts/test_population_atlas_breeding.js
# Expected: "[bp-smoke] ALL CHECKS PASSED" (24 checks)

# Cross-atlas hash routing (10 checks across both atlases)
NODE_PATH=/path/to/node_modules node scripts/test_hash_routing.js
# Expected: "[hash-routing] ALL SUITES PASSED"
```

The smoke test exercises:
- DOM presence (page8 div, "7 breeding" tab, "8 about" renumber)
- All breeding-page helpers exposed on `window`
- Empty-state message before any data loads
- Each of the 5 auto-derivation hooks (Diversity F_ROH, Diversity θπ,
  Relatedness, Burden, Inversions, Markers)
- Source counter cards (correct counts per category)
- Status sorting (red → yellow → green)
- Simple (4 cols) vs detailed (6 cols) toggle
- User JSON upload (new + overlay on existing auto row)
- Filters (status + source) isolating subset rows
- TSV upload parsing
- CSV export click handler
- Reset clears user rows but restores auto-derived defaults
- Status badge transitions planned → ready
- Empty state restores when upstream data drops

## Manuscript-ready methods sentence

> Sample highlights are surfaced from the cohort's per-sample population
> statistics (F_ROH, mean θπ, mean heterozygosity, deleterious burden) and
> from the inversion atlas's karyotype assignments + marker-control set.
> Top-1% thresholds are computed from the 226-sample empirical distribution
> per metric. Highlights are a review aid, not a breeding-decision table;
> final pairing decisions should be made by the breeding committee with
> reference to pedigree records, the full Diversity Atlas drill-downs, and
> the Inversion Atlas marker validation status.

## Cross-atlas hash routing (v4 turn 79)

The breeding page links to multiple sister atlases:
- Marker rows → `Inversion_atlas.html#page18` (marker readiness panel)
- F_ROH / θπ rows → `Diversity_atlas.html#page5` / `#page2`
- Family-cluster rows → Population Atlas's own `#page3`

For these links to land on the right tab (rather than whatever was last
visited via localStorage), both atlases now have **hash routing**:

1. On initial load, `location.hash` takes precedence over the saved
   `localStorage.activeTab` when the hash matches a real `data-page` value.
2. `hashchange` events fire correctly for in-page navigation.
3. Bogus or missing hashes fall back gracefully (no navigation change).

**Inversion Atlas** marker panel methods note now also includes a
back-link to `Population_atlas.html#page8` so users can navigate from the
marker panel directly to the breeding-review highlights for those control
samples.

The `test_hash_routing.js` smoke test exercises both directions (10
checks across 3 suites — Inversion Atlas hash routing, Population Atlas
hash routing, empty-hash fallback).

---

# Hatchery health page (page 9, visible label "7 hatchery health")

Operationalises the manuscript's **F_ROH × H plane** framework — the novel
analytic contribution introduced in the Discussion as
**diversity-contextualised F_ROH** (F_ROH | H).

## What the framework does

F_ROH alone is ambiguous across studies. A cohort with F_ROH = 0.30 on a
diversity-replete background reflects **recent close-kin inbreeding only**
and is manageable by outcrossing. The same F_ROH = 0.30 layered on top of
already-eroded heterozygosity reflects **sustained demographic decline** —
outcrossing within the stock cannot recover the missing diversity because
the non-ROH genome itself is depleted.

Reporting F_ROH and H jointly resolves this ambiguity. Every stock becomes
a point in 2D space; the position tells you which biological process is
driving the inbreeding signal and which management response is appropriate.

## Quadrants of the plane (verdict map)

| F_ROH | H / H_ref | Verdict (colored text) | Interpretation |
|-------|-----------|------------------------|----------------|
| low   | high      | **Healthy** (green)    | baseline state; preserve |
| high  | high      | **Recent inbreeding only** (amber) | recent close-kin events; manageable by pair selection |
| low   | low       | **Diversity erosion only** (amber) | old founder bottleneck; current pairings ok but limited gene pool |
| high  | low       | **Caution — historical erosion + recent inbreeding** (red) | worst case; outcrossing within stock cannot recover diversity |

Quadrant cuts are computed from the **median F_ROH** and **median H/H_ref**
of currently-loaded cohorts (≥2 cohorts required); when only one cohort is
loaded the page falls back to fixed thresholds (F_ROH = 0.10, H/H_ref = 0.5).

## H_ref resolution chain

1. Explicit per-cohort `H_ref` field wins.
2. Otherwise: **species-peer fallback** — the maximum `H_per_site` across
   loaded cohorts of the same species.
3. Otherwise: **global fallback** — the maximum `H_per_site` across all
   loaded cohorts.

This makes the page useful even when only one cohort has H_ref filled in
(the others auto-reference against the highest-H peer).

## Focal cohort

The focal cohort gets a larger point on the SVG scatter, a bold label, and
a **FOCAL** badge in the cohort table.

- If any cohort has `is_focal: true`, that one is focal.
- Otherwise the focal cohort defaults to the **highest-F_ROH** cohort
  (most likely the one under hatchery review).

The focal cohort's verdict drives the breeding-page auto-row.

## Coloring rule (per user spec)

**Only the verdict text is colored.** Cohort names, F_ROH values, H values,
and H_ref values stay in normal ink. The CSS classes are:

- `.fhh-good`    — `color: var(--good)`   (green, weight 600)
- `.fhh-review`  — `color: var(--accent)` (amber, weight 600)
- `.fhh-caution` — `color: var(--bad)`    (red, weight 600)

These are also reused on the SVG scatter (point fill colors) and on the
breeding-page auto-row's `highlight_html`.

## Schema — `hatchery_health.json`

```json
{
  "metadata": { "panel_name": "F_ROH × H plane", "date": "..." },
  "cohorts": [
    { "cohort": "this study (226 broodstock)",
      "species": "C. gariepinus",
      "is_focal": true,
      "F_ROH": 0.277,
      "H_per_site": 0.00455,
      "H_ref": 0.0156,
      "ref_url": "Diversity_atlas.html#page1",
      "notes": "Thai hatchery, 9x WGS, ngsF-HMM ≥1 Mb" },
    { "cohort": "Nigerian wild C. gariepinus",
      "species": "C. gariepinus",
      "F_ROH": 0.05,
      "H_per_site": 0.0156,
      "H_ref": 0.0156,
      "citation": "Sanda 2026" }
  ]
}
```

TSV variant accepted with header
`cohort · species · F_ROH · H_per_site · H_ref` (optional `citation`,
`notes`, `is_focal`).

See `hatchery_health.demo.json` for a complete example with this study +
five comparator cohorts (Nigerian wild and farmed *C. gariepinus*, rainbow
trout, American alligator, Iberian lynx).

## Cross-feed to the breeding page

When `hatchery_health.json` is loaded, the focal cohort's verdict is
surfaced as a **Diversity-source row** on the breeding page table:

- Sample column: cohort name (plain ink)
- Highlight column: colored verdict text via inline `<span class="fhh-...">`
- Value column: `F_ROH X.XXX · H/H_ref Y.YY (H_ref Z.ZZe-N)`
- Source column: `Diversity` with traffic-light dot mapped from verdict
- Status mapping: good→green, review_*→yellow, caution→red
- ref_url: `#page9` (clicking jumps back to the hatchery health page)

This way the breeding committee sees stock-level health alongside
per-sample flags in one table.

## Tests

`test_population_atlas_hatchery_health.js` — 28 checks across:
- Tab renumbering (7 hatchery health · 8 breeding · 9 about)
- Quadrant classification at 4 corners (good / review_inbreeding /
  review_erosion / caution)
- H_ref resolution chain (explicit / species-peer / global fallback)
- Focal-cohort detection (explicit flag / highest-F_ROH fallback)
- Empty state
- JSON upload + plot SVG render + table render + summary cards
- Verdict text uses correct colored CSS class
- Cohort name stays in normal ink (NOT colored)
- Breeding page auto-row appears + carries colored highlight_html
- Breeding page row status maps verdict → traffic-light dot correctly
- Breeding row removed when FHH state is reset
- CSS rules for `.fhh-good / .fhh-review / .fhh-caution` present
- Hash routing — `#page9` lands on hatchery health
- CSV export, TSV upload, reset
