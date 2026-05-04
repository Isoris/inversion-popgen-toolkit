# SPEC — Two-task atlas surface for inversion age (between-inversion AND within-arrangement)

**Status**: forward-looking spec. Not implemented. Drafted turn 117 in
response to a side-by-side comparison with Wang et al. 2023 (Nat. Commun.
14:1567) Fig. 3 / Table 2 (rice pan-genome inversion index). User direction:

> *"This inversion vs that one — which is older? AND within the same
> inversion, which arrangement / group is older relative to others. Two
> very different tasks."*

> *"Generation time we don't have... so just leave it then. Use mu_year
> directly."*

**This spec is atlas-side only.** The producing pipeline,
`STEP_C01f_c_burden_regression.R` v9.7 (1584 lines), already exists and
already emits the per-candidate relative-age + arrangement-specific-load +
Guerrero-2012 divergence-profile-shape outputs. **One small new HPC script
(`STEP_C01f_d_emit_age_json.py`, ~80 lines)** converts those outputs to a
JSON the atlas can load. Atlas-side wiring is the bulk of the work. The
existing C01f_c §4b will need a small modification (one switch from
`age_per_gen` to `age_per_year` framing — see §1.2 below).

**Reading order**: this spec → `STEP_C01f_c_burden_regression.R` (the
producer; sections 4b/4c/4d/4e are the relevant ones) → atlas
`_renderCandidateFocus` (page 3; this is the primary surface) → atlas
`_renderCatalogue` (page 5; secondary surface, sortable column only) →
`SPEC_focal_vs_background_widget.md` (architectural sibling — same
loader-dispatch + `state.<layer>` data-flow pattern).

---

## 0. Two distinct tasks the user wants the atlas to answer

The user's clarification was that **these are two separate questions** and
the atlas needs to answer both, not collapse them:

### Task A — Between-inversion ranking (the catalogue question)

> *"Inversion 1 vs inversion 2 vs inversion 3 — which is older?"*

Surface: catalogue (page 5) sortable column. Data: per-candidate
`age_my_year` derived from between-arrangement dXY in C01f_c §4b. Rank
within chrom AND across chroms both meaningful. Display as `rank` (1 =
oldest), `My` (with CI), or `× chrom median` (relative to chrom baseline).

### Task B — Within-arrangement-group ranking (the candidate-focus question)

> *"Within the same inversion, which arrangement / group is older relative
> to others — HOM_REF or HOM_INV (or H1/H2/H3 for K=6 cases)?"*

Surface: candidate focus panel (page 3). Data: per-group **within-class
biSNP density** (i.e. within-class π / 2μ, computed inside the inversion).
The arrangement that has accumulated more within-class variation is the
older arrangement (longer time since fixation → more accumulated
substitutions on the inverted background that recombination can't purge).

This is **NOT** what C01f_c §4b currently emits — §4b computes between-class
dXY (a property of the *split*, not of either arrangement individually).
The within-class π values are produced by `region_popstats.c` already, but
they're not currently piped into the same per-candidate JSON. The new
emitter script `STEP_C01f_d` joins them.

### What this avoids

- **Generation-time dependency.** Both tasks compute `age_my = dXY / (2 × μ_year)`
  in years directly, never multiplying by generation time. Hatchery *C.
  gariepinus* generation time is 6 months to 1 year (vs. wild assumption of
  2 yr), but μ_year is the *robust* number from the Siluriformes phylogeny
  and doesn't depend on g at all. **Bypass the problem entirely.**
- **Subclade UPGMA inside HOM_INV.** The K=6 atlas mode already exposes
  within-arrangement substructure via local PCA. Don't add a redundant
  second method on the same question.

---

## 1. Why this is reusable from existing infrastructure

### 1.1 Pipeline producer — `STEP_C01f_c_burden_regression.R` v9.7

Lives at `inversion_codebase_v8.5/STEP_C01f_c_burden_regression.R`.
Sections relevant to the atlas:

- **§4b — Relative age via dXY.** Computes per-candidate dXY between
  HOM_REF and HOM_INV pools. Currently emits `age_mya` using μ_per_gen × gen
  time. **Modification needed**: also emit `age_my_year` using μ_year
  directly (`dXY / (2 × μ_year)`), with μ_year = 3 × 10⁻⁹ per site per year
  (Siluriformes-anchored, Liu et al. 2023). This is one new column in the
  catalogue patch TSV, not a re-run.
- **§4c — Per-chromosome chronology.** Ranks inversions on each chromosome
  oldest→youngest. Output `chrom_age_rank` (1 = oldest on chrom). Used
  by Task A.
- **§4d — Arrangement-specific load.** Per-candidate
  `n_inv_private_load`, `n_ref_private_load`, `n_shared_load`,
  `inv_private_fraction`. Berdan-2023 prediction display.
- **§4e — Guerrero 2012 divergence profile shape.** Per-candidate
  `divergence_profile_shape` (FLAT / U_SHAPE / CENTRAL_PEAK /
  BREAKPOINT_PEAKS / MIXED) + raw per-window dXY trace. Surfaces on page 3
  in Task B's panel.

### 1.2 Mutation rate — μ_year = 3 × 10⁻⁹ per site per year

From Liu et al. 2023 (Mar Life Sci Technol, electric catfish *Malapterurus
electricus* genome paper, Methods § Estimates of effective population
size): "the mutation rate per site per year was set at 3 × 10⁻⁹ estimated
by r8s (Sanderson 2003)". This is the value derived from the Siluriformes
phylogeny calibrated against teleost fossils — the **only catfish-anchored
neutral rate** in the literature.

The paper goes on to derive μ_per_gen = 6 × 10⁻⁹ assuming g = 2 yr for
*M. electricus*. **We do NOT use this.** Hatchery *C. gariepinus* generation
time differs (6 mo to 1 yr). Use μ_year directly:

```
age_my_year = dXY / (2 × μ_year × 1e6)
            = dXY / (2 × 3e-9 × 1e6)
            = dXY / 6e-3
            = dXY × 166.67   (My)
```

**Citation footnote on every age display**: "μ = 3 × 10⁻⁹ per site per
year, Liu et al. 2023 *Mar Life Sci Technol* (Siluriformes-anchored,
r8s-estimated). Generation time not assumed."

### 1.3 Within-class biSNP density — `region_popstats.c`

Already computes within-group θπ per window. The new emitter script
queries the existing C01f_c output for per-candidate per-group θπ averaged
over windows inside the inversion. Per-group **within-class age** is then:

```
age_my_year_within_<group> = pi_<group>_inside_inversion / (2 × μ_year × 1e6)
```

This is a per-group scalar. The older arrangement has higher within-class
π (more time accumulated on the same background under recombination
suppression). Displayed on page 3.

### 1.4 Soft-anchor philosophy (kept verbatim)

Three layers of uncertainty the spec must surface, NOT hide:

1. **μ_year is genus-anchored, not lineage-anchored.** ~2× uncertainty on
   absolute My; **ranking is robust** (the user's stated priority).
2. **Ne ≈ 20 hatchery bottleneck noise.** Demographic stochasticity at
   ≤10 My scales is folded into the bootstrap CI from §4b.
3. **Selection vs age confound.** Balancing selection elevates within-class
   π. Manuscript framing: "consistent with old inversions OR balancing
   selection OR both, ranking robust under either".

The atlas displays these caveats inline — never bare My without the
parenthetical CI, every column header / panel caption flags "soft anchor".

---

## 2. Atlas-side data layer

### 2.1 New JSON shape: `inversion_age_v1`

`STEP_C01f_d_emit_age_json.py` (~80 lines, NEW) joins C01f_c's
catalogue-patch TSV with the within-group π values and emits:

```jsonc
{
  "tool": "inversion_age_v1",
  "schema_version": 1,
  "generated_at": "2026-MM-DDTHH:MM:SSZ",
  "cohort": {
    "n_samples": 226,
    "species": "Clarias gariepinus",
    "reference_haplotype": "fClaHyb_Gar_LG.fa"
  },
  "calibration": {
    "mu_per_site_per_year":     3.0e-9,
    "mu_source_citation":       "Liu et al. 2023, Mar Life Sci Technol, doi:10.1007/s42995-023-00197-8",
    "mu_method":                "r8s-estimated, Siluriformes phylogeny, fossil-calibrated",
    "generation_time_assumed":  null,
    "generation_time_note":     "Per-year rate used directly. Generation time NOT assumed.",
    "uncertainty_note":         "mu is a Siluriformes-genus point estimate. Absolute My values carry approx 2x uncertainty; relative ranking robust.",
    "source_script":            "STEP_C01f_c_burden_regression.R v9.7 §4b, with §4b' year-based reframe"
  },
  "genome_baseline": {
    "by_chrom": {
      "C_gar_LG12": { "median_age_my_year": 4.21, "n_callable_kb": 24831 },
      "C_gar_LG28": { "median_age_my_year": 3.89, "n_callable_kb": 22104 },
      ...
    }
  },
  "inversions": [
    {
      "candidate_id":              "LG28_cand_01",
      "chrom":                     "C_gar_LG28",
      "start_bp":                  15115000,
      "end_bp":                    18005000,
      "size_bp":                   2890000,

      // ── TASK A: between-inversion age ────────────────────────────────
      "between_arrangement_age": {
        "dxy_hom_ref_hom_inv":   0.0110,    // average dXY per site between groups
        "age_my_year":           1.84,
        "age_my_year_ci95":      [1.32, 2.41],
        "chrom_age_rank":        1,         // 1 = oldest on this chrom
        "n_on_chrom":            4,
        "fst_hom_ref_hom_inv":   0.308,     // concordance proxy
        "gds_gap":               0.42       // concordance proxy
      },

      // ── TASK B: within-arrangement age (per group) ───────────────────
      "within_arrangement_age": {
        "K": 3,
        "groups": [
          {
            "group_id":           "HOM_REF",
            "n_samples":          60,
            "pi_inside_inversion": 0.00197,
            "age_my_year":        0.33,   // pi / 2μ_year — older arrangement = higher
            "age_my_year_ci95":   [0.27, 0.39],
            "rank_within_inversion": 2    // 1 = oldest of the 3 groups
          },
          {
            "group_id":           "HOM_INV",
            "n_samples":          60,
            "pi_inside_inversion": 0.00582,
            "age_my_year":        0.97,
            "age_my_year_ci95":   [0.84, 1.11],
            "rank_within_inversion": 1    // oldest — what we want to flag
          },
          {
            "group_id":           "HET",
            "n_samples":          106,
            "pi_inside_inversion": 0.00428,
            "age_my_year":        0.71,
            "age_my_year_ci95":   [0.62, 0.81],
            "rank_within_inversion": null,  // HET is a *combination*, not a homozygous arrangement; rank N/A
            "note":               "HET = combination of HOM_REF and HOM_INV haplotypes; pi here reflects the SUM of both arrangements' diversity within heterozygotes, not a separate arrangement. Reported for completeness, ranking N/A."
          }
        ],
        "older_arrangement": "HOM_INV",
        "older_minus_younger_ratio": 2.94,    // 0.97 / 0.33 — interpretive shortcut
        "interpretation": "HOM_INV carries 2.9x more within-class variation than HOM_REF, consistent with HOM_INV being the older / longer-recombination-suppressed arrangement"
      },

      // ── Existing C01f_c outputs (already produced) ────────────────────
      "load": {
        "n_inv_private":          18,
        "n_ref_private":          4,
        "n_shared":               12,
        "inv_private_fraction":   0.529,
        "berdan_prediction_met":  true
      },
      "divergence_profile": {
        "shape":                  "U_SHAPE",
        "shape_confidence":       0.82,
        "trace_window_centers_mb": [15.20, 15.31, ..., 17.93],
        "trace_dxy_per_kb":        [3.14, 4.02, ...,  3.78],
        "interpretation":         "neutral old (recombination-suppressed long enough for U-shape to develop)"
      }
    },
    ...
  ]
}
```

For K=6 candidates, `within_arrangement_age.groups` extends to 6 entries
(H1H1, H1H2, H1H3, H2H2, H2H3, H3H3), with `rank_within_inversion` only
populated for the homozygous diagonal (H1H1, H2H2, H3H3 — these are the
"arrangements" whose internal π is interpretable as age proxy). The
heterozygous off-diagonals report π for completeness but with
`rank_within_inversion: null`.

### 2.2 Loader integration

Mirrors the cross_species + repeat_density + ncRNA loader pattern exactly:

- `_classifyJSONKind` recogniser block: `data.tool === 'inversion_age_v1'`
  AND `Array.isArray(data.inversions)`
- `_storeInversionAge` / `_persistInversionAge` / `_restoreInversionAge` /
  `_clearInversionAge` (~80 lines combined; copy-paste from
  `_storeNcRNADensity` family with field renames)
- Boot-time restore alongside cross_species and dotplot restores
- Stored on `state.inversionAge`, with two indexed lookups built at load:

```javascript
state.inversionAge = {
  schema_version: 1,
  ...,
  by_candidate_id: Map<id, inversion>,    // primary lookup, used by both surfaces
  by_chrom: Map<chrom, inversion[]>,      // used for chrom-rank renumbering on filtered catalogue views
};
```

LocalStorage key: `inversion_atlas.inversion_age` (matches naming convention).

---

## 3. Atlas-side surfaces — page 3 is primary

User direction:

> *"For your previous spec it's mostly like we use the candidate page.
> Because candidate page has the groups for the inversion. It cannot be
> perfect but we almost done with all data to put on paper, it's no need
> to perfect anyway. If not all pages match it's ok."*

So **page 3 is the primary surface** (it has the karyotype-group machinery
already, where Task B's per-group display fits naturally). Page 5 catalogue
is a smaller secondary surface (Task A's column). Page 16 cross-species
overlap indicator is optional polish — defer if time-constrained.

### 3.1 Page 3 — primary panel: "inversion age & arrangement structure"

Sibling to the existing class-summary panel. Renders when
`state.inversionAge` is loaded AND a candidate is active. Four rows:

#### Row A — Task A summary (between-inversion context)

A single horizontal strip showing where this candidate sits in the
chrom-wide age ranking:

```
TASK A: between-inversion ranking on C_gar_LG28
   chrom rank: 1 / 4   |   age: 1.84 My (1.32–2.41)
   timeline: [●─────●─●─●]    ← this candidate is the leftmost (oldest)

   μ_year = 3×10⁻⁹/site/year (Siluriformes-anchored, Liu et al. 2023).
   Generation time NOT assumed. Ranking robust; absolute My ~2× uncertain.
```

#### Row B — Task B core: per-group within-class age

The arrangement-by-arrangement ranking the user explicitly asked for. A
horizontal bar chart, one bar per group, with the older arrangement
flagged at the top:

```
TASK B: within-arrangement age (which group has accumulated more variation
        on its own background — the older arrangement)

  HOM_INV  ████████████████████████  0.97 My (0.84–1.11)   ★ OLDER
  HOM_REF  ████████                  0.33 My (0.27–0.39)
  HET      ██████████████████        0.71 My (0.62–0.81)   (combination — rank N/A)

  HOM_INV is 2.9× older than HOM_REF — consistent with HOM_INV being
  the longer-recombination-suppressed arrangement.
```

For K=6 candidates the chart shows 6 bars, with the three homozygous
diagonal classes (H1H1, H2H2, H3H3) flagged for ranking and the three
heterozygous classes (H1H2, H1H3, H2H3) annotated `rank N/A`. The bar
colors reuse the existing pca_scrubber_v3 palette for consistency with
the rest of the candidate-focus page.

#### Row C — divergence profile shape (Guerrero 2012)

The dXY-along-the-inversion trace as a small SVG (200×80 px). Background
tinted by shape class (FLAT/U_SHAPE/CENTRAL_PEAK/BREAKPOINT_PEAKS/MIXED).
Caption shows the shape label + interpretation:

```
[U_SHAPE plot]    U-shape — neutral old (recombination suppressed long
                  enough for divergence to peak at center, erode toward
                  breakpoints via residual gene flux). Shape conf: 0.82.
```

#### Row D — arrangement-specific load matrix (Berdan 2023)

Small 2×3 table:

```
                   #variants    fraction
INV-private        18           0.53
REF-private         4           0.12
shared             12           0.35
                              ──────────
Berdan prediction (INV > REF private):  ✓ MET
```

Empty-state hint when the candidate isn't in the JSON, mirroring the
catalogue-column behaviour.

### 3.2 Page 5 — secondary: "rel age" column

`_renderCatalogue` gains a sortable column `rel age` (Task A only — Task B
doesn't fit a single-column display because it's per-group). Three display
modes via header chip:

- **`rank`** (default, Wang-2023-style): integer rank within the visible
  catalogue, 1 = oldest. Renumbered when filters change so chrom-only views
  show 1..n_on_chrom.
- **`My`**: numeric My with CI, e.g. `1.84 (1.32–2.41)`. Background tinted
  cool → warm for young → old.
- **`× chrom`**: relative to chrom median, e.g. `2.4× chrom median`.
  Cross-chrom comparison without committing to absolute My.

Sort handler uses the underlying `between_arrangement_age.age_my_year`
across all three display modes. Empty cell `—` with tooltip "no age
estimate (n_HOM_INV<5 or insufficient SNPs)" when:

- `state.inversionAge` not loaded
- Candidate not in the JSON (catalogue grew since last pipeline run)
- HOM_INV carrier count < 5

Column-header tooltip:

> Relative between-inversion age via dXY/2μ_year (μ = 3×10⁻⁹/site/year,
> Siluriformes-anchored). Absolute My uncertain ~2×; ranking robust.
> Click to sort; chip toggles display mode. For per-group within-arrangement
> ranking, see page 3.

The cross-link to page 3 is important — it makes explicit that the
catalogue answers Task A only and the per-group Task B answer lives on
page 3.

### 3.3 Page 16 — cross-species overlap indicator (OPTIONAL)

When a cs_breakpoint coincides with a polymorphic candidate (existing
`bp.candidate_overlap` array, ✪ flag), the page-16 focus header gains a
single-line indicator:

```
✪ coincides with polymorphic candidate LG28_cand_01 — between-arrangement
age 1.84 My (1.32–2.41), chrom rank 1/4. Cross-species event predates
Gar–Mac split (>30 My); coincidence here = recurrent fragility hotspot,
not shared inheritance of the same event.
```

The "predates Gar–Mac split" framing uses the divergence time from Liu
et al. 2023 Fig. 3A (Siluriformes phylogeny shows *I. punctatus* and
*M. electricus* shared ancestor ~30 My ago — *Clarias* would be similarly
deep). This is the testable claim: any polymorphic catfish inversion is
younger than the Gar–Mac species split, so cs-bp + polymorphic-candidate
coincidence is **always** recurrent fragility, never shared inheritance.

Defer this surface if pressed for time — the manuscript can make the
argument from page-3 + page-5 alone.

### 3.4 Help-page row

```html
<tr><td><b>18 inversion age & arrangement structure</b></td><td>Two-task
soft-anchor relative-age layer for polymorphic inversions in the 226-fish
cohort. Loads from <code>inversion_age_v1.json</code>, output of
<code>STEP_C01f_c_burden_regression.R</code> §4b–4e plus the
<code>STEP_C01f_d_emit_age_json.py</code> joiner.

<b>Task A — between-inversion ranking</b>: which inversion is oldest? Surface:
sortable "rel age" column on page 5 with rank / My / × chrom-median display
modes. <b>Task B — within-arrangement ranking</b>: within one inversion,
which arrangement (HOM_REF vs HOM_INV vs heterozygotes; or H1H1/H2H2/H3H3
for K=6 cases) is older? Surface: page-3 candidate-focus panel with
per-group within-class biSNP-density bar chart, divergence-profile shape
(Guerrero 2012), and Berdan-2023 INV-private-vs-REF-private load summary.

<b>Calibration</b>: μ = 3 × 10⁻⁹ per site per year (Liu et al. 2023, Mar Life
Sci Technol, Siluriformes-anchored r8s estimate). Generation time NOT
assumed — per-year rate used directly to bypass hatchery-vs-wild generation-time
ambiguity. Manuscript framing: y-axis labels say "relative biSNP age" or
"My (soft anchor, ranking robust)", never bare My. The Tajima's D-as-
selection-test interpretation is explicitly disabled at Ne≈20; D is only
shown as a Guerrero-2012 relative age proxy.</td></tr>
```

---

## 4. What the atlas does NOT compute (and why)

| Feature | Why not |
|---|---|
| **Subclade UPGMA inside HOM_INV** | The K=6 atlas mode already exposes within-arrangement substructure via local PCA. Adding a second method on the same question is redundant. |
| **Bootstrap CI on age** | C01f_c §4b already emits CI from window-bootstrap. Atlas just displays it. |
| **Generation-time conversion to absolute years** | Bypassed entirely by using μ_year. The hatchery generation time is uncertain (6 mo to 1 yr vs wild 2 yr). Reporting per-year-anchored My sidesteps this. |
| **Per-arrangement Tajima's D as a selection test** | Explicitly ruled out at Ne≈20 (registry caveat already embedded). D shown only as Guerrero-2012 relative-age proxy. |
| **Cross-species (Mac cohort) age comparison** | Day-2 work — needs the Mac cohort processed through C01f_c. Forward-compat extension: optional `species_cohort_id` field, no breaking change. |

---

## 5. Tests

### Pipeline-side (`STEP_C01f_d_emit_age_json.py`)

- Round-trips a known TSV: emit JSON, parse it back, all fields preserved
- Sentinel handling: missing `age_my_year` (n<5 case) emits `null` not 0
- Per-group rank assignment correct: highest-π homozygous group = rank 1
- HET (or other heterozygous K=6 classes) correctly tagged
  `rank_within_inversion: null`
- Schema validation on all required fields

### Atlas-side

- `state.inversionAge` round-trips localStorage
- Page-5 catalogue column sorts numerically across all three display modes
- Page-5 empty-state shows `—` when JSON not loaded
- Page-3 panel re-renders on candidate-switch
- Page-3 panel hides when JSON not loaded (whole panel collapsed, not
  just rows)
- Page-3 Task B chart renders 3 bars for K=3 candidates and 6 for K=6
- Page-3 Task B chart correctly flags "★ OLDER" on the highest-rank-1
  homozygous group
- Page-3 Task B caption shows `older_minus_younger_ratio`
- Page-3 divergence-profile SVG renders with correct shape-class background tint
- Page-16 cs-overlap indicator (if implemented) shows when both
  `state.crossSpecies` AND `state.inversionAge` are loaded
- Soft-anchor caveat tooltip is reachable from every age display

### Voice-discipline check (manual review, not automated)

- Catalogue column header: "rel age" not "age"
- Page-3 Row A label: "TASK A: between-inversion ranking"
- Page-3 Row B label: "TASK B: within-arrangement age"
- Every numeric My display has parenthetical CI
- Calibration footnote on every age display reachable in tooltip

---

## 6. Open questions / explicitly deferred

- **Manuscript-locked Mode A vs live Mode B**: this layer is Mode-A-only.
  Live recomputation in browser is infeasible.
- **Multi-cohort comparison page**: when Mac cohort lands, a per-locus
  Gar↔Mac age-comparison page becomes interesting. Out of scope for day 1.
- **Polishing the cross-species overlap argument** with cs_breakpoint age
  estimates from BUSCO 4D divergence at the breakpoint flanks. Not blocking.
- **Selection-vs-age decomposition**: identifiability problem, not
  attempted here. Manuscript framing: "consistent with X OR Y OR both,
  ranking robust under either".
- **The C01f_c §4b modification**: needs a one-line addition to also emit
  `age_my_year` alongside the existing `age_mya`. Trivially small. Can
  also be done atlas-side at JSON-load time (multiply `age_mya` by
  `gen_time_assumed` to recover dXY, then divide by 2 × μ_year × 1e6) but
  cleaner to fix at source.

---

## 7. Summary

Two distinct user-questions, two atlas surfaces:

- **Task A (between-inversion ranking)** → page 5 catalogue column
  `rel age` with rank / My / × chrom-median display modes
- **Task B (within-arrangement ranking)** → page 3 candidate-focus panel
  with per-group within-class age bar chart + Guerrero divergence-profile
  shape + Berdan load matrix

One small new HPC script (`STEP_C01f_d_emit_age_json.py`, ~80 lines)
joins existing C01f_c outputs with within-group π. One atlas-side data
layer (~80 lines store/persist/restore), two display surfaces (~250 lines
combined for page-3 panel + page-5 column), one help row. Total
atlas-side: ~400 lines. Reuses cross_species / repeat_density / ncRNA
loader pattern verbatim.

**Calibration**: μ_year = 3 × 10⁻⁹ per site per year (Liu et al. 2023,
Siluriformes-anchored). Generation time NOT assumed — per-year rate used
directly. Sidesteps the hatchery-vs-wild g ambiguity (6 mo to 1 yr vs
2 yr) entirely.

**Voice discipline**: every age display flagged "soft anchor, ranking
robust"; CI parenthetical on every numeric My; Tajima's D explicitly
shown only as Guerrero-2012 age proxy not selection test.

End of spec.
