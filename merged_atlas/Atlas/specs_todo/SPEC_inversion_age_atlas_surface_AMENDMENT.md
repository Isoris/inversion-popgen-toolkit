# SPEC AMENDMENT — `inversion_age_atlas_surface.md` (drop absolute-My)

**Status**: amendment, drafted 2026-05-06.
**Amends**: `Atlas/specs_todo/SPEC_inversion_age_atlas_surface.md` (turn 117).
**Depends on**: `SPEC_busco_4d_age_brackets.md` (sibling, this batch).
**Reason**: After agreement that all dating methods are unreliable at
hatchery Ne ≈ 20 with sweepstakes reproduction (millions of eggs per
spawn), absolute-My displays were retired in favor of: (a) ranking
voice for relative methods 1 and 2, and (b) explicit three-μ bracketing
for the absolute method (Method 3, BUSCO 4D — see sibling spec).

---

## 1. The decision

The original spec carried absolute-My displays (1.84 My with
parenthetical CI) on every age surface. After honest assessment:

- Cohort Ne ≈ 20 (226 fish from few founders)
- Sweepstakes reproduction (millions of eggs / spawn → wildly skewed
  coalescent tree, sweepstakes-noise dominates standard SFS theory)
- Single published catfish μ (Liu 2023 Siluriformes, ~2× anchor uncertainty)
- Generation-time genuinely uncertain (6 mo vs 4 yr, hatchery vs wild,
  historical vs current)

…all dating methods give per-locus absolute numbers with CIs that are
honest only when they're so wide as to be useless ("1.84 My (0.6–5.5)").
Relative comparisons are robust; absolute single-method point estimates
are not. Manuscript voice must reflect this.

## 2. Three-method age framework (final)

The atlas displays **three** age methods alongside each other on the
candidate-focus page (page 3). All three are visible at once so the
reader / reviewer can triangulate. Manuscript text leads with rankings
and uses the bracketed Method 3 for absolute-time framing.

| Method | What it measures | Output type | Surface |
|---|---|---|---|
| **Method 1** — relative between-inversions | dXY between arrangements compared *across* inversions | `chrom_age_rank` (1 = oldest), ratio vs chrom-median | Page 5 catalogue column |
| **Method 2** — relative within-arrangement | within-class π per arrangement, *inside one inversion* | `rank_within_inversion`, ratio (HOM_INV is 2.94× HOM_REF) | Page 3 panel Row B |
| **Method 3** — absolute via BUSCO 4D + three-μ brackets | dXY at 4-fold-degenerate sites within BUSCO genes inside the inversion, displayed at three μ | `age_my_under_mu_low` / `_mid` / `_high`, each with bootstrap CI | Page 3 panel Row C |

**Methods 1 + 2 are dimensionless** (ranks and ratios). **Method 3 is
absolute** but reported as three numbers spanning the plausible μ range,
not a single point estimate.

The cheat30 GDS panel (already shipped turn 119) provides the **fourth**
piece of evidence on the same page: relative-age binning (OLD /
INTERMEDIATE / YOUNG) + origin classification (single / recurrent).
Together: four readings of the same underlying biology.

## 3. Removed surfaces (vs original spec)

The following were specified in the turn-117 draft and are now **NOT**
to be implemented:

- **Page 5 column display modes "My" and "× chrom-median".** Only
  `rank` and `× chrom-median` survive. The "My" mode is removed.
  Column header tooltip drops the absolute-My language.
- **Page 3 Row A "age: 1.84 My (1.32–2.41)" inline display.** Replaced
  with rank + ratio: "rank 1 of 4 on chrom; dXY 2.4× chrom median".
- **Page 3 Row B per-group My bars.** Replaced with per-group π values
  and ranks: "HOM_INV: π=0.00582 (rank 1 oldest); HOM_REF: π=0.00197
  (rank 2)" + "ratio: 2.94× HOM_INV/HOM_REF".

Everything else from the original spec is preserved: state.inversionAge
storage, JSON loader pattern, Berdan load matrix, Guerrero divergence
profile shape, cross-species overlap indicator, page-5 sortable column
infrastructure.

## 4. Added surfaces (vs original spec)

### 4.1 Page 3 Row C — Method 3 (BUSCO 4D + three-μ brackets)

Renders only when `state.inversionAge[candidate].busco_4d_age_brackets`
is present. Empty state otherwise (one-line "Method 3 not yet computed
on this chromosome — run STEP_C01f_e_emit_busco_4d_age.py").

Layout (~120px tall when populated, ~24px empty):

```
TASK 3: absolute age (BUSCO 4D-neutral sites, three-μ bracketing)
  μ_low  = 1×10⁻⁹/site/year      5.5 My (3.9–7.2)   slow molecular clock
  μ_mid  = 3×10⁻⁹/site/year      1.84 My (1.32–2.41)  Liu 2023 Siluriformes
  μ_high = 9×10⁻⁹/site/year      0.61 My (0.44–0.80)  fast molecular clock

  n BUSCO 4D sites used: 2,847
  ⓘ Catfish-specific μ unknown. Three brackets span plausible teleost
  range. The truth is in there. Ranking across inversions robust
  under any μ.
```

### 4.2 Page 3 Row D — supporting population-genomic stats

Renders when `state.regionPopstats` (see §6 below) has data for this
candidate. Displays the per-window traces produced by
`region_popstats.c` over the inversion interval:

- **dXY trace** between HOM_REF and HOM_INV (line plot, x = bp window
  center, y = per-site dXY)
- **θπ trace per arrangement** (multi-line: HOM_REF, HOM_INV, HET if
  applicable)
- **Tajima's D shape** (Guerrero 2012 age-proxy classes: U_SHAPE,
  CENTRAL_PEAK, BREAKPOINT_PEAKS, FLAT, MIXED) with the per-window D
  trace below the shape badge

These stats are **diagnostic / triangulation aids**, not age estimates.
They support the Methods 1–3 readings (e.g. U-shape Tajima's D plus
high HOM_INV π plus Method-3 μ_mid age in the My range = "old neutral
inversion under recombination suppression").

## 5. JSON-shape addition: `busco_4d_age_brackets`

Extension to the existing `inversion_age_v1` schema. Every inversion
gets an optional `busco_4d_age_brackets` block (omitted if not yet
computed):

```json
{
  "inversions": [
    {
      "candidate_id": "LG28_cand_01",

      // Existing blocks kept verbatim:
      "between_arrangement_age": { ... },     // Method 1 ranks live here
      "within_arrangement_age": { ... },      // Method 2 ranks live here
      "load": { ... },
      "divergence_profile": { ... },

      // New block (see SPEC_busco_4d_age_brackets.md for full schema):
      "busco_4d_age_brackets": {
        "n_busco_genes_in_inversion": 14,
        "n_4d_sites_used": 2847,
        "dxy_4d_between_arrangements": 0.0078,
        "mu_low":  { "value": 1.0e-9, "age_my": 3.90, "age_my_ci95": [2.78, 5.12], "rationale": "slow teleost clock" },
        "mu_mid":  { "value": 3.0e-9, "age_my": 1.30, "age_my_ci95": [0.93, 1.71], "rationale": "Liu 2023 Siluriformes" },
        "mu_high": { "value": 9.0e-9, "age_my": 0.43, "age_my_ci95": [0.31, 0.57], "rationale": "fast teleost clock" }
      }
    }
  ]
}
```

The three μ values are **fixed defaults**. They are not user-tuneable
in the atlas (changing μ would invalidate every age estimate; if a new
catfish μ becomes available, regenerate the JSON LANTA-side).

## 6. JSON-shape addition: `region_popstats_v1`

NEW per-candidate JSON layer for the supporting stats (Row D). Output
of `region_popstats.c` reformatted by a small joiner script
(`STEP_C01f_f_emit_region_stats_per_candidate.py`, ~80 lines, NEW —
spec'd as a sibling task in `SPEC_busco_4d_age_brackets.md` §10).

```json
{
  "tool": "region_popstats_v1",
  "schema_version": 1,
  "generated_at": "2026-MM-DDTHH:MM:SSZ",
  "per_candidate": {
    "LG28_cand_01": {
      "windows": [
        {
          "start_bp": 15100000, "end_bp": 15600000,
          "dxy_ref_inv": 0.011,
          "fst_ref_inv": 0.31,
          "pi_ref": 0.0019, "pi_inv": 0.0058, "pi_het": 0.0042,
          "tajima_d_ref": -0.42, "tajima_d_inv": -0.15, "tajima_d_het": -0.30,
          "n_callable_kb": 487
        },
        ...
      ],
      "tajima_shape_class": "U_SHAPE",
      "tajima_shape_confidence": 0.82
    }
  }
}
```

State slot: `state.regionPopstats`, same per-chrom keying as
`state.cheat30Results`. LocalStorage key: `inversion_atlas.region_popstats`.

## 7. Voice discipline (final)

| Surface | Voice |
|---|---|
| Page 5 "rel age" column header | "rel age" — never "age" |
| Page 3 Row A label | "TASK 1: between-inversion ranking" |
| Page 3 Row B label | "TASK 2: within-arrangement ranking" |
| Page 3 Row C label | "TASK 3: absolute age (BUSCO 4D, three-μ)" |
| Page 3 Row D label | "supporting population-genomic stats" |
| Method 1 numeric display | rank + ratio, NEVER years |
| Method 2 numeric display | rank + π ratio, NEVER years |
| Method 3 numeric display | three rows, three μ, three Mys, ALWAYS as a triple |
| Method 3 caveat (always reachable) | "Catfish-specific μ unknown. Three brackets span plausible teleost range." |
| Manuscript bundle | Methods 1+2 ranks first, Method 3 brackets second, support stats third |
| Tajima's D | only as Guerrero-2012 age-proxy class, **never** as selection test |

## 8. Implementation order (next chat)

1. Loader + state slot for `region_popstats_v1` JSON (~50 LOC) —
   simplest, mirrors existing `cross_species` loader.
2. Page-3 Row D supporting-stats panel (dXY trace + θπ multi-line
   trace + Tajima D shape badge) (~150 LOC).
3. Loader + state slot for `inversion_age_v1` JSON (~80 LOC).
4. Page-3 Row A (Method 1 ranking strip) (~40 LOC).
5. Page-3 Row B (Method 2 per-arrangement ranking + ratio bars) (~60 LOC).
6. Page-3 Row C (Method 3 three-μ-bracket display) (~80 LOC) —
   gated by `busco_4d_age_brackets` presence; empty state otherwise.
7. Page-5 catalogue "rel age" column with rank / × chrom-median
   display modes (~80 LOC).
8. Manuscript bundle entry combining all three methods + support
   stats (~80 LOC).
9. Help-page row (~10 LOC).
10. Tests (~250 LOC, ~35 assertions, all four loaders + all four panel
    rows + page-5 column).

**Total atlas-side**: ~880 LOC. Larger than the original 670 estimate
because Row D (region_popstats) is new scope.

## 9. What this is NOT

- Not changing the LANTA-side `STEP_C01f_c_burden_regression.R`
  pipeline. The `between_arrangement_age` and `within_arrangement_age`
  blocks of `inversion_age_v1.json` are unchanged. Only the **atlas
  display** stops showing absolute My from those blocks.
- Not requiring `busco_4d_age_brackets` to ship the rest. Methods 1 + 2
  + supporting stats can land before BUSCO 4D is computed. The Row C
  panel just shows its empty state.
- Not asking the atlas to pick a μ. The three μ values are pipeline-side
  decisions baked into the JSON. The atlas displays whatever it gets.
- Not displaying dXY age in years from the existing
  `between_arrangement_age.age_my_year` field. That field can stay in
  the JSON (for backward compatibility, manuscript bundle math, etc.)
  but is not surfaced visually on page 3 or page 5.

## 10. Tests (atlas-side, list)

- `state.regionPopstats` round-trips localStorage
- `state.inversionAge` round-trips localStorage
- Page-3 Row A shows rank, never My
- Page-3 Row B shows π ratios, never My
- Page-3 Row C shows three rows when busco_4d block present
- Page-3 Row C shows empty state when busco_4d block absent
- Page-3 Row D dXY trace renders for K=3 candidates
- Page-3 Row D dXY trace renders for K=6 candidates
- Page-3 Row D Tajima D shape class badge matches JSON value
- Page-5 column sorts by rank in ascending mode, descending mode
- Page-5 column tooltip never says "My" without three-μ bracket framing
- Manuscript bundle TSV includes all three methods + Tajima shape
- Manuscript bundle markdown leads with rank, never with single-My

End of amendment.
