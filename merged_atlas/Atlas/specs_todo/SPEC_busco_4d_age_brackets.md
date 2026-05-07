# SPEC — `busco_4d_age_brackets`: BUSCO 4D-neutral sites + three-μ bracketing for absolute inversion age

**Status**: drafted 2026-05-06.
**Sibling spec**: `SPEC_inversion_age_atlas_surface_AMENDMENT.md` (atlas display).
**Producer**: NEW — `STEP_C01f_e_emit_busco_4d_age.py` on LANTA, runs after
`STEP_C01f_c_burden_regression.R`. Does not exist yet — this spec defines it.
**Reads**: BUSCO ortholog table for *C. gariepinus*, GFF3 of CDS coordinates,
BEAGLE dosage matrix, per-candidate karyotype calls.
**Writes**: extension block in `inversion_age_v1.json` per inversion (see §3).

---

## 1. Why this method

Methods 1 + 2 (between-inversion and within-arrangement ranking) give
**relative** age — robust under any μ, but no absolute years. The
manuscript / reviewer / breeding-program reader will ask "but how old?"

Method 3 answers in absolute years, but only inside an explicit
**bracket** of μ uncertainty. The bracket is the honesty: "the truth is
in [age_low, age_high], we don't know exactly where, ranking is
preserved under any choice." This is more defensible than a single μ
point estimate with a tight statistical CI that hides the calibration
uncertainty.

The BUSCO 4D restriction matters because:

- 4-fold-degenerate sites in BUSCO genes are the closest the genome
  comes to neutral evolution
- They're highly conserved orthologs across teleosts → BUSCO's
  cross-species set is reliable
- Restricting to 4D within BUSCO removes purifying-selection bias
  (1st/2nd codon positions, UTRs, intergenic enhancers) and balancing-
  selection bias (the inversion itself may carry locally selected
  variants)
- The resulting dXY estimator is the cleanest "molecular clock" the
  data allows

## 2. Pipeline — `STEP_C01f_e_emit_busco_4d_age.py`

### 2.1 Inputs

| Input | Source | File | Format |
|---|---|---|---|
| BUSCO ortholog set | run_BUSCO previously | `BUSCO_runs/Cgar_v2/full_table.tsv` | BUSCO TSV |
| CDS coordinates | gffread / GFF3 | `annotation/Cgar_v2.cds.gff3` | GFF3 |
| Dosage matrix | BEAGLE 5.4 | `phase_5_followup/dosage/Cgar_LG{N}.dosage.tsv.gz` | matrix |
| Karyotype calls | C01f_c §3 | `inversion_codebase_v8.5/per_candidate/karyotype_calls.tsv` | TSV |
| Inversion intervals | C01f_c catalogue | `inversion_catalogue_v9.7.tsv` | TSV |
| Sample list | NAToRA-pruned | `n81_unrelated.list` | one sample per line |

### 2.2 Algorithm

For each inversion in the catalogue:

1. **Find BUSCO genes overlapping the inversion interval.**
   Match BUSCO `start..end` (chromosome, start_bp, end_bp) against
   inversion `[start_bp, end_bp]`. Require ≥ 80% of the BUSCO gene
   inside the inversion (avoid partial-overlap noise from genes
   straddling breakpoints).

2. **Extract 4-fold-degenerate sites within those BUSCO genes.**
   For each CDS feature in the GFF3 belonging to a matched BUSCO:
   - Read the codon frame (CDS phase + start_bp)
   - Identify each codon's third position
   - Keep only third positions where all four nucleotide substitutions
     give the same amino acid (4D sites: codon-table lookup)
   - Filter: site must be biallelic in the cohort (rare alleles or
     multiallelics dropped), MAF ≥ 0.01, callable in ≥ 80% of samples

3. **Compute dXY at those 4D sites between HOM_REF and HOM_INV.**
   Use Hudson 1992 unbiased dXY formula (already implemented in
   `region_popstats.c`). Pool: 4D sites only, restricted to BUSCO genes,
   restricted to the inversion interval. Output: `dxy_4d_between_arrangements`
   (per-site mean across the 4D site set).

4. **Bootstrap CI by resampling 4D sites.**
   1000 bootstrap replicates, sample with replacement at the *site*
   level (not the gene level — site-level bootstrap is the standard
   for dXY CI). Report 2.5% / 97.5% quantiles.

5. **Apply three μ values to compute the bracket.**
   ```
   age_my = dxy_4d / (2 × μ × 1e6)
   ```
   for μ ∈ {μ_low, μ_mid, μ_high}. Each bootstrap replicate is divided
   by all three μ to produce three CI brackets.

### 2.3 The three μ values

| Tag | Value | Rationale | Reference |
|---|---|---|---|
| `mu_low` | 1.0×10⁻⁹/site/year | Slow molecular clock — long-lived large teleosts (e.g. cichlid radiations have been estimated at 0.7–1.5×10⁻⁹/yr) | typical lower-bound from teleost phylogenomics |
| `mu_mid` | 3.0×10⁻⁹/site/year | Liu et al. 2023 *Mar Life Sci Technol*, electric catfish (*Malapterurus electricus*), r8s-estimated from Siluriformes phylogeny calibrated against teleost fossils | the only catfish-anchored μ in the literature |
| `mu_high` | 9.0×10⁻⁹/site/year | Fast molecular clock — short-generation teleosts (some pufferfish / killifish lineages estimated at 5–10×10⁻⁹/yr) | typical upper-bound |

**The brackets are fixed.** Do not parameterize. Changing them invalidates
all comparisons across catalogue entries; if the field publishes a
catfish-specific μ later, regenerate the JSON with new values and bump
the schema version.

### 2.4 Site-count gate

If `n_4d_sites_used < 200` for an inversion, do NOT emit the
`busco_4d_age_brackets` block. Set it to `null` (or omit). The atlas
empty state will say "Method 3 not computable: too few BUSCO 4D sites
in this inversion (need ≥ 200)."

This matters because small inversions may not contain enough BUSCO
genes (BUSCO gene density ≈ 1 per ~200 kb in catfish), and forcing a
calculation on a handful of 4D sites would produce noise, not signal.

## 3. JSON output schema

Extension to the existing `inversion_age_v1.json` (no schema-version
bump — additive only).

```json
{
  "tool": "inversion_age_v1",
  "schema_version": 1,
  "calibration": {
    "mu_per_site_per_year": 3.0e-9,
    "mu_source_citation": "Liu et al. 2023, Mar Life Sci Technol, doi:10.1007/s42995-023-00197-8",
    "mu_method": "r8s-estimated, Siluriformes phylogeny, fossil-calibrated",
    "mu_brackets_added_in": "STEP_C01f_e_emit_busco_4d_age.py",
    "mu_brackets_rationale": "Catfish-specific mu unknown. Three brackets span plausible teleost range."
  },
  "inversions": [
    {
      "candidate_id": "LG28_cand_01",
      "between_arrangement_age": { ... existing ... },
      "within_arrangement_age":  { ... existing ... },
      "load":                    { ... existing ... },
      "divergence_profile":      { ... existing ... },

      "busco_4d_age_brackets": {
        "n_busco_genes_in_inversion": 14,
        "n_busco_genes_with_4d_sites": 13,
        "n_4d_sites_used": 2847,
        "n_4d_sites_callable_filter_failed": 612,
        "n_bootstrap_replicates": 1000,

        "dxy_4d_between_arrangements": 0.0078,
        "dxy_4d_ci95": [0.0061, 0.0095],

        "mu_low": {
          "value": 1.0e-9,
          "rationale": "slow teleost clock",
          "age_my": 3.90,
          "age_my_ci95": [3.05, 4.75]
        },
        "mu_mid": {
          "value": 3.0e-9,
          "rationale": "Liu 2023 Siluriformes (electric catfish, r8s)",
          "age_my": 1.30,
          "age_my_ci95": [1.02, 1.58]
        },
        "mu_high": {
          "value": 9.0e-9,
          "rationale": "fast teleost clock",
          "age_my": 0.43,
          "age_my_ci95": [0.34, 0.53]
        }
      }
    }
  ]
}
```

When `n_4d_sites_used < 200`:

```json
"busco_4d_age_brackets": null
```

## 4. Atlas display (Page 3 Row C)

Per `SPEC_inversion_age_atlas_surface_AMENDMENT.md` §4.1:

```
TASK 3: absolute age (BUSCO 4D-neutral sites, three-μ bracketing)
──────────────────────────────────────────────────────────────
   μ_low  = 1×10⁻⁹/site/year      3.90 My (3.05–4.75)   slow molecular clock
   μ_mid  = 3×10⁻⁹/site/year      1.30 My (1.02–1.58)   Liu 2023 Siluriformes
   μ_high = 9×10⁻⁹/site/year      0.43 My (0.34–0.53)   fast molecular clock

   n BUSCO genes: 14    n 4D sites: 2,847    bootstrap reps: 1,000
   ⓘ Catfish-specific μ unknown. Three brackets span plausible teleost
   range. The truth is in there. Ranking across inversions robust
   under any μ.
```

Empty state (when `busco_4d_age_brackets === null` or absent):

```
TASK 3: absolute age (BUSCO 4D, three-μ)        [not computed]
──────────────────────────────────────────────────────────────
   Method 3 not computed for this candidate.
   Either BUSCO 4D pipeline hasn't run yet on this chrom, or this
   inversion has too few BUSCO 4D sites (need ≥ 200).
```

## 5. Manuscript-bundle integration

Markdown block (`_bundleMethod3Block(c)`):

```markdown
**Absolute age via BUSCO 4D-neutral sites (three-μ bracketing)**

| μ assumption | Value | Age (My) | 95% CI |
|---|---|---|---|
| slow clock | 1×10⁻⁹/site/year | 3.90 | 3.05–4.75 |
| Siluriformes anchor | 3×10⁻⁹/site/year | 1.30 | 1.02–1.58 |
| fast clock | 9×10⁻⁹/site/year | 0.43 | 0.34–0.53 |

Computed from dXY = 0.0078 across 2,847 BUSCO 4D sites (14 genes,
1,000 site-bootstrap replicates). Catfish-specific μ unknown; three
brackets span the plausible teleost range. Inter-inversion ranking
robust under any μ.
```

TSV scalar columns added to `_candidateTSV`:

- `busco4d_n_4d_sites`
- `busco4d_dxy_4d`
- `busco4d_age_my_under_mu_low`
- `busco4d_age_my_ci95_low_under_mu_low` (and `_high_`)
- `busco4d_age_my_under_mu_mid`
- `busco4d_age_my_ci95_low_under_mu_mid` (and `_high_`)
- `busco4d_age_my_under_mu_high`
- `busco4d_age_my_ci95_low_under_mu_high` (and `_high_`)

All nullable. Empty when the block is null.

## 6. Voice-discipline rules

- **Never quote a single absolute My from Method 3 without the bracket.**
  Forbidden: "the inversion is 1.30 My old". Required: "1.30 My under
  the Liu 2023 anchor; 0.43–3.90 My across plausible μ".
- **Never use Method 3 numbers to drive ranking.** Ranking comes from
  Methods 1 + 2 (which are dimensionless). Method 3 is for absolute
  framing only.
- **Never let reviewers see μ_mid alone.** All three rows must be
  visible together. The display refuses to render one without the
  others.

## 7. Pipeline-side tests

- `STEP_C01f_e_emit_busco_4d_age.py` round-trips a known TSV
- 4D-site identification correct on a hand-checked codon set
- BUSCO-gene overlap filter rejects partial overlaps below 80%
- Site-count gate emits `null` when `n_4d_sites_used < 200`
- Bootstrap CI excludes the point estimate roughly 5% of the time
  (sanity check on the bootstrap)
- All three μ values produce ages such that
  `age_low > age_mid > age_high` (math sanity: lower μ → older)
- Schema validation: every required field present in non-null blocks

## 8. Cross-references

- `SPEC_inversion_age_atlas_surface_AMENDMENT.md` — atlas display layer
- `SPEC_busco_anchors.md` — existing BUSCO usage for synteny anchors
  (different purpose, same BUSCO ortholog set)
- `STEP_C01f_c_burden_regression.R` v9.7 — sibling producer for
  Methods 1 + 2 outputs
- `region_popstats.c` — provides the Hudson dXY implementation; this
  pipeline restricts the input site set to BUSCO 4D and otherwise calls
  the same engine

## 9. What this is NOT

- Not a replacement for the existing Methods 1 + 2. Three methods
  reported alongside each other.
- Not a fastsimcoal2 / Relate / GEVA / msprime simulation-based dater.
  Method 3 is point-estimate dXY divided by 2μ, plus three μ values.
  No demographic model, no coalescent simulator. Cheap to compute.
- Not parameterizable in the atlas. The three μ values are pipeline-side
  constants. If the field publishes a new catfish μ, regenerate the
  JSON LANTA-side; do not let the atlas multiply.
- Not applied to between-species comparisons. The dXY here is between
  *arrangements* of the same species (HOM_REF vs HOM_INV); cross-species
  μ-anchored dating uses different infrastructure (Liu 2023 r8s on
  whole genomes, calibrated to fossils, not at this resolution).
- Not generation-time-dependent. The μ is per-year. Generation time is
  not assumed.

## 10. Sibling pipeline: `STEP_C01f_f_emit_region_stats_per_candidate.py`

Per the AMENDMENT spec §6, the supporting-stats Row D needs a
per-candidate `region_popstats_v1.json`. Pipeline-side this is a
**simple joiner** (~80 lines):

- Read the existing per-window region_popstats output
- For each candidate in the catalogue, slice the windows in
  `[start_bp, end_bp]`
- Emit a per-candidate JSON with the windows array + Tajima shape
  classification

The Tajima-shape classifier is already in C01f_c §4e — re-use that.

Total NEW pipeline work: two small Python scripts (`C01f_e` for
Method 3 + `C01f_f` for region-stats joiner). Both are cheap and
side-effect-free; can run sequentially after `C01f_c`.

End of spec.
