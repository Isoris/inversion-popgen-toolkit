# `phase_4_catalog/` — catalog birth + per-candidate blocks

This is where the candidate catalog is born. Phase 4 consumes phase_2
discovery output (staircase, precomp, boundaries, SV prior, seeded
regions) and produces `candidate_scores.tsv.gz` — the single source of
truth for every downstream phase. It also writes four structured
registry blocks per candidate (regime, local-structure, distance,
boundary-refined) that feed phase_9 characterization.

No groups are proposed here — `group_validation = NONE` for every
script in this folder. Group proposal happens in phase_7.

## Pipeline DAG

```
phase_2_discovery/01c_landscape      → boundary_catalog_<chr>.tsv.gz
phase_2_discovery/2d_candidate_det   → scoring_table_<chr>.tsv
phase_2_discovery/2c_precomp         → <chr>.precomp.rds (C01a)
                                     → seeded_regions_summary (C01b_1)
                                     → sv_prior_<chr>.rds (C00)
                            │
                            ▼
                ┌──── STEP_C01g_boundary_catalog
                │       │
                │       └──> boundary_catalog_unified.tsv.gz
                │            (feeds C01d D11)
                │
                ▼
         STEP_C01d_candidate_scoring     ← CATALOG BIRTH
                │
                └──> candidate_scores.tsv.gz
                     (12 scoring dims + Tier 1-4 + popgen annotation)
                     = THE CATALOG for all downstream phases
                            │
        ┌───────────────────┼───────────────────┬───────────────────┐
        ▼                   ▼                   ▼                   ▼
 STEP_C01e           STEP_C01j           STEP_C01l           STEP_C01m
 candidate_figures   regime_              local_              distance_
 (diagnostic PDFs;   compatibility_      structure_           concordance
  no registry        engine              segments            (per-chrom
  block)             │                   │                    blocks)
                     ▼                   ▼                    ▼
                  registry:          registry:             registry:
                  regime_segments    local_structure_      distance_
                  block              segments block        concordance
                                                           block
                            │
                            ▼
               phase_5_qc_triage (QC triage, q_qc_shelf_*)
                            │
                            ▼
               phase_6_breakpoint_refinement (bp bp1/bp2 consensus)
                            │
                            ▼
               phase_7_karyotype_groups (proposal + validation)
                            │
                            ▼
               phase_8_evidence_biology (Q4/Q5/Q6/Q7 mechanism, age, robustness)
                            │
                            ▼
               phase_9_classification (characterize_candidate.R)
```

C01g runs before C01d; everything else runs after C01d and can parallelize.

## What runs here

| Step | File | Role | Writes |
|---|---|---|---|
| C01g | `STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` | **Boundary unification.** Merges five boundary sources (01C landscape, staircase, seeded regions, SV prior, blue-cross inners) by proximity and runs cheats 4/10/11 (Fst step, depth anomaly, clipped-read pileup) on each. | `boundary_catalog_unified.tsv.gz`, `boundary_concordance.tsv.gz`, `boundary_summary.tsv` |
| C01d | `STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | **Catalog birth.** 12-dimensional scoring against staircase scoring tables + precomp + optional boundary/hypothesis input. Assigns Tier 1-4 per candidate. Adds popgen annotation (Fst/θπ/θW/D) via Engine B. | `candidate_scores.tsv.gz` (THE catalog), `boundary_refined` registry block, `q1_layer_a_*` flat keys |
| C01j | `STEP_C01j_regime_compatibility_engine.R` | **Compatibility-based regime engine.** At each sliding window, groups samples into compatible site-vector classes (inspired by hifiasm). Detects regime transitions = boundary candidates. Labels segments as background_soup / structured_inversion / nested_regime / recombinant_zone / transition. | `regime_segments` registry block (per-candidate) |
| C01l | `STEP_C01l_local_structure_segments.R` | **Local structure quantification.** For each Tier 1/2 candidate, measures ancestry concentration across 5 spatial segments (left flank / inversion halves / core / right flank). Uses dosage (not NGSadmix Q) for local k-means + PCA → Δ₁₂, entropy, ENA per segment. | `local_structure_segments` registry block (per-candidate) |
| C01m | `STEP_C01m_distance_concordance.R` | **Multi-scale concordance matrix.** Sample × sample co-grouping at distances d = {80, 160, 320, 640} windows. Family LD decays with d; inversions persist → contrast matrix highlights inversion-supporting pairs. | `distance_concordance` registry block (per-chromosome) |
| C01e | `STEP_C01e_candidate_figures.R` | **Per-candidate manuscript figures.** Diagnostic PDFs only — no registry writes. Panels A-F: ideogram / regional PCA / genotype heatmap / popgen signals / breakpoint evidence / gene-repeat context. | `<outdir>/candidate_<chr>_<start>_<end>/panel_*.png` + `data/*.tsv` |

C01d is where the catalog comes to exist. C01g must run first (its
output feeds C01d's D11 boundary-concordance dimension). C01e/j/l/m
read the catalog and can run in parallel.

## Catalog-birth flow

At the moment C01d runs, it:

1. Reads the staircase blocks from `phase_2/2d_candidate_detection/`
   (the primary boundary detector, via the bridge-converted
   `triangle_intervals.tsv.gz`).
2. Optionally consumes, as evidence dimensions feeding the final score:
   - `--cores_dir` → seeded regions from `phase_2/2c_precomp/STEP_C01b_1`
     (populates D2 band + D5 sub-regime)
   - `--boundary_dir` → `boundary_catalog_unified.tsv.gz` from C01g
     (populates D11 boundary concordance)
   - `--hyp_dir` → hypothesis verdicts from
     `phase_7/validation/STEP_C01f_hypothesis_tests` (populates D8
     peel-or-hypothesis; only available in 2-pass mode)
3. Computes 12 scoring dimensions, assigns Tier 1/2/3/4, writes
   `candidate_scores.tsv.gz` and the `boundary_refined` registry block.

Downstream phase_5..9 all read the flat catalog
(`candidate_scores.tsv.gz`) directly. There is no separate pre-C01d
merge or assignment step — seeded regions and boundary layers are
consumed by C01d at the same time it creates the catalog. See
`phase_2/2c_precomp/README.md` "Architectural note" for why the old
merge step was dropped.

### C01d input contract

| CLI flag | Points at | File C01d reads | From |
|---|---|---|---|
| positional `<detector_dir>` | staircase output dir | `scoring_table_<chr>.tsv` (+ block/peel/consensus tables) | `phase_2/2d_candidate_detection/` |
| positional `<outdir>` | where to write | `candidate_scores.tsv.gz` + per-candidate folders | (output) |
| `--precomp_dir` | C01a precompute cache | `<precomp_dir>/precomp/<chr>.precomp.rds` | `phase_2/2c_precomp/STEP_C01a` |
| `--cores_dir` | seeded-region output dir | `seeded_regions_summary_<chr>.tsv.gz` | `phase_2/2c_precomp/STEP_C01b_1` |
| `--boundary_dir` | unified boundary catalog dir | `boundary_catalog_unified.tsv.gz` | `STEP_C01g` in this folder |
| `--hyp_dir` | hypothesis verdicts dir | `hypothesis_verdicts.tsv` | `phase_7/validation/STEP_C01f` (second pass) |

The `--cores_dir` reader has a backward-compat fallback that also
accepts the legacy `snake1_core_regions_<chr>.tsv.gz` filename and the
legacy columns `core_family` / `snake_id` / `cheat26_status`. New runs
should use the renamed C01b_1 output; the fallback exists so output
directories produced before the rename still work.

### C01g flag summary

C01g reads `--precomp` (required), `--landscape`, `--staircase`,
`--cores`, `--sv_prior`, `--samples`, `--bam_manifest`, `--mosdepth_dir`,
`--repeatmasker`, and `--chrom`.

### Deprecated flags (accepted with warning)

Pass 22 (2026-04-24) replaced silent-ignore with deprecation warnings
for three obsolete flags. They're still accepted for back-compat with
existing launcher scripts, but running them prints a one-line notice:

| Script | Flag | Status | Notice |
|---|---|---|---|
| C01d | `--sv_prior_dir` | Accepted + ignored + warns | D7 (sv_breakpoint) reads `sv_overlap_pct` / `n_sv_hits` from the staircase scoring_table directly (populated by `phase_2/2d/STEP_D06_sv_overlap.R`). Safe to drop. |
| C01g | `--ref_fasta` | Accepted + ignored + warns | Never read anywhere. Safe to drop. |
| C01g | `--scores` | Accepted + ignored + warns | Fed the Cheat 17 fossil detection archived on 2026-04-17 to `_archive_superseded/cheat17_fossil_detection/`. Safe to drop. |

### Cheat 17 fossil detection — archived 2026-04-17

C01g used to run a "Cheat 17" block that classified boundaries outside
any candidate as FOSSIL_CANDIDATE / MISSED_INVERSION / FALSE_BOUNDARY.
Implementation required reading C01d's `candidate_scores.tsv.gz` via
`--scores` — but C01d runs *after* C01g in the canonical flow (C01d
reads C01g's `boundary_catalog_unified.tsv.gz` for D11). Getting both
blocks fully populated required a two-pass orchestration that was
never wired, so in every single-pass run Cheat 17 silently produced
nothing.

The feature has been archived to `_archive_superseded/cheat17_fossil_detection/`.
Downstream columns `cheat17_class` and `cheat17_inv_likeness` are
still emitted as all-NA (for schema stability); `boundary_activity`
is now `"active"` unconditionally. If the fossil distinction ever
becomes manuscript-relevant, the archive README describes how to bring
it back as a post-hoc annotation script.

## Layer sources

Phase 4a integrates evidence from:

| Layer | Source | Produced by |
|---|---|---|
| **A** — genotype covariance | `phase_2_discovery/2c_precomp/` | local PCA → MDS → seeded region-growing |
| **B** — SV concordance (primary) | `phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R` → `sv_prior/` | DELLY / Manta INV calls parsed into per-chromosome `sv_prior` |
| **B** — SV concordance (BND rescue) | `phase_3_refine/STEP_B06_bnd_rescue.py` → registry blocks | Paired CT=3to3/INV3 + CT=5to5/INV5 junctions; recovers inversions the strict INV callers missed (see `Modules/MODULE_4E_delly_bnd/WIKI.md` for CT semantics) |
| **C** — GHSL haplotype contrast | `phase_2_discovery/2e_ghsl/` | Clair3 phased genotypes → per-window sample partitioning |
| **D** — OR association | `phase_3_refine/STEP_D03_statistical_tests_and_seeds.py` → registry blocks | Fisher + Armitage test linking PCA groups to physical breakpoint evidence |

> **NOTE 2026-04-17 (chat 5 audit):** the Layer B/D attribution was
> corrected in this session. An earlier revision of this README claimed
> `phase_3_refine/` fed Layer B directly as `DELLY / Manta inversion
> calls` — it doesn't; that's C00's job. Phase 3 contributes Layer B
> **supplementary evidence** (BND rescue via STEP_B06) and owns Layer D
> entirely (OR test via STEP_D03). Both phase 3 contributions reach phase
> 4a/4e via registry blocks (`existence_layer_d`,
> `existence_layer_b_bnd_rescue`) whose `keys_extracted` directives
> auto-materialise `q7_layer_d_*` and `q7b_bnd_*` flat keys.

If any layer is missing, its score is NA and the combined
`existence_score` is computed from the layers that ARE available.

## Contract with phase 4b (and subsequent phases)

Phase 4a must write these flat keys into the evidence registry before
phase 4b starts:

- `q1_layer_a_pass` (boolean)
- `q1_layer_b_pass` (boolean)
- `q1_layer_c_pass` (boolean)
- `q1_existence_score` (0-100)

And the corresponding Tier-2 blocks:

- `existence_layer_a.json`
- `existence_layer_b.json`
- `existence_layer_c.json`
- `boundary.json` / `boundary_refined.json`

Per-candidate blocks produced by C01j/C01l/C01m are consumed by phase_9
characterization; they don't gate phase_5 entry.

## What's NOT here

- **No group proposals.** No Fst/dXY/theta computed with respect to
  groups. No recombinant detection. Those all belong to phase_7
  (proposal / validation) or phase_8 (evidence biology) depending on
  whether they require validated groups.
- **No per-candidate folder materialisation.** The `create_candidate_folders.sh`
  pattern (`01_detection/`, `02_genotypes/`, ...) from older chat
  transcripts is not wired in this deployment. Downstream phases 5..9
  read the flat catalog directly.

## Launchers

- `LAUNCH_C01g_boundary.sh` — runs C01g before C01d in the canonical flow
- `LAUNCH_C01d_scoring_pass1.sh` — runs C01d cold (no `--hyp_dir`)
- `../phase_9_classification/LAUNCH_C01d_scoring_pass2.sh` — reruns C01d
  with `--hyp_dir` populated (after phase_7/validation/C01f)

C01e/j/l/m launchers are not provided in phase_4_catalog itself; they
run from `LAUNCH_characterize_classify.sh` in phase_9 or from
per-candidate drivers.

## See also

- `../phase_2_discovery/README.md` — upstream data sources
- `../phase_7_karyotype_groups/` — where proposed groups get validated
- `../phase_9_classification/SPEC_VS_REALITY.md` — the 317-key spec and
  which phase_4 scripts contribute to each question
- `../../Modules/MODULE_4E_delly_bnd/WIKI.md` — DELLY BND CT semantics,
  relevant to Layer B BND rescue via phase_3/STEP_B06
- `../../docs/PHASE4_ARCHITECTURE.md` — full architectural notes and
  historical rationale
