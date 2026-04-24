# `phase_4_catalog/` — existence evidence (no groups required)

Collects the three existence layers per candidate. Runs first in phase
4. Writes evidence blocks that downstream sub-phases consume, but no
groups are proposed here — group_validation = NONE.

## What runs here

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | aggregate evidence score per candidate; **catalog birth** | staircase blocks from `phase_2/2d_candidate_detection/` (via `triangle_intervals.tsv.gz`), precomp RDS, optionally `--cores_dir`, `--boundary_dir`, `--hyp_dir` | `candidate_scores.tsv.gz` (the catalog), Layer A block |
| `STEP_C01e_candidate_figures.R` | per-candidate diagnostic figures | candidate scoring output | PNG / PDF per candidate |
| `STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` | unified boundary catalog from five sources | PHASE_01C landscape, staircase blocks, seeded regions, SV prior, blue-cross inner boundaries | `boundary_catalog_unified.tsv.gz` |

### Catalog-birth flow

C01d is where the candidate catalog is born. At the moment it runs, it:

1. Reads the staircase blocks from `phase_2/2d_candidate_detection/` (the
   primary boundary detector, via the bridge-converted
   `triangle_intervals.tsv.gz`)
2. Optionally consumes, as evidence dimensions feeding the final score:
   - `--cores_dir` → seeded regions from `phase_2/2c_precomp/STEP_C01b_1`
     (D2 bands, D5 sub-regime)
   - `--boundary_dir` → `boundary_catalog_unified.tsv.gz` from C01g
     (D11 boundary concordance)
   - `--hyp_dir` → hypothesis verdicts from `phase_4/4c/STEP_C01f`
     (D8 peel-or-hypothesis), when available
3. Computes 12 scoring dimensions, assigns Tier 1/2/3/4, writes
   `candidate_scores.tsv.gz`

**Not currently done in this deployment**: per-candidate folder
materialisation (`01_detection/`, `02_genotypes/`, `04_breakpoints/`,
`05_mechanism/`, `06_evolution/`) via a `create_candidate_folders.sh`
— that script exists in older chat transcripts but is not wired into
this deployment. Downstream phase 4b..4e still work because they read
the flat catalog (`candidate_scores.tsv.gz`) directly.

This is the "assignment" step: seeded regions and boundary layers are
consumed as evidence by C01d at the same time it creates the catalog —
there is no separate pre-C01d merge or assignment step. See
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
| `--hyp_dir` | hypothesis verdicts dir | `hypothesis_verdicts.tsv` | `phase_4/4c/STEP_C01f` (second pass) |

The `--cores_dir` reader has a backward-compat fallback that also
accepts the legacy `snake1_core_regions_<chr>.tsv.gz` filename and the
legacy columns `core_family` / `snake_id` / `cheat26_status`. New runs
should use the renamed C01b_1 output; the fallback exists so output
directories produced before the rename still work.

**Removed flags (2026-04-17):** `--flashlight_dir` was previously
documented here but was dead code — C01d gets its SV info from the
staircase scoring table's `sv_overlap_pct` / `n_sv_hits` columns, not
from the SV prior file. The flag is still accepted (silently ignored)
for back-compat with existing launchers, but the row was removed from
this table.

### C01g flag summary

C01g reads `--precomp` (required), `--landscape`, `--staircase`,
`--cores`, `--flashlight`, `--samples`, `--bam_manifest`,
`--mosdepth_dir`, `--repeatmasker`, and `--chrom`.

**Removed flags (2026-04-17):** `--ref_fasta` (never read) and
`--scores` (fossil detection archived — see "Cheat 17 archive" below).

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
| **B** — SV concordance (primary) | `phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R` → `sv_prior/` | DELLY / Manta INV calls parsed into per-chromosome sv_prior |
| **B** — SV concordance (BND rescue) | `phase_3_refine/STEP_B06_bnd_rescue.py` → registry blocks | Paired CT=3to3/INV3 + CT=5to5/INV5 junctions; recovers inversions the strict INV callers missed |
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

## Contract with phase 4b

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
- `boundary.json` (left + right, refined)

## What's NOT here

No group proposals. No Fst, dXY, or theta computed with respect to
groups. No recombinant detection. Those all belong to phase 4b, 4c,
or 4d depending on whether they require validated groups.

## Upload-batch note

When you upload C01d / C01e / C01g, also upload:

- The current `STEP_C01d_candidate_scoring.R` source
- One example `scoring_table_<chr>.tsv` from a real run
- The `triangle_intervals.tsv.gz` output of phase 2d on the same chr

This lets me verify the bridge contract end-to-end and apply any
needed renames.
