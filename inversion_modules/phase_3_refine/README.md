# phase_3_refine — Breakpoint validation + BND rescue + Layer D writer

**Layout note (2026-04-17):** previously nested as
`phase_3_refine/MODULE_5A2_breakpoint_validation/`. The wrapper directory
was unnecessary — phase 3 *is* this module — and was flattened. Scripts
now live directly in `phase_3_refine/`. Comments, echos, and doc headers
that still say "MODULE_5A2" are cosmetic.

## Position in the pipeline (4-layer evidence model)

The inversion catalogue is built from four **independent** evidence layers
that converge per candidate in the evidence registry. Phase 3 is the
writer for **Layer D** (statistical association between PCA genotypes
and physical breakpoint evidence) and for **supplementary Layer B**
(BND-rescued candidates that strict INV callers missed):

| Layer | What | Producer | Registry keys |
|---|---|---|---|
| **A** — PCA sim_mat | Local genotype-covariance structure | phase 2c (C01a) + phase 2d (inv_detect) | `q1_d*`, `q7_layer_a_*` |
| **B** — SV callers (primary) | DELLY2 / Manta INV concordance | phase 2c (`STEP_C00_build_sv_prior.R`) | `q7_layer_b_*` |
| **B** — SV callers (BND rescue) | Paired CT=3to3/INV3 + CT=5to5/INV5 junctions the strict callers missed | **phase 3 STEP06** | `q7b_bnd_*` |
| **C** — GHSL haplotypes | Clair3 phased within-sample haplotype divergence | phase 2e (`STEP_C04_snake3_ghsl_v6.R` + `STEP_C04b_snake3_ghsl_classify.R`) | `q7_layer_c_*` |
| **D** — OR association | Fisher / Armitage test linking PCA groups to read evidence | **phase 3 STEP03** | `q7_layer_d_*` |

The gate that matters is in phase 4c `compute_group_validation()`:
> **if `q7_layer_d_fisher_p < 0.05` AND `q7_layer_d_fisher_or > 5` → group validation = VALIDATED**

No other step writes `q7_layer_d_*`. If phase 3 STEP03 doesn't run for
a candidate, that candidate cannot reach the VALIDATED tier. This is
by design: Layer D is the statistical link between the structural layers
(A/B/C) and the biology (do PCA-INV carriers actually carry the
rearrangement?).

## Why BND rescue exists

At 5–9× coverage DELLY's paired `CT=3to3`+`CT=5to5` rule fails on many
true inversions — the weaker junction misses the threshold and DELLY
emits only the surviving junction as a BND. Manta has the equivalent
with `INV3`/`INV5` flags on pre-conversion BND records. If we stop at
the INV catalogs alone, we miss those inversions. STEP06 rescues them
by:

1. Filtering DELLY BND catalog to `CT∈{3to3,5to5}` (inversion orientation).
2. Filtering Manta's **raw pre-conversion** merged VCF to records whose
   ALT-bracket pattern is `]p]t` (3to3) or `t[p[` (5to5). The
   post-conversion BND catalog is useless here — convertInversion
   already moved all intrachromosomal inversion BNDs into `<INV>`.
3. Greedy pairing left-to-right within 5 Mb on the same chromosome;
   each junction used at most once.
4. Cross-referencing paired positions against the DELLY+Manta INV
   catalogs (±1 kb). Pairs with no catalog match are **orphan rescue
   candidates** — these are the inversions we were missing.
5. Writing `existence_layer_b_bnd_rescue` blocks to the evidence
   registry so phase 4a/4e Layer B scoring sees the rescue evidence.

## Pipeline steps

```
STEP 01 — Extract DELLY+Manta INV calls, match to phase 2d staircase candidates
            → matched_inv_candidates.tsv
STEP 02 — Per-sample pysam evidence extraction at each bp ± 300 bp
            + REF/HET/INV group assignment (Snake2 bands → local PCA fallback)
            → {inv_id}_evidence.tsv + {inv_id}_group_assignments.tsv
STEP 03 — Fisher exact + Chi-square + Cochran–Armitage trend tests + seeds
          + WRITE existence_layer_d block per candidate (Layer D)
            → {inv_id}_tests.tsv + all_candidates_tests.tsv + Layer D registry blocks
STEP 04 — Publication figures (OR forest, contingency table, evidence heatmap)
STEP 05 — DELLY × Manta concordance table + manuscript sentences
STEP 06 — BND orphan rescue (DELLY CT + Manta raw INV3/INV5)
          + WRITE existence_layer_b_bnd_rescue block per paired junction
            → bnd_inv_crossref.tsv + bnd_inv_rescued_candidates.tsv + Layer B rescue blocks
```

## Key outputs

Flat files (for manuscript & inspection):

| File | Description |
|------|-------------|
| `01_per_sample_evidence/{inv_id}_evidence.tsv` | Per-sample support + read-type breakdown |
| `03_statistical_tests/all_candidates_tests.tsv` | All test results across candidates |
| `03_statistical_tests/manuscript_sentences.txt` | Copy-paste text for Methods/Results |
| `04_deconvolution_seeds/{inv_id}_seeds.tsv` | Confirmed REF/HET/INV sample seeds for C01i |
| `04_deconvolution_seeds/seed_summary.tsv` | Which candidates qualify for seeding |
| `05_plots/{inv_id}_validation.pdf` | Per-candidate figure (OR forest + contingency + heatmap) |
| `06_bnd_signal/bnd_inv_paired.tsv` | Paired BND junctions |
| `06_bnd_signal/bnd_inv_rescued_candidates.tsv` | Orphan rescues (no catalog match) |

Registry blocks (for phase 4):

| Block type | Written by | Populates flat keys |
|---|---|---|
| `existence_layer_d` | STEP03 (per candidate) | `q7_layer_d_fisher_or`, `q7_layer_d_fisher_p`, `q7_layer_d_fisher_ci_lower/upper`, `q7_layer_d_armitage_z/p`, `q7_layer_d_n_inv_with_support`, `q7_layer_d_n_inv_total` |
| `existence_layer_b_bnd_rescue` | STEP06 (per paired junction) | `q7b_bnd_rescued`, `q7b_bnd_rescue_source`, `q7b_bnd_pair_bp1/bp2`, `q7b_bnd_pair_size_bp`, `q7b_bnd_match_type`, `q7b_bnd_matched_inv_id`, `q7b_bnd_left/right_pe`, `q7b_bnd_left/right_sr` |

Schemas in `registries/schemas/structured_block_schemas/`:
- `existence_layer_d.schema.json`
- `existence_layer_b_bnd_rescue.schema.json`

## Downstream consumers

- **phase 4c `STEP_C01f_hypothesis_tests.R`** — reads `q7_layer_d_fisher_p`
  and `q7_layer_d_fisher_or` in `compute_group_validation()` for the
  VALIDATED promotion gate (L525–526, also in patches 01/02/03).
- **phase 4c `group_validation_gate.R`** — reads `q7_layer_d_fisher_p`
  (L50) for the OR-pass gate.
- **phase 4e `compute_candidate_status.R`** — reads the full set
  `q7_layer_d_{tested,fisher_or,fisher_p,armitage_z,armitage_p,concordance}`
  (L188–190) for tier classification.
- **phase 4b `STEP_C01i_decompose.R`** — the seed files
  `04_deconvolution_seeds/{inv_id}_seeds.tsv` are intended as seeded
  initialization for k-means clustering. Wiring this replaces unsupervised
  k=3 on PC1 with prior-informed decomposition. **(Not yet wired —
  future work, flagged in chat 5 audit log as deferred.)**

## Registry wiring usage

Registry writes are on by default; phase 3 no-ops gracefully if
`REGISTRIES_ROOT` is unset or points to a non-existent directory.

```bash
# Standard run (registry writes enabled)
cd /scratch/.../inversion_codebase_v8.5/inversion_modules/phase_3_refine
bash run_breakpoint_validation.sh

# Standalone mode (no registry writes, flat TSVs only)
REGISTRIES_ROOT="" bash run_breakpoint_validation.sh

# With candidate_id mapping (recommended when phase 4b/C01i already ran)
CANDIDATE_MAP=/path/to/inv_id_to_candidate_id.tsv \
  bash run_breakpoint_validation.sh
```

## Biological caveats — manuscript-worth flagging

1. **One-sided BNDs can't be rescued.** If only one of the two inversion
   junctions passes BND emission, pairing fails. These inversions are
   invisible to STEP06 and must be validated via phase 2d's PCA signal +
   STEP02's per-sample BAM extraction at the *inferred* breakpoint from
   the PCA candidate boundary.
2. **Repeat-rich breakpoints fail twice.** A breakpoint inside a repeat
   fails INV-typing (lands in BND) AND often fails BND emission too
   (split reads can't cluster). Effectively invisible to both SV callers
   and the rescue pathway.
3. **Fixed inversions (all 226 carriers)** show zero PCA contrast →
   invisible to Layer A → D1 = 0 even when Layers B/C/D all pass.
   Flag as "SV-only, fixed" candidates.
4. **Ambiguous groups.** STEP02 Snake2 band assignments use majority
   vote across overlapping windows. If Snake2 covers <50% of samples,
   STEP02 falls back to local k-means PCA (k=3 on PC1). Both are
   documented per-candidate in `{inv_id}_group_assignments.tsv`.

## Dependencies

- Python 3.9+ with pysam, numpy, matplotlib
- bcftools (STEP01 candidate extraction)
- DELLY2 INV + BND catalogs (from MODULE_4D + MODULE_4F)
- Manta INV PASS + raw pre-conversion catalogs (from MODULE_4H)
- Markdup BAMs (from MODULE_4B)
- Snake 2 band assignments (optional, for STEP02 group source A)
- Dosage files (optional, for STEP02 group source B fallback)
- Registry API at `registries/api/python/registry_loader.py` (optional —
  gracefully skipped if absent, script remains runnable standalone)
