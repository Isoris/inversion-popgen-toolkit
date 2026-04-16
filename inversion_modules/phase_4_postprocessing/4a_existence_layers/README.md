# `4a_existence_layers/` — existence evidence (no groups required)

Collects the three existence layers per candidate. Runs first in phase
4. Writes evidence blocks that downstream sub-phases consume, but no
groups are proposed here — group_validation = NONE.

## What runs here

These scripts are **not yet in this delivery**. Upload them here when
ready:

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| `STEP_C01d_candidate_scoring.R` (pass-1) | aggregate evidence score per candidate | `triangle_intervals.tsv.gz` from phase 2d, precomp RDS | scoring_table + Layer A block |
| `STEP_C01e_candidate_figures.R` | per-candidate diagnostic figures | Layer A block | PNG / PDF per candidate |
| `STEP_C01g_boundary_refinement.R` | refine boundaries against Phase 3 output | phase 3 breakpoint validation | boundary.schema.json block |

## Layer sources

Phase 4a integrates evidence from:

| Layer | Source | Produced by |
|---|---|---|
| **A** — genotype covariance | `phase_2_discovery/2c_precomp/` | local PCA → MDS → seeded region-growing |
| **B** — SV concordance | `phase_3_refine/MODULE_5A2_breakpoint_validation/` | DELLY / Manta inversion calls |
| **C** — GHSL haplotype contrast | `phase_2_discovery/2e_ghsl/` *(not yet built)* | Clair3 phased genotypes → per-window sample partitioning |

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
