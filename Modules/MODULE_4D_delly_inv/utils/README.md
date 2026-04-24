# `MODULE_4D_delly_inv/utils/`

Shared utilities for SV-caller post-processing.

## `population_regenotype.py` — cross-caller, one canonical copy here

This script lives at `MODULE_4D_delly_inv/utils/` as a pass-15 relocation
from `phase_4_postprocessing/4f_group_dependent/` (where it didn't belong —
it runs before candidate inversions exist, not on them).

**It is not DELLY-specific.** Despite its home, it serves three SV-caller
modules:

| Module              | SV type            | Caller             | FORMAT fields read |
|---------------------|--------------------|--------------------|--------------------|
| `MODULE_4D_delly_inv` | INV              | DELLY2             | DR/DV + RR/RV       |
| `MODULE_4E_delly_bnd` | BND              | DELLY2             | DR/DV + RR/RV       |
| `MODULE_4G_manta_sv`  | INV (and others) | Manta              | PR + SR             |

The corresponding SLURM launcher at `MODULE_4D_delly_inv/slurm/SLURM_A03b_population_regenotype.sh`
auto-detects which module it was invoked from by looking for
`00_module4d_config.sh`, `00_module4e_config.sh`, or `00_module4g_config.sh`
in the parent directory, and sets `CALLER=delly` or `CALLER=manta`
accordingly. To run it for MODULE_4E or MODULE_4G, copy the SLURM script
into that module's `slurm/` directory — the script will find this canonical
Python copy via its `../utils/population_regenotype.py` fallback.

## Pipeline slot

Sits between `SLURM_A03_merge_genotype.sh` (produces the merged cohort VCF
with naive per-sample GTs) and `SLURM_A04_annotation_layers.sh` (which
would otherwise use the dropout-biased GTs for strict filtering).

See the script's own docstring for the full method description (Bayesian
population-prior rescue of low-coverage dropout at ~9× WGS).

## `plot_INV_results.R` + `plot_INV_results_v1.R`

Standard MODULE_4D post-processing plots. DELLY-specific, INV-specific.
