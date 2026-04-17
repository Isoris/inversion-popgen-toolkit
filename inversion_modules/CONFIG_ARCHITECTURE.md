# CONFIG_ARCHITECTURE.md — how configuration works across the pipeline

## TL;DR

**Keep `00_inversion_config.sh` where it is, inside the inversion codebase.
Do NOT add a separate repo-root master config.** Per-module configs (like
`00_breakpoint_validation_config.sh`) source the inversion config as their
parent. This is a proven pattern in the codebase. A repo-root master config
would add a layer that does not pay for itself.

## The existing pattern

The codebase already has a clean two-tier pattern:

```
inversion_codebase_v8.5/
├── 00_inversion_config.sh                         # MASTER for anything inversion
│
├── phase_3_refine/
│   └── 00_breakpoint_validation_config.sh         # sources master, adds module-specific
│                                                  # (flat layout — was MODULE_5A2_breakpoint_validation/)
│
├── utils/
│   ├── load_bridge.R
│   ├── sample_map.R
│   ├── sample_registry.R
│   └── theme_systems_plate.R
│
├── phase_2_discovery/
│   ├── 2c_precomp/STEP_C00_build_sv_prior.R       # sources master via env
│   └── ...
│
└── ...
```

A module's config does two things, in this order:

1. Source the master inversion config to inherit all project-wide paths and
   variables (`BASE`, `INVDIR`, `SAMPLES_IND`, `RSCRIPT_BIN`, helpers like
   `inv_init_dirs`, ...).
2. Define everything module-specific (per-module output dirs, per-module
   parameters, per-module helper functions, per-module SBATCH defaults).

The `breakpoint_validation` config demonstrates this exactly:

```bash
# 00_breakpoint_validation_config.sh, lines 37-42
SCRIPT_DIR_BPV="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INV_CONFIG="${SCRIPT_DIR_BPV}/../00_inversion_config.sh"
[[ -f "${INV_CONFIG}" ]] || INV_CONFIG="${BASE:-/scratch/.../inversion_codebase_v8.5}/00_inversion_config.sh"
[[ -f "${INV_CONFIG}" ]] || { echo "ERROR: Cannot find 00_inversion_config.sh" >&2; exit 1; }
source "${INV_CONFIG}"
```

Then it adds its own variables: `BPV_ROOT`, `BPV_EVIDENCE`, `MIN_DELLY_QUAL`,
`bpv_log`, `bpv_init_dirs`, and so on.

## What the master config does

`00_inversion_config.sh` is project-wide. It defines:

- Project root (`BASE`), results root (`INVDIR`), codebase root (`CODEBASE`)
- Stable external inputs: reference FASTA (`REF`), sample list
  (`SAMPLES_IND`), BAM list (`BAMLIST`), ancestry config (`ANCESTRY_CONFIG`)
- Results subdirectories for all phases (`DOSAGE_DIR`, `LPCA_DIR`, `MDS_DIR`,
  `BRIDGE_DIR`, `REGIONAL_DIR`, `FOLLOWUP_RESULTS`, ...)
- Shared defaults (`NPC`, `WINSIZE`, `WINSTEP`, `MDS_MODE_DEFAULT`,
  `MDS_SEED_DEFAULT`, `Z_THRESH_DEFAULT`, ...)
- Shared helpers (`inv_log`, `inv_err`, `inv_die`, `inv_check_file`,
  `inv_init_dirs`)
- Per-module directory aliases when a module's output dir belongs in the
  "results" tree (`SNAKE1_DIR`, `SNAKE2_DIR`, etc.)

## What a module config adds

A per-module config:

- Sources the master
- Defines module output root (typically under `INVDIR/<module_number>_<name>/`)
- Defines module-specific input paths (caller VCFs for that module, etc.)
- Defines module-specific parameters (filtering thresholds, seed criteria,
  ...)
- Defines module-specific helpers (`<prefix>_log`, `<prefix>_init_dirs`,
  ...)
- May override defaults for the specific module (rare — only when the
  module genuinely needs a different value)

## Scripts source the config, not their launchers

R scripts that need config values read them from the environment:

```r
BASE <- Sys.getenv("BASE", "/scratch/.../Quentin_project_KEEP_2026-02-04")
INVDIR <- Sys.getenv("INVDIR", file.path(BASE, "inversion_localpca_v7"))
```

with sensible fallbacks when the env is not set (for interactive runs).
Launchers handle the `source` + `set -a` so that the subprocess
inherits everything:

```bash
set -a
source "${CONFIG}"
set +a
```

This is the pattern every launcher already uses. Do not change it.

## Why NOT a repo-root master config

A repo-root `config.sh` that source-chains into the inversion config has
been suggested. It would not pay for itself:

1. **The inversion pipeline is the repo.** The "toolkit" aspect is about
   the inversion pipeline as a unit, not a collection of unrelated tools
   that happen to share a config. There is no cross-project config to
   centralise.
2. **The master is already at the top of the pipeline.**
   `00_inversion_config.sh` sits at the root of the codebase directory.
   Adding a layer above it just means an extra `source` hop for nothing.
3. **It introduces a circular lookup.** If the repo-root config points at
   the inversion config, and the inversion config defines `BASE` which
   points at a LANTA path separate from the repo, then the repo-root
   config needs its own copy of `BASE` (because the inversion config has
   not been sourced yet). That duplication is worse than the problem it
   solves.
4. **Per-module configs already handle it.** If a module ever needs a
   different project root, its config can override `BASE` before sourcing
   the master. That is one line in one file.

A repo-root doc (**this file**) pointing new developers at the existing
pattern is strictly better than adding a config file that has to be
maintained.

## Recommended checklist for a new module

When adding a new module (say `MODULE_6X_new_module/`):

1. Create `MODULE_6X_new_module/00_<module_name>_config.sh` that sources
   the master via the relative path.
2. Define module-specific paths: `MODX_ROOT`, `MODX_<subdir>_DIR`, etc.
3. Define module-specific parameters.
4. Define `modx_log`, `modx_err`, `modx_die`, `modx_init_dirs` helpers.
5. Each launcher in the module sources the module config (which in turn
   sources the master):

   ```bash
   CONFIG="${SCRIPT_DIR}/00_<module_name>_config.sh"
   set -a; source "${CONFIG}"; set +a
   modx_init_dirs
   ```

6. Each script inside the module reads from env with fallbacks.

## When to add a variable to the master vs the module config

**Master** if the variable describes a project-wide fact that multiple
modules need:

- Paths to stable inputs (`REF`, `SAMPLES_IND`, `BAMLIST`)
- Project-wide output roots (`INVDIR`, `DOSAGE_DIR`, `MDS_DIR`)
- Shared defaults for algorithms that run across modules
  (`MDS_MODE_DEFAULT`, `NPC`, `WINSIZE`)

**Module** if the variable is specific to that module:

- SV-caller paths that only that module uses
- Thresholds that only make sense in that module's context
- Output subdirectories specific to that module

When in doubt, put it in the module first. You can promote it to the
master later if a second module needs it. Demoting is harder because it
requires updating every reference.

## Current state

Already compliant with this pattern:

- `00_inversion_config.sh` — master
- `phase_3_refine/00_breakpoint_validation_config.sh` —
  module config, sources master

Suggested next modules to add a config for:

- `MODULE_4B_DEL_Delly/` — `00_delly_del_config.sh`
  (currently sources master implicitly via env vars set at runtime)
- `MODULE_4D_INV_Delly/` — `00_delly_inv_config.sh`
- `MODULE_4F_BND_Delly/` — `00_delly_bnd_config.sh`
- `MODULE_4H_ALL_Manta/` — `00_manta_config.sh`

Each of these would centralise the per-caller VCF paths, genotype-field
expectations (DELLY DR/DV/RR/RV vs Manta PR/SR), and caller-specific
filter thresholds.

## See also

- `00_inversion_config.sh` — master project config
- `phase_3_refine/00_breakpoint_validation_config.sh` —
  reference implementation of a module config
- `phase_2_discovery/*/README.md` — per-phase documentation
- `phase_2_discovery/2c_precomp/RENAMING.md` — terminology migration
