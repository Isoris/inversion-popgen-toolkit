# PATH_WIRING.md — DELLY and Manta catalog paths across modules

## Problem

Two scripts in this codebase parse the same SV catalog VCFs with different
hardcoded default paths. Environment variables override the defaults, so
whether a given script gets the right data depends entirely on what the
launcher exports. The defaults themselves are inconsistent.

## What each script expects

### `STEP_C00_build_sv_prior.R` defaults (old "delly_sv_*" names)

```
DELLY_INV_VCF  $BASE/delly_sv_INV/07_final_catalogs/catalog_226.INV.vcf.gz
DELLY_DEL_VCF  $BASE/delly_sv/07_final_catalogs/catalog_226.DEL.vcf.gz
DELLY_BND_VCF  $BASE/delly_sv_BND/07_final_catalogs/catalog_226.BND.vcf.gz
MANTA_INV_VCF  $BASE/MODULE_4H_ALL_Manta/05_final_catalogs/catalog_226.INV.PASS.vcf.gz
MANTA_DEL_VCF  $BASE/MODULE_4H_ALL_Manta/05_final_catalogs/catalog_226.DEL.PASS.vcf.gz
POP_CONF_DIR   $BASE/MODULE_4H_ALL_Manta/10_population_confidence
```

### `phase_3_refine` config defaults (new MODULE_4* names)

```
DELLY_BASE      $BASE/MODULE_4B_DEL_Delly
DELLY_INV_BASE  $BASE/MODULE_4D_INV_Delly
DELLY_BND_BASE  $BASE/MODULE_4F_BND_Delly
MANTA_BASE      $BASE/MODULE_4H_ALL_Manta

DELLY_INV_VCF   $DELLY_INV_BASE/07_final_catalogs/catalog_226.INV.vcf.gz
DELLY_BND_VCF   $DELLY_BND_BASE/07_final_catalogs/catalog_226.BND.vcf.gz
MANTA_INV_VCF   $MANTA_BASE/05_final_catalogs/catalog_226.INV.PASS.vcf.gz
```

## Why this is fragile

These are the same data in both cases. If both directory names exist on
LANTA (via symlinks), everything works silently. If one is real and the
other is a broken/missing symlink, whichever script uses that default
will fail. Worse, whichever script is not failing could be reading stale
catalog data from an older run of the SV-calling pipeline without
anyone noticing.

## Proposed fix

Add the four DELLY/Manta base paths to the master inversion config, then
have every script use those variables. One canonical place, no defaults
scattered across scripts.

### 1. Edit `00_inversion_config.sh`

Add this block near the "Stable external inputs" section, before the
"Results subdirectories" section:

```diff
 # ── Stable external inputs ───────────────────────────────────────────────────
 HETDIR="${BASE}/het_roh"
 REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
 REF_FAI="${REF}.fai"
 BAMLIST="${HETDIR}/01_inputs_check/bamlist_qcpass.txt"
 SAMPLES_IND="${HETDIR}/01_inputs_check/samples.ind"

+# ── SV caller module bases (canonical) ──────────────────────────────────────
+# Every script that reads a SV catalog VCF derives its path from these.
+# Do not hardcode delly_sv_* or MODULE_4B/4D/4F/4H paths in any other script.
+DELLY_DEL_BASE="${BASE}/MODULE_4B_DEL_Delly"
+DELLY_INV_BASE="${BASE}/MODULE_4D_INV_Delly"
+DELLY_BND_BASE="${BASE}/MODULE_4F_BND_Delly"
+MANTA_BASE="${BASE}/MODULE_4H_ALL_Manta"
+
+# Canonical VCFs (override via environment if running against alternate builds).
+export DELLY_INV_VCF="${DELLY_INV_BASE}/07_final_catalogs/catalog_226.INV.vcf.gz"
+export DELLY_DEL_VCF="${DELLY_DEL_BASE}/07_final_catalogs/catalog_226.DEL.vcf.gz"
+export DELLY_BND_VCF="${DELLY_BND_BASE}/07_final_catalogs/catalog_226.BND.vcf.gz"
+export MANTA_INV_VCF="${MANTA_BASE}/05_final_catalogs/catalog_226.INV.PASS.vcf.gz"
+export MANTA_DEL_VCF="${MANTA_BASE}/05_final_catalogs/catalog_226.DEL.PASS.vcf.gz"
+export MANTA_BND_POST_VCF="${MANTA_BASE}/05_final_catalogs/catalog_226.BND.PASS.vcf.gz"
+export MANTA_RAW_MERGED_VCF="${MANTA_BASE}/02_merged_cohort/cohort_226.ALL.raw.vcf.gz"
+export POP_CONF_DIR="${MANTA_BASE}/10_population_confidence"
```

### 2. Update `phase_3_refine/00_phase3_config.sh`

Remove the duplicated `DELLY_INV_VCF` etc. definitions; keep only the
module-specific ones (like `DELLY_INV_BED`, `DELLY_INV_GT`,
`MANTA_BND_POST_VCF` — anything not already in the master). Inherit
everything else from the master via the existing `source` line.

### 3. Update `STEP_C00_build_sv_prior.R`

Drop the `Sys.getenv()` defaults since the env is guaranteed to be set
by the master config. Keep the `Sys.getenv()` lookups themselves (they
still provide the abstraction layer and allow per-run overrides):

```r
# C00 line 97 and onward, after:
# SV catalogs
DELLY_INV_VCF <- Sys.getenv("DELLY_INV_VCF")
DELLY_DEL_VCF <- Sys.getenv("DELLY_DEL_VCF")
MANTA_INV_VCF <- Sys.getenv("MANTA_INV_VCF")
MANTA_DEL_VCF <- Sys.getenv("MANTA_DEL_VCF")
DELLY_BND_VCF <- Sys.getenv("DELLY_BND_VCF")
POP_CONF_DIR  <- Sys.getenv("POP_CONF_DIR")
for (v in c("DELLY_INV_VCF","DELLY_DEL_VCF","DELLY_BND_VCF",
            "MANTA_INV_VCF","MANTA_DEL_VCF","POP_CONF_DIR")) {
  if (!nzchar(Sys.getenv(v))) {
    stop("[C00] Env var ", v, " not set — did you source 00_inversion_config.sh?")
  }
}
```

This makes it impossible to silently use a stale default.

## Which directory name is canonical?

**This needs to be confirmed against LANTA.** The old (`delly_sv_INV`,
`delly_sv`, `delly_sv_BND`) names predate the `MODULE_4B/4D/4F`
renumbering. Right now:

| Module | New name (canonical by convention) | Old name (may be symlink) |
|---|---|---|
| DELLY DEL | `MODULE_4B_DEL_Delly` | `delly_sv` |
| DELLY INV | `MODULE_4D_INV_Delly` | `delly_sv_INV` |
| DELLY BND | `MODULE_4F_BND_Delly` | `delly_sv_BND` |
| Manta (all SV types) | `MODULE_4H_ALL_Manta` | (no old name — Manta was always under MODULE_4H) |

Run this on LANTA to verify which is real and which is a symlink:

```bash
for d in delly_sv delly_sv_INV delly_sv_BND MODULE_4B_DEL_Delly MODULE_4D_INV_Delly MODULE_4F_BND_Delly MODULE_4H_ALL_Manta; do
  path="$BASE/$d"
  if [ -d "$path" ]; then
    echo "$(readlink -f "$path")  <=  $path"
  else
    echo "MISSING: $path"
  fi
done
```

If the old names are symlinks pointing at the MODULE_4* names, no cleanup
needed — just adopt the MODULE_4* naming in the master config as shown
above, and the old paths continue to resolve. If the old names are the
real directories, reverse the symlinks so `MODULE_4*` is the real name
(the newer convention).

## Precedent elsewhere in the pipeline

The master inversion config already exports several "stable external
inputs" this way (`BAMLIST`, `SAMPLES_IND`, `REF`). The SV catalog VCFs
are stable external inputs to phase 2 and phase 3. They belong in the
same block.

## What this changes

- One place to update when SV catalog paths move
- No silent divergence between C00 and breakpoint_validation
- Launchers no longer need to `export` SV paths individually
- A script running against a module that has not been built yet gets a
  clear error instead of reading a stale or missing catalog
