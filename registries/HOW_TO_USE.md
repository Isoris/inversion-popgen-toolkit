# How to use the registry

**The single line you need at the top of every R script:**

```r
source("utils/registry_bridge.R")
```

That's it. After sourcing, every script has:

- `reg` — the full 4-table registry (samples / intervals / evidence / results)
- `smap` — Ind↔CGA sample-name translator
- `get_Q()`, `get_Q_summary()` — live ancestry Q matrices
- `get_region_stats()` — unified popstats dispatcher (FST, dxy, theta, MI)
- `BRIDGE_PATHS` — named list of resolved absolute paths

Idempotent: safe to source multiple times. If a script is called from
another script that already sourced the bridge, the second source is a
no-op.

---

## Why there's a separate `utils/registry_bridge.R` instead of `registries/load.R`

It's historical. The bridge file used to be called `utils/load_bridge.R`
and was thought of as "pipeline glue" before the 4-table registry
existed. Chat-18 renamed it to `utils/registry_bridge.R` once it became
obvious that the registry IS the main thing it loads, but the file stays
in `utils/` because that's where other shared R libraries live
(`sample_map.R`, `lib_ghsl_panel.R`, `theme_systems_plate.R`).

The old name `utils/load_bridge.R` still works as a 3-line shim.
Suppress the deprecation notice with `QUIET_LOAD_BRIDGE_SHIM=1`.

---

## Short examples

### Read karyotype groups for a candidate

```r
source("utils/registry_bridge.R")

cid <- "LG12_17"
grs <- reg$samples$get_groups_for_candidate(cid)
ref_samples <- grs$HOM_REF
het_samples <- grs$HET
inv_samples <- grs$HOM_INV
```

### Check what's already been computed

```r
source("utils/registry_bridge.R")

have <- reg$results$ask_what_for_candidate("LG12_17")
# data.table with one row per existing result; empty if nothing computed yet
```

### Compute pairwise FST if missing (compute-if-missing pattern)

```r
source("utils/registry_bridge.R")

for (cid in reg$evidence$list_candidates()) {
  existing <- reg$results$ask_what_for_candidate(cid)
  if (!nrow(existing[kind == "candidate_q" & K == 8])) {
    reg$compute$ancestry_q_and_f_for_candidate(cid, K = 8, persist = TRUE)
  }
}
```

### Prove reproducibility before submitting

```r
source("utils/registry_bridge.R")
ic <- reg$results$integrity_check()
stopifnot(attr(ic, "all_pass"))
# All 7 checks passed — ready to zip and upload
```

---

## SLURM / bash side

```bash
source utils/pipeline_bridge.sh

# $REGISTRY_BRIDGE is now set to the correct absolute path
Rscript -e "source('${REGISTRY_BRIDGE}'); reg\$status()"
```

The old env var `$LOAD_BRIDGE` is preserved as an alias for back-compat.

---

## Full API reference

See `API_CHEATSHEET.md` in this folder for the complete method listing,
`DATABASE_DESIGN.md` for the schema and design rationale.
