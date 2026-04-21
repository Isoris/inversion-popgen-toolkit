#!/usr/bin/env Rscript

# =============================================================================
# utils/sample_registry.R — COMPATIBILITY SHIM
#
# This file has moved. The canonical location is now:
#     registries/api/R/sample_registry.R
#
# Historical callers that source("utils/sample_registry.R") continue to work
# through this shim, which simply forwards to the canonical file. The exposed
# function `load_registry()` behaves exactly as before (returns the OLD flat
# sample-only API with reg$has_group / reg$add_group / reg$list_groups).
#
# NEW CODE should NOT source this file. Source the full registry loader
# instead, which provides the 4-table nested API (reg$samples, reg$intervals,
# reg$evidence, reg$results, reg$query, reg$compute, reg$ask):
#
#     source("registries/api/R/registry_loader.R")
#     reg <- load_registry()              # full nested object
#     reg$samples$has_group("all_226")    # new API
#     reg$has_group("all_226")            # still works (compat alias)
#
# Migration path: once every caller of utils/sample_registry.R has been
# audited and switched to registry_loader.R (which internally sources the
# canonical sample_registry.R via its fallback chain), this shim can be
# deleted. Target: chat-19 or later.
#
# Chat-18 change: moved from utils/ to registries/api/R/.
# =============================================================================

# Resolve a path to the canonical file. Search order:
#   1. Relative to this shim (works when sourced from codebase root)
#   2. $BASE/registries/api/R/sample_registry.R (works from SLURM)
#   3. Relative to cwd
.shim_find_canonical <- function() {
  # Try self-location first
  here <- tryCatch({
    if (sys.nframe() > 0) normalizePath(sys.frame(1)$ofile, mustWork = FALSE)
    else NULL
  }, error = function(e) NULL)

  cands <- character(0)
  if (!is.null(here) && file.exists(here)) {
    # here is .../utils/sample_registry.R → up one, into registries/api/R
    cands <- c(cands, file.path(dirname(dirname(here)),
                                "registries", "api", "R", "sample_registry.R"))
  }

  base <- Sys.getenv("BASE", "")
  if (nzchar(base)) {
    cands <- c(cands,
      file.path(base, "registries", "api", "R", "sample_registry.R"),
      file.path(base, "inversion-popgen-toolkit", "registries", "api", "R", "sample_registry.R"))
  }

  # Relative to cwd
  cands <- c(cands, "registries/api/R/sample_registry.R",
                    "../registries/api/R/sample_registry.R")

  for (p in cands) {
    if (file.exists(p)) return(p)
  }
  return(NA_character_)
}

.canonical_path <- .shim_find_canonical()
if (is.na(.canonical_path)) {
  stop("[utils/sample_registry.R SHIM] Cannot find canonical file at ",
       "registries/api/R/sample_registry.R. If you moved the repo, set $BASE ",
       "to the repo root before sourcing.")
}

# One-shot deprecation notice (suppress by setting QUIET_SAMPLE_REGISTRY_SHIM=1)
if (!nzchar(Sys.getenv("QUIET_SAMPLE_REGISTRY_SHIM", ""))) {
  message("[utils/sample_registry.R] NOTE: this path is a shim. Canonical: ",
          "registries/api/R/sample_registry.R. ",
          "New scripts should use registries/api/R/registry_loader.R instead.")
}

source(.canonical_path)

# Tidy up shim internals so sourced code sees a clean environment
rm(.shim_find_canonical, .canonical_path)
