#!/usr/bin/env Rscript

# =============================================================================
# utils/load_bridge.R — COMPATIBILITY SHIM (chat-18)
#
# This file has been renamed. The canonical location is now:
#     utils/registry_bridge.R
#
# Historical callers that source("utils/load_bridge.R") keep working via
# this shim, which simply forwards to the canonical file. The exposed
# symbols (reg, smap, get_Q, get_region_stats, BRIDGE_PATHS, etc.) are
# identical.
#
# NEW CODE should source utils/registry_bridge.R directly. The new name
# reflects what the file actually does (it bridges the registry + ancestry
# engines into a single sourceable entry point). The old name "load_bridge"
# was opaque about its purpose.
#
# Migration target: once every caller has been updated (grep-friendly; the
# old path appears in a small number of scripts + docs), this shim can be
# deleted. Target removal: chat-20 or later.
# =============================================================================

.shim_find_canonical <- function() {
  here <- tryCatch({
    if (sys.nframe() > 0) normalizePath(sys.frame(1)$ofile, mustWork = FALSE)
    else NULL
  }, error = function(e) NULL)

  cands <- character(0)
  if (!is.null(here) && file.exists(here)) {
    # here is .../utils/load_bridge.R → registry_bridge.R in same dir
    cands <- c(cands, file.path(dirname(here), "registry_bridge.R"))
  }

  base <- Sys.getenv("BASE", "")
  if (nzchar(base)) {
    cands <- c(cands,
      file.path(base, "utils", "registry_bridge.R"),
      file.path(base, "inversion-popgen-toolkit", "utils", "registry_bridge.R"))
  }

  cands <- c(cands, "utils/registry_bridge.R", "../utils/registry_bridge.R")

  for (p in cands) {
    if (file.exists(p)) return(p)
  }
  return(NA_character_)
}

.canonical_path <- .shim_find_canonical()
if (is.na(.canonical_path)) {
  stop("[utils/load_bridge.R SHIM] Cannot find canonical file at ",
       "utils/registry_bridge.R. If you moved the repo, set $BASE ",
       "to the repo root before sourcing.")
}

if (!nzchar(Sys.getenv("QUIET_LOAD_BRIDGE_SHIM", ""))) {
  message("[utils/load_bridge.R] NOTE: this path is a shim. Canonical: ",
          "utils/registry_bridge.R. New scripts should source that directly.")
}

source(.canonical_path, local = FALSE)

# Clean up shim temporaries from the sourced environment
rm(.shim_find_canonical, .canonical_path)
