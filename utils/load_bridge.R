#!/usr/bin/env Rscript
# =============================================================================
# load_bridge.R — Universal loader for cross-module wiring
#
# Source this at the top of any R script. It provides:
#   - smap:  sample_map object (Ind↔CGA conversion)
#   - reg:   sample_registry object (group management)
#   - get_Q, get_Q_summary, get_region_stats (from unified ancestry, if avail)
#   - BRIDGE_PATHS: named list of resolved paths
#
# Design:
#   - Auto-detects its own location → resolves all paths relative to BASE
#   - Reads env vars set by pipeline_bridge.sh (SLURM mode)
#   - Falls back to scanning known locations if env vars are unset
#   - Idempotent: safe to source multiple times
#
# Usage in scripts:
#   source("utils/load_bridge.R")           # from codebase root
#   — or —
#   source(Sys.getenv("LOAD_BRIDGE"))       # from SLURM (set by pipeline_bridge.sh)
#
# After sourcing:
#   smap$to_real("Ind0")                    # → "CGA009"
#   reg$add_group("my_grp", c("CGA009"...))
#   q <- get_Q("C_gar_LG01", 1e6, 5e6)     # instant Q (if configured)
#   stats <- get_region_stats(...)          # unified stats dispatcher
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# Guard: only run once per session
if (exists(".LOAD_BRIDGE_DONE", envir = .GlobalEnv) &&
    isTRUE(get(".LOAD_BRIDGE_DONE", envir = .GlobalEnv))) {
  # Already loaded — skip silently
  invisible(NULL)
} else {

  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a

  # =========================================================================
  # STEP 1: Resolve BASE and key paths
  # =========================================================================

  BASE <- Sys.getenv("BASE", "") %||%
    "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"

  # Try to find our own directory for relative path resolution
  .bridge_self <- tryCatch({
    if (sys.nframe() > 0) normalizePath(sys.frame(1)$ofile, mustWork = FALSE)
    else NULL
  }, error = function(e) NULL)

  .bridge_dir <- if (!is.null(.bridge_self) && file.exists(.bridge_self)) {
    dirname(.bridge_self)
  } else {
    # Fallback: scan likely locations
    cands <- c(
      file.path(BASE, "inversion_codebase_v8.5", "utils"),
      "utils",
      "."
    )
    found <- NULL
    for (d in cands) {
      if (file.exists(file.path(d, "load_bridge.R"))) { found <- d; break }
      if (file.exists(file.path(d, "sample_map.R")))   { found <- d; break }
    }
    found %||% "."
  }

  # =========================================================================
  # STEP 2: Build path registry
  # =========================================================================

  BRIDGE_PATHS <- list(
    BASE              = BASE,
    SAMPLES_IND       = Sys.getenv("SAMPLES_IND", "") %||%
                        file.path(BASE, "het_roh/01_inputs_check/samples.ind"),
    PRUNED_LIST       = Sys.getenv("PRUNED_LIST", "") %||%
                        file.path(BASE, "popstruct_thin/06_relatedness/pruned_samples.txt"),
    SAMPLE_LIST       = Sys.getenv("SAMPLE_LIST", "") %||%
                        file.path(BASE, "popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"),

    # Registry
    REGISTRY_DIR      = Sys.getenv("REGISTRY_DIR", "") %||%
                        file.path(BASE, "sample_registry"),

    # Unified ancestry
    ANCESTRY_CONFIG   = Sys.getenv("ANCESTRY_CONFIG", "") %||%
                        file.path(BASE, "unified_ancestry/00_ancestry_config.sh"),
    INSTANT_Q_R       = Sys.getenv("INSTANT_Q_R", "") %||%
                        file.path(BASE, "unified_ancestry/wrappers/instant_q.R"),
    DISPATCHER_R      = Sys.getenv("DISPATCHER_R", "") %||%
                        file.path(BASE, "unified_ancestry/dispatchers/region_stats_dispatcher.R"),

    # Inversion pipeline
    INVERSION_CONFIG  = Sys.getenv("INVERSION_CONFIG", "") %||%
                        file.path(BASE, "inversion_codebase_v8.5/00_inversion_config.sh"),

    # MODULE_2B
    MODULE2B_CONFIG   = Sys.getenv("MODULE2B_CONFIG", "") %||%
                        file.path(BASE, "MODULE_2B_Structure/00_module2b_config.sh"),

    # BEAGLE
    BEAGLE_DIR        = Sys.getenv("BEAGLE_DIR", "") %||%
                        file.path(BASE, "popstruct_thin/04_beagle_byRF_majmin"),

    # Relatedness
    RELATEDNESS_RES   = Sys.getenv("RELATEDNESS_RES", "") %||%
                        file.path(BASE, "popstruct_thin/06_relatedness",
                                  paste0("catfish_226_relatedness.res")),

    # Theme
    THEME_FILE        = Sys.getenv("THEME_FILE", "") %||%
                        file.path(BASE, "unified_ancestry/utils/theme_systems_plate.R")
  )

  # =========================================================================
  # STEP 3: Source sample_map.R
  # =========================================================================

  .smap_file <- NULL
  for (p in c(
    file.path(.bridge_dir, "sample_map.R"),
    file.path(BASE, "inversion_codebase_v8.5/utils/sample_map.R"),
    "utils/sample_map.R"
  )) {
    if (file.exists(p)) { .smap_file <- p; break }
  }

  if (!is.null(.smap_file)) {
    source(.smap_file, local = FALSE)
    smap <- load_sample_map(BRIDGE_PATHS$SAMPLES_IND)
  } else {
    message("[load_bridge] WARNING: sample_map.R not found")
    smap <- NULL
  }

  # =========================================================================
  # STEP 4: Source sample_registry.R
  # =========================================================================

  .reg_file <- NULL
  for (p in c(
    file.path(.bridge_dir, "sample_registry.R"),
    file.path(BASE, "inversion_codebase_v8.5/utils/sample_registry.R"),
    "utils/sample_registry.R"
  )) {
    if (file.exists(p)) { .reg_file <- p; break }
  }

  if (!is.null(.reg_file)) {
    source(.reg_file, local = FALSE)
    reg <- load_registry(BRIDGE_PATHS$REGISTRY_DIR)
  } else {
    message("[load_bridge] WARNING: sample_registry.R not found")
    reg <- NULL
  }

  # =========================================================================
  # STEP 5: Source unified ancestry (instant_q.R → get_Q, get_Q_summary)
  # =========================================================================

  .iq_loaded <- FALSE
  if (file.exists(BRIDGE_PATHS$INSTANT_Q_R)) {
    tryCatch({
      source(BRIDGE_PATHS$INSTANT_Q_R, local = FALSE)
      if (file.exists(BRIDGE_PATHS$ANCESTRY_CONFIG)) {
        configure_instant_q(config_file = BRIDGE_PATHS$ANCESTRY_CONFIG)
      }
      .iq_loaded <- TRUE
    }, error = function(e) {
      message("[load_bridge] WARNING: instant_q.R load failed: ", conditionMessage(e))
    })
  }

  # =========================================================================
  # STEP 6: Source dispatcher (get_region_stats)
  # =========================================================================

  .disp_loaded <- FALSE
  if (file.exists(BRIDGE_PATHS$DISPATCHER_R)) {
    tryCatch({
      source(BRIDGE_PATHS$DISPATCHER_R, local = FALSE)
      if (file.exists(BRIDGE_PATHS$ANCESTRY_CONFIG)) {
        configure_dispatcher(config_file = BRIDGE_PATHS$ANCESTRY_CONFIG)
      }
      .disp_loaded <- TRUE
    }, error = function(e) {
      message("[load_bridge] WARNING: region_stats_dispatcher.R load failed: ", conditionMessage(e))
    })
  }

  # =========================================================================
  # STEP 7: Convenience functions
  # =========================================================================

  #' Rename all Ind-style columns in a data.table to CGA names.
  #' Handles common prefixes: PC_1_, PC_2_, dosage_, gl_, etc.
  #' @param dt data.table to modify (in-place)
  #' @param prefixes Character vector of column prefixes to scan
  #' @return dt (invisibly)
  rename_ind_to_cga <- function(dt, prefixes = c("PC_1_", "PC_2_")) {
    if (is.null(smap)) {
      message("[load_bridge] smap not loaded; cannot rename")
      return(invisible(dt))
    }
    for (pfx in prefixes) {
      smap$rename_dt_columns(dt, from = "ind", prefix = pfx)
    }
    invisible(dt)
  }

  #' Get real CGA sample names (226 in BAM order)
  get_sample_ids <- function() {
    if (!is.null(smap)) return(smap$real_names)
    # Fallback
    sf <- BRIDGE_PATHS$SAMPLES_IND
    if (file.exists(sf)) {
      x <- trimws(readLines(sf)); x[nzchar(x)]
    } else {
      stop("[load_bridge] Cannot resolve sample names")
    }
  }

  #' Get pruned sample names (81 unrelated)
  get_pruned_ids <- function() {
    pf <- BRIDGE_PATHS$PRUNED_LIST
    if (file.exists(pf)) {
      x <- trimws(readLines(pf)); x[nzchar(x)]
    } else {
      message("[load_bridge] Pruned list not found: ", pf)
      character(0)
    }
  }

  #' Register a sample group (convenience wrapper)
  register_group <- function(group_id, sample_ids, ...) {
    if (is.null(reg)) {
      message("[load_bridge] Registry not loaded; skipping group registration")
      return(invisible(FALSE))
    }
    reg$add_group(group_id, sample_ids, ...)
  }

  # =========================================================================
  # DONE
  # =========================================================================

  assign(".LOAD_BRIDGE_DONE", TRUE, envir = .GlobalEnv)

  message("[load_bridge] Ready: ",
          if (!is.null(smap)) paste0("smap(", smap$n, ") ") else "NO_SMAP ",
          if (!is.null(reg)) paste0("registry(", nrow(reg$list_groups()), " groups) ") else "NO_REG ",
          if (.iq_loaded) "instant_q " else "",
          if (.disp_loaded) "dispatcher " else "")

  # Clean up temporaries
  rm(.smap_file, .reg_file, .iq_loaded, .disp_loaded,
     .bridge_self, .bridge_dir, envir = environment())
}
