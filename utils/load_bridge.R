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
    # Fallback: scan likely locations. Flattened layout first, then legacy v8.5.
    cands <- c(
      file.path(BASE, "utils"),
      file.path(BASE, "inversion_modules", "utils"),
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
    # 2026-04-17: ancestry cache moved out of the code tree.
    LOCAL_Q_DIR       = Sys.getenv("LOCAL_Q_DIR", "") %||%
                        file.path(BASE, "ancestry_cache"),
    CANONICAL_K       = Sys.getenv("CANONICAL_K", "") %||%
                        Sys.getenv("DEFAULT_K", "") %||% "8",
    # Stats cache — chat-15 new store for persisting FST / dxy / pairwise results
    STATS_CACHE_DIR   = Sys.getenv("STATS_CACHE_DIR", "") %||%
                        file.path(BASE, "registries", "data", "stats_cache"),

    # Inversion pipeline — search flattened layout first, legacy v8.5 fallback
    INVERSION_CONFIG  = Sys.getenv("INVERSION_CONFIG", "") %||% {
                          cands <- c(
                            file.path(BASE, "00_inversion_config.sh"),
                            file.path(BASE, "inversion_modules", "00_inversion_config.sh"),
                            file.path(BASE, "inversion_codebase_v8.5", "00_inversion_config.sh")
                          )
                          hit <- cands[file.exists(cands)][1]
                          hit %||% cands[1]
                        },

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
    file.path(BASE, "utils/sample_map.R"),
    file.path(BASE, "inversion_modules/utils/sample_map.R"),
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
    file.path(BASE, "utils/sample_registry.R"),
    file.path(BASE, "inversion_modules/utils/sample_registry.R"),
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
  # STEP 6.5: Auto-register ancestry-derived groups (chat-15)
  # =========================================================================
  # The inversion sample_registry is the source of truth for "who is in which
  # group." Historically the ancestry module referenced group names like
  # `ancestry_K8_Q3` or `all_226` without registering them. This step closes
  # that gap by registering:
  #
  #   all_226        — every sample in sample_master (full cohort)
  #   unrelated_81   — NAToRA-pruned subset from PRUNED_LIST
  #   ancestry_K<K>_Q1 .. Q<K>  — samples whose canonical-K assigned_pop == k
  #                               (read from the canonical-K local_Q cache's
  #                               per-sample summary; fallback to an empty
  #                               group if the cache hasn't been populated yet)
  #
  # Idempotent: sample_registry's add_group() skips if the group already
  # exists (overwrite=FALSE by default). Safe to source load_bridge.R
  # multiple times.
  .anc_groups_registered <- 0L
  if (!is.null(reg)) {
    tryCatch({
      # all_226
      if (!is.null(smap) && length(smap$real_names) > 0 &&
          !reg$has_group("all_226")) {
        reg$add_group("all_226", smap$real_names,
                       description = "Full cohort (226 samples)")
        .anc_groups_registered <- .anc_groups_registered + 1L
      }

      # unrelated_81
      pf <- BRIDGE_PATHS$PRUNED_LIST
      if (nzchar(pf) && file.exists(pf) && !reg$has_group("unrelated_81")) {
        pruned_ids <- trimws(readLines(pf, warn = FALSE))
        pruned_ids <- pruned_ids[nzchar(pruned_ids)]
        if (length(pruned_ids) > 0) {
          reg$add_group("unrelated_81", pruned_ids,
                         description = "NAToRA-pruned unrelated subset")
          .anc_groups_registered <- .anc_groups_registered + 1L
        }
      }

      # ancestry_K<K>_Q<k> — read from canonical-K samples cache, any chrom.
      # The assigned_pop column partitions samples into K clusters; we collapse
      # over chromosomes by majority vote (most-frequent assigned_pop per
      # sample). This yields a stable cohort-level ancestry-cluster membership
      # even when per-window assignments wobble.
      canon_K <- as.integer(BRIDGE_PATHS$CANONICAL_K %||% 8L)
      lqd <- BRIDGE_PATHS$LOCAL_Q_DIR
      k_shard <- file.path(lqd, sprintf("K%02d", canon_K))
      # Look for a samples.tsv.gz in the K-sharded dir; fall back to legacy
      # flat layout at the canonical K.
      samples_files <- character(0)
      if (dir.exists(k_shard)) {
        samples_files <- list.files(k_shard,
                                     pattern = "\\.local_Q_samples\\.tsv\\.gz$",
                                     full.names = TRUE)
      }
      if (length(samples_files) == 0 && dir.exists(lqd)) {
        # Legacy flat — treat as canonical K
        samples_files <- list.files(lqd,
                                     pattern = "\\.local_Q_samples\\.tsv\\.gz$",
                                     full.names = TRUE)
      }
      if (length(samples_files) > 0) {
        # Read up to a few files to stabilise assigned_pop; one file is usually
        # enough but reading 2-3 makes the majority vote less brittle to any
        # single chrom being weird.
        sf_use <- head(samples_files, 3L)
        parts <- lapply(sf_use, function(f) {
          tryCatch(data.table::fread(f, select = c("sample_id", "assigned_pop")),
                    error = function(e) NULL)
        })
        parts <- Filter(Negate(is.null), parts)
        if (length(parts) > 0) {
          all_rows <- data.table::rbindlist(parts, fill = TRUE)
          if ("sample_id" %in% names(all_rows) &&
              "assigned_pop" %in% names(all_rows) &&
              nrow(all_rows) > 0) {
            mode_pop <- all_rows[,
              list(pop = names(sort(table(assigned_pop), decreasing = TRUE))[1]),
              by = sample_id]
            for (k_val in seq_len(canon_K)) {
              gid <- sprintf("ancestry_K%d_Q%d", canon_K, k_val)
              if (!reg$has_group(gid)) {
                members <- mode_pop[as.integer(pop) == k_val, sample_id]
                reg$add_group(gid, as.character(members),
                               description = sprintf(
                                 "Canonical-K=%d ancestry cluster %d (majority-vote assigned_pop)",
                                 canon_K, k_val))
                .anc_groups_registered <- .anc_groups_registered + 1L
              }
            }
          }
        }
      }
    }, error = function(e) {
      message("[load_bridge] ancestry group registration failed: ",
              conditionMessage(e), " — continuing")
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
          if (.disp_loaded) "dispatcher " else "",
          if (.anc_groups_registered > 0L)
            paste0("ancestry_groups_new=", .anc_groups_registered, " ") else "")

  # Clean up temporaries
  rm(.smap_file, .reg_file, .iq_loaded, .disp_loaded, .anc_groups_registered,
     .bridge_self, .bridge_dir, envir = environment())
}
