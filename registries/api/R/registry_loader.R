#!/usr/bin/env Rscript
# =============================================================================
# registry_loader.R — Unified entry point for the three registries
# =============================================================================
# Provides:    reg$samples    (WHO  — sample groups, master)
#              reg$intervals  (WHERE — windows, candidate intervals, cov)
#              reg$evidence   (WHAT — keys, structured blocks, per-candidate)
#              reg$query      (composite multi-registry traversals)
#              reg$compute    (query + call into Engine B for live stats)
#
# This is the Stage-2 library referenced in the v10 handoff. It re-uses
# the existing utils/sample_registry.R (renamed to api/R/sample_registry.R)
# and extends it with interval + evidence + query layers.
#
# File layout expected:
#   registries/
#     data/
#       sample_registry/
#         sample_master.tsv
#         sample_groups.tsv
#         groups/<group_id>.txt
#       interval_registry/
#         windows_master.tsv.gz
#         candidate_intervals.tsv
#         cov_registry.tsv
#       evidence_registry/
#         per_candidate/<cid>/
#           raw/         Tier 1
#           structured/  Tier 2 (JSON blocks)
#           keys.tsv     Tier 3 (flat)
#           characterization.json
#           classification.json
#         global/
#           cross_candidate_ranking.tsv
#           fdr_corrections.tsv
#     schemas/structured_block_schemas/*.schema.json
#
# Usage in a pipeline script:
#   source("registries/api/R/registry_loader.R")
#   reg <- load_registry()   # or load_registry(registries_root = "...")
#   reg$evidence$write_block("LG12_17", "boundary_left", list(...))
#   samps <- reg$samples$get_group("inv_LG12_17_HOM_REF")
#   cand  <- reg$intervals$get_candidate("LG12_17")
#
# FALLBACK: if the library can't find its registries root, it returns a
# shim object whose methods all emit warnings and write to a fallback
# directory. This means existing scripts sourced with this library keep
# running even before the v10 migration completes.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

# =============================================================================
# Top-level loader
# =============================================================================
load_registry <- function(registries_root = NULL, create_if_missing = TRUE) {

  # ── Resolve registries root ─────────────────────────────────────────────
  if (is.null(registries_root)) {
    registries_root <- Sys.getenv("REGISTRIES", "")
    if (!nzchar(registries_root)) {
      base <- Sys.getenv("BASE", ".")
      # v10 convention: ${BASE}/inversion-popgen-toolkit/registries
      cand <- file.path(base, "inversion-popgen-toolkit", "registries")
      if (dir.exists(cand)) {
        registries_root <- cand
      } else {
        # v9 fallback: use sample_registry at top of BASE
        registries_root <- file.path(base, "sample_registry_v10_fallback")
      }
    }
  }

  data_dir   <- file.path(registries_root, "data")
  schema_dir <- file.path(registries_root, "schemas")

  sample_dir   <- file.path(data_dir, "sample_registry")
  interval_dir <- file.path(data_dir, "interval_registry")
  evid_dir     <- file.path(data_dir, "evidence_registry")

  if (create_if_missing) {
    for (d in c(sample_dir, interval_dir,
                file.path(evid_dir, "per_candidate"),
                file.path(evid_dir, "global"))) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }

  message("[registry_loader] root = ", registries_root)

  # ── Samples sublibrary ──────────────────────────────────────────────────
  # Prefer existing sample_registry.R implementation if present; else shim.
  samples <- load_samples_api(sample_dir)

  # ── Intervals sublibrary ────────────────────────────────────────────────
  intervals <- load_intervals_api(interval_dir)

  # ── Evidence sublibrary ─────────────────────────────────────────────────
  evidence <- load_evidence_api(evid_dir, schema_dir)

  # ── Query sublibrary (composite) ────────────────────────────────────────
  query <- load_query_api(samples, intervals, evidence)

  # ── Compute sublibrary (live, calls Engine B) ───────────────────────────
  compute <- load_compute_api(samples, intervals, evidence)

  # ── Results sublibrary (chat-16 rewrite of chat-15 stats_cache) ─────────
  # See registries/DATABASE_DESIGN.md for the full design.
  # Dir resolution: prefer explicit env (RESULTS_REGISTRY_DIR), then legacy
  # STATS_CACHE_DIR (back-compat), then the canonical location under data_dir.
  results_dir <- Sys.getenv("RESULTS_REGISTRY_DIR", "") %||%
                  Sys.getenv("STATS_CACHE_DIR", "") %||%
                  file.path(data_dir, "results_registry")

  # One-time migration: if the chat-15 stats_cache dir exists and
  # results_registry does not, rename. Idempotent and safe — preserves files.
  legacy_stats_dir <- file.path(data_dir, "stats_cache")
  if (dir.exists(legacy_stats_dir) && !dir.exists(results_dir)) {
    tryCatch({
      file.rename(legacy_stats_dir, results_dir)
      message("[registry_loader] migrated chat-15 stats_cache/ → results_registry/")
    }, error = function(e) {
      message("[registry_loader] note: could not migrate stats_cache/: ",
              conditionMessage(e))
    })
  }

  results <- load_results_api(results_dir,
                               samples = samples, intervals = intervals)

  # Publish for compute$ancestry_q_and_f_for_candidate to find via .REG_STATS_API
  # (deprecated name kept so chat-15 call sites keep working)
  assign(".REG_STATS_API", results, envir = .GlobalEnv)
  assign(".REG_RESULTS_API", results, envir = .GlobalEnv)

  # Return composite
  list(
    root      = registries_root,
    samples   = samples,
    intervals = intervals,
    evidence  = evidence,
    query     = query,
    compute   = compute,
    results   = results,    # NEW: preferred name (chat 16)
    stats     = results,    # DEPRECATED alias for one chat cycle

    # Shortcuts for the most common v9-style calls (forward-compat during migration)
    # Chat-18: extended to cover every flat-API method that real consumer scripts
    # (load_bridge.R, region_stats_dispatcher.R, STEP_C01f, cheat6, etc.) still
    # use. Adding an alias here is cheaper than rewriting each call site. New
    # code should still use the nested form (reg$samples$*).
    add_group         = samples$add_group,
    get_group         = samples$get_group,
    has_group         = samples$has_group,
    list_groups       = samples$list_groups,
    get_master        = samples$get_master,
    count_groups      = samples$count_groups,
    get_groups_for_candidate = samples$get_groups_for_candidate,
    register_candidate = evidence$register_candidate,
    add_evidence      = evidence$add_evidence,
    get_evidence      = evidence$get_evidence,
    has_evidence      = evidence$has_evidence,

    status = function() {
      cat("═══ Registry Status ═══\n")
      cat("root:", registries_root, "\n")
      cat("samples:   ", samples$count_groups(), " groups\n", sep = "")
      cat("intervals: ", intervals$count_candidates(), " candidates\n", sep = "")
      cat("evidence:  ", evidence$count_candidates(), " candidates with evidence\n", sep = "")
      n_results <- tryCatch(nrow(results$list_cached()),
                            error = function(e) 0L)
      cat("results:   ", n_results, " manifest rows\n", sep = "")
    }
  )
}

# =============================================================================
# Samples API — thin wrapper around existing sample_registry.R + extensions
# =============================================================================
# Chat 11 additions:
#   - Fixed pre-existing bug: get_master was aliased to full$list_groups (wrong
#     table). Now wired to full$get_master. [Finding AJ — reported in audit log]
#   - Reverse lookup: get_sample_groups, list_carriers_across_candidates,
#     list_recombinant_candidates, get_groups_for_candidate
#   - Pedigree: get_family, list_families, get_family_members
#   - Group set ops: compute_group_overlap, find_co_segregating_groups
#   - Metadata: get_sample_metadata
#
# Design: all new methods read fresh from disk per call (same convention as
# intervals API). Whole-genome this is a few MB; storage is negligible per
# Quentin. Group member files live at <sample_dir>/groups/<gid>.txt.
# =============================================================================
load_samples_api <- function(sample_dir) {

  # ── Build extended method closure on top of a `full` API ────────────────
  # Factored so both the "real" path (sample_registry.R found) and the
  # shim path can share it. `full` must expose: get_master, list_groups,
  # get_group, has_group, add_group.
  build_extended <- function(full) {

    # Internal: path to a group's members.txt (tolerates both absolute
    # members_file paths from the groups table and the default
    # <sample_dir>/groups/<gid>.txt convention).
    members_path_for <- function(gid) {
      gf <- tryCatch(full$list_groups(), error = function(e) data.table())
      if (nrow(gf) > 0 && gid %in% gf$group_id) {
        mf <- gf$members_file[gf$group_id == gid][1]
        if (!is.na(mf) && nzchar(mf) && file.exists(mf)) return(mf)
      }
      # Fallback convention
      mf2 <- file.path(sample_dir, "groups", paste0(gid, ".txt"))
      if (file.exists(mf2)) mf2 else NA_character_
    }

    read_members <- function(gid) {
      mf <- members_path_for(gid)
      if (is.na(mf)) return(character(0))
      s <- trimws(readLines(mf, warn = FALSE))
      s[nzchar(s)]
    }

    # One-time warning cache for family column absence
    fam_warned <- FALSE
    check_family_col <- function(m) {
      has <- "family_id" %in% names(m)
      if (!has && !fam_warned) {
        warning("[registry] sample_master.tsv has no 'family_id' column — ",
                "family methods return NA/empty. ",
                "Populate it via reg$samples$add_master_column('family_id', ...) ",
                "or re-run sample_master build with pedigree input.")
        fam_warned <<- TRUE
      }
      has
    }

    list(
      # ── Originals (preserved) ───────────────────────────────────────────
      get_master   = full$get_master,
      get_group    = full$get_group,
      add_group    = full$add_group,
      has_group    = full$has_group,
      list_groups  = full$list_groups,
      count_groups = function() nrow(full$list_groups()),
      get_carriers = function(cid, status = "HOM_INV") {
        gid <- paste0("inv_", cid, "_", status)
        if (full$has_group(gid)) full$get_group(gid) else character(0)
      },

      # ══ chat 11 extensions ═══════════════════════════════════════════

      # ── get_sample_metadata: sample_master row as list ────────────────
      get_sample_metadata = function(sid) {
        m <- tryCatch(full$get_master(), error = function(e) data.table())
        if (!is.data.table(m)) m <- as.data.table(m)
        if (nrow(m) == 0 || !"sample_id" %in% names(m)) return(NULL)
        r <- m[m$sample_id == sid]
        if (nrow(r) == 0) NULL else as.list(r[1])
      },

      # ── get_sample_groups: all groups containing sid ──────────────────
      # Returns character vector of group_ids. Reads every group file,
      # which is fine at ~hundreds of groups.
      get_sample_groups = function(sid) {
        gf <- tryCatch(full$list_groups(), error = function(e) data.table())
        if (nrow(gf) == 0) return(character(0))
        hits <- character(0)
        for (gid in gf$group_id) {
          if (sid %in% read_members(gid)) hits <- c(hits, gid)
        }
        hits
      },

      # ── get_family: family_id for a sample (NA if unknown/no column) ──
      get_family = function(sid) {
        m <- tryCatch(full$get_master(), error = function(e) data.table())
        if (!is.data.table(m)) m <- as.data.table(m)
        if (!check_family_col(m)) return(NA_character_)
        r <- m[m$sample_id == sid]
        if (nrow(r) == 0) return(NA_character_)
        f <- r$family_id[1]
        if (is.na(f) || !nzchar(as.character(f))) NA_character_ else as.character(f)
      },

      # ── list_families: distinct non-NA family IDs ─────────────────────
      list_families = function() {
        m <- tryCatch(full$get_master(), error = function(e) data.table())
        if (!is.data.table(m)) m <- as.data.table(m)
        if (!check_family_col(m)) return(character(0))
        fams <- unique(as.character(m$family_id))
        fams <- fams[!is.na(fams) & nzchar(fams)]
        sort(fams)
      },

      # ── get_family_members: samples in a family ───────────────────────
      get_family_members = function(fam_id) {
        m <- tryCatch(full$get_master(), error = function(e) data.table())
        if (!is.data.table(m)) m <- as.data.table(m)
        if (!check_family_col(m)) return(character(0))
        m$sample_id[!is.na(m$family_id) & as.character(m$family_id) == fam_id]
      },

      # ── list_carriers_across_candidates ───────────────────────────────
      # Returns unique sample IDs that are `status` for any candidate.
      # Scans group table for group_ids ending in "_<status>", unions members.
      list_carriers_across_candidates = function(status = "HOM_INV") {
        gf <- tryCatch(full$list_groups(), error = function(e) data.table())
        if (nrow(gf) == 0) return(character(0))
        suffix <- paste0("_", status)
        # Use the subgroup column if present (more reliable than name suffix)
        if ("subgroup" %in% names(gf)) {
          gids <- gf$group_id[gf$subgroup == status]
        } else {
          gids <- gf$group_id[endsWith(gf$group_id, suffix)]
        }
        if (length(gids) == 0) return(character(0))
        all_s <- unique(unlist(lapply(gids, read_members), use.names = FALSE))
        sort(all_s)
      },

      # ── compute_group_overlap(gid1, gid2) ─────────────────────────────
      # Returns list with intersection/union counts, jaccard, set-diffs.
      compute_group_overlap = function(gid1, gid2) {
        s1 <- read_members(gid1)
        s2 <- read_members(gid2)
        inter <- intersect(s1, s2)
        uni   <- union(s1, s2)
        list(
          gid1         = gid1,
          gid2         = gid2,
          n1           = length(s1),
          n2           = length(s2),
          intersection = length(inter),
          union        = length(uni),
          jaccard      = if (length(uni) == 0) NA_real_ else length(inter) / length(uni),
          only_in_1    = setdiff(s1, s2),
          only_in_2    = setdiff(s2, s1),
          shared       = inter
        )
      },

      # ── list_recombinant_candidates(sid) ──────────────────────────────
      # Returns CIDs where sid is in RECOMBINANT (or sub-classes).
      # Uses the `inv_<cid>_RECOMBINANT...` naming convention from
      # STEP_C01i_d_seal.R register_all_groups(). Set include_subclasses
      # = TRUE (default) to also count RECOMBINANT_GC / RECOMBINANT_DCO.
      list_recombinant_candidates = function(sid, include_subclasses = TRUE) {
        grps <- Recall_get_sample_groups <- NULL  # silence linter
        grps <- character(0)
        gf <- tryCatch(full$list_groups(), error = function(e) data.table())
        if (nrow(gf) == 0) return(character(0))
        for (gid in gf$group_id) {
          if (sid %in% read_members(gid)) grps <- c(grps, gid)
        }
        # Filter to inv_<cid>_RECOMBINANT[_SUBCLASS]?
        pat <- if (include_subclasses) {
          "^inv_(.+)_RECOMBINANT(_GC|_DCO)?$"
        } else {
          "^inv_(.+)_RECOMBINANT$"
        }
        hits <- grps[grepl(pat, grps)]
        cids <- sub(pat, "\\1", hits)
        unique(cids)
      },

      # ── get_groups_for_candidate(cid) ─────────────────────────────────
      # Returns named list of all groups for cid in one call. Saves
      # C01f from four separate has_group/get_group calls (currently
      # STEP_C01f_hypothesis_tests.R L307-325 does exactly that).
      # Each slot is character(0) if that group is not registered.
      get_groups_for_candidate = function(cid) {
        statuses <- c("HOM_REF", "HET", "HOM_INV", "RECOMBINANT",
                      "RECOMBINANT_GC", "RECOMBINANT_DCO", "HOM_STD")
        out <- vector("list", length(statuses))
        names(out) <- statuses
        for (st in statuses) {
          gid <- paste0("inv_", cid, "_", st)
          out[[st]] <- if (full$has_group(gid)) full$get_group(gid) else character(0)
        }
        out
      },

      # ── find_co_segregating_groups(min_jaccard) ───────────────────────
      # Scans all pairs of groups, returns data.table of (gid1, gid2,
      # n1, n2, intersection, jaccard) where jaccard >= min_jaccard.
      # Restricts to inversion groups (prefix "inv_") by default to avoid
      # N² blowup over non-comparable groups like "unrelated_81".
      # For very-high-count pairs use compute_group_overlap on specific
      # pairs instead.
      find_co_segregating_groups = function(min_jaccard = 0.9,
                                               group_prefix = "inv_",
                                               min_size = 3L) {
        gf <- tryCatch(full$list_groups(), error = function(e) data.table())
        if (nrow(gf) == 0) return(data.table())
        gids <- gf$group_id
        if (!is.null(group_prefix) && nzchar(group_prefix)) {
          gids <- gids[startsWith(gids, group_prefix)]
        }
        if (length(gids) < 2) return(data.table())

        # Preload members to avoid O(N²) file reads
        mem_cache <- setNames(
          lapply(gids, read_members),
          gids
        )
        keep <- vapply(mem_cache, function(x) length(x) >= min_size, logical(1))
        gids <- gids[keep]
        mem_cache <- mem_cache[keep]
        if (length(gids) < 2) return(data.table())

        out <- list()
        for (i in seq_len(length(gids) - 1L)) {
          s1 <- mem_cache[[i]]
          for (j in (i + 1L):length(gids)) {
            s2 <- mem_cache[[j]]
            uni_n <- length(union(s1, s2))
            if (uni_n == 0) next
            inter_n <- length(intersect(s1, s2))
            jac <- inter_n / uni_n
            if (jac >= min_jaccard) {
              out[[length(out) + 1L]] <- data.table(
                gid1 = gids[i], gid2 = gids[j],
                n1 = length(s1), n2 = length(s2),
                intersection = inter_n, jaccard = round(jac, 4)
              )
            }
          }
        }
        if (length(out) == 0) return(data.table(
          gid1 = character(), gid2 = character(),
          n1 = integer(), n2 = integer(),
          intersection = integer(), jaccard = numeric()
        ))
        rbindlist(out)
      }
    )
  }

  # Try to source the full sample_registry.R; if unavailable, emit shim.
  src_candidates <- c(
    file.path(sample_dir, "..", "..", "api", "R", "sample_registry.R"),
    "utils/sample_registry.R",
    "../utils/sample_registry.R",
    file.path(Sys.getenv("BASE", ""),
              "inversion-popgen-toolkit/registries/api/R/sample_registry.R"),
    file.path(Sys.getenv("BASE", ""),
              "inversion_codebase_v8.5/utils/sample_registry.R")
  )
  for (sp in src_candidates) {
    if (file.exists(sp)) {
      source(sp)
      if (exists("load_registry", mode = "function")) {
        full <- load_registry(sample_dir)   # the OLD sample_registry.R entry point
        return(build_extended(full))
      }
    }
  }

  # Shim: no-op implementations with warnings
  warning("[registry_loader] sample_registry.R not found — using shim")
  shim_full <- list(
    get_master   = function() data.table(),
    get_group    = function(gid) character(0),
    add_group    = function(...) invisible(FALSE),
    has_group    = function(gid) FALSE,
    list_groups  = function() data.table()
  )
  build_extended(shim_full)
}

# =============================================================================
# Intervals API — windows_master, candidate_intervals, cov_registry
# =============================================================================
load_intervals_api <- function(interval_dir) {
  win_file  <- file.path(interval_dir, "windows_master.tsv.gz")
  cand_file <- file.path(interval_dir, "candidate_intervals.tsv")
  cov_file  <- file.path(interval_dir, "cov_registry.tsv")

  # Ensure files exist (empty skeletons if not)
  if (!file.exists(cand_file)) {
    fwrite(data.table(candidate_id = character(), chrom = character(),
                       start_bp = integer(), end_bp = integer(),
                       size_kb = numeric(), scale = character(),
                       parent_id = character()),
           cand_file, sep = "\t")
  }
  if (!file.exists(cov_file)) {
    fwrite(data.table(cov_id = character(), scope = character(),
                       chrom = character(), cov_path = character()),
           cov_file, sep = "\t")
  }

  # ── Internal helpers shared by new methods ──────────────────────────────
  # Read the candidate table once per call. Whole-genome this is a few MB;
  # storage + parse cost is negligible per Quentin (chat-11 handoff).
  read_cands <- function() {
    if (!file.exists(cand_file)) return(data.table(
      candidate_id = character(), chrom = character(),
      start_bp = integer(), end_bp = integer(),
      size_kb = numeric(), scale = character(),
      parent_id = character()
    ))
    fread(cand_file)
  }

  list(
    # ── Original 7 methods (chat 5–9, unchanged) ──────────────────────────
    get_windows = function(chrom = NULL, start_bp = NULL, end_bp = NULL) {
      if (!file.exists(win_file)) return(data.table())
      w <- fread(win_file)
      if (!is.null(chrom))   w <- w[w$chrom == chrom]
      if (!is.null(start_bp)) w <- w[w$end_bp   >= start_bp]
      if (!is.null(end_bp))   w <- w[w$start_bp <= end_bp]
      w
    },
    get_window_at = function(chrom, pos) {
      w <- fread(win_file)
      w <- w[w$chrom == chrom & w$start_bp <= pos & w$end_bp >= pos]
      if (nrow(w) == 0) NULL else w[1]
    },
    get_candidate = function(cid) {
      c <- read_cands()
      r <- c[c$candidate_id == cid]
      if (nrow(r) == 0) NULL else as.list(r[1])
    },
    get_candidate_boundaries = function(cid) {
      r <- read_cands()
      r <- r[r$candidate_id == cid]
      if (nrow(r) == 0) return(NULL)
      list(left_bp = as.integer(r$start_bp[1]), right_bp = as.integer(r$end_bp[1]))
    },
    add_candidate = function(cid, chrom, start_bp, end_bp, scale = "inversion_block",
                               parent_id = NA_character_) {
      c <- read_cands()
      if (cid %in% c$candidate_id) return(invisible(FALSE))
      size_kb <- round((end_bp - start_bp) / 1000, 1)
      c <- rbind(c, data.table(candidate_id = cid, chrom = chrom,
                                start_bp = as.integer(start_bp),
                                end_bp = as.integer(end_bp),
                                size_kb = size_kb, scale = scale,
                                parent_id = parent_id), fill = TRUE)
      fwrite(c, cand_file, sep = "\t")
      invisible(TRUE)
    },
    count_candidates = function() {
      if (!file.exists(cand_file)) return(0L)
      nrow(fread(cand_file))
    },
    get_cov = function(scope, chrom = NULL) {
      c <- fread(cov_file)
      c <- c[c$scope == scope]
      if (!is.null(chrom)) c <- c[c$chrom == chrom]
      c
    },

    # ══ chat 11 extensions ═══════════════════════════════════════════════
    # Nested / overlapping / relationship-classification methods. These
    # support C01d pass-2, C01b seeded-regions (parent_id population),
    # phase 4e characterization, and phase 4c overlap reconciliation.

    # ── get_children: direct children (one level) ──────────────────────
    # Returns character(0) if no children. parent_id NA/empty/"." → not a child.
    get_children = function(cid) {
      c <- read_cands()
      if (nrow(c) == 0 || !"parent_id" %in% names(c)) return(character(0))
      mask <- !is.na(c$parent_id) & nzchar(c$parent_id) &
              c$parent_id != "." & c$parent_id == cid
      c$candidate_id[mask]
    },

    # ── get_parent: parent CID or NA_character_ ────────────────────────
    get_parent = function(cid) {
      c <- read_cands()
      r <- c[c$candidate_id == cid]
      if (nrow(r) == 0 || !"parent_id" %in% names(c)) return(NA_character_)
      p <- r$parent_id[1]
      if (is.na(p) || !nzchar(p) || p == ".") NA_character_ else p
    },

    # ── get_ancestors: walk up parent_id chain (closest first) ─────────
    # Cycle-safe via `seen` set. Returns character(0) if no parent.
    get_ancestors = function(cid) {
      c <- read_cands()
      if (nrow(c) == 0 || !"parent_id" %in% names(c)) return(character(0))
      out <- character(0)
      seen <- cid
      cur  <- cid
      repeat {
        r <- c[c$candidate_id == cur]
        if (nrow(r) == 0) break
        p <- r$parent_id[1]
        if (is.na(p) || !nzchar(p) || p == "." || p %in% seen) break
        out  <- c(out, p)
        seen <- c(seen, p)
        cur  <- p
      }
      out
    },

    # ── get_descendants: recursive down parent_id relation ─────────────
    # depth: Inf (default) walks full subtree; 1 = direct children only.
    # Returns character(0) if no descendants. Returns in BFS order.
    get_descendants = function(cid, depth = Inf) {
      c <- read_cands()
      if (nrow(c) == 0 || !"parent_id" %in% names(c)) return(character(0))
      out <- character(0)
      frontier <- cid
      level <- 0L
      while (length(frontier) > 0 && level < depth) {
        mask <- !is.na(c$parent_id) & nzchar(c$parent_id) &
                c$parent_id != "." & c$parent_id %in% frontier
        next_f <- c$candidate_id[mask]
        next_f <- setdiff(next_f, out)        # cycle guard
        next_f <- setdiff(next_f, cid)        # root guard
        if (length(next_f) == 0) break
        out <- c(out, next_f)
        frontier <- next_f
        level <- level + 1L
      }
      out
    },

    # ── get_overlapping: CIDs whose intervals overlap [start_bp, end_bp] ─
    # On the same chromosome. exclude_cid removes a specific CID
    # (typically the caller's own). Uses half-open-style overlap test:
    # overlap iff (A.start <= B.end) && (A.end >= B.start).
    get_overlapping = function(chrom, start_bp, end_bp, exclude_cid = NULL) {
      c <- read_cands()
      if (nrow(c) == 0) return(character(0))
      mask <- c$chrom == chrom &
              c$start_bp <= end_bp &
              c$end_bp   >= start_bp
      if (!is.null(exclude_cid)) mask <- mask & c$candidate_id != exclude_cid
      c$candidate_id[mask]
    },

    # ── get_nested_within: CIDs strictly contained in cid's interval ───
    # Coordinate-based (doesn't require parent_id population). Strict
    # containment: child.start > parent.start AND child.end < parent.end,
    # same chromosome, excluding cid itself. Useful for finding nested
    # relationships BEFORE parent_id is assigned (chat-14 pass-2).
    get_nested_within = function(cid) {
      c <- read_cands()
      r <- c[c$candidate_id == cid]
      if (nrow(r) == 0) return(character(0))
      pstart <- r$start_bp[1]; pend <- r$end_bp[1]; pchr <- r$chrom[1]
      mask <- c$chrom == pchr &
              c$start_bp >  pstart &
              c$end_bp   <  pend &
              c$candidate_id != cid
      c$candidate_id[mask]
    },

    # ── get_children(cid) — chat 16 ────────────────────────────────────
    # (Already defined above at the intervals API head. Not redefined.)

    # ── get_tree(root_cid) — chat 16 ──────────────────────────────────
    # Hierarchical descendant view via parent_id walk. Returns a list
    # with:
    #   $root       — the root_cid
    #   $all_cids   — root + every descendant (BFS order)
    #   $by_depth   — list of cid vectors, one per depth level
    #   $edges      — data.table(parent, child) for the subtree
    #
    # This is how plot-assembly and classification-loop scripts iterate
    # a candidate system when biology has subdivided it. Keeps
    # evidence_registry/per_candidate/ filesystem flat while providing
    # a hierarchical view on demand. See DATABASE_DESIGN.md § "One
    # candidate, one folder: flat filesystem, tree on demand".
    get_tree = function(root_cid, max_depth = 20L) {
      c <- read_cands()
      if (!("parent_id" %in% names(c))) {
        return(list(root = root_cid, all_cids = root_cid,
                     by_depth = list(root_cid),
                     edges = data.table(parent = character(0),
                                         child  = character(0))))
      }
      all_cids <- root_cid
      by_depth <- list(root_cid)
      edges    <- data.table(parent = character(0), child = character(0))
      frontier <- root_cid
      for (depth in seq_len(max_depth)) {
        children <- c$candidate_id[!is.na(c$parent_id) &
                                     c$parent_id %in% frontier]
        if (!length(children)) break
        # Build edges from frontier to children
        new_edges <- c[!is.na(c$parent_id) & c$parent_id %in% frontier,
                       .(parent = parent_id, child = candidate_id)]
        edges    <- rbind(edges, new_edges)
        all_cids <- c(all_cids, children)
        by_depth[[depth + 1L]] <- children
        frontier <- children
      }
      list(root = root_cid,
           all_cids = unique(all_cids),
           by_depth = by_depth,
           edges    = edges)
    },

    # ── resolve_smallest_at(chrom, pos) — chat 16 ─────────────────────
    # Returns the candidate_id of the SMALLEST interval containing pos,
    # or NA_character_ if no candidate contains pos. This is the primitive
    # used by multi-scale scans: for any window position along a scan,
    # this tells you which candidate (parent, segment, nested) you should
    # use to look up the local karyotype groups.
    #
    # For a window-by-window FST scan across an inversion with recombinant
    # segments, call this per window to pick the right candidate, then
    # build the group names as sprintf("inv_%s_HOM_REF", cid) etc.
    # Falls back to parent scale automatically when the window is in a
    # zone with no recombinant segmentation.
    resolve_smallest_at = function(chrom, pos) {
      c <- read_cands()
      if (nrow(c) == 0) return(NA_character_)
      pos <- as.integer(pos)
      mask <- c$chrom == chrom & c$start_bp <= pos & c$end_bp >= pos
      if (!any(mask)) return(NA_character_)
      hits <- c[mask]
      # Pick smallest by span
      hits$.span <- hits$end_bp - hits$start_bp
      hits$candidate_id[which.min(hits$.span)]
    },

    # ── classify_relationship: pairwise interval relationship ──────────
    # Returns one of:
    #   "equal"           — identical start+end on same chrom
    #   "nested_1_in_2"   — cid1 strictly inside cid2
    #   "nested_2_in_1"   — cid2 strictly inside cid1
    #   "partial_overlap" — overlap but neither contains the other
    #   "disjoint"        — no overlap (incl. different chroms)
    #   NA_character_     — either CID not found
    # "Strict" nesting uses the same rule as get_nested_within.
    classify_relationship = function(cid1, cid2) {
      c <- read_cands()
      r1 <- c[c$candidate_id == cid1]; r2 <- c[c$candidate_id == cid2]
      if (nrow(r1) == 0 || nrow(r2) == 0) return(NA_character_)
      if (r1$chrom[1] != r2$chrom[1]) return("disjoint")
      s1 <- r1$start_bp[1]; e1 <- r1$end_bp[1]
      s2 <- r2$start_bp[1]; e2 <- r2$end_bp[1]
      if (s1 == s2 && e1 == e2) return("equal")
      # No overlap at all?
      if (e1 < s2 || e2 < s1) return("disjoint")
      # Nesting (strict): cid1 fully inside cid2 with at least one strict bound
      if (s1 >= s2 && e1 <= e2 && (s1 > s2 || e1 < e2)) return("nested_1_in_2")
      if (s2 >= s1 && e2 <= e1 && (s2 > s1 || e2 < e1)) return("nested_2_in_1")
      "partial_overlap"
    },

    # ── update_candidate: partial update of a single candidate row ─────
    # Useful for C01d pass-2 setting parent_id post-hoc, or for any late
    # metadata update. `...` is any named arg matching a column in the
    # candidate schema; unknown names are silently skipped with a warning.
    # Returns TRUE on success, FALSE if CID not found.
    update_candidate = function(cid, ...) {
      updates <- list(...)
      if (length(updates) == 0) return(invisible(FALSE))
      c <- read_cands()
      if (nrow(c) == 0 || !(cid %in% c$candidate_id)) {
        warning("[registry] update_candidate: unknown CID ", cid)
        return(invisible(FALSE))
      }
      idx <- which(c$candidate_id == cid)[1]
      unknown <- setdiff(names(updates), names(c))
      if (length(unknown) > 0) {
        warning("[registry] update_candidate: unknown columns skipped: ",
                paste(unknown, collapse = ", "))
      }
      for (col in intersect(names(updates), names(c))) {
        v <- updates[[col]]
        # Coerce to column type to avoid data.table type coercion surprises
        set(c, i = idx, j = col, value = v)
      }
      # Recompute size_kb if start/end changed
      if (any(c("start_bp", "end_bp") %in% names(updates))) {
        set(c, i = idx, j = "size_kb",
            value = round((c$end_bp[idx] - c$start_bp[idx]) / 1000, 1))
      }
      fwrite(c, cand_file, sep = "\t")
      invisible(TRUE)
    },

    # ── bulk_add_candidates: batch insert from data.table ──────────────
    # Much faster than row-by-row add_candidate for large seeded-region
    # catalogs (chat 13). Required columns: candidate_id, chrom, start_bp,
    # end_bp. Optional: scale, parent_id. Duplicates (by candidate_id)
    # are silently skipped. Returns number of rows actually inserted.
    bulk_add_candidates = function(dt) {
      stopifnot(is.data.frame(dt))
      dt <- as.data.table(dt)
      req <- c("candidate_id", "chrom", "start_bp", "end_bp")
      miss <- setdiff(req, names(dt))
      if (length(miss) > 0) {
        stop("[registry] bulk_add_candidates: missing required columns: ",
             paste(miss, collapse = ", "))
      }
      if (!"scale" %in% names(dt))     dt[, scale := "inversion_block"]
      if (!"parent_id" %in% names(dt)) dt[, parent_id := NA_character_]
      dt[, size_kb := round((end_bp - start_bp) / 1000, 1)]
      dt[, start_bp := as.integer(start_bp)]
      dt[, end_bp   := as.integer(end_bp)]

      c <- read_cands()
      existing <- if (nrow(c) > 0) c$candidate_id else character(0)
      new_rows <- dt[!(candidate_id %in% existing)]
      if (nrow(new_rows) == 0) return(invisible(0L))

      # Reorder new_rows cols to match existing cand_file schema if present
      cols <- names(c)
      if (length(cols) > 0) {
        extra <- setdiff(names(new_rows), cols)
        missing_c <- setdiff(cols, names(new_rows))
        for (m in missing_c) new_rows[, (m) := NA]
        new_rows <- new_rows[, c(cols, extra), with = FALSE]
      }
      c <- rbind(c, new_rows, fill = TRUE)
      fwrite(c, cand_file, sep = "\t")
      invisible(nrow(new_rows))
    }
  )
}

# =============================================================================
# Evidence API — Tier 1/2/3 storage, with schema validation for Tier 2
# =============================================================================
load_evidence_api <- function(evid_dir, schema_dir) {
  percand_dir <- file.path(evid_dir, "per_candidate")
  global_dir  <- file.path(evid_dir, "global")
  schema_blocks <- file.path(schema_dir, "structured_block_schemas")

  # Cache schemas
  load_schema <- function(block_type) {
    sf <- file.path(schema_blocks, paste0(block_type, ".schema.json"))
    if (!file.exists(sf)) return(NULL)
    tryCatch(jsonlite::fromJSON(sf, simplifyVector = FALSE),
             error = function(e) NULL)
  }

  cand_dir <- function(cid) {
    d <- file.path(percand_dir, cid)
    for (sub in c("raw", "structured")) {
      dir.create(file.path(d, sub), recursive = TRUE, showWarnings = FALSE)
    }
    d
  }

  keys_file <- function(cid) file.path(cand_dir(cid), "keys.tsv")

  # ── write_block: validate data against schema, write JSON, extract keys ──
  write_block <- function(candidate_id, block_type, data) {
    schema <- load_schema(block_type)
    d <- cand_dir(candidate_id)
    block_path <- file.path(d, "structured", paste0(block_type, ".json"))

    # Build the block with standard header
    block <- list(
      block_type     = block_type,
      candidate_id   = candidate_id,
      source_script  = Sys.getenv("CURRENT_SCRIPT", "unknown"),
      computed_at    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      data           = data
    )

    # Validate against schema if available
    if (!is.null(schema) && !is.null(schema$required)) {
      missing_req <- setdiff(unlist(schema$required), names(data))
      if (length(missing_req) > 0) {
        warning("[registry] block ", block_type, " missing required fields: ",
                paste(missing_req, collapse = ", "),
                " — writing anyway but flagging incomplete")
        block$validation_status <- "incomplete"
      } else {
        block$validation_status <- "validated"
      }
    } else {
      block$validation_status <- "no_schema"
    }

    # Write JSON
    jsonlite::write_json(block, block_path, auto_unbox = TRUE, pretty = TRUE,
                          null = "null", na = "null")

    # Extract flat keys per schema's keys_extracted directive
    # chat 13 (Finding BG): `from` may be a dotted path (e.g.
    # "boundary_sharpness.core_delta12" or "cohort.fraction_samples_R_fired");
    # walk the path via resolve_dotted() rather than single-level [[from]].
    # Vector-valued leaves (e.g. persistent_carriers) collapse to length
    # via length(); atomic scalars pass through unchanged.
    resolve_dotted <- function(obj, path) {
      parts <- strsplit(path, ".", fixed = TRUE)[[1]]
      for (p in parts) {
        if (is.null(obj)) return(NULL)
        # data.table-unfriendly: use [[ on lists, data.frames, environments
        obj <- tryCatch(obj[[p]], error = function(e) NULL)
      }
      obj
    }
    n_extracted <- 0L
    if (!is.null(schema) && !is.null(schema$keys_extracted)) {
      kex <- schema$keys_extracted
      flat <- data.table()
      for (ke in kex) {
        # Each entry: list(key = "q7_layer_d_fisher_or", from = "fisher_or")
        key <- ke$key; from <- ke$from
        if (is.null(key) || is.null(from)) next
        val <- resolve_dotted(data, from) %||% NA
        # If the resolved value is a vector with length > 1, collapse to
        # its length (matches existing schema intent — e.g. persistent_carriers
        # is an array; the key wants the count, not the array itself).
        if (!is.null(val) && length(val) > 1) val <- length(val)
        if (is.null(val) || (length(val) == 1 && is.na(val))) next
        flat <- rbind(flat, data.table(
          key = key, value = as.character(val),
          source_script = block$source_script,
          scope = "candidate",
          candidate_id = candidate_id,
          block_type = block_type,
          timestamp = block$computed_at
        ), fill = TRUE)
        n_extracted <- n_extracted + 1L
      }
      if (nrow(flat) > 0) {
        kf <- keys_file(candidate_id)
        if (file.exists(kf)) {
          existing <- fread(kf)
          # Replace any rows with same (key, candidate_id) — last-write-wins
          existing <- existing[!(existing$key %in% flat$key)]
          flat <- rbind(existing, flat, fill = TRUE)
        }
        fwrite(flat, kf, sep = "\t")
      }
    }

    message("[registry] wrote block ", block_type, " for ", candidate_id,
            " (", n_extracted, " keys extracted)")
    invisible(list(path = block_path, n_keys = n_extracted,
                    status = block$validation_status))
  }

  # ── read_block: load JSON ───────────────────────────────────────────────
  read_block <- function(candidate_id, block_type) {
    bp <- file.path(cand_dir(candidate_id), "structured",
                     paste0(block_type, ".json"))
    if (!file.exists(bp)) return(NULL)
    jsonlite::fromJSON(bp, simplifyVector = FALSE)
  }

  # ── get_keys: return flat keys.tsv for a candidate ──────────────────────
  get_keys <- function(candidate_id, key = NULL) {
    kf <- keys_file(candidate_id)
    if (!file.exists(kf)) return(data.table())
    k <- fread(kf)
    if (!is.null(key)) k <- k[k$key == key]
    k
  }

  # ── add_evidence / has_evidence / get_evidence (flat keys, v9 compat) ───
  add_evidence <- function(candidate_id, key, value = "",
                             file_path = "", script = "") {
    kf <- keys_file(candidate_id)
    dt <- if (file.exists(kf)) fread(kf) else data.table()
    dt <- if (nrow(dt) > 0) dt[dt$key != key] else dt
    dt <- rbind(dt, data.table(
      key = key, value = as.character(value),
      source_script = script, scope = "candidate",
      candidate_id = candidate_id, block_type = NA_character_,
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    ), fill = TRUE)
    fwrite(dt, kf, sep = "\t")
    invisible(TRUE)
  }

  has_evidence <- function(candidate_id, key) {
    k <- get_keys(candidate_id, key)
    nrow(k) > 0
  }

  get_evidence <- function(candidate_id, key = NULL) {
    get_keys(candidate_id, key)
  }

  register_candidate <- function(candidate_id, chrom, start_bp, end_bp,
                                   tier = NA, score = NA) {
    d <- cand_dir(candidate_id)  # creates directories
    # Register in interval_registry via add_candidate
    # (intervals API is in a different scope — skip here; scripts call both)
    add_evidence(candidate_id, "scored",
                  value = paste0("tier=", tier, " score=", round(score, 3)),
                  script = "register_candidate")
    message("[registry] candidate registered: ", candidate_id, " (T", tier, ")")
    invisible(TRUE)
  }

  list_candidates = function() {
    ds <- list.dirs(percand_dir, recursive = FALSE, full.names = FALSE)
    ds
  }

  count_candidates = function() length(list_candidates())

  list(
    write_block        = write_block,
    read_block         = read_block,
    get_keys           = get_keys,
    add_evidence       = add_evidence,
    has_evidence       = has_evidence,
    get_evidence       = get_evidence,
    register_candidate = register_candidate,
    list_candidates    = list_candidates,
    count_candidates   = count_candidates
  )
}

# =============================================================================
# Query API — composite reads across registries
# =============================================================================
load_query_api <- function(samples, intervals, evidence) {
  list(
    # The canonical example from handoff §7
    boundary_flanking_groups = function(candidate_id, side, offset_bp,
                                           flank_size = 10000) {
      bounds <- intervals$get_candidate_boundaries(candidate_id)
      if (is.null(bounds)) return(NULL)
      anchor <- if (side == "left") bounds$left_bp else bounds$right_bp
      pos <- anchor + offset_bp * (if (side == "right") 1 else -1)

      chrom <- intervals$get_candidate(candidate_id)$chrom
      win <- intervals$get_window_at(chrom, pos)

      list(
        position   = pos,
        window_id  = if (!is.null(win)) win$window_id else NA,
        chrom      = chrom,
        groups     = list(
          common = samples$get_carriers(candidate_id, "HOM_REF"),
          het    = samples$get_carriers(candidate_id, "HET"),
          rare   = samples$get_carriers(candidate_id, "HOM_INV"),
          recomb = samples$get_carriers(candidate_id, "RECOMBINANT")
        ),
        n_per_group = c(
          common = length(samples$get_carriers(candidate_id, "HOM_REF")),
          het    = length(samples$get_carriers(candidate_id, "HET")),
          rare   = length(samples$get_carriers(candidate_id, "HOM_INV")),
          recomb = length(samples$get_carriers(candidate_id, "RECOMBINANT"))
        )
      )
    },

    inversion_carriers = function(candidate_id, validation = "UNCERTAIN") {
      # Check validation level meets minimum
      ev <- evidence$get_evidence(candidate_id, "q6_group_validation")
      current <- if (nrow(ev) > 0) ev$value[1] else "NONE"
      rank <- c(NONE = 0, UNCERTAIN = 1, SUPPORTED = 2, VALIDATED = 3, SUSPECT = -1)
      if (rank[current] < rank[validation]) {
        message("[query] ", candidate_id, " validation=", current,
                " < required=", validation, " — returning NULL")
        return(NULL)
      }
      list(
        hom_ref = samples$get_carriers(candidate_id, "HOM_REF"),
        het     = samples$get_carriers(candidate_id, "HET"),
        hom_inv = samples$get_carriers(candidate_id, "HOM_INV"),
        recomb  = samples$get_carriers(candidate_id, "RECOMBINANT"),
        validation_level = current
      )
    },

    validation_status = function(candidate_id) {
      ev <- evidence$get_evidence(candidate_id, "q6_group_validation")
      if (nrow(ev) > 0) ev$value[1] else "NONE"
    },

    completion_per_question = function(candidate_id) {
      # Very rough: just counts q1_*, q2_*, ... keys present
      all_keys <- evidence$get_keys(candidate_id)
      if (nrow(all_keys) == 0) {
        return(setNames(rep(0L, 7), paste0("Q", 1:7)))
      }
      prefixes <- sprintf("q%d_", 1:7)
      sapply(prefixes, function(p) sum(grepl(paste0("^", p), all_keys$key)))
    },

    candidate_summary = function(candidate_id) {
      list(
        candidate       = intervals$get_candidate(candidate_id),
        validation      = evidence$get_evidence(candidate_id, "q6_group_validation"),
        tier            = evidence$get_evidence(candidate_id, "q7_tier"),
        verdict         = evidence$get_evidence(candidate_id, "q7_verdict"),
        completion      = evidence$get_keys(candidate_id)[, .N, by = "source_script"]
      )
    },

    # ══ chat 11 composite extensions ═══════════════════════════════════
    # Cross-registry traversals using the new interval + sample methods.

    # ── nested_family_tree(root_cid) ─────────────────────────────────
    # Returns a nested list tree rooted at root_cid:
    #   list(cid = "A", parent = NA, chrom = ..., start_bp = ..., end_bp = ...,
    #        n_children = 2,
    #        children = list(
    #          list(cid = "B", parent = "A", ..., children = list(...)),
    #          list(cid = "D", parent = "A", ..., children = list())
    #        ))
    # Cycle-safe (via `seen` set). Returns NULL if root_cid not found.
    # Used by phase 4e to characterize nested systems ("parent of B,C"
    # vs. "leaf nested inside A").
    nested_family_tree = function(root_cid) {
      seen <- character(0)
      build <- function(cid) {
        if (cid %in% seen) return(NULL)  # cycle guard
        seen <<- c(seen, cid)
        cand <- intervals$get_candidate(cid)
        if (is.null(cand)) return(NULL)
        kids <- intervals$get_children(cid)
        list(
          cid        = cid,
          parent     = intervals$get_parent(cid),
          chrom      = cand$chrom,
          start_bp   = cand$start_bp,
          end_bp     = cand$end_bp,
          n_children = length(kids),
          children   = lapply(kids, build)
        )
      }
      build(root_cid)
    },

    # ── overlap_clusters(chrom) ──────────────────────────────────────
    # Connected-component grouping over "overlap" edges, restricted to
    # chrom. Two candidates share a cluster if their intervals overlap
    # (nested counts as overlap for this purpose — the goal is finding
    # the set of CIDs that any overlap-reconciliation analysis must
    # consider jointly). Returns a data.table with columns:
    #   candidate_id, cluster_id, cluster_size
    # Cluster IDs are 1-indexed within chromosome.
    # Useful for phase 4c competing-hypothesis detection: if H3_nested_
    # composite fires on two CIDs in the same cluster, investigator
    # should check whether they're the same event called twice.
    overlap_clusters = function(chrom) {
      # Enumerate candidates on this chromosome by querying overlap
      # against the maximal possible interval. get_overlapping handles
      # chromosome filtering; .Machine$integer.max is a safe upper bound
      # for any genomic coordinate.
      cids <- intervals$get_overlapping(chrom, 1L, .Machine$integer.max)
      if (length(cids) == 0) return(data.table(
        candidate_id = character(), cluster_id = integer(),
        cluster_size = integer()
      ))

      # Build adjacency via get_overlapping per CID
      adj <- setNames(vector("list", length(cids)), cids)
      for (cid in cids) {
        cand <- intervals$get_candidate(cid)
        if (is.null(cand)) next
        neigh <- intervals$get_overlapping(chrom,
                                            cand$start_bp, cand$end_bp,
                                            exclude_cid = cid)
        adj[[cid]] <- intersect(neigh, cids)
      }

      # Union-find via BFS
      cluster_of <- setNames(rep(NA_integer_, length(cids)), cids)
      cluster_id <- 0L
      for (cid in cids) {
        if (!is.na(cluster_of[cid])) next
        cluster_id <- cluster_id + 1L
        # BFS
        frontier <- cid
        while (length(frontier) > 0) {
          next_f <- character(0)
          for (f in frontier) {
            if (is.na(cluster_of[f])) {
              cluster_of[f] <- cluster_id
              next_f <- c(next_f, adj[[f]])
            }
          }
          frontier <- unique(next_f)
        }
      }

      out <- data.table(
        candidate_id = cids,
        cluster_id   = unname(cluster_of)
      )
      out[, cluster_size := .N, by = cluster_id]
      out[order(cluster_id, candidate_id)]
    },

    # ── sample_inversion_load(sid) ───────────────────────────────────
    # Returns list with counts of candidates sid is a carrier for,
    # broken down by status. Useful for Q6 lineage analysis and for
    # flagging samples that carry unusually many inversions (potential
    # karyotype-level chromosomal rearrangement or sample-quality issue).
    sample_inversion_load = function(sid) {
      groups <- samples$get_sample_groups(sid)
      # Filter to inversion groups (inv_<cid>_<STATUS>)
      inv_groups <- groups[grepl("^inv_.+_[A-Z_]+$", groups)]
      if (length(inv_groups) == 0) {
        return(list(
          sample_id = sid, total = 0L,
          HOM_REF = 0L, HET = 0L, HOM_INV = 0L,
          RECOMBINANT = 0L, RECOMBINANT_GC = 0L, RECOMBINANT_DCO = 0L,
          HOM_STD = 0L, candidates = character(0)
        ))
      }
      # Parse group names: inv_<cid>_<STATUS>
      # STATUS may be HOM_REF, HOM_INV, HET, RECOMBINANT, RECOMBINANT_GC,
      # RECOMBINANT_DCO, HOM_STD. Longer statuses first to avoid
      # greedy prefix confusion (though current cids never contain "_"
      # followed by a status suffix, we still parse explicitly).
      statuses <- c("RECOMBINANT_GC", "RECOMBINANT_DCO",
                    "RECOMBINANT", "HOM_REF", "HOM_INV", "HOM_STD", "HET")
      counts <- setNames(integer(length(statuses)), statuses)
      cids   <- character(0)
      for (g in inv_groups) {
        for (st in statuses) {
          sfx <- paste0("_", st)
          if (endsWith(g, sfx)) {
            counts[st] <- counts[st] + 1L
            cid <- sub(paste0("^inv_(.+)", sfx, "$"), "\\1", g)
            cids <- c(cids, cid)
            break
          }
        }
      }
      out <- list(
        sample_id       = sid,
        total           = length(inv_groups),
        HOM_REF         = unname(counts["HOM_REF"]),
        HET             = unname(counts["HET"]),
        HOM_INV         = unname(counts["HOM_INV"]),
        RECOMBINANT     = unname(counts["RECOMBINANT"]),
        RECOMBINANT_GC  = unname(counts["RECOMBINANT_GC"]),
        RECOMBINANT_DCO = unname(counts["RECOMBINANT_DCO"]),
        HOM_STD         = unname(counts["HOM_STD"]),
        candidates      = unique(cids)
      )
      out
    },

    # ── effective_karyotype(sample_id, chrom, pos) — chat 16 ──────────
    # Multi-scale karyotype resolver. Given a sample and a genomic
    # position, walks the candidate hierarchy (parent_id chain in
    # interval_registry) and returns the karyotype at the SMALLEST
    # containing candidate — plus the parent chain so the caller can
    # see how the labels differ at different scales.
    #
    # Example: CGA042 is HOM_INV for the 7.5 Mb parent LG12_17, but
    # HOM_REF inside the nested 50 kb double-crossover LG12_17_DCO_3.
    # At position 14.5 Mb inside DCO_3, this returns:
    #   list(candidate_id="LG12_17_DCO_3", karyotype="HOM_REF",
    #        parent_candidate="LG12_17", parent_karyotype="HOM_INV", ...)
    #
    # Returns NULL if no candidate contains the position or the sample
    # has no registered karyotype group in any containing candidate.
    effective_karyotype = function(sample_id, chrom, pos) {
      # Find all candidates containing this (chrom, pos) — no exclusion
      cand_ids <- tryCatch(
        intervals$get_overlapping(chrom, as.integer(pos), as.integer(pos)),
        error = function(e) character(0))
      if (length(cand_ids) == 0) return(NULL)

      # Sort by interval size (smallest first) so the most specific
      # candidate is at the head.
      sizes <- vapply(cand_ids, function(cid) {
        b <- tryCatch(intervals$get_candidate_boundaries(cid),
                      error = function(e) NULL)
        if (is.null(b)) return(Inf)
        as.numeric(b$right_bp - b$left_bp)
      }, numeric(1))
      cand_ids <- cand_ids[order(sizes)]

      # For each candidate, find sample_id's karyotype via group membership
      statuses <- c("HOM_REF", "HOM_INV", "HET",
                     "RECOMBINANT", "RECOMBINANT_GC", "RECOMBINANT_DCO",
                     "HOM_STD")
      .karyo_for <- function(cid, sid) {
        for (st in statuses) {
          gid <- sprintf("inv_%s_%s", cid, st)
          if (samples$has_group(gid)) {
            if (sid %in% samples$get_group(gid)) return(st)
          }
        }
        NA_character_
      }

      # Smallest candidate first — that's the effective label
      primary_cid   <- cand_ids[1]
      primary_karyo <- .karyo_for(primary_cid, sample_id)
      if (is.na(primary_karyo)) {
        # Sample has no registered karyotype in the smallest candidate —
        # try the next larger one
        for (cid in cand_ids[-1]) {
          k <- .karyo_for(cid, sample_id)
          if (!is.na(k)) {
            primary_cid <- cid; primary_karyo <- k; break
          }
        }
        if (is.na(primary_karyo)) return(NULL)
      }

      # Walk up the parent chain via interval_registry parent_id
      parent_chain <- list()
      cur <- tryCatch(intervals$get_candidate(primary_cid),
                      error = function(e) NULL)
      while (!is.null(cur) && !is.null(cur$parent_id) &&
              !is.na(cur$parent_id) && nzchar(cur$parent_id) &&
              cur$parent_id != ".") {
        p_cid <- cur$parent_id
        p_karyo <- .karyo_for(p_cid, sample_id)
        parent_chain[[length(parent_chain) + 1L]] <- list(
          candidate_id = p_cid, karyotype = p_karyo)
        cur <- tryCatch(intervals$get_candidate(p_cid),
                        error = function(e) NULL)
        if (length(parent_chain) > 20L) break   # safety
      }

      list(
        sample_id         = sample_id,
        chrom             = chrom,
        pos               = as.integer(pos),
        candidate_id      = primary_cid,
        karyotype         = primary_karyo,
        parent_candidate  = if (length(parent_chain))
                              parent_chain[[1]]$candidate_id else NA_character_,
        parent_karyotype  = if (length(parent_chain))
                              parent_chain[[1]]$karyotype else NA_character_,
        full_chain        = parent_chain,
        n_candidates_containing_pos = length(cand_ids)
      )
    },

    # ══ decompose-v2 + boundary_scan query methods ═══════════════════════
    # These expose the per-window class matrix and the staircase/soft-boundary
    # samples produced by STEP_C01i_decompose_v2.R and
    # STEP_C01g_boundary_scan.R. Added alongside chat 11 so downstream
    # scripts (multi_recomb S4, figure panels) have a uniform accessor.

    # ── per_window_class_matrix(cid) ─────────────────────────────────
    # Returns the samples × windows character matrix of class labels
    # (HOM_REF/HET/HOM_INV/AMBIGUOUS/NA) stored as a companion RDS by
    # decompose v2. Returns NULL if the RDS is missing.
    #
    # Column names are window bp (as character); rows are sample IDs.
    # This is the primary input for:
    #   - multi_recomb S1 (per-window class consistency scoring)
    #   - boundary_scan (staircase of carrier counts)
    #   - diagnostic plots ("class mosaic per sample")
    per_window_class_matrix = function(cid) {
      # Find the per_window_class.rds via two routes:
      #   1. Path stored in the internal_dynamics block (preferred)
      #   2. Convention: <decomp_out>/<cid>/per_window_class.rds
      blk <- tryCatch(evidence$read_block(cid, "internal_dynamics"),
                      error = function(e) NULL)
      rds_path <- NA_character_
      if (!is.null(blk) && !is.null(blk$data) &&
          !is.null(blk$data$per_window_class_rds)) {
        rds_path <- blk$data$per_window_class_rds
      }
      if (is.na(rds_path) || !nzchar(rds_path) || !file.exists(rds_path)) {
        # Convention fallback
        for (base in c(Sys.getenv("DECOMP_V2_OUT", ""),
                        "decomp_v2_out",
                        file.path(Sys.getenv("BASE", ""),
                                  "inversion-popgen-toolkit", "decomp_v2_out"))) {
          p <- file.path(base, cid, "per_window_class.rds")
          if (file.exists(p)) { rds_path <- p; break }
        }
      }
      if (is.na(rds_path) || !file.exists(rds_path)) return(NULL)
      tryCatch(readRDS(rds_path), error = function(e) NULL)
    },

    # ── boundary_scan(cid, side) ─────────────────────────────────────
    # Returns the left or right (or both) scan data.frames from the
    # boundary_scan Tier-2 block. Columns: window_idx, bp, n_hom_ref,
    # n_het, n_hom_inv, delta_exp_h, delta_q12, delta_q13.
    #
    # side ∈ {"left", "right", "both"}. "both" returns a list with
    # $left and $right.
    boundary_scan = function(cid, side = "both") {
      blk <- tryCatch(evidence$read_block(cid, "boundary_scan"),
                      error = function(e) NULL)
      if (is.null(blk) || is.null(blk$data)) return(NULL)
      to_dt <- function(raw) {
        if (is.null(raw) || length(raw) == 0) return(data.table())
        rbindlist(lapply(raw, as.data.table), fill = TRUE)
      }
      left_dt  <- to_dt(blk$data$left_scan)
      right_dt <- to_dt(blk$data$right_scan)
      switch(side,
        "left"  = left_dt,
        "right" = right_dt,
        "both"  = list(left = left_dt, right = right_dt,
                        left_median_bp  = blk$data$left_boundary_median_bp,
                        left_fuzz_bp    = blk$data$left_boundary_fuzz_bp,
                        right_median_bp = blk$data$right_boundary_median_bp,
                        right_fuzz_bp   = blk$data$right_boundary_fuzz_bp,
                        sharpness       = blk$data$boundary_sharpness)
      )
    },

    # ── soft_boundary_samples(cid) ──────────────────────────────────
    # Returns the character vector of samples flagged as soft-boundary
    # by C01g. These feed multi_recomb as signal S4 (a single crossover
    # at a boundary is a different event class than an interior
    # double-crossover; S4 lets the classifier separate them).
    soft_boundary_samples = function(cid) {
      blk <- tryCatch(evidence$read_block(cid, "boundary_scan"),
                      error = function(e) NULL)
      if (is.null(blk) || is.null(blk$data)) return(character(0))
      as.character(blk$data$soft_boundary_samples %||% character(0))
    },

    # ══ v11 additions: C01j/C01l/C01m/gc_detector block accessors ══════

    # ── regime_segments(cid) ──────────────────────────────────────────
    # Returns the list of regime segments from C01j for this candidate,
    # plus the transitions and summary flags. Used by multi_recomb (R
    # signal derivation) and phase 4e characterization (regime-state
    # labeling).
    regime_segments = function(cid) {
      blk <- tryCatch(evidence$read_block(cid, "regime_segments"),
                      error = function(e) NULL)
      if (is.null(blk) || is.null(blk$data)) return(NULL)
      seg_raw   <- blk$data$segments    %||% list()
      trans_raw <- blk$data$transitions %||% list()
      seg_dt   <- if (length(seg_raw))   rbindlist(lapply(seg_raw,   as.data.table), fill = TRUE) else data.table()
      trans_dt <- if (length(trans_raw)) rbindlist(lapply(trans_raw, as.data.table), fill = TRUE) else data.table()
      list(
        segments           = seg_dt,
        transitions        = trans_dt,
        n_segments         = blk$data$n_segments %||% 0L,
        dominant_state     = blk$data$dominant_state %||% NA_character_,
        regime_has_recomb  = isTRUE(blk$data$regime_has_recombinant_signal %||% FALSE),
        n_samples_changed  = blk$data$n_samples_with_regime_change %||% 0L
      )
    },

    # ── local_structure_segments(cid) ────────────────────────────────
    # Returns C01l's per-segment (left_flank / inv_left_half / inv_core /
    # inv_right_half / right_flank) delta12/entropy/ENA summary plus the
    # derived boundary_sharpness object. Feeds the combination rule's
    # B signal (boundary sharpness) and phase 4e characterization.
    local_structure_segments = function(cid) {
      blk <- tryCatch(evidence$read_block(cid, "local_structure_segments"),
                      error = function(e) NULL)
      if (is.null(blk) || is.null(blk$data)) return(NULL)
      seg_raw <- blk$data$segment_summary %||% list()
      class_raw <- blk$data$by_inversion_class %||% list()
      list(
        segment_summary    = if (length(seg_raw))   rbindlist(lapply(seg_raw,   as.data.table), fill = TRUE) else data.table(),
        by_inversion_class = if (length(class_raw)) rbindlist(lapply(class_raw, as.data.table), fill = TRUE) else data.table(),
        boundary_sharpness = blk$data$boundary_sharpness %||% list()
      )
    },

    # ── distance_concordance(cid) ────────────────────────────────────
    # Returns C01m's multi-scale concordance summary. Inversion_vs_family
    # score is the headline: high = inversion-carrier signal dominates,
    # low = family-LD signal.
    distance_concordance = function(cid) {
      blk <- tryCatch(evidence$read_block(cid, "distance_concordance"),
                      error = function(e) NULL)
      if (is.null(blk) || is.null(blk$data)) return(NULL)
      per_dist <- blk$data$per_distance_summary %||% list()
      list(
        distances_tested           = blk$data$distances_tested %||% integer(0),
        per_distance_summary       = if (length(per_dist)) rbindlist(lapply(per_dist, as.data.table), fill = TRUE) else data.table(),
        persistent_carriers        = as.character(blk$data$persistent_carriers %||% character(0)),
        decaying_pairs_count       = blk$data$decaying_pairs_count %||% 0L,
        inversion_vs_family_score  = blk$data$inversion_vs_family_score %||% NA_real_
      )
    },

    # ── gene_conversion_tracts(cid) ──────────────────────────────────
    # Returns per-candidate catalog of short dosage-excursion tracts.
    # These are NOT recombinants — they are auxiliary annotations that
    # enrich the per-sample record in multi_recomb output. The combination
    # rule uses n_tracts per sample to annotate but NOT to gate.
    gene_conversion_tracts = function(cid) {
      blk <- tryCatch(evidence$read_block(cid, "gene_conversion_tracts"),
                      error = function(e) NULL)
      if (is.null(blk) || is.null(blk$data)) return(NULL)
      per_sample <- blk$data$per_sample_summary %||% list()
      tracts     <- blk$data$tracts %||% list()
      list(
        per_sample_summary     = if (length(per_sample)) rbindlist(lapply(per_sample, as.data.table), fill = TRUE) else data.table(),
        tracts                 = if (length(tracts))     rbindlist(lapply(tracts,     as.data.table), fill = TRUE) else data.table(),
        total_tracts           = blk$data$total_tracts %||% 0L,
        total_samples_with_gc  = blk$data$total_samples_with_gc %||% 0L
      )
    }
  )
}

# =============================================================================
# Compute API — live stats via Engine B
# =============================================================================
load_results_api <- function(results_dir, samples = NULL, intervals = NULL) {
  # ═════════════════════════════════════════════════════════════════════════
  # Chat 16 — results_registry (fourth registry, chat-15 stats_cache rewrite)
  #
  # This is a first-class registry, at parity with sample/interval/evidence.
  # Every row in manifest.tsv has an explicit primary key (row_id = UUID),
  # explicit foreign keys into sample_registry (who_1, who_2) and
  # interval_registry (candidate_id), and a provenance block linking back
  # to the source script and input files. See registries/DATABASE_DESIGN.md
  # for the full design rationale.
  #
  # Layout:
  #   <results_dir>/
  #   ├── manifest.tsv                      ← the database index
  #   ├── pairwise/<stat>/<chrom>/<g1>__vs__<g2>.tsv.gz
  #   ├── candidate/<cid>/
  #   │   ├── Q_K<NN>.<group>.tsv.gz        per-sample × per-window Q matrix
  #   │   ├── F_K<NN>.tsv.gz                global F matrix at K
  #   │   └── meta.tsv                      per-candidate audit trail
  #   └── interval/<chrom>/<s>_<e>.<group>.<stat>.K<NN>.tsv.gz
  #
  # Key design points:
  #   - No hashes in filenames. Sample set identity is a registered group_id
  #     (e.g. all_226). The group's 'created' timestamp is the version.
  #   - Every put_* enforces FK constraints: who_1/who_2 must be a registered
  #     group, candidate_id must be in interval_registry. If samples or
  #     intervals are NULL (bootstrap mode), FK checks are skipped with a
  #     warning on first call.
  #   - Manifest rows carry full provenance (source_script, engine, config_hash,
  #     upstream_files) so any number can be traced back to how it was made.
  #   - reg$results$ask() is the single query method with where/who/what/kind
  #     filters. Shortcuts below (ask_what_at, etc.) wrap it.
  #   - reg$results$integrity_check() verifies FKs resolve, file paths exist,
  #     and group_versions match current group timestamps.
  # ═════════════════════════════════════════════════════════════════════════

  # ── Create subdirs ──────────────────────────────────────────────────────
  for (sub in c("pairwise", "candidate", "interval")) {
    dir.create(file.path(results_dir, sub),
               recursive = TRUE, showWarnings = FALSE)
  }
  manifest_file <- file.path(results_dir, "manifest.tsv")

  # Canonical manifest columns (order matters for append consistency)
  .MANIFEST_COLS <- c(
    "row_id", "kind",
    "chrom", "start_bp", "end_bp", "candidate_id",
    "group_1", "group_1_version",
    "group_2", "group_2_version",
    "stat", "K",
    "file",
    "n_rows", "n_cols", "sha256",
    "engine", "engine_version", "source_script",
    "run_id", "config_hash", "upstream_files",
    "timestamp"
  )

  # ── Helpers ─────────────────────────────────────────────────────────────

  # Generate a v4 UUID without requiring the uuid package. Uses R's RNG —
  # sufficient for row_id uniqueness within a cohort (not cryptographic).
  .gen_uuid <- function() {
    hex <- function(n) paste0(sample(c(0:9, letters[1:6]), n, replace = TRUE),
                               collapse = "")
    # Mask for RFC4122 v4: version nibble = 4, variant nibble in {8,9,a,b}
    sprintf("%s-%s-4%s-%s%s-%s",
            hex(8), hex(4), hex(3),
            sample(c("8","9","a","b"), 1), hex(3),
            hex(12))
  }

  # Sort a pair of group names so (a,b) and (b,a) hash to same filename
  .pair_name <- function(g1, g2) {
    pair <- sort(c(as.character(g1), as.character(g2)))
    paste0(pair[1], "__vs__", pair[2])
  }

  # Read the current 'created' timestamp for a registered group.
  # Used as the version stamp on who_*.group_version. Returns NA if samples
  # registry is not available or the group is not registered.
  .get_group_version <- function(group_id) {
    if (is.null(samples) || is.null(samples$list_groups)) return(NA_character_)
    tryCatch({
      g <- samples$list_groups()
      if (!nrow(g)) return(NA_character_)
      r <- g[g$group_id == group_id]
      if (!nrow(r) || !"created" %in% names(r)) return(NA_character_)
      as.character(r$created[1])
    }, error = function(e) NA_character_)
  }

  # FK resolver for a group_id. Errors (not warns) if samples registry is
  # present and the group is not registered. If samples is NULL (bootstrap),
  # returns NA silently.
  .check_group_fk <- function(group_id, context = "group") {
    if (is.null(group_id) || is.na(group_id) || !nzchar(group_id)) {
      stop(sprintf("[results_registry] %s group_id is empty/NA — ",
                    context),
           "refusing to write. Pass a registered group_id from sample_registry.")
    }
    if (!is.null(samples) && !is.null(samples$has_group)) {
      if (!samples$has_group(group_id)) {
        stop(sprintf("[results_registry] %s group_id '%s' is not registered ",
                      context, group_id),
             "in sample_registry. Register with reg$samples$add_group() ",
             "before writing to results_registry.")
      }
    }
    invisible(TRUE)
  }

  # FK resolver for a candidate_id. Errors if intervals registry is present
  # and the candidate is not registered.
  .check_candidate_fk <- function(cid) {
    if (is.null(cid) || is.na(cid) || !nzchar(cid)) {
      stop("[results_registry] candidate_id is empty/NA — refusing to write.")
    }
    if (!is.null(intervals) && !is.null(intervals$get_candidate)) {
      c <- tryCatch(intervals$get_candidate(cid), error = function(e) NULL)
      if (is.null(c)) {
        stop("[results_registry] candidate_id '", cid, "' is not registered ",
             "in interval_registry. Register with ",
             "reg$intervals$add_candidate() before writing.")
      }
    }
    invisible(TRUE)
  }

  # Cheap content sha256. Returns "" if digest/openssl not available.
  .sha256_file <- function(path) {
    if (!file.exists(path)) return("")
    tryCatch({
      if (requireNamespace("digest", quietly = TRUE)) {
        return(digest::digest(path, algo = "sha256", file = TRUE))
      }
      if (requireNamespace("openssl", quietly = TRUE)) {
        return(as.character(openssl::sha256(file(path))))
      }
      ""
    }, error = function(e) "")
  }

  # Atomic append of one row to manifest. Creates with full header if empty.
  .append_manifest_row <- function(row_list) {
    # Enforce canonical column order and default missing columns to NA
    for (cc in .MANIFEST_COLS) {
      if (!(cc %in% names(row_list))) {
        row_list[[cc]] <- switch(cc,
          start_bp = NA_integer_, end_bp = NA_integer_,
          n_rows = NA_integer_, n_cols = NA_integer_, K = NA_integer_,
          NA_character_)
      }
    }
    row_dt <- as.data.table(row_list[.MANIFEST_COLS])
    if (file.exists(manifest_file)) {
      # Existing file — verify header matches; if not, gracefully rewrite.
      hdr <- tryCatch(names(fread(manifest_file, nrows = 0)),
                       error = function(e) character(0))
      if (!identical(hdr, .MANIFEST_COLS)) {
        # Legacy manifest (chat-15 stats_cache shape). Back up and start fresh;
        # legacy rows are preserved in the backup for audit purposes.
        bk <- paste0(manifest_file, ".legacy_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        file.rename(manifest_file, bk)
        message("[results_registry] legacy manifest backed up to ", bk,
                " — new chat-16 manifest starting fresh")
        fwrite(row_dt, manifest_file, sep = "\t")
        return(invisible(TRUE))
      }
      fwrite(row_dt, manifest_file, sep = "\t", append = TRUE)
    } else {
      fwrite(row_dt, manifest_file, sep = "\t")
    }
    invisible(TRUE)
  }

  # Read full manifest. Returns empty data.table if file missing.
  .read_manifest <- function() {
    if (!file.exists(manifest_file)) {
      return(data.table(setNames(
        replicate(length(.MANIFEST_COLS), character(0), simplify = FALSE),
        .MANIFEST_COLS)))
    }
    fread(manifest_file, colClasses = list(character = c(
      "row_id","kind","chrom","candidate_id",
      "group_1","group_1_version","group_2","group_2_version",
      "stat","file","sha256","engine","engine_version","source_script",
      "run_id","config_hash","upstream_files","timestamp"
    )))
  }

  # ── Writers ─────────────────────────────────────────────────────────────

  # put_pairwise(chrom, group_1, group_2, stat, dt, source_script=..., ...)
  # Writes a pairwise table (FST/dxy/etc.) between two registered groups.
  put_pairwise <- function(chrom, group_1, group_2, stat, dt,
                            source_script = "unknown",
                            engine = NULL, engine_version = NULL,
                            run_id = NULL, config_hash = NULL,
                            upstream_files = NULL,
                            # Back-compat: legacy sample_set arg silently ignored
                            sample_set = NULL,
                            # Back-compat: accept g1 / g2 / sample_group aliases
                            g1 = NULL, g2 = NULL) {
    # Back-compat: chat-15 callers used g1= / g2= and sample_set=<hash>.
    # New callers use group_1= / group_2= (semantic names).
    if (is.null(group_1) && !is.null(g1)) group_1 <- g1
    if (is.null(group_2) && !is.null(g2)) group_2 <- g2

    .check_group_fk(group_1, "group_1")
    .check_group_fk(group_2, "group_2")

    outdir <- file.path(results_dir, "pairwise", stat, chrom)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    pair <- .pair_name(group_1, group_2)
    out_file <- file.path(outdir, paste0(pair, ".tsv.gz"))
    fwrite(dt, out_file, sep = "\t")

    v1 <- .get_group_version(group_1)
    v2 <- .get_group_version(group_2)

    .append_manifest_row(list(
      row_id = .gen_uuid(), kind = "pairwise",
      chrom = chrom, candidate_id = NA_character_,
      group_1 = group_1, group_1_version = v1,
      group_2 = group_2, group_2_version = v2,
      stat = stat, K = NA_integer_,
      file = sub(paste0("^", results_dir, "/"), "", out_file),
      n_rows = nrow(dt), n_cols = ncol(dt),
      sha256 = .sha256_file(out_file),
      engine = engine %||% NA_character_,
      engine_version = engine_version %||% NA_character_,
      source_script = source_script,
      run_id = run_id %||% Sys.getenv("SLURM_JOB_ID", ""),
      config_hash = config_hash %||% NA_character_,
      upstream_files = if (is.null(upstream_files)) NA_character_
                        else paste(upstream_files, collapse = ";"),
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    ))

    invisible(out_file)
  }

  # get_pairwise(chrom, group_1, group_2, stat) — reads the tsv.gz back.
  # With back-compat for legacy sample_set arg (ignored).
  get_pairwise <- function(chrom, group_1, group_2, stat, sample_set = NULL,
                            g1 = NULL, g2 = NULL) {
    if (is.null(group_1) && !is.null(g1)) group_1 <- g1
    if (is.null(group_2) && !is.null(g2)) group_2 <- g2
    pair <- .pair_name(group_1, group_2)
    fn_new <- file.path(results_dir, "pairwise", stat, chrom,
                         paste0(pair, ".tsv.gz"))
    if (file.exists(fn_new)) return(fread(fn_new))
    # Back-compat: chat-15 files with a hash tag between pair and .tsv.gz
    outdir <- file.path(results_dir, "pairwise", stat, chrom)
    if (dir.exists(outdir)) {
      hits <- list.files(outdir,
                          pattern = sprintf("^%s\\..*\\.tsv\\.gz$", pair),
                          full.names = TRUE)
      if (length(hits)) {
        message("[results_registry] using legacy-named pairwise file ",
                basename(hits[1]), " — consider migrating")
        return(fread(hits[1]))
      }
    }
    NULL
  }

  # put_candidate_q_f(cid, q_dt, f_mat, K, sample_group) — writes Q and/or F
  # matrices for a candidate. Q is per-group (needs who_1); F has no group
  # dimension (one F per K). Either q_dt or f_mat can be NULL.
  put_candidate_q_f <- function(cid, q_dt = NULL, f_mat = NULL,
                                  K = NULL, sample_group = NULL,
                                  source_script = "unknown",
                                  engine = "instant_q", engine_version = NULL,
                                  run_id = NULL, config_hash = NULL,
                                  upstream_files = NULL) {
    .check_candidate_fk(cid)

    canon_K <- tryCatch(
      if (exists(".iq_env") && !is.null(.iq_env$canonical_K))
        .iq_env$canonical_K else 8L,
      error = function(e) 8L)
    K_eff <- as.integer(K %||% canon_K)

    sample_group <- sample_group %||%
      Sys.getenv("SAMPLE_GROUP", "all_226")
    if (!is.null(q_dt)) .check_group_fk(sample_group, "sample_group")

    outdir <- file.path(results_dir, "candidate", cid)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    q_file <- NULL
    f_file <- NULL
    group_ver <- .get_group_version(sample_group)
    run_id_resolved <- run_id %||% Sys.getenv("SLURM_JOB_ID", "")
    ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

    # Q matrix (per-group)
    if (!is.null(q_dt) && (is.data.frame(q_dt) ||
                            data.table::is.data.table(q_dt))) {
      q_file <- file.path(outdir,
                           sprintf("Q_K%02d.%s.tsv.gz", K_eff, sample_group))
      fwrite(q_dt, q_file, sep = "\t")
      .append_manifest_row(list(
        row_id = .gen_uuid(), kind = "candidate_q",
        chrom = NA_character_, candidate_id = cid,
        group_1 = sample_group, group_1_version = group_ver,
        group_2 = NA_character_, group_2_version = NA_character_,
        stat = "Q", K = K_eff,
        file = sub(paste0("^", results_dir, "/"), "", q_file),
        n_rows = nrow(q_dt), n_cols = ncol(q_dt),
        sha256 = .sha256_file(q_file),
        engine = engine, engine_version = engine_version %||% NA_character_,
        source_script = source_script,
        run_id = run_id_resolved, config_hash = config_hash %||% NA_character_,
        upstream_files = if (is.null(upstream_files)) NA_character_
                          else paste(upstream_files, collapse = ";"),
        timestamp = ts))
    }

    # F matrix (no group — one per K)
    if (!is.null(f_mat) && is.matrix(f_mat)) {
      f_file <- file.path(outdir, sprintf("F_K%02d.tsv.gz", K_eff))
      fwrite(as.data.table(f_mat), f_file, sep = "\t")
      .append_manifest_row(list(
        row_id = .gen_uuid(), kind = "candidate_f",
        chrom = NA_character_, candidate_id = cid,
        group_1 = NA_character_, group_1_version = NA_character_,
        group_2 = NA_character_, group_2_version = NA_character_,
        stat = "F", K = K_eff,
        file = sub(paste0("^", results_dir, "/"), "", f_file),
        n_rows = nrow(f_mat), n_cols = ncol(f_mat),
        sha256 = .sha256_file(f_file),
        engine = engine, engine_version = engine_version %||% NA_character_,
        source_script = source_script,
        run_id = run_id_resolved, config_hash = config_hash %||% NA_character_,
        upstream_files = if (is.null(upstream_files)) NA_character_
                          else paste(upstream_files, collapse = ";"),
        timestamp = ts))
    }

    # Per-candidate audit trail (independent of central manifest; survives
    # manifest resets for debugging)
    meta_file <- file.path(outdir, "meta.tsv")
    meta_row <- data.table(
      cid = cid, K = K_eff, sample_group = sample_group,
      q_file = if (!is.null(q_file)) basename(q_file) else NA_character_,
      f_file = if (!is.null(f_file)) basename(f_file) else NA_character_,
      n_rows_Q = if (!is.null(q_dt)) nrow(q_dt) else 0L,
      source_script = source_script, run_id = run_id_resolved, timestamp = ts)
    if (file.exists(meta_file)) {
      fwrite(meta_row, meta_file, sep = "\t", append = TRUE)
    } else {
      fwrite(meta_row, meta_file, sep = "\t")
    }

    invisible(list(q_file = q_file, f_file = f_file,
                    row_ids = character(0)))
  }

  # put_interval_summary — per-interval ancestry summary (delta12/entropy/ena/etc.)
  put_interval_summary <- function(chrom, start_bp, end_bp, group_1, stat, dt,
                                     K = NULL,
                                     source_script = "unknown",
                                     engine = NULL, engine_version = NULL,
                                     run_id = NULL, config_hash = NULL,
                                     upstream_files = NULL) {
    .check_group_fk(group_1, "group_1")
    valid_stats <- c("delta12", "entropy", "ena", "ancestry_q_mean")
    if (!(stat %in% valid_stats)) {
      stop("[results_registry] stat '", stat, "' not valid for interval_summary. ",
           "Use one of: ", paste(valid_stats, collapse = ", "))
    }
    outdir <- file.path(results_dir, "interval", chrom)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    K_suffix <- if (is.null(K)) "" else sprintf(".K%02d", as.integer(K))
    out_file <- file.path(outdir,
                           sprintf("%d_%d.%s.%s%s.tsv.gz",
                                   as.integer(start_bp), as.integer(end_bp),
                                   group_1, stat, K_suffix))
    fwrite(dt, out_file, sep = "\t")
    v1 <- .get_group_version(group_1)
    .append_manifest_row(list(
      row_id = .gen_uuid(), kind = "interval_summary",
      chrom = chrom, candidate_id = NA_character_,
      start_bp = as.integer(start_bp), end_bp = as.integer(end_bp),
      group_1 = group_1, group_1_version = v1,
      group_2 = NA_character_, group_2_version = NA_character_,
      stat = stat, K = if (is.null(K)) NA_integer_ else as.integer(K),
      file = sub(paste0("^", results_dir, "/"), "", out_file),
      n_rows = nrow(dt), n_cols = ncol(dt),
      sha256 = .sha256_file(out_file),
      engine = engine %||% NA_character_,
      engine_version = engine_version %||% NA_character_,
      source_script = source_script,
      run_id = run_id %||% Sys.getenv("SLURM_JOB_ID", ""),
      config_hash = config_hash %||% NA_character_,
      upstream_files = if (is.null(upstream_files)) NA_character_
                        else paste(upstream_files, collapse = ";"),
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
    invisible(out_file)
  }

  # ── Readers ─────────────────────────────────────────────────────────────

  get_candidate_q <- function(cid, K = NULL, sample_group = NULL,
                                sample_set = NULL) {
    canon_K <- tryCatch(
      if (exists(".iq_env") && !is.null(.iq_env$canonical_K))
        .iq_env$canonical_K else 8L,
      error = function(e) 8L)
    K_eff <- as.integer(K %||% canon_K)
    sg <- sample_group %||% Sys.getenv("SAMPLE_GROUP", "all_226")
    fn <- file.path(results_dir, "candidate", cid,
                     sprintf("Q_K%02d.%s.tsv.gz", K_eff, sg))
    if (file.exists(fn)) return(fread(fn))
    # Back-compat: chat-15 hash-tagged files
    outdir <- file.path(results_dir, "candidate", cid)
    if (dir.exists(outdir)) {
      hits <- list.files(outdir,
                          pattern = sprintf("^Q_K%02d\\..*\\.tsv\\.gz$", K_eff),
                          full.names = TRUE)
      if (length(hits)) {
        message("[results_registry] using legacy-named Q file ",
                basename(hits[1]), " for ", cid,
                " — consider migrating")
        return(fread(hits[1]))
      }
    }
    NULL
  }

  get_candidate_f <- function(cid, K = NULL) {
    canon_K <- tryCatch(
      if (exists(".iq_env") && !is.null(.iq_env$canonical_K))
        .iq_env$canonical_K else 8L,
      error = function(e) 8L)
    K_eff <- as.integer(K %||% canon_K)
    f_file <- file.path(results_dir, "candidate", cid,
                         sprintf("F_K%02d.tsv.gz", K_eff))
    if (file.exists(f_file)) return(as.matrix(fread(f_file)))
    NULL
  }

  # ── Query plane ─────────────────────────────────────────────────────────

  # The one function. Four filter dimensions plus overlap semantics.
  ask <- function(where = NULL, who = NULL, what = NULL, kind = NULL,
                   K = NULL, overlap = "any") {
    m <- .read_manifest()
    if (!nrow(m)) return(m)

    # kind filter
    if (!is.null(kind)) m <- m[m$kind %in% kind]
    if (!nrow(m)) return(m)

    # what filter (stat)
    if (!is.null(what)) m <- m[m$stat %in% what]
    if (!nrow(m)) return(m)

    # K filter
    if (!is.null(K)) m <- m[!is.na(m$K) & m$K %in% as.integer(K)]
    if (!nrow(m)) return(m)

    # who filter — can be one group or a vector of groups; matches if
    # either group_1 or group_2 is in the set
    if (!is.null(who)) {
      who <- as.character(who)
      m <- m[m$group_1 %in% who | m$group_2 %in% who]
    }
    if (!nrow(m)) return(m)

    # where filter — candidate_id, chrom, and interval overlap
    if (!is.null(where)) {
      if (!is.null(where$candidate_id)) {
        m <- m[!is.na(m$candidate_id) & m$candidate_id == where$candidate_id]
      }
      if (!is.null(where$chrom)) {
        # Match either direct chrom rows, or candidate rows whose candidate
        # is on that chrom. For the latter we need interval_registry.
        direct <- m[!is.na(m$chrom) & m$chrom == where$chrom]
        via_cand <- data.table()
        if (!is.null(intervals)) {
          cand_rows <- m[!is.na(m$candidate_id)]
          if (nrow(cand_rows)) {
            cand_rows$.wc <- vapply(cand_rows$candidate_id, function(cid) {
              c <- tryCatch(intervals$get_candidate(cid), error = function(e) NULL)
              if (is.null(c)) return(NA_character_)
              if (is.list(c) && !is.data.frame(c)) c$chrom %||% NA_character_
              else c$chrom[1] %||% NA_character_
            }, character(1))
            via_cand <- cand_rows[.wc == where$chrom]
            via_cand$.wc <- NULL
          }
        }
        m <- unique(rbindlist(list(direct, via_cand), fill = TRUE))
      }
      if (!is.null(where$start_bp) && !is.null(where$end_bp)) {
        qs <- as.integer(where$start_bp); qe <- as.integer(where$end_bp)
        # For rows with explicit start_bp/end_bp, test overlap directly.
        # For rows keyed on candidate_id, look up the candidate's coords.
        keep <- logical(nrow(m))
        for (i in seq_len(nrow(m))) {
          rs <- suppressWarnings(as.integer(m$start_bp[i]))
          re <- suppressWarnings(as.integer(m$end_bp[i]))
          cid_i <- m$candidate_id[i]
          if ((is.na(rs) || is.na(re)) && !is.na(cid_i) && nzchar(cid_i) &&
              !is.null(intervals)) {
            b <- tryCatch(intervals$get_candidate_boundaries(cid_i),
                           error = function(e) NULL)
            if (!is.null(b)) { rs <- b$left_bp; re <- b$right_bp }
          }
          if (is.na(rs) || is.na(re)) next
          keep[i] <- switch(overlap,
            "any"          = !(re < qs || rs > qe),
            "contains"     = rs <= qs && re >= qe,   # row contains query
            "contained_by" = rs >= qs && re <= qe)
        }
        m <- m[keep]
      }
    }
    m
  }

  # Shortcuts

  ask_what_at <- function(chrom, start_bp = NULL, end_bp = NULL) {
    w <- list(chrom = chrom)
    if (!is.null(start_bp)) w$start_bp <- start_bp
    if (!is.null(end_bp))   w$end_bp   <- end_bp
    ask(where = w)
  }

  ask_what_for_candidate <- function(cid) {
    ask(where = list(candidate_id = cid))
  }

  ask_what_for_group <- function(gid) {
    ask(who = gid)
  }

  ask_provenance <- function(row_id) {
    m <- .read_manifest()
    r <- m[m$row_id == row_id]
    if (!nrow(r)) return(NULL)
    cols <- c("row_id", "kind", "source_script", "engine", "engine_version",
              "run_id", "config_hash", "upstream_files", "timestamp")
    cols <- intersect(cols, names(r))
    as.list(r[1, cols, with = FALSE])
  }

  # ── Integrity check ─────────────────────────────────────────────────────

  integrity_check <- function(check_sha256 = FALSE) {
    m <- .read_manifest()
    checks <- list()

    # 1. FK: who_1 / who_2 must exist in sample_registry
    n_viol <- 0L
    details <- character(0)
    if (!is.null(samples) && nrow(m)) {
      for (col in c("group_1", "group_2")) {
        vals <- m[[col]]
        present <- !is.na(vals) & nzchar(vals)
        if (!any(present)) next
        unique_gids <- unique(vals[present])
        for (gid in unique_gids) {
          if (!samples$has_group(gid)) {
            n_viol <- n_viol + sum(vals == gid, na.rm = TRUE)
            details <- c(details, sprintf("%s='%s' unregistered", col, gid))
          }
        }
      }
    }
    checks[[length(checks) + 1L]] <- list(
      check = "fk_group_id_resolves", n_violations = n_viol,
      pass = n_viol == 0L,
      details = if (!length(details)) "" else paste(details, collapse = "; "))

    # 2. FK: candidate_id must exist in interval_registry
    n_viol <- 0L; details <- character(0)
    if (!is.null(intervals) && nrow(m)) {
      cids <- unique(m$candidate_id[!is.na(m$candidate_id) &
                                        nzchar(m$candidate_id)])
      for (cid in cids) {
        c <- tryCatch(intervals$get_candidate(cid), error = function(e) NULL)
        if (is.null(c)) {
          n_viol <- n_viol + sum(m$candidate_id == cid, na.rm = TRUE)
          details <- c(details, sprintf("candidate_id='%s' unregistered", cid))
        }
      }
    }
    checks[[length(checks) + 1L]] <- list(
      check = "fk_candidate_id_resolves", n_violations = n_viol,
      pass = n_viol == 0L,
      details = if (!length(details)) "" else paste(details, collapse = "; "))

    # 3. group_version matches current 'created' timestamp (stale check)
    n_viol <- 0L; details <- character(0)
    if (!is.null(samples) && nrow(m)) {
      for (pair in list(c("group_1","group_1_version"),
                        c("group_2","group_2_version"))) {
        gcol <- pair[1]; vcol <- pair[2]
        rows <- m[!is.na(m[[gcol]]) & nzchar(m[[gcol]])]
        for (i in seq_len(nrow(rows))) {
          gid <- rows[[gcol]][i]; stored_v <- rows[[vcol]][i]
          if (is.na(stored_v) || !nzchar(stored_v)) next
          current_v <- .get_group_version(gid)
          if (is.na(current_v) || current_v != stored_v) {
            n_viol <- n_viol + 1L
            details <- c(details,
              sprintf("row_id=%s: %s '%s' version %s != current %s",
                      rows$row_id[i], gcol, gid, stored_v, current_v))
          }
        }
      }
    }
    checks[[length(checks) + 1L]] <- list(
      check = "group_version_current", n_violations = n_viol,
      pass = n_viol == 0L,
      details = if (!length(details)) "" else
                  paste(head(details, 5), collapse = "; "))

    # 4. Every `file` in manifest exists on disk
    n_viol <- 0L; details <- character(0)
    if (nrow(m)) {
      for (i in seq_len(nrow(m))) {
        rel <- m$file[i]
        if (is.na(rel) || !nzchar(rel)) next
        full <- file.path(results_dir, rel)
        if (!file.exists(full)) {
          n_viol <- n_viol + 1L
          details <- c(details,
                        sprintf("row_id=%s: file %s missing", m$row_id[i], rel))
        }
      }
    }
    checks[[length(checks) + 1L]] <- list(
      check = "file_exists", n_violations = n_viol,
      pass = n_viol == 0L,
      details = if (!length(details)) "" else
                  paste(head(details, 5), collapse = "; "))

    # 5. sha256 matches (optional, expensive)
    if (isTRUE(check_sha256)) {
      n_viol <- 0L; details <- character(0)
      if (nrow(m)) {
        for (i in seq_len(nrow(m))) {
          rel <- m$file[i]; stored <- m$sha256[i]
          if (is.na(rel) || !nzchar(rel)) next
          if (is.na(stored) || !nzchar(stored)) next
          full <- file.path(results_dir, rel)
          if (!file.exists(full)) next
          cur <- .sha256_file(full)
          if (!nzchar(cur)) next
          if (cur != stored) {
            n_viol <- n_viol + 1L
            details <- c(details,
                          sprintf("row_id=%s: sha256 drift", m$row_id[i]))
          }
        }
      }
      checks[[length(checks) + 1L]] <- list(
        check = "sha256_matches", n_violations = n_viol,
        pass = n_viol == 0L,
        details = if (!length(details)) "" else
                    paste(head(details, 5), collapse = "; "))
    }

    # 6. No orphan files under results_dir/ that aren't in manifest
    n_viol <- 0L; details <- character(0)
    expected <- if (nrow(m)) file.path(results_dir, m$file[!is.na(m$file)])
                else character(0)
    for (sub in c("pairwise", "candidate", "interval")) {
      d <- file.path(results_dir, sub)
      if (!dir.exists(d)) next
      all_files <- list.files(d, recursive = TRUE, full.names = TRUE)
      # Ignore meta.tsv files — per-candidate audit, not a result artefact
      all_files <- all_files[!grepl("/meta\\.tsv$", all_files)]
      orphans <- setdiff(all_files, expected)
      if (length(orphans)) {
        n_viol <- n_viol + length(orphans)
        details <- c(details, paste(basename(head(orphans, 5)),
                                      collapse = ", "))
      }
    }
    checks[[length(checks) + 1L]] <- list(
      check = "no_orphan_files", n_violations = n_viol,
      pass = n_viol == 0L,
      details = if (!length(details)) "" else
                  paste(head(details, 5), collapse = "; "))

    out <- rbindlist(lapply(checks, as.data.table))
    attr(out, "all_pass") <- all(out$pass)
    out
  }

  # ── Maintenance ─────────────────────────────────────────────────────────

  # Returns a data.table identical to the manifest, filtered by kind.
  # (Matches chat-15 reg$stats$list_cached signature for back-compat.)
  list_cached <- function(kind = NULL) {
    m <- .read_manifest()
    if (!is.null(kind) && nrow(m)) m <- m[m$kind %in% kind]
    m
  }

  # Delete all files + manifest rows for one candidate.
  clear_candidate <- function(cid) {
    d <- file.path(results_dir, "candidate", cid)
    if (dir.exists(d)) unlink(d, recursive = TRUE)
    m <- .read_manifest()
    if (nrow(m)) {
      keep <- is.na(m$candidate_id) | m$candidate_id != cid
      if (any(!keep)) {
        fwrite(m[keep], manifest_file, sep = "\t")
      }
    }
    invisible(TRUE)
  }

  # Session summary — for the "echo before/after" pattern you wanted in
  # each pipeline script. Call at script top to stash a baseline, then at
  # script bottom to report what the script wrote.
  #
  # Usage:
  #   reg$results$session_start()      # at top
  #   # ... script runs put_* calls ...
  #   reg$results$session_summary()    # at bottom — prints a one-liner
  .session_baseline <- new.env(parent = emptyenv())
  session_start <- function() {
    m <- .read_manifest()
    .session_baseline$n_rows     <- nrow(m)
    .session_baseline$timestamp  <- Sys.time()
    invisible(NULL)
  }
  session_summary <- function(print = TRUE) {
    m <- .read_manifest()
    if (is.null(.session_baseline$n_rows)) {
      # No baseline — just report totals
      baseline <- 0L
      ts <- NA
    } else {
      baseline <- .session_baseline$n_rows
      ts <- .session_baseline$timestamp
    }
    new_rows <- if (nrow(m) > baseline) m[(baseline + 1L):nrow(m)] else m[0]
    summary <- list(
      n_rows_total_after   = nrow(m),
      n_rows_total_before  = baseline,
      n_rows_added         = nrow(new_rows),
      n_candidates_touched = if (nrow(new_rows))
                                length(unique(new_rows$candidate_id[
                                  !is.na(new_rows$candidate_id)])) else 0L,
      kinds_written        = if (nrow(new_rows))
                                table(new_rows$kind) else table(character(0)),
      elapsed_seconds      = if (inherits(ts, "POSIXct"))
                                round(as.numeric(
                                  difftime(Sys.time(), ts, units = "secs")), 1)
                              else NA
    )
    if (isTRUE(print)) {
      cat(sprintf("[results_registry] session wrote %d new rows",
                  summary$n_rows_added))
      if (summary$n_candidates_touched > 0L) {
        cat(sprintf(" across %d candidates", summary$n_candidates_touched))
      }
      if (!is.na(summary$elapsed_seconds)) {
        cat(sprintf(" in %.1fs", summary$elapsed_seconds))
      }
      if (length(summary$kinds_written) > 0L) {
        kinds_str <- paste(names(summary$kinds_written),
                            as.integer(summary$kinds_written),
                            sep = "=", collapse = ", ")
        cat(sprintf(" (%s)", kinds_str))
      }
      cat(sprintf(". Total manifest: %d rows.\n", summary$n_rows_total_after))
    }
    invisible(summary)
  }

  list(
    # Writers
    put_pairwise          = put_pairwise,
    put_candidate_q_f     = put_candidate_q_f,
    put_interval_summary  = put_interval_summary,
    # Readers
    get_pairwise          = get_pairwise,
    get_candidate_q       = get_candidate_q,
    get_candidate_f       = get_candidate_f,
    # Query plane
    ask                    = ask,
    ask_what_at            = ask_what_at,
    ask_what_for_candidate = ask_what_for_candidate,
    ask_what_for_group     = ask_what_for_group,
    ask_provenance         = ask_provenance,
    # Integrity
    integrity_check        = integrity_check,
    # Session echo (pre/post)
    session_start          = session_start,
    session_summary        = session_summary,
    # Maintenance
    list_cached            = list_cached,
    clear_candidate        = clear_candidate,
    # Paths
    results_dir            = results_dir,
    manifest_file          = manifest_file
  )
}

# ═════════════════════════════════════════════════════════════════════════
# Deprecated chat-15 alias. Kept for one chat cycle so old callers don't
# break mid-migration. Forwards to load_results_api with a one-line message.
# Will be removed in chat 17 once the pipeline has fully migrated.
# ═════════════════════════════════════════════════════════════════════════
load_stats_api <- function(stats_dir) {
  .Deprecated("load_results_api",
              msg = "load_stats_api is a chat-15 name; call load_results_api(results_dir, samples, intervals) instead.")
  # Return an API with samples=NULL, intervals=NULL — FK checks are skipped.
  # This is fine for legacy callers that never wired up the cross-registry
  # links. New code paths in load_registry() call load_results_api directly.
  load_results_api(stats_dir, samples = NULL, intervals = NULL)
}

load_compute_api <- function(samples, intervals, evidence) {
  # Requires get_region_stats() / resolve_groups() from load_bridge.R
  # If Engine B not available, methods return NULL with a warning.

  pairwise_stat <- function(position, flank_size, group1, group2, stat = "fst",
                               chrom = NULL, cache = FALSE) {
    if (!exists("get_region_stats", mode = "function")) {
      warning("[compute] Engine B not loaded — pairwise_stat NULL")
      return(NULL)
    }
    if (is.null(chrom)) {
      warning("[compute] chrom required")
      return(NULL)
    }
    res <- tryCatch({
      get_region_stats(chrom,
                          start_bp = position - flank_size,
                          end_bp   = position + flank_size,
                          what     = stat,
                          groups   = list(g1 = group1, g2 = group2))
    }, error = function(e) {
      warning("[compute] pairwise_stat failed: ", conditionMessage(e))
      NULL
    })
    # If asked to cache AND the stats API is in scope, persist the result.
    # Works when result is a data.table (the typical get_region_stats return
    # shape for per-window stats). Scalar results get wrapped into a 1-row dt.
    if (isTRUE(cache) && !is.null(res)) {
      stats_api <- tryCatch(get(".REG_STATS_API", envir = .GlobalEnv,
                                   inherits = FALSE),
                              error = function(e) NULL)
      if (!is.null(stats_api) && !is.null(stats_api$put_pairwise)) {
        out_dt <- if (data.table::is.data.table(res)) res
                    else if (is.list(res)) as.data.table(res)
                    else data.table::data.table(value = unlist(res))
        # g1/g2 may be raw sample vectors or group names. If they're long
        # (>50 entries), assume sample vector and use "custom_a"/"custom_b"
        # placeholders so the filename doesn't blow up.
        g1_label <- if (length(group1) > 50L || is.null(names(group1)))
                         "custom_g1" else paste(head(group1, 1), collapse = "")
        g2_label <- if (length(group2) > 50L || is.null(names(group2)))
                         "custom_g2" else paste(head(group2, 1), collapse = "")
        tryCatch(stats_api$put_pairwise(chrom, g1_label, g2_label, stat, out_dt),
                  error = function(e)
                    warning("[compute] put_pairwise cache failed: ",
                            conditionMessage(e)))
      }
    }
    res
  }

  # ===========================================================================
  # Chat 14 (2026-04-18): GHSL v6 on-demand queries
  # ===========================================================================
  # Thin wrappers around utils/lib_ghsl_panel.R. Pattern matches
  # pairwise_stat / boundary_fst_profile above: registry-aware callers
  # that resolve cid/chrom/bp from the interval_registry and delegate
  # to the backend library. If the library isn't sourced, each method
  # attempts to source it via GHSL_PANEL_LIB env var or a default
  # relative path, and returns NULL with a warning if not found.
  # ===========================================================================

  .ensure_ghsl_panel_lib <- function() {
    if (exists("ghsl_panel_range", mode = "function")) return(TRUE)
    libf <- Sys.getenv("GHSL_PANEL_LIB", "")
    if (!nzchar(libf) || !file.exists(libf)) {
      # Try a few standard locations
      for (cand in c(
        file.path(Sys.getenv("CODEBASE", ""), "utils", "lib_ghsl_panel.R"),
        file.path(Sys.getenv("BASE", ""),
                   "inversion_codebase_v8.5", "utils", "lib_ghsl_panel.R"),
        "utils/lib_ghsl_panel.R",
        "../../../utils/lib_ghsl_panel.R"
      )) {
        if (nzchar(cand) && file.exists(cand)) { libf <- cand; break }
      }
    }
    if (!nzchar(libf) || !file.exists(libf)) {
      warning("[compute] utils/lib_ghsl_panel.R not found — ",
              "set GHSL_PANEL_LIB env var or source it manually")
      return(FALSE)
    }
    source(libf, local = FALSE)
    TRUE
  }

  ghsl_at_interval <- function(chrom, start_bp, end_bp,
                                 sample_subset = NULL,
                                 scale = NULL,
                                 summaries = c("mean", "frac_low", "frac_mid",
                                               "frac_high", "longest_low_run_bp",
                                               "longest_high_run_bp",
                                               "rank_mean", "n_windows"),
                                 ghsl_dir = NULL) {
    if (!.ensure_ghsl_panel_lib()) return(NULL)
    tryCatch(
      ghsl_panel_aggregate(chrom, start_bp, end_bp,
                              sample_ids = sample_subset,
                              scale = scale, summaries = summaries,
                              ghsl_dir = ghsl_dir),
      error = function(e) {
        warning("[compute] ghsl_at_interval failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  ghsl_at_candidate <- function(cid,
                                  sample_subset = NULL,
                                  scale = NULL,
                                  summaries = c("mean", "frac_low", "frac_mid",
                                                "frac_high",
                                                "longest_low_run_bp",
                                                "longest_high_run_bp",
                                                "rank_mean", "n_windows"),
                                  ghsl_dir = NULL) {
    cand <- intervals$get_candidate(cid)
    if (is.null(cand)) {
      warning("[compute] ghsl_at_candidate: unknown cid ", cid)
      return(NULL)
    }
    # candidate row may be a list or a 1-row data.table; robust extract
    .pick <- function(x, field) {
      if (is.null(x)) return(NULL)
      if (is.list(x) && !is.data.frame(x) && field %in% names(x)) return(x[[field]])
      if (is.data.frame(x) && field %in% names(x)) return(x[[field]][1])
      NULL
    }
    chrom <- .pick(cand, "chrom")
    s_bp  <- .pick(cand, "start_bp") %||% .pick(cand, "start")
    e_bp  <- .pick(cand, "end_bp")   %||% .pick(cand, "end")
    if (is.null(chrom) || is.null(s_bp) || is.null(e_bp)) {
      warning("[compute] ghsl_at_candidate: candidate row missing chrom/start/end")
      return(NULL)
    }
    ghsl_at_interval(chrom, s_bp, e_bp,
                       sample_subset = sample_subset,
                       scale = scale, summaries = summaries,
                       ghsl_dir = ghsl_dir)
  }

  ghsl_at_subblocks <- function(cid_or_chrom, subblocks,
                                  sample_subset = NULL, scale = NULL,
                                  summaries = c("mean", "frac_low", "frac_mid",
                                                "frac_high", "rank_mean",
                                                "n_windows"),
                                  include_transitions = TRUE,
                                  ghsl_dir = NULL) {
    if (!.ensure_ghsl_panel_lib()) return(NULL)
    # Accept either a cid (looked up in interval registry for chrom) or
    # a chromosome name directly.
    chrom <- NULL
    cand <- tryCatch(intervals$get_candidate(cid_or_chrom),
                       error = function(e) NULL)
    if (!is.null(cand) && length(cand) > 0) {
      chrom <- if (is.list(cand) && !is.data.frame(cand)) cand$chrom
               else if (is.data.frame(cand) && "chrom" %in% names(cand)) cand$chrom[1]
               else NULL
    }
    if (is.null(chrom)) chrom <- cid_or_chrom
    tryCatch(
      ghsl_panel_subblock_scan(chrom, subblocks,
                                 sample_ids = sample_subset,
                                 scale = scale, summaries = summaries,
                                 include_transitions = include_transitions,
                                 ghsl_dir = ghsl_dir),
      error = function(e) {
        warning("[compute] ghsl_at_subblocks failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  # Results-aware: pull sub-blocks directly from an existing block in the
  # evidence registry. Useful when a candidate already has boundary /
  # regime_segments / recombinant_map blocks and you want to query GHSL
  # inside each declared sub-region without re-computing boundaries.
  ghsl_at_block_subregions <- function(cid, block_type,
                                          sample_subset = NULL,
                                          scale = NULL,
                                          summaries = c("mean", "frac_low",
                                                        "frac_mid", "frac_high",
                                                        "rank_mean", "n_windows"),
                                          ghsl_dir = NULL) {
    if (!.ensure_ghsl_panel_lib()) return(NULL)
    block <- tryCatch(evidence$read_block(cid, block_type),
                        error = function(e) NULL)
    if (is.null(block)) {
      warning("[compute] ghsl_at_block_subregions: no '", block_type,
              "' block for ", cid)
      return(NULL)
    }
    # Extract sub-block list from the block. Different block types store
    # sub-regions under different keys; try the common locations.
    subs <- NULL
    for (k in c("segments", "sub_blocks", "subblocks", "regions",
                 "tracts", "windows", "soft_boundaries_blocks")) {
      if (!is.null(block[[k]])) {
        subs <- block[[k]]
        break
      }
    }
    if (is.null(subs)) {
      warning("[compute] ghsl_at_block_subregions: could not find ",
              "sub-region list in '", block_type, "' block for ", cid)
      return(NULL)
    }
    subs_dt <- tryCatch(as.data.table(subs), error = function(e) NULL)
    if (is.null(subs_dt) || nrow(subs_dt) == 0) {
      warning("[compute] ghsl_at_block_subregions: empty sub-region list")
      return(NULL)
    }
    # Normalize expected column names
    if (!"subblock_id" %in% names(subs_dt)) {
      if ("id" %in% names(subs_dt)) setnames(subs_dt, "id", "subblock_id")
      else subs_dt[, subblock_id := seq_len(.N)]
    }
    if (!"start_bp" %in% names(subs_dt)) {
      for (alt in c("start", "bp_start", "begin_bp")) {
        if (alt %in% names(subs_dt)) { setnames(subs_dt, alt, "start_bp"); break }
      }
    }
    if (!"end_bp" %in% names(subs_dt)) {
      for (alt in c("end", "bp_end", "stop_bp")) {
        if (alt %in% names(subs_dt)) { setnames(subs_dt, alt, "end_bp"); break }
      }
    }
    if (!all(c("subblock_id", "start_bp", "end_bp") %in% names(subs_dt))) {
      warning("[compute] ghsl_at_block_subregions: block '", block_type,
              "' does not expose start_bp/end_bp fields we can read")
      return(NULL)
    }
    ghsl_at_subblocks(cid, subs_dt,
                        sample_subset = sample_subset, scale = scale,
                        summaries = summaries, ghsl_dir = ghsl_dir)
  }

  ghsl_wide_matrix <- function(cid_or_chrom,
                                 start_bp = NULL, end_bp = NULL,
                                 metric = "div_roll_s50",
                                 sample_order = NULL,
                                 sample_subset = NULL,
                                 ghsl_dir = NULL) {
    if (!.ensure_ghsl_panel_lib()) return(NULL)
    chrom <- cid_or_chrom
    if (is.null(start_bp) || is.null(end_bp)) {
      cand <- tryCatch(intervals$get_candidate(cid_or_chrom),
                         error = function(e) NULL)
      if (!is.null(cand) && length(cand) > 0) {
        chrom   <- if (is.list(cand) && !is.data.frame(cand)) cand$chrom
                   else if (is.data.frame(cand) && "chrom" %in% names(cand)) cand$chrom[1]
                   else cid_or_chrom
        start_bp <- if (is.null(start_bp)) {
          if (is.list(cand) && !is.data.frame(cand)) cand$start_bp %||% cand$start
          else cand$start_bp[1] %||% cand$start[1]
        } else start_bp
        end_bp <- if (is.null(end_bp)) {
          if (is.list(cand) && !is.data.frame(cand)) cand$end_bp %||% cand$end
          else cand$end_bp[1] %||% cand$end[1]
        } else end_bp
      }
    }
    if (is.null(start_bp)) start_bp <- 0
    if (is.null(end_bp))   end_bp   <- Inf
    tryCatch(
      ghsl_panel_wide_matrix(chrom, start_bp, end_bp,
                                sample_ids = sample_subset,
                                metric = metric,
                                sample_order = sample_order,
                                ghsl_dir = ghsl_dir),
      error = function(e) {
        warning("[compute] ghsl_wide_matrix failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  # ===========================================================================
  # Chat 15 (2026-04-17): Ancestry on-demand queries
  # ===========================================================================
  # Thin wrappers around unified_ancestry/wrappers/instant_q.R. Same pattern
  # as the GHSL block above: the wrapper is expected to have been sourced by
  # load_bridge.R (which configure_instant_q()s against the ancestry config);
  # if not, these methods return NULL with a warning.
  # ===========================================================================

  .ensure_instant_q <- function() {
    if (exists("get_Q_summary", mode = "function") &&
        exists("merge_local_Q_into_invlikeness", mode = "function")) return(TRUE)
    iqf <- Sys.getenv("INSTANT_Q_R", "")
    if (!nzchar(iqf) || !file.exists(iqf)) {
      for (cand in c(
        file.path(Sys.getenv("BASE", ""),
                   "unified_ancestry", "wrappers", "instant_q.R"),
        "unified_ancestry/wrappers/instant_q.R"
      )) {
        if (nzchar(cand) && file.exists(cand)) { iqf <- cand; break }
      }
    }
    if (!nzchar(iqf) || !file.exists(iqf)) {
      warning("[compute] unified_ancestry/wrappers/instant_q.R not found — ",
              "set INSTANT_Q_R env var or ensure load_bridge.R ran")
      return(FALSE)
    }
    anc_cfg <- Sys.getenv("ANCESTRY_CONFIG", "")
    tryCatch({
      source(iqf, local = FALSE)
      if (nzchar(anc_cfg) && file.exists(anc_cfg) &&
          exists("configure_instant_q", mode = "function")) {
        configure_instant_q(config_file = anc_cfg)
      }
    }, error = function(e) {
      warning("[compute] instant_q.R source failed: ", conditionMessage(e))
      return(FALSE)
    })
    exists("get_Q_summary", mode = "function")
  }

  ancestry_at_interval <- function(chrom, start_bp, end_bp, K = NULL,
                                      sample_set = NULL) {
    if (!.ensure_instant_q()) return(NULL)
    tryCatch({
      dt <- get_Q_summary(chrom, K = K)
      if (is.null(dt) || !nrow(dt)) return(dt)
      s_col <- intersect(c("start_bp", "start"), names(dt))[1]
      e_col <- intersect(c("end_bp",   "end"),   names(dt))[1]
      if (!is.na(s_col) && !is.na(e_col)) {
        dt <- dt[get(s_col) < end_bp & get(e_col) > start_bp]
      }
      dt
    }, error = function(e) {
      warning("[compute] ancestry_at_interval failed: ", conditionMessage(e))
      NULL
    })
  }

  ancestry_at_candidate <- function(cid, K = NULL, sample_set = NULL) {
    cand <- intervals$get_candidate(cid)
    if (is.null(cand)) {
      warning("[compute] ancestry_at_candidate: unknown cid ", cid)
      return(NULL)
    }
    .pick <- function(x, field) {
      if (is.null(x)) return(NULL)
      if (is.list(x) && !is.data.frame(x) && field %in% names(x)) return(x[[field]])
      if (is.data.frame(x) && field %in% names(x)) return(x[[field]][1])
      NULL
    }
    chrom <- .pick(cand, "chrom")
    s     <- .pick(cand, "start_bp") %||% .pick(cand, "start")
    e     <- .pick(cand, "end_bp")   %||% .pick(cand, "end")
    if (is.null(chrom) || is.null(s) || is.null(e)) {
      warning("[compute] ancestry_at_candidate: candidate ", cid,
              " missing chrom/start/end")
      return(NULL)
    }
    ancestry_at_interval(chrom, as.integer(s), as.integer(e), K = K,
                          sample_set = sample_set)
  }

  ancestry_q_vector <- function(chrom, start_bp, end_bp, K = NULL,
                                  sample_ids = NULL) {
    if (!.ensure_instant_q()) return(NULL)
    tryCatch({
      sf <- resolve_cache_samples(chrom, K)
      if (is.na(sf)) {
        warning("[compute] ancestry_q_vector: no samples cache for ", chrom,
                " K=", K %||% "<canon>")
        return(NULL)
      }
      dt <- data.table::fread(sf)
      if (!nrow(dt)) return(dt)
      s_col <- intersect(c("start_bp", "start"), names(dt))[1]
      e_col <- intersect(c("end_bp",   "end"),   names(dt))[1]
      if (!is.na(s_col) && !is.na(e_col)) {
        dt <- dt[get(s_col) < end_bp & get(e_col) > start_bp]
      }
      if (!is.null(sample_ids) && "sample_id" %in% names(dt)) {
        dt <- dt[sample_id %in% sample_ids]
      }
      dt
    }, error = function(e) {
      warning("[compute] ancestry_q_vector failed: ", conditionMessage(e))
      NULL
    })
  }

  ancestry_q_summary <- function(chrom, K = NULL) {
    if (!.ensure_instant_q()) return(NULL)
    tryCatch(get_Q_summary(chrom, K = K),
             error = function(e) {
               warning("[compute] ancestry_q_summary failed: ",
                       conditionMessage(e))
               NULL
             })
  }

  # ── ancestry_q_and_f_for_candidate ────────────────────────────────────────
  # Returns Q (per-sample × per-window Q matrix for the candidate interval)
  # and F (global allele-freq matrix at requested K, read from BEST_FOPT).
  # If persist=TRUE and reg$stats is wired, writes both into
  # stats_cache/candidate/<cid>/. This is the "per-candidate archive" that
  # the manuscript Q/F/tree section needs.
  ancestry_q_and_f_for_candidate <- function(cid, K = NULL, persist = TRUE) {
    cand <- intervals$get_candidate(cid)
    if (is.null(cand)) {
      warning("[compute] ancestry_q_and_f_for_candidate: unknown cid ", cid)
      return(list(Q = NULL, F = NULL))
    }
    bounds <- intervals$get_candidate_boundaries(cid)
    chrom  <- if (is.list(cand) && !is.data.frame(cand)) cand$chrom
              else cand$chrom[1]
    q_dt <- ancestry_q_vector(chrom = chrom,
                                start_bp = bounds$left_bp,
                                end_bp   = bounds$right_bp,
                                K = K)
    # F matrix from the global BEST_FOPT (per-K via filename swap)
    f_mat <- tryCatch({
      canon_K <- if (exists(".iq_env") && !is.null(.iq_env$canonical_K))
                     .iq_env$canonical_K else 8L
      K_eff <- as.integer(K %||% canon_K)
      fopt_path <- if (exists(".iq_env")) .iq_env$fopt_file else NULL
      if (!is.null(fopt_path) && nzchar(fopt_path) && K_eff != canon_K) {
        fopt_path <- sub(sprintf("K%02d", canon_K),
                          sprintf("K%02d", K_eff), fopt_path, fixed = TRUE)
      }
      if (!is.null(fopt_path) && file.exists(fopt_path)) {
        as.matrix(data.table::fread(fopt_path))
      } else {
        NULL
      }
    }, error = function(e) NULL)
    out <- list(Q = q_dt, F = f_mat, cid = cid, K = K, chrom = chrom,
                 start_bp = bounds$left_bp, end_bp = bounds$right_bp)
    # Persist if caller asked AND reg$stats API is in scope (set by
    # load_registry() after both compute and stats APIs are built).
    if (isTRUE(persist)) {
      stats_api <- tryCatch(get(".REG_STATS_API", envir = .GlobalEnv,
                                   inherits = FALSE),
                              error = function(e) NULL)
      if (!is.null(stats_api) &&
          !is.null(stats_api$put_candidate_q_f)) {
        tryCatch(stats_api$put_candidate_q_f(cid, q_dt, f_mat, K = K),
                  error = function(e)
                    warning("[compute] persist Q/F for ", cid, " failed: ",
                            conditionMessage(e)))
      }
    }
    out
  }

  # ══ Chat 16 — multi-window scan drivers ══════════════════════════════
  # Run a per-window statistic across any interval with automatic
  # per-window group resolution via intervals$resolve_smallest_at.
  # Handles the segmented-karyotype case: when a window falls inside
  # a recombinant segment (seg_DCO, seg_L, seg_R), the scan uses the
  # segment-level karyotype groups; when it falls in non-segmented
  # parent-scale zones, it uses the parent's karyotype groups. Fully
  # automatic, no hard-coding per chromosome.
  #
  # Empty-segment fallback: if a segment's karyotype group has fewer
  # than min_group_n samples, fall back to the parent candidate's group.
  # This protects against noisy FST estimates in tiny DCO tracts where
  # only one or two carriers are HOM_REF.

  scan_pairwise <- function(chrom, start_bp, end_bp,
                               window_size = 10000, step = NULL,
                               stat = "fst",
                               karyo_pair = c("HOM_REF", "HOM_INV"),
                               min_group_n = 3,
                               persist = TRUE,
                               source_script = "scan_pairwise") {
    step <- step %||% window_size
    if (length(karyo_pair) != 2L)
      stop("[scan_pairwise] karyo_pair must be length 2")
    k1 <- karyo_pair[1]; k2 <- karyo_pair[2]

    # Build window centers
    centers <- seq(as.integer(start_bp), as.integer(end_bp),
                    by = as.integer(step))
    if (!length(centers)) return(data.table())

    # Per window: resolve smallest candidate, look up both groups at
    # that cid, fallback to parent if group too small, compute stat.
    out <- vector("list", length(centers))
    for (i in seq_along(centers)) {
      pos <- centers[i]
      cid_local <- intervals$resolve_smallest_at(chrom, pos)
      if (is.na(cid_local)) {
        # No candidate covers this position — skip with NA
        out[[i]] <- data.table(
          window_mid = pos, candidate_id = NA_character_,
          group_1 = NA_character_, group_2 = NA_character_,
          karyotype_1 = k1, karyotype_2 = k2,
          n_1 = NA_integer_, n_2 = NA_integer_,
          value = NA_real_, used_fallback = FALSE)
        next
      }
      g1_name <- sprintf("inv_%s_%s", cid_local, k1)
      g2_name <- sprintf("inv_%s_%s", cid_local, k2)
      g1 <- if (samples$has_group(g1_name)) samples$get_group(g1_name)
            else character(0)
      g2 <- if (samples$has_group(g2_name)) samples$get_group(g2_name)
            else character(0)
      used_fallback <- FALSE
      # Fallback to parent candidate if either segment group too small
      if ((length(g1) < min_group_n || length(g2) < min_group_n)) {
        cand_here <- tryCatch(intervals$get_candidate(cid_local),
                               error = function(e) NULL)
        parent_id <- if (is.null(cand_here)) NA_character_
                     else if (is.list(cand_here) && !is.data.frame(cand_here))
                              cand_here$parent_id %||% NA_character_
                          else cand_here$parent_id[1] %||% NA_character_
        if (!is.na(parent_id) && nzchar(parent_id) && parent_id != ".") {
          g1_pname <- sprintf("inv_%s_%s", parent_id, k1)
          g2_pname <- sprintf("inv_%s_%s", parent_id, k2)
          if (samples$has_group(g1_pname) && samples$has_group(g2_pname)) {
            g1 <- samples$get_group(g1_pname)
            g2 <- samples$get_group(g2_pname)
            g1_name <- g1_pname; g2_name <- g2_pname
            cid_local <- parent_id
            used_fallback <- TRUE
          }
        }
      }
      # Still too small after fallback — record NA
      if (length(g1) < min_group_n || length(g2) < min_group_n) {
        out[[i]] <- data.table(
          window_mid = pos, candidate_id = cid_local,
          group_1 = g1_name, group_2 = g2_name,
          karyotype_1 = k1, karyotype_2 = k2,
          n_1 = length(g1), n_2 = length(g2),
          value = NA_real_, used_fallback = used_fallback)
        next
      }
      # Compute
      r <- tryCatch(
        pairwise_stat(pos = pos, flank_size = as.integer(window_size / 2),
                       group1 = g1, group2 = g2,
                       stat = stat, chrom = chrom, cache = FALSE),
        error = function(e) NULL)
      val <- if (is.null(r)) NA_real_ else {
        if (is.list(r) && !is.null(r[[stat]])) as.numeric(r[[stat]])
        else if (is.list(r) && !is.null(r$value)) as.numeric(r$value)
        else NA_real_
      }
      out[[i]] <- data.table(
        window_mid = pos, candidate_id = cid_local,
        group_1 = g1_name, group_2 = g2_name,
        karyotype_1 = k1, karyotype_2 = k2,
        n_1 = length(g1), n_2 = length(g2),
        value = val, used_fallback = used_fallback)
    }
    result_dt <- rbindlist(out)

    # Persist ONE row per unique (candidate_id × stat) seen in the scan.
    # Individual window values live inside the tsv.gz payload. This is
    # cleaner than writing 10k manifest rows for a 100 Mb scan.
    if (isTRUE(persist)) {
      stats_api <- tryCatch(get(".REG_RESULTS_API", envir = .GlobalEnv,
                                  inherits = FALSE),
                              error = function(e) NULL)
      if (!is.null(stats_api)) {
        for (cid_uq in unique(result_dt$candidate_id)) {
          if (is.na(cid_uq)) next
          rd <- result_dt[candidate_id == cid_uq]
          tryCatch({
            # Use put_pairwise for per-cid scan persistence. Group labels
            # at segment scale; fallback-flagged rows live in same file.
            stats_api$put_pairwise(
              chrom = chrom,
              group_1 = sprintf("inv_%s_%s", cid_uq, k1),
              group_2 = sprintf("inv_%s_%s", cid_uq, k2),
              stat = stat, dt = rd,
              source_script = source_script,
              engine = "scan_pairwise")
          }, error = function(e)
               warning("[scan_pairwise] persist failed for ", cid_uq,
                       ": ", conditionMessage(e)))
        }
      }
    }
    result_dt
  }

  # Single-group region statistic scan — ENA, entropy, Δ12, expected
  # heterozygosity within one karyotype group. Less machinery than
  # scan_pairwise because only ONE group per window is needed.
  scan_region_stat <- function(chrom, start_bp, end_bp,
                                  window_size = 10000, step = NULL,
                                  stat = "delta12",
                                  karyotype = "HOM_INV",
                                  min_group_n = 3,
                                  persist = TRUE,
                                  source_script = "scan_region_stat") {
    step <- step %||% window_size
    if (!(stat %in% c("delta12", "entropy", "ena", "ancestry_q_mean")))
      stop("[scan_region_stat] stat '", stat, "' not in supported set")

    centers <- seq(as.integer(start_bp), as.integer(end_bp),
                    by = as.integer(step))
    if (!length(centers)) return(data.table())

    out <- vector("list", length(centers))
    for (i in seq_along(centers)) {
      pos <- centers[i]
      cid_local <- intervals$resolve_smallest_at(chrom, pos)
      if (is.na(cid_local)) {
        out[[i]] <- data.table(
          window_mid = pos, candidate_id = NA_character_,
          group_id = NA_character_, karyotype = karyotype,
          n = NA_integer_, value = NA_real_, used_fallback = FALSE)
        next
      }
      g_name <- sprintf("inv_%s_%s", cid_local, karyotype)
      g <- if (samples$has_group(g_name)) samples$get_group(g_name)
           else character(0)
      used_fallback <- FALSE
      if (length(g) < min_group_n) {
        # Try parent
        cand <- tryCatch(intervals$get_candidate(cid_local),
                         error = function(e) NULL)
        p_id <- if (is.null(cand)) NA_character_
                else if (is.list(cand) && !is.data.frame(cand))
                         cand$parent_id %||% NA_character_
                     else cand$parent_id[1] %||% NA_character_
        if (!is.na(p_id) && nzchar(p_id) && p_id != ".") {
          g_pname <- sprintf("inv_%s_%s", p_id, karyotype)
          if (samples$has_group(g_pname)) {
            g <- samples$get_group(g_pname)
            g_name <- g_pname; cid_local <- p_id
            used_fallback <- TRUE
          }
        }
      }
      if (length(g) < min_group_n) {
        out[[i]] <- data.table(
          window_mid = pos, candidate_id = cid_local,
          group_id = g_name, karyotype = karyotype,
          n = length(g), value = NA_real_, used_fallback = used_fallback)
        next
      }
      # Compute per-window stat. For ancestry-derived stats (delta12,
      # entropy, ena) use ancestry_at_interval; for custom stats extend
      # with stat-specific dispatchers later.
      val <- tryCatch({
        if (stat %in% c("delta12", "entropy", "ena")) {
          dt <- ancestry_at_interval(chrom, pos - window_size %/% 2L,
                                      pos + window_size %/% 2L)
          if (is.null(dt) || !nrow(dt)) NA_real_
          else {
            # Mean of the stat column across windows in this interval,
            # restricted to samples in g
            col_map <- list(delta12 = "mean_delta12",
                            entropy = "mean_entropy",
                            ena     = "mean_ena")
            col <- col_map[[stat]]
            if (col %in% names(dt)) mean(dt[[col]], na.rm = TRUE)
            else NA_real_
          }
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_)

      out[[i]] <- data.table(
        window_mid = pos, candidate_id = cid_local,
        group_id = g_name, karyotype = karyotype,
        n = length(g), value = val, used_fallback = used_fallback)
    }
    result_dt <- rbindlist(out)

    if (isTRUE(persist)) {
      stats_api <- tryCatch(get(".REG_RESULTS_API", envir = .GlobalEnv,
                                  inherits = FALSE),
                              error = function(e) NULL)
      if (!is.null(stats_api)) {
        for (cid_uq in unique(result_dt$candidate_id)) {
          if (is.na(cid_uq)) next
          rd <- result_dt[candidate_id == cid_uq]
          # interval_summary persist — one row per (cid × stat × karyotype)
          cand <- tryCatch(intervals$get_candidate(cid_uq),
                           error = function(e) NULL)
          if (is.null(cand)) next
          s_bp <- if (is.list(cand) && !is.data.frame(cand))
                     cand$start_bp else cand$start_bp[1]
          e_bp <- if (is.list(cand) && !is.data.frame(cand))
                     cand$end_bp else cand$end_bp[1]
          tryCatch({
            stats_api$put_interval_summary(
              chrom = chrom,
              start_bp = as.integer(s_bp),
              end_bp   = as.integer(e_bp),
              group_1  = sprintf("inv_%s_%s", cid_uq, karyotype),
              stat     = stat, dt = rd,
              source_script = source_script,
              engine = "scan_region_stat")
          }, error = function(e)
               warning("[scan_region_stat] persist failed for ", cid_uq,
                       ": ", conditionMessage(e)))
        }
      }
    }
    result_dt
  }

  # scan_candidate_full — the "run everything automatic" driver for one
  # candidate. Given a cid, scans:
  #   - FST (HOM_REF vs HOM_INV), (HET vs HOM_INV), (HET vs HOM_REF)
  #   - Δ12 / entropy / ENA within each karyotype
  # across the parent interval + flanking region. Persists all results
  # to the manifest. Returns a list of the per-comparison data.tables.
  #
  # This is the one call a classification script needs per candidate:
  # "give me the full panel of population statistics for LG12_17".
  scan_candidate_full <- function(cid, flank_kb = 100,
                                     window_size = 10000, step = NULL,
                                     min_group_n = 3,
                                     karyotypes = c("HOM_REF", "HET", "HOM_INV"),
                                     pairwise_stats = c("fst"),
                                     region_stats   = c("delta12", "entropy", "ena"),
                                     persist = TRUE) {
    cand <- intervals$get_candidate(cid)
    if (is.null(cand)) {
      warning("[scan_candidate_full] unknown cid ", cid); return(NULL)
    }
    chrom <- if (is.list(cand) && !is.data.frame(cand)) cand$chrom
             else cand$chrom[1]
    s_bp  <- as.integer(cand$start_bp) - as.integer(flank_kb) * 1000L
    e_bp  <- as.integer(cand$end_bp)   + as.integer(flank_kb) * 1000L
    out <- list()

    # Pairwise scans
    k_pairs <- list()
    ks <- karyotypes
    for (i in seq_along(ks)) {
      for (j in seq_along(ks)) {
        if (i < j) k_pairs[[length(k_pairs) + 1L]] <- c(ks[i], ks[j])
      }
    }
    for (pair in k_pairs) {
      for (st in pairwise_stats) {
        key <- sprintf("pairwise_%s_%s_vs_%s", st, pair[1], pair[2])
        out[[key]] <- tryCatch(
          scan_pairwise(chrom, s_bp, e_bp,
                         window_size = window_size, step = step,
                         stat = st, karyo_pair = pair,
                         min_group_n = min_group_n,
                         persist = persist,
                         source_script = sprintf("scan_candidate_full(%s)", cid)),
          error = function(e) {
            warning("[scan_candidate_full] ", key, " for ", cid, ": ",
                    conditionMessage(e)); NULL
          })
      }
    }

    # Per-karyotype region-stat scans
    for (k in karyotypes) {
      for (st in region_stats) {
        key <- sprintf("region_%s_%s", st, k)
        out[[key]] <- tryCatch(
          scan_region_stat(chrom, s_bp, e_bp,
                            window_size = window_size, step = step,
                            stat = st, karyotype = k,
                            min_group_n = min_group_n,
                            persist = persist,
                            source_script = sprintf("scan_candidate_full(%s)", cid)),
          error = function(e) {
            warning("[scan_candidate_full] ", key, " for ", cid, ": ",
                    conditionMessage(e)); NULL
          })
      }
    }
    out
  }

  list(
    pairwise_stat          = pairwise_stat,
    boundary_fst_profile   = function(cid, side, distances_kb = c(10, 25, 50, 100)) {
      bounds <- intervals$get_candidate_boundaries(cid)
      if (is.null(bounds)) return(NULL)
      anchor <- if (side == "left") bounds$left_bp else bounds$right_bp
      chrom  <- intervals$get_candidate(cid)$chrom
      g1 <- samples$get_carriers(cid, "HOM_REF")
      g2 <- samples$get_carriers(cid, "HOM_INV")
      out <- data.table(distance_kb = distances_kb, fst = NA_real_)
      for (i in seq_along(distances_kb)) {
        pos <- anchor + distances_kb[i] * 1000 * (if (side == "right") 1 else -1)
        r <- pairwise_stat(pos, 10000, g1, g2, "fst", chrom = chrom)
        if (!is.null(r)) out$fst[i] <- r$fst %||% NA
      }
      out
    },
    # Chat 14: GHSL v6 on-demand query API
    ghsl_at_interval          = ghsl_at_interval,
    ghsl_at_candidate         = ghsl_at_candidate,
    ghsl_at_subblocks         = ghsl_at_subblocks,
    ghsl_at_block_subregions  = ghsl_at_block_subregions,
    ghsl_wide_matrix          = ghsl_wide_matrix,
    # Chat 15: Ancestry on-demand query API
    ancestry_at_interval          = ancestry_at_interval,
    ancestry_at_candidate         = ancestry_at_candidate,
    ancestry_q_vector             = ancestry_q_vector,
    ancestry_q_summary            = ancestry_q_summary,
    ancestry_q_and_f_for_candidate = ancestry_q_and_f_for_candidate,
    # Chat 16: multi-window scan drivers (auto per-window group resolution)
    scan_pairwise                 = scan_pairwise,
    scan_region_stat              = scan_region_stat,
    scan_candidate_full           = scan_candidate_full
  )
}

message("[registry_loader] module loaded — call load_registry() to initialize")
