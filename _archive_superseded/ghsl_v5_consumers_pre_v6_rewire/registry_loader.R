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

  # Return composite
  list(
    root      = registries_root,
    samples   = samples,
    intervals = intervals,
    evidence  = evidence,
    query     = query,
    compute   = compute,

    # Shortcuts for the most common v9-style calls (forward-compat during migration)
    add_group         = samples$add_group,
    get_group         = samples$get_group,
    has_group         = samples$has_group,
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
load_compute_api <- function(samples, intervals, evidence) {
  # Requires get_region_stats() / resolve_groups() from load_bridge.R
  # If Engine B not available, methods return NULL with a warning.

  pairwise_stat <- function(position, flank_size, group1, group2, stat = "fst",
                               chrom = NULL) {
    if (!exists("get_region_stats", mode = "function")) {
      warning("[compute] Engine B not loaded — pairwise_stat NULL")
      return(NULL)
    }
    if (is.null(chrom)) {
      warning("[compute] chrom required")
      return(NULL)
    }
    tryCatch({
      res <- get_region_stats(chrom,
                                start_bp = position - flank_size,
                                end_bp   = position + flank_size,
                                what     = stat,
                                groups   = list(g1 = group1, g2 = group2))
      res
    }, error = function(e) {
      warning("[compute] pairwise_stat failed: ", conditionMessage(e))
      NULL
    })
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
    }
  )
}

message("[registry_loader] module loaded — call load_registry() to initialize")
