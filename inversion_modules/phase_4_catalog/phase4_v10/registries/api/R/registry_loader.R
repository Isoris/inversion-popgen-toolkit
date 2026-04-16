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
# Samples API — thin wrapper around existing sample_registry.R
# =============================================================================
load_samples_api <- function(sample_dir) {
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
        return(list(
          get_master   = full$list_groups,
          get_group    = full$get_group,
          add_group    = full$add_group,
          has_group    = full$has_group,
          list_groups  = full$list_groups,
          count_groups = function() nrow(full$list_groups()),
          # Query convenience: get samples by carrier status for a candidate
          get_carriers = function(cid, status = "HOM_INV") {
            gid <- paste0("inv_", cid, "_", status)
            if (full$has_group(gid)) full$get_group(gid) else character(0)
          }
        ))
      }
    }
  }

  # Shim: no-op implementations with warnings
  warning("[registry_loader] sample_registry.R not found — using shim")
  list(
    get_master   = function() data.table(),
    get_group    = function(gid) character(0),
    add_group    = function(...) invisible(FALSE),
    has_group    = function(gid) FALSE,
    list_groups  = function() data.table(),
    count_groups = function() 0L,
    get_carriers = function(cid, status) character(0)
  )
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

  list(
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
      c <- fread(cand_file)
      r <- c[c$candidate_id == cid]
      if (nrow(r) == 0) NULL else as.list(r[1])
    },
    get_candidate_boundaries = function(cid) {
      r <- fread(cand_file)
      r <- r[r$candidate_id == cid]
      if (nrow(r) == 0) return(NULL)
      list(left_bp = as.integer(r$start_bp[1]), right_bp = as.integer(r$end_bp[1]))
    },
    add_candidate = function(cid, chrom, start_bp, end_bp, scale = "inversion_block",
                               parent_id = NA_character_) {
      c <- fread(cand_file)
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
    n_extracted <- 0L
    if (!is.null(schema) && !is.null(schema$keys_extracted)) {
      kex <- schema$keys_extracted
      flat <- data.table()
      for (ke in kex) {
        # Each entry: list(key = "q7_layer_d_fisher_or", from = "fisher_or")
        key <- ke$key; from <- ke$from
        if (is.null(key) || is.null(from)) next
        val <- data[[from]] %||% NA
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
