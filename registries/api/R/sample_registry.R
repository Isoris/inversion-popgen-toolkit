#!/usr/bin/env Rscript

# =============================================================================
# sample_registry.R — Append-only sample group registry
#
# Manages two tables:
#   sample_master.tsv  — one row per fish, static
#   sample_groups.tsv  — one row per group, append-only (grows over time)
#
# SAFETY:
#   - Never deletes rows (append-only)
#   - Auto-backup before every write (.bak.YYYYMMDD_HHMMSS)
#   - Detects duplicate group_id (skips, doesn't overwrite)
#   - Validates members exist in sample_master
#
# Usage:
#   source("utils/sample_registry.R")
#   reg <- load_registry()                # loads or creates both tables
#   reg$add_group("unrelated_81", samples, description = "NAToRA pruned")
#   reg$get_group("unrelated_81")         # returns sample_id vector
#   reg$list_groups()                     # shows all registered groups
#   reg$has_group("unrelated_81")         # TRUE/FALSE
#   reg$get_master()                      # full sample_master table
#   reg$sample_id_from_ind("Ind0")        # → "CGA009"
#
# =============================================================================

suppressPackageStartupMessages(library(data.table))

load_registry <- function(registry_dir = NULL) {

  # ── Find registry directory ──
  if (is.null(registry_dir)) {
    registry_dir <- Sys.getenv("REGISTRY_DIR", "")
    if (!nzchar(registry_dir)) {
      base <- Sys.getenv("BASE", ".")
      registry_dir <- file.path(base, "sample_registry")
    }
  }
  dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(registry_dir, "groups"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(registry_dir, "backups"), recursive = TRUE, showWarnings = FALSE)

  master_file <- file.path(registry_dir, "sample_master.tsv")
  groups_file <- file.path(registry_dir, "sample_groups.tsv")

  # ── Safe backup before write ──
  safe_backup <- function(f) {
    if (file.exists(f)) {
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      bak <- file.path(registry_dir, "backups",
                        paste0(basename(f), ".bak.", ts))
      file.copy(f, bak)
    }
  }

  # ── Initialize master table if needed ──
  if (!file.exists(master_file)) {
    # Try to build from sample list
    source_utils <- file.path(dirname(registry_dir), "inversion_codebase_v8.5", "utils", "sample_map.R")
    if (!file.exists(source_utils)) {
      source_utils <- Sys.glob(file.path(dirname(registry_dir), "inversion_codebase_*/utils/sample_map.R"))[1]
    }

    sample_list <- NULL
    for (f in c(
      Sys.getenv("SAMPLES_IND", ""),
      file.path(Sys.getenv("BASE", "."), "het_roh/01_inputs_check/samples.ind"),
      "het_roh/01_inputs_check/samples.ind",
      file.path(Sys.getenv("BASE", "."), "popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv")
    )) {
      if (nzchar(f) && file.exists(f)) { sample_list <- f; break }
    }

    if (!is.null(sample_list)) {
      real <- trimws(readLines(sample_list))
      real <- real[nzchar(real)]
      n <- length(real)
      master <- data.table(
        sample_id = real,
        ind_id_226 = paste0("Ind", seq_len(n) - 1L),
        sample_index = seq_len(n) - 1L
      )
      fwrite(master, master_file, sep = "\t")
      message("[registry] Created sample_master.tsv: ", n, " samples")
    } else {
      master <- data.table(sample_id = character(), ind_id_226 = character(),
                            sample_index = integer())
      fwrite(master, master_file, sep = "\t")
      message("[registry] WARNING: no sample list found, empty master created")
    }
  }

  # ── Initialize groups table if needed ──
  if (!file.exists(groups_file)) {
    groups <- data.table(
      group_id = character(), n = integer(), chrom = character(),
      inv_id = character(), subgroup = character(),
      description = character(), members_file = character(),
      created = character()
    )
    fwrite(groups, groups_file, sep = "\t")
    message("[registry] Created sample_groups.tsv")
  }

  # ── Load ──
  master <- fread(master_file)
  groups <- fread(groups_file)

  message("[registry] Loaded: ", nrow(master), " samples, ", nrow(groups), " groups")

  # ── Public API ──

  get_master <- function() {
    fread(master_file)  # always read fresh
  }

  list_groups <- function() {
    fread(groups_file)
  }

  has_group <- function(gid) {
    g <- fread(groups_file)
    gid %in% g$group_id
  }

  get_group <- function(gid) {
    g <- fread(groups_file)
    row <- g[group_id == gid]
    if (nrow(row) == 0) {
      message("[registry] Group not found: ", gid)
      return(character())
    }
    mf <- row$members_file[1]
    if (file.exists(mf)) {
      trimws(readLines(mf))
    } else {
      # Try relative to registry
      mf2 <- file.path(registry_dir, "groups", basename(mf))
      if (file.exists(mf2)) trimws(readLines(mf2))
      else { message("[registry] Members file not found: ", mf); character() }
    }
  }

  add_group <- function(group_id, sample_ids, chrom = ".", inv_id = ".",
                         subgroup = ".", description = "", overwrite = FALSE) {

    # Check duplicate
    # FIX 2026-04-17 (Finding AK, chat 11): force `created` to character on
    # read. fread auto-detects the "%Y-%m-%d %H:%M:%S" format and coerces to
    # POSIXct, which then class-mismatches the character value from
    # format(Sys.time()) on rbind — fatal after the first add_group in a
    # session. Bulk-writer modules (chat 13 seeded regions) hit this hard.
    g <- fread(groups_file, colClasses = c(created = "character"))
    if (group_id %in% g$group_id) {
      if (!overwrite) {
        message("[registry] Group already exists: ", group_id, " (", 
                g[group_id == (group_id)]$n[1], " samples). Skipping.")
        return(invisible(FALSE))
      }
      message("[registry] Overwriting group: ", group_id)
      safe_backup(groups_file)
      g <- g[group_id != (group_id)]  # remove old entry
    } else {
      safe_backup(groups_file)
    }

    # Validate samples exist in master
    m <- fread(master_file)
    unknown <- setdiff(sample_ids, m$sample_id)
    if (length(unknown) > 0) {
      message("[registry] WARNING: ", length(unknown), " sample_ids not in master: ",
              paste(head(unknown, 5), collapse = ", "))
    }

    # Write members file
    members_file <- file.path(registry_dir, "groups",
                               paste0(group_id, ".txt"))
    writeLines(sample_ids, members_file)

    # Append to groups table
    new_row <- data.table(
      group_id = group_id, n = length(sample_ids),
      chrom = chrom, inv_id = inv_id, subgroup = subgroup,
      description = description, members_file = members_file,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    g <- rbind(g, new_row)
    fwrite(g, groups_file, sep = "\t")

    message("[registry] Added group: ", group_id, " (", length(sample_ids),
            " samples) -> ", members_file)
    invisible(TRUE)
  }

  # ── Convenience: Ind → real name ──
  sample_id_from_ind <- function(ind_id) {
    m <- fread(master_file)
    idx <- match(ind_id, m$ind_id_226)
    ifelse(is.na(idx), ind_id, m$sample_id[idx])
  }

  ind_from_sample_id <- function(sid) {
    m <- fread(master_file)
    idx <- match(sid, m$sample_id)
    ifelse(is.na(idx), sid, m$ind_id_226[idx])
  }

  # ── Add metadata column to master ──
  add_master_column <- function(col_name, values, by_sample_id = NULL) {
    safe_backup(master_file)
    m <- fread(master_file)
    if (!is.null(by_sample_id)) {
      # values is a named vector or data.table with sample_id + col
      if (is.data.table(values)) {
        m <- merge(m, values, by = "sample_id", all.x = TRUE)
      } else {
        m[, (col_name) := values[match(sample_id, by_sample_id)]]
      }
    } else {
      if (length(values) == nrow(m)) {
        m[, (col_name) := values]
      } else {
        message("[registry] Length mismatch: ", length(values), " vs ", nrow(m))
        return(invisible(FALSE))
      }
    }
    fwrite(m, master_file, sep = "\t")
    message("[registry] Added column '", col_name, "' to sample_master")
    invisible(TRUE)
  }

  # Chat-16: get all groups associated with a candidate, labeled by karyotype
  # and dimension. Helps classification scripts discover sub-clusters without
  # string-parsing group names themselves. See DATABASE_DESIGN.md §
  # "Sample group naming convention" for the suffix conventions.
  #
  # Returns a data.table: group_id / karyotype / dimension / parent_group / n
  # where dimension is derived from the group name (or the stored `dimension`
  # column if present).
  get_subgroups <- function(cid) {
    all_g <- fread(groups_file, colClasses = c(created = "character"))
    prefix <- sprintf("^inv_%s_", cid)
    hits <- all_g[grepl(prefix, all_g$group_id)]
    if (!nrow(hits)) return(data.table(group_id = character(0),
                                          karyotype = character(0),
                                          dimension = character(0),
                                          parent_group = character(0),
                                          n = integer(0)))
    # Parse karyotype + dimension from group_id. Conventions:
    #   inv_<cid>_<KARYO>                    → top-level karyotype
    #   inv_<cid>_<KARYO>__sub<N>            → karyotype_subcluster
    #   inv_<cid>_<KARYO>__ancestry<K>_<k>   → ancestry
    #   inv_<cid>_<KARYO>__ghsl_<band>       → ghsl
    #   inv_<cid>_<KARYO>__family_<fid>      → family
    body <- sub(sprintf("^inv_%s_", cid), "", hits$group_id)
    has_suffix <- grepl("__", body)
    karyo <- ifelse(has_suffix, sub("__.*", "", body), body)
    suffix <- ifelse(has_suffix, sub(".*?__", "", body), "")
    dim <- rep(".", nrow(hits))
    dim[!has_suffix] <- "karyotype"
    dim[has_suffix & grepl("^sub[0-9]+$",   suffix)] <- "karyotype_subcluster"
    dim[has_suffix & grepl("^noise$",       suffix)] <- "karyotype_subcluster_noise"
    dim[has_suffix & grepl("^ancestry",     suffix)] <- "ancestry"
    dim[has_suffix & grepl("^ghsl_",        suffix)] <- "ghsl"
    dim[has_suffix & grepl("^family_",      suffix)] <- "family"
    # Override with stored dimension column if the groups_file has one
    if ("dimension" %in% names(hits)) {
      stored <- as.character(hits$dimension)
      dim <- ifelse(!is.na(stored) & nzchar(stored) & stored != ".", stored, dim)
    }
    parent <- ifelse(has_suffix, sprintf("inv_%s_%s", cid, karyo), NA_character_)
    data.table(
      group_id     = hits$group_id,
      karyotype    = karyo,
      dimension    = dim,
      parent_group = parent,
      n            = hits$n
    )
  }

  list(
    dir = registry_dir,
    master_file = master_file,
    groups_file = groups_file,
    get_master = get_master,
    list_groups = list_groups,
    has_group = has_group,
    get_group = get_group,
    add_group = add_group,
    get_subgroups = get_subgroups,
    sample_id_from_ind = sample_id_from_ind,
    ind_from_sample_id = ind_from_sample_id,
    add_master_column = add_master_column
  )
}

# =============================================================================
# Standalone: initialize registry
# =============================================================================
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  reg_dir <- if (length(args) >= 1) args[1] else NULL
  reg <- load_registry(reg_dir)
  cat("\nSample master:", reg$master_file, "\n")
  cat("Groups file:", reg$groups_file, "\n")
  cat("Samples:", nrow(reg$get_master()), "\n")
  cat("Groups:", nrow(reg$list_groups()), "\n")
}
