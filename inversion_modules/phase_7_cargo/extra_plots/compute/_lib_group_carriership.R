#!/usr/bin/env Rscript

# =============================================================================
# _lib_group_carriership.R
#
# Helper sourced by PLOT_06/07/08. Returns a long-format data.table:
#   var_id | class | group | n_group | n_carriers_in_group | freq_in_group
#
# Resolves group membership from --groups-from <tsv> if provided, otherwise
# from NGSADMIX_Q_FILE argmax cluster.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

resolve_groups <- function(groups_from = NULL, sample_list = NULL,
                            ngsadmix_q = NULL, canonical_k = 8) {
  if (!is.null(groups_from) && file.exists(groups_from)) {
    gd <- fread(groups_from)
    return(setNames(gd$group_id, gd$sample))
  }
  if (!is.null(ngsadmix_q) && file.exists(ngsadmix_q) &&
      !is.null(sample_list) && file.exists(sample_list)) {
    Q <- as.matrix(fread(ngsadmix_q, header = FALSE))
    snames <- as.character(fread(sample_list, header = FALSE)[[1]])
    n <- min(length(snames), nrow(Q))
    am <- apply(Q[seq_len(n), , drop = FALSE], 1, which.max)
    return(setNames(paste0("K", canonical_k, "_Q", am), snames[seq_len(n)]))
  }
  NULL
}

classify_variants <- function(v) {
  cls <- rep("OTHER", nrow(v))
  ref_len <- if ("ref" %in% names(v)) nchar(v$ref) else NA
  alt_len <- if ("alt" %in% names(v)) nchar(v$alt) else NA
  if (!is.na(ref_len[1])) {
    cls[ref_len == 1 & alt_len == 1] <- "SNP"
    cls[ref_len != alt_len & pmax(ref_len, alt_len) <= 50] <- "INDEL"
  }
  if ("sv_type" %in% names(v)) {
    sv <- toupper(v$sv_type)
    for (k in c("DEL","DUP","INV","BND","INS")) cls[sv == k] <- k
  }
  if ("variant_class" %in% names(v))
    cls[toupper(v$variant_class) == "PAV"] <- "PAV"
  cls
}

build_carriership_long <- function(variant_master_path, sample_to_group) {
  v <- fread(variant_master_path)
  v[, class := classify_variants(.SD)]
  v[, var_id := if ("var_key" %in% names(.SD)) var_key else paste0(
      if ("chr" %in% names(.SD)) chr else chrom, ":", pos)]

  if (!"samples_with_alt" %in% names(v)) {
    # Fall back to per-sample columns
    sample_cols <- intersect(names(sample_to_group), names(v))
    if (length(sample_cols) < 5) return(NULL)
    # Vectorised melt in chunks (saves time vs per-sample subset loop)
    chunk_n <- 30L
    rows <- list()
    for (chunk_start in seq(1, length(sample_cols), by = chunk_n)) {
      chunk_end <- min(chunk_start + chunk_n - 1L, length(sample_cols))
      cols_chunk <- sample_cols[chunk_start:chunk_end]
      mlt <- melt(v[, c("var_id", "class", cols_chunk), with = FALSE],
                   id.vars = c("var_id", "class"),
                   variable.name = "sample",
                   value.name = "dosage")
      rows[[length(rows) + 1]] <- mlt[dosage > 0,
                                        .(var_id, class, sample = as.character(sample))]
    }
    long <- rbindlist(rows)
  } else {
    long <- v[, .(sample = unlist(strsplit(samples_with_alt, ";"))),
              by = .(var_id, class)]
    long <- long[nzchar(sample)]
  }
  long[, group := sample_to_group[sample]]
  long <- long[!is.na(group)]
  groups_avail <- sort(unique(unname(sample_to_group)))
  group_size_dt <- data.table(group = groups_avail,
                                n_group = as.integer(table(sample_to_group)[groups_avail]))
  out <- long[, .(n_carriers_in_group = uniqueN(sample)),
              by = .(var_id, class, group)]
  out <- merge(out, group_size_dt, by = "group", all.x = TRUE)
  out[, freq_in_group := n_carriers_in_group / n_group]
  out[]
}
