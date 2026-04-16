#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(paste(
    "Usage: Rscript STEP10b_parse_pestPG_to_sample_theta_windows.R",
    "<pestPG_dir> <samples_file> <outprefix>",
    "[scale_pattern=win50000.step10000] [chr_filter=]"
  ))
}

pestPG_dir    <- args[1]
samples_file  <- args[2]
outprefix     <- args[3]
scale_pattern <- if (length(args) >= 4) args[4] else "win50000.step10000"
chr_filter    <- if (length(args) >= 5 && nchar(args[5]) > 0) args[5] else ""

if (!file.exists(samples_file)) stop("Samples file not found: ", samples_file)

samp_raw <- fread(samples_file, header = FALSE, sep = "\t", fill = TRUE)
sample_names <- samp_raw[[1]]
sample_names <- sample_names[nchar(sample_names) > 0]

message("[INFO] Samples loaded: ", length(sample_names))

all_files <- list.files(pestPG_dir, pattern = "\\.pestPG$", full.names = TRUE, recursive = TRUE)
if (length(all_files) == 0) stop("No .pestPG files found in: ", pestPG_dir)

if (nchar(scale_pattern) > 0) {
  all_files <- all_files[grepl(scale_pattern, basename(all_files), fixed = TRUE)]
}
if (length(all_files) == 0) stop("No .pestPG files matching pattern '", scale_pattern, "' found")

message("[INFO] Found ", length(all_files), " .pestPG files matching '", scale_pattern, "'")

parse_one_pestPG <- function(filepath, sample_id, chr_filter = "") {
  lines <- readLines(filepath, warn = FALSE)
  if (length(lines) < 2) return(NULL)

  hdr_idx <- grep("^#", lines)
  if (length(hdr_idx) == 0) return(NULL)

  data_lines <- lines[(hdr_idx[1] + 1):length(lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

  out <- vector("list", length(data_lines))
  keep_n <- 0L

  for (line in data_lines) {
    paren_end <- regexpr("\\)\t", line)
    if (paren_end < 0) next

    paren_part <- substr(line, 1, paren_end)
    rest_part  <- substr(line, paren_end + 2, nchar(line))

    rest_fields <- strsplit(rest_part, "\t", fixed = TRUE)[[1]]
    if (length(rest_fields) < 5) next

    chrom_now <- rest_fields[1]

    # IMPORTANT: skip immediately if not the chromosome we want
    if (nchar(chr_filter) > 0 && chrom_now != chr_filter) next

    paren_matches <- gregexpr("\\(([^)]+)\\)", paren_part)
    paren_contents <- regmatches(paren_part, paren_matches)[[1]]
    if (length(paren_contents) < 3) next

    win_pair <- gsub("[()]", "", paren_contents[3])
    win_parts <- strsplit(win_pair, ",", fixed = TRUE)[[1]]
    if (length(win_parts) != 2) next

    keep_n <- keep_n + 1L
    out[[keep_n]] <- data.table(
      chrom     = chrom_now,
      WinStart  = as.numeric(win_parts[1]),
      WinStop   = as.numeric(win_parts[2]),
      WinCenter = as.numeric(rest_fields[2]),
      tP        = as.numeric(rest_fields[4]),
      nSites    = as.integer(rest_fields[length(rest_fields)]),
      sample    = sample_id
    )
  }

  if (keep_n == 0L) return(NULL)

  res <- rbindlist(out[seq_len(keep_n)], fill = TRUE)
  res <- res[is.finite(tP)]
  if (nrow(res) == 0) return(NULL)

  res
}

all_results <- list()

for (sid in sample_names) {
  matching <- all_files[grepl(sid, basename(all_files), fixed = TRUE)]

  if (length(matching) == 0) {
    message("[WARN] No .pestPG file found for sample: ", sid)
    next
  }

  if (length(matching) > 1) {
    exact <- matching[grepl(paste0("^", sid, "\\."), basename(matching))]
    if (length(exact) >= 1) matching <- exact[1] else matching <- matching[1]
  } else {
    matching <- matching[1]
  }

  parsed <- parse_one_pestPG(matching, sid, chr_filter)
  if (!is.null(parsed) && nrow(parsed) > 0) {
    all_results[[length(all_results) + 1]] <- parsed
  }
}

if (length(all_results) == 0) stop("No .pestPG data successfully parsed")

merged <- rbindlist(all_results, fill = TRUE)

message("[INFO] Parsed ", nrow(merged), " window-sample records across ", uniqueN(merged$sample), " samples")
if (nchar(chr_filter) > 0) {
  message("[INFO] Chromosome filter: ", chr_filter)
}

tsv_out <- paste0(outprefix, ".sample_tP_windows.tsv.gz")
rds_out <- paste0(outprefix, ".sample_tP_windows.rds")

fwrite(merged, tsv_out, sep = "\t")
saveRDS(merged, rds_out)

message("[DONE] Wrote:")
message("  ", tsv_out, "  (", nrow(merged), " rows, ", uniqueN(merged$sample), " samples, ", uniqueN(merged$chrom), " chromosomes)")
message("  ", rds_out)
