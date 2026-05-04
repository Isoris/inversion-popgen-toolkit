#!/usr/bin/env Rscript
# =============================================================================
# q05_aggregate_theta.R
#
# Parse all per-sample ANGSD thetaStat .pestPG files in a directory (filtered
# by scale pattern) and write one long-format TSV for a given chromosome.
#
# pestPG format (14 cols):
#   #(indexStart,indexStop)(firstPos,lastPos)(WinStart,WinStop) Chr WinCenter
#     tW tP tF tH tL Tajima fuf fud fayh zeng nSites
#
# Sample name is taken from the file basename:
#   <SAMPLE>.win<W>.step<S>.pestPG     -> SAMPLE
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
THETA_DIR <- get_arg("--theta_dir")
SCALE     <- get_arg("--scale", "win500000.step500000")
CHROM     <- get_arg("--chrom")
OUT       <- get_arg("--out")
MIN_SITES <- as.integer(get_arg("--min_sites", "100"))
stopifnot(!is.null(THETA_DIR), !is.null(CHROM), !is.null(OUT))

PESTPG_COLS <- c("WinInfo", "Chr", "WinCenter",
                 "tW", "tP", "tF", "tH", "tL",
                 "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

pattern <- paste0("\\.", SCALE, "\\.pestPG$")
files <- list.files(THETA_DIR, pattern = pattern, full.names = TRUE)
message("[q05] ", length(files), " pestPG files at scale ", SCALE)
if (length(files) == 0) stop("No files match pattern ", pattern, " in ", THETA_DIR)

# Extract sample ID from filename
rx <- paste0("^(.+)\\.", SCALE, "\\.pestPG$")
sample_ids <- sub(rx, "\\1", basename(files))

parse_one <- function(f, sample) {
  # Read header first to validate column count
  hdr_line <- readLines(f, n = 1)
  hdr <- strsplit(sub("^#", "", hdr_line), "\t", fixed = TRUE)[[1]]
  if (length(hdr) != length(PESTPG_COLS)) {
    warning(sample, ": header has ", length(hdr), " cols, expected ",
            length(PESTPG_COLS), " — skipping")
    return(NULL)
  }
  dt <- tryCatch(
    fread(f, sep = "\t", skip = 1, header = FALSE, col.names = PESTPG_COLS),
    error = function(e) { warning(sample, ": ", conditionMessage(e)); NULL })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt <- dt[Chr == CHROM & nSites >= MIN_SITES]
  if (nrow(dt) == 0) return(NULL)
  data.table(
    sample             = sample,
    chrom              = dt$Chr,
    win_center_bp      = as.integer(dt$WinCenter),
    tP                 = as.numeric(dt$tP),
    nSites             = as.integer(dt$nSites),
    theta_pi_persite   = as.numeric(dt$tP) / pmax(1L, as.integer(dt$nSites))
  )
}

message("[q05] Parsing ", length(files), " files for ", CHROM, " ...")
all_rows <- vector("list", length(files))
for (i in seq_along(files)) {
  all_rows[[i]] <- parse_one(files[i], sample_ids[i])
  if (i %% 50 == 0) message("[q05]   ", i, " / ", length(files))
}
dt <- rbindlist(all_rows)
message("[q05] Combined rows: ", nrow(dt),
        "   samples with data: ", length(unique(dt$sample)))

if (nrow(dt) == 0) stop("No data for ", CHROM)

# Sort for reproducibility and efficient downstream use
setorder(dt, sample, win_center_bp)
fwrite(dt, OUT, sep = "\t", quote = FALSE, compress = "gzip")
message("[q05] Wrote ", OUT)
