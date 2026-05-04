#!/usr/bin/env Rscript
# =============================================================================
# aggregate_ncrna_density.R
# -----------------------------------------------------------------------------
# Sibling to aggregate_repeat_density.R. Reads the canonical ncRNA GFF3s
# (tRNAscan-SE / barrnap / Rfam) and emits one per-chromosome
# <chrom>_ncrna_density.json file matching the schema the atlas's loader
# expects (tool = "ncrna_density_v1", chromosomes[].by_class shape mirroring
# the TEfull layer).
#
# Usage:
#   Rscript aggregate_ncrna_density.R \
#     --trna-gff3   /path/to/fClaHyb_Gar.tRNA.gff3 \
#     --rrna-gff3   /path/to/fClaHyb_Gar.rRNA.gff3 \
#     --ncrna-gff3  /path/to/fClaHyb_Gar.ncRNA.gff3 \
#     --windows-tsv /path/to/Gar.LG28.windows.tsv \
#     --species     "C. gariepinus" \
#     --out-dir     /path/to/out/
#
# The windows TSV is the same per-chrom window grid used by the
# repeat-density layer and the precomp JSONs (columns: chrom,
# start_bp, end_bp, center_mb). One JSON output per chrom present in
# both the windows TSV and at least one of the GFF3s.
#
# Sub-categories produced (preserving GFF3 attribute distinctions):
#   tRNA family:  tRNA_all, tRNA_HC, tRNA_pseudo, tRNA_intronic, tRNA_Sec
#   rRNA family:  rRNA_all, rRNA_5S, rRNA_5_8S, rRNA_18S, rRNA_28S, rRNA_partial
#   ncRNA family: ncRNA_all, ncRNA_snRNA, ncRNA_snoRNA, ncRNA_miRNA, ncRNA_other
# Density unit: loci/Mb (matching the atlas's TE-layer unit-consistency).
#
# Sub-class assignment heuristics:
#   tRNA_HC        : tRNAscan-SE high-confidence (no `pseudo=true` or `note=Pseudo`)
#   tRNA_pseudo    : tRNAscan-SE pseudogene flag
#   tRNA_intronic  : feature_type or attribute mentions "intron"
#   tRNA_Sec       : Name=Sec or anticodon=TCA
#   rRNA_5S/5_8S/18S/28S : barrnap Name= attribute (5S_rRNA / 5_8S_rRNA / etc.)
#   rRNA_partial   : barrnap "(partial)" tag in note/Name
#   ncRNA_snRNA    : Rfam type column = snRNA OR Name starts with U[0-9]
#   ncRNA_snoRNA   : Rfam type column = snoRNA OR scaRNA
#   ncRNA_miRNA    : Rfam type column = miRNA
#   ncRNA_other    : everything else (lncRNA, RNase, etc.)
#
# Defensive fallback: if heuristics can't classify a feature, it still goes
# into the *_all bucket (so total counts always match the GFF3 row count for
# each family).
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(data.table)
})

option_list <- list(
  make_option("--trna-gff3",  type="character", default=NULL, help="tRNAscan-SE GFF3"),
  make_option("--rrna-gff3",  type="character", default=NULL, help="barrnap GFF3"),
  make_option("--ncrna-gff3", type="character", default=NULL, help="Rfam/Infernal GFF3 (ex-tRNA/rRNA)"),
  make_option("--windows-tsv", type="character", default=NULL, help="Per-window grid TSV (chrom, start_bp, end_bp, center_mb)"),
  make_option("--species",    type="character", default="unknown", help="Species name (free text, written into JSON metadata)"),
  make_option("--out-dir",    type="character", default=".",     help="Output directory (one JSON per chrom)")
)
opt <- parse_args(OptionParser(option_list=option_list))

stopifnot(!is.null(opt$`windows-tsv`))
if (is.null(opt$`trna-gff3`) && is.null(opt$`rrna-gff3`) && is.null(opt$`ncrna-gff3`)) {
  stop("at least one of --trna-gff3 / --rrna-gff3 / --ncrna-gff3 must be provided")
}

# -----------------------------------------------------------------------------
# Read GFF3 — minimal parser; keeps cols 1 (seqid), 4-5 (start/end), 9 (attrs)
# -----------------------------------------------------------------------------
read_gff3 <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  message("[aggregate_ncrna] reading: ", path)
  dt <- fread(
    path, sep="\t", header=FALSE,
    col.names=c("seqid","source","type","start","end","score","strand","phase","attrs"),
    skip=0, comment.char="#"
  )
  # Keep only feature rows (some GFF3s repeat the header inside the body)
  dt <- dt[!is.na(seqid) & seqid != "" & !grepl("^#", seqid, fixed=FALSE)]
  dt[, mid_bp := as.integer((start + end) / 2L)]
  dt
}

# Attribute parser — KEY=VALUE pairs separated by ';'. Returns NA if not found.
attr_get <- function(attrs, key) {
  pat <- paste0("(^|;)", key, "=([^;]*)")
  m <- regmatches(attrs, regexec(pat, attrs))
  vapply(m, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))
}

# -----------------------------------------------------------------------------
# Sub-class assignment
# -----------------------------------------------------------------------------
classify_trna <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(list())
  pseudo <- grepl("pseudo=true", dt$attrs, ignore.case=TRUE) |
            grepl("Pseudo", dt$attrs)
  intron <- grepl("intron", dt$attrs, ignore.case=TRUE)
  sec    <- grepl("Name=Sec", dt$attrs) | grepl("anticodon=TCA", dt$attrs, ignore.case=TRUE)
  list(
    tRNA_all       = dt,
    tRNA_HC        = dt[!pseudo & !intron & !sec],
    tRNA_pseudo    = dt[pseudo],
    tRNA_intronic  = dt[intron],
    tRNA_Sec       = dt[sec]
  )
}

classify_rrna <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(list())
  name <- attr_get(dt$attrs, "Name")
  partial <- grepl("partial", dt$attrs, ignore.case=TRUE)
  list(
    rRNA_all     = dt,
    rRNA_5S      = dt[grepl("(^|_)5S(_|$)",       name)],
    rRNA_5_8S    = dt[grepl("5[._]8S",            name)],
    rRNA_18S     = dt[grepl("(^|_)18S(_|$)",      name)],
    rRNA_28S     = dt[grepl("(^|_)28S(_|$)",      name)],
    rRNA_partial = dt[partial]
  )
}

classify_ncrna <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(list())
  ftype <- dt$type
  name  <- attr_get(dt$attrs, "Name")
  is_sn  <- grepl("snRNA",   ftype, ignore.case=TRUE) |
            grepl("snRNA",   name,  ignore.case=TRUE) |
            grepl("^U[0-9]", name)
  is_sno <- grepl("snoRNA",  ftype, ignore.case=TRUE) |
            grepl("scaRNA",  ftype, ignore.case=TRUE) |
            grepl("snoRNA",  name,  ignore.case=TRUE) |
            grepl("scaRNA",  name,  ignore.case=TRUE)
  is_mi  <- grepl("miRNA",   ftype, ignore.case=TRUE) |
            grepl("miRNA",   name,  ignore.case=TRUE) |
            grepl("MIR-",    name,  ignore.case=TRUE)
  list(
    ncRNA_all    = dt,
    ncRNA_snRNA  = dt[is_sn],
    ncRNA_snoRNA = dt[is_sno],
    ncRNA_miRNA  = dt[is_mi],
    ncRNA_other  = dt[!is_sn & !is_sno & !is_mi]
  )
}

# -----------------------------------------------------------------------------
# Per-window density: number of feature midpoints landing in [start_bp, end_bp]
# divided by window size in Mb. Vectorized via cut() into the window grid.
# -----------------------------------------------------------------------------
densify <- function(features_dt, win_dt) {
  if (is.null(features_dt) || nrow(features_dt) == 0) {
    return(rep(0.0, nrow(win_dt)))
  }
  # Per-chrom: assign each feature to a window via binary search on start_bp
  feat <- features_dt[seqid == win_dt$chrom[1L]]
  if (nrow(feat) == 0) return(rep(0.0, nrow(win_dt)))
  win_starts <- win_dt$start_bp
  win_ends   <- win_dt$end_bp
  win_size_mb <- (win_ends - win_starts) / 1e6
  win_size_mb[win_size_mb <= 0] <- 1e-9     # guard
  counts <- integer(nrow(win_dt))
  for (i in seq_len(nrow(feat))) {
    idx <- findInterval(feat$mid_bp[i], win_starts)
    if (idx >= 1 && idx <= length(counts) && feat$mid_bp[i] <= win_ends[idx]) {
      counts[idx] <- counts[idx] + 1L
    }
  }
  counts / win_size_mb     # loci/Mb
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
trna_dt  <- read_gff3(opt$`trna-gff3`)
rrna_dt  <- read_gff3(opt$`rrna-gff3`)
ncrna_dt <- read_gff3(opt$`ncrna-gff3`)

trna_buckets  <- classify_trna(trna_dt)
rrna_buckets  <- classify_rrna(rrna_dt)
ncrna_buckets <- classify_ncrna(ncrna_dt)
all_buckets   <- c(trna_buckets, rrna_buckets, ncrna_buckets)

windows <- fread(opt$`windows-tsv`)
required_cols <- c("chrom", "start_bp", "end_bp", "center_mb")
missing <- setdiff(required_cols, colnames(windows))
if (length(missing) > 0) stop("windows TSV missing columns: ", paste(missing, collapse=","))

dir.create(opt$`out-dir`, showWarnings=FALSE, recursive=TRUE)
chroms <- unique(windows$chrom)
for (ch in chroms) {
  win_ch <- windows[chrom == ch][order(start_bp)]
  if (nrow(win_ch) == 0) next
  by_class <- list()
  for (nm in names(all_buckets)) {
    feats <- all_buckets[[nm]]
    by_class[[nm]] <- list(
      densities = round(densify(feats, win_ch), 4)
    )
  }
  obj <- list(
    tool            = "ncrna_density_v1",
    schema_version  = 1L,
    species         = opt$species,
    n_chromosomes   = 1L,
    n_classes       = length(by_class),
    classes         = names(by_class),
    default_class   = "tRNA_all",
    generated_at    = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
    chromosomes     = list(list(
      chrom              = ch,
      n_windows          = nrow(win_ch),
      window_centers_mb  = win_ch$center_mb,
      window_start_bp    = as.integer(win_ch$start_bp),
      window_end_bp      = as.integer(win_ch$end_bp),
      by_class           = by_class
    ))
  )
  out_path <- file.path(opt$`out-dir`, paste0(ch, "_ncrna_density.json"))
  write_json(obj, out_path, auto_unbox=TRUE, digits=NA)
  message("[aggregate_ncrna] wrote ", out_path,
          "  (", nrow(win_ch), " windows, ", length(by_class), " classes)")
}

message("[aggregate_ncrna] done")
