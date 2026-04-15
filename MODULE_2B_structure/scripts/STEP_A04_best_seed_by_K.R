#!/usr/bin/env Rscript
# =============================================================================
# STEP_A04_best_seed_by_K.R — Select best NGSadmix seed per K value
#
# v9.0 REWIRED:
#   - Auto-registers ancestry cluster groups to sample_registry
#   - Uses sample_id (CGA names) in all outputs
#   - Exports ancestry tables with sample_id column
#
# Usage:
#   Rscript scripts/STEP_A04_best_seed_by_K.R \
#     --run-dir  ${MODULE2B_RESULTS}/wholegenome_thin500_all226 \
#     --eval-dir ${MODULE2B_EVALADMIX} \
#     --sample-file ${SAMPLE_LIST} \
#     --file-prefix wholegenome_thin500_all226 \
#     --k-min 2 --k-max 20 --seeds 1,2,3,4,5
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(optparse) })

option_list <- list(
  make_option("--run-dir",      type = "character", help = "Dir with .qopt/.fopt.gz/.log"),
  make_option("--eval-dir",     type = "character", default = NULL),
  make_option("--sample-file",  type = "character", help = "One-sample-per-line canonical list"),
  make_option("--file-prefix",  type = "character", help = "e.g. wholegenome_thin500_all226"),
  make_option("--k-min",        type = "integer", default = 2),
  make_option("--k-max",        type = "integer", default = 20),
  make_option("--seeds",        type = "character", default = "1,2,3,4,5"),
  make_option("--palette-name", type = "character", default = "catfish_ngsadmix_v1"),
  make_option("--best-k",       type = "integer", default = 8,
              help = "K value to use for registry group registration"),
  make_option("--bridge",       type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

run_dir  <- opt[["run-dir"]];   stopifnot(!is.null(run_dir), dir.exists(run_dir))
eval_dir <- opt[["eval-dir"]]
sfp      <- opt[["sample-file"]]; stopifnot(!is.null(sfp), file.exists(sfp))
fp       <- opt[["file-prefix"]]; stopifnot(!is.null(fp))
K_vals   <- seq(opt[["k-min"]], opt[["k-max"]])
seeds    <- as.integer(strsplit(opt$seeds, ",")[[1]])
pal_name <- opt[["palette-name"]]
BEST_K   <- opt[["best-k"]]

# ── Load bridge (optional but recommended for registry) ─────────────────────
bridge <- opt$bridge %||% Sys.getenv("LOAD_BRIDGE", "")
has_bridge <- FALSE
if (nzchar(bridge) && file.exists(bridge)) {
  source(bridge)
  has_bridge <- TRUE
} else {
  # Try auto-detect
  for (p in c("utils/load_bridge.R", "../utils/load_bridge.R")) {
    if (file.exists(p)) { source(p); has_bridge <- TRUE; break }
  }
}

# ── Palette ──────────────────────────────────────────────────────────────────
pal_k20 <- c("#3B6FA0","#CF6839","#C44E52","#6BA08E","#5A8F4A",
             "#C9A83E","#8B6DAD","#E8919C","#5C7A3A","#B07850",
             "#4A8C9F","#9E6B8A","#7A8B3C","#C47A5E","#5E7FAA",
             "#A06B4F","#6B9E7A","#B0855A","#7E6EA0","#A0856B")
get_pal <- function(K) if (K <= 20) pal_k20[1:K] else colorRampPalette(pal_k20)(K)

# ── Samples ──────────────────────────────────────────────────────────────────
samples <- trimws(readLines(sfp, warn = FALSE)); samples <- samples[nzchar(samples)]
n <- length(samples); cat("[INFO] Samples:", n, "\n")

# ── Helpers ──────────────────────────────────────────────────────────────────
stem <- function(K, seed) sprintf("%s_K%02d_seed%d", fp, K, seed)
best_stem <- function(K) sprintf("%s_K%02d_best", fp, K)

ll_from_log <- function(f) {
  if (!file.exists(f)) return(NA_real_)
  x <- grep("like", readLines(f, warn = FALSE), ignore.case = TRUE, value = TRUE)
  if (!length(x)) return(NA_real_)
  nums <- regmatches(tail(x,1), gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", tail(x,1)))[[1]]
  if (!length(nums)) return(NA_real_)
  as.numeric(tail(nums, 1))
}

eval_metrics <- function(f) {
  null <- list(mar=NA_real_, max_ar=NA_real_, mr=NA_real_)
  if (is.null(f) || !file.exists(f)) return(null)
  M <- tryCatch(as.matrix(read.table(f, header=FALSE)), error=function(e) NULL)
  if (is.null(M) || nrow(M)!=ncol(M)) return(null)
  diag(M) <- NA
  list(mar=mean(abs(M),na.rm=TRUE), max_ar=max(abs(M),na.rm=TRUE), mr=mean(M,na.rm=TRUE))
}

read_qopt <- function(f, K) {
  ln <- trimws(readLines(f, warn=FALSE)); ln <- ln[nzchar(ln)]
  do.call(rbind, lapply(strsplit(ln, "[[:space:]]+"), as.numeric))
}

# ── Scan all runs ────────────────────────────────────────────────────────────
all <- rbindlist(lapply(K_vals, function(K) rbindlist(lapply(seeds, function(s) {
  st <- stem(K, s)
  qf <- file.path(run_dir, paste0(st, ".qopt"))
  ff <- file.path(run_dir, paste0(st, ".fopt.gz"))
  lf <- file.path(run_dir, paste0(st, ".log"))
  cf <- if (!is.null(eval_dir)) file.path(eval_dir, paste0(st, ".corres.txt")) else NULL
  em <- eval_metrics(cf)
  nr <- if (file.exists(qf)) length(readLines(qf, warn=FALSE)) else NA_integer_
  data.table(K=K, seed=s, stem=st,
    qopt_ok=file.exists(qf), fopt_ok=file.exists(ff), log_ok=file.exists(lf),
    cor_ok=!is.null(cf)&&file.exists(cf), nrows=nr, rows_match=!is.na(nr)&nr==n,
    loglik=ll_from_log(lf), mean_abs_resid=em$mar, max_abs_resid=em$max_ar,
    mean_raw_resid=em$mr, qf=qf, ff=ff, lf=lf, cf=ifelse(is.null(cf),NA_character_,cf))
}))))

all[, usable := qopt_ok & fopt_ok & log_ok & rows_match]

# ── Rank + select ────────────────────────────────────────────────────────────
r <- copy(all)
r[is.na(loglik), ll_r := -Inf]; r[!is.na(loglik), ll_r := loglik]
r[is.na(mean_abs_resid), mar_r := Inf]; r[!is.na(mean_abs_resid), mar_r := mean_abs_resid]
setorder(r, K, -ll_r, mar_r, seed)
best <- r[usable==TRUE, .SD[1], by=K]
if (!nrow(best)) stop("No usable runs found")

# ── Write outputs ────────────────────────────────────────────────────────────
fwrite(all, file.path(run_dir, paste0(fp, "_all_metrics.tsv")), sep="\t", quote=FALSE, na="NA")
fwrite(best[, .(K,seed,stem,loglik,mean_abs_resid,max_abs_resid,mean_raw_resid,nrows,qf,ff,lf,cf)],
       file.path(run_dir, paste0(fp, "_best_seed_by_K.tsv")), sep="\t", quote=FALSE, na="NA")

cat("\nBest seed by K:\n"); print(best[, .(K, seed, loglik, mean_abs_resid)])

# ── Symlink best files ───────────────────────────────────────────────────────
for (i in seq_len(nrow(best))) {
  K <- best$K[i]; src <- best$stem[i]; dst <- best_stem(K)
  for (ext in c(".qopt", ".fopt.gz", ".log")) {
    sf <- file.path(run_dir, paste0(src, ext)); df <- file.path(run_dir, paste0(dst, ext))
    if (file.exists(sf)) { if (file.exists(df)) file.remove(df); file.symlink(normalizePath(sf), df) }
  }
  if (!is.null(eval_dir) && !is.na(best$cf[i]) && file.exists(best$cf[i])) {
    sf <- best$cf[i]; df <- file.path(eval_dir, paste0(dst, ".corres.txt"))
    if (file.exists(df)) file.remove(df); file.symlink(normalizePath(sf), df)
  }
}

# ── Palette table ────────────────────────────────────────────────────────────
pal_dt <- rbindlist(lapply(unique(best$K), function(K) {
  cols <- get_pal(K)
  data.table(K=K, cluster_index=seq_len(K), cluster_label=paste0("Q",1:K),
             color_hex=cols, palette_name=pal_name)
}))
fwrite(pal_dt, file.path(run_dir, paste0(fp, "_palette.tsv")), sep="\t", quote=FALSE)

# ── Sample order (with sample_id = CGA names) ───────────────────────────────
fwrite(data.table(sample_index=seq_len(n), sample_id=samples),
       file.path(run_dir, paste0(fp, "_sample_order.tsv")), sep="\t", quote=FALSE)

# ── Per-sample ancestry ─────────────────────────────────────────────────────
anc <- rbindlist(lapply(seq_len(nrow(best)), function(i) {
  K <- best$K[i]; qf <- best$qf[i]
  if (!file.exists(qf)) return(NULL)
  qm <- read_qopt(qf, K); stopifnot(nrow(qm)==n)
  dt <- data.table(K=K, sample_index=1:n, sample_id=samples)
  for (j in 1:K) set(dt, j=paste0("Q",j), value=qm[,j])
  qc <- paste0("Q", 1:K)
  dt[, qmax := apply(.SD,1,max), .SDcols=qc]
  dt[, cluster_index := max.col(as.matrix(.SD), ties.method="first"), .SDcols=qc]
  dt[, cluster_label := paste0("Q", cluster_index)]
  dt <- merge(dt, pal_dt[pal_dt$K==K, .(cluster_index, color_hex)],
              by="cluster_index", all.x=TRUE, sort=FALSE)
  dt
}), fill=TRUE)

fwrite(anc, file.path(run_dir, paste0(fp, "_sample_ancestry.tsv")),
       sep="\t", quote=FALSE, na="NA")

# =============================================================================
# v9.0: AUTO-REGISTER ANCESTRY CLUSTER GROUPS TO REGISTRY
# =============================================================================

if (has_bridge && !is.null(reg) && BEST_K %in% best$K) {
  cat("[A04] Registering ancestry cluster groups for K=", BEST_K, "...\n")

  anc_bestK <- anc[K == BEST_K]
  if (nrow(anc_bestK) > 0 && "cluster_index" %in% names(anc_bestK)) {
    # Register each ancestry cluster
    for (ci in sort(unique(anc_bestK$cluster_index))) {
      members <- anc_bestK[cluster_index == ci, sample_id]
      if (length(members) >= 2) {
        grp_id <- paste0("ancestry_K", BEST_K, "_Q", ci)
        register_group(
          group_id    = grp_id,
          sample_ids  = members,
          description = paste0("NGSadmix K=", BEST_K, " Q", ci,
                                " major cluster (n=", length(members), ")")
        )
      }
    }

    # Register the complete ancestry assignment as metadata on master
    if (!is.null(reg)) {
      tryCatch({
        reg$add_master_column(
          paste0("Q_major_K", BEST_K),
          anc_bestK$cluster_label,
          by_sample_id = anc_bestK$sample_id
        )
        cat("[A04] Added Q_major_K", BEST_K, " column to sample_master\n")
      }, error = function(e) {
        cat("[A04] Could not add column to master:", conditionMessage(e), "\n")
      })
    }
  }
}

cat("\n[DONE] STEP_A04 outputs in:", run_dir, "\n")
