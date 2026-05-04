#!/usr/bin/env Rscript
# ============================================================================
# STEP_D15_karyogram.R — blocks-only horizontal karyogram (v1)
# ============================================================================
#
# PURPOSE
# -------
# A blocks-only diagnostic plot showing WHERE on the chromosome each D09
# tree node sits, alongside independent priors (SV calls, 01C registry).
# No sim_mat, no heatmap, no recursion. Answers:
#   "do my detections agree with Manta/DELLY priors geographically?"
#   "which D09 INVERSIONs fall inside 01C-flagged regions?"
#
# Distinct from D13: D13 shows blocks overlaid on the sim_mat heatmap
# (answers: "does the box sit on a bright square?"). The karyogram is
# thin, linear, multi-lane, and trivially vectorized.
#
# DESIGN
# ------
# Horizontal strip per chromosome, five lanes stacked vertically:
#
#   Lane 1 (top)     chromosome axis ruler, Mb
#   Lane 2           D09 INVERSION nodes (span >= inv_min_span_mb), filled
#                     red with per-node label, colored by nn_birth
#   Lane 3           D09 CANDIDATE nodes, filled orange, no labels
#   Lane 4           SV prior (Manta/DELLY INV), thin bars colored by caller
#   Lane 5 (bottom)  01C band boundaries (if registry available)
#
# USAGE
# -----
#   Rscript STEP_D15_karyogram.R \
#     --tree       04_output/nn_tree_C_gar_LG28.tsv \
#     --precomp    01_input/C_gar_LG28.precomp.slim.rds \
#     --sv_prior   01_input/sv_prior_C_gar_LG28.rds \
#     --blocks_01c 01_input/block_registry_C_gar_LG28.tsv.gz \
#     --outdir     04_output \
#     --chr        C_gar_LG28 \
#     --min_span   0.5
#
# --sv_prior and --blocks_01c are optional. --min_span defaults to 0.5 Mb.
# If --tree already has span_mb, we use it; else compute from end_mb -
# start_mb.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------- CLI parsing ----------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

tree_f     <- get_arg("--tree")
precomp_f  <- get_arg("--precomp")
sv_f       <- get_arg("--sv_prior")
blocks_01c_f <- get_arg("--blocks_01c")
outdir     <- get_arg("--outdir", ".")
chr_label  <- get_arg("--chr", "chr")
min_span   <- as.numeric(get_arg("--min_span", "0.5"))

if (is.na(tree_f) || !file.exists(tree_f))
  stop("--tree is required and must exist")
if (is.na(precomp_f) || !file.exists(precomp_f))
  stop("--precomp is required and must exist (used for chromosome length)")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("[karyogram] tree:    ", tree_f,     "\n")
cat("[karyogram] precomp: ", precomp_f,  "\n")
cat("[karyogram] sv_prior:", sv_f,       "\n")
cat("[karyogram] 01c:     ", blocks_01c_f, "\n")
cat("[karyogram] min_span_mb =", min_span, "\n")

# ---------- Load tree ----------
tree <- fread(tree_f)
if (!"span_mb" %in% names(tree)) tree[, span_mb := end_mb - start_mb]
cat("[karyogram] tree nodes loaded:", nrow(tree), "\n")

# Chromosome length from precomp dt
pc <- readRDS(precomp_f)
chr_len_mb <- max(pc$dt$end_bp, na.rm = TRUE) / 1e6
cat("[karyogram] chromosome length:", round(chr_len_mb, 3), "Mb\n")

# ---------- Filter for plotting (B) ----------
# INVERSION  : classification == INVERSION AND span >= min_span
# CANDIDATE  : classification == CANDIDATE (we show these as context, no label)
inv  <- tree[classification == "INVERSION" & span_mb >= min_span]
cand <- tree[classification == "CANDIDATE" & span_mb >= min_span / 2]

cat("[karyogram] INVERSIONs shown (span >=", min_span, "Mb):", nrow(inv), "\n")
cat("[karyogram] CANDIDATEs shown (span >=", min_span / 2, "Mb):", nrow(cand), "\n")

# ---------- Load SV prior (optional) ----------
sv_dt <- data.table()
if (!is.na(sv_f) && file.exists(sv_f)) {
  sv_raw <- tryCatch(readRDS(sv_f), error = function(e) {
    warning("failed to read sv_prior: ", conditionMessage(e)); NULL
  })
  if (!is.null(sv_raw)) {
    # Expected shape: a data.table with $inv (INV calls) or similar.
    # Be defensive — pull whatever looks like start/end/caller.
    sv_candidates <- list()
    if (is.list(sv_raw) && "inv" %in% names(sv_raw)) sv_candidates <- list(sv_raw$inv)
    else if (is.data.frame(sv_raw)) sv_candidates <- list(as.data.table(sv_raw))
    else if (is.list(sv_raw)) sv_candidates <- Filter(is.data.frame, sv_raw)

    for (tbl in sv_candidates) {
      if (!is.null(tbl) && nrow(tbl) > 0) {
        tbl <- as.data.table(tbl)
        # Normalise column names
        if ("POS"   %in% names(tbl) && !"start_bp" %in% names(tbl)) tbl[, start_bp := POS]
        if ("START" %in% names(tbl) && !"start_bp" %in% names(tbl)) tbl[, start_bp := START]
        if ("END"   %in% names(tbl) && !"end_bp"   %in% names(tbl)) tbl[, end_bp   := END]
        if ("CALLER" %in% names(tbl) && !"caller" %in% names(tbl))  tbl[, caller := CALLER]
        if (!"caller" %in% names(tbl)) tbl[, caller := "SV"]
        if ("start_bp" %in% names(tbl) && "end_bp" %in% names(tbl)) {
          sv_dt <- rbind(sv_dt, tbl[, .(start_bp, end_bp, caller)], fill = TRUE)
        }
      }
    }
    if (nrow(sv_dt) > 0) {
      sv_dt[, `:=`(start_mb = start_bp / 1e6, end_mb = end_bp / 1e6)]
      cat("[karyogram] SV prior rows:", nrow(sv_dt),
          " | callers:", paste(unique(sv_dt$caller), collapse = ", "), "\n")
    }
  }
}

# ---------- Load 01C block registry (optional) ----------
blk01c <- data.table()
if (!is.na(blocks_01c_f) && file.exists(blocks_01c_f)) {
  blk01c <- tryCatch(fread(blocks_01c_f), error = function(e) {
    warning("failed to read 01C registry: ", conditionMessage(e)); data.table()
  })
  cat("[karyogram] 01C registry rows:", nrow(blk01c), "\n")
  # Expect start_bp/end_bp or start_mb/end_mb
  if (nrow(blk01c) > 0 && !"start_mb" %in% names(blk01c)) {
    if ("start_bp" %in% names(blk01c)) blk01c[, start_mb := start_bp / 1e6]
    if ("end_bp"   %in% names(blk01c)) blk01c[, end_mb   := end_bp   / 1e6]
  }
}

# ---------- Lane Y-coordinates ----------
# Lanes stack top-down. Higher y = higher lane.
LANE_RULER <- 4.0
LANE_INV   <- 3.0
LANE_CAND  <- 2.2
LANE_SV    <- 1.4
LANE_01C   <- 0.6
LANE_H     <- 0.5   # half-height of each lane

# ---------- Assemble rectangles ----------
# INVERSION (lane 3): colored by nn_birth
inv_rect <- if (nrow(inv) > 0) data.table(
  xmin = inv$start_mb, xmax = inv$end_mb,
  ymin = LANE_INV - LANE_H, ymax = LANE_INV + LANE_H,
  nn_birth = inv$nn_birth, span_mb = inv$span_mb,
  label = paste0("B", inv$node_id)
) else data.table()

# CANDIDATE (lane 2)
cand_rect <- if (nrow(cand) > 0) data.table(
  xmin = cand$start_mb, xmax = cand$end_mb,
  ymin = LANE_CAND - LANE_H, ymax = LANE_CAND + LANE_H
) else data.table()

# SV (lane 4)
sv_rect <- if (nrow(sv_dt) > 0) data.table(
  xmin = sv_dt$start_mb, xmax = sv_dt$end_mb,
  ymin = LANE_SV - LANE_H * 0.4, ymax = LANE_SV + LANE_H * 0.4,
  caller = sv_dt$caller
) else data.table()

# 01C (lane 5)
blk01c_rect <- if (nrow(blk01c) > 0 && all(c("start_mb", "end_mb") %in% names(blk01c)))
  data.table(
    xmin = blk01c$start_mb, xmax = blk01c$end_mb,
    ymin = LANE_01C - LANE_H * 0.5, ymax = LANE_01C + LANE_H * 0.5
  ) else data.table()

# ---------- Plot ----------
p <- ggplot() +
  # Chromosome ruler lane (top)
  geom_segment(
    data = data.table(x = 0, xend = chr_len_mb,
                       y = LANE_RULER, yend = LANE_RULER),
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "#111827", linewidth = 0.6
  ) +
  annotate("text",
    x = chr_len_mb, y = LANE_RULER + 0.25,
    label = paste0(chr_label, ": ", round(chr_len_mb, 2), " Mb"),
    hjust = 1, vjust = 0, fontface = "bold", size = 3.2, color = "#111827"
  )

# Lane 5: 01C
if (nrow(blk01c_rect) > 0) {
  p <- p + geom_rect(
    data = blk01c_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#0d9488", alpha = 0.3, color = "#0d9488", linewidth = 0.3
  )
}

# Lane 4: SV prior — colored by caller
if (nrow(sv_rect) > 0) {
  # Manta / DELLY / other — assign colors manually
  sv_rect[, fill_c := fifelse(tolower(caller) %in% c("manta", "mantainv"), "#f97316",
                       fifelse(tolower(caller) %in% c("delly", "dellyinv"), "#2563eb",
                        "#78350f"))]
  p <- p + geom_rect(
    data = sv_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_c),
    alpha = 0.8, color = NA
  ) + scale_fill_identity()
}

# Lane 3: INVERSION — colored by nn_birth (deeper = redder)
if (nrow(inv_rect) > 0) {
  p <- p + geom_rect(
    data = inv_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = nn_birth),
    color = "#7f1d1d", linewidth = 0.3, alpha = 0.85
  ) + scale_fill_gradient(
    name = "nn_birth",
    low = "#fecaca", high = "#7f1d1d",
    limits = c(min(40, min(inv_rect$nn_birth, na.rm = TRUE)),
               max(320, max(inv_rect$nn_birth, na.rm = TRUE)))
  ) +
  geom_text(
    data = inv_rect,
    aes(x = (xmin + xmax) / 2, y = ymax + 0.15, label = label),
    size = 2.2, color = "#7f1d1d", angle = 45, hjust = 0, vjust = 0.5
  )
}

# Lane 2: CANDIDATE — pale orange
if (nrow(cand_rect) > 0) {
  p <- p + geom_rect(
    data = cand_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#fed7aa", alpha = 0.7,
    color = "#ea580c", linewidth = 0.2
  )
}

# Lane labels on the left margin
lane_labels <- data.table(
  y = c(LANE_RULER, LANE_INV, LANE_CAND, LANE_SV, LANE_01C),
  label = c("chr", "INV (span>=min)", "CAND", "SV prior", "01C bands")
)
p <- p + geom_text(
  data = lane_labels,
  aes(x = -chr_len_mb * 0.01, y = y, label = label),
  hjust = 1, vjust = 0.5, size = 2.6, fontface = "italic", color = "grey30"
)

# Cosmetics
p <- p +
  coord_cartesian(xlim = c(-chr_len_mb * 0.08, chr_len_mb * 1.02),
                   ylim = c(0, 5), clip = "off") +
  labs(
    x = paste0(chr_label, " (Mb)"),
    y = NULL,
    title = paste0(chr_label, " -- Karyogram of D09 INVERSIONs, SV priors, 01C bands"),
    subtitle = paste0(
      nrow(inv), " INV (span>=", min_span, "Mb) | ",
      nrow(cand), " CAND | ",
      if (nrow(sv_dt) > 0) paste0(nrow(sv_dt), " SV calls | ") else "",
      if (nrow(blk01c) > 0) paste0(nrow(blk01c), " 01C bands") else "no 01C"
    ),
    caption = paste(
      "Red = INVERSION (nn_birth + span filter); Orange = CANDIDATE;",
      "Blue bars = DELLY; Orange bars = Manta; Teal = 01C bands.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "#FAFAFA", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(color = "grey40", size = 9),
    plot.caption      = element_text(size = 7, hjust = 0, color = "grey50"),
    legend.position   = "right",
    legend.key.size   = unit(0.35, "cm")
  )

# ---------- Save ----------
out_base <- paste0(chr_label, "_karyogram_D09")
for (ext in c("png", "pdf")) {
  f <- file.path(outdir, paste0(out_base, ".", ext))
  tryCatch({
    ggsave(f, p, width = 14, height = 4.5, dpi = 300)
    cat("[SAVED]", f, "\n")
  }, error = function(e) cat("[FAIL]", ext, ":", conditionMessage(e), "\n"))
}

cat("[karyogram] done.\n")
