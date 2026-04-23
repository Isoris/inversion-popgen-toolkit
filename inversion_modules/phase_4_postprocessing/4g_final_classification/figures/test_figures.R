#!/usr/bin/env Rscript
# test_figures.R — smoke test for the three figures added this turn:
#   fig4e_gene_count_by_biotype.R
#   fig4e_class_along_inversion.R
#   fig4e_cumulative_burden_per_group.R
#
# Creates synthetic inputs under /tmp/fig_smoke/, invokes each script
# via system2(), and checks that the expected output files materialize.
# Does NOT verify plot correctness — this is a pipeline wiring test
# only. Eye-check the PNGs after.
#
# Usage:
#   cd inversion_modules/phase_4_postprocessing/4g_final_classification/figures/
#   Rscript test_figures.R

suppressPackageStartupMessages({ library(data.table) })

set.seed(1)
base <- "/tmp/fig_smoke"
dir.create(base, recursive = TRUE, showWarnings = FALSE)

PASS <- TRUE
report <- function(ok, msg) {
  cat(if (ok) "[OK]   " else "[FAIL] ", msg, "\n")
  if (!ok) PASS <<- FALSE
}

# ── 1. Synthetic gene map ──
gmap <- data.table(
  gene_id    = sprintf("g%04d", 1:500),
  chr        = sample(c("LG01", "LG02", "LG03"), 500, replace = TRUE),
  gene_start = sample(1:10000000, 500),
  biotype    = sample(c("protein_coding", "pseudogene", "lncRNA", "miRNA",
                        "snoRNA", "tRNA"), 500, replace = TRUE,
                      prob = c(0.6, 0.15, 0.12, 0.05, 0.05, 0.03))
)
gmap[, gene_end := gene_start + sample(500:50000, 500, replace = TRUE)]
gene_map_path <- file.path(base, "gene_map.tsv")
fwrite(gmap, gene_map_path, sep = "\t")

# ── 2. Synthetic regions ──
regions <- data.table(
  region_id = c("ROH_A", "ROH_B", "INV_A"),
  chrom     = c("LG01", "LG02", "LG03"),
  start     = c(1000000, 500000, 2000000),
  end       = c(5000000, 3000000, 7000000),
  category  = c("ROH", "ROH", "inversion"),
  label     = c("ROH\nLG01\n4.0 Mb", "ROH\nLG02\n2.5 Mb", "INV\nLG03\n5.0 Mb")
)
regions_path <- file.path(base, "regions.tsv")
fwrite(regions, regions_path, sep = "\t")

# ── 3. Synthetic inversion candidate + per-chrom windows + karyo calls ──
cand_bed <- data.table(chrom = "LG01", start = 1000000, end = 3000000)
cand_bed_path <- file.path(base, "candidate.bed")
fwrite(cand_bed, cand_bed_path, sep = "\t", col.names = FALSE)

# 300 windows on LG01, 10 kb each
win_coords <- data.table(
  chrom            = "LG01",
  global_window_id = 1:300,
  start_bp         = seq(0, by = 10000, length.out = 300),
  end_bp           = seq(10000, by = 10000, length.out = 300)
)
win_coords_path <- file.path(base, "window_coords.tsv")
fwrite(win_coords, win_coords_path, sep = "\t")

# Karyotype calls: 20 samples, some INV_INV runs, some INV_nonINV
n_samples <- 20
karyo_rows <- list()
for (si in seq_len(n_samples)) {
  sid <- sprintf("S%02d", si)
  # 1-3 runs per sample in the inversion window region (100..300)
  n_runs <- sample(1:3, 1)
  for (r in seq_len(n_runs)) {
    ws <- sample(80:280, 1)
    we <- min(300L, ws + sample(10:40, 1))
    karyo_rows[[length(karyo_rows) + 1]] <- data.table(
      sample_id    = sid,
      window_start = ws,
      window_end   = we,
      n_windows    = we - ws + 1L,
      call         = sample(c("INV_INV", "INV_nonINV"), 1),
      mean_rank    = runif(1, 0.05, 0.95)
    )
  }
}
karyo <- rbindlist(karyo_rows)
karyo_path <- file.path(base, "karyo_calls.tsv.gz")
fwrite(karyo, karyo_path, sep = "\t")

# ── 4. Synthetic burden-per-ROH + sample groups ──
n_roh_samples <- 30
roh_rows <- list()
for (si in seq_len(n_roh_samples)) {
  sid <- sprintf("F%02d", si)
  n_roh <- sample(3:10, 1)
  for (r in seq_len(n_roh)) {
    chr <- sample(c("LG01", "LG02", "LG03"), 1)
    rs <- sample(1:8000000, 1)
    rlen <- sample(c(1e6, 2e6, 5e6, 10e6), 1)
    re <- rs + rlen
    n_del <- sample(0:15, 1)
    n_genes <- sample(0:8, 1)
    genes <- if (n_genes == 0) "none" else {
      paste(sprintf("g%04d", sample(1:500, n_genes)), collapse = ";")
    }
    roh_rows[[length(roh_rows) + 1]] <- data.table(
      sample_id = sid, chr = chr,
      roh_start = rs, roh_end = re, roh_length = rlen,
      n_del_variants = n_del,
      n_del_hom = sample(0:n_del, 1),
      sum_priority_score = round(runif(1, 0, 100), 1),
      n_classA = sample(0:n_del, 1),
      n_classB = sample(0:n_del, 1),
      genes_hit = genes
    )
  }
}
burden <- rbindlist(roh_rows)
burden_path <- file.path(base, "burden_per_roh.tsv")
fwrite(burden, burden_path, sep = "\t")

# Sample groups: 4 groups, loosely stratified
group_ids <- c("GroupA", "GroupB", "GroupC", "GroupD")
grp <- data.table(
  sample_id = sprintf("F%02d", 1:n_roh_samples),
  group_id  = sample(group_ids, n_roh_samples, replace = TRUE)
)
grp_path <- file.path(base, "sample_groups.tsv")
fwrite(grp, grp_path, sep = "\t", col.names = FALSE)

# ── Run each script and check outputs ──

run_script <- function(script, args, expected_outputs) {
  cat("\n--- running", script, "---\n")
  status <- system2("Rscript", c(script, args))
  report(status == 0, sprintf("%s exited 0 (got %d)", script, status))
  for (f in expected_outputs) {
    report(file.exists(f), sprintf("%s produced %s", script, f))
  }
}

# fig4e_gene_count_by_biotype
out1 <- file.path(base, "biotype_table.pdf")
run_script(
  "fig4e_gene_count_by_biotype.R",
  c("--gene_map", gene_map_path,
    "--regions",  regions_path,
    "--out",      out1),
  c(out1,
    paste0(tools::file_path_sans_ext(out1), ".counts.tsv"))
)

# fig4e_class_along_inversion
out2 <- file.path(base, "class_along.pdf")
run_script(
  "fig4e_class_along_inversion.R",
  c("--candidate_bed", cand_bed_path,
    "--karyo_calls",   karyo_path,
    "--window_coords", win_coords_path,
    "--out",           out2,
    "--flank_kb",      "200"),
  c(out2,
    paste0(tools::file_path_sans_ext(out2), ".png"),
    paste0(tools::file_path_sans_ext(out2), ".class_proportions.tsv"))
)

# fig4e_cumulative_burden_per_group
out3 <- file.path(base, "cumburden.pdf")
run_script(
  "fig4e_cumulative_burden_per_group.R",
  c("--burden_roh",    burden_path,
    "--sample_groups", grp_path,
    "--out",           out3,
    "--sort_by",       "length",
    "--attribution",   "first"),
  c(out3,
    paste0(tools::file_path_sans_ext(out3), ".png"),
    paste0(tools::file_path_sans_ext(out3), ".cumulative.tsv"))
)

cat("\n===========================\n")
cat(if (PASS) "ALL CHECKS PASSED" else "FAILURES — see above", "\n")
cat("Outputs in:", base, "\n")
cat("Files:\n")
print(list.files(base))
quit(save = "no", status = if (PASS) 0 else 1)
