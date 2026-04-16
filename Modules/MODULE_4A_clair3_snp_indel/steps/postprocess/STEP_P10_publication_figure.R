#!/usr/bin/env Rscript
# ============================================================================
# STEP10_publication_figure.R
#
# Composite figure: Rescuing biologically supported variants from
# Clair3-filtered calls.
#
# Layout (matching mockup):
#   Row 1: (a) Retention funnel   (b) Counts by final class per chromosome
#   Row 2: (c) Stacked bar        (d) Proportions  (e) Rescue composition  (f) Violins
#   Row 3: (g) Weak-indel support  (h) Genome position track
#
# Usage:
#   Rscript STEP10_publication_figure.R <classified_tsv> <weak_tsv> <outdir> \
#       [<chr_length_bp>] [<sample_label>]
#
# Example:
#   Rscript STEP10_publication_figure.R \
#       final_variant_classification.tsv \
#       weak_indel_candidates.tsv \
#       plots/ \
#       33300000 \
#       "CGA009"
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
  library(gridExtra)
  library(forcats)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript STEP10_publication_figure.R <classified.tsv> <weak.tsv> <outdir> [chr_len_bp] [sample_label]")
}

classified_file <- args[1]
weak_file       <- args[2]
outdir          <- args[3]
chr_len_bp      <- if (length(args) >= 4) as.numeric(args[4]) else 33300000
sample_label    <- if (length(args) >= 5) args[5] else "CGA009"
chrom_label     <- "C_gar_LG28"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ────────────────────────────────────────────────────────────────

df <- fread(classified_file)
weak <- if (file.exists(weak_file)) fread(weak_file) else data.table()

# ── Colour palette (publication: dark-blue family) ───────────────────────────

pal <- c(
  "Strict PASS"           = "#1a2744",
  "Rescued strong SNP"    = "#3a6fa0",
  "Rescued strong indel"  = "#5b9bd5",
  "Rescued relaxed indel" = "#a3c4e0",
  "Pending weak"          = "#c8c8c8"
)

# map FINAL_CLASS + RESCUE_REASON to display labels
df <- df %>%
  mutate(
    DISPLAY_CLASS = case_when(
      FINAL_CLASS == "STRICT_PASS"                                      ~ "Strict PASS",
      FINAL_CLASS == "RESCUED_STRONG_SINGLE_SAMPLE" & grepl("snp", RESCUE_REASON) ~ "Rescued strong SNP",
      FINAL_CLASS == "RESCUED_STRONG_SINGLE_SAMPLE" & RESCUE_REASON == "strong_indel" ~ "Rescued strong indel",
      FINAL_CLASS == "RESCUED_STRONG_SINGLE_SAMPLE" & RESCUE_REASON == "relaxed_indel" ~ "Rescued relaxed indel",
      FINAL_CLASS == "RESCUED_POPULATION_SIGNATURE" ~ "Rescued strong indel",
      TRUE ~ "Pending weak"
    ),
    DISPLAY_CLASS = factor(DISPLAY_CLASS, levels = names(pal))
  )

# basic counts
n_raw <- nrow(df)
n_pass <- sum(df$FINAL_CLASS == "STRICT_PASS")
n_filtered <- sum(df$STATUS == "FILTERED")
n_rescued <- sum(grepl("RESCUED", df$FINAL_CLASS))
n_pending <- n_raw - n_pass - n_rescued

# rescued breakdown
rescued_df <- df %>% filter(grepl("RESCUED", FINAL_CLASS))
n_rescued_snp <- sum(rescued_df$VAR_TYPE == "SNP")
n_rescued_del <- sum(rescued_df$VAR_TYPE == "DEL")
n_rescued_ins <- sum(rescued_df$VAR_TYPE == "INS")
n_strong_indel <- sum(rescued_df$RESCUE_REASON == "strong_indel")
n_relaxed_indel <- sum(rescued_df$RESCUE_REASON == "relaxed_indel")
n_strong_snp <- sum(grepl("snp", rescued_df$RESCUE_REASON))

# ── base theme ───────────────────────────────────────────────────────────────

theme_pub <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 5)
  )

fmt_k <- function(x) {
  ifelse(x >= 1000, paste0(round(x/1000, 0), "k"), as.character(x))
}

# ============================================================================
# Panel (a): Retention funnel (built as a table/annotation, not a standard plot)
# ============================================================================

funnel_data <- data.frame(
  step = c("Raw discovered", "Strict PASS", "Filtered out", "Rescued from filtered", "Remaining pending"),
  subtitle = c("All Clair3 candidates", "Clair3 hard filters", "Did not meet hard filters",
                "Biologically supported", "Low/no support"),
  count = c(n_raw, n_pass, n_filtered, n_rescued, n_pending),
  pct = c(100, round(100*n_pass/n_raw, 1), round(100*n_filtered/n_raw, 1),
          round(100*n_rescued/n_raw, 1), round(100*n_pending/n_raw, 1)),
  y = 5:1,
  stringsAsFactors = FALSE
)

# widths proportional to count (for funnel shape)
funnel_data$width <- funnel_data$count / max(funnel_data$count) * 4

pa <- ggplot(funnel_data, aes(x = 0, y = y)) +
  geom_tile(aes(width = width, height = 0.7),
            fill = c("#1a2744", "#3a6fa0", "#8899aa", "#5b9bd5", "#c8c8c8"),
            color = NA) +
  geom_text(aes(x = -2.5, label = step), hjust = 0, fontface = "bold", size = 3) +
  geom_text(aes(x = -2.5, y = y - 0.25, label = subtitle),
            hjust = 0, size = 2.3, color = "grey50") +
  geom_text(aes(x = 2.8, label = paste0(format(count, big.mark = ","),
                                          ifelse(y < 5, paste0("\n(", pct, "%)"), ""))),
            hjust = 1, size = 3, fontface = "plain") +
  scale_y_continuous(limits = c(0.3, 5.7)) +
  scale_x_continuous(limits = c(-2.8, 3)) +
  labs(title = "a   Retention funnel") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
        plot.margin = margin(5, 5, 5, 5))

# ============================================================================
# Panel (b): Counts by final class – grouped bar (chromosome-level)
# ============================================================================

# For pilot single-chromosome, show as a single group
class_counts <- df %>%
  count(DISPLAY_CLASS, name = "N") %>%
  mutate(DISPLAY_CLASS = factor(DISPLAY_CLASS, levels = rev(names(pal))))

pb <- ggplot(class_counts, aes(x = chrom_label, y = N, fill = DISPLAY_CLASS)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "b   Counts by final class",
       subtitle = paste0("Pilot chromosome (", chrom_label, ")"),
       x = NULL, y = "Number of variants") +
  theme_pub +
  theme(legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))

# ============================================================================
# Panel (c): Stacked bar summary
# ============================================================================

stacked_data <- data.frame(
  category = factor(c("Strict PASS", "Rescued SNP\n(strong)", "Rescued indel",
                       "Relaxed indel", "Pending weak"),
                     levels = c("Pending weak", "Relaxed indel", "Rescued indel",
                                "Rescued SNP\n(strong)", "Strict PASS")),
  count = c(n_pass, n_strong_snp, n_strong_indel, n_relaxed_indel, n_pending),
  fill_col = c("Strict PASS", "Rescued strong SNP", "Rescued strong indel",
                "Rescued relaxed indel", "Pending weak"),
  stringsAsFactors = FALSE
)

pc <- ggplot(stacked_data, aes(x = 1, y = count, fill = fill_col)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = paste0(fmt_k(count))),
            position = position_stack(vjust = 0.5), size = 2.5, color = "white") +
  annotate("text", x = 1, y = n_raw + 2000,
           label = paste0(format(n_raw, big.mark = ","), "\n(raw)"),
           size = 3, fontface = "bold") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "c   One-bar summary",
       subtitle = chrom_label, x = NULL, y = "Number of variants") +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ============================================================================
# Panel (d): Proportions
# ============================================================================

prop_data <- data.frame(
  category = factor(c("Strict PASS", "Rescued SNP", "Rescued indel",
                       "Rescued relaxed", "Pending weak"),
                     levels = rev(c("Strict PASS", "Rescued SNP", "Rescued indel",
                                     "Rescued relaxed", "Pending weak"))),
  pct = c(n_pass, n_strong_snp, n_strong_indel + n_relaxed_indel, 0, n_pending) / n_raw * 100,
  fill_col = c("Strict PASS", "Rescued strong SNP", "Rescued strong indel",
                "Rescued relaxed indel", "Pending weak"),
  stringsAsFactors = FALSE
)
# recalculate to match the actual breakdown
prop_data$pct <- c(100*n_pass/n_raw, 100*n_strong_snp/n_raw,
                    100*(n_strong_indel + n_relaxed_indel)/n_raw,
                    0, 100*n_pending/n_raw)

pd_plot <- ggplot(prop_data %>% filter(pct > 0),
                   aes(x = 1, y = pct, fill = fill_col)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = paste0(round(pct, 1), "%")),
            position = position_stack(vjust = 0.5), size = 2.5, color = "white") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "d   Proportions", x = NULL, y = "Proportion of raw calls") +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ============================================================================
# Panel (e): Rescue composition by variant type
# ============================================================================

rescue_comp <- rescued_df %>%
  count(VAR_TYPE, name = "N") %>%
  mutate(pct = round(100 * N / sum(N), 1))

pe <- ggplot(rescue_comp, aes(x = VAR_TYPE, y = N, fill = VAR_TYPE)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(format(N, big.mark = ","), "\n(", pct, "%)")),
            vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("SNP" = "#3a6fa0", "DEL" = "#5b9bd5", "INS" = "#a3c4e0")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = label_comma()) +
  labs(title = "e   What did rescue add?",
       subtitle = paste0("Composition of the ", format(n_rescued, big.mark = ","), " rescued"),
       x = NULL, y = "Number of rescued") +
  theme_pub

# ============================================================================
# Panel (f): Violin distributions of key metrics by class
# ============================================================================

violin_df <- df %>%
  filter(DISPLAY_CLASS != "Pending weak") %>%
  select(DISPLAY_CLASS, QUAL, DP, AF1, AD_ALT) %>%
  pivot_longer(cols = c(QUAL, DP, AF1, AD_ALT), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

# add pending weak as background
pending_violin <- df %>%
  filter(DISPLAY_CLASS == "Pending weak") %>%
  select(DISPLAY_CLASS, QUAL, DP, AF1, AD_ALT) %>%
  pivot_longer(cols = c(QUAL, DP, AF1, AD_ALT), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

violin_all <- bind_rows(violin_df, pending_violin)

pf <- ggplot(violin_all, aes(x = DISPLAY_CLASS, y = value, fill = DISPLAY_CLASS)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.3) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal) +
  labs(title = "f   Distributions of key metrics", x = NULL, y = NULL) +
  theme_pub +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face = "bold", size = 9))

# ============================================================================
# Panel (g): Weak-indel read support histogram
# ============================================================================

if (nrow(weak) > 0 && "N_READS_SUPPORT_INDEL" in names(weak)) {
  weak_hist <- weak %>%
    mutate(support_bin = case_when(
      N_READS_SUPPORT_INDEL == 0 ~ "0",
      N_READS_SUPPORT_INDEL == 1 ~ "1",
      N_READS_SUPPORT_INDEL == 2 ~ "2",
      N_READS_SUPPORT_INDEL == 3 ~ "3",
      TRUE ~ ">3"
    )) %>%
    mutate(support_bin = factor(support_bin, levels = c("0", "1", "2", "3", ">3"))) %>%
    count(support_bin, name = "N") %>%
    mutate(pct = round(100 * N / sum(N), 1))

  pg <- ggplot(weak_hist, aes(x = support_bin, y = N)) +
    geom_col(fill = "#5b9bd5", width = 0.6) +
    geom_text(aes(label = paste0(format(N, big.mark = ","), "\n(", pct, "%)")),
              vjust = -0.2, size = 2.5) +
    annotate("text", x = 4.5, y = max(weak_hist$N) * 0.85,
             label = "Many filtered indels still\ncarry 1-3 reads of support.",
             size = 2.5, hjust = 1, fontface = "italic", color = "grey40") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = label_comma()) +
    labs(title = "g   Weak indel read support",
         x = "N_READS_SUPPORT_INDEL", y = "Number of weak indels") +
    theme_pub
} else {
  pg <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No weak indel data", size = 4) +
    labs(title = "g   Weak indel read support") +
    theme_void()
}

# ============================================================================
# Panel (h): Genome-position density track
# ============================================================================

bin_size <- 100000  # 100 kb bins

density_df <- df %>%
  mutate(bin = floor((POS - 1) / bin_size) * bin_size + bin_size / 2) %>%
  mutate(track = case_when(
    DISPLAY_CLASS == "Strict PASS" ~ "Strict PASS",
    grepl("Rescued", DISPLAY_CLASS) ~ "Rescued",
    TRUE ~ "Unresolved weak"
  )) %>%
  count(bin, track, name = "N")

density_df$track <- factor(density_df$track,
                            levels = c("Strict PASS", "Rescued", "Unresolved weak"))

track_colors <- c("Strict PASS" = "#1a2744", "Rescued" = "#5b9bd5",
                   "Unresolved weak" = "#c8c8c8")

ph <- ggplot(density_df, aes(x = bin / 1e6, y = N, fill = track)) +
  geom_col(width = bin_size / 1e6 * 0.9) +
  facet_wrap(~track, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = track_colors) +
  scale_x_continuous(breaks = seq(0, chr_len_bp / 1e6, 5),
                     labels = function(x) paste0(x, " Mb"),
                     expand = c(0, 0)) +
  labs(title = paste0("h   Genome-position track (", chrom_label, ")"),
       x = chrom_label, y = "Variants per 100 kb bin") +
  theme_pub +
  theme(strip.text = element_text(face = "bold", size = 8),
        strip.placement = "outside")

# ============================================================================
# Composite layout with patchwork
# ============================================================================

# Title
title_grob <- textGrob(
  paste0("Rescuing biologically supported variants from Clair3-filtered calls\n",
         "Chromosome ", chrom_label, " (", sample_label, ")"),
  gp = gpar(fontsize = 14, fontface = "bold"),
  hjust = 0, x = 0.02
)

# Legend text
caption_text <- paste0(
  "Figure. Rescuing Clair3-filtered calls with biologically informed criteria on chromosome ",
  chrom_label, " (", sample_label, "). ",
  "(a) Retention funnel from raw candidates to strict PASS, filtered, rescued, and remaining. ",
  "(b) Stacked counts by final class. ",
  "(c) One-bar partition for ", chrom_label, ". ",
  "(d) Proportions. ",
  "(e) Rescued variants by class. ",
  "(f) Distributions of QUAL, depth (DP), and allele metrics. ",
  "(g) Weak indel read support. ",
  "(h) Genome-wide distribution of strict PASS, rescued, and unresolved weak calls along ",
  chrom_label, " (100 kb bins)."
)

caption_grob <- textGrob(
  str_wrap(caption_text, width = 140),
  gp = gpar(fontsize = 7, fontface = "plain"),
  hjust = 0, x = 0.02, vjust = 1
)

# Row 1: a + b
row1 <- pa + pb + plot_layout(widths = c(1, 1.5))

# Row 2: c + d + e + f
row2 <- pc + pd_plot + pe + pf + plot_layout(widths = c(0.7, 0.7, 1, 2.5))

# Row 3: g + h
row3 <- pg + ph + plot_layout(widths = c(1, 2.5))

# Full composite
composite <- (row1 / row2 / row3) +
  plot_layout(heights = c(1, 1, 1.2))

# Save
out_png <- file.path(outdir, "Figure_rescue_composite.png")
ggsave(out_png, composite, width = 18, height = 14, dpi = 300, bg = "white")
message("[FIG] Composite figure → ", out_png)

out_pdf <- file.path(outdir, "Figure_rescue_composite.pdf")
ggsave(out_pdf, composite, width = 18, height = 14, bg = "white")
message("[FIG] PDF → ", out_pdf)

# ============================================================================
# Individual panels (for flexibility)
# ============================================================================

panels <- list(a = pa, b = pb, c = pc, d = pd_plot, e = pe, f = pf, g = pg, h = ph)
for (nm in names(panels)) {
  fn <- file.path(outdir, paste0("panel_", nm, ".png"))
  ggsave(fn, panels[[nm]], width = 6, height = 5, dpi = 300, bg = "white")
}
message("[FIG] Individual panels saved")

# ============================================================================
# Summary tables for supplementary
# ============================================================================

# rescue performance table
rescue_perf <- data.frame(
  metric = c("Total raw", "Strict PASS", "Total filtered", "Rescued from filtered",
             "Rescued SNP (QUAL>=12)", "Rescued strong indel", "Rescued relaxed indel",
             "Pending weak", "Rescue rate (of filtered)"),
  count = c(n_raw, n_pass, n_filtered, n_rescued,
            n_strong_snp, n_strong_indel, n_relaxed_indel,
            n_pending, NA),
  percentage = c(100, round(100*n_pass/n_raw, 1), round(100*n_filtered/n_raw, 1),
                 round(100*n_rescued/n_raw, 1),
                 round(100*n_strong_snp/n_raw, 1),
                 round(100*n_strong_indel/n_raw, 1),
                 round(100*n_relaxed_indel/n_raw, 1),
                 round(100*n_pending/n_raw, 1),
                 round(100*n_rescued/n_filtered, 1)),
  stringsAsFactors = FALSE
)
fwrite(rescue_perf, file.path(outdir, "table_rescue_performance.tsv"), sep = "\t")

message("[DONE] All figures and tables written to: ", outdir)
