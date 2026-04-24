#!/usr/bin/env Rscript

# =============================================================================
# STEP_C01e_candidate_figures.R  (v8.4)
#
# PER-CANDIDATE DEEP-DIVE FIGURES for manuscript.
#
# For each Tier 1/2 candidate, generates a multi-panel figure with:
#   A: Chromosome overview ideogram with candidate highlighted
#   B: Regional PCA (samples colored by inferred karyotype REF/HET/INV)
#   C: Genotype heatmap (samples x SNPs, ordered by PC1)
#   D: Population-genomic signals (theta, het, Fst per karyotype group)
#   E: Breakpoint evidence (from DELLY/Manta if available)
#   F: Gene/repeat annotation context
#
# Inputs:
#   --scores <candidate_scores.tsv.gz>     -- from C01d
#   --triangles <triangle_dir>             -- from C01c (composition)
#   --precomp <precomp_dir>                -- from C01a (PCA data)
#   --samples <sample_list>                -- sample name mapping
#   [optional] --vcf_dir <dir>             -- Clair3 per-chr VCFs
#   [optional] --sv_dir <dir>              -- DELLY/Manta SV catalogs
#   [optional] --gff <annotation.gff3>     -- gene annotation
#   [optional] --repeats <repeats.bed>     -- repeat annotation
#   [optional] --het_dir <dir>             -- ACCEPTED BUT IGNORED (2026-04-17)
#   [optional] --repeats <file>             -- ACCEPTED BUT IGNORED (2026-04-17)
#
# Output per candidate:
#   <outdir>/candidate_<chr>_<start>_<end>/panel_A_ideogram.png
#   <outdir>/candidate_<chr>_<start>_<end>/panel_B_pca.png
#   <outdir>/candidate_<chr>_<start>_<end>/panel_I_systems.png   (multi-system only)
#   <outdir>/candidate_<chr>_<start>_<end>/data/
#     — sample_karyotypes.tsv, regional_pca.tsv, candidate_metadata.tsv
#     — per-sample group lists (*_samples.txt)
#     — executable extraction scripts (*_cmd.sh) for VCFtools / bcftools / ngsLD
#
# Note: panels C/D/E/F/G are emitted as shell-script extraction commands
# (in data/) rather than inline PNGs — the PNGs are intended to be built
# by those shell scripts in a follow-up step. Composite assembly
# (figure_composite.png) is done externally with Inkscape / figrid and
# is NOT produced by this script.
#
# Usage:
#   Rscript STEP_C01e_candidate_figures.R \
#     --scores candidates/candidate_scores.tsv.gz \
#     --triangles triangles/ \
#     --precomp precomp/ \
#     --samples sample_list.tsv \
#     --outdir figures/ \
#     [--tier_max 2] [--candidate_id 5]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Parse args
args <- commandArgs(trailingOnly = TRUE)
scores_file <- NULL; triangle_dir <- NULL; precomp_dir <- NULL
samples_file <- NULL; outdir <- "candidate_figures"
vcf_dir <- NULL; sv_dir <- NULL; gff_file <- NULL
coseg_dir <- NULL; regime_dir <- NULL
tier_max <- 2L; target_id <- NULL

# BUGFIX 2026-04-17: --repeats and --het_dir were parsed but never read
# anywhere in this script. They're accepted (silently ignored) for
# backward compatibility with existing launchers.
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--scores" && i < length(args))    { scores_file <- args[i+1]; i <- i+2 }
  else if (a == "--triangles" && i < length(args)) { triangle_dir <- args[i+1]; i <- i+2 }
  else if (a == "--precomp" && i < length(args))   { precomp_dir <- args[i+1]; i <- i+2 }
  else if (a == "--samples" && i < length(args))   { samples_file <- args[i+1]; i <- i+2 }
  else if (a == "--outdir" && i < length(args))    { outdir <- args[i+1]; i <- i+2 }
  else if (a == "--vcf_dir" && i < length(args))   { vcf_dir <- args[i+1]; i <- i+2 }
  else if (a == "--sv_dir" && i < length(args))    { sv_dir <- args[i+1]; i <- i+2 }
  else if (a == "--gff" && i < length(args))       { gff_file <- args[i+1]; i <- i+2 }
  else if (a == "--repeats" && i < length(args))   { i <- i+2 }  # accepted but unused
  else if (a == "--het_dir" && i < length(args))   { i <- i+2 }  # accepted but unused
  else if (a == "--coseg_dir" && i < length(args)) { coseg_dir <- args[i+1]; i <- i+2 }
  else if (a == "--regime_dir" && i < length(args)) { regime_dir <- args[i+1]; i <- i+2 }
  else if (a == "--tier_max" && i < length(args))  { tier_max <- as.integer(args[i+1]); i <- i+2 }
  else if (a == "--candidate_id" && i < length(args)) { target_id <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

if (is.null(scores_file)) stop("--scores required")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

DPI <- 400
THEME_BASE <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        plot.caption = element_text(size = 6, color = "grey60", hjust = 0))

KARYOTYPE_PAL <- c("REF" = "royalblue", "HET" = "grey50", "INV" = "red3")

# =============================================================================
# LOAD DATA
# =============================================================================

message("[C01e] Loading candidate scores...")
cand_dt <- fread(scores_file)
if (!is.null(target_id)) {
  cand_dt <- cand_dt[interval_id == target_id]
} else {
  cand_dt <- cand_dt[tier <= tier_max]
}
cand_dt <- cand_dt[order(tier, -final_score)]
message("[C01e] Candidates to plot: ", nrow(cand_dt))

# Composition
# BUGFIX 2026-04-17: The --triangles flag used to read
# triangle_sample_composition.tsv.gz, which was a pre-v9.3 C01c artifact.
# That script was retired when the v9.3 pipeline replaced the C01c
# triangle flow with the inv_detect v9.3 scoring_table track. The file
# doesn't exist on current runs, so comp_dt was always empty and
# Panel B fell through to "everyone is HET" (all samples grey).
#
# Fix: compute bands on demand from precomp$dt PC_1_<sample> columns
# using the same k-means(3) logic C01d uses at L651-680. The triangle
# file is still read if provided (legacy override) but is no longer
# required. compute_bands_for_candidate() returns a data.table with
# columns: sample, band, karyotype, pc1_avg.
compute_bands_for_candidate <- function(pc, start_mb, end_mb,
                                        real_names = NULL, min_samples = 20L) {
  if (is.null(pc) || is.null(pc$dt)) return(data.table())
  dt_chr <- pc$dt
  s_bp <- start_mb * 1e6; e_bp <- end_mb * 1e6

  # Windows fully inside the candidate
  inner_w <- which(dt_chr$start_bp >= s_bp & dt_chr$end_bp <= e_bp)
  if (length(inner_w) < 3) return(data.table())

  pc1_cols <- grep("^PC_1_", names(dt_chr), value = TRUE)
  if (length(pc1_cols) < min_samples) return(data.table())

  avg_pc1 <- colMeans(as.matrix(dt_chr[inner_w, ..pc1_cols]), na.rm = TRUE)
  valid <- is.finite(avg_pc1)
  if (sum(valid) < min_samples) return(data.table())

  vals <- avg_pc1[valid]
  km <- tryCatch(kmeans(vals, centers = 3, nstart = 10),
                 error = function(e) NULL)
  if (is.null(km)) return(data.table())

  co <- order(km$centers[, 1])  # low PC1 → band 1 (REF), high PC1 → band 3 (INV)
  bands <- integer(length(vals))
  for (b in 1:3) bands[km$cluster == co[b]] <- b

  snames <- sub("^PC_1_", "", names(vals))

  # Optional: map IndN → real names if samples_file was provided
  if (!is.null(real_names) && length(real_names) > 0 &&
      length(snames) >= 1 && grepl("^Ind[0-9]", snames[1])) {
    ind_to_real <- setNames(real_names, paste0("Ind", seq_along(real_names) - 1))
    mapped <- ind_to_real[snames]
    snames[!is.na(mapped)] <- mapped[!is.na(mapped)]
  }

  data.table(
    sample    = snames,
    band      = paste0("band", bands),
    karyotype = fifelse(bands == 1L, "REF",
                fifelse(bands == 2L, "HET", "INV")),
    pc1_avg   = round(as.numeric(vals), 4)
  )
}

# Legacy composition file (still accepted if provided, not required)
comp_dt <- data.table()
if (!is.null(triangle_dir)) {
  f <- file.path(triangle_dir, "triangle_sample_composition.tsv.gz")
  if (file.exists(f)) {
    comp_dt <- fread(f)
    message("[C01e] Legacy composition file loaded (", nrow(comp_dt), " rows)")
  } else {
    message("[C01e] No legacy composition file at ", f,
            " — will compute bands from precomp on the fly")
  }
}

# Sample name mapping
name_map <- NULL; real_names <- NULL
if (!is.null(samples_file) && file.exists(samples_file)) {
  slist <- as.character(fread(samples_file, header = FALSE)[[1]])
  real_names <- slist
}

# Precomp (lazy load per chr)
precomp_cache <- list()
load_precomp <- function(chr) {
  if (!is.null(precomp_cache[[chr]])) return(precomp_cache[[chr]])
  if (is.null(precomp_dir)) return(NULL)
  f <- list.files(precomp_dir, pattern = paste0(chr, "\\.precomp\\.rds$"), full.names = TRUE)
  if (length(f) == 0) return(NULL)
  obj <- readRDS(f[1])
  precomp_cache[[chr]] <<- obj
  obj
}

# =============================================================================
# PER-CANDIDATE FIGURE GENERATION
# =============================================================================

for (ci in seq_len(nrow(cand_dt))) {
  cand <- cand_dt[ci]
  chr <- cand$chrom
  start_mb <- cand$start_mb; end_mb <- cand$end_mb
  span_mb <- end_mb - start_mb

  cand_id <- paste0(chr, "_", round(start_mb, 1), "_", round(end_mb, 1))
  cand_dir <- file.path(outdir, paste0("candidate_", cand_id))
  dir.create(cand_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cand_dir, "data"), showWarnings = FALSE)

  message("\n[C01e] === Candidate ", ci, "/", nrow(cand_dt), ": ",
          cand_id, " (Tier ", cand$tier, ", ", cand$pattern, ") ===")

  # --- Get composition (band assignments) ---
  # Prefer legacy file if it has this candidate; otherwise compute from precomp.
  # BUGFIX 2026-04-17: previously only legacy file was consulted, so
  # Panel B on current runs always fell through to "everyone HET" (grey).
  iv_comp <- data.table()
  if (nrow(comp_dt) > 0) {
    iv_comp <- comp_dt[chrom == chr & interval_id == cand$interval_id]
  }
  if (nrow(iv_comp) == 0) {
    pc_for_bands <- load_precomp(chr)
    iv_comp <- compute_bands_for_candidate(pc_for_bands, start_mb, end_mb,
                                           real_names = real_names)
    if (nrow(iv_comp) > 0) {
      iv_comp[, `:=`(chrom = chr, interval_id = cand$interval_id)]
      message("  [bands] computed on the fly: ",
              sum(iv_comp$karyotype == "REF"), " REF / ",
              sum(iv_comp$karyotype == "HET"), " HET / ",
              sum(iv_comp$karyotype == "INV"), " INV")
    }
  }
  if (nrow(iv_comp) > 0) {
    # Ensure karyotype column is present (legacy file uses `band` column)
    if (!"karyotype" %in% names(iv_comp) && "band" %in% names(iv_comp)) {
      iv_comp[, karyotype := fifelse(band == "band1", "REF",
                             fifelse(band == "band2", "HET", "INV"))]
    }
    fwrite(iv_comp, file.path(cand_dir, "data", "sample_karyotypes.tsv"), sep = "\t")
  }

  # --- Panel B: Regional PCA ---
  pc <- load_precomp(chr)
  pB <- NULL
  if (!is.null(pc) && nrow(iv_comp) > 0) {
    dt <- pc$dt
    # Get PC1 columns for the windows in this interval
    start_bp <- start_mb * 1e6; end_bp <- end_mb * 1e6
    win_idx <- which(dt$start_bp >= start_bp & dt$end_bp <= end_bp)
    if (length(win_idx) >= 5) {
      pc1_cols <- grep("^PC_1_", names(dt), value = TRUE)
      pc2_cols <- grep("^PC_2_", names(dt), value = TRUE)

      if (length(pc1_cols) > 0 && length(pc2_cols) > 0) {
        pc1_vals <- colMeans(as.matrix(dt[win_idx, ..pc1_cols]), na.rm = TRUE)
        pc2_vals <- colMeans(as.matrix(dt[win_idx, ..pc2_cols]), na.rm = TRUE)

        snames <- sub("^PC_1_", "", pc1_cols)
        # Map to real names
        if (!is.null(real_names) && length(real_names) == length(snames) &&
            grepl("^Ind[0-9]", snames[1])) {
          snames <- real_names
        }

        pca_dt <- data.table(sample = snames, PC1 = pc1_vals, PC2 = pc2_vals)
        # Merge karyotype
        pca_dt <- merge(pca_dt, iv_comp[, .(sample, karyotype)],
                         by = "sample", all.x = TRUE)
        pca_dt[is.na(karyotype), karyotype := "HET"]

        # Band counts
        n_ref <- sum(pca_dt$karyotype == "REF")
        n_het <- sum(pca_dt$karyotype == "HET")
        n_inv <- sum(pca_dt$karyotype == "INV")

        pB <- ggplot(pca_dt, aes(x = PC1, y = PC2, color = karyotype)) +
          geom_point(size = 2, alpha = 0.7) +
          scale_color_manual(values = KARYOTYPE_PAL,
                             labels = c(paste0("REF (n=", n_ref, ")"),
                                       paste0("HET (n=", n_het, ")"),
                                       paste0("INV (n=", n_inv, ")")),
                             name = "Inferred\nkaryotype") +
          labs(title = "Regional PCA (samples)",
               subtitle = "Three-band pattern consistent\nwith a segregating inversion",
               x = "PC1", y = "PC2") +
          THEME_BASE

        tryCatch(ggsave(file.path(cand_dir, "panel_B_pca.png"),
                         pB, width = 6, height = 5, dpi = DPI),
                 error = function(e) message("  [B] ", e$message))

        fwrite(pca_dt, file.path(cand_dir, "data", "regional_pca.tsv"), sep = "\t")
      }
    }
  }

  # --- Panel A: Chromosome overview ---
  # Simple ideogram with candidate highlighted
  chr_len <- if (!is.null(pc)) max(pc$dt$end_bp) / 1e6 else end_mb * 1.2

  pA <- ggplot() +
    # Chromosome backbone
    annotate("segment", x = 0, xend = chr_len, y = 0, yend = 0,
             color = "grey70", linewidth = 3) +
    # Candidate region
    annotate("rect", xmin = start_mb, xmax = end_mb, ymin = -0.3, ymax = 0.3,
             fill = "red3", alpha = 0.3, color = "red3", linewidth = 0.5) +
    # Label
    annotate("text", x = (start_mb + end_mb) / 2, y = 0.5,
             label = "Candidate inversion", color = "red3", size = 3, fontface = "italic") +
    # Zoom lines
    annotate("segment", x = start_mb, xend = start_mb, y = 0.3, yend = 0.7,
             color = "red3", linewidth = 0.3, linetype = "dashed") +
    annotate("segment", x = end_mb, xend = end_mb, y = 0.3, yend = 0.7,
             color = "red3", linewidth = 0.3, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, chr_len, by = 10)) +
    coord_cartesian(ylim = c(-0.8, 0.8)) +
    labs(title = paste0(chr, " overview"),
         x = paste0(chr, " (Mb)"), y = NULL) +
    THEME_BASE +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

  tryCatch(ggsave(file.path(cand_dir, "panel_A_ideogram.png"),
                   pA, width = 10, height = 2, dpi = DPI),
           error = function(e) message("  [A] ", e$message))

  # --- Panel C: Genotype heatmap placeholder ---
  # This needs VCF data -- create a template that can be filled
  # when VCF directory is provided
  if (!is.null(vcf_dir)) {
    message("  [C] VCF genotype heatmap: TODO (requires bcftools query)")
    # Would need: bcftools query -r chr:start-end -f '%POS\t[%GT\t]\n' vcf
    # For now, save the command
    cmd <- paste0("bcftools query -r ", chr, ":", round(start_mb*1e6), "-", round(end_mb*1e6),
                  " -f '%POS\\t[%GT\\t]\\n' ", vcf_dir, "/", chr, ".vcf.gz")
    writeLines(cmd, file.path(cand_dir, "data", "genotype_extraction_cmd.sh"))
  }

  # --- Panel D: Fst between groups across the region ---
  # Uses C01i core package sample classifications
  # Computes per-window Fst(HOMO_REF vs HOMO_INV) from VCF genotypes
  if (!is.null(vcf_dir)) {
    coseg_markers <- NULL; coseg_samples <- NULL
    if (!is.null(coseg_dir)) {
      # Try to load C01i outputs for this candidate
      mf <- file.path(coseg_dir, "for_C01h",
                       paste0(chr, "_I", cand$interval_id, "_S1_markers.tsv.gz"))
      sf <- file.path(coseg_dir, "for_C01h",
                       paste0(chr, "_I", cand$interval_id, "_S1_samples.tsv.gz"))
      if (file.exists(mf)) coseg_markers <- fread(mf)
      if (file.exists(sf)) coseg_samples <- fread(sf)
    }

    # Generate Fst extraction command (per-window Fst between groups)
    if (!is.null(coseg_samples)) {
      ref_list <- coseg_samples[group == "HOMO_REF"]$sample
      inv_list <- coseg_samples[group %in% c("HOMO_INV", "RARE_INV")]$sample
      het_list <- coseg_samples[group %in% c("HET", "RARE_HET")]$sample

      # Write sample lists for ANGSD or vcftools Fst
      writeLines(ref_list, file.path(cand_dir, "data", "ref_samples.txt"))
      writeLines(inv_list, file.path(cand_dir, "data", "inv_samples.txt"))
      if (length(het_list) > 0)
        writeLines(het_list, file.path(cand_dir, "data", "het_samples.txt"))

      # VCFtools Fst command
      vcf_f <- file.path(vcf_dir, paste0(chr, ".vcf.gz"))
      cmd_fst <- paste0(
        "vcftools --gzvcf ", vcf_f, " \\\n",
        "  --weir-fst-pop ", file.path(cand_dir, "data", "ref_samples.txt"), " \\\n",
        "  --weir-fst-pop ", file.path(cand_dir, "data", "inv_samples.txt"), " \\\n",
        "  --fst-window-size 50000 --fst-window-step 10000 \\\n",
        "  --chr ", chr, " \\\n",
        "  --from-bp ", round(start_mb * 1e6), " --to-bp ", round(end_mb * 1e6), " \\\n",
        "  --out ", file.path(cand_dir, "data", "fst_ref_vs_inv"))
      writeLines(cmd_fst, file.path(cand_dir, "data", "fst_cmd.sh"))
      message("  [D] Fst command saved (", length(ref_list), " REF vs ",
              length(inv_list), " INV)")

      # Per-group het extraction command
      for (grp_name in c("ref", "inv", "het")) {
        grp_file <- file.path(cand_dir, "data", paste0(grp_name, "_samples.txt"))
        if (file.exists(grp_file)) {
          cmd_het <- paste0(
            "vcftools --gzvcf ", vcf_f, " \\\n",
            "  --keep ", grp_file, " \\\n",
            "  --het \\\n",
            "  --chr ", chr, " \\\n",
            "  --from-bp ", round(start_mb * 1e6), " --to-bp ", round(end_mb * 1e6), " \\\n",
            "  --out ", file.path(cand_dir, "data", paste0("het_", grp_name)))
          writeLines(cmd_het, file.path(cand_dir, "data", paste0("het_", grp_name, "_cmd.sh")))
        }
      }
    } else {
      message("  [D] No C01i sample classification -- run C01i first")
    }

    # Genotype heatmap extraction command
    cmd_gt <- paste0(
      "bcftools query -r ", chr, ":", round(start_mb*1e6), "-", round(end_mb*1e6),
      " -f '%POS\\t[%GT\\t]\\n' ", file.path(vcf_dir, paste0(chr, ".vcf.gz")),
      " > ", file.path(cand_dir, "data", "genotypes.tsv"))
    writeLines(cmd_gt, file.path(cand_dir, "data", "genotype_extraction_cmd.sh"))
  }

  # --- Panel E: Breakpoint evidence from DELLY/Manta ---
  if (!is.null(sv_dir)) {
    # Extract INV calls overlapping the candidate region
    for (sv_type in c("INV", "BND")) {
      sv_f <- file.path(sv_dir, paste0("delly_", sv_type, "_final.vcf.gz"))
      if (!file.exists(sv_f)) sv_f <- file.path(sv_dir, paste0("manta_", sv_type, ".vcf.gz"))
      if (file.exists(sv_f)) {
        cmd_sv <- paste0("bcftools view -r ", chr, ":", round((start_mb-0.2)*1e6), "-",
                          round((end_mb+0.2)*1e6), " ", sv_f,
                          " > ", file.path(cand_dir, "data", paste0("sv_", sv_type, ".vcf")))
        writeLines(cmd_sv, file.path(cand_dir, "data", paste0("sv_", sv_type, "_cmd.sh")))
      }
    }
    message("  [E] SV extraction commands saved")
  }

  # --- Panel F: LD confirmation (Arctic cod check) ---
  # LD in all samples should show a block. LD in HOMO_REF only should NOT.
  if (!is.null(vcf_dir) && !is.null(coseg_samples)) {
    # BUGFIX 2026-04-17: gt_data was never constructed anywhere in this
    # script (`%||%` only saves NULLs, not missing bindings). Safe default
    # is the total sample count from coseg_samples (includes all groups).
    n_total_ind <- if (!is.null(coseg_samples) && "sample" %in% names(coseg_samples))
                     length(unique(coseg_samples$sample)) else 226L
    # ref_list may not have been defined if the Panel D branch above
    # didn't enter its HOMO_REF-populating sub-branch. Guard it.
    n_ref_ind <- if (exists("ref_list")) length(ref_list) else 0L

    # ngsLD commands for all-sample vs ref-only
    cmd_ld_all <- paste0(
      "# LD in ALL samples (should show inversion block)\n",
      "ngsLD --geno <beagle_file> --pos <pos_file> --probs --n_ind ", n_total_ind,
      " --max_kb_dist ", round((end_mb - start_mb) * 1000 + 500),
      " --out ", file.path(cand_dir, "data", "ld_all.tsv"))
    cmd_ld_ref <- paste0(
      "# LD in HOMO_REF only (block should disappear)\n",
      "# Subset beagle to ref_samples.txt first\n",
      "ngsLD --geno <beagle_ref_only> --pos <pos_file> --probs --n_ind ", n_ref_ind,
      " --max_kb_dist ", round((end_mb - start_mb) * 1000 + 500),
      " --out ", file.path(cand_dir, "data", "ld_ref_only.tsv"))
    writeLines(c(cmd_ld_all, "", cmd_ld_ref),
               file.path(cand_dir, "data", "ld_confirmation_cmd.sh"))
    message("  [F] LD confirmation commands saved")
  }

  # --- Panel G: Gene annotation context ---
  if (!is.null(gff_file) && file.exists(gff_file)) {
    cmd_gff <- paste0(
      "awk '$1==\"", chr, "\" && $4>=", round(start_mb*1e6),
      " && $5<=", round(end_mb*1e6), "' ", gff_file,
      " > ", file.path(cand_dir, "data", "genes_in_region.gff3"))
    writeLines(cmd_gff, file.path(cand_dir, "data", "gene_extraction_cmd.sh"))
    message("  [G] Gene extraction command saved")
  }

  # --- Panel H: Regime profile from C01j ---
  if (!is.null(regime_dir)) {
    regime_f <- file.path(regime_dir, "regime_windows.tsv.gz")
    if (file.exists(regime_f)) {
      reg_dt <- fread(regime_f)
      reg_cand <- reg_dt[chrom == chr & interval_id == cand$interval_id]
      if (nrow(reg_cand) > 0) {
        fwrite(reg_cand, file.path(cand_dir, "data", "regime_profile.tsv"), sep = "\t")
        message("  [H] Regime profile: ", nrow(reg_cand), " windows")
      }
    }
  }

  # --- Panel I: Multi-system confirmation (Ips typographus style) ---
  # For candidates with multiple overlapping inversion systems (from C01i),
  # generate: system span bars + per-system-pair Fst + LD overlay
  if (!is.null(coseg_dir)) {
    sys_f <- file.path(coseg_dir, "marker_coseg_blocks.tsv.gz")
    samp_f <- file.path(coseg_dir, "marker_coseg_samples.tsv.gz")
    if (file.exists(sys_f) && file.exists(samp_f)) {
      sys_dt <- fread(sys_f)[chrom == chr & interval_id == cand$interval_id]
      samp_all <- fread(samp_f)[chrom == chr & interval_id == cand$interval_id]

      if (nrow(sys_dt) > 0) {
        n_sys <- nrow(sys_dt)
        message("  [I] Multi-system: ", n_sys, " systems detected")

        # System color palette
        sys_colors <- c("steelblue", "darkorange", "green4", "red3",
                         "gold3", "purple3", "brown", "cyan4",
                         "pink3", "grey40")

        # --- I.1: System span bar plot ---
        sys_dt[, sys_label := paste0("S", system_id)]
        sys_dt[, color := sys_colors[pmin(system_id, length(sys_colors))]]
        sys_dt[, start_mb_sys := span_start / 1e6]
        sys_dt[, end_mb_sys := span_end / 1e6]

        pI1 <- ggplot(sys_dt, aes(xmin = start_mb_sys, xmax = end_mb_sys,
                                   ymin = system_id - 0.4, ymax = system_id + 0.4,
                                   fill = sys_label)) +
          geom_rect(alpha = 0.7, color = "grey30", linewidth = 0.3) +
          geom_text(aes(x = (start_mb_sys + end_mb_sys) / 2, y = system_id,
                        label = paste0(sys_label, " (", round(frequency * 100), "%)")),
                    size = 2.5, fontface = "bold") +
          scale_fill_manual(values = setNames(sys_dt$color, sys_dt$sys_label),
                             guide = "none") +
          labs(title = paste0(chr, " I", cand$interval_id, " -- Inversion Systems"),
               subtitle = paste0(n_sys, " co-segregating blocks"),
               x = paste0(chr, " (Mb)"), y = "System") +
          theme_minimal(base_size = 9) +
          theme(axis.text.y = element_blank())

        tryCatch(ggsave(file.path(cand_dir, "panel_I_systems.png"),
                         pI1, width = 12, height = max(2, n_sys * 0.6 + 1), dpi = DPI),
                 error = function(e) message("  [I1] ", e$message))

        # --- I.2: Per-system-pair Fst commands ---
        # For each pair of systems, write sample lists and Fst command
        # Also: within each system, write MJ-MN, MJ-R, MN-R pairs
        vcf_f <- file.path(vcf_dir, paste0(chr, ".vcf.gz"))

        fst_cmds <- character(0)
        for (si in seq_len(n_sys)) {
          sys_samp <- samp_all[system_id == si]
          ref_s <- sys_samp[group == "HOMO_REF"]$sample
          inv_s <- sys_samp[group %in% c("HOMO_INV", "RARE_INV")]$sample
          het_s <- sys_samp[group %in% c("HET", "RARE_HET")]$sample

          # Write per-system sample lists
          prefix <- paste0("S", si)
          if (length(ref_s) >= 3)
            writeLines(ref_s, file.path(cand_dir, "data", paste0(prefix, "_ref.txt")))
          if (length(inv_s) >= 2)
            writeLines(inv_s, file.path(cand_dir, "data", paste0(prefix, "_inv.txt")))
          if (length(het_s) >= 2)
            writeLines(het_s, file.path(cand_dir, "data", paste0(prefix, "_het.txt")))

          # Fst: REF vs INV (main signal)
          if (length(ref_s) >= 3 && length(inv_s) >= 2) {
            sys_start <- sys_dt[system_id == si]$span_start
            sys_end <- sys_dt[system_id == si]$span_end
            flank <- (sys_end - sys_start) * 0.3

            cmd <- paste0(
              "# System ", si, " (", sys_dt[system_id == si]$sys_label,
              "): REF vs INV, freq=", round(sys_dt[system_id == si]$frequency * 100), "%\n",
              "vcftools --gzvcf ", vcf_f, " \\\n",
              "  --weir-fst-pop ", file.path(cand_dir, "data", paste0(prefix, "_ref.txt")), " \\\n",
              "  --weir-fst-pop ", file.path(cand_dir, "data", paste0(prefix, "_inv.txt")), " \\\n",
              "  --fst-window-size 50000 --fst-window-step 10000 \\\n",
              "  --chr ", chr, " \\\n",
              "  --from-bp ", max(1, sys_start - flank),
              " --to-bp ", sys_end + flank, " \\\n",
              "  --out ", file.path(cand_dir, "data", paste0("fst_", prefix, "_ref_vs_inv")))
            fst_cmds <- c(fst_cmds, cmd)
          }

          # Fst: MJ-R and MN-R (recombinant confirmation)
          # Find recombinant samples from C01h if available
          recomb_f <- file.path(coseg_dir, "for_C01h",
                                 paste0(chr, "_I", cand$interval_id, "_S", si, "_het_list.txt"))
          if (file.exists(recomb_f)) {
            recomb_samples <- readLines(recomb_f)
            if (length(recomb_samples) >= 2) {
              writeLines(recomb_samples,
                         file.path(cand_dir, "data", paste0(prefix, "_recomb.txt")))
              # MJ-R Fst
              if (length(ref_s) >= 3) {
                cmd_mjr <- paste0(
                  "# System ", si, ": REF vs RECOMBINANT (MJ-R)\n",
                  "vcftools --gzvcf ", vcf_f, " \\\n",
                  "  --weir-fst-pop ", file.path(cand_dir, "data", paste0(prefix, "_ref.txt")), " \\\n",
                  "  --weir-fst-pop ", file.path(cand_dir, "data", paste0(prefix, "_recomb.txt")), " \\\n",
                  "  --fst-window-size 50000 --fst-window-step 10000 \\\n",
                  "  --chr ", chr, " \\\n",
                  "  --from-bp ", max(1, sys_start - flank),
                  " --to-bp ", sys_end + flank, " \\\n",
                  "  --out ", file.path(cand_dir, "data", paste0("fst_", prefix, "_ref_vs_recomb")))
                fst_cmds <- c(fst_cmds, cmd_mjr)
              }
              # MN-R Fst
              if (length(inv_s) >= 2) {
                cmd_mnr <- paste0(
                  "# System ", si, ": INV vs RECOMBINANT (MN-R)\n",
                  "vcftools --gzvcf ", vcf_f, " \\\n",
                  "  --weir-fst-pop ", file.path(cand_dir, "data", paste0(prefix, "_inv.txt")), " \\\n",
                  "  --weir-fst-pop ", file.path(cand_dir, "data", paste0(prefix, "_recomb.txt")), " \\\n",
                  "  --fst-window-size 50000 --fst-window-step 10000 \\\n",
                  "  --chr ", chr, " \\\n",
                  "  --from-bp ", max(1, sys_start - flank),
                  " --to-bp ", sys_end + flank, " \\\n",
                  "  --out ", file.path(cand_dir, "data", paste0("fst_", prefix, "_inv_vs_recomb")))
                fst_cmds <- c(fst_cmds, cmd_mnr)
              }
            }
          }
        }

        # Cross-system Fst (between systems)
        if (n_sys >= 2) {
          for (s1 in seq_len(n_sys - 1)) {
            for (s2 in (s1 + 1):n_sys) {
              inv1 <- samp_all[system_id == s1 & group %in% c("HOMO_INV", "RARE_INV")]$sample
              inv2 <- samp_all[system_id == s2 & group %in% c("HOMO_INV", "RARE_INV")]$sample
              if (length(inv1) >= 2 && length(inv2) >= 2) {
                writeLines(inv1, file.path(cand_dir, "data", paste0("S", s1, "_inv_cross.txt")))
                writeLines(inv2, file.path(cand_dir, "data", paste0("S", s2, "_inv_cross.txt")))
                cmd_cross <- paste0(
                  "# Cross-system: S", s1, " INV vs S", s2, " INV\n",
                  "vcftools --gzvcf ", vcf_f, " \\\n",
                  "  --weir-fst-pop ", file.path(cand_dir, "data", paste0("S", s1, "_inv_cross.txt")), " \\\n",
                  "  --weir-fst-pop ", file.path(cand_dir, "data", paste0("S", s2, "_inv_cross.txt")), " \\\n",
                  "  --fst-window-size 50000 --fst-window-step 10000 \\\n",
                  "  --chr ", chr, " \\\n",
                  "  --out ", file.path(cand_dir, "data", paste0("fst_S", s1, "_vs_S", s2)))
                fst_cmds <- c(fst_cmds, cmd_cross)
              }
            }
          }
        }

        # Write all Fst commands to one script
        writeLines(c("#!/bin/bash", "set -euo pipefail", "",
                      paste0("# Multi-system Fst commands for ", chr, " I", cand$interval_id),
                      paste0("# ", n_sys, " systems, colored by: ",
                             paste(sys_dt$sys_label, sys_dt$color, sep = "=", collapse = ", ")),
                      "", fst_cmds),
                   file.path(cand_dir, "data", "fst_multi_system.sh"))
        message("  [I2] Multi-system Fst: ", length(fst_cmds), " commands (",
                n_sys, " within-system + cross-system pairs)")

        # Save system metadata for plotting
        fwrite(sys_dt, file.path(cand_dir, "data", "system_spans.tsv"), sep = "\t")
      }
    }
  }

  # --- Save candidate metadata ---
  # BUGFIX 2026-04-17: column name set aligned with current v9.3 C01d output.
  # The old list (d1_triangle..d10_ghsl, tube_stage_d/e, band1_n..band3_n,
  # band_symmetry, interval_type) was from the pre-v9.3 C01c triangle
  # pipeline — none of those columns exist in current candidate_scores.tsv.gz.
  # intersect() silently discarded all of them, leaving candidate_metadata.tsv
  # with only the basic fields (chrom, interval_id, start_mb, end_mb,
  # span_mb, final_score, tier, pattern, hyp_verdict).
  meta_cols <- intersect(names(cand), c(
    # Identity + position (unchanged from v8.4)
    "chrom", "interval_id", "start_mb", "end_mb", "span_mb",
    # Summary scores
    "final_score", "dim_positive", "tier", "pattern",
    # v9.3 scoring dimensions (12)
    "d1_block_strength", "d2_block_shape", "d3_nn_persistence",
    "d4_decay_flatness", "d5_interior_quality", "d6_consensus",
    "d7_sv_breakpoint", "d8_peel_or_hyp", "d9_pca_clusters",
    "d10_partition", "d11_boundary_concordance", "d12_snake_concordance",
    # Source evidence + classifications
    "shape_class", "landscape_category", "hyp_verdict",
    "l1b_peel", "l2_peel", "survives_nn40", "nn_birth",
    "n_variants", "n_sv_hits",
    # v9.3.2 Cheat 25 viability
    "cheat25_status",
    # v9.3.2 boundary / snake concordance annotations
    "snake_overlap", "boundary_verdict_left", "boundary_verdict_right",
    "n_boundary_cheats",
    # v9.3.4 popgen annotations (Engine B)
    "popgen_fst_b1b3", "popgen_theta_pi", "popgen_theta_W",
    "popgen_Tajima_D", "popgen_method",
    # v9.3.4 family Fst ratio (test_05)
    "cheat5_family_fst_ratio",
    # v9 morphology PA summary
    "pa_flat_inv_mean", "pa_spiky_inv_mean", "pa_frag_mean",
    "pa_family_like_mean", "pa_jaggedness_mean",
    "pa_block_compactness_mean"
  ))
  fwrite(cand[, ..meta_cols], file.path(cand_dir, "data", "candidate_metadata.tsv"), sep = "\t")

  message("  Output: ", cand_dir)
}

message("\n[DONE] Generated ", nrow(cand_dt), " candidate figure directories -> ", outdir)
message("Run individual panel scripts or use Inkscape to assemble composites.")
