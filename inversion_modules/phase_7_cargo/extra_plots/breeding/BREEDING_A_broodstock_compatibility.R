#!/usr/bin/env Rscript

# =============================================================================
# BREEDING_A_broodstock_compatibility.R
#
# Hatchery-facing tool: given the inversion karyotype of every available
# broodstock fish, score either (a) a specific pair, or (b) all possible pairs,
# for inversion-haplotype diversity and recombination-block compatibility in
# expected F1 offspring.
#
# Why this matters: in a closed hatchery, picking parents that share inversion
# arrangements at every locus produces F1 with the same homozygous blocks as
# their parents — no new haplotype combinations possible inside inversions.
# Picking parents with complementary arrangements (one HOM_REF × one HOM_INV
# at the same locus → all F1 are HET) maximises haplotype diversity in F1
# but also maximises gamete loss in the F2 (because F1 are obligate
# heterokaryotypes there). The right answer depends on whether the breeder
# wants F1 production volume or long-term genetic diversity.
#
# This script reports both, so the breeder can pick.
#
# Modes:
#   --pair SAMPLE1 SAMPLE2   per-candidate karyotype grid + scores for one pair
#   --all                    all-vs-all matrix, top-N pairs by combined score
#
# Inputs:
#   ${SNAKE_CAND_FILE}                       — candidate intervals
#   ${SAMPLE_REGISTRY}/groups/inv_*.txt      — karyotype memberships
#   ${VARIANT_MASTER} (optional)             — for deleterious load
#   ${VESM_SCORES} (optional)                — for VESM-weighted load
#
# Outputs:
#   pair mode → ${EXTRAS_TBL_DIR}/BREEDING_A_pair__<S1>__x__<S2>.tsv
#                ${EXTRAS_TBL_DIR}/BREEDING_A_pair__<S1>__x__<S2>__summary.txt
#   all  mode → ${EXTRAS_TBL_DIR}/BREEDING_A_all_pairs_ranked.tsv
#                ${EXTRAS_TBL_DIR}/BREEDING_A_top_recommendations.tsv
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

SNAKE_CAND_FILE <- Sys.getenv("SNAKE_CAND_FILE")
SAMPLE_REGISTRY <- Sys.getenv("SAMPLE_REGISTRY")
EXTRAS_TBL_DIR  <- Sys.getenv("EXTRAS_TBL_DIR")
SAMPLE_LIST     <- Sys.getenv("SAMPLE_LIST")
TOPN_PAIRS      <- as.integer(Sys.getenv("BREEDING_TOPN_PAIRS", "50"))

stopifnot(nzchar(SNAKE_CAND_FILE), nzchar(SAMPLE_REGISTRY), nzchar(EXTRAS_TBL_DIR))

args <- commandArgs(trailingOnly = TRUE)
mode <- "all"
pair_s1 <- NULL; pair_s2 <- NULL
if (any(args == "--pair")) {
  i <- which(args == "--pair")
  mode <- "pair"
  if (i + 2 > length(args)) {
    stop("--pair requires two sample IDs: --pair SAMPLE1 SAMPLE2")
  }
  pair_s1 <- args[i + 1]; pair_s2 <- args[i + 2]
  stopifnot(nzchar(pair_s1), nzchar(pair_s2))
}

# ── Load candidates and per-sample karyotype matrix ──────────────────────────
cands <- fread(cmd = paste0("zcat ", shQuote(SNAKE_CAND_FILE)))

read_grp <- function(gid) {
  fp <- file.path(SAMPLE_REGISTRY, "groups", paste0(gid, ".txt"))
  if (!file.exists(fp)) return(character())
  trimws(readLines(fp))
}

# Build sample × candidate matrix of karyotypes (R / H / I / NA)
all_samples <- if (file.exists(SAMPLE_LIST)) {
  as.character(fread(SAMPLE_LIST, header = FALSE)[[1]])
} else {
  unique(unlist(lapply(cands$candidate_id, function(c)
    c(read_grp(paste0("inv_", c, "_HOM_REF")),
      read_grp(paste0("inv_", c, "_HET")),
      read_grp(paste0("inv_", c, "_HOM_INV"))))))
}
n_samples <- length(all_samples)
n_cands   <- nrow(cands)

K <- matrix(NA_character_, n_samples, n_cands,
            dimnames = list(all_samples, cands$candidate_id))
for (i in seq_len(n_cands)) {
  cid <- cands$candidate_id[i]
  K[read_grp(paste0("inv_", cid, "_HOM_REF")), cid] <- "R"
  K[read_grp(paste0("inv_", cid, "_HET")),     cid] <- "H"
  K[read_grp(paste0("inv_", cid, "_HOM_INV")), cid] <- "I"
}
covered <- apply(K, 1, function(r) sum(!is.na(r)))
if (sum(covered > 0) < 2) {
  message("[BREEDING_A] [skip] <2 samples have any karyotype calls"); quit(status = 0)
}
message("[BREEDING_A] ", sum(covered > 0), " samples × ", n_cands, " candidates")

# ── Per-pair-per-candidate cross outcome ─────────────────────────────────────
# Returns a single character: F1 expected karyotype distribution as a code.
#   R+R → R (all F1 homozygous standard)         "concord_R"
#   I+I → I (all F1 homozygous inverted)          "concord_I"
#   R+I → H (all F1 heterokaryotype)              "complementary"   — max F1 diversity but F2 will lose gametes
#   R+H → 50% R / 50% H                            "het_segregating"
#   I+H → 50% I / 50% H                            "het_segregating"
#   H+H → 25% R / 50% H / 25% I                    "het_x_het"      — ½ of F1 are HET, gamete loss in their offspring
cross_outcome <- function(k1, k2) {
  if (is.na(k1) || is.na(k2)) return(NA_character_)
  pair <- paste(sort(c(k1, k2)), collapse = "")
  switch(pair,
         "RR" = "concord_R",
         "II" = "concord_I",
         "IR" = "complementary",
         "HR" = "het_segregating",
         "HI" = "het_segregating",
         "HH" = "het_x_het",
         NA_character_)
}

# ── Score a pair ─────────────────────────────────────────────────────────────
# Two competing scores:
#   diversity_score  — high = F1 has more haplotype diversity
#                      complementary=2, het_segregating=1, het_x_het=1.5,
#                      concord=0
#   load_relief_score — high = F1 expected to have lower deleterious burden
#                       complementary=1 (always HET, masks recessives), others=0.5
# combined = diversity + load_relief
score_pair <- function(p1, p2) {
  outs <- vapply(seq_len(n_cands),
                 function(j) cross_outcome(K[p1, j], K[p2, j]),
                 character(1))
  div <- vapply(outs, function(o) {
    switch(o, "complementary" = 2, "het_x_het" = 1.5, "het_segregating" = 1,
              "concord_R" = 0, "concord_I" = 0, NA_real_)
  }, numeric(1))
  load_relief <- vapply(outs, function(o) {
    switch(o, "complementary" = 1, "het_segregating" = 0.5,
              "het_x_het" = 0.5, "concord_R" = 0, "concord_I" = 0,
              NA_real_)
  }, numeric(1))
  list(outcomes = outs,
       diversity = sum(div, na.rm = TRUE),
       load_relief = sum(load_relief, na.rm = TRUE),
       n_complementary = sum(outs == "complementary", na.rm = TRUE),
       n_het_x_het = sum(outs == "het_x_het", na.rm = TRUE),
       n_concord = sum(outs %in% c("concord_R", "concord_I"), na.rm = TRUE),
       n_evaluable = sum(!is.na(outs)))
}

# ── PAIR MODE ────────────────────────────────────────────────────────────────
if (mode == "pair") {
  if (!(pair_s1 %in% rownames(K)) || !(pair_s2 %in% rownames(K))) {
    stop("Sample(s) not found in karyotype matrix: ",
         setdiff(c(pair_s1, pair_s2), rownames(K)))
  }
  s <- score_pair(pair_s1, pair_s2)
  out <- data.table(
    candidate_id = cands$candidate_id,
    chrom = cands$chrom,
    length_mb = round((cands$end_bp - cands$start_bp) / 1e6, 3),
    karyotype_p1 = K[pair_s1, ],
    karyotype_p2 = K[pair_s2, ],
    cross_outcome = s$outcomes)
  out_path <- file.path(EXTRAS_TBL_DIR,
                          sprintf("BREEDING_A_pair__%s__x__%s.tsv",
                                   pair_s1, pair_s2))
  fwrite(out, out_path, sep = "\t")

  summary_path <- sub("\\.tsv$", "__summary.txt", out_path)
  writeLines(c(
    sprintf("Broodstock compatibility — %s × %s", pair_s1, pair_s2),
    sprintf("  Evaluable candidates:           %d / %d", s$n_evaluable, n_cands),
    sprintf("  Complementary loci:             %d  (HOM_REF × HOM_INV → all F1 = HET)",
            s$n_complementary),
    sprintf("  HET × HET loci:                 %d  (¼ R / ½ H / ¼ I in F1, F2 will lose gametes)",
            s$n_het_x_het),
    sprintf("  Concordant loci:                %d  (no diversity in F1)", s$n_concord),
    "",
    sprintf("  diversity_score:                %.1f", s$diversity),
    sprintf("  load_relief_score:              %.1f  (recessive masking via HET F1)",
            s$load_relief),
    sprintf("  combined_score:                 %.1f", s$diversity + s$load_relief)
  ), summary_path)
  message("[BREEDING_A pair] → ", out_path)
  message("[BREEDING_A pair] → ", summary_path)
  quit(status = 0)
}

# ── ALL-vs-ALL MODE ──────────────────────────────────────────────────────────
elig <- rownames(K)[covered > 0]
if (length(elig) > 200) {
  message("[BREEDING_A] WARNING: ", length(elig),
          " eligible samples → ", choose(length(elig), 2), " pairs.")
  message("              This may take several minutes.")
}

n_e <- length(elig)
pairs_out <- vector("list", choose(n_e, 2))
counter <- 0L
for (i in seq_len(n_e - 1)) {
  for (j in seq.int(i + 1, n_e)) {
    counter <- counter + 1L
    s <- score_pair(elig[i], elig[j])
    pairs_out[[counter]] <- data.table(
      sample_1 = elig[i], sample_2 = elig[j],
      n_evaluable = s$n_evaluable,
      n_complementary = s$n_complementary,
      n_het_x_het = s$n_het_x_het,
      n_concord = s$n_concord,
      diversity_score = s$diversity,
      load_relief_score = s$load_relief,
      combined_score = s$diversity + s$load_relief)
  }
}
all_pairs <- rbindlist(pairs_out)
all_pairs <- all_pairs[order(-combined_score)]
fwrite(all_pairs, file.path(EXTRAS_TBL_DIR, "BREEDING_A_all_pairs_ranked.tsv"),
       sep = "\t")

# Top recommendations: pick high-diversity pairs without reusing the same fish too often
# (greedy non-overlapping selection — produces a rough breeding-plan suggestion)
used <- character(0)
recs <- list()
for (k in seq_len(nrow(all_pairs))) {
  s1 <- all_pairs$sample_1[k]; s2 <- all_pairs$sample_2[k]
  if (s1 %in% used || s2 %in% used) next
  recs[[length(recs) + 1]] <- all_pairs[k]
  used <- c(used, s1, s2)
  if (length(recs) >= TOPN_PAIRS) break
}
recs_dt <- rbindlist(recs)
fwrite(recs_dt, file.path(EXTRAS_TBL_DIR, "BREEDING_A_top_recommendations.tsv"),
       sep = "\t")

message("[BREEDING_A all] ", nrow(all_pairs), " pairs ranked → BREEDING_A_all_pairs_ranked.tsv")
message("[BREEDING_A all] greedy top-", nrow(recs_dt),
        " non-overlapping pairs → BREEDING_A_top_recommendations.tsv")
