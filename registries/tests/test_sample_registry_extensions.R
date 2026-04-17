#!/usr/bin/env Rscript
# =============================================================================
# test_sample_registry_extensions.R — chat 11 test fixture
# =============================================================================
# Exercises the 10 new sample-registry methods added in chat 11:
#   get_sample_metadata, get_sample_groups, get_family, list_families,
#   get_family_members, list_carriers_across_candidates,
#   compute_group_overlap, list_recombinant_candidates,
#   get_groups_for_candidate, find_co_segregating_groups
#
# Plus the 3 composite query methods:
#   query$nested_family_tree, query$overlap_clusters, query$sample_inversion_load
#
# Synthetic cohort: 10 samples (CGA001–CGA010) in 3 families (F1, F2, F3).
# Synthetic inversions:
#   LG12_A polymorphic (HOM_REF/HET/HOM_INV/RECOMBINANT registered)
#   LG12_B nested in A — shares HOM_INV membership with A
#   LG17_X independent — one sample is RECOMBINANT
#
# Run from repo root:
#   Rscript registries/tests/test_sample_registry_extensions.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Resolve paths relative to this script
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0) script_path <- "registries/tests/test_sample_registry_extensions.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))

loader <- c(
  file.path(script_dir, "..", "api", "R", "registry_loader.R"),
  "registries/api/R/registry_loader.R"
)
loader <- loader[file.exists(loader)][1]
if (is.na(loader)) stop("Could not find registry_loader.R")

# Build an isolated cwd with a utils/sample_registry.R shim so the samples
# API hits the "real" (not shim) code path. Uses the fixed file from the
# repo if available, else the candidate search in load_samples_api will try.
sample_reg <- c(
  file.path(script_dir, "..", "..", "utils", "sample_registry.R"),
  "utils/sample_registry.R"
)
sample_reg <- sample_reg[file.exists(sample_reg)][1]

tmp_cwd <- tempfile("samp_ext_test_")
dir.create(tmp_cwd, recursive = TRUE)
if (!is.na(sample_reg)) {
  dir.create(file.path(tmp_cwd, "utils"), showWarnings = FALSE)
  file.copy(sample_reg, file.path(tmp_cwd, "utils/sample_registry.R"))
}
setwd(tmp_cwd)
source(loader)

# Seed master BEFORE load_registry() so sample_registry.R picks it up
root <- file.path(tmp_cwd, "registries")
sample_dir <- file.path(root, "data", "sample_registry")
dir.create(sample_dir, recursive = TRUE)
master_dt <- data.table(
  sample_id    = paste0("CGA", sprintf("%03d", 1:10)),
  ind_id_226   = paste0("Ind", 0:9),
  sample_index = 0:9,
  family_id    = c("F1","F1","F1","F2","F2","F2","F2","F3","F3", NA_character_)
)
fwrite(master_dt, file.path(sample_dir, "sample_master.tsv"), sep = "\t")

reg <- load_registry(root, create_if_missing = TRUE)
sm  <- reg$samples
iv  <- reg$intervals

assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg)
  cat("  ok: ", msg, "\n", sep = "")
}

cat("──── Seed: candidates + groups ────\n")
# Intervals: nested + independent
iv$add_candidate("LG12_A", "chr12", 1e6, 10e6)
iv$add_candidate("LG12_B", "chr12", 3e6,  6e6, parent_id = "LG12_A")
iv$add_candidate("LG17_X", "chr17", 1e6,  5e6)
iv$add_candidate("LG12_E", "chr12", 20e6, 25e6)  # standalone, no groups

# Groups mirroring C01i_d_seal.R conventions (inv_<cid>_<STATUS>)
sm$add_group("inv_LG12_A_HOM_REF",     c("CGA001","CGA002","CGA003","CGA004"),
             chrom="chr12", inv_id="LG12_A", subgroup="HOM_REF")
sm$add_group("inv_LG12_A_HET",         c("CGA005","CGA006"),
             chrom="chr12", inv_id="LG12_A", subgroup="HET")
sm$add_group("inv_LG12_A_HOM_INV",     c("CGA007","CGA008","CGA009"),
             chrom="chr12", inv_id="LG12_A", subgroup="HOM_INV")
sm$add_group("inv_LG12_A_RECOMBINANT", c("CGA010","CGA006"),
             chrom="chr12", inv_id="LG12_A", subgroup="RECOMBINANT")
sm$add_group("inv_LG12_B_HOM_INV",     c("CGA007","CGA008"),
             chrom="chr12", inv_id="LG12_B", subgroup="HOM_INV")
sm$add_group("inv_LG17_X_HOM_REF",     c("CGA002","CGA003","CGA004","CGA005"),
             chrom="chr17", inv_id="LG17_X", subgroup="HOM_REF")
sm$add_group("inv_LG17_X_HOM_INV",     c("CGA006","CGA007","CGA008"),
             chrom="chr17", inv_id="LG17_X", subgroup="HOM_INV")
sm$add_group("inv_LG17_X_RECOMBINANT", c("CGA001"),
             chrom="chr17", inv_id="LG17_X", subgroup="RECOMBINANT")
assert(sm$count_groups() >= 8L, "≥8 groups registered")

cat("\n──── get_sample_metadata ────\n")
md <- sm$get_sample_metadata("CGA003")
assert(!is.null(md),                          "CGA003 metadata found")
assert(identical(md$family_id, "F1"),         "family_id pulled from master")
assert(is.null(sm$get_sample_metadata("GHOST")), "unknown → NULL")

cat("\n──── get_sample_groups (reverse lookup) ────\n")
assert(setequal(sm$get_sample_groups("CGA007"),
                c("inv_LG12_A_HOM_INV","inv_LG12_B_HOM_INV","inv_LG17_X_HOM_INV")),
       "CGA007 in 3 HOM_INV groups")
assert(setequal(sm$get_sample_groups("CGA006"),
                c("inv_LG12_A_HET","inv_LG12_A_RECOMBINANT","inv_LG17_X_HOM_INV")),
       "CGA006: cross-status for LG12_A (HET+RECOMB), HOM_INV for LG17_X")
assert(length(sm$get_sample_groups("GHOST")) == 0, "unknown sample → empty")

cat("\n──── Pedigree methods ────\n")
assert(identical(sm$get_family("CGA001"), "F1"),       "CGA001 → F1")
assert(is.na(sm$get_family("CGA010")),                 "NA in master → NA")
assert(is.na(sm$get_family("GHOST")),                  "unknown → NA")
assert(identical(sm$list_families(),
                 c("F1","F2","F3")),                   "list_families sorted")
assert(identical(sm$get_family_members("F1"),
                 c("CGA001","CGA002","CGA003")),       "F1 members")
assert(length(sm$get_family_members("F99")) == 0,      "unknown family empty")

cat("\n──── list_carriers_across_candidates ────\n")
assert(setequal(sm$list_carriers_across_candidates("HOM_INV"),
                c("CGA007","CGA008","CGA009","CGA006")),
       "HOM_INV across all inversions: {CGA006-9}")
assert(setequal(sm$list_carriers_across_candidates("RECOMBINANT"),
                c("CGA010","CGA006","CGA001")),
       "RECOMBINANT across all inversions: {CGA001,6,10}")

cat("\n──── compute_group_overlap ────\n")
ov <- sm$compute_group_overlap("inv_LG12_A_HOM_INV", "inv_LG12_B_HOM_INV")
assert(ov$intersection == 2L,                          "A∩B HOM_INV = {7,8}")
assert(ov$jaccard == 2/3,                              "jaccard = 2/3")
assert(setequal(ov$shared, c("CGA007","CGA008")),      "shared set correct")
ov2 <- sm$compute_group_overlap("inv_LG12_A_HOM_REF", "inv_LG17_X_HOM_REF")
assert(ov2$intersection == 3L,                         "HOM_REF overlap = 3")

cat("\n──── list_recombinant_candidates ────\n")
assert(setequal(sm$list_recombinant_candidates("CGA006"),
                "LG12_A"),                             "CGA006 → LG12_A")
assert(setequal(sm$list_recombinant_candidates("CGA001"),
                "LG17_X"),                             "CGA001 → LG17_X")
assert(length(sm$list_recombinant_candidates("CGA003")) == 0,
       "non-recombinant sample → empty")

cat("\n──── get_groups_for_candidate ────\n")
gfc <- sm$get_groups_for_candidate("LG12_A")
assert(length(gfc$HOM_REF) == 4L,                      "HOM_REF n=4")
assert(length(gfc$HET) == 2L,                          "HET n=2")
assert(length(gfc$HOM_INV) == 3L,                      "HOM_INV n=3")
assert(length(gfc$RECOMBINANT) == 2L,                  "RECOMBINANT n=2")
assert(length(gfc$RECOMBINANT_GC) == 0L,               "RECOMBINANT_GC empty (none registered)")

cat("\n──── find_co_segregating_groups ────\n")
coseg <- sm$find_co_segregating_groups(min_jaccard = 0.5, min_size = 2L)
assert(nrow(coseg) >= 1L, "at least one co-seg pair at jaccard≥0.5")
# inv_LG12_A_HOM_INV (7,8,9) vs inv_LG17_X_HOM_INV (6,7,8) → jaccard 2/4 = 0.5
has_ax_hom <- any(
  (coseg$gid1 == "inv_LG12_A_HOM_INV" & coseg$gid2 == "inv_LG17_X_HOM_INV") |
  (coseg$gid1 == "inv_LG17_X_HOM_INV" & coseg$gid2 == "inv_LG12_A_HOM_INV")
)
assert(isTRUE(has_ax_hom), "LG12_A vs LG17_X HOM_INV co-seg detected")

# Regression: preserve get_carriers (v9 call site)
assert(setequal(sm$get_carriers("LG12_A", "HOM_INV"),
                c("CGA007","CGA008","CGA009")),
       "regression: get_carriers still works")

cat("\n──── Query API composites ────\n")
qr <- reg$query

tr <- qr$nested_family_tree("LG12_A")
assert(!is.null(tr),                                   "tree(LG12_A) built")
assert(tr$n_children == 1L,                            "A has 1 child (B)")
assert(tr$children[[1]]$cid == "LG12_B",               "child[1] = B")
assert(is.null(qr$nested_family_tree("GHOST")),        "unknown root → NULL")

oc12 <- qr$overlap_clusters("chr12")
# chr12 candidates: LG12_A, LG12_B, LG12_E. A⊃B shares cluster; E disjoint.
assert(nrow(oc12) == 3L,                               "3 candidates on chr12")
a_cl <- oc12$cluster_id[oc12$candidate_id == "LG12_A"]
b_cl <- oc12$cluster_id[oc12$candidate_id == "LG12_B"]
e_cl <- oc12$cluster_id[oc12$candidate_id == "LG12_E"]
assert(a_cl == b_cl,                                   "A,B same cluster")
assert(a_cl != e_cl,                                   "E in separate cluster")

oc99 <- qr$overlap_clusters("chr99")
assert(nrow(oc99) == 0L,                               "unknown chrom empty")

ld <- qr$sample_inversion_load("CGA007")
assert(ld$total == 3L,                                 "CGA007 in 3 inversions")
assert(ld$HOM_INV == 3L,                               "all 3 as HOM_INV")
assert(setequal(ld$candidates, c("LG12_A","LG12_B","LG17_X")),
       "CGA007 candidates enumerated")

ld1 <- qr$sample_inversion_load("CGA001")
assert(ld1$HOM_REF == 1L,                              "CGA001: 1 HOM_REF (LG12_A)")
assert(ld1$RECOMBINANT == 1L,                          "CGA001: 1 RECOMBINANT (LG17_X)")

ldg <- qr$sample_inversion_load("GHOST")
assert(ldg$total == 0L,                                "unknown sample → 0 load")

cat("\nAll sample-registry extension tests (and composite query tests) passed.\n")
cat("Temp root (inspect if needed): ", root, "\n", sep = "")
