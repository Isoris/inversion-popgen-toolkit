#!/usr/bin/env Rscript
# =============================================================================
# test_interval_registry_extensions.R — chat 11 test fixture
# =============================================================================
# Exercises the 9 new interval-registry methods added in chat 11:
#   get_children, get_parent, get_ancestors, get_descendants,
#   get_overlapping, get_nested_within, classify_relationship,
#   update_candidate, bulk_add_candidates
#
# Synthetic topology (chr12 + chr17):
#   LG12_A (outer)          [ 1Mb, 10Mb]
#     └─ LG12_B (nested)    [ 3Mb,  6Mb]  parent_id=LG12_A
#         └─ LG12_C (deep)  [ 4Mb,  5Mb]  parent_id=LG12_B
#   LG12_D (partial overlap [ 8Mb, 12Mb]
#           with A, not nested in A)
#   LG12_E (disjoint)       [20Mb, 25Mb]
#   LG17_X (other chrom)    [ 1Mb,  5Mb]
#
# Run from repo root:
#   Rscript registries/tests/test_interval_registry_extensions.R
#
# Exit 0 on all-pass, exit 1 on first failure (with stop()).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# Resolve registry loader relative to this script
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0) script_path <- "registries/tests/test_interval_registry_extensions.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))

loader <- c(
  file.path(script_dir, "..", "api", "R", "registry_loader.R"),
  "registries/api/R/registry_loader.R"
)
loader <- loader[file.exists(loader)][1]
if (is.na(loader)) stop("Could not find registry_loader.R")
source(loader)

# Use a temp root — no dependencies on BASE or real data
root <- tempfile("int_ext_test_")
dir.create(root, recursive = TRUE)
reg <- load_registry(root, create_if_missing = TRUE)
iv  <- reg$intervals

assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop("FAIL: ", msg)
  cat("  ok: ", msg, "\n", sep = "")
}

cat("──── Setup: synthetic 6-candidate topology ────\n")
stopifnot(iv$add_candidate("LG12_A", "chr12",  1e6, 10e6))
stopifnot(iv$add_candidate("LG12_B", "chr12",  3e6,  6e6, parent_id = "LG12_A"))
stopifnot(iv$add_candidate("LG12_C", "chr12",  4e6,  5e6, parent_id = "LG12_B"))
stopifnot(iv$add_candidate("LG12_D", "chr12",  8e6, 12e6))
stopifnot(iv$add_candidate("LG12_E", "chr12", 20e6, 25e6))
stopifnot(iv$add_candidate("LG17_X", "chr17",  1e6,  5e6))
assert(iv$count_candidates() == 6L, "6 candidates registered")

cat("\n──── get_children ────\n")
assert(identical(iv$get_children("LG12_A"), "LG12_B"), "A has 1 child: B")
assert(identical(iv$get_children("LG12_B"), "LG12_C"), "B has 1 child: C")
assert(length(iv$get_children("LG12_C")) == 0,         "C is a leaf")
assert(length(iv$get_children("GHOST"))  == 0,         "unknown CID returns empty")

cat("\n──── get_parent ────\n")
assert(identical(iv$get_parent("LG12_C"), "LG12_B"),   "C parent = B")
assert(is.na(iv$get_parent("LG12_A")),                 "A is root (NA)")
assert(is.na(iv$get_parent("GHOST")),                  "unknown returns NA")

cat("\n──── get_ancestors ────\n")
assert(identical(iv$get_ancestors("LG12_C"),
                 c("LG12_B", "LG12_A")),               "C ancestors: B→A (closest first)")
assert(length(iv$get_ancestors("LG12_A")) == 0,        "A has no ancestors")

cat("\n──── get_descendants ────\n")
assert(identical(sort(iv$get_descendants("LG12_A")),
                 c("LG12_B", "LG12_C")),               "A has 2 descendants (full)")
assert(identical(iv$get_descendants("LG12_A", depth = 1L),
                 "LG12_B"),                            "A depth=1 → just B")
assert(length(iv$get_descendants("LG12_C")) == 0,      "C leaf → no descendants")

cat("\n──── get_overlapping ────\n")
ov <- iv$get_overlapping("chr12", 5e6, 9e6)
assert(setequal(ov, c("LG12_A","LG12_B","LG12_C","LG12_D")),
       "overlap [5M,9M] hits A,B,C,D")
ov2 <- iv$get_overlapping("chr12", 5e6, 9e6, exclude_cid = "LG12_A")
assert(setequal(ov2, c("LG12_B","LG12_C","LG12_D")), "exclude A drops A")
assert(length(iv$get_overlapping("chr99", 1, 2)) == 0, "unknown chrom empty")

cat("\n──── get_nested_within (coord-based, ignores parent_id) ────\n")
assert(setequal(iv$get_nested_within("LG12_A"),
                c("LG12_B","LG12_C")),                 "A strictly contains B,C")
assert(identical(iv$get_nested_within("LG12_B"), "LG12_C"), "B contains C")
assert(length(iv$get_nested_within("LG12_D")) == 0,    "D (partial overlap) contains nothing")

cat("\n──── classify_relationship ────\n")
cases <- list(
  list("LG12_A", "LG12_B", "nested_2_in_1"),
  list("LG12_B", "LG12_A", "nested_1_in_2"),
  list("LG12_A", "LG12_D", "partial_overlap"),
  list("LG12_A", "LG12_E", "disjoint"),
  list("LG12_A", "LG17_X", "disjoint"),
  list("LG12_A", "LG12_A", "equal")
)
for (cc in cases) {
  got <- iv$classify_relationship(cc[[1]], cc[[2]])
  assert(identical(got, cc[[3]]),
         sprintf("%s vs %s → %s", cc[[1]], cc[[2]], cc[[3]]))
}
assert(is.na(iv$classify_relationship("LG12_A", "GHOST")), "unknown CID → NA")

cat("\n──── update_candidate ────\n")
assert(iv$update_candidate("LG12_D", parent_id = "LG12_A"),
       "update: set D.parent_id = A")
assert(identical(iv$get_parent("LG12_D"), "LG12_A"),   "D parent now = A")
assert(setequal(iv$get_children("LG12_A"),
                c("LG12_B","LG12_D")),                 "A children now include D")
# Unknown column → warn, not fatal
w <- tryCatch(
  withCallingHandlers(iv$update_candidate("LG12_A", bogus = 1),
                      warning = function(w) { invokeRestart("muffleWarning") }),
  error = function(e) FALSE
)
assert(isTRUE(w), "unknown column warns but succeeds")
# Unknown CID → FALSE
assert(isFALSE(suppressWarnings(iv$update_candidate("GHOST", parent_id = "X"))),
       "unknown CID returns FALSE")

cat("\n──── bulk_add_candidates ────\n")
new_dt <- data.table(
  candidate_id = c("LG03_A", "LG03_B", "LG12_A"),   # LG12_A already exists
  chrom        = c("chr3",   "chr3",   "chr12"),
  start_bp     = c(1e6,      2e6,      999),
  end_bp       = c(5e6,      4e6,      1001),
  parent_id    = c(NA_character_, "LG03_A", NA_character_)
)
n_ins <- iv$bulk_add_candidates(new_dt)
assert(n_ins == 2L,                      "bulk: 2 inserted (1 dup skipped)")
assert(iv$count_candidates() == 8L,      "total now 8")
assert(identical(iv$get_parent("LG03_B"), "LG03_A"),
       "bulk-inserted parent_id stored")

bad_dt <- data.table(candidate_id = "X", chrom = "chr1")  # no start/end
ok <- tryCatch({ iv$bulk_add_candidates(bad_dt); FALSE },
               error = function(e) TRUE)
assert(isTRUE(ok), "bulk: missing required cols → error")

cat("\nAll interval-registry extension tests passed.\n")
cat("Temp root (inspect if needed): ", root, "\n", sep = "")
