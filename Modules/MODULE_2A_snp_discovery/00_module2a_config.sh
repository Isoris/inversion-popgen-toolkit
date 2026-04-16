#!/usr/bin/env bash
###############################################################################
# 00_module2a_config.sh — MODULE_2A / STEP_1 centralized configuration
#
# ALL paths, ANGSD parameters, tool settings, thinning distances, K ranges,
# seed lists, and output conventions live here. No other script in STEP_1
# should hardcode these values.
#
# Usage:
#   source 00_module2a_config.sh
#   # Then all variables below are available in the calling script.
###############################################################################

# =============================================================================
# PROJECT PATHS
# =============================================================================
export BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
export REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
export FAI="${REF}.fai"
export BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"
export SAMPLE_LIST="${BASE}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"

# Callable mask prefix (from S00/S01)
export MASK_PREFIX="${BASE}/fClaHyb_Gar_LG.mask_regions"
export CALLABLE_SITES="${MASK_PREFIX}.normalACGT.renamed.1-based.angsd"

# Output roots
export GLOBAL_DIR="${BASE}/popstruct_global"
export THIN_DIR="${BASE}/popstruct_thin"

# =============================================================================
# GENOME
# =============================================================================
export GENOME_SIZE=963905721
export N_SAMPLES=226

# =============================================================================
# ANGSD SHARED FILTERS (used in SAF, biSNP, BEAGLE steps)
# =============================================================================
export ANGSD_GL=1               # SAMtools GL model
export ANGSD_MINQ=25            # base quality
export ANGSD_MINMAPQ=25         # mapping quality
export ANGSD_BAQ=1              # BAQ recalibration
export ANGSD_C=50               # excessive mismatch downweight
export ANGSD_MINDEPTHIND=3      # min depth per individual
export ANGSD_MAXDEPTHIND=57     # max depth per individual
export ANGSD_MININD=200         # min individuals with data
export ANGSD_REMOVE_BADS=1
export ANGSD_UNIQUE_ONLY=1
export ANGSD_ONLY_PROPER_PAIRS=1

# =============================================================================
# SNP DISCOVERY (biSNP step)
# =============================================================================
export SNP_PVAL="1e-6"
export MIN_MAF=0.05

# =============================================================================
# SFS
# =============================================================================
export SFS_BOOTSTRAPS=100

# =============================================================================
# THINNING
# =============================================================================
# Fine thinning (4-col with maj/min alleles for BEAGLE -doMajorMinor 3)
export THIN_FINE=(200 500 1000)
# Broad thinning (2-col for whole-genome structure)
export THIN_BROAD=(5000 10000 25000)
# All thinning distances
export THIN_ALL=(200 500 1000 5000 10000 25000)

# =============================================================================
# ANCESTRY / STRUCTURE
# =============================================================================
export K_MIN=2
export K_MAX=12
export SEEDS=(1 2 3)
export PCANGSD_EIG=6
export PCANGSD_MAF=0.05
export PCANGSD_ITER=250
export NGSADMIX_MINMAF=0.05

# =============================================================================
# EVALADMIX
# =============================================================================
export EVALADMIX_BIN="${BASE}/evalAdmix"
export EVALADMIX_NITS=5
export EVALADMIX_MISTOL=0.05
export EVALADMIX_MINMAF=0.05

# =============================================================================
# RELATEDNESS
# =============================================================================
# Default thinning for ngsRelate (thin-500 is standard)
export RELATE_THIN=500
# KING-based theta cutoffs (single threshold for classification)
export THETA_DUP_MZ=0.354
export THETA_FIRST_DEGREE=0.177
export THETA_SECOND_DEGREE=0.0884
export THETA_THIRD_DEGREE=0.0442

# NAToRA low/high cutoff pairs per relatedness class
export NATORA_DUP_MZ_LOW=0.3540
export NATORA_DUP_MZ_HIGH=0.5000
export NATORA_FIRST_LOW=0.1770
export NATORA_FIRST_HIGH=0.3540
export NATORA_SECOND_LOW=0.0884
export NATORA_SECOND_HIGH=0.1770
export NATORA_THIRD_LOW=0.0442
export NATORA_THIRD_HIGH=0.0884

# =============================================================================
# BEST-SEED SELECTION
# =============================================================================
export PALETTE_NAME="catfish_ngsadmix_v1"
# Selection rule: highest loglik -> lowest mean |evalAdmix residual| -> lowest seed
export BEST_SEED_RULE="highest_loglik_then_lowest_mean_abs_resid"

# =============================================================================
# SLURM DEFAULTS (can be overridden per-step)
# =============================================================================
export DEFAULT_PARTITION="compute"
export DEFAULT_THREADS="${SLURM_CPUS_PER_TASK:-80}"

# =============================================================================
# CONVENIENCE: ANGSD common args as a bash array
# =============================================================================
angsd_common_args() {
  echo "-GL ${ANGSD_GL}" \
       "-minQ ${ANGSD_MINQ} -minMapQ ${ANGSD_MINMAPQ}" \
       "-baq ${ANGSD_BAQ} -C ${ANGSD_C}" \
       "-setMinDepthInd ${ANGSD_MINDEPTHIND} -setMaxDepthInd ${ANGSD_MAXDEPTHIND}" \
       "-minInd ${ANGSD_MININD}" \
       "-remove_bads ${ANGSD_REMOVE_BADS} -uniqueOnly ${ANGSD_UNIQUE_ONLY}" \
       "-only_proper_pairs ${ANGSD_ONLY_PROPER_PAIRS}" \
       "-doCounts 1 -doDepth 1"
}
export -f angsd_common_args
