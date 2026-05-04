#!/usr/bin/env bash
# =============================================================================
# 00_phase3_config.sh — v2: dual-caller (DELLY2 + Manta)
# =============================================================================
#
# v2 changes:
#   - Adds Manta INV candidates as a second source alongside DELLY2
#   - Manta raw pre-conversion VCF for BND inversion signal (STEP_B06)
#   - Caller-aware evidence field documentation
#   - Unified candidate table with caller provenance column
#
# CALLER FIELD DIFFERENCES (for reference):
#   ┌──────────────┬────────────────────┬─────────────────────────┐
#   │              │ DELLY2             │ Manta                   │
#   ├──────────────┼────────────────────┼─────────────────────────┤
#   │ PE support   │ INFO/PE (integer)  │ FORMAT/PR (ref,alt)     │
#   │ SR support   │ INFO/SR (integer)  │ FORMAT/SR (ref,alt)     │
#   │ Per-sample   │ DR/DV (span ref/   │ PR:ref,alt SR:ref,alt   │
#   │   evidence   │  alt), RR/RV       │                         │
#   │              │ (junction ref/alt) │                         │
#   │ Connection   │ INFO/CT (3to3,     │ INFO flags INV3/INV5    │
#   │   type       │  5to5, 3to5, 5to3) │ (on converted INV recs) │
#   │ PRECISE      │ SR>0 & SRQ≥0.8    │ assembly succeeded      │
#   │ Quality      │ QUAL (phred)       │ QUAL (phred)            │
#   │ CIPOS/CIEND  │ same               │ same                    │
#   │ IMPRECISE    │ same flag name     │ same flag name          │
#   │ Breakpoint   │ INFO/CONSENSUS     │ (not available)         │
#   │   consensus  │                    │                         │
#   └──────────────┴────────────────────┴─────────────────────────┘
#
# The BAM-level evidence extraction (STEP_A02) is caller-agnostic — it reads
# BAMs directly with pysam looking for discordant pairs, split reads, and
# soft clips. The caller-specific fields are only used in STEP_A01 (candidate
# extraction) and STEP_B05/06 (cross-caller concordance).
# =============================================================================

# ── Source parent inversion config ───────────────────────────────────────────
SCRIPT_DIR_BPV="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INV_CONFIG="${SCRIPT_DIR_BPV}/../00_inversion_config.sh"
[[ -f "${INV_CONFIG}" ]] || INV_CONFIG="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}/inversion_codebase_v8.5/00_inversion_config.sh"
[[ -f "${INV_CONFIG}" ]] || { echo "ERROR: Cannot find 00_inversion_config.sh" >&2; exit 1; }
source "${INV_CONFIG}"

# ── DELLY2 SV paths ─────────────────────────────────────────────────────────
DELLY_BASE="${BASE}/MODULE_4B_DEL_Delly"
DELLY_INV_BASE="${BASE}/MODULE_4D_INV_Delly"
DELLY_BND_BASE="${BASE}/MODULE_4F_BND_Delly"
DELLY_MARKDUP_DIR="${DELLY_BASE}/00_markdup"

# DELLY2 INV catalog (from MODULE_4D)
DELLY_INV_VCF="${DELLY_INV_BASE}/07_final_catalogs/catalog_226.INV.vcf.gz"
DELLY_INV_BED="${DELLY_INV_BASE}/07_final_catalogs/catalog_226.INV.bed"
DELLY_INV_GT="${DELLY_INV_BASE}/07_final_catalogs/catalog_226.INV.GT_matrix.tsv"

# DELLY2 BND catalog (for BND inversion signal in STEP_B06)
DELLY_BND_VCF="${DELLY_BND_BASE}/07_final_catalogs/catalog_226.BND.vcf.gz"

# ── Manta SV paths (from MODULE_4H) ─────────────────────────────────────────
MANTA_BASE="${BASE}/MODULE_4H_ALL_Manta"

# Manta INV catalog (converted from BND via convertInversion_py3.py)
# These records carry INV3/INV5 flags and EVENT tags.
MANTA_INV_VCF="${MANTA_BASE}/05_final_catalogs/catalog_226.INV.PASS.vcf.gz"

# Manta post-conversion BND catalog — inter-chromosomal translocations ONLY.
# NO INV3/INV5 tags remain (all inversion BNDs were converted to <INV>).
MANTA_BND_POST_VCF="${MANTA_BASE}/05_final_catalogs/catalog_226.BND.PASS.vcf.gz"

# Manta RAW pre-conversion merged VCF — still has BND records with inversion
# orientation (for STEP_B06 BND inversion signal analysis).
MANTA_RAW_MERGED_VCF="${MANTA_BASE}/02_merged_cohort/cohort_226.ALL.raw.vcf.gz"

# ── Inversion pipeline outputs to consume ────────────────────────────────────
SNAKE2_BANDS="${SNAKE2_DIR}/snake2_band_assignments.tsv.gz"
SNAKE_SCORES="${SNAKE1_DIR}/snake_candidate_regions.tsv.gz"
DECOMP_DIR="${DECOMPOSITION_DIR}"

# ── MODULE_5A2 output directories ────────────────────────────────────────────
BPV_ROOT="${INVDIR}/22_breakpoint_validation"
BPV_EVIDENCE="${BPV_ROOT}/01_per_sample_evidence"
BPV_CONCORDANCE="${BPV_ROOT}/02_concordance"
BPV_STATS="${BPV_ROOT}/03_statistical_tests"
BPV_SEEDS="${BPV_ROOT}/04_deconvolution_seeds"
BPV_PLOTS="${BPV_ROOT}/05_plots"
BPV_LOGS="${BPV_ROOT}/logs"

# ── Validation parameters ────────────────────────────────────────────────────
BP_WINDOW=300          # ±bp around breakpoint for evidence search
MIN_MAPQ=20            # minimum MAPQ for supporting reads
MIN_CLIP_LEN=10        # minimum soft-clip length
MIN_DISCORDANT=1       # min discordant pairs for support=yes

# DELLY2-specific filters
MIN_DELLY_QUAL=300     # only validate INVs above this QUAL
MIN_DELLY_PE=3         # only validate INVs with PE >= this

# Manta-specific filters
MIN_MANTA_QUAL=50      # only validate INVs above this QUAL
MIN_MANTA_TOTAL_ALT=6  # sum(PR_alt + SR_alt) across carriers >= this

# Shared size filters
MAX_INV_SIZE_MB=20     # skip INVs larger than this
MIN_INV_SIZE_BP=5000   # skip INVs smaller than this

# ── Statistical test parameters ──────────────────────────────────────────────
FISHER_ALPHA=0.05
BONFERRONI=1

# ── Seed generation parameters ───────────────────────────────────────────────
SEED_MIN_CONCORDANCE=0.80
SEED_MIN_SUPPORT_FRAC=0.60
SEED_MIN_SAMPLES=5

# ── Sample lists ─────────────────────────────────────────────────────────────
SAMPLES_ALL="${BASE}/MODULE_4B_DEL_Delly/samples_all_226.txt"
SAMPLES_UNRELATED="${BASE}/MODULE_4B_DEL_Delly/samples_unrelated_81.txt"

# ── Evidence registry wiring (FIX 29 v2, 2026-04-17) ─────────────────────────
# STEP_D03's OR test writes existence_layer_d blocks per candidate, and
# STEP_B06's BND rescue writes bnd_rescue_* keys into existence_layer_b.
# These feed the 4-layer evidence model (A=PCA, B=SV callers, C=GHSL,
# D=OR association) consumed by C01f's compute_group_validation() and
# C01d's Layer B/D accounting in 4a.
#
# Set REGISTRIES_ROOT="" to run phase_3 standalone without registry writes.
REGISTRIES_ROOT="${REGISTRIES_ROOT:-${BASE}/inversion_codebase_v8.5/registries}"

# Optional inv_id → candidate_id translation TSV. When unset/missing,
# phase_3's inv_id (e.g. "C_gar_LG12:52400000-54600000:delly") is used
# verbatim as the registry candidate_id. When running alongside C01i
# which already registered candidates as LG12_17-style IDs, supply a
# map so Layer D writes land in the right per-candidate folder.
CANDIDATE_MAP="${CANDIDATE_MAP:-}"

# ── Helper ───────────────────────────────────────────────────────────────────
bpv_log() { echo "[$(date '+%F %T')] [BPV] $*"; }
bpv_err() { echo "[$(date '+%F %T')] [BPV] [ERROR] $*" >&2; }
bpv_die() { bpv_err "$@"; exit 1; }

bpv_init_dirs() {
  mkdir -p \
    "${BPV_ROOT}" "${BPV_EVIDENCE}" "${BPV_CONCORDANCE}" \
    "${BPV_STATS}" "${BPV_SEEDS}" "${BPV_PLOTS}" "${BPV_LOGS}"
}
