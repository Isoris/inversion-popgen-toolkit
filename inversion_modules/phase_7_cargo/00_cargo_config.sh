#!/usr/bin/env bash
# =============================================================================
# 00_cargo_config.sh — MODULE_6_Cargo (phase 7) configuration
#
# Sources the main inversion config and layers cargo-specific paths on top.
# All cargo outputs live under ${INVDIR}/30_cargo/ — outside the code tree,
# matching the convention established by chat 15.
#
# Usage in any cargo script:
#   source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/00_cargo_config.sh"
# =============================================================================
set -u

# ── Source main inversion config ─────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Walk up to find 00_inversion_config.sh — try direct parent first, then
# inversion_modules/, then fall back to env-provided path.
for _try in \
    "${SCRIPT_DIR}/../00_inversion_config.sh" \
    "${SCRIPT_DIR}/../../inversion_modules/00_inversion_config.sh" \
    "${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}/inversion_modules/00_inversion_config.sh"; do
  if [[ -f "$_try" ]]; then
    source "$_try"
    break
  fi
done
unset _try

# Sanity check — the main config must have set BASE / INVDIR / REF_GFF3
: "${BASE:?00_inversion_config.sh did not set BASE}"
: "${INVDIR:?00_inversion_config.sh did not set INVDIR}"
: "${REF_GFF3:?00_inversion_config.sh did not set REF_GFF3}"

# ── Cargo output root ────────────────────────────────────────────────────────
export CARGO_DIR="${CARGO_DIR:-${INVDIR}/30_cargo}"
mkdir -p "${CARGO_DIR}"

# ── Sub-directories ──────────────────────────────────────────────────────────
export CARGO_INVENTORY_DIR="${CARGO_DIR}/inventory"          # 6A outputs
export CARGO_BURDEN_DIR="${CARGO_DIR}/per_arrangement_burden" # 6C Level 1
export CARGO_CONFIG_DIR="${CARGO_DIR}/config_spectrum"        # 6C Level 2
export CARGO_ENRICH_DIR="${CARGO_DIR}/functional_enrichment"  # 6B
export CARGO_SYNTH_DIR="${CARGO_DIR}/cross_inversion"         # 6E synthesis
export CARGO_FIG_DIR="${CARGO_DIR}/figures"                   # 6E figures
export CARGO_LOGS_DIR="${CARGO_DIR}/logs"

mkdir -p "${CARGO_INVENTORY_DIR}" "${CARGO_BURDEN_DIR}" \
         "${CARGO_CONFIG_DIR}" "${CARGO_ENRICH_DIR}" \
         "${CARGO_SYNTH_DIR}" "${CARGO_FIG_DIR}" "${CARGO_LOGS_DIR}"

# ── Gene annotation derived files (built by STEP_60) ─────────────────────────
export GENE_BED="${CARGO_DIR}/gene_intervals.bed.gz"
export GENE_BED_TBI="${GENE_BED}.tbi"
export GENE_LENGTH_TSV="${CARGO_DIR}/gene_lengths.tsv"     # gene_id  length_bp  cds_length_bp
export TRANSCRIPT_GENE_MAP="${CARGO_DIR}/transcript_gene_map.tsv"  # transcript_id  gene_id

# ── eggNOG annotation (built by STEP_61, optional) ───────────────────────────
export EGGNOG_OUT="${CARGO_DIR}/eggnog_annotations.tsv"
# Resolved per-gene table: gene_id  best_OG  preferred_name  go_terms  kegg_pathway  family
export GENE_FUNCTION_TSV="${CARGO_DIR}/gene_function.tsv"

# ── Repeat annotation (optional; if unavailable, repeat columns are NA) ──────
# Set REPEAT_BED in the environment if you have RepeatMasker / EarlGrey output.
export REPEAT_BED="${REPEAT_BED:-}"

# ── MODULE_CONSERVATION outputs we consume ───────────────────────────────────
# MODCONS env var is expected to point to the conservation pipeline base.
# Default mirrors the convention used by STEP_14/STEP_15.
export MODCONS="${MODCONS:-${BASE}/conservation_core}"
export VARIANT_MASTER="${MODCONS}/16_merged_variant_tables/variant_master_scored.tsv"
export VESM_SCORES="${MODCONS}/05B_vesm/vesm_variant_scores.tsv"
export NORMALIZED_VCF_DIR="${MODCONS}/03_variants/normalized"

# ── Diagnostic gate thresholds ───────────────────────────────────────────────
# Level 1 (descriptive per-arrangement burden, NS/SS, FST, SFS)
export CARGO_LEVEL1_MIN_HOMO="${CARGO_LEVEL1_MIN_HOMO:-5}"
# Level 2 (configuration spectrum)
export CARGO_LEVEL2_MIN_HOMO="${CARGO_LEVEL2_MIN_HOMO:-10}"
export CARGO_LEVEL2_MIN_NS_PER_GENE="${CARGO_LEVEL2_MIN_NS_PER_GENE:-5}"

# ── Configuration spectrum settings ──────────────────────────────────────────
export CARGO_PERM_N="${CARGO_PERM_N:-1000}"      # permutations for null
export CARGO_MIN_DP="${CARGO_MIN_DP:-3}"          # genotype min depth (matches STEP_15)
export CARGO_MAF_FLOOR="${CARGO_MAF_FLOOR:-0.0}"  # set >0 to drop singletons in MI
export CARGO_THREADS="${CARGO_THREADS:-8}"

# ── Functional enrichment settings ───────────────────────────────────────────
export CARGO_ENRICH_FDR="${CARGO_ENRICH_FDR:-0.05}"
export CARGO_ENRICH_MIN_GENES="${CARGO_ENRICH_MIN_GENES:-3}"
export CARGO_MATCHED_BG_N="${CARGO_MATCHED_BG_N:-100}"  # number of matched collinear regions to draw

# ── Registry hooks ───────────────────────────────────────────────────────────
# These are all already exported by the upstream config; re-asserted here for clarity.
: "${RESULTS_REGISTRY_DIR:?results_registry not configured}"
: "${SAMPLE_REGISTRY:?sample_registry not configured}"
: "${REGISTRIES_DATA_DIR:?registries data root not configured}"
export REGISTRY_LOADER_PY="${REGISTRY_LOADER_PY:-${BASE}/registries/api/python/registry_loader.py}"
export REGISTRY_LOADER_R="${REGISTRY_LOADER_R:-${BASE}/registries/api/R/registry_loader.R}"

# ── Active sample group ──────────────────────────────────────────────────────
export SAMPLE_GROUP="${SAMPLE_GROUP:-all_226}"

# ── Pretty banner ────────────────────────────────────────────────────────────
echo "[cargo_config] CARGO_DIR=${CARGO_DIR}"
echo "[cargo_config] MODCONS=${MODCONS}"
echo "[cargo_config] level1_min_homo=${CARGO_LEVEL1_MIN_HOMO} level2_min_homo=${CARGO_LEVEL2_MIN_HOMO} perm_n=${CARGO_PERM_N}"
