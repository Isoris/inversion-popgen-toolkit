#!/usr/bin/env bash
# =============================================================================
# run_phase4_v5_new_evidence.sh
# =============================================================================
# Runs the v5 evidence scripts + v6 cross-species bridge + v7 bp_pipeline
# bridge + v7 final synthesis to produce a final_catalog.tsv.
#
# Order:
#   1. cheat29b_assembled_junction.py  (assembled-junction evidence)
#   2. bnd_sided_support.py            (single-sided BND at known candidates)
#   3. cross_species_bridge_v6.py         (synteny_dollo + flank_coherence)
#   4. assign_structural_class_v7.py   (final label synthesis)
#
# Inputs (env vars with defaults):
#   CANDIDATES        — TSV with candidate_id, chrom, left_bp, right_bp
#   DELLY_INV_VCF     — MODULE_4D strict catalog (required for #1)
#   DELLY_BND_VCF     — MODULE_4E strict catalog (required for #1 and #2)
#   MANTA_INV_VCF     — MODULE_4G INV catalog (optional for #1)
#   MANTA_RAW_VCF     — MODULE_4G raw pre-conversion VCF (optional for #2)
#   BETWEEN_SPECIES_BED  — catfish-synteny-toolkit STEP_02 output (optional)
#   FLANK_COHERENCE_TSV  — catfish-synteny-toolkit STEP_09c output (optional)
#   POLARIZED_TSV     — catfish-synteny-toolkit STEP_11 output (optional)
#   REGISTRIES_ROOT   — path to inversion-popgen-toolkit/registries/ (optional)
#   KEYS_DIR          — path to per-candidate keys.tsv dir (required for #4)
#   OUTDIR            — where to write per-candidate blocks + final catalog
#
# Script location conventions:
#   SCRIPT_DIR defaults to $(dirname this script) — must contain the four .py files.
# =============================================================================
set -euo pipefail

: "${CANDIDATES:?need CANDIDATES=path/to/candidates.tsv}"
: "${OUTDIR:?need OUTDIR=path/to/output}"
: "${KEYS_DIR:=$OUTDIR/keys_dummy}"   # must be populated by upstream 4e
DELLY_INV_VCF="${DELLY_INV_VCF:-}"
DELLY_BND_VCF="${DELLY_BND_VCF:-}"
MANTA_INV_VCF="${MANTA_INV_VCF:-}"
MANTA_RAW_VCF="${MANTA_RAW_VCF:-}"
BETWEEN_SPECIES_BED="${BETWEEN_SPECIES_BED:-}"
FLANK_COHERENCE_TSV="${FLANK_COHERENCE_TSV:-}"
POLARIZED_TSV="${POLARIZED_TSV:-}"
DOLLO_TSV="${DOLLO_TSV:-}"
BP_CONSENSUS_TSV="${BP_CONSENSUS_TSV:-}"
BP_FRAGMENTS_TSV="${BP_FRAGMENTS_TSV:-}"
BP_PER_METHOD_TSV="${BP_PER_METHOD_TSV:-}"
REGISTRIES_ROOT="${REGISTRIES_ROOT:-}"

SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
V5_BLOCKS="$OUTDIR/v5_blocks"
FINAL="$OUTDIR/final"
mkdir -p "$V5_BLOCKS" "$FINAL"

echo "=============================================================================="
echo "v5 pipeline: $SCRIPT_DIR"
echo "  candidates  = $CANDIDATES"
echo "  outdir      = $OUTDIR"
echo "=============================================================================="

# -------------------------------------------------------------
# 1. Assembled-junction evidence
# -------------------------------------------------------------
echo "[v5] STEP 1: assembled-junction forensics"
python3 "$SCRIPT_DIR/cheat29b_assembled_junction.py" \
    --candidates "$CANDIDATES" \
    --delly_inv_vcf "$DELLY_INV_VCF" \
    --delly_bnd_vcf "$DELLY_BND_VCF" \
    --manta_inv_vcf "$MANTA_INV_VCF" \
    --outdir "$V5_BLOCKS" \
    ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}

# -------------------------------------------------------------
# 2. Single-sided BND support
# -------------------------------------------------------------
echo "[v5] STEP 2: single-sided BND support at known candidates"
python3 "$SCRIPT_DIR/bnd_sided_support.py" \
    --candidates "$CANDIDATES" \
    --delly_bnd_vcf "$DELLY_BND_VCF" \
    --manta_raw_vcf "$MANTA_RAW_VCF" \
    --outdir "$V5_BLOCKS" \
    ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}

# -------------------------------------------------------------
# 3. Cross-species bridge
# -------------------------------------------------------------
echo "[v5] STEP 3: cross-species synteny bridge (synteny_dollo + flanks)"
python3 "$SCRIPT_DIR/cross_species_bridge_v6.py" \
    --candidates "$CANDIDATES" \
    --between_species_bed "$BETWEEN_SPECIES_BED" \
    --flank_coherence_tsv "$FLANK_COHERENCE_TSV" \
    --polarized_tsv "$POLARIZED_TSV" \
    ${DOLLO_TSV:+--dollo_tsv "$DOLLO_TSV"} \
    --outdir "$V5_BLOCKS" \
    ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}


# -------------------------------------------------------------
# 3b. bp_pipeline bridge (v7 addition — reads breakpoint_pipeline TSVs)
# -------------------------------------------------------------
if [[ -n "$BP_CONSENSUS_TSV" && -n "$BP_FRAGMENTS_TSV" ]]; then
    echo "[v5] STEP 3b: bp_pipeline bridge (fragment_distribution + consensus)"
    python3 "$SCRIPT_DIR/bp_pipeline_bridge.py" \
        --consensus_tsv "$BP_CONSENSUS_TSV" \
        --fragments_tsv "$BP_FRAGMENTS_TSV" \
        ${BP_PER_METHOD_TSV:+--per_method_tsv "$BP_PER_METHOD_TSV"} \
        --outdir "$V5_BLOCKS" \
        ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}
else
    echo "[v5] STEP 3b: bp_pipeline bridge — SKIPPED (set BP_CONSENSUS_TSV + BP_FRAGMENTS_TSV to enable)"
fi

# -------------------------------------------------------------
# 4. Final structural-class synthesis
# -------------------------------------------------------------
echo "[v5] STEP 4: final structural-class assignment"
python3 "$SCRIPT_DIR/assign_structural_class_v7.py" \
    --candidates "$CANDIDATES" \
    --keys_dir "$KEYS_DIR" \
    --v5_blocks_dir "$V5_BLOCKS" \
    --outdir "$FINAL"

echo "=============================================================================="
echo "[v5] done. Final catalog: $FINAL/final_catalog.tsv"
echo "=============================================================================="
