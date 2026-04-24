#!/usr/bin/env bash
# =============================================================================
# run_evidence_biology.sh — cohort-level phase_8 → phase_9 synthesis driver
# =============================================================================
# History:
#   pass 15 (2026-04-24): renamed from run_phase4_v5_new_evidence.sh.
#   pass 18 (2026-04-24): step naming regularized to STEP_NN_* (letter sub-
#                         steps for parallel annotation). Relationships
#                         between the five scripts made explicit.
#
# =============================================================================
# THE PIPELINE — a single directed chain with one parallel sub-step
# =============================================================================
#
#   ┌───────────────────────────────────────────────────────────────────┐
#   │  STEP 1: assembled-junction forensics (Q4 mechanism)              │
#   │    q4_mechanism/STEP_01_assembled_junction.py                      │
#   │    in:  DELLY INV+BND VCFs, Manta INV VCF, candidates.tsv          │
#   │    out: <outdir>/<cid>/mechanism_assembled.json                    │
#   │         flat keys: q4b_asm_precise_record_available,               │
#   │                    q4b_asm_junction_class, q4b_asm_source          │
#   └───────────────────────────────────────────────────────────────────┘
#                               │
#                               ▼
#   ┌───────────────────────────────────────────────────────────────────┐
#   │  STEP 2: single-sided BND rescue scoring (Q7 existence)           │
#   │    q7_existence_audit/STEP_02_bnd_sided_support.py                 │
#   │    in:  DELLY BND VCF, Manta raw VCF, candidates.tsv               │
#   │    out: <outdir>/<cid>/bnd_sided_support.json                      │
#   │         flat keys: q7b_bnd_left_support, q7b_bnd_right_support     │
#   └───────────────────────────────────────────────────────────────────┘
#                               │
#                               ▼
#   ┌───────────────────────┬───────────────────────────────────────────┐
#   │  STEP 3A (parallel)   │  STEP 3B (parallel)                       │
#   │  cross-species        │  bp-pipeline bridge                       │
#   │                       │                                           │
#   │  cross_species/       │  bp_bridge/                               │
#   │    STEP_03A_*.py      │    STEP_03B_*.py                          │
#   │                       │                                           │
#   │  in: catfish-synteny- │  in: phase_6 TSVs                         │
#   │      toolkit outputs  │     (candidate_breakpoints_consensus.tsv, │
#   │  out: synteny_v6.json │      candidate_ancestral_fragments.tsv,   │
#   │  flat keys:           │      candidate_breakpoints_per_method.tsv)│
#   │   q5_bs_*,            │  out: fragment_distribution.json          │
#   │   q5_conservation_*,  │  flat keys:                               │
#   │   q5_tree_polarization│   q3_final_{left,right}_bp,               │
#   │   q3_*_flank_*        │   q3_breakpoint_precision_class,          │
#   │                       │   q3_n_methods_agreeing_{left,right},     │
#   │                       │   q2_interior_class,                      │
#   │                       │   q2_fragment_{left,right}_*,             │
#   │                       │   q2_n_fragment_carriers                  │
#   └───────────────────────┴───────────────────────────────────────────┘
#                               │
#                               ▼
#   ┌───────────────────────────────────────────────────────────────────┐
#   │  STEP 4: final structural-class assignment  (phase_9)             │
#   │    phase_9_classification/STEP_04_assign_structural_class.py       │
#   │    in:  <outdir>/<cid>/mechanism_assembled.json      (from STEP 1) │
#   │         <outdir>/<cid>/bnd_sided_support.json        (from STEP 2) │
#   │         <outdir>/<cid>/synteny_v6.json               (from STEP 3A)│
#   │         <outdir>/<cid>/fragment_distribution.json    (from STEP 3B)│
#   │         <keys_dir>/<cid>/keys.tsv          (per-candidate flat keys)│
#   │    out: <outdir>/final/final_catalog.tsv                           │
#   │         <outdir>/<cid>/final_label.json                            │
#   │         decision: one of 16 structural class labels                │
#   └───────────────────────────────────────────────────────────────────┘
#
# Why STEP 3 is split into 3A + 3B:
#   Both run BEFORE STEP 4 and AFTER STEP 1+2, but neither depends on
#   the other — they consume different inputs and produce independent
#   evidence blocks. They're grouped as "step 3" because they share
#   the same position in the DAG (all-annotation-before-synthesis),
#   not because one is a sub-step of the other.
#
# STEP 3B is the phase_6 → registry bridge. Phase_6 produces three TSVs;
# this script promotes them into registry Tier-2 blocks so STEP 4 can
# read them uniformly with the other evidence blocks.
#
# =============================================================================
# INPUTS (env vars with defaults)
# =============================================================================
#   CANDIDATES          — TSV with candidate_id, chrom, left_bp, right_bp
#   OUTDIR              — where to write per-candidate blocks + final catalog
#   KEYS_DIR            — per-candidate keys.tsv dir (required for STEP 4)
#
#   DELLY_INV_VCF       — MODULE_4D strict catalog (for STEP 1)
#   DELLY_BND_VCF       — MODULE_4E strict catalog (for STEP 1 + 2)
#   MANTA_INV_VCF       — MODULE_4G INV catalog (for STEP 1)
#   MANTA_RAW_VCF       — MODULE_4G raw pre-conversion VCF (for STEP 2)
#
#   BETWEEN_SPECIES_BED — catfish-synteny-toolkit STEP_02 output (for STEP 3A)
#   FLANK_COHERENCE_TSV — catfish-synteny-toolkit STEP_09c output (for STEP 3A)
#   POLARIZED_TSV       — catfish-synteny-toolkit STEP_11 output (for STEP 3A)
#   DOLLO_TSV           — Dollo reconstruction output (optional, for STEP 3A)
#
#   BP_CONSENSUS_TSV    — phase_6 candidate_breakpoints_consensus.tsv      (STEP 3B)
#   BP_FRAGMENTS_TSV    — phase_6 candidate_ancestral_fragments.tsv        (STEP 3B)
#   BP_PER_METHOD_TSV   — phase_6 candidate_breakpoints_per_method.tsv     (STEP 3B, optional)
#
#   REGISTRIES_ROOT     — inversion-popgen-toolkit/registries/ (optional)
#
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
# STEP 4 lives in phase_9_classification (moved pass 15).
PHASE9_DIR="${PHASE9_DIR:-$(cd "${SCRIPT_DIR}/../phase_9_classification" && pwd)}"
EVIDENCE_BLOCKS="$OUTDIR/evidence_blocks"
FINAL="$OUTDIR/final"
mkdir -p "$EVIDENCE_BLOCKS" "$FINAL"

echo "=============================================================================="
echo "phase_8 → phase_9 evidence synthesis"
echo "  candidates       = $CANDIDATES"
echo "  outdir           = $OUTDIR"
echo "  evidence blocks  = $EVIDENCE_BLOCKS"
echo "  final catalog    = $FINAL"
echo "=============================================================================="

# -------------------------------------------------------------
# STEP 1 — assembled-junction forensics
# -------------------------------------------------------------
echo ""
echo "[STEP 1] assembled-junction forensics (q4_mechanism)"
python3 "$SCRIPT_DIR/q4_mechanism/STEP_01_assembled_junction.py" \
    --candidates "$CANDIDATES" \
    --delly_inv_vcf "$DELLY_INV_VCF" \
    --delly_bnd_vcf "$DELLY_BND_VCF" \
    --manta_inv_vcf "$MANTA_INV_VCF" \
    --outdir "$EVIDENCE_BLOCKS" \
    ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}

# -------------------------------------------------------------
# STEP 2 — single-sided BND support
# -------------------------------------------------------------
echo ""
echo "[STEP 2] single-sided BND support (q7_existence_audit)"
python3 "$SCRIPT_DIR/q7_existence_audit/STEP_02_bnd_sided_support.py" \
    --candidates "$CANDIDATES" \
    --delly_bnd_vcf "$DELLY_BND_VCF" \
    --manta_raw_vcf "$MANTA_RAW_VCF" \
    --outdir "$EVIDENCE_BLOCKS" \
    ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}

# -------------------------------------------------------------
# STEP 3A — cross-species conservation bridge
# -------------------------------------------------------------
# Writes Q5 conservation + polarization keys AND Q3 flank coherence keys.
# The Q3 flank keys folded in from the superseded flank_coherence schema
# on 2026-04-24 (v6 consolidation).
echo ""
echo "[STEP 3A] cross-species synteny bridge (cross_species)"
python3 "$SCRIPT_DIR/cross_species/STEP_03A_cross_species_bridge.py" \
    --candidates "$CANDIDATES" \
    --between_species_bed "$BETWEEN_SPECIES_BED" \
    --flank_coherence_tsv "$FLANK_COHERENCE_TSV" \
    --polarized_tsv "$POLARIZED_TSV" \
    ${DOLLO_TSV:+--dollo_tsv "$DOLLO_TSV"} \
    --outdir "$EVIDENCE_BLOCKS" \
    ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}

# -------------------------------------------------------------
# STEP 3B — phase_6 breakpoint-pipeline → registry bridge
# -------------------------------------------------------------
# Reads the three TSVs produced by phase_6_breakpoint_refinement's
# run_pipeline.sh and writes derived flat keys (Q2 interior class +
# Q3 breakpoint precision) into a fragment_distribution Tier-2 block.
# STEP 4 consumes q2_interior_class to decide suffix modifiers and
# terminal gating (complex_rearrangement_out_of_scope).
if [[ -n "$BP_CONSENSUS_TSV" && -n "$BP_FRAGMENTS_TSV" ]]; then
    echo ""
    echo "[STEP 3B] bp-pipeline bridge (bp_bridge → fragment_distribution)"
    python3 "$SCRIPT_DIR/bp_bridge/STEP_03B_bp_pipeline_bridge.py" \
        --consensus_tsv "$BP_CONSENSUS_TSV" \
        --fragments_tsv "$BP_FRAGMENTS_TSV" \
        ${BP_PER_METHOD_TSV:+--per_method_tsv "$BP_PER_METHOD_TSV"} \
        --outdir "$EVIDENCE_BLOCKS" \
        ${REGISTRIES_ROOT:+--registries_root "$REGISTRIES_ROOT"}
else
    echo ""
    echo "[STEP 3B] SKIPPED — set BP_CONSENSUS_TSV + BP_FRAGMENTS_TSV to enable"
    echo "          Without this, STEP 4 will report q2_interior_class=degenerate"
    echo "          and never apply interior suffixes."
fi

# -------------------------------------------------------------
# STEP 4 — final structural-class synthesis (phase_9)
# -------------------------------------------------------------
# Reads all four blocks written by STEPs 1/2/3A/3B plus the per-candidate
# keys.tsv, decides one of 16 structural class labels, writes per-candidate
# final_label.json and cohort final_catalog.tsv.
echo ""
echo "[STEP 4] final structural-class assignment (phase_9)"
python3 "$PHASE9_DIR/STEP_04_assign_structural_class.py" \
    --candidates "$CANDIDATES" \
    --keys_dir "$KEYS_DIR" \
    --evidence_blocks_dir "$EVIDENCE_BLOCKS" \
    --outdir "$FINAL"

echo ""
echo "=============================================================================="
echo "[done] Final catalog: $FINAL/final_catalog.tsv"
echo "=============================================================================="
