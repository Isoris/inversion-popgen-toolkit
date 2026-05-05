#!/bin/bash
#SBATCH --job-name=sv_evidence_pipeline
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/sv_evidence_%j.out
#SBATCH --error=logs/sv_evidence_%j.err
#
# Atlas SV evidence producers — LANTA SLURM job
# =============================================
# Produces all three layers for one candidate at a time:
#   1) sv_genotype_counts.json         (STEP_SV_GT_AGG)
#   2) sv_evidence_combinations.json   (STEP_SV_EVID_COMB)
#   3) sv_support_by_sample.json       (STEP_SV_SUPPORT)  — step-6 heatmap
# Output goes into the per-candidate folder under data/<chrom>/candidates/<id>/
# which the local atlas page can drag-drop directly.
#
# Submit one candidate at a time:
#   sbatch run_sv_evidence_pipeline.slurm INV_LG28_003
# Or batch-submit a whole chromosome's candidates:
#   for cid in INV_LG28_001 INV_LG28_002 INV_LG28_003 ...; do
#     sbatch run_sv_evidence_pipeline.slurm "$cid"
#   done

set -euo pipefail

# -----------------------------------------------------------------------------
# Paths — edit if the project moves
# -----------------------------------------------------------------------------
PROJECT_ROOT="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
PRODUCERS="${PROJECT_ROOT}/inversion_codebase_v8.5/atlas_producers/sv_evidence"
INPUTS="${PROJECT_ROOT}/MS_Inversions_inputs"
OUT_ROOT="${PROJECT_ROOT}/atlas_outputs/sv_evidence"

# -----------------------------------------------------------------------------
# Per-candidate inputs
# -----------------------------------------------------------------------------
CANDIDATE_ID="${1:?usage: sbatch run_sv_evidence_pipeline.slurm <candidate_id>}"
CANDIDATE_JSON="${INPUTS}/candidates/${CANDIDATE_ID}.json"
KARYOTYPE_TSV="${INPUTS}/karyotype/${CANDIDATE_ID}.locked_labels.tsv"
SV_VCF="${INPUTS}/sv_calls/${CANDIDATE_ID}.merged.vcf"
SV_TSV="${INPUTS}/sv_calls/${CANDIDATE_ID}.calls.tsv"
EVIDENCE_TSV="${INPUTS}/bam_evidence/${CANDIDATE_ID}.evidence_matrix.tsv"
EVIDENCE_TYPES_JSON="${INPUTS}/bam_evidence/evidence_types.json"

# Pull chrom from the candidate JSON for the output folder path
CHROM=$(python3 -c "import json,sys; print(json.load(open('${CANDIDATE_JSON}'))['chrom'])")

echo "[$(date)] Starting SV evidence pipeline for ${CANDIDATE_ID} on ${CHROM}"

# -----------------------------------------------------------------------------
# Environment
# -----------------------------------------------------------------------------
module load Python/3.11
# Quentin's pinned local env (matches the rest of the codebase)
source "${PROJECT_ROOT}/.venv/bin/activate" 2>/dev/null || true

mkdir -p logs "${OUT_ROOT}"

# -----------------------------------------------------------------------------
# Step 1: genotype counts
# -----------------------------------------------------------------------------
python3 "${PRODUCERS}/STEP_SV_GT_AGG_aggregate_genotype_counts.py" \
    --vcf       "${SV_VCF}" \
    --candidate "${CANDIDATE_JSON}" \
    --karyotype "${KARYOTYPE_TSV}" \
    --out-root  "${OUT_ROOT}" \
    --fdr-cutoff 0.05

# -----------------------------------------------------------------------------
# Step 2: evidence combinations (UpSet panel)
# -----------------------------------------------------------------------------
if [[ -f "${EVIDENCE_TSV}" ]]; then
    python3 "${PRODUCERS}/STEP_SV_EVID_COMB_emit_combinations.py" \
        --evidence       "${EVIDENCE_TSV}" \
        --evidence-types "${EVIDENCE_TYPES_JSON}" \
        --candidate-id   "${CANDIDATE_ID}" \
        --chrom          "${CHROM}" \
        --out-root       "${OUT_ROOT}" \
        --top-n          20
else
    echo "[$(date)] no evidence matrix found at ${EVIDENCE_TSV} — skipping STEP_SV_EVID_COMB"
fi

# -----------------------------------------------------------------------------
# Step 3: per-sample × per-SV dosage matrix (sv_support_by_sample_v1)
# -----------------------------------------------------------------------------
# Powers the step-6 heatmap and unlocks three deferred features:
#   - selection-scoped genotype-count recount in the SV table
#   - dim non-carriers in locus glyphs when an UpSet bar is active
#   - per-cell hover tooltip on the heatmap
if [[ -f "${SV_TSV}" ]]; then
    python3 "${PRODUCERS}/STEP_SV_SUPPORT_emit_support_by_sample.py" \
        --tsv        "${SV_TSV}" \
        --candidate  "${CANDIDATE_JSON}" \
        --karyotype  "${KARYOTYPE_TSV}" \
        --out-root   "${OUT_ROOT}"
else
    echo "[$(date)] no calls TSV found at ${SV_TSV} — skipping STEP_SV_SUPPORT"
fi

echo "[$(date)] pipeline complete; output at ${OUT_ROOT}/${CHROM}/candidates/${CANDIDATE_ID}/"
ls -la "${OUT_ROOT}/${CHROM}/candidates/${CANDIDATE_ID}/" || true
