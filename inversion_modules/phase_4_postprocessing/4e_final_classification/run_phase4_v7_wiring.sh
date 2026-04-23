#!/bin/bash
# =============================================================================
# run_phase4_v7_wiring.sh — Phase 4 v7 master orchestrator
# -----------------------------------------------------------------------------
# Sets the environment variables and invokes run_phase4_v5_new_evidence.sh,
# which runs the full v5+v6+v7 chain:
#
#   STEP 1   cheat29b_assembled_junction  → mechanism_assembled.json
#   STEP 2   bnd_sided_support            → bnd_sided_support.json
#   STEP 3   cross_species_bridge_v6      → synteny_v6.json + flank_coherence.json
#   STEP 3b  bp_pipeline_bridge (v7)      → fragment_distribution.json + breakpoint keys
#   STEP 4   assign_structural_class_v7   → final_label.json + final_catalog.tsv
#
# All I/O paths are controlled by env vars; defaults use the Quentin_project_KEEP
# layout. Override any before sourcing.
#
# Usage:
#   sbatch run_phase4_v7_wiring.sh
#   bash   run_phase4_v7_wiring.sh              # interactive / login node
#
# Override example:
#   CANDIDATES=/path/to/my_cands.tsv \
#   OUTDIR=/path/to/v7_out \
#     bash run_phase4_v7_wiring.sh
# =============================================================================
#SBATCH --job-name=phase4_v7_wiring
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase4_v7_%j.out
#SBATCH --error=logs/phase4_v7_%j.err

set -euo pipefail

# ---- Anchor + project config -----------------------------------------------
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
TOOLKIT="${TOOLKIT:-${BASE}/inversion-popgen-toolkit}"

if [[ -f "${TOOLKIT}/00_inversion_config.sh" ]]; then
  # shellcheck source=/dev/null
  source "${TOOLKIT}/00_inversion_config.sh"
fi

# ---- Paths the v5 launcher expects -----------------------------------------
# Required (will hard-fail if unset):
export CANDIDATES="${CANDIDATES:-${BASE}/CANDIDATE_TABLE.tsv}"
export OUTDIR="${OUTDIR:-${BASE}/phase4_v7_blocks}"
export KEYS_DIR="${KEYS_DIR:-${BASE}/registries/evidence/keys}"

# Optional inputs — SV VCFs
export DELLY_INV_VCF="${DELLY_INV_VCF:-${BASE}/MODULE_4D_sv_delly_inv/delly_inv_strict.vcf.gz}"
export DELLY_BND_VCF="${DELLY_BND_VCF:-${BASE}/MODULE_4E_sv_delly_bnd/delly_bnd_strict.vcf.gz}"
export MANTA_INV_VCF="${MANTA_INV_VCF:-${BASE}/MODULE_4G_sv_manta/manta_inv.vcf.gz}"
export MANTA_RAW_VCF="${MANTA_RAW_VCF:-${BASE}/MODULE_4G_sv_manta/manta_raw.vcf.gz}"

# Optional inputs — synteny toolkit outputs
SYNTENY_DIR="${SYNTENY_DIR:-${BASE}/catfish-synteny-toolkit/output}"
export BETWEEN_SPECIES_BED="${BETWEEN_SPECIES_BED:-${SYNTENY_DIR}/step_02/events.bed}"
export FLANK_COHERENCE_TSV="${FLANK_COHERENCE_TSV:-${SYNTENY_DIR}/step_09c/flank_coherence.tsv}"
export POLARIZED_TSV="${POLARIZED_TSV:-${SYNTENY_DIR}/step_11/polarized.tsv}"
export DOLLO_TSV="${DOLLO_TSV:-${BASE}/dollo_output.tsv}"

# Optional inputs — breakpoint_pipeline outputs (v7 addition)
BP_OUT="${BP_OUT:-${BASE}/breakpoint_pipeline_output}"
export BP_CONSENSUS_TSV="${BP_CONSENSUS_TSV:-${BP_OUT}/candidate_breakpoints_consensus.tsv}"
export BP_FRAGMENTS_TSV="${BP_FRAGMENTS_TSV:-${BP_OUT}/candidate_ancestral_fragments.tsv}"
export BP_PER_METHOD_TSV="${BP_PER_METHOD_TSV:-${BP_OUT}/candidate_breakpoints_per_method.tsv}"

# Registry root
export REGISTRIES_ROOT="${REGISTRIES_ROOT:-${TOOLKIT}/registries}"

# ---- Resolve launcher location ---------------------------------------------
P4="${TOOLKIT}/inversion_modules/phase_4_postprocessing"
LAUNCHER="${P4}/4d_group_dependent/run_phase4_v5_new_evidence.sh"

if [[ ! -f "${LAUNCHER}" ]]; then
  echo "ERROR: launcher not found: ${LAUNCHER}" >&2
  exit 1
fi

mkdir -p "${OUTDIR}" logs

# ---- Pre-flight ------------------------------------------------------------
echo "================================================================"
echo "  Phase 4 v7 wiring — delegating to v5 launcher"
echo "  Started:        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  BASE:           ${BASE}"
echo "  CANDIDATES:     ${CANDIDATES}"
echo "  OUTDIR:         ${OUTDIR}"
echo "  Launcher:       ${LAUNCHER}"
echo ""
echo "  Optional inputs (yes = file exists):"
printf "    %-22s %s\n" "DELLY_INV_VCF:"       "$(test -f "${DELLY_INV_VCF}"       && echo yes || echo no)"
printf "    %-22s %s\n" "DELLY_BND_VCF:"       "$(test -f "${DELLY_BND_VCF}"       && echo yes || echo no)"
printf "    %-22s %s\n" "MANTA_INV_VCF:"       "$(test -f "${MANTA_INV_VCF}"       && echo yes || echo no)"
printf "    %-22s %s\n" "MANTA_RAW_VCF:"       "$(test -f "${MANTA_RAW_VCF}"       && echo yes || echo no)"
printf "    %-22s %s\n" "BETWEEN_SPECIES_BED:" "$(test -f "${BETWEEN_SPECIES_BED}" && echo yes || echo no)"
printf "    %-22s %s\n" "FLANK_COHERENCE_TSV:" "$(test -f "${FLANK_COHERENCE_TSV}" && echo yes || echo no)"
printf "    %-22s %s\n" "POLARIZED_TSV:"       "$(test -f "${POLARIZED_TSV}"       && echo yes || echo no)"
printf "    %-22s %s\n" "DOLLO_TSV:"           "$(test -f "${DOLLO_TSV}"           && echo yes || echo no)"
printf "    %-22s %s\n" "BP_CONSENSUS_TSV:"    "$(test -f "${BP_CONSENSUS_TSV}"    && echo yes || echo no)"
printf "    %-22s %s\n" "BP_FRAGMENTS_TSV:"    "$(test -f "${BP_FRAGMENTS_TSV}"    && echo yes || echo no)"
echo "================================================================"
echo ""

# ---- Run the launcher ------------------------------------------------------
bash "${LAUNCHER}"

# ---- Summary ---------------------------------------------------------------
echo ""
echo "================================================================"
echo "  Phase 4 v7 wiring complete"
echo "  Finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  Final catalog: ${OUTDIR}/final/final_catalog.tsv"
echo "================================================================"

if [[ -f "${OUTDIR}/final/final_catalog.tsv" ]]; then
  n_cands=$(( $(wc -l < "${OUTDIR}/final/final_catalog.tsv") - 1 ))
  echo ""
  echo "Candidates: ${n_cands}"
  echo ""
  echo "Label distribution:"
  awk -F'\t' 'NR > 1 { print $2 }' "${OUTDIR}/final/final_catalog.tsv" | sort | uniq -c | sort -rn
fi
