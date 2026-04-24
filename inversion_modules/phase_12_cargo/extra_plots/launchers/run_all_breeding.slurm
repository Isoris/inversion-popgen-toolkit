#!/usr/bin/env bash
#SBATCH --job-name=cargo_breeding
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/breeding_%j.out
#SBATCH --error=logs/breeding_%j.err
set -uo pipefail
mkdir -p logs

# =============================================================================
# run_all_breeding.slurm — driver for phase_12_cargo/extra_plots/breeding
#
# Runs all 3 BREEDING_* scripts. BREEDING_A defaults to --all mode (all-vs-all
# pair ranking). To run a specific pair instead:
#
#   sbatch run_all_breeding.slurm --pair SAMPLE_X SAMPLE_Y
#
# Forwards extra args (e.g. --groups-from) to scripts that accept them.
# =============================================================================

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_extra_plots_config.sh"

EXTRA_ARGS=("$@")
BR="${SCRIPT_DIR}/../breeding"

run_breed() {
  local label="$1"; shift
  local script="$1"; shift
  echo ""
  echo "============================================================"
  echo "[$(date)] $label"
  echo "============================================================"
  if [[ ! -f "$script" ]]; then echo "  [skip] $script not found"; return 0; fi
  "${RSCRIPT_BIN}" "$script" "$@" || echo "  [warn] $label exited non-zero"
}

run_breed "BREEDING_A broodstock compatibility"  "${BR}/BREEDING_A_broodstock_compatibility.R" "${EXTRA_ARGS[@]}"
run_breed "BREEDING_C founder haplotype tracker" "${BR}/BREEDING_C_founder_haplotype_tracker.R" "${EXTRA_ARGS[@]}"
run_breed "BREEDING_D recombination atlas"       "${BR}/BREEDING_D_recombination_atlas.R"

echo ""; echo "[$(date)] Done — outputs under ${EXTRAS_TBL_DIR}/ and ${EXTRAS_FIG_DIR}/"
ls -1 "${EXTRAS_TBL_DIR}"/BREEDING_*.tsv 2>/dev/null | sed 's/^/  /' || true
ls -1 "${EXTRAS_FIG_DIR}"/BREEDING_*.pdf 2>/dev/null | sed 's/^/  /' || true
