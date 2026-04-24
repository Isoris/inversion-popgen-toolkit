#!/usr/bin/env bash
#SBATCH --job-name=cargo_extra_tables
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/extra_tables_%j.out
#SBATCH --error=logs/extra_tables_%j.err
set -uo pipefail
mkdir -p logs

# =============================================================================
# run_all_extra_tables.slurm — driver for phase_12_cargo/extra_plots/tables
#
# Runs all 6 TABLE_* scripts in sequence, skipping any with missing inputs.
#
# Usage:
#   sbatch launchers/run_all_extra_tables.slurm
#   sbatch launchers/run_all_extra_tables.slurm --groups-from /path/to/groups.tsv
# =============================================================================

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_extra_plots_config.sh"

EXTRA_ARGS=("$@")
TBL="${SCRIPT_DIR}/../tables"

run_table() {
  local label="$1"; shift
  local script="$1"; shift
  echo ""
  echo "============================================================"
  echo "[$(date)] $label"
  echo "============================================================"
  if [[ ! -f "$script" ]]; then
    echo "  [skip] $script not found"; return 0
  fi
  "${RSCRIPT_BIN}" "$script" "$@" || {
    echo "  [warn] $label exited non-zero (continuing)"
  }
}

run_table "TABLE_01 top recurrent w/ group AF"  "${TBL}/TABLE_01_top_recurrent_with_group_AF.R" "${EXTRA_ARGS[@]}"
run_table "TABLE_02 top private singletons"     "${TBL}/TABLE_02_top_private_singletons.R"      "${EXTRA_ARGS[@]}"
run_table "TABLE_03 variant categories summary" "${TBL}/TABLE_03_variant_categories_summary.R"
run_table "TABLE_04 group exclusivity summary"  "${TBL}/TABLE_04_group_exclusivity_summary.R"   "${EXTRA_ARGS[@]}"
run_table "TABLE_05 repeat class enrichment"    "${TBL}/TABLE_05_repeat_class_enrichment.R"
run_table "TABLE_06 top inversion carriers"     "${TBL}/TABLE_06_top_inversion_carriers.R"

echo ""
echo "============================================================"
echo "[$(date)] Done — outputs under ${EXTRAS_TBL_DIR}"
echo "============================================================"
ls -1 "${EXTRAS_TBL_DIR}"/TABLE_*.tsv 2>/dev/null | sed 's/^/  /' || true
