#!/usr/bin/env bash
#SBATCH --job-name=cargo_extra_plots
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/extra_plots_%j.out
#SBATCH --error=logs/extra_plots_%j.err
set -uo pipefail   # NB: no -e; we want to continue past plots that fail/skip
mkdir -p logs

# =============================================================================
# run_all_extra_plots.slurm — driver for phase_12_cargo/extra_plots
#
# Runs all 12 PLOT_* scripts in sequence. Each script self-detects missing
# inputs and exits 0 with a [skip] message, so the whole run finishes even
# when only a subset of upstream pipelines have completed.
#
# Two of the per-candidate plots (PLOT_01) iterate over all candidates by
# default. Override with PLOT01_TARGET=<cid> to run a single one for testing.
#
# Usage:
#   sbatch launchers/run_all_extra_plots.slurm
#   sbatch launchers/run_all_extra_plots.slurm --groups-from /path/to/groups.tsv
# =============================================================================

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_extra_plots_config.sh"

EXTRA_ARGS=("$@")  # forwarded to plots that accept --groups-from

run_plot() {
  local label="$1"; shift
  local script="$1"; shift
  echo ""
  echo "============================================================"
  echo "[$(date)] $label"
  echo "============================================================"
  if [[ ! -f "$script" ]]; then
    echo "  [skip] $script not found"
    return 0
  fi
  "${RSCRIPT_BIN}" "$script" "$@" || {
    echo "  [warn] $label exited non-zero (continuing with rest)"
  }
}

CMP="${SCRIPT_DIR}/../compute"

# ── Per-candidate plot — PLOT_01 ───────────────────────────────────────────
# By default runs all candidates; override with PLOT01_TARGET to do one.
run_plot "PLOT_01 karyotype network"        "${CMP}/PLOT_01_karyotype_network.R" \
    "${PLOT01_TARGET:-all}"

# ── Cross-candidate plot — PLOT_02 ─────────────────────────────────────────
run_plot "PLOT_02 inversion sharing across groups" \
    "${CMP}/PLOT_02_inversion_sharing_across_groups.R" "${EXTRA_ARGS[@]}"

# ── Variant landscape plots — PLOT_03 to PLOT_12 ───────────────────────────
run_plot "PLOT_03 variant burden dashboard" "${CMP}/PLOT_03_variant_burden_dashboard.R" "${EXTRA_ARGS[@]}"
run_plot "PLOT_04 cohort SFS"               "${CMP}/PLOT_04_cohort_sfs.R"
run_plot "PLOT_05 top recurrent variants"   "${CMP}/PLOT_05_top_recurrent_variants.R"
run_plot "PLOT_06 private vs shared by group" "${CMP}/PLOT_06_private_vs_shared_by_group.R" "${EXTRA_ARGS[@]}"
run_plot "PLOT_07 group exclusivity heatmap" "${CMP}/PLOT_07_group_exclusivity_heatmap.R" "${EXTRA_ARGS[@]}"
run_plot "PLOT_08 upset variant sharing"    "${CMP}/PLOT_08_upset_variant_sharing.R" "${EXTRA_ARGS[@]}"
run_plot "PLOT_09 genomic distribution"     "${CMP}/PLOT_09_genomic_distribution.R"
run_plot "PLOT_10 context enrichment"       "${CMP}/PLOT_10_context_enrichment.R"
run_plot "PLOT_11 SVs in inversions"        "${CMP}/PLOT_11_sv_in_inversions.R"
run_plot "PLOT_12 inversion population frequencies" "${CMP}/PLOT_12_inversion_pop_frequencies.R"

echo ""
echo "============================================================"
echo "[$(date)] Done — outputs under ${EXTRAS_FIG_DIR}"
echo "============================================================"
ls -1 "${EXTRAS_FIG_DIR}/" 2>/dev/null | sed 's/^/  /' || true
