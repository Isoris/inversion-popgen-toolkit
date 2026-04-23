#!/usr/bin/env bash
#SBATCH --job-name=cargo_plot03_burden
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/plot03_burden_%j.out
#SBATCH --error=logs/plot03_burden_%j.err
set -euo pipefail
mkdir -p logs

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_extra_plots_config.sh"

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/../compute/PLOT_03_variant_burden_dashboard.R" "$@"

echo "[DONE] PLOT_03"
