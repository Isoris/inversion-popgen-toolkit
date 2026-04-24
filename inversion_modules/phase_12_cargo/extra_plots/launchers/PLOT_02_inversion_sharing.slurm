#!/usr/bin/env bash
#SBATCH --job-name=cargo_plot02_inv_sharing
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/plot02_inv_sharing_%j.out
#SBATCH --error=logs/plot02_inv_sharing_%j.err
set -euo pipefail
mkdir -p logs

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_extra_plots_config.sh"

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/../compute/PLOT_02_inversion_sharing_across_groups.R" "$@"

echo "[DONE] PLOT_02"
