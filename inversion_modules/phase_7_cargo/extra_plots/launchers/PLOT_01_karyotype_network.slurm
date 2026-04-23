#!/usr/bin/env bash
#SBATCH --job-name=cargo_plot01_kary_net
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/plot01_kary_net_%j.out
#SBATCH --error=logs/plot01_kary_net_%j.err
set -euo pipefail
mkdir -p logs

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_extra_plots_config.sh"

TARGET="${1:-all}"

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/../compute/PLOT_01_karyotype_network.R" "${TARGET}"

echo "[DONE] PLOT_01"
