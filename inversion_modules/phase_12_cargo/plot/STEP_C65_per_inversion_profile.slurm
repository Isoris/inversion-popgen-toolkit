#!/usr/bin/env bash
#SBATCH --job-name=cargo_C65_fig
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/C65_fig_%j.out
#SBATCH --error=logs/C65_fig_%j.err
set -euo pipefail
mkdir -p logs

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_cargo_config.sh"

TARGET="${1:-all}"

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_C65_per_inversion_profile.R" "${TARGET}"

echo "[DONE] STEP_C65"
