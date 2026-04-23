#!/usr/bin/env bash
#SBATCH --job-name=cargo_C63_enrich
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/C63_enrich_%j.out
#SBATCH --error=logs/C63_enrich_%j.err
set -euo pipefail
mkdir -p logs

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_cargo_config.sh"

TARGET="${1:-all}"

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_C63_functional_enrichment.R" "${TARGET}"

echo "[DONE] STEP_C63"
