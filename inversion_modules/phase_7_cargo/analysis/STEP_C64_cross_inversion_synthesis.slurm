#!/usr/bin/env bash
#SBATCH --job-name=cargo_C64_synth
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/C64_synth_%j.out
#SBATCH --error=logs/C64_synth_%j.err
set -euo pipefail
mkdir -p logs

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_cargo_config.sh"

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/STEP_C64_cross_inversion_synthesis.R"

echo "[DONE] STEP_C64"
