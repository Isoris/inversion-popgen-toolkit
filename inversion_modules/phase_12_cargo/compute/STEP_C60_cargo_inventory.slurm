#!/usr/bin/env bash
#SBATCH --job-name=cargo_C60_inv
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/C60_inventory_%j.out
#SBATCH --error=logs/C60_inventory_%j.err
set -euo pipefail
mkdir -p logs

source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../00_cargo_config.sh"

# Optional: pass --candidate-id <CID> through "$@"
python3 "${SCRIPT_DIR}/STEP_C60_cargo_inventory.py" "$@"

echo "[DONE] STEP_C60"
echo ""
echo "Diagnostic table summary:"
awk -F'\t' '
  NR==1 {next}
  {n_total++}
  $13=="1" {n_l1++}
  $14=="1" {n_l2++}
  END {
    printf "  Total candidates:    %d\n", n_total
    printf "  Level 1 (burden):    %d\n", n_l1
    printf "  Level 2 (config MI): %d\n", n_l2
  }
' "${CARGO_INVENTORY_DIR}/diagnostic_table.tsv" 2>/dev/null || true
