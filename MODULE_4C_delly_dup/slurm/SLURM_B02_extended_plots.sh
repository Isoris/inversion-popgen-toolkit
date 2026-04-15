#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J delly_dup_plots
#SBATCH -o logs/08_extended_plots.%j.out
#SBATCH -e logs/08_extended_plots.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4c_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

DIR_SUMMARY="${OUTDIR}/11_summary"
DIR_PLOTS="${OUTDIR}/12_plots"
mkdir -p "${DIR_PLOTS}"

RSCRIPT_BIN="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript"

"${RSCRIPT_BIN}" "${SCRIPT_DIR}/11_plot_DUP_extended.R" \
  --summary_dir "${DIR_SUMMARY}" \
  --plot_dir "${DIR_PLOTS}" \
  --ref_fai "${REF_FAI}" \
  --samples_unrelated "${SAMPLES_UNRELATED}"
