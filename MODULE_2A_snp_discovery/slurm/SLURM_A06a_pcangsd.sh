#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -J pcangsd
#SBATCH -o pcangsd.%A_%a.out
#SBATCH -e pcangsd.%A_%a.err
set -euo pipefail
source ~/.bashrc; mamba activate assembly
source "$(dirname "$0")/../config.sh"

OUTBASE="${THIN_DIR}/05_pcangsd_byLG"
P="${SLURM_CPUS_PER_TASK:-32}"

L200="${OUTBASE}/beagle_LG_thin_200.list"
L500="${OUTBASE}/beagle_LG_thin_500.list"
L1000="${OUTBASE}/beagle_LG_thin_1000.list"

KS=($(seq $K_MIN $K_MAX)); NK=${#KS[@]}
N200=$(wc -l < "$L200"); N500=$(wc -l < "$L500"); N1000=$(wc -l < "$L1000")

task=${SLURM_ARRAY_TASK_ID}; idx=$((task-1))
file_gi=$(( idx / NK )); k_i=$(( idx % NK )); K=${KS[$k_i]}

if (( file_gi < N200 )); then thin=200; fi=$file_gi; LIST="$L200"
elif (( file_gi < N200+N500 )); then thin=500; fi=$((file_gi-N200)); LIST="$L500"
else thin=1000; fi=$((file_gi-N200-N500)); LIST="$L1000"; fi

IN=$(sed -n "$((fi+1))p" "$LIST")
[[ -n "$IN" && -s "$IN" ]] || exit 1

bn=$(basename "$IN"); tag="${bn%.beagle.gz}"
LG=$(echo "$bn" | grep -oE 'LG[0-9]+' | head -n1); [[ -n "$LG" ]] || LG="LGNA"

OUTDIR="${OUTBASE}/thin_${thin}/${LG}/K${K}"; mkdir -p "$OUTDIR"
OUT="${OUTDIR}/${tag}.pcangsd"

pcangsd --beagle "$IN" --eig ${PCANGSD_EIG} --threads "$P" \
  --maf ${PCANGSD_MAF} --iter ${PCANGSD_ITER} \
  --admix --admix-K "$K" --admix-seed 1 \
  --tree --tree-samples "${SAMPLE_LIST}" --out "$OUT"
