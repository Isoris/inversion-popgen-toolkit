#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=128G
#SBATCH -t 1-00:00:00
#SBATCH -J manta_filter_all
#SBATCH -o logs/filter_all_sv_types.%j.out
#SBATCH -e logs/filter_all_sv_types.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_4G_ALL_Manta
bash slurm/SLURM_A04_annotation_filter.sh
