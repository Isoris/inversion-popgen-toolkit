#!/usr/bin/env bash
# =============================================================================
# 00_module1_config.sh — Configuration for MODULE_1 (reads + BAMs)
# =============================================================================
# This module covers: fastp QC, species verification (Mash + meryl),
# minimap2 alignment, merge/markdup/clipOverlap, TLEN stats, popgen filter,
# depth QC, provenance manifest.
#
# Two tracks:
#   Track A (STEP_A01..A08 + SLURM_A01/A04) — read-level QC and species check
#   Track B (STEP_B01..B08 + SLURM_B01/B02/B05/B06) — BAM processing
#
# Source from any pipeline script:
#   source "$(dirname "$0")/../00_module1_config.sh"
# =============================================================================

# ── Project root ────────────────────────────────────────────────────────────
PROJECT="${PROJECT:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
BASE="${PROJECT}"

# ── Reference genome ────────────────────────────────────────────────────────
REF="${PROJECT}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"

# ── Module directories ──────────────────────────────────────────────────────
MODULE_DIR="${PROJECT}/Modules/MODULE_1_read_prep"
STEPS_DIR="${MODULE_DIR}/steps"
SLURM_DIR="${MODULE_DIR}/slurm"
LAUNCHERS_DIR="${MODULE_DIR}/launchers"

# ── Working directories ─────────────────────────────────────────────────────
FASTQ_DIR="${PROJECT}/01_raw_reads"
FASTP_DIR="${PROJECT}/02_fastp_trimmed"
BAM_RAW_DIR="${PROJECT}/03_bam_minimap2_raw"
BAM_MERGED_DIR="${PROJECT}/04_bam_merged_markdup"
BAM_POPGEN_DIR="${PROJECT}/05_bam_popgen_filtered"
QC_DIR="${PROJECT}/06_qc_reports"

# ── Key parameters ──────────────────────────────────────────────────────────
# Read QC (fastp)
FASTP_MIN_QUAL=20           # Phred ≥ 20
FASTP_TRIM_FRONT=5          # 5 bp trimmed from front of both reads
FASTP_TRIM_TAIL=5           # 5 bp trimmed from tail of both reads
FASTP_MIN_LEN=30            # minimum read length after trimming

# Species verification (Mash)
MASH_K=31
MASH_SKETCH=50000
MASH_MARGIN=0.002

# Alignment (minimap2)
MM2_PRESET="-ax sr -Y"

# Insert size cutoff (empirical p99, this library)
TLEN_P99=514

# PopGen BAM filter
POPGEN_MAPQ_MIN=60
POPGEN_FLAG_REQUIRED="0x2"    # proper pair
POPGEN_FLAG_FORBIDDEN="0xF0C" # no unmapped/secondary/dup/supp
POPGEN_TLEN_MAX="${TLEN_P99}"

# Depth QC (mosdepth)
DEPTH_THRESHOLDS="1,5,10"
CALLABLE_GENOME_BP=963905721

# ── Manifest output paths ───────────────────────────────────────────────────
SAMPLE_MANIFEST="${PROJECT}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
FASTP_QC_TABLE="${QC_DIR}/fastp_qc_table.tsv"
BAM_QC_TABLE="${QC_DIR}/all_samples.bam_qc.tsv"
COVERAGE_CONTEXT_QC="${QC_DIR}/all_samples.coverage_context_qc.tsv"

# ── SLURM defaults ──────────────────────────────────────────────────────────
SLURM_ACCOUNT="${SLURM_ACCOUNT:-lt200308}"
SLURM_PARTITION="${SLURM_PARTITION:-compute}"

export PROJECT BASE REF REF_FAI
export MODULE_DIR STEPS_DIR SLURM_DIR LAUNCHERS_DIR
export FASTQ_DIR FASTP_DIR BAM_RAW_DIR BAM_MERGED_DIR BAM_POPGEN_DIR QC_DIR
export FASTP_MIN_QUAL FASTP_TRIM_FRONT FASTP_TRIM_TAIL FASTP_MIN_LEN
export MASH_K MASH_SKETCH MASH_MARGIN MM2_PRESET TLEN_P99
export POPGEN_MAPQ_MIN POPGEN_FLAG_REQUIRED POPGEN_FLAG_FORBIDDEN POPGEN_TLEN_MAX
export DEPTH_THRESHOLDS CALLABLE_GENOME_BP
export SAMPLE_MANIFEST FASTP_QC_TABLE BAM_QC_TABLE COVERAGE_CONTEXT_QC
export SLURM_ACCOUNT SLURM_PARTITION
