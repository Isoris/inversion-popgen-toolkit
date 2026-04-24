#!/bin/bash
# =============================================================================
# cheat14_self_align.sh — minimap2 self-alignment for repeat architecture
#
# Usage: bash cheat14_self_align.sh <ref_fasta> <chr> <bp1> <bp2> <outdir>
#   Extracts ±50 kb flanks around each breakpoint, runs self-alignment,
#   and produces PAF + optional dotplot.
# =============================================================================
set -euo pipefail

REF_FASTA="${1:?Usage: $0 <ref_fasta> <chr> <bp1> <bp2> <outdir>}"
CHR="${2:?}"
BP1="${3:?}"
BP2="${4:?}"
OUTDIR="${5:-.}"
FLANK=50000

mkdir -p "${OUTDIR}"

# ── Extract flanks ──────────────────────────────────────────────
LEFT1=$((BP1 > FLANK ? BP1 - FLANK : 1))
RIGHT1=$((BP1 + FLANK))
LEFT2=$((BP2 > FLANK ? BP2 - FLANK : 1))
RIGHT2=$((BP2 + FLANK))

FLANK_FA="${OUTDIR}/flanks_${CHR}_${BP1}_${BP2}.fa"
samtools faidx "${REF_FASTA}" \
  "${CHR}:${LEFT1}-${BP1}" \
  "${CHR}:${BP1}-${RIGHT1}" \
  "${CHR}:${LEFT2}-${BP2}" \
  "${CHR}:${BP2}-${RIGHT2}" \
  > "${FLANK_FA}"

echo "[cheat14_sh] Extracted flanks → ${FLANK_FA}"

# ── Self-alignment ─────────────────────────────────────────────
PAF="${OUTDIR}/self_align_${CHR}_${BP1}_${BP2}.paf"
minimap2 -X -c --eqx "${FLANK_FA}" "${FLANK_FA}" > "${PAF}" 2>/dev/null

N_HITS=$(wc -l < "${PAF}")
N_INV=$(awk '$5 == "-"' "${PAF}" | wc -l)
echo "[cheat14_sh] Alignment: ${N_HITS} hits (${N_INV} inverted)"

# ── Cross-alignment (left_bp flanks vs right_bp flanks) ────────
CROSS_FA_Q="${OUTDIR}/cross_q_${CHR}.fa"
CROSS_FA_T="${OUTDIR}/cross_t_${CHR}.fa"
samtools faidx "${REF_FASTA}" \
  "${CHR}:${LEFT1}-${RIGHT1}" > "${CROSS_FA_Q}"
samtools faidx "${REF_FASTA}" \
  "${CHR}:${LEFT2}-${RIGHT2}" > "${CROSS_FA_T}"

CROSS_PAF="${OUTDIR}/cross_align_${CHR}_${BP1}_${BP2}.paf"
minimap2 -c --eqx "${CROSS_FA_T}" "${CROSS_FA_Q}" > "${CROSS_PAF}" 2>/dev/null

N_CROSS=$(wc -l < "${CROSS_PAF}")
N_CROSS_INV=$(awk '$5 == "-"' "${CROSS_PAF}" | wc -l)
echo "[cheat14_sh] Cross-alignment: ${N_CROSS} hits (${N_CROSS_INV} inverted)"

# ── Summary ────────────────────────────────────────────────────
SUMMARY="${OUTDIR}/repeat_summary_${CHR}_${BP1}_${BP2}.tsv"
echo -e "chr\tbp1\tbp2\tself_hits\tself_inv\tcross_hits\tcross_inv" > "${SUMMARY}"
echo -e "${CHR}\t${BP1}\t${BP2}\t${N_HITS}\t${N_INV}\t${N_CROSS}\t${N_CROSS_INV}" >> "${SUMMARY}"

echo "[cheat14_sh] Done → ${SUMMARY}"
