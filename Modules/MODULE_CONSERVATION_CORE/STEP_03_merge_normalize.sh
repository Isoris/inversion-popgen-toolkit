#!/usr/bin/env bash
#SBATCH --job-name=CONS_03_merge_norm
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/03_merge_norm_%j.out
#SBATCH --error=slurm_logs/03_merge_norm_%j.err
set -euo pipefail

# =============================================================================
# STEP 03 — Merge & Normalize All Variant Sources
# =============================================================================
# Merges Clair3 SNPs/INDELs + DELLY SVs + Manta SVs into unified VCFs.
# Normalizes small variants with bcftools norm.
# SVs kept separate (not normalizable the same way).
# Produces per-chromosome VCFs ready for annotation.
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

OUTDIR="${DIR_VARIANTS}"
mkdir -p "${OUTDIR}/clair3_merged" "${OUTDIR}/delly_merged" "${OUTDIR}/manta_merged" \
         "${OUTDIR}/normalized" "${OUTDIR}/sv_combined" "${DIR_LOGS}"

# ─────────────────────────────────────────
# A. CLAIR3: merge per-sample VCFs into multi-sample per-chrom
# ─────────────────────────────────────────
echo "=== A. Merging Clair3 per-sample VCFs ==="

for CHR in "${CHROMS[@]}"; do
    MERGED="${OUTDIR}/clair3_merged/${CHR}.clair3.merged.vcf.gz"
    NORM="${OUTDIR}/normalized/${CHR}.clair3.norm.vcf.gz"

    if [[ -f "${NORM}" ]]; then
        echo "  [SKIP] ${CHR} — normalized VCF exists"
        continue
    fi

    echo "  [MERGE] ${CHR}..."

    # Collect per-sample VCFs for this chromosome
    VCF_LIST="${OUTDIR}/clair3_merged/${CHR}.vcf_list.txt"
    > "${VCF_LIST}"
    while read -r SAMPLE; do
        vcf="${CLAIR3_DIR}/vcf/${CHR}/${SAMPLE}.${CHR}.clair3.discovery.vcf.gz"
        if [[ -f "${vcf}" ]]; then
            echo "${vcf}" >> "${VCF_LIST}"
        fi
    done < "${SAMPLES_ALL}"

    N_VCF=$(wc -l < "${VCF_LIST}")
    echo "    Found ${N_VCF} sample VCFs"

    if [[ "${N_VCF}" -eq 0 ]]; then
        echo "    WARNING: No VCFs for ${CHR} — skipping"
        continue
    fi

    # Merge with bcftools (handles overlapping records)
    bcftools merge \
        --file-list "${VCF_LIST}" \
        --merge none \
        --threads 4 \
        -Oz -o "${MERGED}"
    bcftools index -t "${MERGED}"

    # Normalize: split multiallelics + left-align
    echo "  [NORM] ${CHR}..."
    bcftools norm \
        -m -both \
        -f "${REF_FASTA}" \
        --threads 4 \
        -Oz -o "${NORM}" \
        "${MERGED}"
    bcftools index -t "${NORM}"

    N_VAR=$(bcftools view -H "${NORM}" | wc -l)
    echo "    ${CHR}: ${N_VAR} normalized variants"
done

# ─────────────────────────────────────────
# B. DELLY: combine SV types into single per-chrom VCF
# ─────────────────────────────────────────
echo ""
echo "=== B. Combining DELLY SV VCFs ==="

# Find DELLY merged/filtered VCFs across SV types
for SV_TYPE in DEL DUP INV BND; do
    SV_DIR="${BASE}/MODULE_5_DELLY_${SV_TYPE}"
    # Look for the merged/filtered population VCF
    SV_VCF=$(find "${SV_DIR}" -name "*merged*filtered*.vcf.gz" -o -name "*population*.vcf.gz" 2>/dev/null | head -1)
    if [[ -n "${SV_VCF}" && -f "${SV_VCF}" ]]; then
        echo "  Found DELLY ${SV_TYPE}: ${SV_VCF}"
        cp "${SV_VCF}" "${OUTDIR}/delly_merged/delly_${SV_TYPE}.vcf.gz"
        bcftools index -t "${OUTDIR}/delly_merged/delly_${SV_TYPE}.vcf.gz" 2>/dev/null || true
    else
        echo "  WARNING: No merged DELLY ${SV_TYPE} VCF found in ${SV_DIR}"
    fi
done

# INS failed per handoff — skip
echo "  [NOTE] DELLY INS skipped (failed per project notes)"

# Concat all DELLY SV types
DELLY_FILES=()
for SV_TYPE in DEL DUP INV BND; do
    f="${OUTDIR}/delly_merged/delly_${SV_TYPE}.vcf.gz"
    [[ -f "$f" ]] && DELLY_FILES+=("$f")
done

if [[ ${#DELLY_FILES[@]} -gt 0 ]]; then
    echo "  Concatenating ${#DELLY_FILES[@]} DELLY VCFs..."
    bcftools concat \
        "${DELLY_FILES[@]}" \
        --allow-overlaps \
        -Oz -o "${OUTDIR}/delly_merged/delly_all_sv.vcf.gz"
    bcftools sort \
        -Oz -o "${OUTDIR}/sv_combined/delly_all_sv.sorted.vcf.gz" \
        "${OUTDIR}/delly_merged/delly_all_sv.vcf.gz"
    bcftools index -t "${OUTDIR}/sv_combined/delly_all_sv.sorted.vcf.gz"
fi

# ─────────────────────────────────────────
# C. MANTA: combine into single VCF
# ─────────────────────────────────────────
echo ""
echo "=== C. Processing Manta SVs ==="

if [[ -d "${MANTA_DIR}" ]]; then
    MANTA_VCF=$(find "${MANTA_DIR}" -name "*diploidSV*vcf.gz" -o -name "*candidateSV*vcf.gz" 2>/dev/null | head -1)
    if [[ -n "${MANTA_VCF}" && -f "${MANTA_VCF}" ]]; then
        echo "  Found Manta VCF: ${MANTA_VCF}"
        bcftools sort \
            -Oz -o "${OUTDIR}/manta_merged/manta_all_sv.sorted.vcf.gz" \
            "${MANTA_VCF}"
        bcftools index -t "${OUTDIR}/manta_merged/manta_all_sv.sorted.vcf.gz"
    fi
else
    echo "  WARNING: Manta directory not found: ${MANTA_DIR}"
fi

# ─────────────────────────────────────────
# D. Combined SV set (DELLY + Manta, deduplicated)
# ─────────────────────────────────────────
echo ""
echo "=== D. Combining DELLY + Manta SVs ==="

SV_FILES=()
[[ -f "${OUTDIR}/sv_combined/delly_all_sv.sorted.vcf.gz" ]] && \
    SV_FILES+=("${OUTDIR}/sv_combined/delly_all_sv.sorted.vcf.gz")
[[ -f "${OUTDIR}/manta_merged/manta_all_sv.sorted.vcf.gz" ]] && \
    SV_FILES+=("${OUTDIR}/manta_merged/manta_all_sv.sorted.vcf.gz")

if [[ ${#SV_FILES[@]} -gt 0 ]]; then
    bcftools concat \
        "${SV_FILES[@]}" \
        --allow-overlaps \
        -Oz -o "${OUTDIR}/sv_combined/all_sv.merged.vcf.gz"
    bcftools sort \
        -Oz -o "${OUTDIR}/sv_combined/all_sv.sorted.vcf.gz" \
        "${OUTDIR}/sv_combined/all_sv.merged.vcf.gz"
    bcftools index -t "${OUTDIR}/sv_combined/all_sv.sorted.vcf.gz"
    echo "  Combined SVs: $(bcftools view -H "${OUTDIR}/sv_combined/all_sv.sorted.vcf.gz" | wc -l) records"
fi

# ─────────────────────────────────────────
# E. Summary
# ─────────────────────────────────────────
echo ""
echo "=== STEP 03 Summary ==="
echo "  Clair3 normalized VCFs: ${OUTDIR}/normalized/"
echo "  DELLY combined SVs:     ${OUTDIR}/sv_combined/delly_all_sv.sorted.vcf.gz"
echo "  Manta SVs:              ${OUTDIR}/manta_merged/"
echo "  All SVs combined:       ${OUTDIR}/sv_combined/all_sv.sorted.vcf.gz"

for CHR in "${CHROMS[@]}"; do
    f="${OUTDIR}/normalized/${CHR}.clair3.norm.vcf.gz"
    if [[ -f "$f" ]]; then
        echo "    ${CHR}: $(bcftools view -H "$f" | wc -l) variants"
    fi
done
