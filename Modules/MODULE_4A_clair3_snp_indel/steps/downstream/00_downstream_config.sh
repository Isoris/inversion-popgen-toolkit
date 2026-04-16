#!/usr/bin/env bash
# ============================================================
# 00_config.sh — Central config for Clair3 downstream analysis
#
# Source this from every script:
#   source "$(dirname "$0")/00_config.sh"
# ============================================================

# ── Project roots ──
export BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
export C3_ROOT="${BASE}/MODULE_4A_SNP_INDEL50_Clair3"
export PP_RESULTS="${C3_ROOT}/postprocess_results"
export PP_SCRIPTS="${C3_ROOT}/postprocess_scripts"

# ── Reference genome ──
export REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
export REF_FAI="${REF}.fai"

# ── Sample manifests ──
export BAM_MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
export SAMPLES_ALL="${BASE}/pa_roary_results/00_manifests/sample_list_226.txt"
export SAMPLES_UNRELATED="${BASE}/pa_roary_results/00_manifests/samples_unrelated_81.txt"

# ── Ancestry (NGSadmix) ──
export PA_NGSADMIX_DIR="${BASE}/pa_roary_results/03_ngsadmix"
export PA_BEST_SEED_TABLE="${PA_NGSADMIX_DIR}/best_seed_per_K.tsv"

# ── Annotation BEDs (from DELLY pipeline or shared) ──
export ANNOT_BEDS="${BASE}/MODULE_5_DELLY_DEL/10_annotation/beds"

# ── Downstream output root ──
export DS_ROOT="${C3_ROOT}/downstream_results"

# ── Chromosome list ──
export CHROM_LIST="${C3_ROOT}/meta/chromosome_list.txt"

# ── Helper functions ──
ds_log() { echo "[$(date '+%H:%M:%S')] $*"; }

ds_init_dirs() {
    local chrom="${1:-}"
    if [[ -n "$chrom" ]]; then
        export DS_CHROM_DIR="${DS_ROOT}/${chrom}"
        mkdir -p "${DS_CHROM_DIR}"/{catalogs,matrices,annotation,per_sample,distances,markers,gene_tables,figures,phase_prep}
    fi
    mkdir -p "${DS_ROOT}"/{_all_chroms,_figures,_tables,logs}
}

ds_resolve_ancestry() {
    # Find best K ancestry and set QOPT + QOPT_SAMPLES
    export QOPT="" QOPT_SAMPLES=""
    if [[ -f "${PA_BEST_SEED_TABLE}" ]]; then
        local best_k best_prefix
        best_k=$(awk -F'\t' 'NR>1 { if ($5 > max || NR==2) { max=$5; k=$1 } } END { print k }' \
                 "${PA_BEST_SEED_TABLE}")
        best_prefix=$(awk -F'\t' -v k="${best_k}" 'NR>1 && $1==k { print $4 }' "${PA_BEST_SEED_TABLE}")
        local qopt="${PA_NGSADMIX_DIR}/${best_prefix}.qopt"
        if [[ -f "${qopt}" ]]; then
            export QOPT="${qopt}"
            export QOPT_SAMPLES="${SAMPLES_ALL}"
            ds_log "Ancestry: K=${best_k} from ${qopt}"
        fi
    fi
}
