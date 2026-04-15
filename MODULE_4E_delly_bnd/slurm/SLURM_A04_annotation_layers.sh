#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 0-04:00:00
#SBATCH -J bnd_annot
#SBATCH -o logs/04_annotation.%j.out
#SBATCH -e logs/04_annotation.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4e_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== ANNOTATION LAYERS (BND) ==="

SV_BED_226="${DIR_FINAL}/catalog_226.BND.bed"
dv_check_file "${SV_BED_226}" "226 BND BED"
SV_BED_81=$(ls "${DIR_FINAL}"/catalog_81.BND.*.PASS.bed 2>/dev/null | head -1)
ANNOT_BEDS="${DIR_ANNOT}/beds"
export SV_BED_226

# ── LAYER 1: Gene/Exon/CDS ─────────────────────────────────────────────────
dv_log "--- LAYER 1: Functional overlap ---"

annotate_functional() {
  local sv_bed="$1" label="$2" outpfx="${DIR_ANNOT}/${2}"

  for feat in GENE EXON CDS; do
    local fb="${ANNOT_BEDS}/${feat}.sorted.bed"
    [[ -f "${fb}" ]] && bedtools intersect -a "${sv_bed}" -b "${fb}" -wa -wb \
      > "${outpfx}.${feat}_overlap.tsv" 2>/dev/null || true
  done

  python3 << PYEOF
from collections import defaultdict
sv_bed="${sv_bed}"; outpfx="${outpfx}"
svs={}
with open(sv_bed) as f:
    for line in f:
        p=line.strip().split('\t'); key=f"{p[0]}:{p[1]}-{p[2]}"
        svs[key]={'chrom':p[0],'start':p[1],'end':p[2],'id':p[3] if len(p)>3 else '.','extra':p[4] if len(p)>4 else '.','genes':set(),'exons':set(),'cds':set()}
def load(path,field):
    try:
        with open(path) as f:
            for l in f:
                p=l.strip().split('\t'); k=f"{p[0]}:{p[1]}-{p[2]}"
                if k in svs: svs[k][field].add(p[8] if len(p)>8 else p[7] if len(p)>7 else '.')
    except: pass
load(f"{outpfx}.GENE_overlap.tsv",'genes'); load(f"{outpfx}.EXON_overlap.tsv",'exons'); load(f"{outpfx}.CDS_overlap.tsv",'cds')
with open(f"{outpfx}.functional_class.tsv",'w') as o:
    o.write("sv_key\tchrom\tstart\tend\tid\textra\tn_genes\tn_exons\tn_cds\tclass\tgene_names\n")
    for k,d in svs.items():
        ng,ne,nc=len(d['genes']),len(d['exons']),len(d['cds'])
        c="CDS_overlap" if nc>0 else "exon_overlap" if ne>0 else "intronic" if ng>0 else "intergenic"
        g=','.join(sorted(d['genes'])) if d['genes'] else '.'
        o.write(f"{k}\t{d['chrom']}\t{d['start']}\t{d['end']}\t{d['id']}\t{d['extra']}\t{ng}\t{ne}\t{nc}\t{c}\t{g}\n")
PYEOF

  dv_log "  ${label}:"; tail -n +2 "${outpfx}.functional_class.tsv" | cut -f10 | sort | uniq -c | sort -rn | while read c l; do dv_log "    ${l}: ${c}"; done
}

annotate_functional "${SV_BED_226}" "catalog_226"
[[ -n "${SV_BED_81:-}" && -f "${SV_BED_81}" ]] && annotate_functional "${SV_BED_81}" "catalog_81"

# ── LAYER 2: Repeat overlap ────────────────────────────────────────────────
dv_log "--- LAYER 2: Repeat overlap ---"
if [[ -f "${ANNOT_BEDS}/REPEATS.sorted.bed" ]]; then
  for label_bed in "catalog_226:${SV_BED_226}" "catalog_81:${SV_BED_81:-}"; do
    label="${label_bed%%:*}"; bed="${label_bed#*:}"
    [[ -n "$bed" && -f "$bed" ]] || continue
    bedtools intersect -a "$bed" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -u > "${DIR_ANNOT}/${label}.BNDs_in_repeats.bed" 2>/dev/null || true
    bedtools intersect -a "$bed" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -v > "${DIR_ANNOT}/${label}.BNDs_not_in_repeats.bed" 2>/dev/null || true
    dv_log "  ${label}: $(wc -l < "${DIR_ANNOT}/${label}.BNDs_in_repeats.bed") in repeats"
  done
fi

# ── No SVLEN QC for BND (no meaningful SVLEN) ──────────────────────────
dv_log "--- Skipping SVLEN QC (not applicable for BND) ---"

dv_log "=== ANNOTATION COMPLETE ==="
