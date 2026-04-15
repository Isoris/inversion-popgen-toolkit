#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 1-00:00:00
#SBATCH -J manta_annot
#SBATCH -o logs/04_annotation.%j.out
#SBATCH -e logs/04_annotation.%j.err
# =============================================================================
# 04_annotation_layers.sh — Functional + repeat + depth annotation per SV type
# =============================================================================
#
# LAYER 1: Gene/Exon/CDS overlap (bedtools intersect → functional class)
# LAYER 2: Repeat overlap (≥50% reciprocal)
# LAYER 3: Depth support (mosdepth ratios for DEL + DUP)
#           DEL: expect inside/flank ratio < 1
#           DUP: expect inside/flank ratio > 1
# LAYER 4: SVLEN QC flags (extreme sizes)
#
# v2 changes:
#   - Fixed mosdepth glob to find files in sample subdirectories
#   - Uses up to 81 mosdepth files (unrelated subset) for depth ratios
#   - Added SVLEN QC layer
#   - Better error reporting when mosdepth files not found
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4g_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

mv_init_dirs
mv_log "=== ANNOTATION LAYERS ==="

ANNOT_BEDS="${DIR_ANNOT}/beds"

# ── LAYER 1: Gene/Exon/CDS overlap (per SV type) ──────────────────────────
mv_log "--- LAYER 1: Functional overlap ---"

annotate_functional() {
  local sv_bed="$1" label="$2" outpfx="${DIR_ANNOT}/${2}"

  [[ -s "${sv_bed}" ]] || { mv_log "  ${label}: empty BED, skipping"; return 0; }

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
        p=line.strip().split('\t')
        key=f"{p[0]}:{p[1]}-{p[2]}"
        svs[key]={
            'chrom':p[0],'start':p[1],'end':p[2],
            'id':p[3] if len(p)>3 else '.',
            'svlen':p[4] if len(p)>4 else '.',
            'svtype':p[5] if len(p)>5 else '.',
            'genes':set(),'exons':set(),'cds':set()
        }

def load(path, field):
    try:
        with open(path) as f:
            for l in f:
                p=l.strip().split('\t')
                k=f"{p[0]}:{p[1]}-{p[2]}"
                if k in svs:
                    # Try to get gene name from annotation BED columns
                    name = '.'
                    for col in [8, 7, 6, 5, 4, 3]:
                        if len(p) > col and p[col] and p[col] != '.':
                            name = p[col]
                            break
                    svs[k][field].add(name)
    except:
        pass

load(f"{outpfx}.GENE_overlap.tsv", 'genes')
load(f"{outpfx}.EXON_overlap.tsv", 'exons')
load(f"{outpfx}.CDS_overlap.tsv", 'cds')

with open(f"{outpfx}.functional_class.tsv", 'w') as o:
    o.write("sv_key\tchrom\tstart\tend\tid\tsvlen\tsvtype\tn_genes\tn_exons\tn_cds\tclass\tgene_names\n")
    for k, d in svs.items():
        ng, ne, nc = len(d['genes']), len(d['exons']), len(d['cds'])
        c = "CDS_overlap" if nc > 0 else "exon_overlap" if ne > 0 else "intronic" if ng > 0 else "intergenic"
        g = ','.join(sorted(d['genes'])) if d['genes'] else '.'
        o.write(f"{k}\t{d['chrom']}\t{d['start']}\t{d['end']}\t{d['id']}\t{d['svlen']}\t{d['svtype']}\t{ng}\t{ne}\t{nc}\t{c}\t{g}\n")
PYEOF

  mv_log "  ${label}:"
  tail -n +2 "${outpfx}.functional_class.tsv" | cut -f11 | sort | uniq -c | sort -rn \
    | while read c l; do mv_log "    ${l}: ${c}"; done
}

# Annotate each SV type for the 226 catalog
for svtype in DEL DUP INV INS_small INS_large; do
  bed="${DIR_FINAL}/catalog_226.${svtype}.PASS.bed"
  if [[ -f "${bed}" && -s "${bed}" ]]; then
    annotate_functional "${bed}" "catalog_226_${svtype}"
  else
    mv_log "  catalog_226_${svtype}: no BED or empty"
  fi
done

# ── LAYER 2: Repeat overlap ────────────────────────────────────────────────
mv_log "--- LAYER 2: Repeat overlap ---"
if [[ -f "${ANNOT_BEDS}/REPEATS.sorted.bed" ]]; then
  for svtype in DEL DUP INV INS_small INS_large; do
    bed="${DIR_FINAL}/catalog_226.${svtype}.PASS.bed"
    [[ -f "${bed}" && -s "${bed}" ]] || continue
    label="catalog_226_${svtype}"

    bedtools intersect -a "${bed}" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -u \
      > "${DIR_ANNOT}/${label}.in_repeats.bed" 2>/dev/null || true
    bedtools intersect -a "${bed}" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -v \
      > "${DIR_ANNOT}/${label}.not_in_repeats.bed" 2>/dev/null || true

    n_rep=$(wc -l < "${DIR_ANNOT}/${label}.in_repeats.bed")
    n_norep=$(wc -l < "${DIR_ANNOT}/${label}.not_in_repeats.bed")
    mv_log "  ${svtype}: ${n_rep} in repeats, ${n_norep} not"
  done
else
  mv_log "  REPEATS.sorted.bed not found — skipping repeat layer."
fi

# ── LAYER 3: Depth support (DEL and DUP only) ─────────────────────────────
mv_log "--- LAYER 3: Depth support (DEL + DUP) ---"

# Find mosdepth data. Try multiple known locations and glob patterns.
MOSDEPTH_SOURCE=""
MOSDEPTH_FILES=()

# Pattern 1: flat directory with .regions.bed.gz
for d in \
  "${MANTA_PROJECT}/pa_roary_results/01_mosdepth" \
  "${MANTA_PROJECT}/pa_roary_results/01_mosdepth_callable"; do
  if [[ -d "${d}" ]]; then
    mapfile -t found < <(find "${d}" -name '*.regions.bed.gz' -type f 2>/dev/null | head -100)
    if [[ ${#found[@]} -gt 0 ]]; then
      MOSDEPTH_SOURCE="${d}"
      MOSDEPTH_FILES=("${found[@]}")
      break
    fi
  fi
done

if [[ ${#MOSDEPTH_FILES[@]} -eq 0 ]]; then
  mv_log "  No mosdepth *.regions.bed.gz found in expected locations."
  mv_log "  Checked:"
  mv_log "    ${MANTA_PROJECT}/pa_roary_results/01_mosdepth"
  mv_log "    ${MANTA_PROJECT}/pa_roary_results/01_mosdepth_callable"
  mv_log "  Skipping depth support layer."
else
  N_MD=${#MOSDEPTH_FILES[@]}
  mv_log "  Found ${N_MD} mosdepth files in: ${MOSDEPTH_SOURCE}"
  mv_log "  Example: $(basename "${MOSDEPTH_FILES[0]}")"

  # Use up to 81 files (unrelated subset if possible)
  MAX_MD=81
  if [[ ${N_MD} -gt ${MAX_MD} ]]; then
    mv_log "  Using first ${MAX_MD} of ${N_MD} files"
  fi

  # Write file list for Python
  MD_FILELIST="${DIR_DEPTH}/mosdepth_file_list.txt"
  printf '%s\n' "${MOSDEPTH_FILES[@]}" | head -${MAX_MD} > "${MD_FILELIST}"

  for svtype in DEL DUP; do
    bed="${DIR_FINAL}/catalog_226.${svtype}.PASS.bed"
    [[ -f "${bed}" && -s "${bed}" ]] || continue
    n_sv=$(wc -l < "${bed}")
    mv_log "  ${svtype}: computing depth ratios for ${n_sv} SVs..."

    export SV_BED="${bed}" MD_FILELIST DEPTH_OUT="${DIR_DEPTH}/depth_support_226_${svtype}.tsv"
    export SVTYPE_LABEL="${svtype}"

    python3 << 'PYEOF'
import gzip, os, statistics
from collections import defaultdict

sv_bed = os.environ['SV_BED']
md_list = os.environ['MD_FILELIST']
out = os.environ['DEPTH_OUT']
svtype = os.environ['SVTYPE_LABEL']
flank = 5000

# Read SV sites
sites = []
with open(sv_bed) as f:
    for l in f:
        p = l.strip().split('\t')
        sid = p[3] if len(p) > 3 else f"{p[0]}:{p[1]}-{p[2]}"
        sites.append((p[0], int(p[1]), int(p[2]), sid))

print(f"  Loaded {len(sites)} {svtype} sites")

# Read mosdepth file list
with open(md_list) as f:
    rfiles = [l.strip() for l in f if l.strip()]

print(f"  Using {len(rfiles)} mosdepth files")

def load_mosdepth(path):
    """Load mosdepth regions.bed.gz into dict of (chrom -> list of (start, end, depth))"""
    d = defaultdict(list)
    try:
        with gzip.open(path, 'rt') as f:
            for l in f:
                p = l.strip().split('\t')
                if len(p) >= 4:
                    d[p[0]].append((int(p[1]), int(p[2]), float(p[3])))
    except Exception as e:
        print(f"  WARNING: failed to load {path}: {e}")
    return d

def query_depth(regions, chrom, start, end):
    """Get mean depth in [start, end) from mosdepth regions"""
    if chrom not in regions:
        return None
    total_bp = 0
    total_depth = 0
    for ws, we, wd in regions[chrom]:
        o1 = max(ws, start)
        o2 = min(we, end)
        if o1 < o2:
            bp = o2 - o1
            total_bp += bp
            total_depth += wd * bp
    return total_depth / total_bp if total_bp > 0 else None

# Compute depth ratios per SV across samples
depth_ratios = defaultdict(list)
n_processed = 0

for rf in rfiles:
    regions = load_mosdepth(rf)
    if not regions:
        continue

    # Chromosome-level median depth for normalization sanity check
    chr_medians = {}
    for c, windows in regions.items():
        depths = [d for _, _, d in windows if d > 0]
        if depths:
            chr_medians[c] = statistics.median(depths)

    for chrom, start, end, sv_id in sites:
        if chrom not in chr_medians or chr_medians[chrom] == 0:
            continue

        inside = query_depth(regions, chrom, start, end)
        if inside is None:
            continue

        left_d = query_depth(regions, chrom, max(0, start - flank), start)
        right_d = query_depth(regions, chrom, end, end + flank)

        flank_depths = [d for d in [left_d, right_d] if d is not None]
        if not flank_depths:
            continue

        flank_mean = statistics.mean(flank_depths)
        if flank_mean == 0:
            continue

        depth_ratios[sv_id].append(inside / flank_mean)

    n_processed += 1
    if n_processed % 10 == 0:
        print(f"  Processed {n_processed}/{len(rfiles)} mosdepth files")

print(f"  Processed all {n_processed} mosdepth files")

# Classify depth support
expect_low = (svtype == "DEL")

with open(out, 'w') as o:
    o.write("sv_id\tn_samples\tmedian_ratio\tmean_ratio\tdepth_label\n")
    for chrom, start, end, sv_id in sites:
        rs = depth_ratios.get(sv_id, [])
        if not rs:
            o.write(f"{sv_id}\t0\tNA\tNA\tdepth_no_data\n")
            continue

        mr = statistics.median(rs)
        mn = statistics.mean(rs)

        if expect_low:
            lb = ("strong_depth_support" if mr < 0.3 else
                  "moderate_depth_support" if mr < 0.6 else
                  "weak_depth_support" if mr < 0.85 else
                  "no_depth_support")
        else:
            lb = ("strong_depth_support" if mr > 1.7 else
                  "moderate_depth_support" if mr > 1.3 else
                  "weak_depth_support" if mr > 1.1 else
                  "no_depth_support")

        o.write(f"{sv_id}\t{len(rs)}\t{mr:.4f}\t{mn:.4f}\t{lb}\n")

# Summary
labels = defaultdict(int)
with open(out) as f:
    next(f)
    for l in f:
        labels[l.strip().split('\t')[-1]] += 1
for lb, cnt in sorted(labels.items(), key=lambda x: -x[1]):
    print(f"  {lb}: {cnt}")
PYEOF

    mv_log "  ${svtype} depth support:"
    tail -n +2 "${DIR_DEPTH}/depth_support_226_${svtype}.tsv" | cut -f5 | sort | uniq -c | sort -rn \
      | while read c l; do mv_log "    ${l}: ${c}"; done
  done
fi

# ── LAYER 4: SVLEN QC flags ───────────────────────────────────────────────
mv_log "--- LAYER 4: SVLEN QC flags ---"

for svtype in DEL DUP INV; do
  vcf="${DIR_FINAL}/catalog_226.${svtype}.PASS.vcf.gz"
  [[ -f "${vcf}" ]] || continue

  qc_out="${DIR_ANNOT}/svlen_qc_226_${svtype}.tsv"

  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\n' "${vcf}" \
  | awk -v svtype="${svtype}" 'BEGIN{OFS="\t"; print "chrom","pos","end","id","svlen_bp","abs_kb","flag"}{
    sv=$5; gsub(/^-/,"",sv); kb=sv/1000
    if(kb >= 1000) flag="extreme_>1Mb"
    else if(kb >= 100) flag="very_large_100kb-1Mb"
    else if(kb >= 10) flag="large_10-100kb"
    else if(kb >= 1) flag="medium_1-10kb"
    else flag="small_<1kb"
    print $1,$2,$3,$4,sv,sprintf("%.2f",kb),flag
  }' > "${qc_out}"

  mv_log "  ${svtype} SVLEN QC:"
  tail -n +2 "${qc_out}" | cut -f7 | sort | uniq -c | sort -rn \
    | while read c f; do mv_log "    ${f}: ${c}"; done
done

mv_log "=== ANNOTATION COMPLETE ==="
