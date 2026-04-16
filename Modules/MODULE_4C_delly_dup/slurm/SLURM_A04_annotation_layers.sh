#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 1-00:00:00
#SBATCH -J delly_dup_annot
#SBATCH -o logs/04_annotation.%j.out
#SBATCH -e logs/04_annotation.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/../00_module4c_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }

set -a
source "${CONFIG}"
set +a

# DUP-specific annotation directories
DIR_ANNOT="${OUTDIR}/08_annotation"
DIR_DEPTH="${OUTDIR}/09_depth_support"
DIR_MATDIST="${OUTDIR}/10_mate_distance_qc"

mkdir -p "${DIR_ANNOT}/beds" "${DIR_DEPTH}" "${DIR_MATDIST}"

dv_log "=== DUP ANNOTATION LAYERS ==="

ANNOT_DIR="${DELLY_PROJECT}/pa_roary_results/annotation_flag_beds"
ANNOT_BEDS="${DIR_ANNOT}/beds"
REPEAT_BED="${DELLY_PROJECT}/fClaHyb_Gar_LG.mask_regions.softacgt.renamed.bed"

# Prepare sorted annotation beds
for feat in GENE EXON CDS; do
  src="${ANNOT_DIR}/features/${feat}.bed"
  if [[ -f "${src}" ]]; then
    bedtools sort -i "${src}" > "${ANNOT_BEDS}/${feat}.sorted.bed"
    dv_log "Prepared ${feat} BED"
  else
    dv_log "WARNING: missing ${src}"
  fi
done

if [[ -f "${REPEAT_BED}" ]]; then
  bedtools sort -i "${REPEAT_BED}" > "${ANNOT_BEDS}/REPEATS.sorted.bed"
fi

DUP_BED_226="${DIR_FINAL}/catalog_226.DUP.bed"
dv_check_file "${DUP_BED_226}" "226 DUP BED"

DUP_BED_81=$(ls "${DIR_FINAL}"/catalog_81.DUP.*.PASS.bed 2>/dev/null | head -1 || true)

# ── LAYER 1: Functional overlap ────────────────────────────────────────────
annotate_functional() {
  local in_bed="$1"
  local label="$2"
  local outpfx="${DIR_ANNOT}/${label}"

  for feat in GENE EXON CDS; do
    local fb="${ANNOT_BEDS}/${feat}.sorted.bed"
    if [[ -f "${fb}" ]]; then
      bedtools intersect -a "${in_bed}" -b "${fb}" -wa -wb > "${outpfx}.${feat}_overlap.tsv" 2>/dev/null || true
    fi
  done

  python3 << PYEOF
from collections import defaultdict

dup_bed = "${in_bed}"
outpfx = "${outpfx}"

events = {}
with open(dup_bed) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        chrom, start, end = p[0], p[1], p[2]
        dup_id = p[3] if len(p) > 3 else f"{chrom}:{start}-{end}"
        span = int(end) - int(start)
        key = dup_id
        events[key] = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "dup_id": dup_id,
            "span": span,
            "genes": set(),
            "exons": set(),
            "cds": set(),
        }

def load_hits(path, field):
    try:
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                dup_id = p[3] if len(p) > 3 else f"{p[0]}:{p[1]}-{p[2]}"
                gene_name = p[8] if len(p) > 8 else (p[7] if len(p) > 7 else ".")
                if dup_id in events:
                    events[dup_id][field].add(gene_name)
    except FileNotFoundError:
        pass

load_hits(f"{outpfx}.GENE_overlap.tsv", "genes")
load_hits(f"{outpfx}.EXON_overlap.tsv", "exons")
load_hits(f"{outpfx}.CDS_overlap.tsv", "cds")

with open(f"{outpfx}.functional_class.tsv", "w") as o:
    o.write("dup_key\tchrom\tstart\tend\tid\tspan_bp\tn_genes\tn_exons\tn_cds\tclass\tgene_names\n")
    for dup_id, d in events.items():
        ng, ne, nc = len(d["genes"]), len(d["exons"]), len(d["cds"])
        if nc > 0:
            cls = "CDS_overlap"
        elif ne > 0:
            cls = "exon_overlap"
        elif ng > 0:
            cls = "intronic"
        else:
            cls = "intergenic"
        genes = ",".join(sorted(d["genes"])) if d["genes"] else "."
        o.write(f"{dup_id}\t{d['chrom']}\t{d['start']}\t{d['end']}\t{d['dup_id']}\t{d['span']}\t{ng}\t{ne}\t{nc}\t{cls}\t{genes}\n")
PYEOF
}

dv_log "--- LAYER 1: Functional overlap ---"
annotate_functional "${DUP_BED_226}" "catalog_226"
if [[ -n "${DUP_BED_81}" && -f "${DUP_BED_81}" ]]; then
  annotate_functional "${DUP_BED_81}" "catalog_81"
fi

# ── LAYER 2: Repeat overlap ────────────────────────────────────────────────
dv_log "--- LAYER 2: Repeat overlap ---"
if [[ -f "${ANNOT_BEDS}/REPEATS.sorted.bed" ]]; then
  for pair in "catalog_226:${DUP_BED_226}" "catalog_81:${DUP_BED_81:-}"; do
    label="${pair%%:*}"
    bed="${pair#*:}"
    [[ -n "${bed}" && -f "${bed}" ]] || continue

    bedtools intersect -a "${bed}" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -u > "${DIR_ANNOT}/${label}.DUPs_in_repeats.bed" 2>/dev/null || true
    bedtools intersect -a "${bed}" -b "${ANNOT_BEDS}/REPEATS.sorted.bed" -wa -f 0.5 -v > "${DIR_ANNOT}/${label}.DUPs_not_in_repeats.bed" 2>/dev/null || true
  done
fi

# ── LAYER 3: Depth support ─────────────────────────────────────────────────
dv_log "--- LAYER 3: Depth support ---"
PA_MOSDEPTH_DIR="${DELLY_PROJECT}/pa_roary_results/01_mosdepth"
MOSDEPTH_SOURCE=""
for d in "${PA_MOSDEPTH_DIR}" "${DELLY_PROJECT}/pa_roary_results/01_mosdepth_callable"; do
  if ls "${d}"/*.regions.bed.gz &>/dev/null 2>&1 || ls "${d}"/*/*.regions.bed.gz &>/dev/null 2>&1; then
    MOSDEPTH_SOURCE="${d}"
    break
  fi
done

if [[ -n "${MOSDEPTH_SOURCE}" ]]; then
  export DUP_BED_226 MOSDEPTH_SOURCE DIR_DEPTH
  python3 << 'PYEOF'
import gzip, os, glob, statistics
from collections import defaultdict

dup_bed = os.environ["DUP_BED_226"]
md = os.environ["MOSDEPTH_SOURCE"]
out = os.environ["DIR_DEPTH"] + "/depth_support_226.tsv"
flank = 5000

sites = []
with open(dup_bed) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        chrom, start, end = p[0], int(p[1]), int(p[2])
        dup_id = p[3] if len(p) > 3 else f"{chrom}:{start}-{end}"
        sites.append((chrom, start, end, dup_id))

rfiles = glob.glob(os.path.join(md, "*.regions.bed.gz")) + glob.glob(os.path.join(md, "*", "*.regions.bed.gz"))
rfiles = rfiles[:40]

def load_regions(path):
    d = defaultdict(list)
    with gzip.open(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            d[p[0]].append((int(p[1]), int(p[2]), float(p[3])))
    return d

def mean_depth_in_interval(regions, chrom, start, end):
    if chrom not in regions:
        return None
    total_bp = 0
    total_depth = 0.0
    for s, e, dep in regions[chrom]:
        a = max(s, start)
        b = min(e, end)
        if a < b:
            bp = b - a
            total_bp += bp
            total_depth += dep * bp
    if total_bp == 0:
        return None
    return total_depth / total_bp

dup_ratios = defaultdict(list)

for rf in rfiles:
    regs = load_regions(rf)
    chrom_med = {}
    for chrom, arr in regs.items():
        vals = [x[2] for x in arr if x[2] > 0]
        if vals:
            chrom_med[chrom] = statistics.median(vals)

    for chrom, start, end, dup_id in sites:
        if chrom not in chrom_med or chrom_med[chrom] == 0:
            continue

        inside = mean_depth_in_interval(regs, chrom, start, end)
        if inside is None:
            continue

        left = mean_depth_in_interval(regs, chrom, max(0, start - flank), start)
        right = mean_depth_in_interval(regs, chrom, end, end + flank)
        flank_vals = [x for x in (left, right) if x is not None]
        if not flank_vals:
            continue

        flank_mean = statistics.mean(flank_vals)
        if flank_mean == 0:
            continue

        ratio = inside / flank_mean
        dup_ratios[dup_id].append(ratio)

with open(out, "w") as o:
    o.write("dup_id\tn_samples\tmedian_ratio\tmean_ratio\tdepth_label\n")
    for chrom, start, end, dup_id in sites:
        vals = dup_ratios.get(dup_id, [])
        if not vals:
            o.write(f"{dup_id}\t0\tNA\tNA\tdepth_no_data\n")
            continue
        med = statistics.median(vals)
        mean = statistics.mean(vals)
        if med >= 1.5:
            label = "strong_DUP_depth_support"
        elif med >= 1.2:
            label = "moderate_DUP_depth_support"
        elif med >= 1.05:
            label = "weak_DUP_depth_support"
        else:
            label = "no_depth_support"
        o.write(f"{dup_id}\t{len(vals)}\t{med:.4f}\t{mean:.4f}\t{label}\n")
PYEOF
fi

# ── LAYER 4: Span QC ───────────────────────────────────────────────────────
dv_log "--- LAYER 4: Span QC ---"
FINAL_226_VCF="${DIR_FINAL}/catalog_226.DUP.vcf.gz"
export FINAL_226_VCF DIR_MATDIST
python3 << 'PYEOF'
import gzip, os

vcf = os.environ["FINAL_226_VCF"]
out = os.environ["DIR_MATDIST"] + "/span_qc_226.tsv"

with gzip.open(vcf, "rt") as f, open(out, "w") as o:
    o.write("chrom\tpos\tend\tid\tspan_bp\tspan_kb\tflag\n")
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        chrom, pos, dup_id = p[0], int(p[1]), p[2]
        info = {}
        for x in p[7].split(";"):
            if "=" in x:
                k, v = x.split("=", 1)
                info[k] = v
        end = int(info.get("END", pos))
        span = abs(end - pos)
        kb = span / 1000.0
        if kb >= 500:
            flag = "extreme_large"
        elif kb >= 100:
            flag = "very_large"
        elif kb >= 20:
            flag = "large"
        else:
            flag = "normal"
        o.write(f"{chrom}\t{pos}\t{end}\t{dup_id}\t{span}\t{kb:.2f}\t{flag}\n")
PYEOF

dv_log "=== DUP ANNOTATION COMPLETE ==="
