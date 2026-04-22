#!/usr/bin/env bash
###############################################################################
# SLURM_B06_qc_depth_table.slurm
#
# Comprehensive depth and coverage QC using mosdepth, across multiple region
# classes (whole-chromosome, repeat, nonrepeat, callable, noncallable).
# Computes per-sample overall + per-chromosome tables, plus summaries.
#
# Region classes come from the mask BED files generated during genome assembly.
#
# Usage:
#   sbatch SLURM_B06_qc_depth_table.slurm
#
# Outputs:
#   ${OUTDIR}/all_samples.coverage_context_qc.tsv
#   ${OUTDIR}/all_samples.per_chrom.coverage_context_qc.tsv
#   ${OUTDIR}/summary.coverage_context_qc.tsv
#   ${OUTDIR}/summary.per_chrom.coverage_context_qc.tsv
#   ${OUTDIR}/SLURM_B06_qc_depth_table.arg
#   ${OUTDIR}/SLURM_B06_qc_depth_table.results
###############################################################################
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J qc_covmap
#SBATCH -o qc_covmap.%j.out
#SBATCH -e qc_covmap.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

STEP="SLURM_B06_qc_depth_table"
timestamp(){ date '+%F %T'; }

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
FAI="${REF}.fai"

PA_REPEAT_BED="${BASE}/fClaHyb_Gar_LG.mask_regions.softacgt.bed"
PA_CALLABLE_BED="${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.0based.bed"

OUTDIR="${BASE}/02-merged_per_sample/results/popgen_qc"
LOGDIR="${OUTDIR}/logs"
TMPROOT="${OUTDIR}/tmp_covmap_ctx"
WORKDIR="${TMPROOT}/work"
PER_SAMPLE_DIR="${TMPROOT}/per_sample"
PER_CHR_DIR="${TMPROOT}/per_chr"

mkdir -p "$OUTDIR" "$LOGDIR" "$TMPROOT" "$WORKDIR" "$PER_SAMPLE_DIR" "$PER_CHR_DIR"

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-64}"
JOBS="${JOBS:-8}"
THREADS_PER_JOB="${THREADS_PER_JOB:-4}"

if (( JOBS * THREADS_PER_JOB > TOTAL_CPUS )); then
  THREADS_PER_JOB=4
  JOBS=$(( TOTAL_CPUS / THREADS_PER_JOB ))
  (( JOBS < 1 )) && JOBS=1
fi

# ---- Tool checks ----
for tool in samtools mosdepth bedtools bgzip tabix; do
  command -v "$tool" >/dev/null 2>&1 || { echo "[ERROR] $tool not found"; exit 1; }
done

# ---- Write .arg ----
ARGFILE="${OUTDIR}/${STEP}.arg"
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "slurm_job_id\t${SLURM_JOB_ID:-NA}"
  echo -e "base\t${BASE}"
  echo -e "manifest\t${MANIFEST}"
  echo -e "reference\t${REF}"
  echo -e "repeat_bed\t${PA_REPEAT_BED}"
  echo -e "callable_bed\t${PA_CALLABLE_BED}"
  echo -e "total_cpus\t${TOTAL_CPUS}"
  echo -e "jobs\t${JOBS}"
  echo -e "threads_per_job\t${THREADS_PER_JOB}"
  echo -e "mosdepth_thresholds\t1,5,10"
  echo -e "region_classes\tchr_all,repeat,nonrepeat,callable,noncallable"
} > "$ARGFILE"

[[ -s "$FAI" ]] || samtools faidx "$REF"

# ---- 1) Build region BEDs ----
CHR_BED="${WORKDIR}/ref.chromosomes.bed"
awk 'BEGIN{OFS="\t"} {print $1,0,$2,$1}' "$FAI" > "$CHR_BED"

REPEAT_SORTED="${WORKDIR}/repeat.sorted.bed"
CALLABLE_SORTED="${WORKDIR}/callable.sorted.bed"
sort -k1,1 -k2,2n "$PA_REPEAT_BED" > "$REPEAT_SORTED"
sort -k1,1 -k2,2n "$PA_CALLABLE_BED" > "$CALLABLE_SORTED"

NONREPEAT_BED="${WORKDIR}/nonrepeat.bed"
NONCALLABLE_BED="${WORKDIR}/noncallable.bed"
bedtools complement -i "$REPEAT_SORTED" -g "$FAI" > "$NONREPEAT_BED"
bedtools complement -i "$CALLABLE_SORTED" -g "$FAI" > "$NONCALLABLE_BED"

# ---- 2) Combined mosdepth regions BED ----
COMBINED_BED="${WORKDIR}/combined_regions.bed"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"chr_all|" $1 "|" NR}' "$CHR_BED" > "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"repeat|" $1 "|" NR}' "$REPEAT_SORTED" >> "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"nonrepeat|" $1 "|" NR}' "$NONREPEAT_BED" >> "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"callable|" $1 "|" NR}' "$CALLABLE_SORTED" >> "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"noncallable|" $1 "|" NR}' "$NONCALLABLE_BED" >> "$COMBINED_BED"
sort -k1,1 -k2,2n "$COMBINED_BED" -o "$COMBINED_BED"

# ---- 3) Sample list ----
SAMPLE_BAMS="${WORKDIR}/sample_bams.tsv"
awk -F'\t' 'NR>1 && $1!="" && $3!="" {print $1 "\t" $3}' "$MANIFEST" > "$SAMPLE_BAMS"
NSAMPLES=$(wc -l < "$SAMPLE_BAMS")
echo "[INFO] Samples: ${NSAMPLES}"
(( NSAMPLES == 0 )) && { echo "[ERROR] No samples in manifest"; exit 1; }

# ---- 4) Worker script ----
WORKER="${WORKDIR}/worker_qc_covmap.sh"
cat > "$WORKER" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail
sample="$1"; bam="$2"; base="$3"; combined_bed="$4"
per_sample_dir="$5"; per_chr_dir="$6"; threads="$7"; workdir="$8"
prefix="${workdir}/${sample}"

[[ -s "$bam" ]] || { echo "[WARN] Missing BAM: ${sample}" >&2; exit 0; }

mosdepth -n -t "$threads" --by "$combined_bed" --thresholds 1,5,10 "$prefix" "$bam" >/dev/null 2>&1

samtools view -@ "$threads" -F 0x904 "$bam" \
| awk 'BEGIN{OFS="\t"} {chr=$3; total++; tot[chr]++; if($5==0){mq0++; mq0c[chr]++} if($5>=30){mq30++; mq30c[chr]++}} END{print "TOTAL",total+0,mq0+0,mq30+0; for(c in tot) print c,tot[c]+0,mq0c[c]+0,mq30c[c]+0}' > "${prefix}.mapq_counts.tsv"

python3 - "$sample" "$bam" "${prefix}.regions.bed.gz" "${prefix}.thresholds.bed.gz" "${prefix}.mapq_counts.tsv" "${per_sample_dir}/${sample}.overall.tsv" "${per_chr_dir}/${sample}.per_chrom.tsv" <<'PY'
import sys, gzip, csv
from collections import defaultdict

sample, bam, regions_gz, thresholds_gz, mapq_tsv, overall_out, perchr_out = sys.argv[1:]

acc = defaultdict(lambda: {"length":0, "depth_sum":0.0, "cov1":0, "cov5":0, "cov10":0})

with gzip.open(regions_gz, "rt") as f:
    for line in f:
        if not line.strip(): continue
        chrom, start, end, name, mean_cov = line.rstrip("\n").split("\t")[:5]
        class_name, chrom_from_name, _ = name.split("|", 2)
        length = int(end) - int(start)
        acc[(class_name, chrom_from_name)]["length"] += length
        acc[(class_name, chrom_from_name)]["depth_sum"] += float(mean_cov) * length

with gzip.open(thresholds_gz, "rt") as f:
    for line in f:
        if not line.strip(): continue
        parts = line.rstrip("\n").split("\t")
        chrom, start, end, name = parts[:4]
        vals = parts[4:]
        class_name, chrom_from_name, _ = name.split("|", 2)
        acc[(class_name, chrom_from_name)]["cov1"] += int(vals[0]) if len(vals) > 0 else 0
        acc[(class_name, chrom_from_name)]["cov5"] += int(vals[1]) if len(vals) > 1 else 0
        acc[(class_name, chrom_from_name)]["cov10"] += int(vals[2]) if len(vals) > 2 else 0

mapq = {}
with open(mapq_tsv) as f:
    for line in f:
        ch, total, mq0, mq30 = line.rstrip("\n").split("\t")
        total, mq0, mq30 = int(total), int(mq0), int(mq30)
        mapq[ch] = {"mapped_primary_reads": total, "mapq0_reads": mq0, "mapq30plus_reads": mq30,
                     "mapq0_frac": mq0/total if total>0 else 0, "mapq30plus_frac": mq30/total if total>0 else 0}

classes = ["chr_all", "repeat", "nonrepeat", "callable", "noncallable"]

def metrics_for(cn, ch):
    d = acc.get((cn, ch))
    if not d or d["length"]==0:
        return {f"{cn}_region_bp":0, f"{cn}_mean_depth_x":"", f"{cn}_breadth_1x":"", f"{cn}_breadth_5x":"", f"{cn}_breadth_10x":""}
    L = d["length"]
    return {f"{cn}_region_bp":L, f"{cn}_mean_depth_x":d["depth_sum"]/L, f"{cn}_breadth_1x":d["cov1"]/L, f"{cn}_breadth_5x":d["cov5"]/L, f"{cn}_breadth_10x":d["cov10"]/L}

chroms = sorted({ch for (_, ch) in acc.keys() if ch != "TOTAL"})

fieldnames = ["sample","bam","chrom","mapped_primary_reads","mapq0_reads","mapq30plus_reads","mapq0_frac","mapq30plus_frac"]
for c in classes:
    fieldnames += [f"{c}_region_bp", f"{c}_mean_depth_x", f"{c}_breadth_1x", f"{c}_breadth_5x", f"{c}_breadth_10x"]

with open(perchr_out, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    for ch in chroms:
        row = {"sample":sample, "bam":bam, "chrom":ch}
        row.update({k: mapq.get(ch,{}).get(k,0) for k in ["mapped_primary_reads","mapq0_reads","mapq30plus_reads","mapq0_frac","mapq30plus_frac"]})
        for c in classes: row.update(metrics_for(c, ch))
        w.writerow(row)

overall_acc = defaultdict(lambda: {"length":0,"depth_sum":0.0,"cov1":0,"cov5":0,"cov10":0})
for (cn, ch), d in acc.items():
    for k in d: overall_acc[cn][k] += d[k] if isinstance(d[k], (int, float)) else 0

overall_row = {"sample":sample, "bam":bam}
overall_row.update({k: mapq.get("TOTAL",{}).get(k,0) for k in ["mapped_primary_reads","mapq0_reads","mapq30plus_reads","mapq0_frac","mapq30plus_frac"]})
for c in classes:
    d = overall_acc[c]
    if d["length"]==0:
        for s in ["region_bp","mean_depth_x","breadth_1x","breadth_5x","breadth_10x"]: overall_row[f"{c}_{s}"] = "" if "region" not in s else 0
    else:
        L = d["length"]
        overall_row[f"{c}_region_bp"] = L
        overall_row[f"{c}_mean_depth_x"] = d["depth_sum"]/L
        for t,k in [(1,"cov1"),(5,"cov5"),(10,"cov10")]: overall_row[f"{c}_breadth_{t}x"] = d[k]/L

with open(overall_out, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=list(overall_row.keys()), delimiter="\t")
    w.writeheader()
    w.writerow(overall_row)
PY
BASH
chmod +x "$WORKER"

# ---- 5) Run workers in parallel ----
running=0
while IFS=$'\t' read -r sample bam; do
  echo "[INFO] Launching ${sample}"
  bash "$WORKER" "$sample" "$bam" "$BASE" "$COMBINED_BED" "$PER_SAMPLE_DIR" "$PER_CHR_DIR" "$THREADS_PER_JOB" "$WORKDIR" &
  running=$((running+1))
  if (( running >= JOBS )); then wait -n; running=$((running-1)); fi
done < "$SAMPLE_BAMS"
wait

# ---- 6) Merge per-sample files ----
ALL_OVERALL="${OUTDIR}/all_samples.coverage_context_qc.tsv"
ALL_PERCHR="${OUTDIR}/all_samples.per_chrom.coverage_context_qc.tsv"

first=1
for f in "$PER_SAMPLE_DIR"/*.overall.tsv; do
  [[ -e "$f" ]] || continue
  if (( first )); then cat "$f" > "$ALL_OVERALL"; first=0; else tail -n +2 "$f" >> "$ALL_OVERALL"; fi
done
first=1
for f in "$PER_CHR_DIR"/*.per_chrom.tsv; do
  [[ -e "$f" ]] || continue
  if (( first )); then cat "$f" > "$ALL_PERCHR"; first=0; else tail -n +2 "$f" >> "$ALL_PERCHR"; fi
done

# Sort stably
python3 - "$ALL_OVERALL" "$ALL_PERCHR" <<'PY'
import sys, csv
def sort_tsv(path, keys):
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t"); rows = list(r); fn = r.fieldnames
    rows.sort(key=lambda x: tuple(x[k] for k in keys))
    with open(path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fn, delimiter="\t"); w.writeheader(); w.writerows(rows)
sort_tsv(sys.argv[1], ["sample"])
sort_tsv(sys.argv[2], ["chrom", "sample"])
PY

# ---- 7) Summary statistics ----
SUMMARY_OVERALL="${OUTDIR}/summary.coverage_context_qc.tsv"
SUMMARY_PERCHR="${OUTDIR}/summary.per_chrom.coverage_context_qc.tsv"

python3 - "$ALL_OVERALL" "$ALL_PERCHR" "$SUMMARY_OVERALL" "$SUMMARY_PERCHR" <<'PY'
import sys, csv
from collections import defaultdict

overall_tsv, perchr_tsv, summary_overall, summary_perchr = sys.argv[1:]

def is_num(x):
    try: float(x); return True
    except: return False

def qstats(vals):
    vals = sorted(vals)
    if not vals: return None
    def q(p):
        k = max(1, int(round(p*len(vals)))); return vals[min(k, len(vals))-1]
    return {"n":len(vals),"min":vals[0],"p10":q(.1),"p25":q(.25),"p50":q(.5),"mean":sum(vals)/len(vals),"p75":q(.75),"p90":q(.9),"max":vals[-1]}

with open(overall_tsv) as f:
    r = csv.DictReader(f, delimiter="\t"); rows = list(r); fields = r.fieldnames
numeric_fields = [x for x in fields if x not in ("sample","bam")]
with open(summary_overall, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["metric","n","min","p10","p25","p50","mean","p75","p90","max"])
    for m in numeric_fields:
        vals = [float(row[m]) for row in rows if is_num(row[m])]
        qs = qstats(vals)
        if qs: w.writerow([m]+[qs[k] for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]])

by_chr = defaultdict(list)
with open(perchr_tsv) as f:
    r = csv.DictReader(f, delimiter="\t"); fields = r.fieldnames
    for row in r:
        for m in fields:
            if m in ("sample","bam","chrom"): continue
            if is_num(row[m]): by_chr[(row["chrom"],m)].append(float(row[m]))
with open(summary_perchr, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["chrom","metric","n","min","p10","p25","p50","mean","p75","p90","max"])
    for (ch,m), vals in sorted(by_chr.items()):
        qs = qstats(vals)
        if qs: w.writerow([ch,m]+[qs[k] for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]])
PY

# ---- Write .results ----
RESULTSFILE="${OUTDIR}/${STEP}.results"
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "all_overall\t${ALL_OVERALL}\tPer-sample overall coverage QC"
  echo -e "all_perchr\t${ALL_PERCHR}\tPer-sample per-chromosome coverage QC"
  echo -e "summary_overall\t${SUMMARY_OVERALL}\tSummary statistics across samples"
  echo -e "summary_perchr\t${SUMMARY_PERCHR}\tPer-chromosome summary statistics"
  echo -e "n_samples\t${NSAMPLES}\tSamples processed"
} > "$RESULTSFILE"

echo "[DONE] $ALL_OVERALL"
echo "[DONE] $SUMMARY_OVERALL"
echo "[DONE] $ALL_PERCHR"
echo "[DONE] $SUMMARY_PERCHR"
