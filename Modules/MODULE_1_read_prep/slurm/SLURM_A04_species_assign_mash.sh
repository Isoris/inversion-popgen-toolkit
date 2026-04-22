#!/usr/bin/env bash
###############################################################################
# SLURM_A04_species_assign_mash.slurm
#
# Builds Mash sketches for Gar and Mac references, sketches all fastp read
# pairs, assigns each run-lane unit by Mash distance in parallel, and writes
# both per-run and collapsed per-CGA species assignment tables.
#
# Usage:
#   sbatch SLURM_A04_species_assign_mash.slurm
#
# Outputs:
#   ${OUTDIR}/results.tsv              Per-run Mash distances + call
#   ${OUTDIR}/results.per_CGA.tsv      Collapsed per-CGA majority call
#   ${OUTDIR}/SLURM_A04_species_assign_mash.arg
#   ${OUTDIR}/SLURM_A04_species_assign_mash.results
#   ${OUTDIR}/SLURM_A04_species_assign_mash.metric.tsv
###############################################################################
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -c 127
#SBATCH --mem=237GB
#SBATCH -t 1-00:00:00
#SBATCH -J mash_assign
#SBATCH -o mash_assign.%j.out
#SBATCH -e mash_assign.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

STEP="SLURM_A04_species_assign_mash"
timestamp(){ date '+%F %T'; }

# ---- Paths ----
BASE="/scratch/lt200308-agbsci/Quentin_project"
FASTP_DIR="${BASE}/00-samples/fastp"
REF_GAR="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_MAC="${BASE}/00-samples/fClaHyb_Mac_LG.fa"

OUTDIR="${BASE}/00-samples/species_assign_mash_fastp"
LOGDIR="${OUTDIR}/logs"
TMPDIR="${OUTDIR}/tmp"
SAMPDIR="${OUTDIR}/mash_samples"
REFDIR="${OUTDIR}/mash_refs"

mkdir -p "$OUTDIR" "$LOGDIR" "$TMPDIR" "$SAMPDIR" "$REFDIR"

# ---- Mash parameters ----
MASH="$(command -v mash)"
K=31
SKETCH=50000
MINCOPIES=2
MARGIN=0.002

# ---- Parallelism ----
TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-127}"
THREADS_PER_SAMPLE="${THREADS_PER_SAMPLE:-8}"
NJOBS=$(( TOTAL_CPUS / THREADS_PER_SAMPLE ))
(( NJOBS < 1 )) && NJOBS=1

# ---- Sidecar paths ----
GAR_MSH="${REFDIR}/gar_ref.msh"
MAC_MSH="${REFDIR}/mac_ref.msh"
IDS="${TMPDIR}/ids.txt"
RESULTS="${OUTDIR}/results.tsv"
PER_CGA="${OUTDIR}/results.per_CGA.tsv"
ARGFILE="${OUTDIR}/${STEP}.arg"
RESULTSFILE="${OUTDIR}/${STEP}.results"
METRICFILE="${OUTDIR}/${STEP}.metric.tsv"

MASH_VERSION="$("$MASH" --version 2>/dev/null | head -n 1 || echo NA)"

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "slurm_job_id\t${SLURM_JOB_ID:-NA}"
  echo -e "tool_name\tmash"
  echo -e "tool_path\t${MASH}"
  echo -e "tool_version\t${MASH_VERSION}"
  echo -e "fastp_dir\t${FASTP_DIR}"
  echo -e "ref_gar\t${REF_GAR}"
  echo -e "ref_mac\t${REF_MAC}"
  echo -e "outdir\t${OUTDIR}"
  echo -e "kmer_size\t${K}"
  echo -e "mash_sketch\t${SKETCH}"
  echo -e "min_copies\t${MINCOPIES}"
  echo -e "margin\t${MARGIN}"
  echo -e "total_cpus\t${TOTAL_CPUS}"
  echo -e "threads_per_sample\t${THREADS_PER_SAMPLE}"
  echo -e "njobs\t${NJOBS}"
} > "$ARGFILE"

echo "[$(timestamp)] [INFO] Host=$(hostname) JOB_ID=${SLURM_JOB_ID:-NA}"
echo "[$(timestamp)] [INFO] MASH=$MASH ($MASH_VERSION)"
echo "[$(timestamp)] [INFO] K=$K SKETCH=$SKETCH MINCOPIES=$MINCOPIES MARGIN=$MARGIN"
echo "[$(timestamp)] [INFO] NJOBS=$NJOBS THREADS_PER_SAMPLE=$THREADS_PER_SAMPLE"

# ---- Build reference sketches ----
if [[ ! -s "$GAR_MSH" ]]; then
  echo "[$(timestamp)] Building Gar reference sketch..."
  "$MASH" sketch -k "$K" -s "$SKETCH" -p "$THREADS_PER_SAMPLE" -o "${REFDIR}/gar_ref" "$REF_GAR" \
    > "${LOGDIR}/mash_ref_gar.log" 2>&1
fi
if [[ ! -s "$MAC_MSH" ]]; then
  echo "[$(timestamp)] Building Mac reference sketch..."
  "$MASH" sketch -k "$K" -s "$SKETCH" -p "$THREADS_PER_SAMPLE" -o "${REFDIR}/mac_ref" "$REF_MAC" \
    > "${LOGDIR}/mash_ref_mac.log" 2>&1
fi

# ---- Discover all fastp pairs ----
: > "$IDS"
while IFS= read -r r1; do
  bn="$(basename "$r1")"
  id="${bn%.R1.fastp.fq.gz}"
  r2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"
  if [[ -s "$r2" ]]; then
    echo "$id" >> "$IDS"
  else
    echo "[$(timestamp)] [WARN] Missing R2 for $id" >&2
  fi
done < <(find "$FASTP_DIR" -type f -name "*.R1.fastp.fq.gz" | sort)

sort -u "$IDS" -o "$IDS"
N=$(wc -l < "$IDS")
echo "[$(timestamp)] Found $N run-lane units with pairs."
[[ "$N" -gt 0 ]] || { echo "[ERROR] No FASTP pairs found" >&2; exit 1; }

# ---- Per-sample worker ----
echo -e "id\tR1\tR2\tdist_gar\tdist_mac\tcall\tabs_diff" > "$RESULTS"

run_one() {
  local id="$1"
  local r1="${FASTP_DIR}/${id}.R1.fastp.fq.gz"
  local r2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"
  local msh="${SAMPDIR}/${id}.msh"
  local slog="${LOGDIR}/mash_${id}.log"

  if [[ ! -s "$msh" ]]; then
    "$MASH" sketch -r -k "$K" -s "$SKETCH" -m "$MINCOPIES" -p "$THREADS_PER_SAMPLE" \
      -o "${SAMPDIR}/${id}" "$r1" "$r2" > "$slog" 2>&1
  fi

  local dist_gar dist_mac
  dist_gar="$("$MASH" dist "$GAR_MSH" "$msh" | awk 'NR==1{print $3}')"
  dist_mac="$("$MASH" dist "$MAC_MSH" "$msh" | awk 'NR==1{print $3}')"

  local call
  call="$(awk -v a="$dist_gar" -v b="$dist_mac" -v m="$MARGIN" 'BEGIN{
    if (a + m < b) print "Gariepinus";
    else if (b + m < a) print "Macrocephalus";
    else print "Ambiguous";
  }')"

  local absdiff
  absdiff="$(awk -v a="$dist_gar" -v b="$dist_mac" 'BEGIN{d=a-b; if(d<0) d=-d; printf "%.6f", d}')"

  echo -e "${id}\t${r1}\t${r2}\t${dist_gar}\t${dist_mac}\t${call}\t${absdiff}"
}
export -f run_one timestamp
export FASTP_DIR SAMPDIR LOGDIR MASH K SKETCH MINCOPIES GAR_MSH MAC_MSH MARGIN THREADS_PER_SAMPLE

echo "[$(timestamp)] Running mash per id in parallel..."
xargs -a "$IDS" -n 1 -P "$NJOBS" bash -lc 'run_one "$0"' >> "$RESULTS"
echo "[$(timestamp)] Wrote per-run table: $RESULTS"

# ---- Collapse to per-CGA ----
python3 - <<PY
import csv, statistics
from collections import defaultdict, Counter

infile = "${RESULTS}"
outfile = "${PER_CGA}"

rows = []
with open(infile, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        rows.append(row)

by_cga = defaultdict(list)
for row in rows:
    cga = row["id"].split(".")[0]
    by_cga[cga].append(row)

def median(xs):
    xs = [float(x) for x in xs if x not in ("NA","")]
    return statistics.median(xs) if xs else float("nan")

with open(outfile, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["CGA_sample","n_runs","median_dist_gar","median_dist_mac","majority_call","calls_breakdown"])
    for cga, lst in sorted(by_cga.items()):
        n = len(lst)
        med_g = median([x["dist_gar"] for x in lst])
        med_m = median([x["dist_mac"] for x in lst])
        calls = [x["call"] for x in lst]
        c = Counter(calls)
        nonamb = [x for x in calls if x != "Ambiguous"]
        maj = Counter(nonamb).most_common(1)[0][0] if nonamb else "Ambiguous"
        breakdown = ",".join([f"{k}:{v}" for k,v in c.most_common()])
        w.writerow([cga, n, f"{med_g:.6f}", f"{med_m:.6f}", maj, breakdown])
PY
echo "[$(timestamp)] Wrote per-CGA table: $PER_CGA"

# ---- Write metric sidecar ----
{
  echo -e "sample\tmodule\tmetric\tvalue\tunit\tderived_from"
  awk -F'\t' 'NR>1{
    print $1 "\tspecies_assign_mash\tdist_gar\t" $4 "\tdistance\tresults.tsv";
    print $1 "\tspecies_assign_mash\tdist_mac\t" $5 "\tdistance\tresults.tsv";
    print $1 "\tspecies_assign_mash\tcall\t" $6 "\tclass\tresults.tsv";
    print $1 "\tspecies_assign_mash\tabs_diff\t" $7 "\tdistance\tresults.tsv";
  }' "$RESULTS"
  awk -F'\t' 'NR>1{
    print $1 "\tspecies_assign_mash_cga\tn_runs\t" $2 "\tcount\tresults.per_CGA.tsv";
    print $1 "\tspecies_assign_mash_cga\tmedian_dist_gar\t" $3 "\tdistance\tresults.per_CGA.tsv";
    print $1 "\tspecies_assign_mash_cga\tmedian_dist_mac\t" $4 "\tdistance\tresults.per_CGA.tsv";
    print $1 "\tspecies_assign_mash_cga\tmajority_call\t" $5 "\tclass\tresults.per_CGA.tsv";
  }' "$PER_CGA"
} > "$METRICFILE"

# ---- Write .results ----
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "per_run_results\t${RESULTS}\tPer-run-lane Mash distances and calls"
  echo -e "per_cga_results\t${PER_CGA}\tPer-CGA collapsed majority calls"
  echo -e "metric_file\t${METRICFILE}\tMetric sidecar"
  echo -e "gar_sketch\t${GAR_MSH}\tGar reference Mash sketch"
  echo -e "mac_sketch\t${MAC_MSH}\tMac reference Mash sketch"
  echo -e "n_run_lane_units\t${N}\tTotal run-lane units processed"
} > "$RESULTSFILE"

echo "[$(timestamp)] [DONE] Summary:"
cut -f6 "$RESULTS" | tail -n +2 | sort | uniq -c || true
cut -f5 "$PER_CGA" | tail -n +2 | sort | uniq -c || true
