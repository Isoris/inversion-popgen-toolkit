#!/usr/bin/env bash
###############################################################################
# S15_summarize_bam_qc.sh
#
# Parses per-sample samtools stats and flagstat outputs from merged BAMs,
# writes a per-sample QC table, and generates a summary report with depth
# distributions and heuristic ANGSD depth cutoff suggestions.
#
# Usage:
#   sbatch S15_summarize_bam_qc.sh
#   # or: bash S15_summarize_bam_qc.sh  (if run interactively)
#
# Outputs:
#   ${OUTDIR}/all_samples.bam_qc.tsv
#   ${OUTDIR}/summary.txt
#   ${OUTDIR}/S15_summarize_bam_qc.arg
#   ${OUTDIR}/S15_summarize_bam_qc.results
###############################################################################
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8GB
#SBATCH -t 01:00:00
#SBATCH -J bam_qc_sum
#SBATCH -o bam_qc_sum.%j.out
#SBATCH -e bam_qc_sum.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly >/dev/null 2>&1 || true

STEP="S15_summarize_bam_qc"
timestamp(){ date '+%F %T'; }

BASE="$(pwd)"
MERGED_DIR="${BASE}/02-merged_per_sample"
GENOME_SIZE=963905721    # fClaHyb_Gar_LG assembly size

OUTDIR="${MERGED_DIR}/results/popgen_qc"
mkdir -p "$OUTDIR"

ARGFILE="${OUTDIR}/${STEP}.arg"
RESULTSFILE="${OUTDIR}/${STEP}.results"

# ---- Write .arg ----
{
  echo -e "key\tvalue"
  echo -e "step\t${STEP}"
  echo -e "script_name\t$(basename "$0")"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "merged_dir\t${MERGED_DIR}"
  echo -e "genome_size\t${GENOME_SIZE}"
  echo -e "outdir\t${OUTDIR}"
  echo -e "depth_metric\tbases_mapped_cigar_div_genome_size"
  echo -e "heuristic_minDepthInd\t3_if_p10_depth_ge_6_else_2"
  echo -e "heuristic_maxDepthInd\tmin(60,round(p90_depth*5))"
} > "$ARGFILE"

python3 - "$MERGED_DIR" "$OUTDIR" "$GENOME_SIZE" <<'PY'
import os, re, sys, glob, statistics
from math import floor

merged_dir, outdir, genome_size = sys.argv[1], sys.argv[2], int(sys.argv[3])

def parse_samtools_stats(path):
    d = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            if not line.startswith("SN\t"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            key = parts[1].strip().rstrip(":")
            val = parts[2].strip()
            mapping = {
                "raw total sequences": ("raw_total_sequences", int),
                "reads properly paired": ("reads_properly_paired", int),
                "reads paired": ("reads_paired", int),
                "reads duplicated": ("reads_duplicated", int),
                "reads MQ0": ("reads_MQ0", int),
                "bases mapped (cigar)": ("bases_mapped_cigar", int),
                "error rate": ("error_rate", float),
            }
            if key in mapping:
                name, typ = mapping[key]
                try: d[name] = typ(val)
                except: pass
    return d

def parse_flagstat(path):
    d = {}
    patterns = [
        (re.compile(r"^(\d+)\s+\+\s+\d+\s+in total"), "total_reads", int, 1),
        (re.compile(r"^(\d+)\s+\+\s+\d+\s+primary$"), "primary", int, 1),
        (re.compile(r"^(\d+)\s+\+\s+\d+\s+supplementary$"), "supplementary", int, 1),
        (re.compile(r"^(\d+)\s+\+\s+\d+\s+duplicates$"), "duplicates", int, 1),
        (re.compile(r"^(\d+)\s+\+\s+\d+\s+mapped\s+\(([\d\.]+)%"), "mapped_pct", float, 2),
        (re.compile(r"^(\d+)\s+\+\s+\d+\s+properly paired\s+\(([\d\.]+)%"), "properly_paired_pct", float, 2),
        (re.compile(r"^(\d+)\s+\+\s+\d+\s+with mate mapped to a different chr$"), "mate_mapped_diff_chr", int, 1),
    ]
    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            for rx, name, typ, grp in patterns:
                m = rx.match(line)
                if m: d[name] = typ(m.group(grp)); break
    return d

def quantiles(xs):
    xs = sorted(xs)
    if not xs: return {}
    def q(p):
        k = max(1, int(round(p*len(xs)))); return xs[min(k, len(xs))-1]
    return {"min":xs[0],"p10":q(.1),"p25":q(.25),"p50":q(.5),"mean":sum(xs)/len(xs),"median":statistics.median(xs),"p75":q(.75),"p90":q(.9),"max":xs[-1],"n":len(xs)}

sample_dirs = sorted([d for d in glob.glob(os.path.join(merged_dir, "CGA*")) if os.path.isdir(d)])
rows, missing = [], []

for sdir in sample_dirs:
    sample = os.path.basename(sdir)
    flagstat_path = os.path.join(sdir, f"{sample}.filtered.flagstat.txt")
    if not os.path.exists(flagstat_path):
        flagstat_path = os.path.join(sdir, f"{sample}.flagstat.txt")
    stats_path = os.path.join(sdir, f"{sample}.samtools.stats.txt")

    if not os.path.exists(stats_path) and not os.path.exists(flagstat_path):
        missing.append(sample); continue

    st = parse_samtools_stats(stats_path) if os.path.exists(stats_path) else {}
    fs = parse_flagstat(flagstat_path) if os.path.exists(flagstat_path) else {}

    bases = st.get("bases_mapped_cigar")
    depth = bases / genome_size if bases else None
    raw = st.get("raw_total_sequences")
    dup = st.get("reads_duplicated")
    mq0 = st.get("reads_MQ0")

    rows.append({
        "sample": sample, "mean_depth_x": depth, "bases_mapped_cigar": bases,
        "raw_total_sequences": raw, "reads_duplicated": dup,
        "dup_frac": dup/raw if raw and dup is not None else None,
        "reads_MQ0": mq0, "mq0_frac": mq0/raw if raw and mq0 is not None else None,
        "error_rate": st.get("error_rate"),
        "filtered_total_reads": fs.get("total_reads"),
        "filtered_properly_paired_pct": fs.get("properly_paired_pct"),
        "filtered_mapped_pct": fs.get("mapped_pct"),
        "filtered_mate_diff_chr": fs.get("mate_mapped_diff_chr"),
    })

# Write per-sample TSV
tsv_path = os.path.join(outdir, "all_samples.bam_qc.tsv")
cols = ["sample","mean_depth_x","bases_mapped_cigar","raw_total_sequences","reads_duplicated",
        "dup_frac","reads_MQ0","mq0_frac","error_rate","filtered_total_reads",
        "filtered_properly_paired_pct","filtered_mapped_pct","filtered_mate_diff_chr"]
with open(tsv_path, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in rows:
        out.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c] if r[c] is not None else "") for c in cols) + "\n")

# Summary + ANGSD depth heuristics
depths = [r["mean_depth_x"] for r in rows if r["mean_depth_x"] is not None]
qd = quantiles(depths)

rec_min = 3 if qd.get("p10") is not None and qd["p10"] >= 6.0 else 2
rec_max = min(60, int(round(qd["p90"] * 5))) if qd.get("p90") else None

summary_path = os.path.join(outdir, "summary.txt")
with open(summary_path, "w") as out:
    out.write(f"Genome size (bp): {genome_size}\n")
    out.write(f"Samples parsed: {len(rows)}\n")
    if missing:
        out.write(f"Missing stats+flagstat: {len(missing)}\n")
    out.write("\n=== Mean depth (X) ===\n")
    if qd:
        for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]:
            out.write(f"{k}\t{qd[k]:.6g}\n" if isinstance(qd[k], float) else f"{k}\t{qd[k]}\n")
    out.write(f"\n=== Suggested ANGSD depth cutoffs (heuristic) ===\n")
    out.write(f"setMinDepthInd\t{rec_min}\n")
    out.write(f"setMaxDepthInd\t{rec_max if rec_max else ''}\n")
    out.write("\nNotes:\n")
    out.write("- Depth from samtools stats 'bases mapped (cigar)' / genome_size.\n")
    out.write("- Properly paired % from filtered.flagstat when available.\n")

print(f"[DONE] {tsv_path}")
print(f"[DONE] {summary_path}")
PY

# ---- Write .results ----
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "bam_qc_table\t${OUTDIR}/all_samples.bam_qc.tsv\tPer-sample BAM QC metrics"
  echo -e "summary\t${OUTDIR}/summary.txt\tSummary statistics and ANGSD depth heuristics"
} > "$RESULTSFILE"

echo "[DONE] $ARGFILE  $RESULTSFILE"
