#!/usr/bin/env bash
###############################################################################
# helpers/06_relatedness.sh — ngsRelate + NAToRA + pruning + plotting
#
# Runs between structure_all and structure_pruned:
#   1) ngsRelate on thin-500 whole-genome BEAGLE
#   2) NAToRA multi-cutoff culling (using config low/high pairs)
#   3) Greedy first-degree pruning
#   4) 3-panel relatedness figure
#
# Produces:
#   - ngsRelate pairwise output
#   - NAToRA cutoff summary table
#   - pruned_samples.txt (input for structure_pruned)
#   - 3-panel relatedness figure
#
# Called by: run_step1.sh relatedness
###############################################################################
set -euo pipefail
source "$(dirname "$0")/../config.sh"

OUTDIR="${THIN_DIR}/06_relatedness"
mkdir -p "${OUTDIR}"

BEAGLE="${THIN_DIR}/03_merged_beagle/catfish.wholegenome.byRF.thin_${RELATE_THIN}.beagle.gz"
MAF_DIR="${GLOBAL_DIR}/02_snps"

[[ -s "$BEAGLE" ]] || { echo "[ERROR] Missing BEAGLE: $BEAGLE" >&2; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "[ERROR] Missing SAMPLE_LIST" >&2; exit 1; }

ARGFILE="${OUTDIR}/06_relatedness.arg"
{
  echo -e "key\tvalue"
  echo -e "step\t06_relatedness"
  echo -e "datetime\t$(timestamp)"
  echo -e "host\t$(hostname)"
  echo -e "beagle\t${BEAGLE}"
  echo -e "sample_list\t${SAMPLE_LIST}"
  echo -e "n_samples\t${N_SAMPLES}"
  echo -e "relate_thin\t${RELATE_THIN}"
  echo -e "natora_dup_mz\t${NATORA_DUP_MZ_LOW},${NATORA_DUP_MZ_HIGH}"
  echo -e "natora_first\t${NATORA_FIRST_LOW},${NATORA_FIRST_HIGH}"
  echo -e "natora_second\t${NATORA_SECOND_LOW},${NATORA_SECOND_HIGH}"
  echo -e "natora_third\t${NATORA_THIRD_LOW},${NATORA_THIRD_HIGH}"
  echo -e "first_degree_prune_theta\t${THETA_FIRST_DEGREE}"
} > "$ARGFILE"

RES="${OUTDIR}/catfish_${N_SAMPLES}_relatedness.res"
P="${DEFAULT_THREADS}"

# ===========================================================================
# 1) ngsRelate
# ===========================================================================
echo "[$(timestamp)] === ngsRelate ==="

if [[ -s "$RES" && "${FORCE:-0}" -eq 0 ]]; then
  echo "[SKIP] ngsRelate output exists: $RES"
else
  # Build frequency file matched to beagle site order
  echo "[$(timestamp)] Building frequency file..."
  zcat ${MAF_DIR}/catfish.*.mafs.gz \
    | awk 'BEGIN{OFS="\t"} $1!="chromo" {print $1"_"$2, $6}' \
    | sort -k1,1 > "${OUTDIR}/all_mafs_freq.tmp"

  zcat "$BEAGLE" | tail -n +2 | cut -f1 > "${OUTDIR}/beagle_sites.tmp"
  NSITES=$(wc -l < "${OUTDIR}/beagle_sites.tmp")

  awk 'NR==FNR {freq[$1]=$2; next} {if($1 in freq) print freq[$1]; else print "NA"}' \
    "${OUTDIR}/all_mafs_freq.tmp" "${OUTDIR}/beagle_sites.tmp" \
    > "${OUTDIR}/freq_for_ngsrelate.txt"

  echo "[$(timestamp)] Running ngsRelate (${N_SAMPLES} samples, ${NSITES} sites)..."
  ngsRelate \
    -G "$BEAGLE" \
    -f "${OUTDIR}/freq_for_ngsrelate.txt" \
    -n "${N_SAMPLES}" \
    -O "$RES" \
    -p "$P" \
    -m 1 \
    -z "$SAMPLE_LIST"

  rm -f "${OUTDIR}/all_mafs_freq.tmp" "${OUTDIR}/beagle_sites.tmp"
  echo "[$(timestamp)] ngsRelate done: $RES"
fi

# ===========================================================================
# 2) NAToRA multi-cutoff
# ===========================================================================
echo ""
echo "[$(timestamp)] === NAToRA multi-cutoff ==="

NATORA_PY="$(dirname "$0")/NAToRA_Public.py"
NATORA_INPUT="${OUTDIR}/catfish_for_natora.txt"
NATORA_SUMMARY="${OUTDIR}/natora_cutoff_summary.tsv"

[[ -f "$NATORA_PY" ]] || { echo "[WARN] NAToRA not found: $NATORA_PY — skipping"; }

if [[ -f "$NATORA_PY" ]]; then
  # Build NAToRA input: col3=ida, col4=idb, col18=theta
  awk 'BEGIN{OFS="\t"} NR>1 {print $3, $4, $18}' "$RES" > "$NATORA_INPUT"
  N_TOTAL=$(awk '{print $1"\n"$2}' "$NATORA_INPUT" | sort -u | wc -l)

  printf "group\tcutoff_label\ttheta_cutoff\tremoved_n\tretained_n\tremoved_pct\toutput_prefix\n" > "$NATORA_SUMMARY"

  run_cutoff() {
    local GRP="$1" LBL="$2" CUT="$3"
    local PFX="${OUTDIR}/natora_${GRP}_${LBL}_c${CUT}"
    python "$NATORA_PY" -i "$NATORA_INPUT" -o "$PFX" -c "$CUT" -v 0.0 -V 0.5 -e 1 2>/dev/null
    local RM=0; [[ -s "${PFX}_toRemove.txt" ]] && RM=$(sort -u "${PFX}_toRemove.txt" | wc -l)
    local RET=$((N_TOTAL - RM))
    local PCT; PCT=$(awk -v r="$RM" -v n="$N_TOTAL" 'BEGIN{printf "%.2f", 100*r/n}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$GRP" "$LBL" "$CUT" "$RM" "$RET" "$PCT" "$PFX" >> "$NATORA_SUMMARY"
  }

  run_cutoff "duplicate_mz"  "low"  "${NATORA_DUP_MZ_LOW}"
  run_cutoff "duplicate_mz"  "high" "${NATORA_DUP_MZ_HIGH}"
  run_cutoff "first_degree"  "low"  "${NATORA_FIRST_LOW}"
  run_cutoff "first_degree"  "high" "${NATORA_FIRST_HIGH}"
  run_cutoff "second_degree" "low"  "${NATORA_SECOND_LOW}"
  run_cutoff "second_degree" "high" "${NATORA_SECOND_HIGH}"
  run_cutoff "third_degree"  "low"  "${NATORA_THIRD_LOW}"
  run_cutoff "third_degree"  "high" "${NATORA_THIRD_HIGH}"

  echo "[$(timestamp)] NAToRA summary: $NATORA_SUMMARY"
  column -t -s $'\t' "$NATORA_SUMMARY" 2>/dev/null || cat "$NATORA_SUMMARY"
fi

# ===========================================================================
# 3) Greedy first-degree pruning
# ===========================================================================
echo ""
echo "[$(timestamp)] === Greedy first-degree pruning ==="

PRUNE_PY="$(dirname "$0")/prune_first_degree_pairs.py"
PRUNE_PREFIX="${OUTDIR}/first_degree_prune"

if [[ -f "$PRUNE_PY" && -s "$NATORA_INPUT" ]]; then
  python3 "$PRUNE_PY" \
    -i "$NATORA_INPUT" \
    -o "$PRUNE_PREFIX" \
    -t "${THETA_FIRST_DEGREE}"

  # Create the canonical pruned_samples.txt for structure_pruned
  PRUNED_SAMPLES="${OUTDIR}/pruned_samples.txt"
  cp "${PRUNE_PREFIX}_toKeep.txt" "$PRUNED_SAMPLES"
  echo "[$(timestamp)] Pruned sample list: $PRUNED_SAMPLES ($(wc -l < "$PRUNED_SAMPLES") samples)"
else
  echo "[WARN] Skipping greedy pruning (missing script or input)"
fi

# ===========================================================================
# 4) Relatedness figure
# ===========================================================================
echo ""
echo "[$(timestamp)] === Relatedness figure ==="

PLOT_R="$(dirname "$0")/plot_relatedness_3panel.R"
if [[ -f "$PLOT_R" && -s "$RES" ]]; then
  Rscript "$PLOT_R" \
    --input "$RES" \
    --outdir "${OUTDIR}/figures" \
    --prefix "catfish_relatedness" \
    2>/dev/null || echo "[WARN] Relatedness plot had issues"
else
  echo "[WARN] Skipping plot (missing R script or ngsRelate output)"
fi

# ===========================================================================
# 5) Write .results
# ===========================================================================
RESULTSFILE="${OUTDIR}/06_relatedness.results"
{
  echo -e "key\tvalue\tdescription"
  echo -e "arg_file\t${ARGFILE}\tParameter record"
  echo -e "ngsrelate_res\t${RES}\tPairwise relatedness output"
  echo -e "natora_summary\t${NATORA_SUMMARY}\tNAToRA cutoff comparison"
  echo -e "pruned_samples\t${OUTDIR}/pruned_samples.txt\tPruned sample list for structure_pruned"
  echo -e "prune_decisions\t${PRUNE_PREFIX}_decisions.tsv\tPer-pair pruning decisions"
  echo -e "prune_summary\t${PRUNE_PREFIX}_summary.tsv\tPruning summary stats"
} > "$RESULTSFILE"

echo ""
echo "[$(timestamp)] [DONE] 06_relatedness"
echo "[$(timestamp)] Next: bash helpers/05_structure.sh --samples ${OUTDIR}/pruned_samples.txt"
