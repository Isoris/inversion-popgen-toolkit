#!/usr/bin/env bash
###############################################################################
# LAUNCH_module2a.sh — MODULE_2A / STEP_1: Variants & Structure Base
#
# ###########################################################################
# MODULE
#   MODULE_2A_ModernDNA_Variants_Ancestry / STEP_1_Variants_StructureBase
#
# PURPOSE
#   Main user-facing wrapper for the modern DNA variant + structure foundation.
#   Centralizes: mask/SFS preparation, biSNP site discovery, panel generation,
#   BEAGLE generation, PCAngsd, NGSadmix, ngsRelate, evalAdmix, relatedness
#   pruning, and canonical best-seed-by-K selection.
#
# ###########################################################################
# AVAILABLE SUBCOMMANDS
# ###########################################################################
#
#   init
#       Validate config, reference, BAM list, sample list, tools.
#
#   masks_sfs
#       Build callable mask, RF/chrom lists, per-chunk SAF, global merged
#       SAF, folded SFS, and mean pest prior.
#       Internally calls: steps/STEP_A01_masks_sfs.sh
#
#   call_bisnps
#       Run chunked ANGSD biSNP discovery using pest prior.
#       Internally calls: steps/STEP_A02_bisnps.sh
#
#   build_panels
#       Merge + thin sites at 200/500/1000 bp and 5/10/25 kb. ANGSD-index.
#       Internally calls: steps/STEP_A03_panels.sh
#
#   make_beagles
#       Generate BEAGLE GL files per-RF and whole-genome for each panel.
#       Merge per-RF BEAGLEs into whole-genome BEAGLEs.
#       Internally calls: steps/STEP_A04_beagles.sh
#
#   structure_all
#       Run PCAngsd + NGSadmix + evalAdmix + best-seed-by-K on ALL samples.
#       Internally calls: steps/STEP_A06_structure.sh --samples all_samples.txt
#
#   relatedness
#       Run ngsRelate, NAToRA multi-cutoff, greedy first-degree pruning,
#       and 3-panel relatedness figure. Produces pruned_samples.txt.
#       Internally calls: steps/STEP_A07_relatedness.sh
#
#   structure_pruned
#       Same as structure_all but on PRUNED sample set (after relatedness).
#       Internally calls: steps/STEP_A06_structure.sh --samples pruned_samples.txt
#
#   merge_summaries
#       Merge best-seed-by-K tables from all + pruned into one combined
#       table with a sample_set column.
#       Internally calls: steps/STEP_A08_merge_structure_summaries.sh
#
#   all
#       Run the full pipeline from current stage (skips completed steps).
#
# ###########################################################################
# CANONICAL STEP 1 OUTPUTS (produced by 'select_best')
# ###########################################################################
#
#   all_seed_metrics_by_K.tsv
#       Every K × seed run with loglik, evalAdmix residuals, file paths,
#       usability flags.
#
#   best_seed_by_K.tsv
#       One row per K: selected best seed, metrics, file paths, selection
#       rule, sample count validation.
#
#   best_seed_copied_files.tsv
#       Provenance: which source files were copied/linked as *_best.
#
#   cluster_palette_by_K.tsv
#       Stable 12-color Tableau-derived palette, one row per K × cluster.
#
#   sample_order_reference.tsv
#       Canonical sample ordering matching BAM list / BEAGLE column order.
#
#   sample_main_ancestry_by_K.tsv
#       Per-sample: dominant cluster, full Q vector, color, for each K.
#       This is the primary table downstream modules consume.
#
# ###########################################################################
# FULL CLI ARGUMENTS
# ###########################################################################
#
#   --config <file>     Config file. Default: ./config.sh
#   --task <name>       Alternative to positional subcommand.
#   --dry-run           Print commands without executing.
#   --force             Rebuild outputs even if existing.
#   --resume            Skip completed outputs.
#   --thin <W>          Override default thinning for ancestry/relatedness.
#   --kmin <int>        Override K_MIN. Default: 2
#   --kmax <int>        Override K_MAX. Default: 12
#   --seeds <list>      Override seeds (comma-separated). Default: 1,2,3
#   --threads <int>     Override thread count.
#   --scope <name>      Analysis scope for select_best: global, chromosome,
#                       or a specific thin label.
#
# ###########################################################################
# DEFAULT DECISION RULES
# ###########################################################################
#
#   Best seed selection (within each scope):
#     1. Highest log-likelihood
#     2. Lowest mean |evalAdmix off-diagonal residual|
#     3. Lowest max |evalAdmix residual|
#     4. Lowest seed number as tie-break
#
#   If evalAdmix outputs are missing for some seeds, selection proceeds
#   on loglik alone with a warning flag.
#
#   If --scope is not given, select_best defaults to global thin-500.
#
# ###########################################################################
# INTERNAL TASK MAP
# ###########################################################################
#
#   masks_sfs    -> steps/STEP_A01_masks_sfs.sh
#                   (mask_regions_from_fasta.py, angsd sites index,
#                    make_chr_rf_chunk_list, angsd -doSaf, realSFS cat,
#                    realSFS -fold -bootstrap)
#
#   call_bisnps  -> steps/STEP_A02_bisnps.sh
#                   (angsd -doMaf -SNP_pval -pest per chunk)
#
#   build_panels -> steps/STEP_A03_panels.sh
#                   (merge mafs.gz, awk thinning, angsd sites index)
#
#   make_beagles -> steps/STEP_A04_beagles.sh
#                   (angsd -doGlf 2 -doMajorMinor 3 per-RF + WG,
#                    merge per-RF to WG)
#
#   run_ancestry -> steps/STEP_A06_structure.sh
#                   (pcangsd per LG×K, NGSadmix per RF×K×seed)
#
#   run_evaladmix -> steps/STEP_A06_structure.sh (evalAdmix)
#                    (evalAdmix per K×seed, summarize_evaladmix.R)
#
#   run_relatedness -> steps/STEP_A07_relatedness.sh
#                      (ngsRelate, NAToRA multi-cutoff,
#                       prune_first_degree_pairs.py,
#                       plot_relatedness_3panel.R)
#
#   select_best  -> utils/select_best_seed_by_K.R
#                   (select_best_seed_by_K.R)
#
# ###########################################################################
# REQUIRED INPUTS
# ###########################################################################
#
#   - config.sh
#   - Reference FASTA + .fai
#   - Final filtered BAM list from MODULE_1_Reads
#   - Sample list (one per line, same order as BAM list)
#
# ###########################################################################
# DEBUG / COPY-PASTE TO CHATGPT
# ###########################################################################
#
#   If something breaks, copy-paste:
#     1. This entire header (LAUNCH_module2a.sh top section)
#     2. Your config.sh
#     3. The exact command you ran
#     4. The error output
#     5. Which subcommand failed
#
#   This header exposes the full control surface so an assistant can see
#   all available knobs and decide what to change without follow-up.
###############################################################################
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
STEPS="${SCRIPT_DIR}/../steps"
CONFIG="${SCRIPT_DIR}/../00_module2a_config.sh"

# ---- Parse args ----
TASK=""
DRY_RUN=0
FORCE=0
RESUME=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)   CONFIG="$2"; shift 2;;
    --task)     TASK="$2"; shift 2;;
    --dry-run)  DRY_RUN=1; shift;;
    --force)    FORCE=1; RESUME=0; shift;;
    --resume)   RESUME=1; shift;;
    --thin)     export RELATE_THIN="$2"; shift 2;;
    --kmin)     export K_MIN="$2"; shift 2;;
    --kmax)     export K_MAX="$2"; shift 2;;
    --seeds)    IFS=',' read -ra SEEDS <<< "$2"; export SEEDS; shift 2;;
    --threads)  export DEFAULT_THREADS="$2"; shift 2;;
    --scope)    export SCOPE="$2"; shift 2;;
    -h|--help)  head -n 120 "$0" | grep -E '^\s*#' | sed 's/^# \?//'; exit 0;;
    *)          TASK="${TASK:-$1}"; shift;;
  esac
done

[[ -n "$TASK" ]] || { echo "Usage: $0 <subcommand> [options]"; echo "Subcommands: init masks_sfs call_bisnps build_panels make_beagles structure_all relatedness structure_pruned merge_summaries all"; exit 1; }

# ---- Load config ----
[[ -f "$CONFIG" ]] || { echo "[ERROR] Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

export DRY_RUN FORCE RESUME SCRIPT_DIR STEPS

timestamp(){ date '+%F %T'; }
export -f timestamp

echo "[$(timestamp)] LAUNCH_module2a.sh task=${TASK}"
echo "[$(timestamp)] config=${CONFIG}"

# ---- Dispatch ----
case "$TASK" in
  init)
    echo "[$(timestamp)] Validating config..."
    [[ -f "$REF" ]]        || { echo "[ERROR] Missing REF: $REF" >&2; exit 1; }
    [[ -f "$BAMLIST" ]]    || { echo "[ERROR] Missing BAMLIST: $BAMLIST" >&2; exit 1; }
    [[ -f "$SAMPLE_LIST" ]] || { echo "[ERROR] Missing SAMPLE_LIST: $SAMPLE_LIST" >&2; exit 1; }
    command -v angsd >/dev/null 2>&1 || echo "[WARN] angsd not in PATH"
    command -v realSFS >/dev/null 2>&1 || echo "[WARN] realSFS not in PATH"
    command -v pcangsd >/dev/null 2>&1 || echo "[WARN] pcangsd not in PATH"
    command -v NGSadmix >/dev/null 2>&1 || echo "[WARN] NGSadmix not in PATH"
    command -v ngsRelate >/dev/null 2>&1 || echo "[WARN] ngsRelate not in PATH"
    N_BAM=$(wc -l < "$BAMLIST")
    N_SAMP=$(wc -l < "$SAMPLE_LIST")
    echo "[OK] REF=$REF"
    echo "[OK] BAMLIST=$BAMLIST ($N_BAM BAMs)"
    echo "[OK] SAMPLE_LIST=$SAMPLE_LIST ($N_SAMP samples)"
    echo "[OK] Config validated."
    ;;

  masks_sfs)
    bash "${STEPS}/STEP_A01_masks_sfs.sh"
    ;;

  call_bisnps)
    bash "${STEPS}/STEP_A02_bisnps.sh"
    ;;

  build_panels)
    bash "${STEPS}/STEP_A03_panels.sh"
    ;;

  make_beagles)
    bash "${STEPS}/STEP_A04_beagles.sh"
    ;;

  structure_all)
    bash "${STEPS}/STEP_A06_structure.sh" --samples "${SAMPLE_LIST}"
    ;;

  relatedness)
    bash "${STEPS}/STEP_A07_relatedness.sh"
    ;;

  structure_pruned)
    PRUNED="${THIN_DIR}/06_relatedness/pruned_samples.txt"
    [[ -s "$PRUNED" ]] || { echo "[ERROR] Run 'relatedness' first to produce pruned_samples.txt" >&2; exit 1; }
    bash "${STEPS}/STEP_A06_structure.sh" --samples "${PRUNED}"
    ;;

  merge_summaries)
    bash "${STEPS}/STEP_A08_merge_structure_summaries.sh"
    ;;

  all)
    echo "[$(timestamp)] Running full pipeline..."
    for step in masks_sfs call_bisnps build_panels make_beagles structure_all relatedness structure_pruned merge_summaries; do
      echo ""
      echo "[$(timestamp)] === ${step} ==="
      "$0" --config "$CONFIG" "${step}"
    done
    echo "[$(timestamp)] Full pipeline complete."
    ;;

  *)
    echo "[ERROR] Unknown task: $TASK" >&2
    echo "Available: init masks_sfs call_bisnps build_panels make_beagles structure_all relatedness structure_pruned merge_summaries all"
    exit 1
    ;;
esac

echo "[$(timestamp)] Done: ${TASK}"
