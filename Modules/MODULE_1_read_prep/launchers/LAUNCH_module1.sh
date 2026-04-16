#!/usr/bin/env bash
# =============================================================================
# LAUNCH_module1.sh — Master runner for MODULE_1 (reads + BAMs)
# =============================================================================
# Two tracks:
#   Track A (read-level): fastp QC, species verification, inventory checks
#     STEP_A01 (SLURM) → A02 → A03 → A04 (SLURM) → A05 → A06 → A07 → A08
#   Track B (BAM-level): align, merge/markdup, TLEN, popgen filter, depth QC
#     SLURM_B01 → SLURM_B02 → B03 → B04 → SLURM_B05 → SLURM_B06 → B07 → B08
#
# This launcher is a skeleton / dispatcher — the heavy work is in each step
# or SLURM wrapper. SLURM steps are expected to be submitted with `sbatch`;
# this launcher echoes the sbatch line rather than submitting, so it's safe
# to run on a login node.
#
# Usage:
#   bash launchers/LAUNCH_module1.sh                # Show the full plan
#   bash launchers/LAUNCH_module1.sh --step A03     # Run a single step
#   bash launchers/LAUNCH_module1.sh --from A05     # Start from step A05
#   bash launchers/LAUNCH_module1.sh --track A      # Track A only
#   bash launchers/LAUNCH_module1.sh --track B      # Track B only
#
# CONFIRM: MODULE_1 did not previously have a config or launcher. This is
# a new skeleton built to match other modules' structure. The step scripts
# themselves may still carry internal self-references to their old S## names
# — that does not affect correctness but is a followup cleanup.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${SCRIPT_DIR}/00_module1_config.sh"
[[ -f "$CONFIG" ]] || { echo "Config not found: $CONFIG" >&2; exit 1; }
source "$CONFIG"

STEPS="${STEPS_DIR}"
SLURM="${SLURM_DIR}"

STEP=""
FROM=""
TRACK=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --step)  STEP="$2"; shift 2 ;;
    --from)  FROM="$2"; shift 2 ;;
    --track) TRACK="$2"; shift 2 ;;
    -h|--help)
      sed -n '3,22p' "$0"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

should_run() {
  local id="$1" track="${id:0:1}"
  [[ -n "$TRACK" && "$track" != "$TRACK" ]] && return 1
  [[ -n "$STEP"  && "$id"    != "$STEP"  ]] && return 1
  [[ -n "$FROM"  && "$id" < "$FROM"      ]] && return 1
  return 0
}

announce() {
  local id="$1" desc="$2"
  echo ""
  echo "================================================================"
  echo "  [$id]  $desc"
  echo "================================================================"
}

echo "=== MODULE_1 pipeline plan ==="
echo "Steps/ : $STEPS"
echo "Slurm/ : $SLURM"
echo ""

# ── Track A: read-level ──────────────────────────────────────────────────────
if should_run A01; then
  announce A01 "fastp trim (SLURM)"
  echo "    sbatch $SLURM/SLURM_A01_fastp_trim_reads.sh"
fi
if should_run A02; then
  announce A02 "fastp EOF / redo table"
  bash "$STEPS/STEP_A02_fastp_eof_redo_table.sh"
fi
if should_run A03; then
  announce A03 "extract fastp QC"
  bash "$STEPS/STEP_A03_extract_fastp_qc.sh"
fi
if should_run A04; then
  announce A04 "species assign (Mash, SLURM)"
  echo "    sbatch $SLURM/SLURM_A04_species_assign_mash.sh"
fi
if should_run A05; then
  announce A05 "species assign (kmers)"
  bash "$STEPS/STEP_A05_species_assign_kmers.sh"
fi
if should_run A06; then
  announce A06 "plot species assignment"
  python3 "$STEPS/STEP_A06_plot_species_assign.py"
fi
if should_run A07; then
  announce A07 "inventory fastp vs BAM"
  bash "$STEPS/STEP_A07_inventory_fastp_vs_bam.sh"
fi
if should_run A08; then
  announce A08 "check BAM/BAI pairs"
  bash "$STEPS/STEP_A08_check_bam_bai_pairs.sh"
fi

# ── Track B: BAM-level ───────────────────────────────────────────────────────
if should_run B01; then
  announce B01 "minimap2 alignment (SLURM)"
  echo "    sbatch --export=ALL,BATCH=00 $SLURM/SLURM_B01_map_minimap2.sh"
fi
if should_run B02; then
  announce B02 "merge / markdup / clipOverlap (SLURM array)"
  echo "    sbatch $SLURM/SLURM_B02_merge_markdup_clip.sh"
fi
if should_run B03; then
  announce B03 "TLEN percentiles"
  bash "$STEPS/STEP_B03_get_tlen_percentiles.sh"
fi
if should_run B04; then
  announce B04 "insert size stats"
  bash "$STEPS/STEP_B04_get_insert_size_stats.sh"
fi
if should_run B05; then
  announce B05 "popgen BAM filter (SLURM)"
  echo "    sbatch $SLURM/SLURM_B05_filter_bam_popgen.sh"
fi
if should_run B06; then
  announce B06 "depth QC (mosdepth, SLURM)"
  echo "    sbatch $SLURM/SLURM_B06_qc_depth_table.sh"
fi
if should_run B07; then
  announce B07 "summarize BAM QC"
  bash "$STEPS/STEP_B07_summarize_bam_qc.sh"
fi
if should_run B08; then
  announce B08 "build provenance manifest"
  bash "$STEPS/STEP_B08_make_bam_provenance.sh"
fi

echo ""
echo "=== MODULE_1 complete ==="
