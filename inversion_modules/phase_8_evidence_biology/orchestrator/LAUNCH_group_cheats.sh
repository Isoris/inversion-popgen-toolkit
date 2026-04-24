#!/usr/bin/env bash
# =============================================================================
# LAUNCH_group_cheats.sh — phase 4d production launcher
# =============================================================================
# Runs cheat30 (GDS), cheat6 (ancestry jackknife), cheat27/28/29 (mechanism),
# and conditionally the Fst age sub-block + burden regression — gating each
# on the corresponding group validation level.
#
# Called by run_phase4.sh after 4c_hypothesis completes.
#
# Usage:
#   sbatch LAUNCH_group_cheats.sh --chroms LG01,LG02,...
#
# FIX 50 (chat 10, 2026-04-17): promoted from LAUNCH_group_cheats_example.sh.
# Pass 15 (2026-04-24): scripts regrouped by question under phase_8_evidence_biology/.
# Pass 17 (2026-04-24): cheat6 moved out of q7_existence_audit (it's a
# group-trust robustness test, not an existence test); cross_species_bridge
# moved out of q5_age_and_origin (it's a cross-species conservation question,
# not an age question, and it writes Q3 flank keys in addition to Q5).
# Pass 18 (2026-04-24): cheat14 (minimap2) + cheat27 (BISER2) merged into
# sd_substrate.R per docs/toolkit_audit.md §D.1. cheat27's launcher call
# was a silent no-op (no CLI parser); now properly wired. cheat14 restored
# from archive and folded in as Angle A. Use NO_MINIMAP2=1 to skip the
# expensive angle.
# Pass 19 (2026-04-24): sd_substrate.R split back into three scripts per
# user request (independent provenance + separate blocks):
#   sd_substrate/sd_substrate_minimap2.R     — Angle A (trusted), BAM output
#   sd_substrate/sd_substrate_biser2.R       — Angle B, catalog lookup
#   sd_substrate/sd_substrate_concordance.R  — joint verdict
# minimap2 preset changed from asm5 to asm10 (catches SDs ≥90% identity).
# Adds --secondary=yes -N 50 -p 0.1 for SD pair detection. Outputs
# sorted/indexed BAM + strand-colored BED9 (red=inv, green=direct,
# blue=breakpoint) for IGV inspection. Everything persists under
# <cid>/sd_substrate/ — no tmp/.
# Pass 21 (2026-04-24): STEP_03_per_read_evidence added — closes the
# Q7B flat-key extraction gap. Row-level BAM ledger (read names + tags +
# 5'/3' clip orientation) at <cid>/per_read_evidence/; VCF sidecar for
# DELLY/Manta aggregates. Summary block exports q7b_bam_n_{split,discord,
# soft}_{left,right}_{INV,HET,REF} keys to the registry.
#
# Path corrections:
#   ORIGINAL EXAMPLE                  PRE-PASS-15             POST-PASS-15            POST-PASS-17/18 (current)
#   --------------------------------  ----------------------  ----------------------  ---------------------------
#   cheats/cheat30_gds_by_genotype.R  4f_group_dependent/     q5_age_and_origin/      q5_age_and_origin/ (same)
#   cheats/cheat6_ancestry_*.R        4f_group_dependent/     q7_existence_audit/     q6_group_robustness/
#   cheats/cheat27_sd_nahr_*.R        (lib-only, no CLI)      q4_mechanism/cheat27_*  q4_mechanism/sd_substrate.R
#   cheats/cheat14_repeat_arch.R      (archive only)          (archive only)          q4_mechanism/sd_substrate.R
#   (cheat28 / 29 — added in FIX 50)                          q4_mechanism/           q4_mechanism/ (same)
#   cross_species_bridge_v6.py        (existed at top level)  q5_age_and_origin/      cross_species/STEP_03A_cross_species_bridge.py
#   scripts/compute_age_fst_subblock  (ASPIRATIONAL)          still aspirational      still aspirational
#   burden/STEP_C01f_c_burden_*       (ASPIRATIONAL)          still aspirational      still aspirational
#
# Also added: cheat27 + cheat28 + cheat29 calls (these 3 exist in q4_mechanism/
# but the example launcher never called them). They require NONE validation
# level (mechanism analysis is group-independent), so no gate is needed.
#
# NOTE: The two aspirational scripts (compute_age_fst_subblock.R,
# STEP_C01f_c_burden_regression.R) are left in place as commented-out calls
# so that when the scripts arrive in the tarball, uncommenting is the only
# action needed. Their validation gates (SUPPORTED / VALIDATED) remain
# active if-guarded via registry_check_validation.
# =============================================================================

#SBATCH -A lt200308
#SBATCH -J 4d_cheats
#SBATCH -t 04:00:00
#SBATCH -c 4
#SBATCH --mem=16G

set -euo pipefail

# ── Parse ────────────────────────────────────────────────────────────────────
CHROMS_CSV=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_CSV="$2"; shift 2 ;;
    *) shift ;;
  esac
done

# ── Config + registry paths ──────────────────────────────────────────────────
BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
MODULE_DIR="${BASE}/inversion-popgen-toolkit"
PHASE8_DIR="${MODULE_DIR}/inversion_modules/phase_8_evidence_biology"

# shellcheck disable=SC1091
source "${MODULE_DIR}/00_inversion_config.sh" 2>/dev/null || true
# shellcheck disable=SC1091
source "${MODULE_DIR}/registries/api/bash/registry_loader.sh"
registry_resolve_paths

echo "[LAUNCH_group_cheats] REGISTRIES=${REGISTRIES}"
echo "[LAUNCH_group_cheats] PHASE8_DIR=${PHASE8_DIR}"
echo "[LAUNCH_group_cheats] chroms=${CHROMS_CSV:-ALL}"

# ── Get candidate list (tier ≤ 3) ────────────────────────────────────────────
CIDS=($(registry_list_candidates_by_tier 3))
echo "[LAUNCH_group_cheats] ${#CIDS[@]} candidates at tier<=3"

# ── Per-candidate dispatch with validation gates ─────────────────────────────
for cid in "${CIDS[@]}"; do
  level=$(registry_get_validation "$cid")
  echo "[LAUNCH_group_cheats] ${cid}: validation=${level}"
  RAW_OUT="${EVIDENCE_REGISTRY}/per_candidate/${cid}/raw/"
  mkdir -p "${RAW_OUT}"

  # ── cheat30 GDS — requires NONE (always runs) ────────────────────────────
  echo "  cheat30 GDS → running"
  Rscript "${PHASE8_DIR}/q5_age_and_origin/cheat30_gds_by_genotype.R" \
    --candidate "$cid" \
    --outdir "${RAW_OUT}" \
    || echo "    cheat30 failed for $cid" >&2

  # ── sd_substrate: Angle A + Angle B + concordance (pass 19) ──────────────
  # Gate: NONE (junction/SD sequence is group-independent).
  # Three independent scripts, three independent registry blocks:
  #   Angle A — sd_substrate_minimap2.R (trusted NAHR call)
  #   Angle B — sd_substrate_biser2.R   (catalog cross-check)
  #   Derived — sd_substrate_concordance.R (joint verdict)
  #
  # Reads breakpoints from the fragment_distribution block written by
  # STEP_03B (bp_pipeline_bridge). If that block isn't present yet,
  # skips with a warning.
  #
  # Env vars:
  #   REF_FASTA       → needed for Angle A
  #   BISER2_TSV      → needed for Angle B
  #   NO_MINIMAP2=1   → skip Angle A (fast BISER2-only sweep)
  #   SD_MM2_PRESET   → override asm10 default (asm5/asm10/asm20)
  #   SD_OVERWRITE=1  → force re-run even if BAM/TSV exist
  echo "  sd_substrate (minimap2 + BISER2 + concordance) → running"
  COORDS_LEFT=$(registry_get_key "$cid" "q3_final_left_bp"  2>/dev/null || echo "")
  COORDS_RIGHT=$(registry_get_key "$cid" "q3_final_right_bp" 2>/dev/null || echo "")
  COORDS_CHR=$(registry_get_key "$cid" "q3_chrom" 2>/dev/null || echo "")
  if [[ -z "$COORDS_LEFT" || -z "$COORDS_RIGHT" || -z "$COORDS_CHR" ]]; then
      echo "    sd_substrate → SKIP (no breakpoint coords in registry for $cid)"
  else
      SD_DIR="${PHASE8_DIR}/q4_mechanism/sd_substrate"
      RAN_ANY=0

      # Angle A — minimap2
      if [[ -z "${NO_MINIMAP2:-}" && -n "${REF_FASTA:-}" ]]; then
          Rscript "${SD_DIR}/sd_substrate_minimap2.R" \
            --candidate "$cid" \
            --chrom "$COORDS_CHR" \
            --left_bp "$COORDS_LEFT" \
            --right_bp "$COORDS_RIGHT" \
            --ref_fasta "$REF_FASTA" \
            ${SD_MM2_PRESET:+--mm2_preset "$SD_MM2_PRESET"} \
            ${SD_OVERWRITE:+--overwrite} \
            ${REGISTRIES:+--registries_root "$REGISTRIES"} \
            --outdir "${RAW_OUT}" \
            && RAN_ANY=1 \
            || echo "    sd_substrate Angle A failed for $cid" >&2
      else
          echo "    sd_substrate Angle A → SKIP (NO_MINIMAP2=1 or no REF_FASTA)"
      fi

      # Angle B — BISER2
      if [[ -n "${BISER2_TSV:-}" ]]; then
          Rscript "${SD_DIR}/sd_substrate_biser2.R" \
            --candidate "$cid" \
            --chrom "$COORDS_CHR" \
            --left_bp "$COORDS_LEFT" \
            --right_bp "$COORDS_RIGHT" \
            --biser2_tsv "$BISER2_TSV" \
            ${SD_OVERWRITE:+--overwrite} \
            ${REGISTRIES:+--registries_root "$REGISTRIES"} \
            --outdir "${RAW_OUT}" \
            && RAN_ANY=1 \
            || echo "    sd_substrate Angle B failed for $cid" >&2
      else
          echo "    sd_substrate Angle B → SKIP (no BISER2_TSV)"
      fi

      # Concordance — only if at least one angle produced a block
      if [[ "$RAN_ANY" -eq 1 ]]; then
          Rscript "${SD_DIR}/sd_substrate_concordance.R" \
            --candidate "$cid" \
            ${REGISTRIES:+--registries_root "$REGISTRIES"} \
            --outdir "${RAW_OUT}" \
            || echo "    sd_substrate concordance failed for $cid" >&2
      fi
  fi

  # ── STEP_03 per-read evidence ledger (pass 21) ──────────────────────────
  # Extracts per-read FF/RR discordants + split reads + soft clips from
  # each sample's BAM, writes:
  #   evidence_reads_bam.tsv.gz    — row-level ledger (29 cols, read names + tags)
  #   evidence_reads_vcf.tsv.gz    — sidecar with DELLY/Manta aggregate counts
  #   per_read_evidence_summary.json — registry block with q7b_bam_* counts
  #
  # Required env vars:
  #   BAM_DIR           directory of <sample>.markdup.bam files
  #   SAMPLES_GROUPS_TSV  TSV (sample_id, group) from Q1 decomposition
  # Optional:
  #   DELLY_INV_VCF     for VCF sidecar
  #   MANTA_INV_VCF     for VCF sidecar
  #   PER_READ_WINDOW   default 2000
  #   PER_READ_MIN_MAPQ default 20
  #   PER_READ_THREADS  default 4
  #   SD_OVERWRITE=1    re-run even if TSVs exist
  if [[ -n "${BAM_DIR:-}" && -n "${SAMPLES_GROUPS_TSV:-}" ]]; then
      echo "  STEP_03 per-read evidence → running"
      python3 "${PHASE8_DIR}/q7_existence_audit/STEP_03_per_read_evidence.py" \
        --candidate "$cid" \
        ${REGISTRIES:+--registries_root "$REGISTRIES"} \
        --bam_dir "$BAM_DIR" \
        --samples_tsv "$SAMPLES_GROUPS_TSV" \
        ${DELLY_INV_VCF:+--delly_inv_vcf "$DELLY_INV_VCF"} \
        ${MANTA_INV_VCF:+--manta_inv_vcf "$MANTA_INV_VCF"} \
        ${PER_READ_WINDOW:+--window "$PER_READ_WINDOW"} \
        ${PER_READ_MIN_MAPQ:+--min_mapq "$PER_READ_MIN_MAPQ"} \
        ${PER_READ_THREADS:+--threads "$PER_READ_THREADS"} \
        ${SD_OVERWRITE:+--overwrite} \
        --outdir "${RAW_OUT}" \
        || echo "    STEP_03 per-read evidence failed for $cid" >&2
  else
      echo "  STEP_03 per-read evidence → SKIP (set BAM_DIR + SAMPLES_GROUPS_TSV to enable)"
  fi

  # ── cheat28 tandem-repeat context — requires NONE ────────────────────────
  echo "  cheat28 tandem-repeat → running"
  Rscript "${PHASE8_DIR}/q4_mechanism/cheat28_tandem_repeat_context.R" \
    --candidate "$cid" \
    --outdir "${RAW_OUT}" \
    || echo "    cheat28 failed for $cid" >&2

  # ── cheat29 junction forensics — requires NONE ───────────────────────────
  echo "  cheat29 junction forensics → running"
  Rscript "${PHASE8_DIR}/q4_mechanism/cheat29_junction_forensics.R" \
    --candidate "$cid" \
    --outdir "${RAW_OUT}" \
    || echo "    cheat29 failed for $cid" >&2

  # ── cheat6 ancestry jackknife — requires UNCERTAIN ───────────────────────
  # Moved 2026-04-24 (pass 17): q7_existence_audit/ → q6_group_robustness/.
  # This is a group-trust diagnostic (robustness to family drop-out), not
  # an existence test. C01f in phase_7/validation has the production T9
  # (self-contained); this standalone exists for per-candidate manuscript
  # plots and audit.
  if registry_check_validation "$cid" UNCERTAIN; then
    echo "  cheat6 jackknife → running"
    Rscript "${PHASE8_DIR}/q6_group_robustness/cheat6_ancestry_jackknife_v934.R" \
      --candidate "$cid" \
      --outdir "${RAW_OUT}" \
      || echo "    cheat6 failed for $cid" >&2
  else
    echo "  cheat6 jackknife → SKIP (validation below UNCERTAIN)"
  fi

  # ── Fst sub-block of age_evidence — requires SUPPORTED ───────────────────
  # ASPIRATIONAL: compute_age_fst_subblock.R not yet in tarball (chat-9
  # Finding W). When it lands, uncomment the call below.
  if registry_check_validation "$cid" SUPPORTED; then
    echo "  age_evidence Fst sub-block → PENDING SCRIPT (SUPPORTED reached)"
    # Rscript "${MODULE_DIR}/scripts/compute_age_fst_subblock.R" \
    #   --candidate "$cid" \
    #   --outdir "${RAW_OUT}" \
    #   || echo "    age_fst failed for $cid" >&2
  else
    echo "  age_evidence Fst sub-block → SKIP (below SUPPORTED)"
  fi

  # ── Burden regression — requires VALIDATED ───────────────────────────────
  # ASPIRATIONAL: STEP_C01f_c_burden_regression.R not yet in tarball.
  if registry_check_validation "$cid" VALIDATED; then
    echo "  burden regression → PENDING SCRIPT (VALIDATED reached)"
    # Rscript "${MODULE_DIR}/burden/STEP_C01f_c_burden_regression.R" \
    #   --candidate "$cid" \
    #   --outdir "${RAW_OUT}" \
    #   || echo "    burden failed for $cid" >&2
  else
    echo "  burden regression → SKIP (below VALIDATED)"
  fi
done

echo "[LAUNCH_group_cheats] done"
