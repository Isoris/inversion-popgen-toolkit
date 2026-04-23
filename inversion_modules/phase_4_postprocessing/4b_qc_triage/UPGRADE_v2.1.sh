#!/bin/bash
# =============================================================================
# UPGRADE_v2.1.sh — phase_qc_shelf v2.0 -> v2.1 incremental upgrade
# =============================================================================
# Adds the dropout-semantics features on top of a working v2.0 install:
#   - No-data background strips in Q04 (decorate())
#   - Two-variant Q04 output (with / without strips)
#   - Two-layer ideogram stipple (zero-SNP + reference-N)
#   - New STEP_Q09_gap_characterization + R/q09_gap_characterization.R
#   - New scripts/build_ref_n_bed.sh
#   - New run_chrom.sh per-chromosome wrapper
#   - STEP_Q04 driver: OUT_SUFFIX + Q04_NO_STRIPS + REF_N_BED env support
#
# Idempotent: re-runnable. Backs everything up with .preUPGRADE_21_<ts>.
#
# Usage:
#   bash UPGRADE_v2.1.sh
#   bash UPGRADE_v2.1.sh --dry-run
# =============================================================================
set -euo pipefail

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

MOD_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${MOD_ROOT}"
[[ -f 00_config.sh ]] || { echo "ERROR: run from phase_qc_shelf/"; exit 1; }

log()  { printf "[UPGRADE21] %s\n" "$*"; }
apply() { if [[ "${DRY_RUN}" == "1" ]]; then log "DRY: would $*"; else log "$*"; eval "$@"; fi; }

ts="$(date +%Y%m%d_%H%M%S)"
backup() {
  local f="$1"
  [[ -f "${f}" ]] || return 0
  local bk="${f}.preUPGRADE_21_${ts}"
  [[ "${DRY_RUN}" == "1" ]] || cp "${f}" "${bk}"
  log "  backup: ${bk}"
}

# -----------------------------------------------------------------------------
# 1. R/q04_compose_plot.R : ensure decorate() uses no-data strips
# -----------------------------------------------------------------------------
log "== 1/5: Q04 decorate() with no-data strips + stipple ideogram =="
if grep -q "detect_nodata_regions" R/q04_compose_plot.R && \
   grep -q "make_stipple" R/q04_compose_plot.R && \
   grep -q "NO_STRIPS" R/q04_compose_plot.R; then
  log "  already patched"
else
  backup R/q04_compose_plot.R
  if [[ "${DRY_RUN}" == "1" ]]; then
    log "  DRY: would replace decorate block and ideogram block"
  else
    python3 <<'PY'
with open("R/q04_compose_plot.R") as f:
    s = f.read()

# Patch A: add NO_STRIPS + REF_N_BED after SMOOTH_WIN parse
if "NO_STRIPS" not in s:
    old_a = 'SMOOTH_WIN <- as.integer(get_arg("--smooth_win", "1"))'
    new_a = '''SMOOTH_WIN <- as.integer(get_arg("--smooth_win", "1"))
# --no_nodata_strips disables the gray "no data here" background bands.
NO_STRIPS <- !is.null(get_arg("--no_nodata_strips", NULL))
# --ref_n_bed: optional BED of reference-N regions for the current chromosome
REF_N_BED <- get_arg("--ref_n_bed", NULL)'''
    s = s.replace(old_a, new_a)
    print("A applied")

# Patch B: ensure span_kb is computed at load time
if 'dt[, span_kb := (end_bp - start_bp) / 1e3]' not in s.split("# =")[0]:
    # Find dt <- as.data.table(pc$dt) block and add span_kb right after
    old_b = 'dt <- as.data.table(pc$dt)\ndt[, mb := (start_bp + end_bp) / 2 / 1e6]\nz_col <- NULL'
    new_b = 'dt <- as.data.table(pc$dt)\ndt[, mb := (start_bp + end_bp) / 2 / 1e6]\ndt[, span_kb := (end_bp - start_bp) / 1e3]\nz_col <- NULL'
    if old_b in s:
        s = s.replace(old_b, new_b)
        print("B applied")

# Patch C: decorate() with no-data strips + honoring NO_STRIPS
# If detect_nodata_regions is absent, insert before the current decorate
if "detect_nodata_regions" not in s:
    import re
    decorate_new = '''# =============================================================================
# No-data regions + ideogram stipple helpers
# =============================================================================
.NODATA_CACHE <- NULL
detect_nodata_regions <- function() {
  if (!is.null(.NODATA_CACHE)) return(.NODATA_CACHE)
  if (!exists("dt", envir = .GlobalEnv)) return(data.table())
  w <- copy(dt)
  w[, is_gap := FALSE]
  if ("n_snps"  %in% names(w)) w[n_snps == 0, is_gap := TRUE]
  if ("span_kb" %in% names(w)) {
    span_hi <- stats::quantile(w$span_kb, 0.995, na.rm = TRUE)
    w[span_kb > max(span_hi, 80), is_gap := TRUE]
  }
  if (!any(w$is_gap)) { .NODATA_CACHE <<- data.table(); return(.NODATA_CACHE) }
  g <- w[is_gap == TRUE][order(start_bp)]
  g[, run_grp := cumsum(c(TRUE, diff(start_bp) > 50e3))]
  out <- g[, .(xmin = min(start_bp) / 1e6, xmax = max(end_bp) / 1e6),
           by = run_grp]
  out[, run_grp := NULL]
  .NODATA_CACHE <<- out
  out
}

make_stipple <- function(xmin, xmax, y_lo = 0.05, y_hi = 0.95,
                         dx = 0.05, dy = 0.15, offset = 0) {
  if (!length(xmin) || !length(xmax)) return(data.frame(x = numeric(0), y = numeric(0)))
  dots <- list()
  for (i in seq_along(xmin)) {
    xs <- seq(xmin[i] + dx / 2 + offset * dx / 2, xmax[i] - dx / 2, by = dx)
    ys <- seq(y_lo + dy / 2, y_hi - dy / 2, by = dy)
    if (!length(xs) || !length(ys)) next
    grid_rows <- list()
    for (j in seq_along(ys)) {
      x_row <- if (j %% 2 == 1) xs else xs + dx / 2
      x_row <- x_row[x_row <= xmax[i] - dx / 4]
      if (length(x_row)) grid_rows[[j]] <- data.frame(x = x_row, y = ys[j])
    }
    if (length(grid_rows)) dots[[i]] <- do.call(rbind, grid_rows)
  }
  if (!length(dots)) return(data.frame(x = numeric(0), y = numeric(0)))
  do.call(rbind, dots)
}

load_ref_n_regions <- function(bed_path, chrom) {
  if (is.null(bed_path) || !file.exists(bed_path)) return(data.table(xmin = numeric(0), xmax = numeric(0)))
  b <- tryCatch(fread(bed_path, header = FALSE, col.names = c("chrom", "start", "end")),
                error = function(e) NULL)
  if (is.null(b) || !nrow(b)) return(data.table(xmin = numeric(0), xmax = numeric(0)))
  b <- b[chrom == CHROM][, .(xmin = start / 1e6, xmax = end / 1e6)]
  b
}

decorate <- function(p) {
  nd <- if (isTRUE(NO_STRIPS)) data.table() else detect_nodata_regions()
  if (nrow(nd) > 0) {
    p <- p + geom_rect(data = nd,
                       aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       inherit.aes = FALSE,
                       fill = "#888888", alpha = 0.18)
  }
  if (is.finite(SHELF_A) && is.finite(SHELF_B)) {
    p <- p + annotate("rect", xmin = SHELF_A, xmax = SHELF_B,
                      ymin = -Inf, ymax = Inf,
                      fill = "#f5a524", alpha = 0.10)
  }
  if (length(BPS) > 0) {
    p <- p + geom_vline(xintercept = BPS, linetype = "dashed",
                        color = "#e0555c", linewidth = 0.35, alpha = 0.8)
  }
  p
}
'''
    # Replace existing decorate block
    pat = re.compile(r'# Shelf shading.*?^decorate <- function\(p\) \{.*?^\}\n', re.S | re.M)
    s2, n = pat.subn(decorate_new, s, count=1)
    if n == 1:
        s = s2
        print("C applied (replaced decorate)")
    else:
        print("C FAIL — decorate pattern not found")

with open("R/q04_compose_plot.R", "w") as f:
    f.write(s)
PY
    log "  q04 patches applied"
  fi
fi

# -----------------------------------------------------------------------------
# 2. STEP_Q04 driver: OUT_SUFFIX + Q04_NO_STRIPS + REF_N_BED
# -----------------------------------------------------------------------------
log "== 2/5: STEP_Q04 driver env-var support =="
if grep -q "OUT_SUFFIX" STEP_Q04_plot_diagnostic.sh && \
   grep -q "Q04_NO_STRIPS" STEP_Q04_plot_diagnostic.sh; then
  log "  already patched"
else
  backup STEP_Q04_plot_diagnostic.sh
  if [[ "${DRY_RUN}" == "1" ]]; then
    log "  DRY: would add OUT_SUFFIX / Q04_NO_STRIPS / REF_N_BED handling"
  else
    python3 <<'PY'
with open("STEP_Q04_plot_diagnostic.sh") as f:
    s = f.read()
# Replace the out_pdf line
s = s.replace(
    'local out_pdf="${QC_FIGS}/diagnostic.${chr}.pdf"',
    'local out_pdf="${QC_FIGS}/diagnostic.${chr}${OUT_SUFFIX:-}.pdf"'
)
# Add --no_nodata_strips and --snp_density_scale_kb and --ref_n_bed before the Rscript call
if "Q04_NO_STRIPS" not in s:
    marker = '  [[ -n "${SMOOTH_WIN:-}"   ]] && args+=( --smooth_win     "${SMOOTH_WIN}"     )'
    extra = marker + '\n' + \
            '  [[ -n "${SNP_DENSITY_SCALE_KB:-}" ]] && args+=( --snp_density_scale_kb "${SNP_DENSITY_SCALE_KB}" )\n' + \
            '  [[ "${Q04_NO_STRIPS:-}" == "1" ]] && args+=( --no_nodata_strips yes )\n' + \
            '  # Reference-N BED (for dark-stipple assembly-gap layer on the ideogram)\n' + \
            '  local ref_n_bed=""\n' + \
            '  if [[ -f "${QC_TRACKS}/ref_n.${chr}.bed" ]]; then\n' + \
            '    ref_n_bed="${QC_TRACKS}/ref_n.${chr}.bed"\n' + \
            '  elif [[ -n "${REF_N_BED:-}" && -f "${REF_N_BED}" ]]; then\n' + \
            '    ref_n_bed="${REF_N_BED}"\n' + \
            '  fi\n' + \
            '  [[ -n "${ref_n_bed}" ]] && args+=( --ref_n_bed "${ref_n_bed}" )'
    s = s.replace(marker, extra)

# Log file and OUT_SUFFIX
s = s.replace(
    '2> "${QC_LOGS}/q04_${chr}.log"',
    '2> "${QC_LOGS}/q04_${chr}${OUT_SUFFIX:-}.log"'
)

with open("STEP_Q04_plot_diagnostic.sh", "w") as f:
    f.write(s)
PY
    log "  Q04 driver patched"
  fi
fi

# -----------------------------------------------------------------------------
# 3. New files from v2.1 bundle
# -----------------------------------------------------------------------------
log "== 3/5: new files (Q09, build_ref_n_bed, run_chrom) =="
for new in STEP_Q09_gap_characterization.sh R/q09_gap_characterization.R \
           scripts/build_ref_n_bed.sh run_chrom.sh; do
  if [[ -f "${new}" ]]; then
    log "  ${new} already present"
  else
    log "  NOTE: ${new} missing — extract from phase_qc_shelf_v2.1.tar.gz"
  fi
done

# -----------------------------------------------------------------------------
# 4. Make new shell scripts executable
# -----------------------------------------------------------------------------
log "== 4/5: chmod +x =="
for f in STEP_Q09_gap_characterization.sh run_chrom.sh \
         scripts/build_ref_n_bed.sh; do
  [[ -f "${f}" ]] && apply "chmod +x ${f}"
done

# -----------------------------------------------------------------------------
# 5. Verify
# -----------------------------------------------------------------------------
log "== 5/5: verify =="
for f in $(find . -name "*.sh"); do
  bash -n "${f}" && log "  bash OK  ${f}" || log "  bash FAIL ${f}"
done

log "=== UPGRADE 2.1 complete ==="
log "Run a test render: "
log "  SHELF_START_MB=15 SHELF_END_MB=18 BP1_MB=15.115 BP2_MB=18.005 \\"
log "    bash run_chrom.sh C_gar_LG28"
