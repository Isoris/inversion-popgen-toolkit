#!/bin/bash
# =============================================================================
# UPGRADE.sh — phase_qc_shelf v2.0 upgrade
# =============================================================================
# Idempotent in-place upgrade of an existing phase_qc_shelf install to the
# consolidated v2.0 state. Re-runnable safely.
#
# What it does:
#   - Adds config.local.sh if missing (anchors INSTANT_Q_BIN, LOCAL_Q_DIR,
#     BEAGLE_DIR to the correct repo-local paths)
#   - Patches install.sh set-e bug (missing=$((missing+1)) in && form)
#   - Patches STEP_Q06_precompute.sh to preserve INSTANT_Q_BIN / LOCAL_Q_DIR
#     when 00_ancestry_config.sh is sourced (otherwise it clobbers them)
#   - Replaces R/q06_ancestry_tracks.R with the Engine-B-native version
#     (reads mean_delta12/mean_entropy/mean_ena directly; handles sample_id,
#     assigned_pop, max_q)
#   - Patches STEP_Q03_coverage_tracks.sh to soft-skip missing mosdepth dir
#   - Patches R/q04_compose_plot.R:
#       * Z-column candidate list (max_abs_z, MDS_z fallbacks)
#       * Correct arg name for ancestry_multi / popstats_invgt
#       * NULL-out placeholder panels so missing data doesn't draw blanks
#       * SNP density panel with configurable scale (--snp_density_scale_kb)
#       * n_sites_used filter for popstats panels (drops fake Fst zeros)
#   - Patches R/q08_shelf_heatmap.R subscript-out-of-bounds bug
#     (named-integer lookup -> environment hash)
#
# Usage:
#   bash UPGRADE.sh
#   bash UPGRADE.sh --dry-run          # show what would change
#
# Location: run from the phase_qc_shelf module root (where 00_config.sh lives).
# =============================================================================
set -euo pipefail

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

MOD_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${MOD_ROOT}"

[[ -f 00_config.sh ]] || { echo "ERROR: 00_config.sh not found. Run from phase_qc_shelf/" >&2; exit 1; }

log()  { printf "[UPGRADE] %s\n" "$*"; }
apply() { [[ "${DRY_RUN}" == "1" ]] && log "DRY: would $*" || { log "$*"; eval "$@"; }; }

ts="$(date +%Y%m%d_%H%M%S)"
backup() {
  local f="$1"
  [[ -f "${f}" ]] || return 0
  local bk="${f}.preUPGRADE_${ts}"
  [[ "${DRY_RUN}" == "1" ]] || cp "${f}" "${bk}"
  log "  backup: ${bk}"
}

# -----------------------------------------------------------------------------
# 1. config.local.sh
# -----------------------------------------------------------------------------
log "== Fix 1/7: config.local.sh =="
if [[ ! -f config.local.sh ]]; then
  log "  Creating config.local.sh"
  if [[ "${DRY_RUN}" != "1" ]]; then
    cat > config.local.sh <<'EOF'
# phase_qc_shelf local config — anchors all paths.
# Source order: 00_config.sh -> config.local.sh -> (optionally) 00_ancestry_config.sh
export BASE="${BASE:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"
REPO_ROOT="${BASE}/inversion-popgen-toolkit"
source "${REPO_ROOT}/inversion_modules/00_inversion_config.sh" 2>/dev/null || true

# Layer our defaults on top of the shared inversion config
export HET_DIR="${HET_DIR:-${HETDIR:-${BASE}/het_roh}}"
export SAMPLE_LIST="${SAMPLE_LIST:-${SAMPLES_IND:-${HET_DIR}/01_inputs_check/samples.ind}}"
export SAMPLE_LIST_POPSTATS="${SAMPLE_LIST_POPSTATS:-${SAMPLE_LIST}}"
export BEAGLE_DIR="${BASE}/inversion_localpca_v7/02_snps_beagle"
export UNIFIED_ANCESTRY_DIR="${REPO_ROOT}/unified_ancestry"
export LOCAL_Q_DIR="${UNIFIED_ANCESTRY_DIR}/local_Q"
export REGION_POPSTATS_BIN="${UNIFIED_ANCESTRY_DIR}/engines/fst_dxy/region_popstats"
export INSTANT_Q_BIN="${UNIFIED_ANCESTRY_DIR}/src/instant_q"
EOF
  fi
else
  # Add anchors if missing (idempotent)
  for var in INSTANT_Q_BIN BEAGLE_DIR UNIFIED_ANCESTRY_DIR LOCAL_Q_DIR; do
    if ! grep -q "^export ${var}=" config.local.sh; then
      case "${var}" in
        INSTANT_Q_BIN)        line='export INSTANT_Q_BIN="${UNIFIED_ANCESTRY_DIR}/src/instant_q"' ;;
        BEAGLE_DIR)           line='export BEAGLE_DIR="${BASE}/inversion_localpca_v7/02_snps_beagle"' ;;
        UNIFIED_ANCESTRY_DIR) line='export UNIFIED_ANCESTRY_DIR="${REPO_ROOT}/unified_ancestry"' ;;
        LOCAL_Q_DIR)          line='export LOCAL_Q_DIR="${UNIFIED_ANCESTRY_DIR}/local_Q"' ;;
      esac
      [[ "${DRY_RUN}" != "1" ]] && echo "${line}" >> config.local.sh
      log "  appended: ${line}"
    fi
  done
fi

# -----------------------------------------------------------------------------
# 2. install.sh set-e bug
# -----------------------------------------------------------------------------
log "== Fix 2/7: install.sh set-e missing=… =="
if grep -q '\[\[ "\${hard}" == "1" \]\] && missing=\$((missing + 1))' install.sh 2>/dev/null; then
  backup install.sh
  [[ "${DRY_RUN}" != "1" ]] && \
    sed -i 's|\[\[ "${hard}" == "1" \]\] && missing=$((missing + 1))|if [[ "${hard}" == "1" ]]; then missing=$((missing + 1)); fi|' install.sh
  log "  patched set-e conditional"
else
  log "  already patched (or not present in this version)"
fi

# -----------------------------------------------------------------------------
# 3. STEP_Q06_precompute.sh — preserve INSTANT_Q_BIN / LOCAL_Q_DIR
# -----------------------------------------------------------------------------
log "== Fix 3/7: STEP_Q06_precompute.sh preserve-and-restore =="
if [[ -f STEP_Q06_precompute.sh ]] && grep -qF 'source "${ANC_CFG}"' STEP_Q06_precompute.sh && \
   ! grep -q '_preserve_bin' STEP_Q06_precompute.sh; then
  backup STEP_Q06_precompute.sh
  if [[ "${DRY_RUN}" != "1" ]]; then
    python3 - <<'PY'
import pathlib, re
p = pathlib.Path("STEP_Q06_precompute.sh")
s = p.read_text()
old = '''ANC_CFG="${UNIFIED_ANCESTRY_DIR}/00_ancestry_config.sh"
if [[ -f "${ANC_CFG}" ]]; then
  # shellcheck disable=SC1090
  source "${ANC_CFG}"
fi'''
new = '''ANC_CFG="${UNIFIED_ANCESTRY_DIR}/00_ancestry_config.sh"
if [[ -f "${ANC_CFG}" ]]; then
  _preserve_bin="${INSTANT_Q_BIN:-}"
  _preserve_local_q="${LOCAL_Q_DIR:-}"
  # shellcheck disable=SC1090
  source "${ANC_CFG}"
  [[ -n "${_preserve_bin}"      ]] && INSTANT_Q_BIN="${_preserve_bin}"
  [[ -n "${_preserve_local_q}"  ]] && LOCAL_Q_DIR="${_preserve_local_q}"
  unset _preserve_bin _preserve_local_q
fi'''
if old in s:
    p.write_text(s.replace(old, new))
    print("  PATCHED STEP_Q06_precompute.sh")
else:
    print("  already patched or pattern drifted")
PY
  fi
else
  log "  already patched or file not present"
fi

# -----------------------------------------------------------------------------
# 4. STEP_Q03_coverage_tracks.sh — soft-skip missing mosdepth dir
# -----------------------------------------------------------------------------
log "== Fix 4/7: STEP_Q03 soft-skip =="
if grep -q 'qc_die "No mosdepth dir' STEP_Q03_coverage_tracks.sh 2>/dev/null; then
  backup STEP_Q03_coverage_tracks.sh
  [[ "${DRY_RUN}" != "1" ]] && \
    sed -i 's|\[\[ -d "${search_dir}" \]\] || qc_die "No mosdepth dir: ${search_dir}"|if [[ ! -d "${search_dir}" ]]; then qc_log "SKIP ${chr}: no mosdepth dir ${search_dir}"; return 0; fi|' STEP_Q03_coverage_tracks.sh
  log "  patched Q03 to skip gracefully"
else
  log "  already patched"
fi

# -----------------------------------------------------------------------------
# 5. R/q06_ancestry_tracks.R — replace with Engine B-native reader
# -----------------------------------------------------------------------------
log "== Fix 5/7: R/q06_ancestry_tracks.R (Engine B native reader) =="
# Mark: the new version has "Engine B native columns" in header
if ! grep -q 'Engine B native columns' R/q06_ancestry_tracks.R 2>/dev/null; then
  backup R/q06_ancestry_tracks.R
  if [[ "${DRY_RUN}" != "1" ]]; then
    cat > R/q06_ancestry_tracks.R <<'RSCRIPT_EOF'
#!/usr/bin/env Rscript
# =============================================================================
# q06_ancestry_tracks.R  (Engine B native columns, K-agnostic)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag); if (!length(i)) return(default); args[i + 1]
}
SUMMARY  <- get_arg("--summary")
SAMPLES  <- get_arg("--samples", NULL)
CHROM    <- get_arg("--chrom")
OUT_WIN  <- get_arg("--out_win")
OUT_SAMP <- get_arg("--out_samp")
stopifnot(!is.null(SUMMARY), !is.null(CHROM), !is.null(OUT_WIN))

message("[q06] Loading window summary: ", SUMMARY)
s <- fread(SUMMARY)
setnames(s, tolower(names(s)))
if ("chrom" %in% names(s)) s <- s[chrom == CHROM]

required <- c("start_bp", "end_bp", "mean_delta12", "mean_entropy", "mean_ena")
missing  <- setdiff(required, names(s))
if (length(missing) > 0) {
  stop("Summary missing required columns: ", paste(missing, collapse = ", "),
       "\nAvailable: ", paste(names(s), collapse = ", "))
}
s[, window_mid_bp := as.integer((start_bp + end_bp) / 2)]

sample_dt <- NULL
if (!is.null(SAMPLES) && file.exists(SAMPLES)) {
  message("[q06] Loading per-sample cache: ", SAMPLES)
  sample_dt <- fread(SAMPLES)
  setnames(sample_dt, tolower(names(sample_dt)))
  if ("chrom" %in% names(sample_dt)) sample_dt <- sample_dt[chrom == CHROM]
  sample_col_aliases <- c("sample", "sample_id", "sampleid", "ind", "individual", "id")
  found_sample <- intersect(sample_col_aliases, names(sample_dt))
  if (length(found_sample) > 0 && found_sample[1] != "sample") {
    setnames(sample_dt, found_sample[1], "sample")
  }
  if (!"maxq_label" %in% names(sample_dt)) {
    if ("assigned_pop" %in% names(sample_dt)) {
      sample_dt[, maxq_label := paste0("K", assigned_pop)]
    } else {
      qcols <- grep("^q[0-9]+$", names(sample_dt), value = TRUE)
      if (length(qcols) >= 2) {
        qm <- as.matrix(sample_dt[, ..qcols])
        sample_dt[, maxq_label := paste0("K", max.col(qm, ties.method = "first"))]
      }
    }
  }
}

if (!is.null(sample_dt) && all(c("window_id", "delta12") %in% names(sample_dt))) {
  cv_dt <- sample_dt[, .(
    cv_delta12_across_samples = {
      x <- delta12; m <- mean(x, na.rm = TRUE)
      if (is.na(m) || m == 0) NA_real_ else stats::sd(x, na.rm = TRUE) / m
    }
  ), by = window_id]
  s <- merge(s, cv_dt, by = "window_id", all.x = TRUE, sort = FALSE)
}
if (!"cv_delta12_across_samples" %in% names(s)) s[, cv_delta12_across_samples := NA_real_]

if (!is.null(sample_dt) && all(c("window_id", "maxq_label") %in% names(sample_dt))) {
  dom <- sample_dt[, .(
    maxQ_label = {
      tt <- table(maxq_label, useNA = "no")
      if (length(tt) == 0) NA_character_ else names(tt)[which.max(tt)]
    }
  ), by = window_id]
  s <- merge(s, dom, by = "window_id", all.x = TRUE, sort = FALSE)
}
if (!"maxQ_label" %in% names(s)) s[, maxQ_label := NA_character_]

out_win_dt <- data.table(
  chrom                     = CHROM,
  window_id                 = s$window_id,
  window_start_bp           = as.integer(s$start_bp),
  window_end_bp             = as.integer(s$end_bp),
  window_mid_mb             = round(s$window_mid_bp / 1e6, 4),
  delta12                   = round(s$mean_delta12, 5),
  entropy                   = round(s$mean_entropy, 5),
  ena                       = round(s$mean_ena, 4),
  maxQ_label                = s$maxQ_label,
  cv_delta12_across_samples = round(s$cv_delta12_across_samples, 4)
)
setorder(out_win_dt, window_start_bp)
fwrite(out_win_dt, OUT_WIN, sep = "\t", quote = FALSE)
message("[q06] Wrote ", nrow(out_win_dt), " window rows -> ", OUT_WIN)

if (!is.null(sample_dt) && "sample" %in% names(sample_dt) &&
    nrow(sample_dt) > 0 && !is.null(OUT_SAMP)) {
  keep <- intersect(
    c("sample", "chrom", "window_id", "max_q", "maxq", "maxq_label",
      "assigned_pop", "secondq", "delta12", "delta13",
      "entropy", "effective_num_ancestries", "ena"),
    names(sample_dt))
  out_samp <- sample_dt[, ..keep]
  rn <- names(out_samp)
  rn <- sub("^max_q$",      "maxQ",       rn)
  rn <- sub("^maxq$",       "maxQ",       rn)
  rn <- sub("^maxq_label$", "maxQ_label", rn)
  rn <- sub("^secondq$",    "secondQ",    rn)
  setnames(out_samp, names(out_samp), rn)
  mids <- s[, .(window_id, window_mid_bp)]
  out_samp <- merge(out_samp, mids, by = "window_id", all.x = TRUE, sort = FALSE)
  setorder(out_samp, sample, window_mid_bp)
  fwrite(out_samp, OUT_SAMP, sep = "\t", quote = FALSE, compress = "gzip")
  message("[q06] Wrote ", nrow(out_samp), " per-sample rows -> ", OUT_SAMP)

  if ("maxQ_label" %in% names(out_samp)) {
    wide <- dcast(out_samp, window_mid_bp ~ sample, value.var = "maxQ_label")
    setorder(wide, window_mid_bp)
    out_wide <- sub("\\.tsv\\.gz$", ".maxQ_wide.tsv.gz", OUT_SAMP)
    fwrite(wide, out_wide, sep = "\t", quote = FALSE, compress = "gzip")
    message("[q06] Wrote maxQ wide (", nrow(wide), "x", ncol(wide)-1, ") -> ", out_wide)
  }
}
message("[q06] DONE.")
RSCRIPT_EOF
  fi
  log "  replaced R/q06_ancestry_tracks.R with Engine-B-native version"
else
  log "  already up-to-date"
fi

# -----------------------------------------------------------------------------
# 6. R/q04_compose_plot.R — multiple patches
# -----------------------------------------------------------------------------
log "== Fix 6/7: R/q04_compose_plot.R =="
needs_patch=0

# 6a. Z column candidate list
if ! grep -q '"max_abs_z"' R/q04_compose_plot.R 2>/dev/null; then needs_patch=1; fi
# 6b. placeholders -> NULL (sign: look for the old placeholder pattern)
if grep -q 'p3 <- placeholder' R/q04_compose_plot.R 2>/dev/null; then needs_patch=1; fi
# 6c. snp density configurable
if ! grep -q 'snp_density_scale_kb' R/q04_compose_plot.R 2>/dev/null; then needs_patch=1; fi
# 6d. n_sites_used filter
if ! grep -q 'Q04_MIN_SITES' R/q04_compose_plot.R 2>/dev/null; then needs_patch=1; fi

if [[ ${needs_patch} -eq 1 ]]; then
  backup R/q04_compose_plot.R
  [[ "${DRY_RUN}" != "1" ]] && python3 - <<'PY'
import pathlib, re
p = pathlib.Path("R/q04_compose_plot.R")
s = p.read_text()

# --- 6a: Z candidate list
s = re.sub(
    r'for \(cand in c\("robust_z", "z_robust", "z", "mds_z_robust", "mds_z1_robust"\)\)',
    'for (cand in c("max_abs_z", "robust_z", "z_robust", "z", "mds_z_robust", "mds_z1_robust", "MDS1_z"))',
    s, count=1)

# --- 6b: NULL-ify placeholder panels
for var in ("p3", "p4", "p5", "p6", "p7"):
    s = re.sub(rf'^{var} <- placeholder\(.*\)\s*$', f'{var} <- NULL', s, flags=re.M)

# Replace the compose block with NULL-filtering
old_compose = '''panels  <- list(ideogram, sim_strip, p1, p2, p3, p4, p5, p6, p7)
heights <- c(0.5,          0.5,       1.3, 0.9, 0.9, 0.9, 0.6, 0.9, 1.2)
if (!is.null(p7b)) { panels <- c(panels, list(p7b)); heights <- c(heights, 0.9) }
if (!is.null(p8))  { panels <- c(panels, list(p8));  heights <- c(heights, 0.9) }
if (!is.null(p9))  { panels <- c(panels, list(p9));  heights <- c(heights, 0.9) }'''
new_compose = '''all_panels <- list(
  list(p = ideogram, h = 0.5),
  list(p = sim_strip, h = 0.5),
  list(p = p1,  h = 1.3),
  list(p = p2,  h = 0.9),
  list(p = p3,  h = 0.9),
  list(p = p4,  h = 0.9),
  list(p = p5,  h = 0.6),
  list(p = p6,  h = 0.9),
  list(p = p7,  h = 1.2),
  list(p = p7b, h = 0.9),
  list(p = p8,  h = 0.9),
  list(p = p9,  h = 0.9)
)
keep_idx <- which(sapply(all_panels, function(x) !is.null(x$p)))
panels  <- lapply(all_panels[keep_idx], `[[`, "p")
heights <- sapply(all_panels[keep_idx], `[[`, "h")'''
if old_compose in s:
    s = s.replace(old_compose, new_compose)

# --- 6c: Configurable SNP density panel
old_p2 = '''# Track 2: n_snps
# =============================================================================
p2 <- decorate(
  ggplot(dt, aes(mb, n_snps)) +
    geom_point(size = 0.22, color = "#2c7a39", alpha = 0.55) +
    geom_hline(yintercept = median(dt$n_snps, na.rm = TRUE),
               linetype = "dashed", color = "#2c3e50", linewidth = 0.25)
) + common_x + labs(y = "SNPs/win") + base_theme +
  edge_labels("dense", "sparse")'''
new_p2 = '''# Track 2: SNP density (per kb / 10kb / 50kb, configurable)
# =============================================================================
sd_scale_kb <- as.numeric(get_arg("--snp_density_scale_kb", "10"))
dt[, span_kb    := (end_bp - start_bp) / 1e3]
dt[, snps_per_x := n_snps * sd_scale_kb / pmax(span_kb, 1e-9)]
if (SMOOTH_WIN > 1) {
  dt[, snps_per_x_sm := rolling_median(snps_per_x, SMOOTH_WIN)]
} else {
  dt[, snps_per_x_sm := snps_per_x]
}
sd_unit <- if (sd_scale_kb == 1) "SNPs/kb" else sprintf("SNPs/%gkb", sd_scale_kb)
p2 <- decorate(
  ggplot(dt, aes(mb, snps_per_x_sm)) +
    geom_line(color = "#ffffff", linewidth = 1.2, alpha = 0.75) +
    geom_line(color = "#2c7a39", linewidth = 0.45) +
    geom_hline(yintercept = median(dt$snps_per_x, na.rm = TRUE),
               linetype = "dashed", color = "#2c3e50", linewidth = 0.25)
) + common_x + labs(y = sd_unit) + base_theme +
  edge_labels("dense", "sparse")'''
if old_p2 in s:
    s = s.replace(old_p2, new_p2)

# --- 6d: n_sites_used filter on popstats
old_ps = '''p8 <- NULL; p9 <- NULL
if (!is.null(POPSTATS_INVGT_FILE) && file.exists(POPSTATS_INVGT_FILE)) {
  ps <- fread(POPSTATS_INVGT_FILE)
  setnames(ps, tolower(names(ps)))
  start_col <- intersect(c("window_start", "start_bp", "start", "ws"), names(ps))[1]
  end_col   <- intersect(c("window_end",   "end_bp",   "end",   "we"), names(ps))[1]'''
new_ps = '''p8 <- NULL; p9 <- NULL
if (!is.null(POPSTATS_INVGT_FILE) && file.exists(POPSTATS_INVGT_FILE)) {
  ps <- fread(POPSTATS_INVGT_FILE)
  setnames(ps, tolower(names(ps)))
  min_sites <- as.integer(Sys.getenv("Q04_MIN_SITES", "20"))
  n_used_col <- intersect(c("n_sites_used", "n_used", "nsites"), names(ps))[1]
  if (!is.null(n_used_col)) {
    ps <- ps[get(n_used_col) >= min_sites]
  }
  start_col <- intersect(c("window_start", "start_bp", "start", "ws"), names(ps))[1]
  end_col   <- intersect(c("window_end",   "end_bp",   "end",   "we"), names(ps))[1]'''
if old_ps in s:
    s = s.replace(old_ps, new_ps)

p.write_text(s)
print("  Q04 patched")
PY
  log "  patched (Z columns, NULL placeholders, SNP density, n_sites filter)"
else
  log "  already up-to-date"
fi

# -----------------------------------------------------------------------------
# 7. R/q08_shelf_heatmap.R — env-hash lookup
# -----------------------------------------------------------------------------
log "== Fix 7/7: R/q08_shelf_heatmap.R subscript bug =="
if grep -q 'target_lookup\[\[key\]\]' R/q08_shelf_heatmap.R 2>/dev/null; then
  backup R/q08_shelf_heatmap.R
  [[ "${DRY_RUN}" != "1" ]] && python3 - <<'PY'
import pathlib
p = pathlib.Path("R/q08_shelf_heatmap.R")
s = p.read_text()
old = '''target_lookup <- setNames(seq_along(shelf_row_idx), as.character(shelf_row_idx))
data_row_idx <- 0L
hits <- 0L
chunk_size <- 20000L
max_target <- max(shelf_row_idx)
repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0L) break
  for (ln in lines) {
    data_row_idx <- data_row_idx + 1L
    key <- as.character(data_row_idx)
    if (!is.null(target_lookup[[key]])) {
      hi <- target_lookup[[key]]'''
new = '''target_env <- new.env(hash = TRUE, size = length(shelf_row_idx) * 2L)
for (i in seq_along(shelf_row_idx)) {
  assign(as.character(shelf_row_idx[i]), i, envir = target_env)
}
data_row_idx <- 0L
hits <- 0L
chunk_size <- 20000L
max_target <- max(shelf_row_idx)
repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0L) break
  for (ln in lines) {
    data_row_idx <- data_row_idx + 1L
    key <- as.character(data_row_idx)
    hi <- target_env[[key]]
    if (!is.null(hi)) {'''
if old in s:
    p.write_text(s.replace(old, new))
    print("  Q08 patched")
else:
    print("  pattern drifted; manual review needed")
PY
  log "  patched env-hash lookup"
else
  log "  already patched"
fi

log ""
log "=========================================="
log "UPGRADE COMPLETE"
log "=========================================="
log "Backup files carry suffix: .preUPGRADE_${ts}"
log ""
log "Verify the upgrade:"
log "  bash -n *.sh                              # bash syntax"
log "  grep -c 'max_abs_z' R/q04_compose_plot.R  # should be >= 1"
log "  grep -c 'target_env' R/q08_shelf_heatmap.R  # should be >= 1"
log "  grep -c 'Engine B native' R/q06_ancestry_tracks.R  # should be = 1"
log ""
log "Then rerun the pipeline:"
log "  SHELF_START_MB=15 SHELF_END_MB=18 bash run_all.sh C_gar_LG28"
