# =============================================================================
# registry_key_helpers.R — per-writer block assembly for the 3 HELPER_MISSING
#                          schemas this file is responsible for.
#
# History
# -------
#   - Originally slated for chat 7 (Finding 9), re-slated for chat 15, never
#     built until chat B pass (2026-04-24 continuation). Meanwhile the three
#     caller scripts all guard with `exists(fn, mode="function")` so they
#     ran as silent no-ops.
#
# What this file does
# -------------------
# Provides six functions called from writer scripts under `exists()` guards:
#
#   STEP_C01g_boundary_catalog_wired_*.R  → register_C01g_boundary(bd, cid, side, outdir)
#                                           store_C01g_boundary   (bd, cid, side, outdir)
#   STEP_C01d_candidate_scoring_wired_*.R → register_C01d_keys    (cd, cid, outdir)
#                                           store_C01d_results    (cd, cid, outdir)
#   STEP_C01f_hypothesis_tests.R          → register_C01f_keys    (vd, cid, outdir)
#                                           store_C01f_results    (vd, cid, outdir)
#
# `register_*` builds a block_data list whose NAMES match each schema's `from`
# fields, then calls write_block() (via registry_bridge's reg$evidence), which
# does the `from` → `key` translation automatically using the schema's
# `keys_extracted` directive. Returns the integer key count.
#
# `store_*` writes a per-candidate diagnostic TSV/RDS under `outdir/per_cand/`
# so a reviewer can audit what went into each block. Non-critical; failures
# are swallowed so they don't block the register_* path.
#
# Design rule: SCHEMA WINS AS API
# -------------------------------
# The field NAMES written into block_data follow each schema's `from` vocabulary,
# not the writer's internal column names. Writer-side column aliases are
# translated here. That way schemas are the stable public contract; writers
# can rename their internals without breaking the key extraction.
#
# Examples:
#   C01d internal `final_score` → schema `composite_score`
#   C01d internal `span_mb`     → schema `span_kb`  (span_mb * 1000)
#   C01f internal `t9_jackknife_verdict` → schema `t9_jackknife_status`
#     (NOTE: assumes _status ≡ _verdict semantically. See TODO §1 below.)
#
# Open TODOs (flagged for Quentin; do NOT guess these — they are semantic
# decisions, not naming)
# ------------------------------------------------------------------
#
# TODO §1  `_status` vs `_verdict` synonymy
#          Currently treated as synonyms: schema `t9_jackknife_status` is
#          filled from writer `t9_jackknife_verdict`. If you meant them
#          differently, this mapping is wrong.
#
# TODO §2  `clip_count` (boundary schema)
#          Writer persists `cheat11_clip_score` (float 0..1),
#          `cheat11_clip_enrichment` (ratio), `cheat11_clip_bimodal` (bool).
#          No integer count is persisted. Schema asks for an integer. Leaving
#          as NA_integer_ until you pick a meaning:
#            (a) rename schema key to clip_score → map cheat11_clip_score
#            (b) modify C01g to persist n_clip on `deduped`
#            (c) drop the key
#
# TODO §3  `depth_anomaly` (boundary schema)
#          Writer has `cheat10_depth_dip` and `cheat10_depth_ratio` (both
#          persisted on `deduped`). Schema asks for one scalar. Leaving NA
#          until you pick one or split the schema into `_depth_dip` +
#          `_depth_ratio`.
#
# TODO §4  `is_fossil` (boundary schema)
#          Writer has only the aggregate `n_fossil` at L1394. We derive a
#          per-boundary boolean via `boundary_activity == "fossil"` when that
#          column is present on the row, NA otherwise. Derivation, not
#          mapping; remove if you disagree.
#
# TODO §5  `fdr_q_value` (existence_layer_a schema)
#          Not computed in C01d — FDR correction is typically a post-pass
#          (BH adjustment over all candidates). Left NA_real_ here; should
#          be populated by a separate pass, not by register_C01d_keys.
#
# TODO §6  `q6_family_linkage` in hypothesis_verdict
#          The schema declares `q6_family_linkage` with source `t9_jackknife_status`,
#          but C01i_d_seal ALSO writes q6_family_linkage (from `family_linkage`,
#          L433). Two writers for the same key with different semantics — whichever
#          runs last wins. This helper still writes it (schema is authoritative API)
#          but we emit a message so the collision is visible. Decision needed on
#          which writer should own the key; suggested action is to drop the entry
#          from hypothesis_verdict.schema.json.
#
# TODO §7  LOADER BUG — {side} template substitution is unimplemented
#          `registry_loader.R:load_schema("boundary_5prime")` returns NULL because
#          only boundary.schema.json exists (no boundary_5prime.schema.json).
#          Consequently write_block() sees schema=NULL and skips keys_extracted
#          entirely — no q3_5prime_* / q3_3prime_* keys land in keys.tsv. This
#          also affects `03_consensus_merge.R` / `boundary_refined_{side}`
#          (currently flagged OK_DIRECT but silently broken the same way).
#          Scoped workaround below: register_C01g_boundary loads boundary.schema.json
#          directly, does its own `{side}` substitution, and writes the flat
#          keys via add_evidence(). A proper fix requires patching load_schema +
#          write_block in registry_loader.R to handle the `{side}` templating.
#
# TODO §8  boundary_refined still uses left/right
#          `boundary_refined.schema.json` and `03_consensus_merge.R` still emit
#          `left`/`right` as side values (L473, L499 of 03_consensus_merge.R).
#          To avoid breaking JSON-schema validation against a mismatched enum,
#          boundary_refined was deliberately NOT updated in this pass. Next
#          vocabulary-sweep pass should: (a) change boundary_refined.schema.json
#          enum to [5prime, 3prime], (b) update 03_consensus_merge.R to emit
#          5prime/3prime, (c) rename the keys_extracted templates from
#          q3_refined_{side}_* (resolving to _left/_right) to the same template
#          resolving to _5prime/_3prime. Same dual-write shim pattern used here
#          would apply there.
#
# Returns and side effects
# ------------------------
#   register_* returns integer key count (for the caller's "[reg] N keys" message).
#   store_*    returns invisible(NULL); writes an RDS under
#              <outdir>/registry_helpers/<cid>[_side]/<tag>.rds.
# =============================================================================

suppressPackageStartupMessages({
  requireNamespace("data.table", quietly = TRUE)
  requireNamespace("jsonlite", quietly = TRUE)
})

# ─── tiny utilities ─────────────────────────────────────────────────────────

.rkh_safe_get <- function(row, col, default = NA) {
  # For a single-row data.table (or list / data.frame), fetch a column value
  # robustly. Returns `default` if the column is missing or the value is
  # NULL. NA values pass through unchanged.
  if (is.null(row)) return(default)
  if (!(col %in% names(row))) return(default)
  v <- tryCatch(row[[col]], error = function(e) NULL)
  if (is.null(v) || length(v) == 0) return(default)
  v[[1]]
}

.rkh_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  suppressWarnings(as.numeric(x))
}

.rkh_int <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_integer_)
  suppressWarnings(as.integer(x))
}

.rkh_chr <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  as.character(x)
}

.rkh_bool <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(NA)
  as.logical(x)
}

.rkh_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

.rkh_store_rds <- function(obj, outdir, cid, tag) {
  # Non-critical: store helper-level diagnostic RDS per candidate. Never stop
  # the caller on failure.
  tryCatch({
    d <- .rkh_ensure_dir(file.path(outdir, "registry_helpers", cid))
    saveRDS(obj, file.path(d, paste0(tag, ".rds")))
  }, error = function(e) {
    message("[rkh] store skipped for ", cid, "/", tag, ": ", conditionMessage(e))
  })
  invisible(NULL)
}

.rkh_write_block <- function(reg, cid, block_type, data, source_script) {
  # Full-reg mode (2026-04-24 migration): no JSON fallback. If the registry
  # isn't loaded, fail loudly rather than silently dumping orphan JSON files.
  if (is.null(reg) || is.null(reg$evidence) ||
      is.null(reg$evidence$write_block)) {
    stop("[rkh] reg$evidence$write_block unavailable for block ", block_type,
         " / ", cid, ". The v10 registry is mandatory (full-reg migration). ",
         "Ensure utils/registry_bridge.R is sourced before calling this helper.")
  }
  Sys.setenv(CURRENT_SCRIPT = source_script)
  res <- tryCatch(
    reg$evidence$write_block(candidate_id = cid, block_type = block_type,
                             data = data),
    error = function(e) {
      stop("[rkh] write_block failed for ", cid, "/", block_type, ": ",
           conditionMessage(e))
    }
  )
  as.integer(res$n_keys %||% 0L)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a


# =============================================================================
# NAMING CONVENTION: 5prime / 3prime
# =============================================================================
# All per-boundary labels in this helper use `5prime`/`3prime`, NOT `left`/`right`.
#
# Convention:
#   5prime = bp1 = lower coordinate on the project reference assembly
#   3prime = bp2 = higher coordinate on the project reference assembly
#
# This is 5'/3' *relative to the project reference assembly orientation*,
# NOT relative to any external reference (e.g., NCBI GCA). The pipeline has
# one working reference and this convention is defined against it. A separate
# conversion module could translate to NCBI coordinates in the future; that
# is OUT OF SCOPE for the inversion manuscript.
#
# Back-compat with legacy `left`/`right` callers:
#   - `register_C01g_boundary(side="left")` is silently translated to "5prime"
#     so the existing C01g writer (which emits bd$side = "left"/"right")
#     keeps working without edits.
#   - A one-time deprecation message is emitted per session.
#
# Dual-write shim for legacy consumer scripts:
#   - For every new `q3_5prime_*` / `q3_3prime_*` key written, the helper
#     ALSO emits the legacy `q3_left_*` / `q3_right_*` key with the same
#     value. This keeps `compute_candidate_status.R` and any other downstream
#     reader working during the migration window.
#   - Controlled by module-level flag RKH_EMIT_LEGACY_SIDE_KEYS (TRUE by default).
#     Once no consumer reads q3_left_*/q3_right_*, flip to FALSE and the legacy
#     keys disappear.
# =============================================================================

# Module-level configuration: keep legacy left/right keys during migration.
# Flip to FALSE once grep finds no q3_left_ / q3_right_ readers in the tree.
if (!exists("RKH_EMIT_LEGACY_SIDE_KEYS")) {
  RKH_EMIT_LEGACY_SIDE_KEYS <- TRUE
}

# Canonical side values
.RKH_SIDES_CANONICAL <- c("5prime", "3prime")

# Legacy-to-canonical translator (one-way).
.rkh_canonicalize_side <- function(side) {
  if (identical(side, "5prime") || identical(side, "3prime")) return(side)
  if (identical(side, "left")) {
    if (!isTRUE(getOption(".rkh_left_right_deprecation_warned", FALSE))) {
      message("[rkh] NOTE: received side='left' — translating to '5prime'. ",
              "Callers should migrate to 5prime/3prime; this message appears once per session.")
      options(.rkh_left_right_deprecation_warned = TRUE)
    }
    return("5prime")
  }
  if (identical(side, "right")) {
    return("3prime")  # already warned above on first call
  }
  stop("[rkh] unknown side value: '", side, "' — expected '5prime', '3prime', 'left', or 'right'")
}

# Canonical → legacy (used by the dual-write shim).
.rkh_legacy_side <- function(side) {
  switch(side, "5prime" = "left", "3prime" = "right",
         stop("[rkh] no legacy mapping for side='", side, "'"))
}


# =============================================================================
# C01g BOUNDARY HELPERS — boundary_5prime / boundary_3prime
# =============================================================================
# Schema:  registries/schemas/structured_block_schemas/boundary.schema.json
# Caller:  phase_4_catalog/STEP_C01g_boundary_catalog_wired_*.R around L1434-1438
# Per-side call: once with side="5prime" and once with side="3prime" per
# confirmed boundary row. Legacy callers passing "left"/"right" are auto-
# translated (see .rkh_canonicalize_side).
#
# Workaround for TODO §7 (loader {side} bug):
#   - Load boundary.schema.json manually.
#   - For each keys_extracted entry, substitute {side} → side in the `key`.
#   - Call reg$evidence$write_block so the structured JSON still lands in the
#     candidate folder (it emits a "schema not found" notice for the variant
#     block_type but writes the block anyway).
#   - Then write the flat keys via reg$add_evidence so they appear in the
#     central store.
#
# Dual-write shim:
#   When RKH_EMIT_LEGACY_SIDE_KEYS is TRUE (the default), each flat key written
#   as q3_5prime_* / q3_3prime_* is ALSO written as q3_left_* / q3_right_*
#   with the same value, so legacy consumers keep working.
# -----------------------------------------------------------------------------

register_C01g_boundary <- function(bd, cid, side, outdir) {
  side <- .rkh_canonicalize_side(side)
  stopifnot(side %in% .RKH_SIDES_CANONICAL)
  block_type <- paste0("boundary_", side)

  # ── build block_data using SCHEMA vocabulary (`from` field names) ────────
  block_data <- list(
    side                   = side,
    boundary_bp            = .rkh_int (.rkh_safe_get(bd, "boundary_bp")),
    hardness               = .rkh_num (.rkh_safe_get(bd, "hardness")),
    sharpness              = .rkh_num (.rkh_safe_get(bd, "sharpness")),
    boundary_type          = .rkh_chr (.rkh_safe_get(bd, "boundary_type")),
    boundary_verdict       = .rkh_chr (.rkh_safe_get(bd, "boundary_verdict")),
    n_cheats_supporting    = .rkh_int (.rkh_safe_get(bd, "n_cheats_supporting")),
    # cheat_*_passed booleans: writer has no direct columns, leave absent
    # (schema lists them as optional). Can be added later if writer starts
    # emitting explicit per-cheat pass flags.

    # TODO §2: clip_count intent ambiguous
    clip_count             = NA_integer_,

    # TODO §3: depth_anomaly — writer has two candidates (_dip, _ratio)
    depth_anomaly          = NA_real_,

    # TODO §4: derived from boundary_activity when present
    is_fossil              = {
      ba <- .rkh_chr(.rkh_safe_get(bd, "boundary_activity"))
      if (is.na(ba)) NA else identical(ba, "fossil")
    }
  )

  # ── emit the structured JSON (best-effort; write_block will warn if the
  #    boundary_{side}.schema.json variant doesn't exist — see TODO §7) ─────
  .rkh_write_block(
    reg = if (exists("reg")) reg else NULL,
    cid = cid, block_type = block_type, data = block_data,
    source_script = "STEP_C01g_boundary_catalog"
  )

  # ── TODO §7 workaround: manual {side} substitution + add_evidence ────────
  # Load the template schema and walk keys_extracted, substituting {side}.
  schema_path <- .rkh_find_schema("boundary")
  n_keys <- 0L
  legacy_side <- if (isTRUE(RKH_EMIT_LEGACY_SIDE_KEYS)) .rkh_legacy_side(side) else NULL

  if (!is.null(schema_path)) {
    schema <- tryCatch(jsonlite::fromJSON(schema_path, simplifyVector = FALSE),
                       error = function(e) NULL)
    if (!is.null(schema) && !is.null(schema$keys_extracted) &&
        exists("reg") && !is.null(reg$add_evidence)) {
      for (ke in schema$keys_extracted) {
        key_tmpl <- ke$key; from <- ke$from
        if (is.null(key_tmpl) || is.null(from)) next
        val <- block_data[[from]]
        if (is.null(val) || (length(val) == 1 && is.na(val))) next

        # canonical key (5prime / 3prime)
        key_canonical <- gsub("{side}", side, key_tmpl, fixed = TRUE)
        .rkh_safe_add_evidence(cid, key_canonical, val, "STEP_C01g_boundary_catalog")
        n_keys <- n_keys + 1L

        # legacy dual-write (left / right)
        if (!is.null(legacy_side)) {
          key_legacy <- gsub("{side}", legacy_side, key_tmpl, fixed = TRUE)
          .rkh_safe_add_evidence(cid, key_legacy, val, "STEP_C01g_boundary_catalog")
          # NOTE: legacy keys are written but NOT counted — n_keys reflects the
          # canonical set so the caller's "[reg] N keys" message matches schema
          # intent. The dual-write is a migration aid, not a new data product.
        }
      }
    }
  }
  as.integer(n_keys)
}

store_C01g_boundary <- function(bd, cid, side, outdir) {
  side <- .rkh_canonicalize_side(side)
  .rkh_store_rds(bd, outdir, cid, paste0("C01g_boundary_", side))
}

# Small wrapper: add_evidence with consistent error handling (used by the
# dual-write path so a failure in the legacy key doesn't kill the canonical
# key, and vice versa).
.rkh_safe_add_evidence <- function(cid, key, val, script) {
  tryCatch({
    reg$add_evidence(
      candidate_id = cid,
      key = key,
      value = as.character(val),
      script = script
    )
  }, error = function(e) {
    message("[rkh] add_evidence failed for ", cid, "/", key,
            ": ", conditionMessage(e))
  })
  invisible(NULL)
}


# =============================================================================
# C01d CANDIDATE-SCORING HELPERS — existence_layer_a + morphology
# =============================================================================
# Schemas: registries/schemas/structured_block_schemas/existence_layer_a.schema.json
#          registries/schemas/structured_block_schemas/morphology.schema.json
# Caller: phase_4_catalog/STEP_C01d_candidate_scoring_wired_*.R around L1043-1049
# Per-candidate: one call per row of cand_dt. Writes TWO blocks:
#   (a) existence_layer_a — 20 keys (12 D-dimensions, scores, shape class, ...)
#   (b) morphology         — 6 live keys (span, n_children, stripe_count,
#                            patchiness, interior_homogeneity, size_class) +
#                            4 TODO-NA keys (aspect_ratio, architecture,
#                            sample_grouping, spatial_consistency).
# Return value is the sum of keys extracted from both blocks.
# -----------------------------------------------------------------------------

register_C01d_keys <- function(cd, cid, outdir) {

  # ── schema-vocabulary block_data for existence_layer_a ───────────────────
  # Translations documented inline: schema_name ← writer_name
  span_mb_val <- .rkh_num(.rkh_safe_get(cd, "span_mb"))
  span_kb_val <- if (is.na(span_mb_val)) NA_real_ else round(span_mb_val * 1000, 1)

  final_score_val <- .rkh_num(.rkh_safe_get(cd, "final_score"))

  block_data <- list(
    # composite_score ← final_score  (C01d L361, L436)
    composite_score         = final_score_val,
    dim_positive            = .rkh_int(.rkh_safe_get(cd, "dim_positive")),
    tier                    = .rkh_int(.rkh_safe_get(cd, "tier")),

    d1_block_strength       = .rkh_num(.rkh_safe_get(cd, "d1_block_strength")),
    d2_block_shape          = .rkh_num(.rkh_safe_get(cd, "d2_block_shape")),
    d3_nn_persistence       = .rkh_num(.rkh_safe_get(cd, "d3_nn_persistence")),
    d4_decay_flatness       = .rkh_num(.rkh_safe_get(cd, "d4_decay_flatness")),
    d5_interior_quality     = .rkh_num(.rkh_safe_get(cd, "d5_interior_quality")),
    d6_consensus            = .rkh_num(.rkh_safe_get(cd, "d6_consensus")),
    d7_sv_breakpoint        = .rkh_num(.rkh_safe_get(cd, "d7_sv_breakpoint")),
    d8_peel_or_hyp          = .rkh_num(.rkh_safe_get(cd, "d8_peel_or_hyp")),
    d9_pca_clusters         = .rkh_num(.rkh_safe_get(cd, "d9_pca_clusters")),
    d10_partition           = .rkh_num(.rkh_safe_get(cd, "d10_partition")),
    d11_boundary_concordance= .rkh_num(.rkh_safe_get(cd, "d11_boundary_concordance")),
    d12_snake_concordance   = .rkh_num(.rkh_safe_get(cd, "d12_snake_concordance")),

    shape_class             = .rkh_chr(.rkh_safe_get(cd, "shape_class")),
    landscape_category      = .rkh_chr(.rkh_safe_get(cd, "landscape_category")),
    n_children              = .rkh_int(.rkh_safe_get(cd, "n_children")),

    # span_kb ← span_mb * 1000  (C01d L414)
    span_kb                 = span_kb_val,

    cheat25_status          = .rkh_chr(.rkh_safe_get(cd, "cheat25_status")),

    # TODO §5: FDR q-value not computed in C01d — left NA, populate in post-pass
    fdr_q_value             = NA_real_,

    # Record which pass produced this row, if available
    pass_number             = .rkh_int(.rkh_safe_get(cd, "pass_number"))
  )

  n_a <- .rkh_write_block(
    reg = if (exists("reg")) reg else NULL,
    cid = cid, block_type = "existence_layer_a", data = block_data,
    source_script = "STEP_C01d_candidate_scoring"
  )

  # ── morphology block (chat B continuation, 2026-04-24) ──────────────────
  # The morphology schema declares 10 keys. Of those:
  #   - span_kb, n_children, stripe_count, patchiness_score,
  #     interior_homogeneity, size_class: MATCH or DRIFT fields that can be
  #     derived from cand_dt (now includes homogeneity/patchiness/stripe_count
  #     as of the 2026-04-24 pass-through added to C01d).
  #   - architecture, sample_grouping, spatial_consistency, aspect_ratio:
  #     MISSING — writer has no computation for these. Left as NA with
  #     per-field TODOs in the block. The schema's `required` includes
  #     `architecture` so the block will emit `validation_status=incomplete`
  #     until a writer for these fields appears.
  morph_block <- .rkh_build_morphology_block(cd)
  n_m <- .rkh_write_block(
    reg = if (exists("reg")) reg else NULL,
    cid = cid, block_type = "morphology", data = morph_block,
    source_script = "STEP_C01d_candidate_scoring"
  )

  as.integer(n_a + n_m)
}

# ─── morphology block builder ───────────────────────────────────────────────
# Separated from register_C01d_keys for clarity. Builds the block_data list
# using the morphology schema's `from` vocabulary. Fields flagged TODO are
# populated as NA until a writer appears; the schema-validation status will
# be "incomplete" because `architecture` is required but NA.
.rkh_build_morphology_block <- function(cd) {
  # Derivable fields
  span_mb_val <- .rkh_num(.rkh_safe_get(cd, "span_mb"))
  span_kb_val <- if (is.na(span_mb_val)) NA_real_ else round(span_mb_val * 1000, 1)
  span_bp_val <- if (is.na(span_mb_val)) NA_integer_ else as.integer(round(span_mb_val * 1e6))

  # interior_homogeneity: use cand_dt's pass-through `homogeneity` column.
  # C01d uses this same value to compute D5_interior_quality (L218, L237).
  interior_homog <- .rkh_num(.rkh_safe_get(cd, "homogeneity"))

  # patchiness_score ← `patchiness` (writer column name)
  patchiness_val <- .rkh_num(.rkh_safe_get(cd, "patchiness"))

  # stripe_count: pass-through
  stripe_val <- .rkh_int(.rkh_safe_get(cd, "stripe_count"))

  # n_children: pass-through
  n_children_val <- .rkh_int(.rkh_safe_get(cd, "n_children"))

  # size_class: derive from span_mb per schema description (small <1Mb,
  # medium 1-5Mb, large >5Mb). Writer has no explicit size_class column;
  # derivation is deterministic so we compute it here.
  size_class_val <- if (is.na(span_mb_val)) NA_character_
                    else if (span_mb_val < 1.0) "small"
                    else if (span_mb_val < 5.0) "medium"
                    else "large"

  list(
    span_kb              = span_kb_val,
    span_bp              = span_bp_val,

    # TODO §A1 (morphology): aspect_ratio — writer has no aspect_ratio
    # column. Could be derived from block dimensions if n_windows and
    # n_samples are known; not computed in C01d today.
    aspect_ratio         = NA_real_,

    interior_homogeneity = interior_homog,
    patchiness_score     = patchiness_val,
    stripe_count         = stripe_val,
    n_children           = n_children_val,

    # TODO §A2: architecture enum (single_block / nested_blocks /
    # multi_system / complex). Would need a classifier on (n_children,
    # shape_class, snake_overlap) — heuristic could be written but not
    # in this pass.
    architecture         = NA_character_,

    size_class           = size_class_val,

    # TODO §A3: sample_grouping enum (simple_polymorphic / nested /
    # overlapping / complex). Needs cross-sample decomposition info not
    # present in cand_dt.
    sample_grouping      = NA_character_,

    # TODO §A4: spatial_consistency enum (uniform / gradient / patchy).
    # Requires windowed homogeneity scan across the block interior;
    # C01d has the data (iv$homogeneity per-window) but doesn't carry
    # the profile through to cand_dt. Requires either an upstream
    # pass-through or a post-pass.
    spatial_consistency  = NA_character_
  )
}

store_C01d_results <- function(cd, cid, outdir) {
  .rkh_store_rds(cd, outdir, cid, "C01d_candidate_scoring")
}


# =============================================================================
# C01f HYPOTHESIS-VERDICT HELPERS — hypothesis_verdict
# =============================================================================
# Schema: registries/schemas/structured_block_schemas/hypothesis_verdict.schema.json
# Caller: phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R L2535-2539
# Per-candidate: one call per row of verd_dt.
#
# Note: frequency block updates happen in a separate helper
# (update_C01f_frequency_linkage, defined below). This function handles only
# the hypothesis_verdict block. The v1/v2 frequency schema collision was
# resolved in frequency.v3.schema.json (2026-04-24 reconciliation).
# -----------------------------------------------------------------------------

register_C01f_keys <- function(vd, cid, outdir) {

  # ── schema-vocabulary block_data for hypothesis_verdict ──────────────────

  # TODO §1: `_status` ≡ `_verdict` assumption. Writer has both possibilities:
  # `t9_jackknife_verdict` and a plain `jackknife_verdict`. Prefer the t9-prefixed
  # one when present (it's the canonical t9 output), else fall back.
  t9_status_val <- .rkh_safe_get(vd, "t9_jackknife_verdict")
  if (is.na(t9_status_val %||% NA) || is.null(t9_status_val)) {
    t9_status_val <- .rkh_safe_get(vd, "jackknife_verdict")
  }
  t9_status_val <- .rkh_chr(t9_status_val)

  # TODO §1: same principle for verdict_confidence — writer has `confidence`
  # (bare). Treating them as the same field pending clarification.
  vc_val <- .rkh_safe_get(vd, "verdict_confidence")
  if (is.na(vc_val %||% NA) || is.null(vc_val)) {
    vc_val <- .rkh_safe_get(vd, "confidence")
  }
  vc_val <- .rkh_chr(vc_val)

  block_data <- list(
    verdict                   = .rkh_chr(.rkh_safe_get(vd, "verdict")),
    verdict_confidence        = vc_val,
    t1_ratio                  = .rkh_num(.rkh_safe_get(vd, "t1_ratio")),
    t2_eff_k                  = .rkh_num(.rkh_safe_get(vd, "t2_eff_k")),
    t3_retention              = .rkh_num(.rkh_safe_get(vd, "t3_retention")),
    t8_concordance            = .rkh_num(.rkh_safe_get(vd, "t8_concordance")),

    # schema: t9_jackknife_status  ←  writer: t9_jackknife_verdict   (TODO §1)
    t9_jackknife_status       = t9_status_val,
    t9_max_delta              = .rkh_num(.rkh_safe_get(vd, "t9_max_delta")),
    t10_theta_concordance     = .rkh_num(.rkh_safe_get(vd, "t10_theta_concordance")),

    t11_has_extended          = .rkh_bool(.rkh_safe_get(vd, "t11_has_extended")),
    t11_extent_kb             = .rkh_num (.rkh_safe_get(vd, "t11_extent_kb")),
    t11_fst_decay_rate        = .rkh_num (.rkh_safe_get(vd, "t11_fst_decay_rate")),

    group_validation_after    = .rkh_chr(.rkh_safe_get(vd, "group_validation_after")),
    group_validation_before   = .rkh_chr(.rkh_safe_get(vd, "group_validation_before")),
    comp_source               = .rkh_chr(.rkh_safe_get(vd, "comp_source"))
  )

  # TODO §6: q6_family_linkage collision flag
  # The schema maps q6_family_linkage ← t9_jackknife_status (same field as
  # q7_t9_jackknife_status). C01i_d_seal also writes q6_family_linkage from a
  # different source. Emit a one-line notice the first time we see this per
  # session so the collision is visible in logs.
  if (!isTRUE(getOption(".rkh_c01f_collision_warned", FALSE))) {
    message("[rkh/C01f] NOTE: hypothesis_verdict writes q6_family_linkage from ",
            "t9_jackknife_status; C01i_d_seal also writes q6_family_linkage ",
            "from `family_linkage`. Last writer wins. See TODO §6.")
    options(.rkh_c01f_collision_warned = TRUE)
  }

  .rkh_write_block(
    reg = if (exists("reg")) reg else NULL,
    cid = cid, block_type = "hypothesis_verdict", data = block_data,
    source_script = "STEP_C01f_hypothesis_tests"
  )
}

store_C01f_results <- function(vd, cid, outdir) {
  .rkh_store_rds(vd, outdir, cid, "C01f_hypothesis_verdict")
}


# ─── schema path resolver ───────────────────────────────────────────────────

.rkh_find_schema <- function(block_type) {
  # Try a few canonical locations for a schema file.
  candidates <- c(
    file.path("registries/schemas/structured_block_schemas",
              paste0(block_type, ".schema.json")),
    file.path("../registries/schemas/structured_block_schemas",
              paste0(block_type, ".schema.json")),
    file.path(Sys.getenv("BASE", ""),
              "registries/schemas/structured_block_schemas",
              paste0(block_type, ".schema.json")),
    # and toolkit-relative path if we're inside inversion-popgen-toolkit
    file.path(Sys.getenv("BASE", ""),
              "inversion-popgen-toolkit/registries/schemas/structured_block_schemas",
              paste0(block_type, ".schema.json"))
  )
  for (p in candidates) {
    if (nzchar(p) && file.exists(p)) return(normalizePath(p, mustWork = FALSE))
  }
  NULL
}

# =============================================================================
# C01i_d_seal + C01f FREQUENCY HELPERS — frequency block (Option A architecture)
# =============================================================================
# Schema: registries/schemas/structured_block_schemas/frequency.v3.schema.json
# Architecture: seal writes the block with what seal knows; C01f updates
#               family_linkage + polymorphism_class after the jackknife runs.
#
# Writers:
#   register_C01i_frequency_block(seal_record, cid, outdir)
#     Called from STEP_C01i_d_seal.R per-candidate seal loop.
#     Writes: freq_inv, freq_class, n_total, class_counts, expected_counts_hwe,
#             hwe_chi2, hwe_p, hwe_verdict, genotype_balance, family_linkage
#             (placeholder "unknown"), polymorphism_class (placeholder
#             "unclassified"), computed_from="C01i_d_seal_initial".
#     Derives HWE from class_counts locally (~20 lines, chi-squared 2df).
#     Derives freq_class from freq_inv.
#
#   update_C01f_frequency_linkage(vd, cid, outdir)
#     Called from STEP_C01f_hypothesis_tests.R per-candidate verdict loop.
#     Reads the existing frequency block via reg$evidence$read_block.
#     Updates: family_linkage (from vd$family_linkage),
#              jackknife_max_delta (from vd$t9_max_delta),
#              jackknife_n_contributing (derived from vd if present),
#              polymorphism_class (derived from family_linkage — helper-side,
#              per the design note in v3 schema; keeps C01f source untouched),
#              computed_from = "C01f_updated".
#     Writes the block back. If no block exists yet (seal hasn't run), writes
#     a minimal block with NA for freq fields and the jackknife-based info
#     only.
#
#   store_C01i_frequency_block / store_C01f_frequency_update
#     Diagnostic RDS writers, per convention.

# ─── freq_class derivation ──────────────────────────────────────────────────
# Matches characterize_candidate.R L555-559 so the schema contract agrees
# with the consumer's threshold ladder.
.rkh_classify_freq <- function(freq) {
  if (is.null(freq) || length(freq) == 0 || !is.finite(freq) || is.na(freq)) {
    return(NA_character_)
  }
  if (freq < 0.05) return("rare")
  if (freq < 0.15) return("low")
  if (freq < 0.50) return("intermediate")
  if (freq < 0.85) return("high")
  "nearly_fixed"
}

# ─── HWE chi-squared test from genotype counts ──────────────────────────────
# Returns list(expected = list(HOM_REF, HET, HOM_INV), chi2, p, verdict).
# verdict: "consistent" if p >= 0.05; "het_excess"/"het_deficit" if p < 0.05
# based on sign of (obs_het - exp_het); "underpowered" if n < 30 or any
# expected cell < 5.
#
# Recombinants excluded from both observed and expected (consistent with
# freq_inv denominator logic).
.rkh_hwe_test <- function(n_ref, n_het, n_inv, n_rec = 0) {
  n_diploid <- n_ref + n_het + n_inv  # recombinants excluded
  out <- list(
    expected = list(HOM_REF = NA_real_, HET = NA_real_, HOM_INV = NA_real_),
    chi2     = NA_real_,
    p        = NA_real_,
    verdict  = "underpowered"
  )
  if (n_diploid < 1) return(out)

  # Inverted-allele frequency (recombinants out of denominator)
  p <- (2 * n_inv + n_het) / (2 * n_diploid)
  q <- 1 - p

  exp_ref <- n_diploid * q * q
  exp_het <- n_diploid * 2 * p * q
  exp_inv <- n_diploid * p * p
  out$expected <- list(HOM_REF = round(exp_ref, 3),
                       HET     = round(exp_het, 3),
                       HOM_INV = round(exp_inv, 3))

  # Underpowered: small n or any expected cell < 5
  if (n_diploid < 30 ||
      exp_ref < 5 || exp_het < 5 || exp_inv < 5) {
    return(out)  # verdict stays "underpowered"
  }

  # chi2, 2 df
  parts <- c(
    if (exp_ref > 0) (n_ref - exp_ref)^2 / exp_ref else 0,
    if (exp_het > 0) (n_het - exp_het)^2 / exp_het else 0,
    if (exp_inv > 0) (n_inv - exp_inv)^2 / exp_inv else 0
  )
  chi2 <- sum(parts)
  pval <- pchisq(chi2, df = 2, lower.tail = FALSE)

  verdict <- if (pval >= 0.05) "consistent"
             else if (n_het > exp_het) "het_excess"
             else "het_deficit"

  out$chi2    <- round(chi2, 4)
  out$p       <- signif(pval, 4)
  out$verdict <- verdict
  out
}

# ─── polymorphism_class derivation (helper-side per v3 schema design) ──────
# Maps family_linkage -> polymorphism_class:
#   multi_family          -> cohort_wide
#   few_family            -> lineage_restricted
#   single_family         -> family_restricted
#   pca_family_confounded -> unclassified
#   uncertain, unknown    -> unclassified
.rkh_derive_polymorphism_class <- function(family_linkage) {
  if (is.null(family_linkage) || is.na(family_linkage) ||
      !nzchar(family_linkage)) {
    return("unclassified")
  }
  switch(family_linkage,
         multi_family          = "cohort_wide",
         few_family            = "lineage_restricted",
         single_family         = "family_restricted",
         "unclassified")  # default: pca_family_confounded, uncertain, unknown
}

# =============================================================================
# register_C01i_frequency_block — seal-side writer
# =============================================================================
# Caller: STEP_C01i_d_seal.R, inside the per-candidate seal loop, after
# class_counts have been computed (currently around L369 after freq_inv is
# derived). The seal script already has all the inputs; just pass in a list
# with the count fields plus cid + outdir.
#
# Expected input (seal_record): a list or single-row data.table/data.frame
# with fields: n_HOM_REF, n_HET, n_HOM_INV, n_RECOMBINANT (optional), n_total,
# freq_inv, family_linkage (optional; seal passes "unknown"), polymorphism_class
# (optional; seal passes "unclassified").
#
# The call site in seal would be:
#   if (exists("register_C01i_frequency_block", mode = "function")) {
#     register_C01i_frequency_block(
#       list(n_HOM_REF=n_ref, n_HET=n_het, n_HOM_INV=n_inv,
#            n_RECOMBINANT=n_rec, n_total=n_tot, freq_inv=freq_inv),
#       cid, outdir
#     )
#   }

register_C01i_frequency_block <- function(seal_record, cid, outdir) {
  n_ref <- .rkh_int(.rkh_safe_get(seal_record, "n_HOM_REF"))
  n_het <- .rkh_int(.rkh_safe_get(seal_record, "n_HET"))
  n_inv <- .rkh_int(.rkh_safe_get(seal_record, "n_HOM_INV"))
  n_rec <- .rkh_int(.rkh_safe_get(seal_record, "n_RECOMBINANT"))
  if (is.na(n_rec)) n_rec <- 0L

  n_tot_given <- .rkh_int(.rkh_safe_get(seal_record, "n_total"))
  n_tot <- if (is.na(n_tot_given)) {
    sum(c(n_ref, n_het, n_inv, n_rec), na.rm = TRUE)
  } else {
    n_tot_given
  }

  freq_inv_given <- .rkh_num(.rkh_safe_get(seal_record, "freq_inv"))
  # Recompute freq_inv if not supplied (defensive).
  n_diploid <- sum(c(n_ref, n_het, n_inv), na.rm = TRUE)
  freq_inv <- if (!is.na(freq_inv_given)) {
    freq_inv_given
  } else if (n_diploid > 0) {
    round((2 * n_inv + n_het) / (2 * n_diploid), 4)
  } else {
    NA_real_
  }

  # HWE
  hwe <- .rkh_hwe_test(n_ref %||% 0L, n_het %||% 0L, n_inv %||% 0L, n_rec)

  # genotype_balance mirrors hwe_verdict minus "underpowered" — forces a call
  genotype_balance <- switch(hwe$verdict,
    "consistent"   = "hwe_consistent",
    "het_excess"   = "het_excess",
    "het_deficit"  = "het_deficit",
    "underpowered" = {
      # Even when underpowered, we still record directional info for
      # descriptive purposes. Sign of (observed - expected) heterozygotes.
      if (!is.na(hwe$expected$HET) && !is.na(n_het)) {
        if (n_het > hwe$expected$HET) "het_excess" else "het_deficit"
      } else {
        "hwe_consistent"  # no data: default
      }
    },
    "hwe_consistent"
  )

  block_data <- list(
    freq_inv         = freq_inv,
    freq_class       = .rkh_classify_freq(freq_inv),
    n_total          = as.integer(n_tot),
    class_counts     = list(
      HOM_REF     = n_ref %||% NA_integer_,
      HET         = n_het %||% NA_integer_,
      HOM_INV     = n_inv %||% NA_integer_,
      RECOMBINANT = n_rec %||% NA_integer_
    ),
    expected_counts_hwe = hwe$expected,
    hwe_method       = "chi_squared",
    hwe_chi2         = hwe$chi2,
    hwe_p            = hwe$p,
    hwe_verdict      = hwe$verdict,
    genotype_balance = genotype_balance,

    # Placeholders — C01f overwrites in update_C01f_frequency_linkage
    family_linkage            = .rkh_chr(.rkh_safe_get(seal_record, "family_linkage"))
                                  %||% "unknown",
    jackknife_max_delta       = NA_real_,
    jackknife_n_contributing  = NA_integer_,
    polymorphism_class        = .rkh_chr(.rkh_safe_get(seal_record, "polymorphism_class"))
                                  %||% "unclassified",

    computed_from   = "C01i_d_seal_initial",
    note_on_hatchery_hwe = paste0(
      "F1 hybrid hatchery samples systematically violate HWE at neutral loci ",
      "due to population admixture. het_excess is expected; het_deficit is ",
      "the more informative deviation."
    )
  )

  .rkh_write_block(
    reg = if (exists("reg")) reg else NULL,
    cid = cid, block_type = "frequency", data = block_data,
    source_script = "STEP_C01i_d_seal"
  )
}

store_C01i_frequency_block <- function(seal_record, cid, outdir) {
  .rkh_store_rds(seal_record, outdir, cid, "C01i_frequency_block")
}

# =============================================================================
# update_C01f_frequency_linkage — C01f-side updater
# =============================================================================
# Caller: STEP_C01f_hypothesis_tests.R, inside the per-candidate verdict loop
# (currently around L2535, same place as register_C01f_keys). Takes the vd
# row and updates family_linkage + derived polymorphism_class in the existing
# frequency block.
#
# If the frequency block doesn't yet exist (seal hasn't run, or seal failed),
# writes a minimal block with NA for freq-side fields, just the jackknife
# info. This is rare in practice (seal runs before C01f in the orchestrator).

update_C01f_frequency_linkage <- function(vd, cid, outdir) {
  # Read existing block (if any)
  existing <- tryCatch(
    if (!is.null(reg) && !is.null(reg$evidence$read_block)) {
      reg$evidence$read_block(cid, "frequency")
    } else {
      NULL
    },
    error = function(e) NULL
  )

  # Extract new values from vd
  family_linkage <- .rkh_chr(.rkh_safe_get(vd, "family_linkage")) %||% "unknown"
  if (!family_linkage %in% c("multi_family", "few_family", "single_family",
                              "pca_family_confounded", "uncertain", "unknown")) {
    # Unrecognized; coerce to unknown rather than violating the enum
    message("[rkh/C01f] NOTE: unexpected family_linkage='", family_linkage,
            "' for ", cid, " — coercing to 'unknown'")
    family_linkage <- "unknown"
  }

  jk_max_delta <- .rkh_num(.rkh_safe_get(vd, "t9_max_delta"))
  # jackknife_n_contributing: derived if available. vd may carry an explicit
  # column (not currently in C01f output); else leave NA_integer_.
  jk_n_contrib <- .rkh_int(.rkh_safe_get(vd, "jackknife_n_contributing"))
  if (is.na(jk_n_contrib)) {
    jk_n_contrib <- .rkh_int(.rkh_safe_get(vd, "t9_n_contributing"))
  }

  polymorphism_class <- .rkh_derive_polymorphism_class(family_linkage)

  if (!is.null(existing) && !is.null(existing$data)) {
    # Merge: keep seal's counts/HWE, overwrite C01f's fields.
    block_data <- existing$data
    block_data$family_linkage            <- family_linkage
    block_data$jackknife_max_delta       <- jk_max_delta
    block_data$jackknife_n_contributing  <- jk_n_contrib
    block_data$polymorphism_class        <- polymorphism_class
    block_data$computed_from             <- "C01f_updated"
  } else {
    # No prior block — write a minimal C01f-only block. freq / HWE fields
    # are NA because seal didn't populate them for this candidate.
    message("[rkh/C01f] WARNING: no existing frequency block for ", cid,
            " — writing minimal C01f-only block (seal may not have run)")
    block_data <- list(
      freq_inv         = NA_real_,
      freq_class       = NA_character_,
      n_total          = NA_integer_,
      class_counts     = list(HOM_REF = NA_integer_, HET = NA_integer_,
                               HOM_INV = NA_integer_, RECOMBINANT = NA_integer_),
      hwe_method       = "chi_squared",
      hwe_chi2         = NA_real_,
      hwe_p            = NA_real_,
      hwe_verdict      = "underpowered",
      genotype_balance = "hwe_consistent",
      family_linkage            = family_linkage,
      jackknife_max_delta       = jk_max_delta,
      jackknife_n_contributing  = jk_n_contrib,
      polymorphism_class        = polymorphism_class,
      computed_from   = "C01f_updated"  # only C01f has written; seal missing
    )
  }

  .rkh_write_block(
    reg = if (exists("reg")) reg else NULL,
    cid = cid, block_type = "frequency", data = block_data,
    source_script = "STEP_C01f_hypothesis_tests"
  )
}

store_C01f_frequency_update <- function(vd, cid, outdir) {
  .rkh_store_rds(
    list(family_linkage = .rkh_safe_get(vd, "family_linkage"),
         t9_max_delta   = .rkh_safe_get(vd, "t9_max_delta")),
    outdir, cid, "C01f_frequency_update"
  )
}


message("[registry_key_helpers] loaded (10 functions: ",
        "register_C01g_boundary, store_C01g_boundary, ",
        "register_C01d_keys, store_C01d_results, ",
        "register_C01f_keys, store_C01f_results, ",
        "register_C01i_frequency_block, store_C01i_frequency_block, ",
        "update_C01f_frequency_linkage, store_C01f_frequency_update)")
