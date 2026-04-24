# HANDOFF — chat B continuation (2026-04-24)

Consolidated handoff for everything done in this session and everything deferred to next.

---

## What was done

### 1. HELPER_MISSING bucket — CLOSED (6 of 6 schemas)

| Schema | Helper pair | Status |
|---|---|---|
| `boundary` | `register_C01g_boundary` / `store_C01g_boundary` | Wired + 5prime/3prime vocabulary + dual-write shim for legacy q3_left_*/q3_right_* consumers |
| `existence_layer_a` | `register_C01d_keys` / `store_C01d_results` | Wired; 20 keys populated |
| `hypothesis_verdict` | `register_C01f_keys` / `store_C01f_results` | Wired; 10 keys populated |
| `frequency` (v3) | `register_C01i_frequency_block` / `store_C01i_frequency_block` | Wired from seal side |
| `frequency` (v3) update path | `update_C01f_frequency_linkage` / `store_C01f_frequency_update` | Wired from C01f side |
| `morphology` | folded into `register_C01d_keys` | Wired; 6 live keys + 4 NA-TODO |

`utils/registry_key_helpers.R` now has 10 public functions + 1 internal `.rkh_build_morphology_block`.

### 2. Frequency schema reconciliation

v1 and v2 were colliding on `block_type="frequency"` with different semantics. Reconciled into `frequency.v3.schema.json`. Merge details in `_archive_superseded/bk_schemas_pre_canonical/frequency_reconciled_2026-04-24/README.md`.

Key decisions baked in:
- **5-class `freq_class`** (`rare/low/intermediate/high/nearly_fixed` at 0.05/0.15/0.50/0.85) — matches `characterize_candidate.R`.
- **6-value `family_linkage`** (added `uncertain` to v2's 5-value enum; writer emits it at C01f L553).
- **polymorphism_class derivation helper-side** in `update_C01f_frequency_linkage` (C01f source untouched).
- **flat class-count keys** (`q6_n_HOM_REF`, `q6_n_HOM_STD`, `q6_n_HET`, `q6_n_HOM_INV`, `q6_n_Recombinant`) — fixes the gap flagged in `SPEC_VS_REALITY.md` line 362.
- **HWE computation** (chi-squared, 2df) done locally in the helper from genotype counts.

Architecture: **Option A** — seal writes partial block + C01f updates `family_linkage`/`polymorphism_class` after jackknife. `computed_from` field records provenance.

### 3. 5prime / 3prime vocabulary

`boundary.schema.json` migrated from `left`/`right` to `5prime`/`3prime`. Convention: 5prime = bp1 = lower coordinate on the **project reference assembly** (NOT NCBI). Explicitly documented in schema description.

Back-compat:
- Caller-side: legacy `"left"/"right"` values auto-translated with one-time deprecation message.
- Consumer-side: dual-write shim emits both canonical `q3_5prime_*` and legacy `q3_left_*` keys. Flag `RKH_EMIT_LEGACY_SIDE_KEYS <- TRUE`.

### 4. Full-reg migration — ALL FALLBACKS REMOVED

Every v9/v10 dual-path and every "JSON fallback if registry unavailable" path gone. List:

| File | What was removed |
|---|---|
| `utils/registry_bridge.R` | v9 flat `sample_registry.R` loader fallback |
| `inversion_modules/phase_7_karyotype_groups/proposal/lib_decompose_helpers.R` | v9 loader path in `try_load_registry`, JSON fallback in `write_block_safe`/`read_block_safe`, dead `outdir_fallback` param (and removed it from all 8 call sites) |
| `unified_ancestry/dispatchers/region_stats_dispatcher.R` | flat-loader fallback |
| `registries/api/R/registry_loader.R` | silent no-op shim in `load_samples_api` |
| `utils/registry_key_helpers.R` | orphan JSON dump in `.rkh_write_block`; removed `outdir` param from signature and 5 internal call sites |

If the v10 registry can't load, scripts now **crash visibly** at startup instead of silently producing orphan JSON files nobody reads.

### 5. Call-site cleanup (Part D)

All 8 call sites passing `outdir_fallback = ...` stripped. Signatures updated. No callers, no deprecated param.

Edited files:
- `inversion_modules/phase_4_catalog/STEP_C01l_local_structure_segments.R` (1)
- `inversion_modules/phase_4_catalog/STEP_C01j_regime_compatibility_engine.R` (1)
- `inversion_modules/phase_4_catalog/STEP_C01m_distance_concordance.R` (1)
- `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R` (3)
- `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R` (2)

### 6. `add_evidence` audit (Part E)

Audited all 19 `reg$add_evidence` / `reg$evidence$add_evidence` call sites. Classification:

- **Bucket 1 — figure path keys** (5 sites): `figure_breakpoint_diagnostic_pdf_path`, `figure_stream_graph_*_path` etc. Legitimate flat keys (Tier-3 pointers to rendered figures). **Kept.**
- **Bucket 2 — stale duplicates of frequency block placeholders** (2 sites in `STEP_C01i_d_seal.R`): `q6_family_linkage = "unknown"`, `q6_polymorphism_class = "unclassified"`. These are now written inside the `frequency` block and extracted as flat keys via `keys_extracted`. **Deleted.**
- **Bucket 3 — genuine flat scalar keys** (3 sites): `q6_group_validation`, `q2_decomp_quality_flags`, `q6_validation_promotion_cap`. Single scalars with no natural block home. **Kept.**
- **Bucket 4 — guard checks** (9 sites): `!is.null(reg$add_evidence)` conditions, not writes. **Kept.**

Net effect: removed 2 redundant writes; the other 17 are legitimate. The "19 sites of duplication" framing was misleading — real duplication was just 2.

### 7. `check_deployment_v10.py`

Updated to reference `frequency.v3.schema.json` instead of v2.

---

## Files changed (complete list)

**Created:**
- `utils/registry_key_helpers.R` (~975 lines, 10 public functions + internal helpers)
- `registries/schemas/structured_block_schemas/frequency.v3.schema.json`
- `_archive_superseded/bk_schemas_pre_canonical/frequency_reconciled_2026-04-24/README.md`
- `docs/HANDOFF_2026-04-24_chat_B.md` (this file)

**Edited:**
- `registries/schemas/structured_block_schemas/boundary.schema.json` (5prime/3prime vocabulary)
- `utils/registry_bridge.R` (v9 fallback stripped)
- `inversion_modules/phase_7_karyotype_groups/proposal/lib_decompose_helpers.R` (fallbacks stripped; signature simplified)
- `unified_ancestry/dispatchers/region_stats_dispatcher.R` (flat fallback stripped)
- `registries/api/R/registry_loader.R` (silent shim → `stop()`)
- `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R` (frequency block call site added; 2 stale placeholders removed)
- `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R` (frequency update call site added)
- `inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` (morphology pass-through columns added to cand_dt)
- `inversion_modules/phase_4_catalog/STEP_C01l_local_structure_segments.R` (outdir_fallback stripped)
- `inversion_modules/phase_4_catalog/STEP_C01j_regime_compatibility_engine.R` (outdir_fallback stripped)
- `inversion_modules/phase_4_catalog/STEP_C01m_distance_concordance.R` (outdir_fallback stripped)
- `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R` (3× outdir_fallback stripped)
- `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R` (2× outdir_fallback stripped)
- `inversion_modules/check_deployment_v10.py` (v2 → v3)

**Archived:**
- `_archive_superseded/bk_schemas_pre_canonical/frequency_reconciled_2026-04-24/frequency.v1.schema.json` (was `frequency.schema.json`)
- `_archive_superseded/bk_schemas_pre_canonical/frequency_reconciled_2026-04-24/frequency.v2.schema.json`

---

## TO DO (next sessions)

### B. PRODUCES_BUT_NOT_WIRED schemas — 5 schemas, 24 keys

Different bucket than HELPER_MISSING. These scripts run and produce output files (RDS / TSV) but never call any registry function, so declared keys never flow to `keys.tsv`.

- `band_composition` (2 keys, `STEP_C01i_decompose.R`)
- `block_detect` (6 keys, `PHASE_01C_block_detect.R`)
- `existence_layer_b` (6 keys, `STEP_C00_build_sv_prior.R`)
- `existence_layer_c` (7 keys, `STEP_C04_snake3_ghsl_v6.R`)
- `flank_coherence` (3 keys, `cross_species_bridge_v6.py`)

Each needs:
1. Decision on "fold into existing block vs standalone" (`band_composition` may want to merge into `internal_dynamics`; `flank_coherence` may want to merge into `synteny_v6`)
2. Writer patch — add `write_block` call at the right spot
3. Field mapping resolution (most have DRIFT/MISSING fields against the writer, similar to HELPER_MISSING)

Estimate: 3–5 hours, medium risk — touches 5 scripts actively used in the manuscript pipeline.

### C. Helper TODOs from first pass

These are semantic decisions, not mechanical work:

1. **§2 `clip_count`** — schema has it, writer has only `cheat11_clip_score` / `cheat11_clip_enrichment` / `cheat11_clip_bimodal`. Pick: map schema to `cheat11_clip_score`, modify writer to persist `n_clip`, or drop the schema key.
2. **§3 `depth_anomaly`** — writer has `cheat10_depth_dip` AND `cheat10_depth_ratio`. Schema has one slot. Pick one, combine, or split schema.
3. **§5 `fdr_q_value`** — not computed in C01d; is a post-pass (BH adjustment over all candidates). Confirm this is a separate pipeline step, not C01d's job.
4. **§6 `q6_family_linkage` collision** — `hypothesis_verdict` schema AND `frequency.v3` schema both define it. Last-write-wins; `frequency.v3` is the canonical source (from the jackknife via `update_C01f_frequency_linkage`). Suggested fix: drop the key from `hypothesis_verdict.schema.json` `keys_extracted`.
5. **§7 LOADER BUG** — `{side}` template substitution unimplemented in `registry_loader.R::load_schema`. Scoped workaround lives in `register_C01g_boundary` (loads `boundary.schema.json` manually, substitutes `{side}`, writes via `add_evidence`). Proper fix patches `load_schema` + `write_block` to handle the substitution. Benefit: retroactively fixes `boundary_refined` too. Estimate: 1 hour + targeted test.
6. **§8 `boundary_refined` vocabulary** — still uses `left`/`right`. `03_consensus_merge.R` hardcodes them at L473, L499. Coordinated pass: (a) update schema enum + templates, (b) patch writer to emit 5prime/3prime. Estimate: 30 min.

### Morphology TODOs (introduced this pass)

Four fields in `morphology.schema.json` are populated as NA because the writer doesn't compute them. They have per-field TODOs in `.rkh_build_morphology_block`:

- **§A1 `aspect_ratio`** — would need block-dimension info (n_windows × n_samples). Not in cand_dt.
- **§A2 `architecture`** — enum classifier on `(n_children, shape_class, snake_overlap)`. Could be written as a heuristic.
- **§A3 `sample_grouping`** — needs cross-sample decomposition info not in cand_dt.
- **§A4 `spatial_consistency`** — needs windowed homogeneity profile. Data exists in `iv` but not carried to `cand_dt`.

Until these are computed, the morphology block will be marked `validation_status = incomplete` by `write_block` because `architecture` is required.

---

## Verification plan (on LANTA)

```bash
# 1. Sanity: helpers load, 10 functions visible
cd $BASE
Rscript -e '
source("utils/registry_bridge.R")
source("utils/registry_key_helpers.R")
fns <- c("register_C01g_boundary","store_C01g_boundary",
         "register_C01d_keys","store_C01d_results",
         "register_C01f_keys","store_C01f_results",
         "register_C01i_frequency_block","store_C01i_frequency_block",
         "update_C01f_frequency_linkage","store_C01f_frequency_update")
stopifnot(all(vapply(fns, exists, logical(1), mode="function")))
cat("OK: all 10 helpers visible\n")
stopifnot(!is.null(reg$evidence$write_block))
cat("OK: v10 registry active\n")
'

# 2. Run C01d on LG28 — expect existence_layer_a + morphology blocks written
sbatch --wrap='Rscript inversion_modules/phase_4_catalog/STEP_C01d_*.R --chr LG28'

# 3. Run C01g on LG28 — expect boundary_5prime + boundary_3prime
# 4. Run C01i_d_seal — expect frequency block (computed_from=C01i_d_seal_initial)
# 5. Run C01f — expect hypothesis_verdict + frequency update (computed_from=C01f_updated)

# 6. Spot-check a known candidate
CID=<pick one>  # e.g., C_gar_LG28_inv_001
cat <cand_dir>/$CID/structured/existence_layer_a.json
cat <cand_dir>/$CID/structured/morphology.json
cat <cand_dir>/$CID/structured/frequency.json       # should have computed_from=C01f_updated
cat <cand_dir>/$CID/structured/boundary_5prime.json
cat <cand_dir>/$CID/structured/boundary_3prime.json

# 7. keys.tsv should contain dual-write rows
grep -E '^q3_(5prime|3prime|left|right)_' <cand_dir>/$CID/keys.tsv | wc -l
# Expect ~18 (9 keys × 2 sides, each × 2 for dual-write = 36 - but last-write-wins
# dedups on identical (key, cid), so actual row count is 18)

# 8. freq keys (new)
grep '^q6_' <cand_dir>/$CID/keys.tsv
# Expect q6_freq_inv, q6_freq_class, q6_n_total, q6_n_HOM_REF, q6_n_HOM_STD,
# q6_n_HET, q6_n_HOM_INV, q6_n_Recombinant, q6_expected_het_hwe, q6_hwe_chi2,
# q6_hwe_p, q6_hwe_verdict, q6_genotype_balance, q6_family_linkage,
# q6_jackknife_max_delta, q6_jackknife_n_contributing, q6_polymorphism_class
```

---

## Expected keys.tsv row count per candidate (rough order of magnitude)

Pre-pass baseline: ~50 rows (boundary, existence_layer_a subset, hypothesis_verdict subset, scattered add_evidence).

Post-pass target:
- `boundary` canonical + legacy: ~18 rows
- `existence_layer_a`: 20 rows
- `morphology`: ~6 rows (4 NA-skipped)
- `hypothesis_verdict`: 10 rows
- `frequency`: 17 rows
- Other blocks + add_evidence: ~10 rows

Rough total: **~80 rows per fully-processed candidate.** Compare against pre-pass on the same candidate to verify ~30 extra rows appeared from the newly-wired blocks.

---

## Risk items for LANTA deployment

1. **Crash on missing bridge.** If `utils/registry_bridge.R` fails to load (missing config, unreachable schema dir, etc.), scripts now hard-fail at startup. Previously they'd run to completion with silent orphan JSON dumps. Watch the first SLURM run's stderr.

2. **`boundary_refined` still on `left`/`right`.** `03_consensus_merge.R` hardcodes these. The `boundary_refined.schema.json` enum still lists `["left", "right"]`. Unchanged — works as before (still silently broken same as prior to this pass, per TODO §7). Will finish in the TODO §7 + §8 combined follow-up.

3. **`{side}` loader bug** still affects both `boundary` and `boundary_refined`. The helper's scoped workaround handles `boundary`; `boundary_refined` is still affected. See TODO §7.

4. **Dual-write doubles some keys.** `q3_5prime_bp` and `q3_left_bp` both exist. Legacy consumers keep working. Flip `RKH_EMIT_LEGACY_SIDE_KEYS <- FALSE` when migration is done.
