# AUDIT_LOG_chat13_2026-04-18.md

Session log for chat 13 (registry wiring pass after chat-12 audit).

**Date:** 2026-04-18
**Preceded by:** chat 12 (phase-4 coherence audit + DAG rewrite).
**Follows into:** chat 14 (first HPC run on LANTA).

---

## Primary deliverables

1. **All four target block writers wired.**
   - `STEP_C01j_regime_compatibility_engine.R` emits
     `write_block_safe(cid, "regime_segments", ...)` per candidate inside
     the main loop. `dominant_state` computed by window-weighted mode;
     `regime_has_recombinant_signal` from structured↔structured
     transitions; `n_samples_with_regime_change` from stability classes.
     Schema fields (`segments`, `transitions`) serialised as list-of-list
     for JSON via a local `.as_block_list()` helper.
   - `STEP_C01l_local_structure_segments.R` emits
     `write_block_safe(cid, "local_structure_segments", ...)` after the
     segment loop. `boundary_sharpness.overall` ∈ {sharp, moderate,
     diffuse, asymmetric} computed from core-vs-flank Δ12 deltas
     (thresholds: diff > 0.30 = sharp, > 0.10 = moderate, else diffuse;
     asymmetric when left and right classifications differ).
   - `STEP_C01m_distance_concordance.R` emits
     `write_block_safe(cid, "distance_concordance", ...)` at the end of
     each chromosome loop iteration. Interval-restricted aggregates:
     mean concordance at each distance; n_pairs_high and n_pairs_decaying
     counters; persistent_carriers list; and `inversion_vs_family_score =
     longest_mean / shortest_mean` clamped [0, 1].
   - `STEP_C01i_b_multi_recomb.R` already emits `regime_sample_dag`
     correctly from chat 12 — confirmed no change needed.

2. **Seal wiring:** `STEP_C01i_d_seal.R` now registers
   `inv_<cid>_RECOMBINANT_GC` and `inv_<cid>_RECOMBINANT_DCO` groups
   using multi_recomb's canonical `recomb_subgroup` field, with
   event-class-based fallback. `recomb_subgroup` threaded through
   `final_records` so `register_all_groups()` can read it.

3. **Characterize Q2 wiring:** `characterize_candidate.R::characterize_q2`
   pulls the 14 new keys (regime_segments 4, local_structure_segments 6,
   distance_concordance 3, regime_sample_dag 5 — with overlap on q2_regime_*)
   as annotation-only evidence lines. Existing status/label decision
   logic unchanged per handoff: "Wire the new blocks as annotation-only
   first; promote to gates after HPC calibration in chat 14."

4. **Threshold + assertion fixes:**
   - **AT:** C01l `flank_bp` scales per candidate:
     `cand_flank_bp = clamp(span_bp, 200000L, 500000L)`.
     `define_segments(..., flank_bp = cand_flank_bp)`.
   - **AV:** decompose emits non-gating `decomp_quality ∈ {"clean",
     "noisy", NA}` based on `silhouette >= 0.40`.
   - **BB:** `lib_ghsl_confirmation.R::resolve_karyo_bp` now has
     `stopifnot(uniqueN(annot_dt$global_window_id) == nrow(annot_dt))`
     at the top with inline rationale.

5. **Doc fixes:**
   - **AO:** C01d header: "10 dimensions from detector" →
     "12 dimensions from detector" with chat-13 AO note.
   - **AY:** seed_loader header annotates the implemented behaviour as
     STRICTER than the chat-9 "priority: flashlight > STEP03" language —
     drop-on-conflict, not priority-pick.

6. **Path verification (resolved as no-op):**
   - **AQ:** C01e Panel H reads `regime_windows.tsv.gz` via `--regime_dir`.
     C01j still writes that filename at the same relative location inside
     `--outdir`; chat-11.5 directory move did not change the written name.
     No hardcoded stale path anywhere in C01e.
   - **BD:** C01k reads `regime_state_labels.tsv.gz` via `--regime_dir`.
     Same story as AQ — sidecar TSV path unchanged. Keeping TSV reads
     rather than migrating to registry block reads is consistent with
     the handoff's "annotation-only first" stance; migration can be a
     chat-15 cleanup.

7. **BC wired-fraction metric:**
   - Registry-wide: **73 / 319 = 22.9%** of spec keys have a schema
     writer producing them.
   - Q2 specifically: **22 / 73 = 30.1%** of q2 spec keys have a schema
     writer.
   - Of those 22 covered Q2 keys, 18 are the chat-13 newly-wired keys.
   - The handoff's aspirational 85% target is for Q2 coverage after ALL
     phase-4b schemas expose `keys_extracted` directives; see Finding BK.
   - The BG dotted-path fix (session 2) unlocked **21 previously-dropped
     keys** across 4 schemas — materially more impact than chat 12
     anticipated when it flagged the issue. Keys now reach keys.tsv:
     5 from regime_sample_dag, 6 from local_structure_segments,
     9 from age_evidence, 1 from gene_conversion_tracts.

---

## Findings closed this session

- **V** / **AO** / **AQ** / **AT** / **AV** / **AY** / **BB** / **BC** /
  **BD** — all handoff items applied or verified.
- **C01i_b regime_sample_dag round-trip** — confirmed correct from
  chat-12 code; no change.

---

## Findings opened this session (BG onwards)

| # | Severity | Owner | Summary |
|---|----------|-------|---------|
| BG | bug | chat 13 | `registry_loader.R::write_block` key extractor used `data[[from]]`, no dotted-path support. Fixed with `resolve_dotted()` helper; vector leaves collapse to `length()`. Unlocked 21 keys across 4 schemas. |
| BH | bug (blocker) | chat 13 | C01j header advertised `regime_memberships.tsv.gz` as an output, but no writer existed. C01i_b requires that file via `--regime_memb`. Pipeline could not run end-to-end until this was fixed. Now writes via `rbindlist(all_memberships, fill = TRUE)` at the end of the main loop. |
| BI | bug | chat 13 | Seal's `RECOMBINANT_GC` classifier matched on `event_class == "gene_conversion"`, but C01i_b (chat-12 rewrite) emits `"gene_conversion_embedded"`. Samples that should have landed in RECOMBINANT_GC were silently falling through to the unspecified-subgroup bucket. Fixed to accept both names; prefer `recomb_subgroup` field directly when present. |
| BJ | spec maintenance | chat 14 | `build_key_spec()::q2` did not list the 18 keys chat-13 is now producing. Added. Any future new block writer must also land its keys in the spec or completion accounting under-reports. |
| BK | spec maintenance | chat 14/15 | Four phase-4b schemas (`internal_dynamics`, `recombinant_map`, `internal_ancestry_composition`, and nested_composition's own output) have NO `keys_extracted` directives even though their block data contains 51 of the Q2 spec keys. Adding these directives is the path to the handoff's 85% Q2 aspiration. Not in chat-13 scope (the keys exist in the JSON; they just don't get flattened to keys.tsv). |

---

## Parse-check status

No R interpreter available in this environment. All edited files were
checked with `_rcheck.py` (custom brace/paren/string balance checker at
the work root). Files that passed:

- `registries/api/R/registry_loader.R` (BG fix)
- `phase_4_postprocessing/4a_existence_layers/STEP_C01j_regime_compatibility_engine.R` (wiring + BH)
- `phase_4_postprocessing/4a_existence_layers/STEP_C01l_local_structure_segments.R` (wiring + AT)
- `phase_4_postprocessing/4a_existence_layers/STEP_C01m_distance_concordance.R` (wiring)
- `phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` (AO)
- `phase_4_postprocessing/4b_group_proposal/STEP_C01i_d_seal.R` (BI + recomb_subgroup threading)
- `phase_4_postprocessing/4b_group_proposal/STEP_C01i_decompose.R` (AV)
- `phase_4_postprocessing/4b_group_proposal/lib_ghsl_confirmation.R` (BB)
- `phase_4_postprocessing/4b_group_proposal/lib_step03_seed_loader.R` (AY doc)
- `phase_4_postprocessing/4e_final_classification/characterize_candidate.R` (Q2 wiring)
- `phase_4_postprocessing/4e_final_classification/compute_candidate_status.R` (BJ)

---

## Test status

**Chat-12 tests NOT re-run** (no R interpreter in this environment). The
50/50 DAG derive_R + 33/33 GC detector test suites from chat 12 cover
compute logic in `lib_recomb_combination.R` and `gene_conversion_detector.R`
— neither touched in chat 13, so no behavioural regression expected.

**Chat 14 MUST re-run both suites** under the HPC R to confirm no
regressions. Parse-check passes on every edited file; tests should be
clean.

---

## Posture for chat 14

Chat 13 is a pure wiring pass. Every compute path from chat 12 is
preserved; chat 13 only added output routing, per-candidate block
writes, and a handful of small fixes/docs. The pipeline should now run
end-to-end for the first time (BH was a blocker — without it, C01i_b's
regime_memb input did not exist).

Chat 14 is the first HPC run on LANTA:
1. LG12 first (known strong inversion). Verify DAG fires R correctly,
   verify boundary_sharpness comes back "sharp", verify
   `inversion_vs_family_score > 0.5`.
2. LG25 second (paralog-heavy control). Verify the null — DAG R should
   not fire broadly, boundary_sharpness should be "diffuse" or
   "asymmetric", inv_vs_family_score should be low.
3. Calibrate **AS** (structure_score cutoffs), **AW** (min_dco_bp),
   **AV** (silhouette threshold for decomp_quality) against observed
   distributions.
4. Re-run 50/50 + 33/33 test suites under HPC R to confirm no
   regressions from chat 13.
