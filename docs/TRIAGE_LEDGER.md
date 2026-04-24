# Schema <-> Writer triage ledger

**Generated:** chat B, post-audit verification.  
**Scope:** every schema in `registries/schemas/structured_block_schemas/` mapped to its intended writer script (per `source_script` field) and the actual state of that writer in the live tree.

## TL;DR

37 schemas, 273 declared keys. Four distinct failure modes, not one:

| Status | Schemas | Keys | Meaning |
|---|---|---|---|
| Wired | 18 | 127 | `OK_DIRECT` + `OK_VIA_SAFE` — keys flow |
| Helper missing | 6 | 74 | chat-7 Finding 9 — `utils/registry_key_helpers.R` never built |
| Script runs, doesn't wire | 5 | 24 | output files exist, no `write_block` call added |
| Script is a function library | 2 | 21 | called by orchestrator but has no CLI entry; effectively dead |
| No script declared | 6 | 27 | forward-declared schemas or `source_script` missing |

The single highest-leverage fix is **building `utils/registry_key_helpers.R`** — unblocks 74 keys (27% of all declared) across 6 schemas, closes a 2-chat-old deferred item, no architectural changes needed.

Separately, the audit's "7 schemas missing `block_type`" is a documentation drift, not a runtime issue: `load_schema()` looks up by filename, not by the internal `block_type` field. 4 of those 7 have working writers that the audit missed because their `source_script` is blank.

## Spec-coverage view

`build_key_spec()` in `phase_9_classification/compute_candidate_status.R` declares **406 live keys + 96 aspirational** across Q1–Q7 + Q_QC_SHELF (paste0-expanded; keys like `q3_left_fst_m200kb` are counted).

- **133 spec keys have a writer source** (schema `keys_extracted`, `add_evidence`, or the env-gated `_qc_shelf_reader.R` appender)
- **185 spec keys are real orphans** (in spec, no writer, not aspirational)
- **164 keys land in `keys.tsv` but nobody registered them in the spec** — read-only via `read_block()`

Per-axis orphan count:

| Axis | Spec | Wired | Orphan | Aspir |
|---|---|---|---|---|
| q1 | 49 | 20 | 29 | 0 |
| q2 | 91 | 55 | 36 | 0 |
| q3 | 61 | 12 | 40 | 9 |
| q4 | 47 | 2 | 27 | 19 |
| q5 | 39 | 2 | 26 | 11 |
| q6 | 31 | 15 | 0 | 16 |
| q7 | 75 | 14 | 27 | 41 |
| q_qc_shelf | 13 | 13 | 0 | 0 |

`q6` and `q_qc_shelf` have zero orphans — every live spec key has a writer. `q4`/`q5` are the worst because their writers are in the HELPER_MISSING + SCRIPT_IS_LIB categories.

---

## Verdict legend

| Verdict | Meaning | Fix posture |
|---|---|---|
| `OK_DIRECT` | writer calls `write_block(cid, "<bt>", ...)` directly | nothing to do |
| `OK_VIA_SAFE` | writer uses `write_block_safe()` helper (passthrough) | nothing to do |
| `HELPER_MISSING` | script calls `register_C01x_*`/`store_C01x_*` but the helper file (`utils/registry_key_helpers.R`) was never built — silent no-op at runtime | create the helper file (original chat-7/chat-11 plan) |
| `SCRIPT_IS_LIB` | script exists but has no CLI entry point; defines functions only, orchestrator calls `Rscript script.R` which just loads definitions and exits | wire a runner (add `commandArgs` + main dispatch) |
| `PRODUCES_BUT_NOT_WIRED` | script runs and produces output files (RDS/TSV/etc.), but never calls any registry function | add `write_block`/`add_evidence` calls for the outputs |
| `SCRIPT_STUB` | script exists, has a main entry, but produces no file outputs | likely abandoned or incomplete — check intent |
| `NO_SCRIPT_DECLARED` | `source_script` is `(none)` and no candidate writer found | schema is forward-declared for future work, or obsolete |

## Headline counts

| Verdict | Schemas | Keys affected |
|---|---|---|
| `OK_DIRECT` | 12 | 80 |
| `HELPER_MISSING` | 6 | 74 |
| `NO_SCRIPT_DECLARED` | 6 | 27 |
| `OK_VIA_SAFE` | 6 | 47 |
| `PRODUCES_BUT_NOT_WIRED` | 5 | 24 |
| `SCRIPT_IS_LIB` | 2 | 21 |
| **TOTAL** | **37** | **273** |

## OK — writers wired and producing keys

| Schema | Q | Keys | Writer(s) |
|---|---|---|---|
| `ancestral_fragments_summary` | Q3 | 6 | `inversion_modules/phase_6_breakpoint_refinement/02_ancestral_fragments.R` |
| `bnd_sided_support` | Q7 | 5 | `inversion_modules/phase_8_evidence_biology/q7_existence_audit/bnd_sided_support.py` |
| `boundary_refined` | Q3 | 6 | `inversion_modules/phase_6_breakpoint_refinement/03_consensus_merge.R` |
| `breakpoints_per_method` | Q3 | 1 | `inversion_modules/phase_6_breakpoint_refinement/03_consensus_merge.R` |
| `distance_concordance` | — | 3 | `inversion_modules/phase_4_catalog/STEP_C01m_distance_concordance.R` |
| `dosage_blocks` | Q3 | 6 | `inversion_modules/phase_6_breakpoint_refinement/01_dosage_signal.R` |
| `encoding_robustness` | Q2 | 3 | `inversion_modules/phase_6_breakpoint_refinement/05_four_encoding_diagnostic.R` |
| `existence_layer_b_bnd_rescue` | Q7 | 11 | `inversion_modules/phase_3_refine/STEP_B06_bnd_rescue.py` |
| `existence_layer_d` | Q7 | 8 | `inversion_modules/phase_3_refine/STEP_D03_statistical_tests_and_seeds.py` |
| `fragment_distribution` | Q2+Q3 | 10 | `inversion_modules/phase_8_evidence_biology/bp_bridge/bp_pipeline_bridge.py` |
| `internal_ancestry_composition` | Q1, Q2 | 12 | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py` |
| `internal_dynamics` | Q1, Q2 | 19 | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R` |
| `local_structure_segments` | — | 6 | `inversion_modules/phase_4_catalog/STEP_C01l_local_structure_segments.R` |
| `mechanism_assembled` | Q4 | 6 | `inversion_modules/phase_8_evidence_biology/q4_mechanism/cheat29b_assembled_junction.py` |
| `recombinant_map` | Q2, Q4 | 10 | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R` |
| `regime_sample_dag` | — | 5 | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R` |
| `regime_segments` | — | 4 | `inversion_modules/phase_4_catalog/STEP_C01j_regime_compatibility_engine.R` |
| `synteny_v6` | Q5 | 6 | `inversion_modules/phase_8_evidence_biology/q5_age_and_origin/cross_species_bridge_v6.py` |

## HELPER_MISSING — blocked on the chat-7 deferred helper file

These scripts use the pattern:

```r
for (.hf in c("utils/registry_key_helpers.R", "../utils/registry_key_helpers.R")) {
  if (file.exists(.hf)) { source(.hf); break }
}
if (exists("register_C01x_keys", mode = "function")) { ... }
```

The helper file was intentionally not built in chat 7 (flagged as Finding 9), was slated for chat 15 per `HANDOFF_PROMPT_chat11`, but chat 15 did BK schemas instead. Nothing has built it since. Silent no-op at runtime.

**Fix scope:** build `utils/registry_key_helpers.R` with ~6 functions: `register_C01d_keys`, `store_C01d_results`, `register_C01f_keys`, `store_C01f_results`, `register_C01g_boundary`, `store_C01g_boundary`. Signatures are already determined by their callers.

| Schema | Q | Keys | Declared source_script | Found scripts |
|---|---|---|---|---|
| `boundary` | Q3 | 9 | STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R | `inversion_modules/phase_4_catalog/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` |
| `existence_layer_a` | Q1, Q7 | 23 | C01a + C01b + C01d | `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`<br>`inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` |
| `frequency` | Q6 | 8 | C01i + C01f | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R`<br>`inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R` |
| `frequency.v2` | Q6 | 10 | C01i_d_seal + C01f | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R`<br>`inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R` |
| `hypothesis_verdict` | Q7 | 14 | STEP_C01f_hypothesis_tests_wired_6_9_12_23_v934_registry.R | `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R` |
| `morphology` | Q1 | 10 | C01a precompute + C01d scoring | `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R`<br>`inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` |

## SCRIPT_IS_LIB — function-library-only, no runner

These scripts define R functions and exit. The orchestrator calls them with `Rscript script.R --candidate X --outdir Y`, but the script has no `commandArgs` handler; R parses the file, registers the function definitions, and exits with no output. Effectively dead. These are good candidates for a thin CLI wrapper.

| Schema | Q | Keys | Declared source_script | Found scripts |
|---|---|---|---|---|
| `age_evidence` | Q5 | 11 | cheat30_gds_by_genotype + unified_ancestry/region_popstats | `inversion_modules/phase_8_evidence_biology/q5_age_and_origin/cheat30_gds_by_genotype.R` |
| `mechanism` | Q4 | 10 | cheat14 + cheat21 + cheat27 + cheat28 + cheat29 | `inversion_modules/phase_8_evidence_biology/q4_mechanism/cheat27_sd_nahr_substrate.R`<br>`inversion_modules/phase_8_evidence_biology/q4_mechanism/cheat28_tandem_repeat_context.R`<br>`inversion_modules/phase_8_evidence_biology/q4_mechanism/cheat29_junction_forensics.R` |

## PRODUCES_BUT_NOT_WIRED — output files exist, never flow to registry

Script runs, writes output files (RDS, TSV, etc.), but never calls `write_block` / `add_evidence`. Keys declared by the schema are not populated.

**2026-04-24 (chat C) resolution pass:**
- `band_composition` — FOLDED into `internal_dynamics`. Archived at `_archive_superseded/bk_schemas_pre_canonical/band_composition_folded_into_internal_dynamics_2026-04-24/`. `q1_k_used` dual-written alongside `q2_pca_k_used`; `q1_ancestry_div_hom_ref_vs_hom_inv` declared NA pending instant_q join inside decompose.
- `flank_coherence` — FOLDED into `synteny_v6`. Archived at `_archive_superseded/bk_schemas_pre_canonical/flank_coherence_folded_into_synteny_v6_2026-04-24/`. `cross_species_bridge_v6.py` already emitted the fields; schema's `keys_extracted` list extended to pick them up.
- `block_detect`, `existence_layer_b`, `existence_layer_c` — remain BLOCKED with documented reasons. Each has a `REGISTRY_CONTRACT` header in its source script explaining the architectural blocker (producers run before candidate_ids exist). See `docs/SCRIPT_CONTRACT.md`.

| Schema | Q | Keys | Declared source_script | Status |
|---|---|---|---|---|
| `band_composition` | Q1 | 2 | C01i + unified_ancestry/instant_q | **FOLDED → internal_dynamics (2026-04-24)** |
| `block_detect` | Q1, Q3 | 6 | PHASE_01C block_detect | BLOCKED_ON_NO_CANDIDATE_JOIN + 5 missing producer fields |
| `existence_layer_b` | Q7 | 6 | C00 build_flashlight + MODULE_4H_ALL_Manta + MODULE_4_DELLY_INV | BLOCKED_ON_NO_CANDIDATE_JOIN |
| `existence_layer_c` | Q7 | 7 | C04 GHSL v5 (within-sample haplotype divergence) + STEP_C01c triangle insulation | BLOCKED_ON_NO_CANDIDATE_JOIN |
| `flank_coherence` | Q3 | 3 | cross_species_bridge.py | **FOLDED → synteny_v6 (2026-04-24)** |

## NO_SCRIPT_DECLARED — schema is orphan / forward-declared

`source_script` is `(none)` and no candidate writer file exists. These schemas might be obsolete or forward-declarations for future work.

| Schema | Q | Keys | Declared source_script | Found scripts |
|---|---|---|---|---|
| `burden` | Q6 | 6 | burden_regression + STEP_C01f_c_burden_regression.R | — |
| `carrier_reconciliation` | Q3, Q6, Q7 | 3 | audit_cores_vs_sv (new) + carrier_audit.R | — |
| `gene_conversion_tracts` | — | 3 | (none) | — |
| `inv_detect_blocks` | Q1, Q3 | 6 | inv_detect_v9.3 + 12_bridge_to_codebase.R | — |
| `peel_diagnostic` | Q7 | 4 | C01n peel / inv_detect Phase 9 | — |
| `synteny_dollo` | Q5 | 5 | synteny_inversion_detection.sh + dollo_parsimony_inversions.py + mashmap | — |

## Ad-hoc `add_evidence` writers — flat keys that bypass schemas

These keys land in `keys.tsv` via `add_evidence()` (single-key writes), not through a schema's `keys_extracted`. Consumers can read them, but they are not formally part of any schema's API.

| Key | Writer script |
|---|---|
| `figure_breakpoint_diagnostic_pdf_path` | `inversion_modules/phase_6_breakpoint_refinement/04_diagnostic_figure.R` |
| `figure_breakpoint_diagnostic_png_path` | `inversion_modules/phase_6_breakpoint_refinement/04_diagnostic_figure.R` |
| `figure_stream_graph_bands_tsv_path` | `inversion_modules/phase_6_breakpoint_refinement/06_regime_stream_graph.R` |
| `figure_stream_graph_pdf_path` | `inversion_modules/phase_6_breakpoint_refinement/06_regime_stream_graph.R` |
| `figure_stream_graph_png_path` | `inversion_modules/phase_6_breakpoint_refinement/06_regime_stream_graph.R` |
| `q2_decomp_quality_flags` | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R` |
| `q6_family_linkage` | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R` |
| `q6_group_validation` | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R` |
| `q6_polymorphism_class` | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R` |
| `q6_validation_promotion_cap` | `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R` |

## Dead helper-function references

These function names are called with `exists(fn, mode="function")` guards, but no definition of them exists anywhere in the tree. The guards silently skip the call. All six live in `utils/registry_key_helpers.R`, which was never built (see HELPER_MISSING section).

| Function | Called by |
|---|---|
| `register_C01d_keys` | `inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` |
| `register_C01f_keys` | `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`, `inversion_modules/phase_7_karyotype_groups/validation/_PRISTINE_v9.3.4.R`, `inversion_modules/phase_9_classification/patches/01_C01f_comp_from_registry.R` |
| `register_C01g_boundary` | `inversion_modules/phase_4_catalog/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` |
| `store_C01d_results` | `inversion_modules/phase_4_catalog/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` |
| `store_C01f_results` | `inversion_modules/phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`, `inversion_modules/phase_7_karyotype_groups/validation/_PRISTINE_v9.3.4.R` |
| `store_C01g_boundary` | `inversion_modules/phase_4_catalog/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` |
