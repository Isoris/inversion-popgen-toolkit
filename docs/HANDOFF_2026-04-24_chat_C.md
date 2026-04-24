# HANDOFF 2026-04-24 — chat C → chat D

## State of the tree

Working on `inversion-popgen-toolkit` (private repo, flat layout).
Baseline is the `inversion-popgen-toolkit.tar` uploaded at the start of
chat B. This handoff's delta sits in `/home/claude/work/delta/` ready to
drop onto the baseline. DELETIONS.txt in that dir lists the 6 files to
remove (see below).

LANTA is still under maintenance; this session was codebase-only.

## What chat C accomplished

### 1. nested_composition dedup (done, verified)

Four byte-identical copies of the ancestry-composition algorithm →
one canonical engine + one inversion wrapper.

- Canonical engine: `unified_ancestry/engines/nested_composition/internal_ancestry_composition.py`
  (renamed from `nested_composition.py`). Has new reusable library
  function `analyze_parent_composition()` extracted from the duplicated
  logic. CLI behavior preserved.
- Wrapper kept: `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py`
  imports engine via `ANCESTRY_ENGINE_DIR` sys.path. Its `analyze_parent`
  is now a thin shim adding `composite_flag` interpretation on top of the
  engine's output.
- Deleted:
  - `inversion_modules/phase_7_karyotype_groups/proposal/nested_composition_core.py`
  - `inversion_modules/phase_2_discovery/_archive/analysis/nested_composition.py`
  - `unified_ancestry/engines/nested_composition/nested_composition.py` (old name)
- Stale refs fixed in: `check_deployment_v10.py`, `tools/_code_field_check.py`,
  `test_phase4b_integration.py`, `docs/ENGINES.md`.
- Dangling config fixed: `unified_ancestry/00_ancestry_config.sh` —
  `NESTED_COMPOSITION_PY` now points at the engine.
- Smoke-tested: homogeneous → clean, mixed 3/1/1 → maybe_composite,
  all 12 schema `from` fields populated.

### 2. SCRIPT_CONTRACT convention introduced (done)

New header block pattern for every registry-participating script.
Scanner catches drift between contract `keys:` list and schema
`keys_extracted`. BLOCKED statuses document exact architectural reasons.

- `docs/SCRIPT_CONTRACT.md` — spec.
- `tools/scan_script_contracts.py` — enforcement tool. Supports:
  - `BLOCKS_WRITTEN: none` for compute engines
  - `{side}` template expansion for boundary
  - Fallback chain `block_type` → `title` → filename in
    `enumerate_schemas()` for schemas missing `block_type`
- `docs/NESTED_VS_COMPOSITE.md` — design note on the nested/composite
  overload + role split + 20% threshold calibration concern at 9×
  coverage.

### 3. Schema folds

Two schemas retired by folding their keys into larger sibling schemas
that shared their writer. Consistent pattern: dissolve redundant schema,
archive standalone with README, extend target schema's `keys_extracted`,
archive with implementation sketch for deferred computation.

**`flank_coherence` → `synteny_v6`** (done, writer already emits):
- Added 5 keys to `synteny_v6.schema.json`: q3_left_flank_coherence_score,
  q3_right_flank_coherence_score, q3_left_flank_n_families,
  q3_right_flank_n_families, q3_flank_coherence_class.
- Writer `cross_species_bridge_v6.py` (L350-371) already emitted the
  fields. No writer change needed.
- Archived at `_archive_superseded/bk_schemas_pre_canonical/flank_coherence_folded_into_synteny_v6_2026-04-24/`.

**`band_composition` → `internal_dynamics`** (done, 1 field NA pending join):
- Quentin picked **option (a)**: fold not new block.
- New keys in `internal_dynamics.schema.json`:
  - `q1_decompose_k_used` ← `k_used` (dual-write with existing
    `q2_pca_k_used`; same K, different q-axis meaning — q2 asks "how
    did decompose cluster samples?", q1 asks "how many bands does the
    block have?")
  - `q1_decompose_ancestry_div_hom_ref_vs_hom_inv` ←
    `ancestry_div_hom_ref_vs_hom_inv` (NEW field, NA until instant_q
    join lands)
- Naming convention: explicit `decompose_` prefix shows provenance in
  keys.tsv. Future ancestry-divergence values from a different script
  would use that script's prefix (e.g. `q1_unified_ancestry_div_*`).
- New property `ancestry_div_hom_ref_vs_hom_inv` added to schema's
  `properties` block, nullable number, fully documented.
- `STEP_C01i_decompose.R` L556-570: `block_data$ancestry_div_hom_ref_vs_hom_inv = NA_real_`
  in seeded_ok path only (no_seeding stubs stay minimal) with TODO
  comment pointing at archive README.
- Archived at `_archive_superseded/bk_schemas_pre_canonical/band_composition_folded_into_internal_dynamics_2026-04-24/`
  with 98-line README covering decision rationale, name-change rationale,
  deferred instant_q join sketch (~30 lines), no-action-required note.
- TRIAGE_LEDGER updated: both folds marked FOLDED with pointers.

### 4. Bonus drift fix caught during fold audit

`internal_dynamics.schema.json` had `flashlight_mode` as the property
name + `from` field, but decompose actually writes `sv_prior_mode` (the
field was renamed upstream but the schema wasn't updated). Result:
`q2_pca_seed_mode` in keys.tsv was silently empty. Renamed the property
and the `from` value. This is the same class of bug to look for in
other schemas — the scanner would now catch it via the contract header.

### 5. Regime-family contracts (4 of 5 scripts contracted)

All 4 writes already existed from chat 13 (pass 15) — just needed
REGISTRY_CONTRACT headers so the scanner can see them.

| Script | Block(s) | Keys | Status |
|---|---|---|---|
| `STEP_C01j_regime_compatibility_engine.R` | regime_segments | 4 | WIRED |
| `STEP_C01l_local_structure_segments.R` | local_structure_segments | 6 (dotted-path) | WIRED |
| `STEP_C01m_distance_concordance.R` | distance_concordance | 3 (array→count via registry loader) | WIRED |
| `STEP_C01i_b_multi_recomb.R` | recombinant_map + regime_sample_dag | 10 + 5 | WIRED |

All key/from pairs cross-checked against live block_data — zero drift.

### 6. Schema block_type consistency fix

4 regime-family schemas used `title` instead of `block_type`:
regime_segments, regime_sample_dag, local_structure_segments,
distance_concordance. Added `block_type` field to each (kept `title`
too). Also hardened scanner's `enumerate_schemas()` with title/$id/
filename fallback chain.

Scanner then surfaced 3 MORE schemas still missing `block_type` —
they were previously invisible:
- `boundary.schema.json`
- `boundary_refined.schema.json`
- `gene_conversion_tracts.schema.json`

These still need `block_type` added. Trivial fix, just wasn't done
this session.

### 7. Removed redundant frequency schemas

- `registries/schemas/structured_block_schemas/frequency.schema.json` (v1) — deleted (chat B had
  archived the content but left the live file).
- `registries/schemas/structured_block_schemas/frequency.v2.schema.json` — same.
- Only `frequency.v3.schema.json` remains live.

## Scanner final state

```
scripts with REGISTRY_CONTRACT: 11     (started session at 6)
schemas in registry:           34     (27 before title-fallback surfaced 7 more)
BLOCKS_WRITTEN entries:
  BLOCKED_ON_NO_CANDIDATE_JOIN  2
  NONE_DECLARED                 1
  WIRED                        13     (started at 6)
UNCOVERED SCHEMAS              19     (28 at start)
PARSE ERRORS                    0
HARD DRIFT                      0
```

## PENDING — next chat (chat D) should do these

### High priority — finish what we started

**[P1] C01i_d_seal contract** — ONLY remaining script from the 5
regime-family batch. Different shape from the others: writes flat keys
via `add_evidence` (not `write_block`) + calls the frequency-block
registration helper. Keys to declare:
  - `q6_group_validation`
  - `q1_composite_flag`
  - `q2_decomp_quality_flags`
  - `q6_validation_promotion_cap`
Plus the frequency block it registers via `register_C01i_frequency_block`.

The contract may need a new status tag or a different section (not
`BLOCKS_WRITTEN` — maybe `KEYS_WRITTEN`) since it's using a different
API. See what makes sense when you look at the file.
Location: `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R`.

**[P2] Bundle + drop delta** — the delta is staged but chat C ran out
of turns before packaging the tarball. After [P1], create a tarball of
`/home/claude/work/delta/` for Quentin to unpack on baseline.

**[P3] Re-run scanner after [P1]** — expect uncovered count to drop
from 19 to 18.

### Medium priority

**[P4] 3 remaining schemas missing block_type**:
`boundary.schema.json`, `boundary_refined.schema.json`,
`gene_conversion_tracts.schema.json`. Add `block_type` field to each.
Then re-enumerate — may reveal even more hidden-from-indexer schemas.

**[P5] band_composition's NA field** — extend `STEP_C01i_decompose.R`
with the instant_q join (~30 lines) so
`q1_decompose_ancestry_div_hom_ref_vs_hom_inv` stops being NA. Sketch
in the archive README at
`_archive_superseded/bk_schemas_pre_canonical/band_composition_folded_into_internal_dynamics_2026-04-24/README.md`.
Load the per-window local-Q TSV for the candidate's chromosome from
`${LOCAL_Q_DIR}/K*/${chr}.local_Q_samples.tsv.gz`, filter to windows
overlapping the interval, average Q per sample, group by
`class_assignment`, average per class → 3 mean Q vectors, L1 distance
of HOM_REF vs HOM_INV.

**[P6] Morphology schema's 4 NA-required fields** in
`utils/registry_key_helpers.R` (TODOs §A1-§A4): `aspect_ratio`,
`architecture`, `sample_grouping`, `spatial_consistency`. Block
`validation_status=incomplete` until architecture classifier written.

### Long-horizon (need HPC or design decisions)

**Bucket B blocked schemas** still need producer→consumer join work:

- `block_detect`: PHASE_01C needs 5 additional columns in
  `block_registry_<chr>.tsv.gz` (confidence_score, inner_boundary_*,
  mean_outside_similarity, contrast_ratio); C01d then joins via
  interval_id=block_id.
- `existence_layer_b`: need `STEP_C00b_attribute_sv_per_candidate.R`
  reading sv_prior_<chr>.rds × candidate regions.
- `existence_layer_c`: symmetric C04b per-candidate attribution step.

**Verification plan on LANTA** (from chat B handoff, not run yet):
sanity Rscript for 10 helpers, LG28 chain C01d→C01g→C01i_d_seal→C01f,
spot-check keys.tsv row count (~80 target up from ~50 baseline).

### Other chat B items still pending (not touched this session)

- §2 clip_count semantics
- §3 depth_anomaly mapping
- §5 fdr_q_value post-pass location
- §6 q6_family_linkage collision (drop from hypothesis_verdict)
- §7 `{side}` loader bug proper fix
- §8 boundary_refined still on left/right

## Naming convention established this session

When a key's value can be sourced by more than one script, put the
producing-script's name (or short-name) as a prefix in the key itself.
This stays visible in keys.tsv without cross-referencing the schema.

Examples:
- `q2_pca_k_used` — "q2 axis, from decompose's PCA clustering"
- `q1_decompose_k_used` — "q1 axis, from decompose" (dual-write to above)
- `q1_decompose_ancestry_div_hom_ref_vs_hom_inv` — "q1 axis, from
  decompose"; future sibling `q1_unified_ancestry_div_*` would come
  from a different script without collision

Apply retroactively if you find another key with the same ambiguity
(e.g. if C01j and C01l both end up producing an `n_segments` value
for different notions of segment, prefix them by script).

## Delta files staged at /home/claude/work/delta/

**Modified** (vs baseline tarball):
```
docs/ENGINES.md
docs/TRIAGE_LEDGER.md
inversion_modules/check_deployment_v10.py
inversion_modules/phase_2_discovery/2c_precomp/PHASE_01C_block_detect.R
inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R
inversion_modules/phase_4_catalog/STEP_C01j_regime_compatibility_engine.R
inversion_modules/phase_4_catalog/STEP_C01l_local_structure_segments.R
inversion_modules/phase_4_catalog/STEP_C01m_distance_concordance.R
inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R
inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py
inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R
inversion_modules/phase_8_evidence_biology/q5_age_and_origin/cross_species_bridge_v6.py
inversion_modules/phase_9_classification/tests/test_phase4b_integration.py
registries/schemas/structured_block_schemas/distance_concordance.schema.json
registries/schemas/structured_block_schemas/internal_dynamics.schema.json
registries/schemas/structured_block_schemas/local_structure_segments.schema.json
registries/schemas/structured_block_schemas/regime_sample_dag.schema.json
registries/schemas/structured_block_schemas/regime_segments.schema.json
registries/schemas/structured_block_schemas/synteny_v6.schema.json
tools/_code_field_check.py
unified_ancestry/00_ancestry_config.sh
utils/registry_key_helpers.R
```

**Created**:
```
docs/NESTED_VS_COMPOSITE.md
docs/SCRIPT_CONTRACT.md
tools/scan_script_contracts.py
unified_ancestry/engines/nested_composition/internal_ancestry_composition.py   (renamed)
_archive_superseded/bk_schemas_pre_canonical/flank_coherence_folded_into_synteny_v6_2026-04-24/README.md
_archive_superseded/bk_schemas_pre_canonical/flank_coherence_folded_into_synteny_v6_2026-04-24/flank_coherence.schema.json
_archive_superseded/bk_schemas_pre_canonical/band_composition_folded_into_internal_dynamics_2026-04-24/README.md
_archive_superseded/bk_schemas_pre_canonical/band_composition_folded_into_internal_dynamics_2026-04-24/band_composition.schema.json
```

**Deleted** (DELETIONS.txt):
```
inversion_modules/phase_7_karyotype_groups/proposal/nested_composition_core.py
inversion_modules/phase_2_discovery/_archive/analysis/nested_composition.py
unified_ancestry/engines/nested_composition/nested_composition.py
registries/schemas/structured_block_schemas/flank_coherence.schema.json
registries/schemas/structured_block_schemas/frequency.schema.json
registries/schemas/structured_block_schemas/frequency.v2.schema.json
```

## Architectural patterns set this session (reusable for future passes)

1. **Fold-into-existing-schema** when two schemas share a writer.
   Archive standalone with README documenting rationale + computation
   sketch.
2. **Provenance-prefixed key names** when same field serves multiple
   q-axes (`q1_decompose_*` vs `q2_pca_*`). Prefix = script short-name.
3. **Dual-write shims** during migrations. Controlled by flags where
   possible; flip off when grep finds no consumers of legacy names.
4. **Contract header + scanner enforcement** for every
   registry-participating script. BLOCKED statuses document why wiring
   is impossible with precise architectural reasons, not handwave.
5. **Three-level schema indexing fallback** (`block_type` → `title` →
   filename) so authoring inconsistencies don't silently hide schemas
   from the scanner.

## Don't re-do

- Reading past chats about band_composition / flank_coherence — both
  folded now, decisions captured in archive READMEs.
- The nested_composition engine/wrapper split — done and verified. If
  anything looks broken, read `docs/NESTED_VS_COMPOSITE.md` first
  (answers 90% of the likely confusion).
- The naming pattern for the band_composition fold — Quentin agreed on
  `q1_decompose_*` prefix. Don't second-guess.

## Contact / context

- Quentin Andres, PhD researcher Kasetsart University Bangkok.
- LG28 inversion: 15.115-18.005 Mb, 60/106/60 karyotype,
  Fst_Hom1_Hom2=0.308 (same as before).
- K clusters = hatchery broodline structure, NOT hybrid/admixture.
- Terminology: "arrangement-discordant IBS tracts consistent with gene
  conversion" (not "gene conversion events").
- Working style: terse, direct, expects autonomous scope management.
  Does not want tests re-run after every change; scale test re-runs to
  change scope. Does not want stale archive reading at session start —
  just current handoff.
