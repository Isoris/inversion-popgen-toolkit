# Toolkit audit — where to merge, combine, simplify

**Scope:** the whole `inversion-popgen-toolkit` as of the Git clone
on this HPC-offline weekend (2026-04-23). Audit is structural (reorganise
to reduce surface area) not logical (keep all logic intact).

**Principle:** preserve every bit of logic. Reduce file count, schema
count, launcher count, and name confusion.

---

## Merges that pay off

### Merge 1 — `cheat27` + `cheat14` → `sd_substrate.R`

**Current state**
- `cheat14_minimap2_inverted_sd.R` (171 lines) — minimap2 self-alignment
  of the boundary region, detects inverted-SD pairs.
- `cheat27_sd_nahr_substrate.R` (219 lines) — loads BISER2 pre-computed
  SD catalog, detects inverted-SD pairs at boundaries, checks
  concordance with cheat14 minimap2 output.

Both measure **the same biological entity** (inverted repeat pair
straddling an inversion boundary) in the same reference genome. They
are two tools detecting SDs; the "concordance" check inside cheat27 is
two tools agreeing on the same feature, which is NOT two independent
lines of evidence.

**Proposed merge**
- Single script `sd_substrate.R` that runs both detectors (minimap2
  + BISER2 if available), outputs a unified block:

```
q4_sd_substrate_class ∈ {
  inverted_both_methods_agree,
  inverted_biser2_only,
  inverted_minimap2_only,
  direct_only,
  none,
  untested_no_biser2
}
q4_sd_best_identity_pct   float (whichever method found the stronger pair)
q4_sd_best_length_bp      int
q4_sd_methods_used        list
```

**Benefit**
- ~400 → ~250 lines (one tool instead of two + concordance glue)
- One schema replaces two
- Eliminates the overclaim risk that cheat27's `concordance_mechanism
  = "NAHR_confirmed"` generates — the new key honestly says
  "two tools detected the same SD", not "mechanism confirmed".
- Downstream (characterize_q4) reads one block, not two.

**Risk**
- Low. Both current scripts have separate schemas consumed only by
  characterize_q4; migration is a schema rename + one reader change.

---

### Merge 2 — `cheat29` + new `cheat29b_assembled_junction.py`
                  → `junction_forensics.py` (dispatcher)

**Current state (after v5)**
- `cheat29_junction_forensics.R` (347 lines) — reads REFERENCE sequence
  at estimated boundary, detects MH / TSD / TE.
- `cheat29b_assembled_junction.py` (v5, 512 lines) — reads
  CONSENSUS/HOMLEN from PRECISE SV records, actual junction forensics.

Currently v5 runs both and emits two separate blocks (`mechanism.json`
from cheat29 + `mechanism_assembled.json` from cheat29b). The
synthesis script then reads both.

**Proposed merge**
- Dispatcher script `junction_forensics.py` that:
  1. First tries cheat29b logic (PRECISE assembled junction).
  2. If not available, falls back to cheat29 reference-based logic.
  3. Writes ONE `mechanism.json` block with a `source` field:
     `source ∈ {assembled_delly_inv, assembled_delly_bnd_pair,
     assembled_delly_bnd_single, assembled_manta_inv, reference_fallback,
     none}`

**Keys**
```
q4_junction_source       from above enum
q4_junction_class        blunt / MH_short / MH_long / TE / insertion /
                          unresolved
q4_junction_homlen       int
q4_junction_class_confidence  high (assembly) / low (reference) / none
```

**Benefit**
- Eliminates dual-schema confusion.
- Forces the correct precedence (assembly beats reference) at
  measurement time, not at synthesis.
- The `_confidence` key propagates the "reference-based is weaker"
  signal directly, so every downstream consumer sees it without having
  to re-derive the distinction.
- ~850 → ~500 lines net.

**Risk**
- Medium. The reference-fallback logic is R; the assembled logic is
  Python. Either port cheat29 to Python (recommended; Python is
  already in the pipeline), or keep R side as a helper called by the
  dispatcher. Python port is ~150 lines.

---

### Merge 3 — `STEP_C01l_local_structure_segments.R`
                 + `STEP_C01m_distance_concordance.R`
                 → `STEP_C01lm_local_structure.R`

**Current state (Phase 4a)**
- Both consume the boundary catalog output of C01g.
- Both produce summary statistics per candidate.
- Both have separate LAUNCH_*.sh scripts (not in the dir listing I saw,
  but implied by the file layout).

**Proposed merge**
- One script with `--mode {segments, distances, all}`.
- Launcher remains a wrapper for each mode.

**Benefit**
- Saves ~200 lines of duplicated I/O and argument-parsing boilerplate.
- One set of tests instead of two.

**Risk**
- Very low. Pure refactor, no logic change.

---

### Merge 4 — five `existence_layer_*.schema.json` into three

**Current state** (in `registries/schemas/structured_block_schemas/`)
- `existence_layer_a.schema.json`
- `existence_layer_b.schema.json`
- `existence_layer_b_bnd_rescue.schema.json`
- `existence_layer_c.schema.json`
- `existence_layer_d.schema.json`

**Proposed merge**
- Keep `_a`, `_c`, `_d` as is.
- Merge `_b` and `_b_bnd_rescue` into one schema with a
  `rescue_source` field. Readers already treat them as the same logical
  layer (Layer B = all SV geometry evidence).

**Benefit**
- 5 → 4 schemas.
- Layer B consumption code becomes one `evidence.read_block("existence_layer_b")`
  call instead of two.

**Risk**
- Low. Migration is a one-liner in each consumer.

---

## Wiring gaps to close (logic exists, keys not exposed)

### Gap 1 — `breakpoint_evidence_audit.py` Q7B keys not read by 4e
- The audit writes `q7b_pca_carriers_strong_sv` etc. to the registry.
- `compute_candidate_status.R::build_key_spec()` doesn't list them as
  wired → they don't contribute to Axis 2 completion %.
- **Fix:** add 8 q7b keys to `q7_wired` in `build_key_spec()`. ~10 lines.

### Gap 2 — `cheat30` GDS keys not consumed by Q5 age
- cheat30 writes `q5_hom_inv_gds_*` but `characterize_q5()` doesn't
  use them to upgrade age_class from "measured" to "answered".
- **Fix:** add GDS read in characterize_q5 and a `gds_supports_age` gate
  in the convergence rules. ~15 lines.

### Gap 3 — v5/v6 `q_overall_structural_class` not in Axis 5
- v5's `assign_structural_class_v5.py` (and v6's v6 equivalent) write
  `final_label.json` in an outdir parallel to the registry.
- `compute_candidate_status.R` does not read it, so the final label
  isn't promoted to an official axis output.
- **Fix:** add Axis 5 to `compute_candidate_status.R`:
  ```r
  axis_5 <- list(
    final_class = read_final_label(cid, keys_dir)$q_overall_structural_class,
    weakest_component = read_final_label(cid, keys_dir)$weakest_component
  )
  ```
  ~15 lines.

### Gap 4 — `gene_conversion_detector.R` outputs not in Q2
- The GC tract detector in 4b writes `gene_conversion_tracts.json`.
- `characterize_q2()` reads `internal_dynamics.json` but not
  `gene_conversion_tracts.json`.
- **Fix:** add the block to Q2 read list, add `q2_gc_tract_count` and
  `q2_gc_tract_total_bp` to wired keys. ~10 lines.

**Total wiring work:** ~50 lines across 3 R files + 1 Python. One PR.

---

## Simplifications worth doing

### Simplify 1 — Drop the 15 v9↔v10 alias keys

**Current state**
- `build_key_spec()` has 15 intentional alias pairs (e.g. `q6_n_HOM_REF`
  + `q6_n_HOM_STD`) maintained for v9 readers.
- The comment in `SPEC_VS_REALITY.md` states this is for migration.

**Proposal**
- v9 writers are gone. Find the last v9 reader, patch it to v10 names,
  drop the aliases.
- Saves 15 lines in spec + eliminates the confusion of two names for
  one concept.

**Risk**
- Medium. Need to grep for every `q6_n_HOM_STD`, `q7_t9_jackknife_verdict`,
  etc. access and confirm v10 alternative works.
- ~1 session of careful grep + test.

### Simplify 2 — Skip 4b steps b and c when `q1_composite_flag = clean`

**Current state**
- 4b always runs: decompose → multi_recomb → nested_composition → seal.
- For candidates with `q1_composite_flag = clean` AND `q1_n_children = 0`,
  steps b and c do no useful work but still execute.

**Proposal**
- Add a gate at the top of the 4b orchestrator launcher:
  ```bash
  if [ "$Q1_COMPOSITE_FLAG" = "clean" ] && [ "$Q1_N_CHILDREN" = "0" ]; then
      # skip steps b and c, go straight to seal
      bash LAUNCH_C01i_decompose.sh
      bash LAUNCH_C01i_d_seal.sh
  else
      # full pipeline
      ...
  fi
  ```
- Saves 60-80% of compute on simple candidates.

**Risk**
- Low if done at orchestrator level. The seal step handles missing
  b/c inputs gracefully (empty recombinant_map and
  internal_ancestry_composition mean "no recombinants, no nesting").

### Simplify 3 — Rename `cheat*` scripts to descriptive names

**Current state**
- 6 scripts in 4d named `cheat6`, `cheat27`, ..., `cheat30`.
- Names are historical (from the "flashlight cheats" framework).
- Make discovery by new readers (reviewers, collaborators) harder.

**Proposal**
- Rename with descriptive prefix, keep cheat number as suffix for
  audit trail:
  - `cheat6_ancestry_jackknife_v934.R` → `ancestry_jackknife.R`
    (comment: `# (formerly cheat6)`)
  - `cheat27_sd_nahr_substrate.R` → merged into `sd_substrate.R`
  - `cheat28_tandem_repeat_context.R` → `tandem_repeat_context.R`
  - `cheat29_junction_forensics.R` → merged into `junction_forensics.py`
  - `cheat30_gds_by_genotype.R` → `gds_by_genotype.R`

**Risk**
- Low. Two-line git mv + update callers. One launcher update each.

### Simplify 4 — Collapse LAUNCH_ wrappers with shared env

**Current state**
- Each 4b step has its own `LAUNCH_*.sh` that sets env vars + runs the
  R/Python script.
- Env vars are repeated across launchers.

**Proposal**
- One `04b_orchestrator.sh` that sources a shared
  `00_module4b_config.sh` (same pattern as MODULE_4D etc.) and runs all
  steps in order.
- Keep individual LAUNCH_*.sh for dev use; mark them "dev_only".

**Risk**
- Very low.

---

## Organisation — don't change this

The `4a_existence_layers/` / `4b_group_proposal/` / `4c_group_validation/`
/ `4d_group_dependent/` / `4e_final_classification/` directory layout
**is good**. It maps directly to the state-machine architecture
(NONE → UNCERTAIN → SUPPORTED → VALIDATED). Keeping it role-based rather
than question-based makes the pipeline order obvious.

**Resist the urge to reorganise by question** (Q1/, Q2/, etc.). A given
question (e.g. Q7) gets evidence from multiple role-layers (4a Layer A,
4a Layer B, 4c Layer D, etc.). Role-by-role is the natural unit.

---

## Where to put the new v6 code

### New in 4b
- `interior_structure_diagnostic.py` → `4b_group_proposal/`
  - It's a group-dependent analysis (requires Q-groups from decompose).
  - Launched after seal step in the 4b orchestrator.

### New in 4d
- `boundary_region_annotations.R` (v6 Kuang-style — priority F4, not yet
  coded) → `4d_group_dependent/`
- `boundary_internal_alignment.sh` (priority F5, not yet coded) →
  `4d_group_dependent/`
- `cross_species_bridge_v6.py` → `4d_group_dependent/`
  (replaces v5 cross_species_bridge.py)

### New in 4e
- `assign_structural_class_v6.py` → `4e_final_classification/`
  (replaces v5 assign_structural_class_v5.py)
- Axis 5 extension to `compute_candidate_status.R` → `4e_final_classification/`

### New schemas
- `interior_structure.schema.json` → `registries/schemas/structured_block_schemas/`
- `synteny_v6.schema.json` → same
- Deprecate `synteny_dollo.schema.json` (keep the file, mark as
  deprecated in its description; v6 readers use `synteny_v6`).

---

## Implementation order for the audit-driven work

1. **Wiring gaps (Gap 1-4):** trivial, 50 lines, 1 session. No risk.
2. **Merge 4 (schemas):** 1 session, low risk.
3. **Simplify 3 (rename cheats):** 1 short session, low risk.
4. **Simplify 2 (4b skip gate):** 1 session, moderate risk, big speed win.
5. **Merge 1 (SD substrate):** 2 sessions — design the merged schema,
   port logic, migrate consumers. Moderate risk.
6. **Merge 2 (junction dispatcher):** 2 sessions — port cheat29 to
   Python, build dispatcher, migrate. Moderate risk.
7. **Simplify 1 (drop v9 aliases):** 1 careful session, medium risk.
8. **Merge 3 (STEP_C01lm):** 1 session, very low risk.

~10 sessions of audit-driven refactoring. None of it is urgent; the
LG28 manuscript can ship without any of it. Do this when the HPC is up
and the backend is stable, as a "clean-up before public release" phase.

---

## What NOT to merge

- DO NOT merge BISER2 and minimap2 substrate detection **methods**
  together internally. They are genuinely different algorithms and
  having both run independently is a methodological strength. What we
  merge is the **reporting layer** (one block, not two), not the
  detection layer.
- DO NOT merge cheat28 TR context with cheat27/cheat14 SD detection.
  TRs and SDs are biologically distinct substrate classes.
- DO NOT merge 4a Layer B and 4d cheat-based SV audits. The 4a layer
  decides existence; the 4d layer audits a known candidate. Different
  roles.
- DO NOT merge Q4 substrate and Q4 junction. They are separate
  sub-verdicts for a reason (v4/v5 framework argument).

---

## Status update — what this chat (2026-04-23) resolved

Items from the audit that this chat's bundle addresses:

**RESOLVED this chat:**
- *Wiring gap 3 (Axis 5 for `q_overall_structural_class`)* — the v7
  assigner writes `final_label.json`; the HANDOFF §6a wires it into
  `compute_candidate_status.R` as Axis 5.

**DEFERRED until wiring session (addressed in HANDOFF.md §8):**
- Wiring gap 1 — Q7B keys from breakpoint_evidence_audit
- Wiring gap 2 — cheat30 GDS into characterize_q5
- Wiring gap 4 — gene_conversion_detector into characterize_q2

**NOT YET ADDRESSED (future sessions):**
- All merges (1, 2, 3, 4)
- All simplifications (1, 2, 3)
- Schema dedup (merge 4)
- Cheat script renaming (simplify 3)

**NEW item added from this chat:**
- *Merge candidate 5 (DELETED, no action):* v6's
  `interior_structure_diagnostic.py` was designed then superseded
  within the same chat by `bp_pipeline_bridge.py`. The v6 file is NOT
  in this bundle. If an earlier session deployed it, delete from
  `4b_group_proposal/` per HANDOFF §4a.

**New wired keys from this bundle (per HANDOFF §6b):**

| Block source | New wired keys | Count |
|---|---|---|
| `fragment_distribution` | q2 + q3 interior/breakpoint keys | 10 |
| `synteny_v6` | q5 cross-species keys | 6 |
| `mechanism_assembled` | q4b assembled-junction keys | 4 |
| `bnd_sided_support` | q7b BND rescue keys | 3 |

Projected wired-key coverage after this bundle: 271/367 → **292/367**
(73.8% → 79.6%).
