# Handoff prompt for chat 11 — FINAL (v3)

**SUPERSEDES** all prior chat-11 handoff drafts. This version reflects
Quentin's direction after chat 10:

- Storage + parse cost for flat keys is negligible (whole-genome
  registries are a few MB). **Do NOT prune the 367-key spec. Keep
  aggressive flat-key registration as-is.**
- Evidence registry is mature; the library's `reg$evidence$*` methods
  are solid. Phase 3 is wired (chat 5 FIX 29 v2) and is the template.
- **Sample registry and interval registry are underbuilt.** Chat 11's
  main job is to beef them up with methods that nested inversions,
  overlapping inversions, recombinants, and Engine B queries will
  actually need downstream.
- After chat 11, upstream writers (precomp, SV priors, block detect,
  GHSL, C01d, C01g) get wired one per chat. Then phase 2e staircase
  work starts. Then the HPC runs happen once, definitively.

## The reality of what exists

The registry library is at the project root: `registries/` is a sibling
of `inversion_modules/`, `unified_ancestry/`, `utils/`, `Modules/`.

```
registries/
├── api/
│   ├── R/registry_loader.R          567 lines  — 3 classes + query + compute + shortcuts
│   ├── python/registry_loader.py    404 lines  — mirror, Layer 1 atomic only
│   └── bash/registry_loader.sh      112 lines  — SLURM launcher helpers
├── data/
│   ├── sample_registry/              sample_master.tsv, sample_groups.tsv, groups/<id>.txt
│   ├── interval_registry/            windows_master.tsv.gz, candidate_intervals.tsv, cov_registry.tsv
│   └── evidence_registry/            per_candidate/<cid>/{raw,structured}/, global/
└── schemas/structured_block_schemas/  20 Tier-2 JSON schemas, each with keys_extracted directive
```

### Three registries, their current semantics, and where the gaps are

**1. Interval registry — WHERE (spatial, Engine B's backend)**

Current methods: `get_windows`, `get_window_at(chrom, pos)`,
`get_candidate(cid)`, `get_candidate_boundaries(cid)`,
`add_candidate(cid, chrom, start_bp, end_bp, scale, parent_id)`,
`count_candidates`, `get_cov(scope, chrom)`.

**Gaps that hurt nested + overlapping inversions:**
- `parent_id` is a schema column but no method traverses it
- No `get_children(cid)` / `get_parent(cid)` / `get_ancestors(cid)` /
  `get_descendants(cid)`
- No `get_overlapping(chrom, start, end, exclude_cid = NULL)` —
  can't answer "what other candidates overlap this region?" without
  a TSV scan every time
- No `classify_relationship(cid1, cid2)` returning
  nested / partial_overlap / disjoint / equal
- No `get_nested_within(cid)` (all candidates whose interval is
  strictly contained in cid's interval)
- `get_cov` is thin — no bulk lookup of cov paths for a group of
  candidates

These matter because C01d already sets `pattern = nested_fixed /
nested_rare / complex_system` (see L387–394) and `n_children` on
candidates with nested structure, but nothing ties those to the
`parent_id` column in the interval registry. The structural
relationship between A (big) and B (nested inside A) is implicit in
coordinates only.

**2. Sample registry — WHO (sample groups, pedigree)**

Current methods: `get_master`, `get_group(gid)`, `add_group(...)`,
`has_group(gid)`, `list_groups`, `count_groups`,
`get_carriers(cid, status)`.

**Gaps that hurt recombinants + family linkage + cross-inversion
analysis:**
- No `get_sample_metadata(sid)` returning sample_master row
  (family, chrom_parent_origin, coverage, etc.)
- No `get_sample_groups(sid)` — reverse lookup, "which groups does
  sample X appear in?" Needed for: "sample CGA045 is HOM_INV for
  LG12_17 AND RECOMBINANT for LG17_42, is it the same event?"
- No `get_family(sid)` / `list_families()` / `get_family_members(fam_id)`
  — critical for jackknife. Currently C01f reads family info from
  a flat pedigree file, not from the registry.
- No `list_carriers_across_candidates(status = "HOM_INV")` — "which
  samples are HOM_INV for any polymorphic inversion?"
- No `compute_group_overlap(gid1, gid2)` — returns Jaccard / overlap
  count. For overlapping inversions A and B, this answers "do the
  HOM_INV groups of A and B share samples?" which tests whether
  they're really the same haplotype.
- No `list_recombinant_candidates(sid)` — "for sample X, which
  candidates flag it as RECOMBINANT?"
- No bulk `get_groups_for_candidate(cid)` returning all 4–6
  per-candidate groups at once (HOM_REF / HET / HOM_INV / RECOMBINANT
  / RECOMBINANT_GC / RECOMBINANT_DCO)

**3. Evidence registry — WHAT (keys and structured blocks)**

Current methods: `write_block(cid, block_type, data)`,
`read_block(cid, block_type)`, `get_keys(cid, key)`,
`add_evidence(cid, key, value, ...)`, `has_evidence(cid, key)`,
`get_evidence(cid, key)`, `register_candidate(...)`,
`list_candidates`, `count_candidates`.

**Mostly solid.** Phase 3 STEP03 is the reference writer (chat 5 FIX
29 v2). The `keys_extracted` auto-materialization works, schema
validation works, v9 fallback works.

**Minor gap that still matters:**
- No bulk reader: `get_all_blocks(cid)` or `get_blocks_by_type(block_type)`
  — for cross-candidate aggregation (e.g. Q7B dropout ranking across
  all candidates, FDR correction globally)
- No query-by-key-value: `find_candidates_where("q6_group_validation",
  "VALIDATED")` — useful for Phase 4d's tier-gated dispatch

## The design problem for nested + overlapping + recombinants

Quentin's concern (verbatim-ish): we need to understand everywhere
groups get used, especially for recombinants and decomposition of
sample groups, for complex/nested inversions and when inversions
overlap.

### What the 4b scripts do today (audited cover-to-cover in chat 11
prep; this is a summary)

**STEP_C01i_decompose.R (4b.1)**
- Input: candidate_id, chrom, start_bp, end_bp, tier ≤ 3
- Method: seeded k-means on mean PC1 (seeded by SV priors if ≥10
  genotyped anchors), else random init with cheat2 het-DEL
  constraint (forces HOM_INV → HET for samples with het-DEL at bp)
- Output: `internal_dynamics` Tier-2 block with per_sample classes
  + per_window_class.rds (too big for JSON)
- Does NOT register groups yet
- Flags: `flashlight_mode` (seeded/random), `discordant_samples`
  (SV ≠ clustering), `cheat2_constrained_samples`

**STEP_C01i_b_multi_recomb.R (4b.2)**
- Input: per_window_class.rds from 4b.1
- Method: three independent signals combined:
  - S1: per-window class consistency (recombinant = <70% consistent)
  - S2: Clair3 phase switches within interval
  - S3: flashlight hemizygous segment with discordant genotype
  - A sample is recombinant if S1 AND (S2 OR S3), or S3 alone
- cheat24 classifies events: gene_conversion / double_crossover /
  suspicious / ambiguous
- Output: `recombinant_map` Tier-2 block with per-sample switchpoints,
  mosaic lengths, event_class, posterior

**STEP_C01i_c_nested_composition.py (4b.3)**
- Input: Engine B local Q cache (per-window × per-sample Q matrix)
- Method: for each candidate, classify internal structure:
  - homogeneous (1 ancestry throughout)
  - dominant_plus_secondary (>80% one, <20% other)
  - two_block_composite (2 ancestries, balanced) ← nested inversion signature
  - continuous_gradient (smooth transition)
  - multi_block_fragmented (≥40% switches) ← often artifact or
    complex rearrangement
  - diffuse_mixed (many ancestries, no dominant)
- Output: `internal_ancestry_composition` Tier-2 block with
  composite_flag (likely_composite / maybe_composite / not_composite /
  unknown_no_engine_b)

**STEP_C01i_d_seal.R (4b.4)**
- Input: reads all three prior blocks
- Rules applied per sample:
  - Rule 1: recomb detected by 4b.2 → FINAL = RECOMBINANT
  - Rule 2: structure = multi_block_fragmented but 4b.2 says no
    recomb → FINAL = RECOMBINANT (weak, event_class = ambiguous)
  - Rule 3: structure = two_block_composite → keep pca_class, set
    `in_composite_region = TRUE`
  - Rule 4: else → FINAL = pca_class
- Registers groups in sample_registry: HOM_REF, HET, HOM_INV,
  RECOMBINANT, plus RECOMBINANT_GC / RECOMBINANT_DCO subgroups
  (if cheat24 resolved event class)
- Writes HOM_STD alias for v9 back-compat
- Determines quality flags and promotion_cap:
  - composite_flag = likely_composite → promotion_cap = UNCERTAIN
    (C01f cannot promote above this)
  - silhouette < 0.25 → quality_flag = "low_silhouette"
  - bic_gap < 0.05 → quality_flag = "weak_k3_fit"
  - phase_concordance < 0.30 → quality_flag = "low_phase"
- Writes flat keys: q6_group_validation=UNCERTAIN,
  q2_decomp_quality_flags, q6_validation_promotion_cap (if capped),
  q6_family_linkage="unknown" (placeholder — C01f overwrites),
  q6_polymorphism_class="unclassified" (placeholder — C01f overwrites)

### The nested-inversion lifecycle as it exists today

Consider A (big outer inversion) and B (nested inside A):

1. Phase 2 discovery produces A and B as separate candidates
2. C01d sets `pattern = nested_fixed` on A, `n_children = 1`; B has
   `pattern = strong_inversion` (B is a normal-looking inversion; it
   just happens to sit inside A)
3. `parent_id` in `candidate_intervals.tsv` is NOT populated — the
   relationship is encoded only in coordinates (B's interval is
   contained in A's)
4. 4b runs decompose + multi_recomb + nested_comp + seal on A and B
   independently
5. A's nested_composition likely returns `two_block_composite`
   (inside-B vs. outside-B within A, two distinct ancestries)
6. A's seal sees composite_flag=likely_composite → sets
   promotion_cap=UNCERTAIN on A. Correct behaviour. A will stay at
   UNCERTAIN even if C01f finds strong Layer D evidence.
7. B's seal sees composite_flag=not_composite or homogeneous → no
   cap. B can be promoted to VALIDATED normally.
8. Sample CGA045 is HOM_INV for A and HOM_REF for B → appears in
   both `inv_A_HOM_INV` and `inv_B_HOM_REF`. Correct, independent
   groups.
9. C01f hypothesis test may verdict A as H3_nested_composite (or
   not — depends on the T4/T5/T6 composition tests)
10. No explicit link between A and B at any registry level.

**What's missing:**
- `reg$intervals$get_children("A") → ["B"]` — would help phase 4e
  characterize A as "nested parent of B" instead of just "composite"
- Ability to propagate RECOMBINANT detection across parent/child:
  if a sample switches class precisely at B's boundary within A,
  that's more likely a real recombination event than noise
- `reg$samples$compute_group_overlap("inv_A_HOM_INV", "inv_B_HOM_INV")`
  — for cases where two inversions co-segregate perfectly (haplotype
  block containing both). Currently impossible to detect without
  hand-written set intersection code.

### The overlapping (non-nested) case

A and B partially overlap but neither contains the other. Happens
rarely but it does happen — e.g. two independent recombination events
producing overlapping but distinct inversions.

Current state: both scored independently, both decomposed
independently, both get groups registered. Nothing says they're
related. If C01f hypothesis H3 fires on both, it's unclear which is
"real" — could be the same event called twice at slightly different
coordinates, or two real overlapping events.

**What's missing:**
- `reg$intervals$get_overlapping(chrom, start, end, exclude = cid)`
  — "what other candidates overlap this region?" Answer lets Phase 4c
  detect competing hypotheses.
- A new `overlap_reconciliation` step (not in chat 11 scope, but the
  sample-registry method to detect it is)

## The chat-11 job

**Build out the sample registry and interval registry APIs with the
methods listed as "gaps" above.** Do NOT touch the evidence registry
(it's mature and the Phase 3 wiring is the reference). Do NOT wire
upstream writers (that's chats 12–17). Do NOT prune the spec
(Quentin: storage is nothing).

### Scope for chat 11

**Part A — Interval registry extension**

Add these methods to `registries/api/R/registry_loader.R` inside
`load_intervals_api()`:

- `get_children(cid)` — returns character vector of child CIDs (those
  with parent_id == cid)
- `get_parent(cid)` — returns parent CID or NA
- `get_ancestors(cid)` — walks up the parent_id chain; returns vector
  (closest first)
- `get_descendants(cid, depth = Inf)` — recursive down
- `get_overlapping(chrom, start_bp, end_bp, exclude_cid = NULL)` —
  returns CIDs whose intervals overlap the range
- `get_nested_within(cid)` — returns CIDs whose intervals are strictly
  contained in cid's interval (synonym-ish of descendants but
  coordinate-based, doesn't require parent_id to be populated)
- `classify_relationship(cid1, cid2)` — returns
  "equal" / "nested_1_in_2" / "nested_2_in_1" / "partial_overlap" /
  "disjoint"
- `update_candidate(cid, ...)` — partial update (e.g. setting parent_id
  post-hoc when C01d discovers a nesting relationship in pass-2)
- `bulk_add_candidates(dt)` — batch insert from a data.table; much
  faster than row-by-row for large catalogs

Mirror all of these in `registries/api/python/registry_loader.py`
(Python is Layer-1 atomic per chat-5 convention — it just needs to
read, not query/compute).

**Part B — Sample registry extension**

Add these methods to `load_samples_api()` (the shim fallback list
AND the real path that sources `utils/sample_registry.R`):

- `get_sample_metadata(sid)` — returns sample_master row as a list
- `get_sample_groups(sid)` — reverse lookup, all groups containing sid
- `get_family(sid)` — returns family ID or NA (requires sample_master
  to have a `family_id` column; add it if missing)
- `list_families()` — distinct family IDs
- `get_family_members(fam_id)` — samples in a family
- `list_carriers_across_candidates(status = "HOM_INV")` — returns
  unique sample IDs that are `status` for any candidate
- `compute_group_overlap(gid1, gid2)` — returns list(
  intersection = int, union = int, jaccard = numeric,
  only_in_1 = char, only_in_2 = char)
- `list_recombinant_candidates(sid)` — returns CIDs where sid is
  RECOMBINANT
- `get_groups_for_candidate(cid)` — returns named list of all groups
  for cid (HOM_REF / HET / HOM_INV / RECOMBINANT / RECOMBINANT_GC /
  RECOMBINANT_DCO) in one call
- `find_co_segregating_groups(min_jaccard = 0.9)` — across all
  registered groups, find pairs whose membership overlap is ≥
  threshold. Useful for flagging inversions that may be the same
  event called twice.

Python: same methods, atomic read-only (`compute_group_overlap` and
`find_co_segregating_groups` can live in R only if Python implementation
is annoying).

**Part C — Query API extensions**

Add composite methods to `load_query_api()` that use the new
interval + sample methods:

- `nested_family_tree(root_cid)` — returns a tree structure:
  `list(cid, parent, children = list(list(cid, parent, children=...)))`
  Useful for phase 4e characterization of complex nested systems.
- `overlap_clusters(chrom)` — groups overlapping candidates into
  connected components; each component is a set of candidates that
  share samples. For resolving competing hypotheses.
- `sample_inversion_load(sid)` — how many candidates is this sample
  a carrier for? Broken down by status (HOM_INV / HET / RECOMBINANT).
  Useful for Q6 lineage analysis.

### Definition of done for chat 11

- `registries/api/R/registry_loader.R` grows with the new methods
  (probably ~200–300 new lines, pushing file to ~800+ lines)
- `registries/api/python/registry_loader.py` mirrors the atomic ones
  (~100–150 new lines)
- `registries/tests/test_sample_registry_extensions.R` — test fixtures
  for the sample registry additions (synthetic groups, reverse
  lookups, family queries, group overlaps)
- `registries/tests/test_interval_registry_extensions.R` — test
  fixtures for nested/overlapping candidates, relationship
  classification, tree traversal
- `registries/API_CHEATSHEET.md` updated or created: a one-page ref
  for all methods (current + new), grouped by registry, with the
  typical caller (which script calls which method)

### Files chat 11 should NOT touch

- `registries/api/*/registry_loader.*` — evidence API class
  (`load_evidence_api`). It's fine. Don't touch.
- `registries/schemas/structured_block_schemas/*.schema.json` —
  schemas are the frozen contract. Don't touch.
- Any phase_4_postprocessing scripts. All 4a/4b/4c/4d/4e work happened
  in chats 5–10; chat 11 is registry library work.
- Anything in phase_2_discovery/. Upstream writers get wired in
  chats 12–17.

## Upstream writers that chats 12–17 will wire

**All at actual paths in the tarball:**

| Chat | Script | Path | Writes block |
|---|---|---|---|
| 12 | precompute | `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R` | existence_layer_a (skeleton — mostly C01d's job; precomp writes precompute_diagnostic or similar) |
| 12 | SV prior | `inversion_modules/phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R` | existence_layer_b |
| 12 | block detect | `inversion_modules/phase_2_discovery/2c_precomp/PHASE_01C_block_detect.R` | block_detect (schema already exists) + populates interval_registry windows + candidates |
| 13 | seeded regions | `inversion_modules/phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R` | possibly a `seeded_regions` block; registers child candidates via `add_candidate(parent_id=...)` |
| 14 | C01d scoring | `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | existence_layer_a (the main layer-A block) |
| 15 | C01g boundaries | `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R` | boundary_left / boundary_right. Requires creating `utils/registry_key_helpers.R` (currently missing — closes chat-7 Finding 9) |
| 16 | GHSL | `inversion_modules/phase_2_discovery/2e_ghsl/STEP_C04_snake3_ghsl_v5.R` | existence_layer_c |
| 17 | phase 2e staircase | (Quentin's actual end-goal) | — |

## Reading order for chat 11 (≈45 min)

**Primary — understand what's there:**

1. `registries/api/R/registry_loader.R` end-to-end (567 lines). Focus
   on `load_samples_api` (L143–187) and `load_intervals_api` (L192–260)
   — these are what you'll extend.
2. `registries/api/python/registry_loader.py` (404 lines) — focus on
   the same classes in Python.
3. `registries/api/bash/registry_loader.sh` (112 lines) — to confirm
   nothing breaks in the bash helpers when R API grows.

**Secondary — understand the callers, so you know what methods to
expose:**

4. `inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_d_seal.R`
   (422 lines) — this is the biggest sample-registry writer. Read the
   `register_all_groups()` function (L221–295) to confirm the groups
   you're indexing match what it registers.
5. `inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_decompose.R`
   L310–448 — to see what `internal_dynamics` block writes, which
   seal then reads.
6. `inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_b_multi_recomb.R`
   L350–437 — the recombinant output structure.
7. `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
   L380–450 — the `pattern` classification (nested_fixed / nested_rare /
   complex_system) that you'll use when implementing `get_children` +
   relationship classification.
8. `inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R`
   — scan for `samples$get_carriers` calls (grep is enough); these
   are the biggest sample-API consumers and their patterns tell you
   what convenience methods would save them code.

**Tertiary — schema context:**

9. `registries/schemas/structured_block_schemas/internal_dynamics.schema.json`
10. `registries/schemas/structured_block_schemas/recombinant_map.schema.json`
11. `registries/schemas/structured_block_schemas/internal_ancestry_composition.schema.json`
12. `registries/schemas/structured_block_schemas/hierarchy.schema.json` or
    `block_index.schema.json` if they exist (check the 20-schema
    directory; a hierarchy block MAY need to be added, but that's a
    chat-12+ decision, not chat-11 scope).

**Skip entirely:**
- Phase 4e scripts (done in chat 10; don't touch)
- Phase 3 scripts (wired in chat 5)
- Evidence API in registry_loader.R (mature; don't touch)
- Schemas (frozen contracts; don't touch)

## Cumulative findings carry-forward

Still open from prior chats:
- Finding T (chat 8): README structure_type drift — doc pass
- Finding U (chat 8): `q2_decomp_quality_flags` forward-declared,
  no consumer — acceptable since seal writes it and nobody currently
  reads it; chat 11 may add a reader if it's useful for the new sample
  API methods (unlikely)
- Finding V (chat 8): cheat24 threshold divergence — design call
  needed
- Finding W (chat 9): `register_C01f_keys` helper verification —
  deferred until the helper is actually written (chat 15 with C01g)
- Finding Y (chat 9): **resolved** by chat 10 FIX 50
- Finding Z (chat 9): 4d README stale — doc pass
- Finding AA (chat 9): `population_regenotype.py` misfiled — defer
- Finding AB (chat 9): 4e key-count cosmetic — defer
- Finding AF (chat 9): **resolved** by chat 10 FIX 45
- Finding AG (chat 9): **resolved** by chat 10
- Finding 8 (chat 7): 4a writes no registry blocks — chat 14
- Finding 9 (chat 7): C01g silent skip on missing helper — chat 15
- Finding AH (chat 10): driver runtime-untested on real data —
  unblocked by chat 17's HPC run
- Finding AI (chat 10): example launcher coexists with prod — cosmetic

Chat 11 may surface new findings while auditing the sample and interval
APIs (e.g. bugs in the existing `get_group` fallback, or missing
fields in `sample_master.tsv`). Flag them; don't fix unless it blocks
the new-method work.

## Parse-check backlog

Same cumulative list as chat 10 (14 items). Chat 11 will add:
- Modified `registry_loader.R` (parse + smoke test)
- Modified `registry_loader.py` (AST + smoke test)
- Two new test fixtures

## Posture

Understand-before-code (chat 9) + smoke-test-after-code (chat 10).
The smoke test for chat 11 is the test fixtures — they exercise every
new method against a synthetic registry with nested candidates,
overlapping candidates, and recombinant samples.

**No three-fix ceiling this chat.** Sample + interval registry API
work is architectural; every new method is defensible against the
"will this be needed downstream for nested/overlap/recomb analysis?"
question, and there are ~15–20 new methods to add total.

**Don't get tempted to wire writers.** Chat 11 is pure library work.
No changes to any script outside `registries/`. If reading C01d or
C01i surfaces a bug, flag it — don't fix.

## End of chat 11

- `AUDIT_LOG_chat11_<date>.md` — document the new methods + any
  findings surfaced
- Append to `SESSION_SUMMARY_<date>.md` and
  `FIXES_APPLIED_<date>.md`
- New `registries/API_CHEATSHEET.md` or update existing
- Repack as `inversion-popgen-toolkit_chat11_registry_extensions_<date>.tar`

## The revised ~6-chat plan

| Chat | Scope | End state |
|---|---|---|
| 11 | Sample + interval registry API buildout | ~15–20 new methods, test fixtures, cheatsheet |
| 12 | Wire precomp + SV prior + block detect | C00, C01a, PHASE_01C all call registry |
| 13 | Wire C01b seeded regions (registers child candidates with parent_id) | Nested relationships materialize in interval registry |
| 14 | Wire C01d scoring (existence_layer_a block) | Layer A evidence flows to registry |
| 15 | Wire C01g boundaries + create `utils/registry_key_helpers.R` | Boundaries flow to registry; Finding 9 closed |
| 16 | Wire GHSL (existence_layer_c) | Layer C evidence flows to registry |
| **17** | **Phase 2e staircase** | Upstream all registry-correct; HPC runs can happen |

After chat 16, all upstream modules are registry-wired. Precomp, SV
priors, block detect, C01b, C01d, C01g, GHSL all run on 28
chromosomes definitively. A day or two of HPC time. They never run
again because their contracts are correct.

Chats 18+: return to downstream. Phase 4e real-data smoke against the
populated registry (the chat-10 driver gets its first real-data
invocation). Q7B breakpoint_evidence_audit wiring. Phase 4a
characterization pass. Eventually: the `final_catalog.tsv` with 14-axis
classification per candidate. Then manuscript work begins.
