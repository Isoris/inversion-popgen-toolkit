# AUDIT LOG — chat 11 of the phase-4 audit series (2026-04-17) — PART H

Scope: sample and interval registry API buildout per chat-11 handoff. 22 new
methods across R + Python bindings, plus two shippable test fixtures and the
API cheatsheet. Three pre-existing bugs surfaced and flagged; one fixed
in-session because it blocked smoke testing of the new methods. Two orphan
duplicate files deleted.

## Reading done this chat

- `HANDOFF_PROMPT_chat11_2026-04-17.md` — full handoff (priority tasking,
  scope, design problem for nested/overlap/recomb, reading order)
- `AUDIT_LOG_chat10_2026-04-17.md` — chat-10 summary, FIX 45–50
- `registries/api/R/registry_loader.R` (567 lines) — all three load_*_api
  functions, the query + compute layers, the shortcuts
- `registries/api/python/registry_loader.py` (404 lines) — SamplesAPI,
  EvidenceAPI, IntervalsAPI
- `registries/api/bash/registry_loader.sh` (112 lines) — path resolution
  and validation gates
- `utils/sample_registry.R` (254 lines) — the delegated implementation
  backing `load_samples_api`
- `inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_d_seal.R`
  — register_all_groups() to confirm group naming (`inv_<cid>_<STATUS>`)
- `inversion_modules/phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R`
  `comp_from_registry()` — the biggest sample-API consumer pattern
- `inversion_modules/phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
  — pattern classification (nested_fixed/nested_rare/complex_system)
  to confirm the relationship between pattern and parent_id

Total ≈ 1600 lines read before coding. Less than chat 10 because the
handoff was dense and the library structure was clear from the first pass.

## Posture

Handoff-prescribed: understand-before-code + smoke-test-after-code.
No three-fix ceiling (chat-11 is architectural library work). Each new
method defended against "needed for nested/overlap/recomb downstream?"
question before writing it.

## Methods added this session (22 total)

### Interval registry extensions (9 methods)

Added to `load_intervals_api()` in `registries/api/R/registry_loader.R`
and mirrored in `IntervalsAPI` in `registries/api/python/registry_loader.py`.

| Method | Purpose | Downstream caller |
|---|---|---|
| `get_children(cid)` | direct children where `parent_id == cid` | phase 4e nested-system characterization |
| `get_parent(cid)` | direct parent CID or NA | characterization of nested leaves |
| `get_ancestors(cid)` | walk up parent_id chain, closest first, cycle-safe | same |
| `get_descendants(cid, depth)` | BFS down, depth-capped | C01d pass-2 recursive propagation |
| `get_overlapping(chrom, s, e, exclude?)` | CIDs overlapping a range | C01b seeded regions, 4c overlap reconciliation |
| `get_nested_within(cid)` | coord-based strict containment (no parent_id) | chat-13 seeded-region child discovery before parent_id assignment |
| `classify_relationship(cid1, cid2)` | equal / nested_N_in_M / partial_overlap / disjoint | 4c competing-hypothesis detection |
| `update_candidate(cid, ...)` | partial column update (e.g. parent_id post-hoc) | C01d pass-2 nesting discovery |
| `bulk_add_candidates(dt)` | batch insert, dup-skip | chat-13 seeded regions bulk writer |

Internal helpers: `read_cands()` (internal closure), `_has_parent()` (Python).
Existing 7 methods preserved verbatim.

### Sample registry extensions (10 methods)

Added to `load_samples_api()` in R and `SamplesAPI` in Python.

| Method | Purpose | Downstream caller |
|---|---|---|
| `get_sample_metadata(sid)` | sample_master row as list | any script needing per-sample context |
| `get_sample_groups(sid)` | reverse lookup: all groups containing sid | cross-inversion sample tracking |
| `get_family(sid)` | family_id (NA if no column; warn-once) | jackknife, Q6 family linkage |
| `list_families()` | distinct family IDs | same |
| `get_family_members(fam_id)` | samples in a family | same |
| `list_carriers_across_candidates(status)` | union of `inv_*_<status>` members | Q6 polymorphism-class assignment |
| `compute_group_overlap(gid1, gid2)` | intersection/union/jaccard/set-diffs | overlapping-inversion same-event detection |
| `list_recombinant_candidates(sid)` | CIDs where sid is RECOMBINANT | sample-level recombinant tracking |
| `get_groups_for_candidate(cid)` | named list of all 7 status groups in one call | C01f collapse (currently does 4 individual has_group+get_group calls per cid) |
| `find_co_segregating_groups(min_jaccard)` | all pairs with membership overlap ≥ threshold | repeat-called-event detection across whole catalog |

Internal helpers: `build_extended()` (refactored factory shared by real-path
and shim-path), `members_path_for()`, `read_members()`, `check_family_col()`
(warn-once cache via closure in R / class attr in Python).

### Query API composite extensions (3 methods, R only)

Added to `load_query_api()`.

| Method | Purpose | Downstream caller |
|---|---|---|
| `nested_family_tree(root_cid)` | recursive tree structure | phase 4e characterization format |
| `overlap_clusters(chrom)` | connected-component grouping by overlap | 4c competing-hypothesis cluster identification |
| `sample_inversion_load(sid)` | per-status counts for one sample | Q6 lineage + sample-level QC |

Python intentionally NOT mirrored for these three per handoff direction
(they live in R because they compose interval + sample API calls).

## Findings surfaced this session

### Finding AJ (fixed in-session)

**Bug:** In `load_samples_api` (real-path branch), `get_master` was aliased
to `full$list_groups` — returning the groups table instead of the master
table. Pre-existing from the initial v10 wiring.

**Why it blocked chat-11 work:** `get_sample_metadata(sid)` needs to read
`sample_master.tsv` rows. The wrong alias would have silently returned
group rows and produced garbage metadata.

**Fix:** Rewired to `full$get_master` (which the old `sample_registry.R`
already exposes). No signature change to consumers.

### Finding AK (fixed in-session)

**Bug:** In `utils/sample_registry.R::add_group()`, line 154 (post-fix 149)
read `groups_file` via `fread(groups_file)` without `colClasses`. `data.table`
auto-detects the `created` column ("%Y-%m-%d %H:%M:%S" format) as `POSIXct`,
but `format(Sys.time(), ...)` produces a `character`. On the second
`add_group()` call in the same R session, `rbind(g, new_row)` fatals with
"Class attribute on column 8 of item 2 does not match".

**Why it hadn't surfaced before:** Real production runs typically call
`add_group` once per R session. Chat-11 smoke testing hit it immediately
(6 groups registered in one session for the fixture). Chat-13's planned
bulk seeded-region writer would have hit it hard.

**Fix:** `fread(groups_file, colClasses = c(created = "character"))`.
One-line fix. Annotated in-place with 2026-04-17 Finding-AK comment.

**Why this is a legitimate "blocks new-method work" fix (per handoff rule):**
The new `compute_group_overlap`, `find_co_segregating_groups`, and test
fixtures cannot be smoke-tested without registering multiple groups in one
session.

### Finding AL (flagged, not fixed — deferred to chat 12)

**Bug:** Function-name collision between
`registries/api/R/registry_loader.R::load_registry` (the unified entry
point) and `utils/sample_registry.R::load_registry` (the legacy sample
groups entry point). When `load_samples_api` sources the legacy file, the
global-env copy of the legacy `load_registry` overwrites the unified one.
Any subsequent call in the same session to `load_registry()` resolves to
the legacy version (different signature: `registry_dir` vs.
`registries_root, create_if_missing`).

**Why not fixed now:** Fix requires renaming one function and updating all
call sites. The legacy `load_registry` is sourced by `utils/load_bridge.R`
and by any script that sources `sample_registry.R` directly. Cross-codebase
rename is chat-12 scope.

**Proposed fix for chat 12:** Rename
`utils/sample_registry.R::load_registry` to `load_sample_groups_api` (or
rename file + function together). Update `load_bridge.R` step 4 and any
direct consumers.

### Finding AM (resolved in chat 11 tarball)

**Bug:** `inversion_modules/utils/sample_registry.R_v254` is byte-identical
to `utils/sample_registry.R` (md5: `3a3b4d1d502db8df20316fc20bee1b94`).
The `_v254` suffix suggests an informal archive marker. No scripts source
this copy (confirmed via grep across the codebase). Orphan duplicate.

**Fix:** Deleted from chat-11 tarball (`--exclude` during rsync stage).

### Finding AN (resolved in chat 11 tarball)

**Bug:** `inversion_modules/utils/sample_map.R` is byte-identical to
`utils/sample_map.R` (md5: `6ddc06fee8e40bed9b1981a10855138d`). No scripts
source this copy. Orphan duplicate.

**Fix:** Deleted from chat-11 tarball.

## Files changed

**Modified:**

- `registries/api/R/registry_loader.R` (566 → 1166 lines; +600 lines)
- `registries/api/python/registry_loader.py` (403 → 812 lines; +409 lines)
- `utils/sample_registry.R` (254 → 261 lines; +7 lines for Finding AK +
  comments)

**Created:**

- `registries/tests/test_interval_registry_extensions.R` (new dir, new file)
- `registries/tests/test_sample_registry_extensions.R`
- `registries/API_CHEATSHEET.md`

**Deleted:**

- `inversion_modules/utils/sample_registry.R_v254` (Finding AM)
- `inversion_modules/utils/sample_map.R` (Finding AN)

**Untouched (explicit non-scope per handoff):**

- `registries/api/R/registry_loader.R::load_evidence_api` — evidence API
  mature; Phase 3 STEP03 wiring is reference writer
- `registries/schemas/structured_block_schemas/*.schema.json` — frozen
  contracts
- All phase_2_*, phase_3_*, phase_4_* scripts — upstream writer wiring
  deferred to chats 12–17
- `registries/api/bash/registry_loader.sh` — no extensions needed; bash
  is SLURM-gate-only

## Smoke tests

### Interval extensions

Synthetic 6-candidate topology covering every relationship type:

```
LG12_A  [ 1Mb, 10Mb]  (outer)
  └─ LG12_B  [ 3Mb,  6Mb]  parent_id=LG12_A  (nested)
      └─ LG12_C  [ 4Mb,  5Mb]  parent_id=LG12_B  (deep nested)
LG12_D  [ 8Mb, 12Mb]  (partial overlap with A, not nested)
LG12_E  [20Mb, 25Mb]  (disjoint same chrom)
LG17_X  [ 1Mb,  5Mb]  (different chrom)
```

All 9 new methods exercised. All assertions pass. Cycle-safety confirmed
by depth-capped descendants. Dup-skip in `bulk_add_candidates` confirmed.

### Sample extensions

Synthetic 10-sample cohort (CGA001–CGA010), 3 families (F1, F2, F3, plus
one no-family sample), multiple inversions with overlapping membership:

- LG12_A polymorphic: HOM_REF/HET/HOM_INV/RECOMBINANT registered
- LG12_B nested: HOM_INV subset of LG12_A HOM_INV
- LG17_X independent: HOM_INV identical to LG12_A HOM_INV (co-segregating
  test case)

All 10 new sample methods + 3 query composites exercised. Co-segregation
detection confirmed. Family-column-absent graceful fallback confirmed.

### Python mirror

Same topology and cohort. All 19 mirrored methods exercised. Parity
confirmed on: `get_children/parent/ancestors/descendants`,
`get_overlapping`, `get_nested_within`, `classify_relationship`,
`update_candidate`, `bulk_add_candidates`, all 10 sample methods.

## Parse-check backlog

Same 14-item cumulative list as chat 10, plus:

- `registries/api/R/registry_loader.R` — parsed (9 top-level exprs)
- `registries/api/python/registry_loader.py` — AST OK
- `utils/sample_registry.R` — parsed (3 top-level exprs, same as pre-fix)
- `registries/tests/test_interval_registry_extensions.R` — ran end-to-end
- `registries/tests/test_sample_registry_extensions.R` — ran end-to-end

## Cumulative findings carry-forward

Still open from prior chats:
- Finding T (chat 8): README structure_type drift — doc pass
- Finding V (chat 8): cheat24 threshold divergence — design call needed
- Finding W (chat 9): `register_C01f_keys` helper verification — deferred
  to chat 15 (when helper is written)
- Finding Z (chat 9): 4d README stale — doc pass
- Finding AA (chat 9): `population_regenotype.py` misfiled — defer
- Finding AB (chat 9): 4e key-count cosmetic — defer
- Finding 8 (chat 7): 4a writes no registry blocks — chat 14
- Finding 9 (chat 7): C01g silent skip on missing helper — chat 15
- Finding AH (chat 10): driver runtime-untested on real data — chat 17 HPC
- Finding AI (chat 10): example launcher coexists with prod — cosmetic

New this chat:
- Finding AJ: resolved in-session (get_master alias bug)
- Finding AK: resolved in-session (column class drift)
- Finding AL: open, deferred to chat 12 (load_registry shadowing)
- Finding AM: resolved (orphan duplicate deleted)
- Finding AN: resolved (orphan duplicate deleted)

## End-of-chat-11 state

- Registry library: 1166 R lines + 812 Python lines. Three new composite
  query methods in R. Parse-clean. Test-fixture-clean end-to-end.
- Two orphan duplicate files removed.
- Cheatsheet + audit log + test fixtures shipped.
- Next chat (12): wire precomp + SV prior + block detect writers. Also
  resolve Finding AL (load_registry rename).

## Architecture clarified for Quentin

Restating the three-registry architecture in the naming form we'll use
going forward:

- **sample_registry** (data/sample_registry/): WHO — samples, groups,
  pedigree
- **interval_registry** (data/interval_registry/): WHERE — windows,
  candidates, cov paths
- **evidence_registry** (data/evidence_registry/): WHAT — per-candidate
  keys and Tier-2 blocks

Three language bindings (not three APIs): R (full), Python (atomic +
writes), bash (SLURM gates). One on-disk data store. R owns queries and
computes; Python and bash are atomic.

Two pre-existing bootstrappers live alongside, not inside, the registry
library: `utils/load_bridge.R` (R) and `utils/pipeline_bridge.sh` (shell).
These are universal script-top sourcers that export paths and source
Engine B. Eventual retirement is chat-17+ scope.
