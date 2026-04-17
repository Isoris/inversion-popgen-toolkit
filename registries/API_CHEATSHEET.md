# Registry API Cheatsheet

One-page reference for the three registries across all three language bindings.

## The three registries (proper names)

| Data dir | Role | What it answers |
|---|---|---|
| `data/sample_registry/` | WHO | Which samples are in which group? Family/pedigree? |
| `data/interval_registry/` | WHERE | What genomic region? What windows? Which cov file? |
| `data/evidence_registry/` | WHAT | What do we know about candidate X? |

## The three language bindings

| Binding | File | Role | Scope |
|---|---|---|---|
| R | `api/R/registry_loader.R` | Full API | atomic + queries + composites + live compute |
| Python | `api/python/registry_loader.py` | Writer + reader | atomic reads + Tier-2 block writes |
| Bash | `api/bash/registry_loader.sh` | SLURM gate-checks | path resolution + validation gates |

**Rule:** Queries and computes live in R. Python and bash are atomic only.

## Entry point

```r
source("registries/api/R/registry_loader.R")
reg <- load_registry()               # auto-detects REGISTRIES env or BASE
reg <- load_registry("/path/to/registries")   # explicit
reg$status()                         # print counts
```

```python
from registry_loader import load_registry
reg = load_registry()                # auto-detects REGISTRIES env or BASE
reg = load_registry("/path/to/registries")
```

```bash
source registries/api/bash/registry_loader.sh
registry_resolve_paths
registry_check_validation LG12_17 SUPPORTED && echo ok || echo not_yet
```

## reg$samples — the sample_registry

| Method | R | Py | Returns | Purpose |
|---|:-:|:-:|---|---|
| `get_master()` | ✓ | ✓ | data.table/list | full sample_master |
| `get_group(gid)` | ✓ | ✓ | chr[] | members of group |
| `has_group(gid)` | ✓ | ✓ | bool | exists? |
| `add_group(gid, sids, ...)` | ✓ | ✓ | bool | register new group |
| `list_groups()` | ✓ | ✓ | data.table/list | all groups |
| `count_groups()` | ✓ | — | int | N groups |
| `get_carriers(cid, status)` | ✓ | ✓ | chr[] | `inv_<cid>_<status>` members |
| **— chat 11 additions —** | | | | |
| `get_sample_metadata(sid)` | ✓ | ✓ | list/dict | sample_master row as list |
| `get_sample_groups(sid)` | ✓ | ✓ | chr[] | reverse lookup: all groups containing sid |
| `get_family(sid)` | ✓ | ✓ | chr/NA | family_id from master (NA if no column) |
| `list_families()` | ✓ | ✓ | chr[] | distinct family_ids |
| `get_family_members(fam_id)` | ✓ | ✓ | chr[] | samples in a family |
| `list_carriers_across_candidates(status)` | ✓ | ✓ | chr[] | union of `inv_*_<status>` members |
| `compute_group_overlap(gid1, gid2)` | ✓ | ✓ | list/dict | inter/union/jaccard/diffs |
| `list_recombinant_candidates(sid)` | ✓ | ✓ | chr[] | CIDs where sid is RECOMBINANT |
| `get_groups_for_candidate(cid)` | ✓ | ✓ | named list | all 7 status groups in one call |
| `find_co_segregating_groups(min_jaccard)` | ✓ | ✓ | data.table/list | pairs with high membership overlap |

**Group naming convention:** `inv_<candidate_id>_<STATUS>` where
STATUS ∈ {`HOM_REF`, `HET`, `HOM_INV`, `RECOMBINANT`, `RECOMBINANT_GC`,
`RECOMBINANT_DCO`, `HOM_STD`}. `HOM_STD` is v9 alias for `HOM_REF`.

## reg$intervals — the interval_registry

| Method | R | Py | Returns | Purpose |
|---|:-:|:-:|---|---|
| `get_windows(chrom?, start?, end?)` | ✓ | — | data.table | window rows |
| `get_window_at(chrom, pos)` | ✓ | — | list/NULL | window containing pos |
| `get_candidate(cid)` | ✓ | ✓ | list/dict | candidate row |
| `get_candidate_boundaries(cid)` | ✓ | ✓ | list | `{left_bp, right_bp}` |
| `add_candidate(cid, chrom, s, e, scale?, parent_id?)` | ✓ | ✓ | bool | register; skips if exists |
| `count_candidates()` | ✓ | — | int | N candidates |
| `get_cov(scope, chrom?)` | ✓ | — | data.table | cov paths |
| **— chat 11 additions —** | | | | |
| `get_children(cid)` | ✓ | ✓ | chr[] | direct children (parent_id-based) |
| `get_parent(cid)` | ✓ | ✓ | chr/NA | direct parent |
| `get_ancestors(cid)` | ✓ | ✓ | chr[] | walk up chain, closest first |
| `get_descendants(cid, depth?)` | ✓ | ✓ | chr[] | BFS down; depth caps levels |
| `get_overlapping(chrom, s, e, exclude?)` | ✓ | ✓ | chr[] | CIDs whose interval overlaps |
| `get_nested_within(cid)` | ✓ | ✓ | chr[] | coord-based, ignores parent_id |
| `classify_relationship(cid1, cid2)` | ✓ | ✓ | chr | equal/nested_N_in_M/partial_overlap/disjoint |
| `update_candidate(cid, ...)` | ✓ | ✓ | bool | partial column update (e.g. set parent_id post-hoc) |
| `bulk_add_candidates(dt)` | ✓ | ✓ | int | batch insert; returns N inserted |

## reg$evidence — the evidence_registry

| Method | R | Py | Returns | Purpose |
|---|:-:|:-:|---|---|
| `write_block(cid, block_type, data)` | ✓ | ✓ | list/dict | write Tier-2 JSON + extract keys per schema |
| `read_block(cid, block_type)` | ✓ | ✓ | list/dict/NULL | read Tier-2 JSON |
| `get_keys(cid, key?)` | ✓ | ✓ | data.table/list | flat keys.tsv rows |
| `add_evidence(cid, key, value, ...)` | ✓ | ✓ | bool | write flat key (v9 compat) |
| `has_evidence(cid, key)` | ✓ | ✓ | bool | flat key present? |
| `get_evidence(cid, key?)` | ✓ | — | data.table | alias of get_keys |
| `register_candidate(cid, chrom, s, e, tier, score)` | ✓ | — | bool | ensure dir + add `scored` key |
| `list_candidates()` | ✓ | ✓ | chr[] | all registered CIDs |
| `count_candidates()` | ✓ | — | int | N candidates |

**Untouched in chat 11 (mature; Phase 3 wiring is the reference writer).**

## reg$query — composite reads (R only)

| Method | Returns | Purpose |
|---|---|---|
| `boundary_flanking_groups(cid, side, offset_bp, flank?)` | list | window at breakpoint + 4 carrier groups |
| `inversion_carriers(cid, validation?)` | list | all carrier groups, gated by validation level |
| `validation_status(cid)` | chr | current q6_group_validation level |
| `completion_per_question(cid)` | named int | N keys per Q1–Q7 |
| `candidate_summary(cid)` | list | candidate row + validation + tier + verdict |
| **— chat 11 additions —** | | |
| `nested_family_tree(root_cid)` | nested list | recursive tree for characterization |
| `overlap_clusters(chrom)` | data.table | connected-component grouping by overlap |
| `sample_inversion_load(sid)` | list | per-status counts for one sample |

## reg$compute — live Engine B stats (R only)

| Method | Returns | Purpose |
|---|---|---|
| `pairwise_stat(pos, flank, g1, g2, stat, chrom, cache=FALSE)` | list/NULL | FST/dxy/etc. via get_region_stats(). `cache=TRUE` auto-persists to reg$stats |
| `boundary_fst_profile(cid, side, distances_kb)` | data.table | FST at multiple distances from breakpoint |
| **— chat 14 additions (GHSL v6) —** | | |
| `ghsl_at_interval(chrom, s, e, ...)` | list | GHSL panel summaries for an interval |
| `ghsl_at_candidate(cid, ...)` | list | same, looked up by cid |
| `ghsl_at_subblocks(cid, ...)` | list | sub-block GHSL |
| `ghsl_at_block_subregions(cid, ...)` | list | sub-region GHSL inside a block |
| `ghsl_wide_matrix(cid_or_chrom, s, e, metric, ...)` | matrix | sample × window GHSL matrix |
| **— chat 15 additions (Ancestry) —** | | |
| `ancestry_at_interval(chrom, s, e, K=NULL)` | data.table | per-window summary (mean_delta12, mean_entropy, mean_ena) |
| `ancestry_at_candidate(cid, K=NULL)` | data.table | same, looked up by cid |
| `ancestry_q_vector(chrom, s, e, K=NULL, sample_ids=NULL)` | data.table | per-sample × per-window Q matrix for a region |
| `ancestry_q_summary(chrom, K=NULL)` | data.table | whole-chrom window summary |
| `ancestry_q_and_f_for_candidate(cid, K=NULL, persist=TRUE)` | list(Q, F, ...) | candidate Q + F; auto-persists to reg$stats |

Requires `utils/load_bridge.R` to have been sourced (or Engine B / instant_q
wired separately). Falls back to NULL with warning if dependencies missing.

## reg$results — persisted numerical results (chat 16 rewrite, R + Python)

The fourth first-class registry alongside sample/interval/evidence. Replaces
chat-15's `reg$stats` (kept as a deprecated alias for one cycle). See
`DATABASE_DESIGN.md` for the full design rationale.

**What changed from chat 15:**
- Filenames use registered `group_id` strings (e.g. `all_226`, `ancestry_K8_Q3`),
  not opaque content hashes. No more R-vs-bash tag mismatches.
- Every manifest row has a UUID `row_id`, full provenance
  (`source_script`/`engine`/`engine_version`/`run_id`/`config_hash`/`upstream_files`),
  and sha256 of its data file.
- `put_*` methods enforce FK constraints: group_id must be in sample_registry,
  candidate_id must be in interval_registry — or the write refuses.
- Group versions (`created` timestamp) are copied into the manifest on write.
  `integrity_check()` flags stale results when a group's membership changes.
- One new `ask()` method with four filter dimensions + interval overlap semantics.

### Writer methods

| Method | Returns | Purpose |
|---|---|---|
| `put_pairwise(chrom, group_1, group_2, stat, dt, source_script, ...)` | path | FST/dxy/etc. between two registered groups |
| `put_candidate_q_f(cid, q_dt, f_mat, K, sample_group, source_script, ...)` | list(q_file, f_file) | candidate Q (per-group) + F (no group) matrices |
| `put_interval_summary(chrom, start_bp, end_bp, group_1, stat, dt, K, source_script, ...)` | path | per-interval ancestry summary (delta12 / entropy / ena / ancestry_q_mean) |

All writers accept `engine`, `engine_version`, `run_id`, `config_hash`,
`upstream_files` as optional provenance kwargs.

### Reader methods (direct file access)

| Method | Returns |
|---|---|
| `get_pairwise(chrom, group_1, group_2, stat)` | data.table/NULL |
| `get_candidate_q(cid, K, sample_group)` | data.table/NULL |
| `get_candidate_f(cid, K)` | matrix/NULL |

### Query plane — one method, four filters

| Method | Purpose |
|---|---|
| `ask(where, who, what, kind, K, overlap="any")` | data.table of manifest rows matching filters. `where` = list(chrom, start_bp, end_bp, candidate_id); `who` = one group or vector; `what` = stat enum; `overlap` ∈ {"any","contains","contained_by"} |
| `ask_what_at(chrom, start_bp=NULL, end_bp=NULL)` | everything at/overlapping an interval |
| `ask_what_for_candidate(cid)` | everything for a candidate |
| `ask_what_for_group(gid)` | everything involving a group |
| `ask_provenance(row_id)` | full provenance block for one row |

### Integrity & session

| Method | Purpose |
|---|---|
| `integrity_check(check_sha256=FALSE)` | 6 FK/version/file-exists/orphan checks. Returns data.table with `attr(x,"all_pass")` |
| `session_start()` / `session_summary()` | bracket a script to echo "wrote N new rows in Ts" at the end |
| `list_cached(kind=NULL)` | alias of manifest read, filterable by kind |
| `clear_candidate(cid)` | delete all files + manifest rows for one candidate |

### Filename convention (no hashes)

```
pairwise/<stat>/<chrom>/<group_1>__vs__<group_2>.tsv.gz      (pair sorted)
candidate/<cid>/Q_K<NN>.<group>.tsv.gz                        per-sample Q
candidate/<cid>/F_K<NN>.tsv.gz                                global F
candidate/<cid>/meta.tsv                                      audit trail
interval/<chrom>/<start_bp>_<end_bp>.<group>.<stat>.K<NN>.tsv.gz
```

### Common use

```r
# 1. Compute + auto-persist a pairwise FST; cache=TRUE routes through reg$results
r <- reg$compute$pairwise_stat(pos, 10000,
       reg$get_group("ancestry_K8_Q3"),
       reg$get_group("ancestry_K8_Q5"),
       stat = "fst", chrom = "C_gar_LG12", cache = TRUE)

# 2. Archive Q + F for a candidate
qf <- reg$compute$ancestry_q_and_f_for_candidate("LG12_17", persist = TRUE)

# 3. Find everything we've computed for a candidate, in one call
reg$results$ask_what_for_candidate("LG12_17")

# 4. Find everything overlapping a region
reg$results$ask_what_at("C_gar_LG12", 1e6, 1.1e6)

# 5. Compute-if-missing pattern (classification loop)
for (cid in reg$evidence$list_candidates()) {
  have <- reg$results$ask_what_for_candidate(cid)
  if (!nrow(have[kind == "candidate_q" & K == 8])) {
    reg$compute$ancestry_q_and_f_for_candidate(cid, K = 8, persist = TRUE)
  }
}

# 6. Before submitting to journal: prove the database is consistent
ic <- reg$results$integrity_check()
stopifnot(attr(ic, "all_pass"))

# 7. Before/after echo (the script-level check)
reg$results$session_start()
# ... script does its put_* calls ...
reg$results$session_summary()
# prints: [results_registry] session wrote 3 new rows across 1 candidate
#         in 42.1s (candidate_q=1, candidate_f=1). Total manifest: 847 rows.
```

### Python parity (atomic writer + read-only query)

```python
from registry_loader import load_registry
reg = load_registry()
reg.results.put_pairwise(chrom="C_gar_LG12", group_1="all_226",
                         group_2="unrelated_81", stat="fst",
                         rows=[...], fieldnames=[...],
                         source_script="my_script.py")
rows = reg.results.ask(where={"chrom": "C_gar_LG12"}, what="fst")
```

## Bash binding

| Function | Purpose |
|---|---|
| `registry_resolve_paths` | export REGISTRIES, SAMPLE_REGISTRY, RESULTS_REGISTRY, ... |
| `registry_list_candidates_by_tier MAX_TIER` | echo CIDs where q7_tier ≤ N |
| `registry_check_validation CID MIN` | exit 0/1 based on q6_group_validation rank |
| `registry_get_validation CID` | echo current level |
| `registry_has_group GID` | exit 0 if member file exists |
| **— chat 16 (results_registry) —** | |
| `registry_results_path <kind> ...` | echo expected path for a result artifact (pairwise / candidate_q / candidate_f / interval_summary) |
| `registry_results_has <kind> ...` | exit 0 if the file exists (gate for "skip if cached") |
| `registry_results_count_by_kind <kind>` | echo row count for a given kind in manifest.tsv |

## Common patterns

### C01f-style: get all carrier bands for a candidate

```r
# Old: 4 has_group + 4 get_group calls
# New: 1 call
grs <- reg$samples$get_groups_for_candidate(cid)
ref_samps <- grs$HOM_REF
het_samps <- grs$HET
inv_samps <- grs$HOM_INV
```

### Phase 4e nested-system characterization

```r
tree <- reg$query$nested_family_tree(cid)
if (tree$n_children >= 2) {
  # complex_system — this candidate has multiple nested children
} else if (tree$n_children == 1) {
  # nested_parent — one inner inversion
}
```

### Chat 4c competing-hypothesis detection

```r
clusters <- reg$query$overlap_clusters(chrom)
# any cluster with size > 1 needs overlap reconciliation
```

### Detecting repeat-called events

```r
coseg <- reg$samples$find_co_segregating_groups(min_jaccard = 0.9)
# pairs with jaccard ≥ 0.9 may be the same event at slightly different coords
```

### Sample-level QC: flag unusually high inversion load

```r
for (sid in sample_ids) {
  load <- reg$query$sample_inversion_load(sid)
  if (load$HOM_INV > threshold) flag(sid, "high_inversion_load")
}
```

## Findings open after chat 11

- **AL** — `load_registry` shadowing between `utils/sample_registry.R` and
  `registries/api/R/registry_loader.R`. Deferred until rename in chat 12.
- **AM, AN** — orphan duplicate files `inversion_modules/utils/sample_*` —
  resolved (deleted) in chat 11 tarball.

See `AUDIT_LOG_chat11_2026-04-17.md` for full list.
