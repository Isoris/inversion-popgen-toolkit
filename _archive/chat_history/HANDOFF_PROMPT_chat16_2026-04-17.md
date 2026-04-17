# Handoff prompt for chat 16

**Chat 16 is the first HPC run on LANTA** (same priority as the
chat-14→15 handoff asked for, because chat 15 couldn't execute — no R
interpreter in the chat-15 environment, and chat 15 did the BK
schema-canonicalisation pass instead). Two test batches are now
waiting on real-data validation: the chat-14 GHSL v6 panel path and
the chat-15 three-block BK key extractions.

**Read at the start:**

1. This file.
2. `SESSION_SUMMARY.md` — current-state snapshot (post-chat-15).
3. `FIXES_APPLIED.md` — chat-15 section (BK, BK-rename, BK-docs) and
   chat-14 section (BL / BM / BN / BO / BP / BQ) — both need to execute
   against HPC R in chat 16.
4. `registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md` —
   per-key biological/computational reference for the 41 chat-15
   wired keys. Read the top-of-file intro + the sections for whichever
   block you're debugging.
5. `inversion_modules/phase_2_discovery/2e_ghsl/README.md` — GHSL v6
   wiring status (chat-14).
6. `utils/lib_ghsl_panel_README.md` — panel query API (chat-14).

**Do NOT read** unless a specific earlier decision needs reconstructing:
- `_archive/chat_history/` (pre-chat-15 handoffs and audit logs)
- `_archive_superseded/` (old versions of files — v5 GHSL consumers,
  pre-BK canonical internal_dynamics, the three phase_4b schema
  drafts, etc.)

---

## What chat 15 delivered (non-HPC only)

Chat 15 had no R interpreter and no LANTA access, so ran no code. It
did two things: the BK schema-canonicalisation pass the chat-14
handoff deferred, then — based on user direction mid-chat — a full
ancestry wiring pass that surfaced and fixed a silent no-op that had
been in place since chat 12.

### Phase A — BK schema canonicalisation

- **Three canonical phase-4b schemas** at
  `registries/schemas/structured_block_schemas/`:
  - `internal_dynamics.schema.json` — rewritten from pre-chat-12
    recombinant-flavor to chat-12 decompose-flavor (silhouette / BIC /
    flashlight / seed counts / per_sample PCA class).
  - `recombinant_map.schema.json` — new; chat-12 R + G gate block from
    `STEP_C01i_b_multi_recomb.R`.
  - `internal_ancestry_composition.schema.json` — new; Engine B
    nested-composition block from `STEP_C01i_c_nested_composition.py`.
- **41 `keys_extracted` directives** total across the three schemas.
- **31 scientific renames** (e.g. `q2_decomp_status` →
  `q2_pca_decomp_status`; `q2_silhouette_score` →
  `q2_pca_cluster_silhouette`; `q2_n_cheat2_reclass` →
  `q2_n_hetDEL_breakpoint_reclass`). Full rename table in
  `BK_KEYS_EXPLAINED.md`.
- **`BK_KEYS_EXPLAINED.md`** — ~500-line per-key reference: provenance
  / biological meaning / downstream consumer for every one of the 41
  keys.
- **`compute_candidate_status.R::build_key_spec()`** updated: q2 gains
  36 BK keys + the seal-written `q2_decomp_quality_flags` (distinct
  from the schema's silhouette-derived flag); q6 gains 3 `_prelim`
  counts; q1 entry renamed; 5 stale pre-chat-12 entries removed.
- **Three validators** in work root: `_schema_check.py`,
  `_code_field_check.py`, `_bk_rename.py`.
- **Archive** — pre-BK `internal_dynamics.schema.json` and the three
  phase_4b drafts under `_archive_superseded/bk_schemas_pre_canonical/`.

Coverage: Q2 30.1% → 64.8% (22/73 → 59/91).

### Phase B — Ancestry full-wiring + stats_cache + config rewrite

User flagged mid-chat that "the ancestry bridge is maybe kind of old —
wire fully to registry and to precomp." Investigation found that
`STEP_C01a_precompute.R`'s local-Q merge block had been silently
no-op'ing since chat 12 (broken fallback path lookup for the old
`ancestry_bridge.R`). User also asked: full K=2..20 sweep, persist
Q+F per candidate, update the inversion config "it's been a while,"
and don't mix data with code.

Deliverables:

1. **Fixed the silent no-op** — C01a now calls the already-in-scope
   `merge_local_Q_into_invlikeness()` from
   `unified_ancestry/wrappers/instant_q.R` (sourced via `load_bridge.R`
   since chat 11). Stamps `localQ_delta12_K08`, `localQ_entropy_K08`,
   `localQ_ena_K08` onto every window. Non-fatal if cache not
   populated — prints a clear "run LAUNCH_instant_q_precompute.slurm
   first" message.

2. **K=2..20 sweep in a K-sharded cache**. The instant_q wrapper is
   now fully K-parameterised. Launcher array is (28 × 19) = 532 tasks.
   Per-K qinit/fopt resolved by filename swap (`thin500_K08_seedN` →
   `thin500_K12_seedM` via `BEST_SEED_BY_K`). Only canonical K (8)
   gets flattened into the precomp RDS; other K levels live only in
   the cache.

3. **Sample-set fingerprinting** — every cache filename now includes
   a tag `N<count>_<6char-sha1>` computed from the sorted SAMPLE_LIST.
   Identical implementation in R (wrapper) and bash (launcher), so a
   subset run cannot silently overwrite a 226-cohort run. Legacy
   untagged files still readable via fallback chain.

4. **Data dirs moved out of the code tree**. `LOCAL_Q_DIR` default
   changed from `$BASE/unified_ancestry/local_Q` → `$BASE/ancestry_cache`.
   Three-plane separation: code under `unified_ancestry/` (git),
   pipeline outputs under `$INVDIR/`, mutable state under
   `$BASE/ancestry_cache/` + `$BASE/sample_registry/` + (new)
   `$BASE/registries/data/stats_cache/`.

5. **Ancestry groups auto-registered** — `load_bridge.R` STEP 6.5
   now registers `all_226` (from sample_master), `unrelated_81` (from
   PRUNED_LIST), and `ancestry_K<K>_Q<k>` for k=1..K (majority-vote
   `assigned_pop` across the first 3 cached chroms). Idempotent —
   safe to source repeatedly.

6. **Five new `reg$compute` ancestry methods** (same chat-14 GHSL
   pattern — thin wrappers, NULL+warning if deps missing):
   - `ancestry_at_interval(chrom, s, e, K=NULL)`
   - `ancestry_at_candidate(cid, K=NULL)`
   - `ancestry_q_vector(chrom, s, e, K=NULL, sample_ids=NULL)`
   - `ancestry_q_summary(chrom, K=NULL)`
   - `ancestry_q_and_f_for_candidate(cid, K=NULL, persist=TRUE)`

7. **NEW `reg$stats` subsystem** — fourth registry for persisting
   numerical results that compute methods used to return then
   discard. Layout:
   ```
   $STATS_CACHE_DIR/
   ├── pairwise/<stat>/<chrom>/<g1>__vs__<g2>.<sample_set>.tsv.gz
   ├── candidate/<cid>/Q_K<NN>.<sample_set>.tsv.gz
   ├── candidate/<cid>/F_K<NN>.tsv.gz
   ├── candidate/<cid>/meta.tsv
   └── manifest.tsv
   ```
   Methods: `put_pairwise`, `get_pairwise`, `put_candidate_q_f`,
   `get_candidate_q`, `get_candidate_f`, `list_cached`,
   `clear_candidate`. `pairwise_stat()` gains `cache=TRUE` flag that
   auto-persists. `ancestry_q_and_f_for_candidate(persist=TRUE)` hits
   `reg$stats$put_candidate_q_f()`.

8. **Inversion config rewrite** (`inversion_modules/00_inversion_config.sh`).
   Auto-discovery for CODEBASE / UTILS_DIR / LOAD_BRIDGE so it works
   under both flattened and legacy v8.5 layouts. Added LOCAL_Q_DIR,
   STATS_CACHE_DIR, CANONICAL_K, K_SWEEP, SV_PRIOR_DIR, GHSL_DIR,
   GHSL_PANEL_LIB, REGISTRIES_DATA_DIR, PHASE2/3/4_DIR. NGSADMIX_DIR
   synced to `05_ngsadmix_global` to match the ancestry config (was
   divergent `05_NGSadmix`). `inv_init_dirs()` mkdirs the new data
   dirs.

9. **Archive + docs** — old `inversion_modules/utils/ancestry_bridge.R`
   archived to `_archive_superseded/ancestry_bridge_pre_unified/` with
   README. Stale `unified_ancestry/README_REWIRING_v9_need_rewrite.md`
   archived (the rewiring it planned is now done). `unified_ancestry/README.md`
   fully rewritten with code-vs-data layout and K-sweep commands.
   `2c_precomp/README.md` gains ancestry precompute section.
   `API_CHEATSHEET.md` extended with chat-14 GHSL, chat-15 ancestry,
   and new reg$stats tables.

What chat 15 did NOT do:

- Any HPC run. The three BK schemas have not executed against real
  candidates. The new ancestry wiring has not been run end-to-end.
  The K-sweep launcher has not submitted 532 tasks.
- Re-run chat-12 unit tests (50/50 DAG + 33/33 GC detector). No R.
- Validate that the GHSL v6 panel path (from chat 14) works on real
  data.
- Tree persistence per candidate — user said "Q and F and the tree or
  if no tree at least we need them for that interval." Q + F
  delivered. Tree deferred (requires picking a tree library and
  format — can reconstruct offline from Q + F as needed).


---

## Posture for chat 16

**Priority order.** Largely identical to what chat-14's handoff asked
of chat 15, plus the chat-15 BK deliverables now need their own
live-run validation.

1. **Validate the BK schemas extract keys on real data.** Run one of:
   - `STEP_C01i_decompose.R` on a chrom + candidate set → inspect
     `<cand_outdir>/structured/internal_dynamics.json` exists AND
     `<cand_outdir>/keys.tsv` contains the 19 expected keys (16 q2_* +
     3 q6_*_prelim).
   - `STEP_C01i_b_multi_recomb.R` → same for `recombinant_map.json`
     and 10 expected keys.
   - `STEP_C01i_c_nested_composition.py` → same for
     `internal_ancestry_composition.json` and 12 expected keys.
   If any key is missing, the most likely cause is that my
   `resolve_dotted()` call (chat-13 BG) didn't find the field at
   runtime — check the block JSON first for the field's presence, then
   trace through `write_block` in `registries/api/R/registry_loader.R`.
2. **Validate the v6 end-to-end GHSL path works** on LG12. Run
   STEP_C04 heavy engine + STEP_C04b classifier. Confirm all three
   per-chrom RDS land
   (`annot/*.ghsl_v6.annot.rds`,
    `annot/*.ghsl_v6.karyotypes.rds`,
    `per_sample/*.ghsl_v6.per_sample.rds`) and that
   `reg$compute$ghsl_at_candidate("LG12_17")` (or any real LG12 cid)
   returns sane per-sample rows. See `BK_KEYS_EXPLAINED.md` for what
   `sane` looks like per key.
3. **Unit tests.** Re-run 50/50 DAG + 33/33 GC detector under HPC R.
   No chat-12 code touched in 13/14/15; no regression expected;
   verification needed.
4. **LG25 paralog control.** Run the same two paths. Expect LG25 panel
   data to look messy — that's the null.
5. **Sweep remaining 26 chroms** once LG12 + LG25 look right.
6. **Thresholds.** v5-inherited defaults stand. If LG12 shows weak
   signal that's a finding, not a miscalibration (user direction from
   chat 15).

---

## Smoke-test runbook for ancestry wiring (chat 15 Phase B)

Before the BK validation block runs, the ancestry cache must be
populated — C01a's new localQ stamps read from it, and so do the new
`reg$compute$ancestry_*` methods.

**One-time setup on LANTA.**

```bash
source "$BASE/utils/pipeline_bridge.sh"         # pulls in all env vars
# Or, if the flattened layout isn't present yet on the HPC box:
#   source "$BASE/inversion_modules/00_inversion_config.sh"
#   source "$BASE/unified_ancestry/00_ancestry_config.sh"

# Verify config resolved the new data dirs correctly
echo "LOCAL_Q_DIR=$LOCAL_Q_DIR"
echo "STATS_CACHE_DIR=$STATS_CACHE_DIR"
echo "CANONICAL_K=$CANONICAL_K"
echo "K_SWEEP=$K_SWEEP"
# Expected:
#   LOCAL_Q_DIR=/scratch/.../Quentin_project_KEEP_2026-02-04/ancestry_cache
#   STATS_CACHE_DIR=/scratch/.../Quentin_project_KEEP_2026-02-04/registries/data/stats_cache
#   CANONICAL_K=8
#   K_SWEEP=2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
```

**Step A — canonical-K only first (28 tasks, sanity check).**

```bash
K_SWEEP=8 sbatch "$BASE/unified_ancestry/launchers/LAUNCH_instant_q_precompute.slurm"
# ~10 min array-parallel. Check:
ls "$LOCAL_Q_DIR/K08/"
# Expect: C_gar_LG01.N226_xxxxxx.local_Q_{samples,summary,meta}.tsv.gz × 28 chroms
cat "$LOCAL_Q_DIR/manifest.tsv" | head
# Expect columns: chrom K sample_set summary_file samples_file timestamp
```

**Step B — verify C01a now stamps localQ columns.**

```bash
Rscript "$BASE/inversion_modules/phase_2_discovery/2c_precomp/STEP_C01a_precompute.R" \
    --step10_out <step10_prefix> \
    --outdir     <workdir>/precomp_LG12 \
    --dosage_dir "$DOSAGE_DIR"

# Inspect the output RDS
Rscript -e '
  x <- readRDS("<workdir>/precomp_LG12/C_gar_LG12.precomp.rds")
  cat("columns with localQ:\n")
  print(grep("^localQ_", names(x), value = TRUE))
  cat("n non-NA delta12_K08:", sum(!is.na(x$localQ_delta12_K08)), "/", nrow(x), "\n")
'
# Expected 4 columns present (delta12, entropy, ena, assigned_pop all with _K08 suffix);
# n non-NA should be the full window count for the chromosome.
```

If localQ columns are absent or all NA, the cache step failed or the
`merge_local_Q_into_invlikeness()` returned empty. Look for
`[instant_q.R] Merging K=8 from ... manifest rows` in the precomp log.

**Step C — full K sweep (if A+B pass).**

```bash
sbatch "$BASE/unified_ancestry/launchers/LAUNCH_instant_q_precompute.slurm"
# 532 tasks, probably 30-60 min elapsed with normal queue concurrency.
# Total disk: ~21 GB under $LOCAL_Q_DIR.
```

**Step D — sanity the registry ancestry query path.**

```r
source(Sys.getenv("LOAD_BRIDGE"))
reg <- load_registry()

# Auto-registered groups
reg$has_group("all_226")       # expect TRUE
reg$has_group("unrelated_81")  # expect TRUE
reg$has_group("ancestry_K8_Q1") # expect TRUE after cache populated

# Query a candidate
reg$compute$ancestry_at_candidate("LG12_17")
# Expect data.table with window_id / chrom / start_bp / end_bp / mean_delta12 / mean_entropy / mean_ena

# Persist Q + F
qf <- reg$compute$ancestry_q_and_f_for_candidate("LG12_17", K = 8, persist = TRUE)
str(qf)
# Expect list with Q (data.table, ~226 × n_windows rows), F (matrix n_sites × 8)
list.files(file.path(Sys.getenv("STATS_CACHE_DIR"), "candidate", "LG12_17"))
# Expect: Q_K08.N226_xxxxxx.tsv.gz  F_K08.tsv.gz  meta.tsv
```

**Step E — sanity a pairwise FST with caching.**

```r
g_ref <- reg$get_group("ancestry_K8_Q3")
g_alt <- reg$get_group("ancestry_K8_Q5")
r <- reg$compute$pairwise_stat(
       position = 20e6, flank_size = 500e3,
       group1 = g_ref, group2 = g_alt,
       stat = "fst", chrom = "C_gar_LG12", cache = TRUE)
# Check cached file lands:
list.files(file.path(Sys.getenv("STATS_CACHE_DIR"), "pairwise", "fst", "C_gar_LG12"))
# Expect a file: ancestry_K8_Q3__vs__ancestry_K8_Q5.N226_xxxxxx.tsv.gz (or a similar canonical-sort pair label)
```

If all five steps pass, the chat-15 Phase B wiring is validated.
Then proceed with the BK schema validation below.

---

## Smoke-test runbook for BK validation

For one candidate on one chrom, the whole BK validation is four
commands + an eyeball. Example for LG12, candidate LG12_17:

```bash
# Prereqs — phase-2 outputs already exist; regime_memberships.tsv.gz
# from C01j; per_window_class.rds from decompose first.

# 1. decompose writes internal_dynamics
Rscript inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_decompose.R \
  --candidates <c01d_out>/candidates.tsv \
  --outdir     <workdir>/decompose_out \
  --seed_tsv   <step03>/04_deconvolution_seeds/LG12_17_seeds.tsv

# 2. multi_recomb writes recombinant_map
Rscript inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_b_multi_recomb.R \
  --candidates <c01d_out>/candidates.tsv \
  --regime_memb <c01j_out>/regime_memberships.tsv.gz \
  --decomp_dir  <workdir>/decompose_out \
  --outdir      <workdir>/multi_recomb_out

# 3. nested_composition writes internal_ancestry_composition
python3 inversion_modules/phase_4_postprocessing/4b_group_proposal/STEP_C01i_c_nested_composition.py \
  --candidates  <c01d_out>/candidates.tsv \
  --q_cache_dir <engineb>/local_Q_cache \
  --outdir      <workdir>/nested_comp_out

# 4. Inspect the per-candidate keys.tsv
cid=LG12_17
cat <registry>/per_candidate/$cid/keys.tsv | cut -f1 | sort -u
```

Expected keys in `keys.tsv` for LG12_17 (41 BK keys + any pre-existing
keys from earlier phase-4a blocks):

From `internal_dynamics`:
```
q2_pca_decomp_status
q2_pca_decomp_skip_reason       (only if skipped; typically absent)
q2_pca_seed_source
q2_n_seed_class_conflicts
q2_pca_cluster_silhouette
q2_pca_cluster_separation_flag
q2_bic_gap_k3_vs_k2
q2_het_phase_support_fraction
q2_n_het_phase_supported
q2_pca_seed_mode
q2_n_seed_HOM_REF
q2_n_seed_HET
q2_n_seed_HOM_INV
q2_n_samples_sv_pca_discordant
q2_n_hetDEL_breakpoint_reclass
q2_pca_k_used
q6_n_HOM_REF_prelim
q6_n_HET_prelim
q6_n_HOM_INV_prelim
```

From `recombinant_map`:
```
q2_n_recombinant_samples
q2_n_recombinant_gc_samples
q2_n_recombinant_dco_samples
q2_n_regime_ghsl_disputed_samples
q2_n_ghsl_split_only_samples
q2_recomb_posterior_source
q2_recomb_gate_rule_version
q2_recomb_min_regime_dev_fraction
q2_recomb_min_regime_dev_bp
q2_recomb_dco_threshold_bp
```

From `internal_ancestry_composition`:
```
q1_composite_flag
q1_ancestry_dominant_pattern
q2_pct_samples_two_block_ancestry
q2_pct_samples_multi_block_ancestry
q2_pct_samples_homogeneous_ancestry
q2_pct_samples_gradient_ancestry
q2_pct_samples_diffuse_ancestry
q2_mean_ancestry_fragmentation
q2_mean_ancestry_entropy_within_sample
q2_mean_ancestry_switches_per_sample
q2_ancestry_K_used
q2_ancestry_n_samples_analyzed
```

**If any key is missing from `keys.tsv`.** The registry loader's
`write_block` skips keys whose values resolve to NULL or NA. This is
intentional — but for the smoke test you want all keys present, so
feed a candidate that exercises every code path. For example,
`q2_pca_decomp_skip_reason` is only populated on skipped candidates;
that's fine to see missing on a successful run. `q2_recomb_posterior_source`
should always be present because it's a fixed string.

---

## Threshold plan (unchanged from chat-14 handoff except posture)

User direction chat 15: defaults stand unless observed distributions
break something. The defaults to NOT tune unless the data argues:

| Threshold | Default | Where | Handoff name |
|---|---|---|---|
| `structure_score` cutoffs in C01j | 0.35 / 0.6 | C01j | AS |
| `min_dco_bp` in multi_recomb | 200 kb | recombinant_map gate_params | AW |
| Silhouette cutoff for `decomp_quality` flag | 0.40 | decompose | AV |
| `KARYO_LO` / `KARYO_HI` quantiles | 0.15 / 0.70 | GHSL classifier | (v5-inherited) |

If LG12 signal looks weak, report as finding; don't retune to force a
positive.

---

## BK-docs maintenance protocol

When a new phase-4 block is added or an existing block writer changes
its output fields:
1. Update the block's JSON schema in
   `registries/schemas/structured_block_schemas/<block>.schema.json`
   (the only path `registry_loader.R::load_schema` reads).
2. Add / update `keys_extracted` with scientific names (see the
   rename rationale in `BK_KEYS_EXPLAINED.md` §Pre-rename / post-rename
   key map).
3. Add a section to `BK_KEYS_EXPLAINED.md` with provenance /
   biological meaning / downstream consumer for each new key.
4. Register new keys in
   `compute_candidate_status.R::build_key_spec()::q1/q2/q3/q4/q5/q6/q7`.
5. Run `_schema_check.py` + `_code_field_check.py` from work root.
6. If you renamed keys, run `_bk_rename.py`-style migration over the R
   files that consume them (grep `"q2_old_name"` and replace).

---

## Known caveats heading into chat 16

1. **`_code_field_check.py` R parser has a known limitation.** Its
   strict comma-splitter under-reports fields in R list() bodies
   where a value is a multi-line if/else expression; a permissive
   fallback (line-level `^<ws><ident> =<not-=>` regex) catches them.
   Not an issue for the chat-15 schemas (all 41 head segments resolve
   under the combined pass), but if a future block's R writer has
   complex embedded function calls as values, verify by grep before
   ruling a field missing.
2. **`q2_decomp_quality_flags` vs `q2_pca_cluster_separation_flag`**
   are two DIFFERENT keys that look similar. The first is seal-written
   (comma-joined aggregate of validation quality_flags); the second is
   the silhouette-derived clean/noisy flag from
   `internal_dynamics.schema.json`. The chat-15 commit comments make
   this explicit, but if they diverge in meaning later, adjust names
   before it becomes a long-term confusion.
3. **`q1_composite_flag` is v10.1 canonical (NOT renamed).** Chat-9
   Finding X removed the duplicate `q1_ancestry_composite_flag` alias.
   Don't re-rename in a later pass.
4. **Three phase_4b drafts now archived**
   (`_archive_superseded/bk_schemas_pre_canonical/phase4b_drafts/`).
   `frequency.v2.schema.json` is still there (out of BK scope — no
   canonical `frequency.schema.json` exists yet). If BC coverage work
   resumes, that's the next natural target.
5. **v5 back-compat branches in `run_all.R` and
   `lib_ghsl_confirmation.R`** (chat 14) — still in. Once the HPC
   sweep completes all 28 chroms with v6, remove.
6. **Panel memoization** — `load_ghsl_panel` caches per-chrom in
   memory; call `panel_clear_cache(chrom)` between sweep iterations
   (chat-14 caveat).
7. **SPLIT detection still routes through run-overlap Tier-3.** The
   panel-based `ghsl_per_sample_panel_in_interval` is additive — it
   does NOT replace the run-overlap logic. Intentional (chat-14
   caveat).

---

## Housekeeping at end of chat 16

```bash
mv HANDOFF_PROMPT_chat16_2026-04-17.md \
   _archive/chat_history/
# plus chat-16's own audit log if one is produced
```

Create `HANDOFF_PROMPT_chat17_<date>.md`. Update `SESSION_SUMMARY.md`
+ `FIXES_APPLIED.md` in place (roll chat-16 section to top). Tarball
as `inversion-popgen-toolkit_chat16_<topic>_<date>.tar`.
