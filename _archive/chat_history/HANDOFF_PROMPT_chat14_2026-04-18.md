# Handoff prompt for chat 14

**Chat 14 is NOT the first HPC run.** Chat 13's tarball went wider than
its mandate — GHSL v5 → v6 was swapped mid-session, and the downstream
wiring around it was only partly fixed. Chat 14 is the "finish what
chat 13 started" session. HPC first run is chat 15 earliest.

**Read at the start:**

1. This file.
2. `SESSION_SUMMARY.md` — current-state snapshot.
3. `AUDIT_LOG_chat13_2026-04-18.md` — findings BG / BH / BI / BJ / BK.
4. `inversion_modules/phase_2_discovery/2e_ghsl/README.md` — the v6
   install notes. It documents what's wired and what isn't.
5. `AUDIT_PHASE4_COHERENCE_2026-04-18.md` §5 findings AS / AW / AV
   only — reference for the eventual threshold calibration.

**Do NOT read** unless a specific earlier decision needs reconstructing:
- `_archive/chat_history/` (pre-chat-13 handoffs and audit logs)
- `_archive_superseded/` (old versions of files — v5 GHSL, old C01f,
  etc.)

---

## Honest notes from chat 13

Chat 13 was long and drifted. Two things the next chat needs to know:

1. **GHSL v5 → v6 swap went wider than scoped.** Chat 13's mandate was
   phase-4 wiring. Near the end of the session, the user correctly
   raised that v5 had been superseded by v6 (rolling-window smoothed,
   heavy/light split). Chat 13 swapped the files in
   `phase_2_discovery/2e_ghsl/` and archived v5 to
   `_archive_superseded/2e_ghsl_v5/`. BUT chat 13 did not update the
   four downstream consumers that still read v5 paths/columns:
   - `phase_4_postprocessing/4b_group_proposal/lib_ghsl_confirmation.R`
     (lines 60–62, 85–86) — reads `<chr>.ghsl_v5.annot.rds` and
     `<chr>.ghsl_v5.karyotypes.rds`. Feeds Tier-3 confirmation in C01i_b.
     **Without this fix, phase-4b Tier-3 returns UNINFORMATIVE for
     every sample — silent.**
   - `phase_2_discovery/2d_candidate_detection/run_all.R` (lines 453–525)
     — reads `<chr>.ghsl_v5.annot.rds` with v5 column names and
     aggregates per-block into the scoring table as `ghsl_v5_score_max`
     etc.
   - `phase_4_postprocessing/4a_existence_layers/STEP_C01d_candidate_scoring_wired_25_v934_registry.R`
     (line 319) — reads `iv$ghsl_v5_score_max` from the scoring table
     (depends on run_all.R).
   - `phase_2_discovery/2d_candidate_detection/STEP_D05_ghsl_stability.R`
     — header says "should match what STEP_C04_snake3_ghsl_v5.R writes."
     Likely reads v5 column names directly.

   The v6 install README flagged this as "wiring gap for chat 14/15"
   but the consumers weren't patched. **This is chat 14's first job.**

2. **"LG12 as known strong inversion" framing in the previous handoff
   was wrong.** LG12 has never been independently confirmed — the user
   described it as "a playtoy mockup." Chat 13 carried this stale
   framing forward from earlier handoff docs. The correct framing for
   an eventual HPC run is:
   - LG12 is the first test target, not a validated positive control.
   - LG25 is paralog-heavy and used as a negative control.
   - If LG12 shows weak signal, that's a valid finding (not evidence
     the pipeline is miscalibrated). Threshold calibration depends on
     genuinely known positives; we don't have one.

3. **Cohort terminology.** 226-sample North African catfish (*Clarias
   gariepinus*), not F1 hybrid. Fixed in SESSION_SUMMARY; check other
   docs still at root and update if you see "F1" or
   "*C. gariepinus × C. macrocephalus*".

---

## Posture for chat 14

Small, careful, biology-grounded. The user explicitly said earlier
Claude versions were better at biology framing than code framing; this
chat should lean into that. If a question is "what does this metric
mean for an INV/INV fish at 9×", answer biology first, then code. If a
question is "is this wired correctly", answer code first.

The primary deliverable is: **GHSL v6 wiring cleanup, properly.**
No mock-HPC launch attempts until that is done.

---

## Primary work: proper GHSL v6 wiring cleanup

User's direction (chat 13 end): **"Proper and put the old in _archive"**
— i.e. do Option B (update all consumers to read v6), not Option A
(v5-compat emit in classifier).

### Plan

**Step 1 — Decide the file interface that v6 offers per-chromosome.**

v5 gave consumers `<chr>.ghsl_v5.annot.rds` (per-window metrics, keyed
by `global_window_id`) and `<chr>.ghsl_v5.karyotypes.rds` (per-sample
stable runs, one row per run).

v6's current outputs (from `STEP_C04b_snake3_ghsl_classify.R`) are
**genome-wide TSVs** — `snake3v6_window_track.tsv.gz` and
`snake3v6_karyotype_calls.tsv.gz`. Consumers that read per-chromosome
would need to subset these.

Cleaner: have the v6 classifier also emit **per-chromosome annot RDS
files** with v6-native column names:
`<outdir>/annot/<chr>.ghsl_v6.annot.rds` and
`<outdir>/annot/<chr>.ghsl_v6.karyotypes.rds`. Same shape as v5 annot
RDS but columns prefixed `ghsl_v6_*`. This is a ~30-line addition at
the end of the classifier's main loop.

**Step 2 — Patch classifier.** Add the per-chrom RDS emit block.

**Step 3 — Archive v5-dependent versions of the four consumers.**
Move pre-patch copies to
`_archive_superseded/ghsl_v5_consumers_pre_v6_rewire/` with a short
README explaining what they were and when they got patched.

**Step 4 — Update each consumer.** Column renames throughout:
`ghsl_v5_score` → `ghsl_v6_score`, `ghsl_v5_score_max` →
`ghsl_v6_score_max`, etc. File-path renames:
`<chr>.ghsl_v5.annot.rds` → `<chr>.ghsl_v6.annot.rds`, same for
karyotypes.

Files to touch:
- `lib_ghsl_confirmation.R` — update `load_ghsl_annot` (line ~83) and
  `load_ghsl_karyo` (line ~58) path lists.
- `run_all.R` phase-2d (lines 453–525) — update annot path and
  scoring-table output column names.
- `STEP_C01d_candidate_scoring_wired_25_v934_registry.R` line 319 —
  rename `iv$ghsl_v5_score_max` → `iv$ghsl_v6_score_max`.
- `STEP_D05_ghsl_stability.R` — update column refs per v6 schema.

**Step 5 — Secondary non-blocking cleanups.**
- `inversion_modules/phase_2_discovery/2d_candidate_detection/00_config.R`
  line 20: `GHSL_DIR <- file.path(CFG$SNAKES, "ghsl_v5")` →
  `"ghsl_v6"` (only if the directory name actually changes; else
  leave alone with a comment).
- `inversion_modules/phase_3_refine/README.md` line 22 — replace
  `STEP_C04_snake3_ghsl_v5.R` with v6 reference.
- `phase_4_postprocessing/4e_final_classification/SPEC_VS_REALITY.md`
  line 52 — update GHSL v5 section to v6.
- `inversion_modules/phase_2_discovery/2c_precomp/patches/README.md` —
  update `patch_C04_ghsl_flashlight.R` target from v5 to v6.

**Step 6 — Test under HPC R.**
- Parse-check every patched file.
- Re-run chat-12 test suites (50/50 + 33/33). Neither should regress;
  GHSL v6 is upstream of them.
- If possible, dry-run the v6 classifier on one chromosome with a
  stub `snake3v6_interval_genotypes.tsv.gz`, then verify
  `load_ghsl_annot` in `lib_ghsl_confirmation.R` finds and reads the
  new `<chr>.ghsl_v6.annot.rds`.

### One design decision chat 14 still needs to make

Column names in the scoring table (`run_all.R` output) that feed C01d:
**rename everywhere to `ghsl_v6_*`** (clean, intellectually honest — no
v5 names masquerading as truth), **or keep `ghsl_v5_*` as a stable
contract** (smaller diff, but confusingly names v6 values with v5
labels). Chat 13 recommended rename-everywhere. Chat 14 can confirm or
flip if there's a good reason.

---

## Additional chat-13 hardening that was partial

### Chunking in v6 heavy engine

`STEP_C04_snake3_ghsl_v6.R` lines 109–110 load every chromosome's
precomp RDS eagerly before applying the `--chrom` filter (line 112).
Wastes ~1–3 GB RAM for single-chromosome runs. Fix: move `readRDS`
into the main loop so only the active chromosome is in memory. Apply
the filter before load. Trivial rewrite, ~5 lines. Not blocking for
single-chrom LG12 dry run; is blocking for the full 28-chrom sweep.

### k=4 and k=5 cluster labels in v6 classifier Part C

`STEP_C04b_snake3_ghsl_classify.R` lines 383–397. At k=2 labels are
`LOW_DIV / HIGH_DIV` (honest). At k=3 labels are `INV_INV / HET /
INV_nonINV` (good). At k=4+ labels fall back to generic `CLASS_1,
CLASS_2, ...` which loses biological meaning.

User's point: "Maybe have many subgroups too." The right fix at k=4
is `INV_INV / INTER_LOW / INTER_HIGH / INV_nonINV` (driven by cluster
center rank). At k=5 add HET in the middle: `INV_INV / INTER_LOW / HET
/ INTER_HIGH / INV_nonINV`.

Open question for chat 14 to think through BEFORE coding: is "many
subgroups" at k=4/5 a divergence-level axis (Part C) or a sub-inversion
boundary axis (Part D — CUSUM changepoint clustering)? If it's both,
the labeling scheme needs two dimensions, not one. Biology question.
Ask before patching.

### "triangle_intervals" naming drift

`STEP_C04b_snake3_ghsl_classify.R` CLI says `--intervals
<triangle_intervals.tsv.gz>`. The file still exists as a compatibility
alias produced by `phase_2_discovery/2d_candidate_detection/STEP_D12_bridge_to_C01d.R`,
which converts `inv_detect_v9.3` output (`blocks_<chr>.tsv`) into the
old triangle schema for downstream back-compat.

There is **no "triangle" script anymore**. The name is historical.
Either rename the CLI arg to `--candidate_intervals` (honest but
breaks launcher scripts), or keep the arg name and fix the docs to
say "intervals from `STEP_D12_bridge_to_C01d.R` — historical name
from `STEP_C01c_triangle_regimes.R` which no longer exists."

Low priority; the intervals themselves are correct.

---

## Chat-13 checklist status (for reference)

All items from the chat-13 handoff closed. New findings opened:

- **BG** (fixed) — `registry_loader.R::write_block` now resolves
  dotted-path keys. Unlocked 21 previously-dropped keys.
- **BH** (fixed) — C01j writes `regime_memberships.tsv.gz` (it was
  advertised but missing). Unblocks C01i_b.
- **BI** (fixed) — seal accepts both canonical
  (`gene_conversion_embedded`) and legacy (`gene_conversion`)
  event_class names; prefers `rec_info$recomb_subgroup` when present.
- **BJ** (fixed) — `build_key_spec()::q2` now lists the 18 new
  chat-13 keys.
- **BK** (open, chat 14/15) — four phase-4b schemas
  (`internal_dynamics`, `recombinant_map`, nested composition, and
  the upstream Q2 producers) have no `keys_extracted` directives.
  51 Q2 spec keys have no schema writer. Adding these is the path to
  the handoff's 85% Q2 coverage aspiration.

**BC metric (as of end of chat 13):** Registry-wide coverage
73/319 = 22.9%; Q2 coverage 22/73 = 30.1%. The 30.1% figure will move
toward 50-60% once chat 14's BK work lands.

---

## Eventually (not in chat 14): first HPC run

After GHSL v6 wiring + BK schema extensions are done, **chat 15 or
later** does the first LANTA run. By then:

- v6 pipeline is wire-complete end to end.
- Tier-3 karyotype confirmation actually fires (currently silent).
- Scoring table has v6 columns from phase-2d feeding D10 in C01d.
- Coverage metric moved to ~50-60%.

The run order and calibration plan from the previous handoff (plot
`structure_score`, `silhouette`, `longest_deviation_bp` distributions
LG12 vs LG25, tune AS / AW / AV) still applies — but **drop the
"known strong inversion" framing**. Run them, plot the distributions,
and interpret what you actually see. If LG12 looks weak, that's
information, not a bug report.

---

## Housekeeping at end of chat 14

```bash
mv HANDOFF_PROMPT_chat14_2026-04-18.md AUDIT_LOG_chat13_2026-04-18.md \
   _archive/chat_history/
```

Create `HANDOFF_PROMPT_chat15_<date>.md`. Update `SESSION_SUMMARY.md`
+ `FIXES_APPLIED.md` in place (roll chat-14 section to top). Tarball
as `inversion-popgen-toolkit_chat14_ghsl_v6_wiring_<date>.tar`.
