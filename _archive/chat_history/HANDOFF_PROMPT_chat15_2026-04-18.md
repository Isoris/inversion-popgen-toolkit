# Handoff prompt for chat 15

**Chat 15 is the first HPC run on LANTA.** The GHSL v6 wiring gap that
chat 13 left and chat 14 closed is now resolved end-to-end. The pipeline
can run from phase 2 through phase 4 with v6 panel outputs feeding
Tier-3 karyotype confirmation, run_all.R 2d block scoring, and C01d
Layer C scoring. Chat 14 also added a per-sample panel backend with
on-demand registry query API (`reg$compute$ghsl_*`). None of it has
been executed on real data yet.

**Read at the start:**

1. This file.
2. `SESSION_SUMMARY.md` — current-state snapshot.
3. `FIXES_APPLIED.md` — chat-14 section (BL / BM / BN / BO / BP / BQ +
   doc fixes). Full chat-14 detail on what was wired.
4. `inversion_modules/phase_2_discovery/2e_ghsl/README.md` — current
   output layout + wiring status.
5. `utils/lib_ghsl_panel_README.md` — panel query API.
6. `AUDIT_PHASE4_COHERENCE_2026-04-18.md` §5 findings AS / AW / AV —
   threshold calibration targets. **Note the honest framing:** LG12 is
   the first test target, not a validated positive control. LG25 is
   paralog-heavy and used as a negative control. If LG12 shows weak
   signal, that's a finding, not evidence of miscalibration.

**Do NOT read** unless a specific earlier decision needs reconstructing:
- `_archive/chat_history/` (pre-chat-14 handoffs and audit logs)
- `_archive_superseded/` (old versions of files — v5 GHSL, pre-v6
  consumer files, etc.)

---

## Honest notes from chat 14

Chat 14 went wider than the handoff's stated scope. The original
handoff asked for a narrow wiring cleanup: update four consumers to
read v6 paths and columns. Discussion with the user (Quentin) during
the session led to a broader design: a per-sample panel backend (dense
long-format RDS per chromosome, multi-scale divergence + rank + band,
joined with Part B stable-run membership) with on-demand query functions
and registry integration. The reasoning is in the chat-14 turn-level
design discussion; the biological framing is captured in
`utils/lib_ghsl_panel_README.md`.

Concretely, chat 14 delivered:

- Classifier emits three per-chrom RDS files (annot, karyotypes,
  per-sample panel) instead of only genome-wide TSVs.
- `utils/lib_ghsl_panel.R` — ~500-line backend library with range /
  aggregate / sub-block / wide-matrix / plot helpers, memoized per-chrom.
- `registries/api/R/registry_loader.R::load_compute_api()` extended
  with `ghsl_at_candidate`, `ghsl_at_interval`, `ghsl_at_subblocks`,
  `ghsl_at_block_subregions` (results-aware — pulls sub-blocks from
  `boundary`, `regime_segments`, `recombinant_map` blocks), and
  `ghsl_wide_matrix`.
- Heavy-engine chunking fix (single-chrom runs no longer load all 28
  chroms' precomp into memory).
- Heavy-engine default scale ladder widened from `20, 50, 100` to
  `10, 20, 30, 40, 50, 100`.
- Part C k=4 / k=5 labels use the divergence-rank scheme
  (`INV_INV / INTER_LOW / [HET] / INTER_HIGH / INV_nonINV`) instead of
  generic `CLASS_N`.
- All four documented v5 → v6 consumer renames, plus the 2e / 3 /
  2c / 4e doc updates.
- Archive of pre-patch copies at
  `_archive_superseded/ghsl_v5_consumers_pre_v6_rewire/`.

**What chat 14 did NOT do:**
- No HPC run. Everything is parse-checked (`_rcheck.py`) but has not
  executed against real data. Chat 15 is the first live test.
- No re-run of chat-12 unit tests (50/50 DAG + 33/33 GC detector).
  No R interpreter in the chat-14 environment. Chat 15 MUST re-run
  both suites under HPC R to confirm no regressions — no chat-12 code
  was touched, so behavioural regression is not expected, but
  verification is needed.
- No BK schema-extraction directives (deferred from chat 13/14 to 15).
- No AS / AW / AV threshold calibration (requires HPC run first).

---

## Posture for chat 15

**Priority order:**

1. **Validate the v6 end-to-end path works** on one chromosome. Run
   STEP_C04 heavy engine + STEP_C04b classifier on LG12. Confirm all
   three per-chrom RDS land (`annot/*.ghsl_v6.annot.rds`,
   `annot/*.ghsl_v6.karyotypes.rds`,
   `per_sample/*.ghsl_v6.per_sample.rds`). Confirm
   `reg$compute$ghsl_at_candidate("LG12_17")` (or whatever a real LG12
   candidate is) returns sane per-sample rows.
2. **Unit tests.** Re-run 50/50 + 33/33 under HPC R.
3. **LG25 paralog control.** Run the same path. Expect LG25 panel data
   to look messy/diffuse — that's the null.
4. **Sweep remaining 26 chroms** once LG12 + LG25 look right.
5. **Threshold calibration** AS / AW / AV against observed
   distributions. If LG12 is weak, report it as a finding, not a
   miscalibration.
6. **BK schema-extraction.** Add `keys_extracted` directives to the
   four 4b schemas (`internal_dynamics`, `recombinant_map`,
   `internal_ancestry_composition`, nested_composition). Target is
   moving Q2 coverage from 30% toward 50-60%.

---

## Primary validation work: v6 panel smoke test

### Step 1 — classifier emits three RDS per chrom

After running STEP_C04 heavy engine on LG12, then STEP_C04b classifier:

```bash
ls -la <ghsl_dir>/annot/C_gar_LG12.ghsl_v6.*.rds
ls -la <ghsl_dir>/per_sample/C_gar_LG12.ghsl_v6.per_sample.rds
```

Sizes should be roughly:
- annot RDS: a few hundred KB (one row per window, ~1500 rows at
  LG12 scale).
- karyotypes RDS: small (a few KB) — one row per stable run,
  typically 0-200 rows per chrom at 9×.
- per-sample panel RDS: **30-80 MB compressed**, depending on chrom
  length. LG12 is mid-sized so expect ~40-60 MB.

If the panel RDS is much smaller, check the classifier logs — one
common failure mode would be `mat_data$rolling` being empty (heavy
engine didn't emit the scale). Another is the panel being only the
base columns (no div_roll_<s> columns) because the scale key loop
found nothing.

### Step 2 — loader + registry smoke test

Start an R session with the toolkit sourced:

```r
source("utils/lib_ghsl_panel.R")
panel_set_default_dir("<ghsl_dir>")

# Direct loader
panel <- load_ghsl_panel("C_gar_LG12")
nrow(panel)                                # ≈ 226 × n_windows
attr(panel, "ghsl_panel_meta")$scales_available
# expect c("s10", "s20", "s30", "s40", "s50", "s100")
grep("^div_roll_",      names(panel), value = TRUE)
grep("^rank_in_cohort_", names(panel), value = TRUE)
grep("^rank_band_",     names(panel), value = TRUE)

# Aggregate across a candidate
agg <- ghsl_panel_aggregate("C_gar_LG12",
                              start_bp = <cand_start>,
                              end_bp   = <cand_end>)
agg[order(mean)][1:10]     # most-INV_INV-like samples
agg[order(-mean)][1:10]    # most-INV_nonINV-like samples

# Registry entry point
reg <- load_registry()
reg$compute$ghsl_at_candidate("LG12_17")
```

### Step 3 — Tier-3 ghsl confirmation

`STEP_C01i_b_multi_recomb.R` should now pick up karyo_dt from the
v6 karyotypes RDS and annot_dt from the v6 annot RDS. Check the
multi_recomb log for Tier-3 activity — samples should no longer be
UNINFORMATIVE en masse. If they still are, check:
- karyotypes RDS has non-zero rows (meaning Part B actually called
  stable runs on LG12)
- `lib_ghsl_confirmation.R::ghsl_per_sample_in_interval` is being
  called with the v6-loaded karyo_dt

### Step 4 — sub-block scan validation

Pick a complex candidate (one with a `boundary` block that has
multiple soft boundaries, or a `regime_segments` block with several
segments):

```r
# Use the results-aware helper
scan <- reg$compute$ghsl_at_block_subregions("LG12_17", "regime_segments")
scan
patterns <- attr(scan, "per_sample_patterns")
patterns[n_transitions >= 1]      # samples with within-candidate band transitions
```

The `patterns` table with `n_transitions >= 1` is the within-candidate
recombinant-detection signal. If the candidate has true recombinants,
those carriers should show here.

---

## Threshold calibration plan

Once LG12 + LG25 + a handful more chroms have real distributions:

1. Plot distributions of `structure_score` (C01j), `silhouette`
   (decompose), `longest_deviation_bp` (multi_recomb), across LG12 vs
   LG25.
2. **Tune AS** — `structure_score` cutoffs currently 0.35/0.6 in
   C01j. Find the knee in the LG25 null and set the lower cutoff
   just above it.
3. **Tune AW** — `min_dco_bp` currently 200 kb in multi_recomb.
   Chat 12 called it an "educated guess"; recalibrate against LG12
   recombinant tract length distribution.
4. **Tune AV** — silhouette cutoff for `decomp_quality` flag
   currently 0.40. Check distribution of silhouette values across
   real data.
5. **Tune karyo quantiles** — `KARYO_LO = 0.15`, `KARYO_HI = 0.70`
   in the v6 classifier. These are v5-inherited. Chat-14 2e README
   flags them as "bias toward low-frequency inversions." Part C
   (k-means, no fixed quantiles) is unbiased and may be a more
   reliable primary signal at high inversion frequency — empirical
   check: correlation of Part B per-sample stable-run calls vs
   Part C interval_class assignments at candidates with frequency
   near 0.5.

---

## BK schema-extraction work

Four phase-4b schemas have no `keys_extracted` directives even though
their block data contains 51 Q2 spec keys:
- `internal_dynamics`
- `recombinant_map`
- `internal_ancestry_composition`
- `nested_composition`

Add `keys_extracted` to each (follow the pattern in
`registries/schemas/structured_block_schemas/regime_sample_dag.schema.json`
and the other chat-13-wired schemas). Use `resolve_dotted()` for any
dotted paths. Confirm with Python inspection that keys reach
`keys.tsv` post-write.

Target: Q2 coverage from 30.1% → 50-60%+.

---

## Known caveats heading into chat 15

1. **v5 back-compat branches in `run_all.R` and
   `lib_ghsl_confirmation.R`** — these accept `.ghsl_v5.*.rds` if
   `.ghsl_v6.*.rds` is missing. Intentional for transitional runs
   where some chroms have been re-classified and others haven't. Once
   the HPC sweep completes all 28 chroms with v6, we can remove the
   v5 fallback branches (cleanup for chat 16).
2. **`reg$compute$ghsl_at_block_subregions` field-name fallback
   list** — currently tries `segments`, `sub_blocks`, `subblocks`,
   `regions`, `tracts`, `windows`, `soft_boundaries_blocks`. If a
   block you care about stores sub-regions under a different key,
   extend the list in `load_compute_api()` or wrap an ad-hoc caller.
3. **Panel memoization** — `load_ghsl_panel` caches per-chrom in
   memory. A script touching 28 chroms sequentially will hold all in
   memory unless `panel_clear_cache(chrom)` is called between iters.
   ~40-60 MB per chrom × 28 chroms ≈ 1.5 GB if all retained. Not a
   problem on LANTA scratch, but worth a `panel_clear_cache` in the
   sweep loop for tidiness.
4. **Rank_in_cohort is full-cohort, not subset-relative.** If you
   query with `sample_subset = recombinants_only`, rank values are
   still relative to all 226 samples. If you want within-subset
   ranks, compute them on the `ghsl_panel_range` output.
5. **SPLIT detection still routes through run-overlap Tier-3.** The
   panel-based `ghsl_per_sample_panel_in_interval` is additive — it
   does NOT replace the run-overlap logic for SPLIT detection
   because within-sample band toggling collapses under an
   interval-mean. multi_recomb's Tier-3 calls
   `ghsl_per_sample_in_interval`, not the panel version. This is
   intentional (documented in `lib_ghsl_confirmation.R` header).
6. **The triangle_intervals CLI arg name is historical**. The file
   is produced by `STEP_D12_bridge_to_C01d.R`, not by any "triangle"
   script. The v6 classifier's `--intervals` CLI arg is still named
   after the historical upstream. Low priority; documented in
   chat-13 handoff.

---

## Chat-14 findings summary (for the record)

- **BL** (fixed) — Four consumers read v5 paths/columns; chat 14
  rewired to v6.
- **BM** (delivered) — New `utils/lib_ghsl_panel.R` backend.
- **BN** (delivered) — Registry `reg$compute$ghsl_*` methods.
- **BO** (fixed) — Heavy engine loaded all chroms eagerly; chat 14
  chunking fix.
- **BP** (done) — Heavy-engine default scales widened.
- **BQ** (done) — Part C k=4/k=5 labels.

All chat-13 items closed in chat 13; chat 14 added the above and
closed them in-session.

**BC coverage metric:** unchanged from end of chat 13
(22 / 73 = 30.1% Q2). Chat 15 BK work is the path to 50-60%.

---

## Housekeeping at end of chat 15

```bash
mv HANDOFF_PROMPT_chat15_2026-04-18.md \
   _archive/chat_history/
# plus chat-15's own audit log if one is produced
```

Create `HANDOFF_PROMPT_chat16_<date>.md`. Update `SESSION_SUMMARY.md`
+ `FIXES_APPLIED.md` in place (roll chat-15 section to top). Tarball
as `inversion-popgen-toolkit_chat15_<topic>_<date>.tar`.
