# MASTER HANDOFF PROMPT — Part B (phase 3 unpack + audit)

Copy everything from here to the end, paste into the next chat with
the tarball attached.

---

I'm continuing a multi-session code audit of my catfish hybrid
population-genomics pipeline. The attached tarball
`inversion-popgen-toolkit_phase2_phase4a_fixes_2026-04-17.tar` is a
working copy. Part A was done in the previous chat — phase 2d + 2e,
the C01f port, plus two critical late-session fixes to the precomp
NN-smoothing pipeline and the 2e→C01d wiring. Part B is this chat.

**THIS CHAT = PART B.** Do NOT do Part A work. Part A is done.

## CRITICAL ACTION ITEM BEFORE RUNNING ANYTHING

**Precomp MUST be rebuilt** for all 28 chromosomes before running
any phase 2d detection. FIX 21 (applied in Part A) changes the
precomp output schema: it now writes
`<precomp_dir>/sim_mats/<chr>.sim_mat_nn<k>.rds` files that the
live precomp was NOT producing before. Without the rebuild,
phase 2d's D02/D09 NN-persistence track stays silently dead, and
C01d's D3 score stays zero. See FIX 21 section of
`AUDIT_LOG_chat4_2026-04-17.md` for full details.

This is not Part B's job to run, but Part B should not try to
smoke-test anything in phase 2d/4a until precomp has been rebuilt.
Part B's phase-3 audit is independent and can proceed.

## Required reading before starting

Read in this order at the repo root of the extracted tarball:

1. `SESSION_SUMMARY_2026-04-17.md` — running summary of the whole
   audit series, including all chat 4 additions (FIX 15-17, 20,
   21, 22 and FINDING 18, 19)
2. `FIXES_APPLIED_2026-04-17.md` — per-fix changelog
3. `AUDIT_LOG_chat4_2026-04-17.md` — full Part A audit log
   including the late-session continuation (FIX 21, FIX 22). Pay
   attention to the "What was audited and is clean" section —
   those modules don't need re-auditing. Also note the "Revised
   understanding of the NN tree" section at the end, which
   documents the correct persistence-barcode semantics of D09
   that the original prompt didn't fully convey.
4. `_archive_superseded/cheat17_fossil_detection/README.md` — the
   circular-dependency pattern from chat 2. Part B may encounter
   similar patterns in phase 3.

Conventions (same as Part A):
- Work tree at `/home/claude/work/fixed/`. Extract tarball there.
- Comment markers: `BUGFIX 2026-04-17` for changes, `ARCHIVED
  2026-04-17` for deletions.
- Severity scale: CRASH / SILENT / STALE / DESIGN. Tag each
  finding.
- Parse-check every modified R file with
  `Rscript -e "parse(file=...)"`. R 4.3.3 and
  `r-cran-data.table` are already installed via apt in the
  container (no internet access for CRAN). Python files:
  `python3 -c "import ast; ast.parse(open('x').read())"`.
- Archive, don't delete. Move dead-limb scripts to
  `_archive_superseded/<n>/` with a README.

## Part A chat 4 summary — what happened in the two batches

**Six fixes applied and parse-checked:**

- **FIX 15 (SILENT, run_all.R Phase 8)** — D01/D08 `contrast`
  column collision. `merge()` produced `contrast.x`/`contrast.y`,
  C01d silently computed D1_block_strength without its contrast
  term. Reproduced in a fresh R session. D01's `contrast` dropped
  before the merge.
- **FIX 16 (DESIGN, run_all.R Phase 8)** — D09 NN sweep tree was
  built but never merged onto scoring table; C01d's D3 was capped
  at 0.6 because `iv$nn_birth` was always NA. Reciprocal-overlap
  join added. **NOTE: this fix was a no-op in practice until FIX
  21 was also applied, because the tree was always NULL upstream.**
- **FIX 17 (SILENT, STEP_D04_interior_cv.R)** — D04's `patchiness`
  (tail fraction) collided with D08's `patchiness` (CV). Renamed
  to `below_median_frac`. D08's CV version — which C01d's L222
  scaling is designed for — now reaches C01d.
- **FIX 20 (SILENT, STEP_C01f_hypothesis_tests.R)** — ported
  C01e's `compute_bands_for_candidate()` as a second-tier
  fallback. Three-tier chain: registry → k-means → legacy file.
  Candidates with no registry entry no longer silently dropped
  from hypothesis testing.
- **FIX 21 (CRASH/DESIGN, 2c_precomp + 00_inversion_config + 2d
  defaults)** — the live `STEP_C01a_precompute.R` did NOT produce
  `sim_mat_nn<k>.rds` files (the NN-smoothing loop existed only
  in archived snake1 precomp). Entire D02/D09 NN-persistence
  track was silently dead. Ported the MDS-space kNN smoothing
  loop; expanded scale ladder to `20,40,80,120,160,200,240,320`
  (was `20,40,80`). ACTION: rebuild precomp before running 2d.
- **FIX 22 (DESIGN, run_all.R Phase 8 + C01d D10)** — wired
  2e/C04 GHSL v5 phased-Clair3 signal into the scoring table.
  Six GHSL columns merged per block. C01d's D10 now blends
  `0.60 * d10_ghsl + 0.40 * d10_simmat` when C04 data exists,
  falls back to sim_mat alone otherwise. Emits `d10_source`
  column. Layer C of the 4-layer framework is now live.

**Two findings documented:**

- **FINDING 18 (DESIGN)** — D14 landscape classifier is orphaned
  from run_all.R. C01d correctly falls through to its shape-class
  fallback. Not a bug, dead contract. Header comment added.
- **FINDING 19 (DESIGN)** — SUPERSEDED by FIX 22. 2e/C04 is now
  wired into C01d via the D10 dual-source patch. 2e README was
  updated to reflect the new status.

**None of these touch phase 3 or MODULE_5A2.** You start with a
clean slate for Part B.

## Part B — phase 3 unpack + audit

### Structural change first (do this before auditing)

Move the contents of
`inversion_modules/phase_3_refine/MODULE_5A2_breakpoint_validation/`
one folder up, so the scripts land directly in
`inversion_modules/phase_3_refine/`. Then delete the now-empty
`MODULE_5A2_breakpoint_validation/` directory. Update:

- `phase_3_refine/README.md` (or create if missing) — describe the
  new flat layout
- Any launcher or SLURM script that references the old path
- Any downstream consumer (check
  `phase_4_postprocessing/4d_group_dependent/` —
  `breakpoint_evidence_audit.py` is a candidate) that imports
  paths from phase 3

Quentin specifically said the MODULE_5A2 wrapper is unnecessary
nesting. The module IS the whole of phase 3. Flatten it.

### Audit priorities

1. **`breakpoint_validator_standalone.py`** — ~1300 lines, the
   centerpiece. Read the whole boundary-refinement logic: reads
   DELLY/Manta VCFs + BAMs, produces bp-resolution boundary
   coordinates with concordance scores. Verify the algorithm AND
   the output schema.

2. **The numbered shell / Python scripts** 01 through 06:
   - `01_extract_inv_candidates.sh`
   - `02_extract_breakpoint_evidence.py`
   - `03_statistical_tests_and_seeds.py`
   - `04_validation_plots.py` (skim — plotting)
   - `05_delly_manta_concordance.py`
   - `06_bnd_inversion_signal.py`

   Verify each reads what the previous one writes. Cross-check
   output schemas of 03/05/06 against what
   `population_regenotype.py` (in 4d_group_dependent) reads.

3. **`annotate_population_confidence.sh`** — this writes the
   `POP_CONF_DIR/226.INV.population_confidence.tsv` that STEP_C00
   (phase 2c) reads. Audit the output column schema against
   STEP_C00 L511's read pattern. This is a cross-phase contract.

4. **`00_breakpoint_validation_config.sh` +
   `run_breakpoint_validation.sh`** — config + master launcher.
   Verify env vars propagate correctly.

### The gate question for phase 3

Does phase 3 actually feed 4a, or only 4d?

- The `phase_2_discovery/2c_precomp/README.md` claims
  `MODULE_5A2` feeds C01g (phase 4a).
- Session notes from chat 1/2 said phase 3 feeds 4d via
  `population_regenotype.py`, and the SV evidence reaching C01g
  comes independently from C00 (in 2c_precomp), not from phase 3.

Reconcile this. Grep for phase 3 output file names in all of
phase 4. If 4a doesn't read anything from phase 3, update both
READMEs (2c_precomp/README.md and the new phase_3_refine/README.md).

### Also in-scope for Part B

`classify_inversions` is mentioned as a downstream consumer in
the 2c_precomp README. Determine whether it exists anywhere in
the repo, or if it's a planned-but-never-written script. Document
the finding either way.

### End-of-chat paperwork

- Write `AUDIT_LOG_chat5_2026-04-17.md` at repo root summarizing
  what this chat found and fixed
- Update `SESSION_SUMMARY_2026-04-17.md` with chat 5 additions
- Update `FIXES_APPLIED_2026-04-17.md` with chat 5 entries
- Re-create the tarball:

    tar --format=pax -cf \
        /mnt/user-data/outputs/inversion-popgen-toolkit_phase2_phase4a_fixes_2026-04-17.tar \
        --exclude='.git' --exclude='__pycache__' \
        -C /home/claude/work fixed/

- `present_files` on the tarball

## Stop rule

Quentin overrode the stop rule for chat 4 on the grounds that the
method is architecturally complete and the remaining work is
"plumbing and thinking so hard to make sure everything is
coherent". Chat 4 proved him right, finding six real bugs
(including the architecturally-broken NN-smoothing loop in FIX 21
and the silent-candidate-drop-out in FIX 20) and wiring Layer C
into the main scoring for the first time. Similar rigour is
warranted for Part B.

After Part B the audit is closed unless a bug surfaces during
manuscript figure generation. No open-ended "let's also audit
phase 5" requests.

## Pipeline state context (updated 2026-04-17 late)

- DELLY2 and Manta SV catalogs: 100% done
- BEAGLE dosages: 100% done
- MDS (phase 2b): done
- Clair3: LG01-LG12 done, LG13-LG16 running, LG17-LG27 and
  LG02-LG06 stragglers over next ~3 days
- **2c_precomp: NEEDS REBUILD** before any new 2d runs (FIX 21
  adds nn-scaled sim_mat outputs)
- 2d_candidate_detection: algorithm fully audited in chat 4,
  FIX 15/16/17/21 applied, smoke-test pending after precomp
  rebuild
- 2e_ghsl: wired to C01d via FIX 22. Can be run per-chromosome
  as Clair3 phasing completes
- 4a: fixed in chat 2, chat 4 added FIX 22 to C01d's D10
- 4b/4c/4d/4e: column-read contract verified against 4a outputs
  in chat 2. Chat 4 added FIX 20 in C01f (4c).
- **Phase 3: NOT YET AUDITED** — this chat's scope.

## Start

After reading the four prior docs, confirm you're Part B and
start with the structural unpack (move
MODULE_5A2_breakpoint_validation/ contents up one level), then
grep for all references to the old path across the repo, then
begin the script audit starting with
`breakpoint_validator_standalone.py`.
