# HANDOFF PROMPT — next chat, starting 2026-04-17 (late)

Paste this at the start of the next chat, with the tarball attached.

---

I'm continuing a multi-session audit of my catfish hybrid
population-genomics pipeline. The attached tarball is a working copy
with 14 fixes from the previous session already applied to
`phase_2_discovery/2c_precomp/` and `phase_4_postprocessing/4a_existence_layers/`.

**Read first, in this order:**

1. `SESSION_SUMMARY_2026-04-17.md` at the repo root — explains the 14 fixes
   from the last session and what was deliberately not touched.
2. `FIXES_APPLIED_2026-04-17.md` — the per-fix change log.
3. `_archive_superseded/cheat17_fossil_detection/README.md` — explains
   why Cheat 17 fossil detection was archived (circular dependency with
   C01d that was never orchestrated as the required two-pass run).

**Scope for this session, in priority order:**

1. **Phase 2e audit** (`phase_2_discovery/2e_ghsl/`) — just 2 files:
   - `STEP_C04_snake3_ghsl_v5.R` (32K, ~700 lines) — GHSL v5 detector
   - `STEP_C04b_snake3_ghsl_figures.R` (21K, ~500 lines) — figures

   Verify: CLI contract, what it reads from upstream (Clair3 phased
   genotypes), what it writes, and whether any downstream consumer
   actually reads its output (4a/C01d currently does NOT — Layer C
   in the 4-layer framework is documented but C01d's D10 reads
   sim_mat-based partition_stability from 2d/STEP_D05 instead).
   This is worth clarifying: is 2e output consumed anywhere, or is
   it a dead limb?

2. **Phase 3 audit** (`phase_3_refine/MODULE_5A2_breakpoint_validation/`) —
   ~11 scripts:
   - `00_breakpoint_validation_config.sh`
   - `01_extract_inv_candidates.sh` through `06_bnd_inversion_signal.py`
   - `breakpoint_validator_standalone.py` (~1300 lines — the main script)
   - `annotate_population_confidence.sh`
   - `run_breakpoint_validation.sh`
   - `README.md`

   Verify: this module feeds 4d (not 4a). Confirm
   `population_regenotype.py` in 4d reads what this module writes.
   Also check if `classify_inversions` (mentioned in 2c_precomp/README
   as a downstream consumer) actually exists anywhere.

3. **Apply fixes** to anything genuinely broken. Same severity scale:
   - CRASH: runtime error, guaranteed
   - SILENT: produces wrong data without failing
   - STALE: dead flag / obsolete doc
   - DESIGN: architectural note

4. **One known leftover from last session** — found at the very end of
   that chat but not fixed: C01f (`phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R`)
   has the same `triangle_sample_composition.tsv.gz` read pattern that
   was fixed in C01e last session (BUG 3). In C01f it's partially
   papered over by a registry-first fallback — the triangle file is
   only consulted if the C01i registry groups are empty. Worth
   porting the same `compute_bands_for_candidate()` helper from C01e
   to C01f as a second fallback. Not a blocker but closes the loop.

**Out of scope / don't touch:**

- **Phase 4b/4c/4d/4e internals beyond what's needed for the phase 2e
  and phase 3 scopes.** I audited the column-read contracts last
  session; diving into their internals is a separate project.
- **C01a internals** (1875 lines of inv_likeness / test_07 Beta /
  test_05 F_ST scan math). Output schema was verified last session.
  Internal science audit is a separate project.
- **The flashlight → sv_prior rename in C01d's output column name**
  (`cheat5_family_fst_ratio` stayed legacy-named to avoid cascading
  changes to 4e/compute_candidate_status.R registry keys and
  test_registry_sanity.py). Leave it legacy-named.

**My actual priority is not more code auditing, it's the manuscript.**
Previous chats pointed this out and I agreed. Push back if I slide into
polishing code when the manuscript is what needs writing. The compute
pipeline is running in parallel — Clair3 LG13–LG16 tonight, then
LG17–LG20, LG21–LG24, LG25–LG27 + LG02–LG06 stragglers over the
following days. C00 sv_prior rebuild is also queued for tomorrow. I'm
using the compute-wait time to finish phase 2e and phase 3 audits, not
to optimize crufty code.

**Pipeline state:**
- DELLY2 SV catalogs: 100% done (all 28 chr)
- Manta SV catalogs: 100% done
- angsd/BEAGLE dosages: 100% done
- MDS (phase 2b): done
- Clair3 phased VCFs: LG01–LG12 done; LG13+ running in background
- 2c_precomp (C00 + C01a + C01b_1 + PHASE_01C): mostly done with
  existing fixes applied last session
- 2d_candidate_detection (staircase): ran previously
- 2e_ghsl: has not run — needs Clair3 complete first
- 4a (C01d + C01e + C01g): fixed last session, smoke test scheduled
  for tomorrow morning

**Workflow convention (same as every session):**
- Extract working copy to `/home/claude/work/fixed/`, apply fixes in
  place, tag each change with `BUGFIX 2026-04-17` or
  `ARCHIVED 2026-04-17` comment markers
- Parse-check every modified R file with `Rscript -e "parse(...)"`
  before declaring done
- At end of session: update SESSION_SUMMARY_2026-04-17.md with the
  new fixes, produce a new tarball
- If finding genuinely dead code or circular dependencies — archive to
  `_archive_superseded/` with a README documenting why

Begin with the phase 2e audit. Start by reading the two script headers
end-to-end before anything else.
