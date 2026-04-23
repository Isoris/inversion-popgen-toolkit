# DROP_README — pass 10: phase_6 MODULE_5E archive

**Pass:** 10
**Date:** 2026-04-24
**Format:** diff-only tarball
**Scope:** Archive `MODULE_5E_Inversion_HOBS` (zero callers, superseded
by `phase_qc_shelf/STEP_Q07b + STEP_Q07c`); remove dead `HOBSDIR`
alias from master config; update 6 README / doc files to reflect the
new Hobs-confirmation home.

Second pass of phase 3-6 reorg. Phase 4b/4c/4d + phase 5 still to come.

---

## ⚠️ Manual commands (required because diff-only)

```bash
cd /mnt/c/Users/quent/Desktop/inversion-popgen-toolkit

# Remove the original MODULE_5E directory (its contents are preserved
# in _archive_superseded/ by this tarball).
rm -rf inversion_modules/phase_6_secondary/MODULE_5E_Inversion_HOBS
```

That's it — one `rm -rf`. The tarball drops the archived copy at
`inversion_modules/_archive_superseded/MODULE_5E_Inversion_HOBS_superseded_by_Q07b/`
(with an added `SUPERSEDED_NOTE.md`), plus the documentation updates.

After `rm -rf` + drag-drop, GitHub Desktop should show: one directory
deletion (phase_6_secondary/MODULE_5E_Inversion_HOBS), one directory
addition (_archive_superseded/MODULE_5E_Inversion_HOBS_superseded_by_Q07b),
and modifications to 6 other files.

---

## Why archive 5E (not keep, not delete)

**Zero live callers.** Grep on 2026-04-24 confirmed no launcher, no
orchestrator, no `sbatch` script references `run_hobs_confirmation_module.R`
anywhere in the live tree. The `HOBSDIR` env var in
`inversion_modules/00_inversion_config.sh` is set but consumed
nowhere.

**Scientifically superseded by Q07b + Q07c.** Per-karyotype-group
Hobs (running ANGSD `-doHWE` separately on Hom1 / Het / Hom2) is the
scientifically correct Mérot approach — `Hexp` at an inversion SNP
should be computed from *within-group* allele frequencies, not a
cohort-pooled estimate. MODULE_5E's single-subset wrapper is a
narrower predecessor of the same idea. Q07b + Q07c additionally runs
a patched ANGSD (`Isoris/angsd_fixed_HWE`) and the C `hobs_windower`
binary at seven scales per group.

**Archive vs delete.** The R logic is a compact reference for the
simple single-subset Hobs approach, in case any future use case needs
a lightweight check without the per-group refinement. The module
identity (HOBSDIR / "5E") is also a recognisable pointer from chat
history and older design docs — preserving by name helps grep against
historical context.

Full side-by-side comparison lives in the new `SUPERSEDED_NOTE.md`.

---

## Why keep MODULE_5C and MODULE_5D

Both are actively wired:

- **MODULE_5C_Inversion_LD** — `run_ld_candidate_suite_v4.slurm`
  invokes `run_ld_candidate_suite_v4.py`; `compute/` has ngsLD-driver
  SLURM scripts; `plot/` has 6 plotting scripts. Self-contained.
- **MODULE_5D_Inversion_FST** — `STEP19_candidate_FST_scan.slurm`
  invokes `STEP19_candidate_FST_scan.R`. Self-contained.

Both modules produce the LD / Fst secondary confirmation signals the
manuscript relies on for §3.7. The `MODULE_5*` folder naming is
legacy, but self-relative path references inside their launchers mean
a rename needs a coordinated sweep — tracked for the HPC-coordinated
Cat 3 RENAMING pass (see `phase_2_discovery/2c_precomp/RENAMING.md`).

---

## File changes

### Moved to archive (contents preserved)

- `inversion_modules/phase_6_secondary/MODULE_5E_Inversion_HOBS/` →
  `inversion_modules/_archive_superseded/MODULE_5E_Inversion_HOBS_superseded_by_Q07b/`
  - `README.md` (original, preserved as-is)
  - `run_hobs_confirmation_module.R` (original, preserved as-is)
  - **`SUPERSEDED_NOTE.md`** (new — side-by-side comparison with Q07b+Q07c)

### Modified

- `inversion_modules/00_inversion_config.sh` — removed
  `HOBSDIR="${CODEBASE}/MODULE_5E_Inversion_HOBS"` dead alias;
  replaced with a comment explaining the archive.
- `inversion_modules/README.md` — tree diagram updated to drop 5E and
  note Hobs moved to phase_qc_shelf.
- `inversion_modules/phase_6_secondary/README.md` — rewritten:
  2 modules instead of 3, new "Hobs confirmation moved" section
  pointing at phase_qc_shelf.
- `inversion_modules/phase_5_followup/README.md` — ASCII DAG comment
  + "Feeds phase 6 secondary" section updated to reflect that
  per-candidate Hobs overlay now consumes phase_qc_shelf Q07b+Q07c
  outputs instead of feeding a phase_6 5E module.
- `Modules/MODULE_2A_snp_discovery/README.md` — "MODULE_5E" → "per-group
  Hobs confirmation (phase_qc_shelf STEP_Q07b + Q07c, previously
  MODULE_5E)".
- `Modules/MODULE_3_heterozygosity_roh/README.md` — same treatment for
  the "MODULE_5E (Hobs confirmation)" mention.
- `docs/MODULE_MAP.md` — phase_6 row, quick-lookup table row for Hobs,
  and numbering-oddity #2 all updated.

---

## Pre-ship verification

- [x] `bash -n inversion_modules/00_inversion_config.sh` passes
- [x] All markdown fences balanced on edited files (phase_6/README,
      SUPERSEDED_NOTE, MODULE_MAP, inversion_modules/README,
      phase_5/README, MODULE_2A README, MODULE_3 README)
- [x] `phase_6_secondary/` now contains only MODULE_5C + MODULE_5D + README
- [x] `_archive_superseded/MODULE_5E_Inversion_HOBS_superseded_by_Q07b/`
      contains the original 2 files + new SUPERSEDED_NOTE.md
- [x] No `MODULE_5E` refs remain in any live `.sh`/`.py`/`.slurm`/`.R`
      (only intentional comment in `00_inversion_config.sh` explaining
      the archive — that's documentation)
- [x] No dangling refs to `HOBSDIR` anywhere (grep confirmed dead env var)
- [x] No cohort-identity regressions (grep F1/hybrid clean on edited files)

---

## Commit message

```
pass 10: archive MODULE_5E (superseded by phase_qc_shelf Q07b+Q07c)

MODULE_5E_Inversion_HOBS has zero live callers in the repo and is
scientifically superseded by the per-karyotype-group Hobs pair
phase_qc_shelf/STEP_Q07b + STEP_Q07c (Engine H).

Moved:
  inversion_modules/phase_6_secondary/MODULE_5E_Inversion_HOBS/
    -> inversion_modules/_archive_superseded/
       MODULE_5E_Inversion_HOBS_superseded_by_Q07b/

Added:
  .../SUPERSEDED_NOTE.md  (side-by-side comparison + rationale)

Modified:
  inversion_modules/00_inversion_config.sh          (removed HOBSDIR)
  inversion_modules/README.md                       (tree diagram)
  inversion_modules/phase_6_secondary/README.md     (rewrite: 3->2)
  inversion_modules/phase_5_followup/README.md      (DAG + feeds note)
  Modules/MODULE_2A_snp_discovery/README.md         (prose)
  Modules/MODULE_3_heterozygosity_roh/README.md     (prose)
  docs/MODULE_MAP.md                                (row + oddity + lookup)

MODULE_5C and MODULE_5D stay as-is — both are actively wired, only
the legacy folder naming is tracked for the HPC-coordinated Cat 3
RENAMING pass.
```

---

## Phase 3-6 reorg progress

- ✅ **Pass 9**: phase_3 Layer-tagged rename
- ✅ **Pass 10** (this pass): phase_6 MODULE_5E archive
- ⏳ **phase_4b** (next candidate): 10 flat files + lib_*/engine/core scattered at top level
- ⏳ **phase_4d**: 14 flat files (cheats + bridges + regenotype + wiring)
- ⏳ **phase_5**: 2-3 chats, big reorg

From `docs/OPEN_TODOS.md`:
- 🟢 T2, T3 — BREEDING + Q7B wiring (needs prior-chat spec)
- 🔵 T7 — RENAMING Cat 2 sed pass
- 🔵 T12 — manuscript MAPQ inconsistency
- 🟡 T5, T6 — HPC-blocked (full v7 run + Engine H validation)

Pick whichever you want next.
