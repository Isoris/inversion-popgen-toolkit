# Superseded snapshots

Old packaging tarballs that were sitting loose in the live tree. Kept
here only as insurance against "did we lose something during the rename?"
— the live tree has newer versions of every file each one contains.

| Tarball | Dated | What | Why it's superseded |
|---|---|---|---|
| `phase_2_discovery_2026-04-18.tar` | 2026-04-18 | Snapshot of `inversion_modules/phase_2_discovery/` taken during chat 17 | Live tree files are newer. Example: `STEP_D01_staircase_boundaries.R` is at v9.6 in live tree, v9.3 in this snapshot. 3 files appear in the tarball under their pre-rename names (`STEP_D01_staircase_boundaries_v9.4_pmin_fix.R`, etc.); equivalents under current names live in `_archive_old_names/` in the active tree. |
| `2d_candidate_detection_2026-04-18.tar` | 2026-04-18 | Snapshot of the staircase-detector subfolder | Live tree has 5 files the tarball doesn't (`STEP_D15_karyogram.R`, `_archive_old_names/`, etc.). |
| `phase_qc_shelf_pre_v3.0_2026-04-19.tar.gz` | 2026-04-19 | Pre-v3.0 snapshot of phase_qc_shelf | Live tree is at v3.1, with Q07b, Q07c, Q09, Q09b, Q10 additions and registry wiring that are absent from this snapshot. |

If you ever want to check "was this script different a week ago?" these
are the reference. Otherwise they can be deleted whenever you like.
