# `phase_6_secondary/` — secondary confirmation tests on candidate inversions

Three orthogonal secondary confirmation tests that each predict a
specific signature of a real balanced inversion polymorphism. These
are **not** discovery modules — they're confirmation layers applied
after `phase_4_postprocessing/` has produced a classified catalog and
`phase_5_followup/` has generated per-candidate data products.

Each module answers a different prediction:

| Module | Prediction tested |
|---|---|
| `MODULE_5C_Inversion_LD/` | Suppressed recombination → elevated interior LD |
| `MODULE_5D_Inversion_FST/` | Balanced polymorphism → elevated Fst between karyotype groups |
| `MODULE_5E_Inversion_HOBS/` | Heterokaryotype carriers → locally elevated Hobs vs HWE |

A candidate that passes all three tests in addition to the phase 3/4
breakpoint + boundary evidence is strongly supported as a real
polymorphic inversion in the 226-sample *C. gariepinus* cohort.

## Layout

```
phase_6_secondary/
├── MODULE_5C_Inversion_LD/           linkage disequilibrium
│   ├── README.md
│   ├── analysis/   compute/   plot/
│   ├── run_ld_candidate_suite_v4.py
│   └── run_ld_candidate_suite_v4.slurm
│
├── MODULE_5D_Inversion_FST/          Hudson Fst scans
│   ├── README.md
│   ├── STEP19_candidate_FST_scan.R
│   └── STEP19_candidate_FST_scan.slurm
│
└── MODULE_5E_Inversion_HOBS/         Hobs/HWE confirmation (Mérot-style)
    ├── README.md
    └── run_hobs_confirmation_module.R
```

## Dependencies

Each module reads:

- Candidate table + karyotype / contrast group assignments from phase 4
  (registry-backed; see `phase_4_postprocessing/4b_group_proposal/`
  and `4c_group_validation/`).
- Per-candidate dosage / BEAGLE / GL data from phase 5 engines where
  relevant (residual dosage for LD; raw dosage for Fst; windowed Hobs
  from phase 5's `analysis/04_candidate_overlay.py`).

None of the phase 6 modules depend on each other — they can run in
parallel.

## Naming note

The `MODULE_5*` prefixes persist from the pre-phase-refactor scheme
where these were direct sub-modules of MODULE_5. The folder paths are
referenced from several wrapper scripts (`run_hobs_confirmation_module.R`,
`run_ld_candidate_suite_v4.py`, etc.) via `source()`/imports; renaming
them requires a coordinated sweep and is tracked in
`inversion_modules/phase_2_discovery/2c_precomp/RENAMING.md`
Category 3 (cross-module path renames, HPC-coordinated).

## Manuscript mapping

All three modules feed into §3.7 (secondary confirmation) of
`MS_Inversions_North_african_catfish`. LD and Hobs tracks appear as
per-candidate supplementary panels; the Fst shelf contrast is the
headline karyotype-vs-flanking statistic (Fst 0.308 vs flanking
~0.032–0.055 on the LG28 test case).
