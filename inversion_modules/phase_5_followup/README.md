# `phase_5_followup/` — per-candidate deep analysis

After `phase_4_postprocessing/` has produced a classified candidate
catalog, phase 5 runs per-candidate deep-dive analyses: residual dosage,
rare-allele sharing, Hobs/HWE overlays, FST/dXY tracks, LD panels, and
the sample-structure-first STEP20–41 suite from the legacy MODULE_5B.

Where phase 4 asks "is this a real inversion and what kind?", phase 5
asks "for each candidate we're keeping, what do its samples look like
from every angle we can throw at it?"

## Layout

```
phase_5_followup/
├── LEGACY_FOLLOWUP_BRIDGE.md         integration map: STEP20–41 → v8.5
├── analysis/                         lightweight python analysis scripts
│   └── 04_candidate_overlay.py       Hobs/HWE × candidate overlay
│
├── engines/                          C binaries (high-performance compute)
│   ├── export_q_residual_dosage.c    Q-corrected residual dosage from BEAGLE
│   └── rare_sfs_pairwise.c           pairwise rare-allele sharing (Bergström-style)
│   + Makefile
│
├── figures/                          R plotting scripts (publication figures)
│   ├── 05_plot_hobs_hwe.R             genome tracks / heatmaps / outlier / candidate stacks
│   ├── 06_plot_fst_dxy_tracks.R       Mérot-style multi-track windowed stats
│   ├── 07_plot_rare_sfs_heatmap.R     Bergström Science Fig-2 style
│   ├── plot_ld_panels.R               6+3 / 7-panel LD figure suites
│   └── plot_local_Q_diagnostics.R     local Q behavior across candidates
│
├── launchers/                        SLURM wrappers for the C engines
│   ├── LAUNCH_q_residual_dosage.slurm
│   └── LAUNCH_rare_sfs_pairwise.slurm
│
├── utils/
│   └── export_module5b.py             legacy MODULE_5B export helpers
│
└── codebase_8.5_experimental/        STEP20–41 suite (sample-structure-first)
    ├── current_followup/              STEP20–41 (live — 25 R scripts)
    ├── legacy_followup/               STEP14–18 (frozen reference)
    ├── config_inversion_followup.R
    ├── run_candidate_followup.slurm
    ├── run_candidate_followup_v6.sh
    ├── README.md                      describes the STEP20–41 method
    └── DIAGNOSIS.md                   bug audit + fixes for the figure suite
```

## Workflow

Phase 5 is **not** a strict DAG — it's a menu of per-candidate
analyses, each with its own inputs and outputs. Run what you need per
candidate:

```
phase_4 candidate catalog (candidate_status.tsv + contrast groups)
        │
        ├─────────────────────────────────────┬──────────────────────┐
        ▼                                     ▼                      ▼
  engines/export_q_residual_dosage    engines/rare_sfs_pairwise   STEP20–41
  → residual BEAGLE per cand          → pairwise sharing matrix   → figure suite
        │                                     │                      │
        ▼                                     ▼                      ▼
  figures/plot_ld_panels.R          figures/07_plot_rare_sfs_heatmap   per-cand
  figures/06_plot_fst_dxy_tracks.R                                    composite
  figures/plot_local_Q_diagnostics                                    figures
                                                                      (STEP38)
        │
        ▼
  analysis/04_candidate_overlay.py
  figures/05_plot_hobs_hwe.R
  (consumed by phase_qc_shelf Q07b+Q07c outputs; the Hobs confirmation
   formerly in phase_6_secondary/MODULE_5E was archived in April 2026)
```

## Relationship to other phases

- **Reads from phase 4**: `candidate_status.tsv` (Axis 1–3 tiering),
  contrast group manifests (STEP17c equivalent), registry-backed
  candidate coordinates.
- **Shares BEAGLE inputs with phase 1**: residual dosage computation
  re-reads the BEAGLE from phase_1 outputs.
- **Feeds phase 6 secondary**: the LD and Fst modules read
  per-candidate data from phase_5 compute products. The Hobs
  confirmation path moved to `phase_qc_shelf/STEP_Q07b + STEP_Q07c`
  (per-group Hobs, superseded `MODULE_5E_Inversion_HOBS` in April 2026).
  `analysis/04_candidate_overlay.py` + `05_plot_hobs_hwe.R` still
  serve as the per-candidate Hobs overlay / plotter against the
  phase_qc_shelf outputs.

## Current vs legacy

- **`current_followup/STEP20–STEP41`**: live publication-figure
  pipeline. STEP38 (composite figure, 537 lines) is the manuscript
  panel generator.
- **`legacy_followup/STEP14–STEP18`**: retained as reference
  implementations. Do not run; `LEGACY_FOLLOWUP_BRIDGE.md` documents
  which legacy step was superseded by which current step.

## Notes

- The `codebase_8.5_experimental/` name is a misnomer — nothing here is
  "experimental" anymore. The name persists from the v8.5 refactor when
  these steps were being triaged. Consider renaming after the
  manuscript is submitted.
- Several filenames carry `snake` / `MODULE_5B` legacy vocabulary.
  These are tracked in `inversion_modules/phase_2_discovery/2c_precomp/RENAMING.md`
  Category 2 (code-identifier renames), pending a coordinated sed pass.
