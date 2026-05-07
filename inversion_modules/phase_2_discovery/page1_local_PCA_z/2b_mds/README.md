# `2b_mds/` — lostruct distances + MDS, per-chromosome

Second block of phase 2. Takes local-PCA window summaries from `2a_local_pca/`
and produces per-chromosome MDS coordinates plus z-outlier candidate
regions. Feeds `2c_precomp/`.

## Layout

```
2b_mds/
├── STEP_B01_legacy_mds_perchr.R       # monolithic per-chr MDS (v7.4)
├── STEP_B01_mds_multimode.R           # full-genome multimode MDS (monolithic)
├── STEP_B01_mds_multimode_stage1.R    # per-focal-chr MDS (parallel stage 1)
├── STEP_B01_mds_multimode_stage2.R    # merge per-chr results → final .mds.rds
└── LAUNCH_B01_*.slurm                 # SLURM wrappers (legacy + stage1 + stage2)
```

## Workflow

```
2a window_pca.rds (or per-chr RDS registry)
        │
        ▼  per-chromosome lostruct-style window distances
        │         (within-chromosome only; cross-chr distances are meaningless
        │          for inversion detection and explode O(n²))
        │
        ▼  cmdscale MDS per chromosome
        │
        ▼  merge → .mds.rds, .window_mds.tsv.gz,
                   .candidate_regions.tsv.gz,
                   .candidate_window_membership.tsv.gz,
                   .mds_mode_metadata.tsv
        │
        ▼  2c_precomp/STEP_C01a_precompute.R (reads .mds.rds)
```

## MDS modes

- `per_chr`: only the focal chromosome's windows (baseline).
- `chunked_Nx`: focal + N·|focal| background windows sampled from other
  chromosomes. Sharpens outlier contrast. `chunked_4x` is the typical
  setting. Background sampling is logged in `mds_background_<chr>.txt`.

## Candidate regions

Outlier windows (z-score ≥ threshold, default 3) are clustered into
candidate regions with `gap_bp` (default 500 kb) and `min_windows`
(default 3). These are **MDS outlier seeds** — not final candidates —
and feed `2c_precomp/STEP_C01b_1_seeded_regions.R` for region growing.

## Legacy vs parallel

- `legacy_mds_perchr.R` / `mds_multimode.R`: monolithic (single job).
- `_stage1.R` + `_stage2.R`: v8.3 SLURM-array parallel. Preferred.
- Output structure is identical across variants — downstream readers
  (`STEP_C01a_precompute.R`, all `snake` / region scripts) expect
  `mds_obj$per_chr` in the merged .mds.rds.
