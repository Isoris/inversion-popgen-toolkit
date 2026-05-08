# `2a_local_pca/` — dosage + per-chromosome local PCA windows

First block of phase 2. Converts ANGSD BEAGLE genotype likelihoods into
per-chromosome dosage files, then computes local-PCA summaries in
sliding windows along each chromosome. Feeds `2b_mds/`.

## Layout

```
2a_local_pca/
├── STEP_A01_beagle_to_dosage_by_chr.py    # BEAGLE → per-chr dosage + sites
├── STEP_A02_local_pca_windows_by_chr.R    # per-chr windowed local PCA (legacy monolithic)
├── STEP_A03_dense_registry_stage1.R       # per-chr dense-window PCA (parallel stage 1)
├── STEP_A03_dense_registry_stage2.R       # merge per-chr outputs into master registry
└── LAUNCH_A0{1,2,3}*.slurm                # SLURM wrappers (array-per-chromosome)
```

## Workflow

```
ANGSD .beagle.gz
        │
        ▼  STEP_A01 (python): dosage = P(AB) + 2·P(BB)
        │         optional --gl-triplet keeps P(AA|AB|BB) vectors
        │
        ▼  $DOSAGE_DIR/<chr>.dosage.tsv.gz + <chr>.sites.tsv.gz
        │
        ├──────────────────────────────┐
        ▼                              ▼
STEP_A02 (legacy monolithic)   STEP_A03 stage1 (SLURM array, 1 chr/task)
per-chr windowed PCA           per-chr windowed PCA → tmp/<chr>.*
$outprefix.window_pca.rds                       │
                                               ▼
                                       STEP_A03 stage2 (merge)
                                       → master window_pca.rds + pca_table.tsv.gz
```

## Parameters (defaults)

- Window: 100 SNPs, step 20 SNPs (dense registry); A02 is similar.
- `--npc 2`: eigenvectors kept per window.
- Covariance matrix is over **individuals**, not sites — window is an n×n matrix.

## Downstream

`2b_mds/` reads `window_pca.rds` (A02) or per-chr RDS dir (A03-stage2) to
compute lostruct distances + MDS.

## Legacy vs parallel

- A02 = monolithic v7.4 (single node, all chromosomes sequentially).
- A03 = v8.3-parallel (SLURM array stage1 + merge stage2). Preferred.
- Both produce the same downstream format.
