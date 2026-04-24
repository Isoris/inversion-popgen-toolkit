# MODULE_5D — Inversion FST

Hudson FST scans from dosage for inversion candidate contrasts.

## Dependencies

- **Upstream:** MODULE_5A STEP08 (dosage + sites), STEP10 (candidate table)
- **Grouping:** MODULE_5B STEP17c (contrast group manifests)

## Script

`STEP19_candidate_FST_scan.R` (v2)

Computes windowed Hudson FST for three contrasts per candidate:
- FULL_A vs FULL_B (main inversion contrast)
- FULL_A vs HALF
- FULL_B vs HALF

Uses the same STEP17c-exported contrast groups as the LD module.

## Usage

```bash
# All candidates
Rscript STEP19_candidate_FST_scan.R <group_dir> <candidate_table> <dosage_dir> <outdir>

# Single candidate
Rscript STEP19_candidate_FST_scan.R <group_dir> <candidate_table> <dosage_dir> <outdir> 42

# Via SLURM
sbatch STEP19_candidate_FST_scan.slurm
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| flank_frac | 0.10 | Fraction of candidate length added as flanking region |
| win_snps | 50 | SNPs per FST window |
| min_per_group | 5 | Minimum samples per group to compute FST |
