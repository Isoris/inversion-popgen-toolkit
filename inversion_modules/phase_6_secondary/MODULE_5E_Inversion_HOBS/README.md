# MODULE_5E — Inversion HOBS

Secondary confirmation module using observed heterozygosity (Hobs) and Hardy-Weinberg equilibrium (HWE) sliding windows.

## Important

This is **NOT primary inversion discovery**. Use only as secondary confirmation / interpretation / support after candidates have been identified by MODULE_5A and characterized by MODULE_5B.

## Dependencies

- **External input:** ANGSD `-doHWE` output (tab-separated with header)
- **Optional:** MODULE_5A STEP10 candidate table for overlay

## Script

`run_hobs_confirmation_module.R`

Adapted from Claire Mérot's Hobs/HWE sliding-window logic, made subset-aware and multiscale for structured hatchery catfish data.

Computes:
1. Site-level Hexp, Hobs, F from ANGSD -doHWE output
2. Multiscale window summaries with mean, median, outlier burden
3. Candidate interval overlays (if candidate table provided)

## Usage

```bash
Rscript STEP_HOBS_confirmation_module.R \
    <hwe_file> <subset_label> <outdir> \
    [candidate_table=NONE] [chrom_col=Chromo] [pos_col=Position] \
    [freq_col=Freq] [f_col=F] [hwefreq_col=hweFreq]
```

## Citation

Hobs/HWE sliding-window logic: Claire Mérot et al.
ANGSD/thetaStat: Korneliussen et al.
