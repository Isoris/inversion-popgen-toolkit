# MODULE 4D — Inversion Structural Variant Discovery

## Methods

### Inversion discovery and genotyping (A01–A03)

Inversions were called using DELLY v1.7.3 (`delly call -t INV`) on the shared duplicate-marked BAMs with the callable-region exclusion BED. Per-sample discovery was parallelized across 30 concurrent calls with 4 threads each. DELLY detects inversions from discordant read-pair orientations and split-read alignments spanning inversion breakpoints. Per-sample BCFs were merged into a unified site list, all 226 samples were regenotyped at merged sites, and the cohort BCF was produced. An 81-sample unrelated subset was extracted and germline-filtered.

### Annotation and summary (A04–A05)

Inversions were annotated with gene, exon, and CDS overlap, repeat content, and depth support. Summary statistics were computed including size distribution, per-sample inversion burden, and sharing spectrum.

### Cross-validation with local PCA inversion detection (downstream)

The DELLY INV catalog provides independent SV-caller-based evidence for inversions that complements the local PCA approach in MODULE_5A. Inversions detected by both methods receive higher confidence in the MODULE_5A hypothesis testing framework, where DELLY INV genotypes serve as "Cheat 9" concordance evidence — an independent confirmation that the three-class genotype structure observed in local PCA is reproducible by a paired-end/split-read method.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| DELLY version | 1.7.3 | A02, A03 |
| SV type | INV | A02 |
| Threads per call | 4 | A02 |
| Parallel samples | 30 | A02 |
| Shared exclude BED | from MODULE_4B | A02 |
| Shared markdup BAMs | from MODULE_4B | A02, A03 |
| Full cohort | 226 samples | A02–A03 |
| Unrelated subset | 81 samples | A03 |
