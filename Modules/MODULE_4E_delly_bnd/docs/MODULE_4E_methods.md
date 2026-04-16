# MODULE 4E — Breakend (Translocation) Structural Variant Discovery

## Methods

### Breakend discovery and genotyping (A01–A03)

Breakend (BND) variants representing translocations and complex rearrangements were called using DELLY v1.7.3 (`delly call -t BND`) on the shared duplicate-marked BAMs with the callable-region exclusion BED. Per-sample discovery was parallelized across 30 concurrent calls with 4 threads each. DELLY detects BND events from discordant read pairs mapping to different chromosomes (inter-chromosomal) or distant same-chromosome locations (intra-chromosomal), as well as split-read evidence spanning breakpoint junctions. Per-sample BCFs were merged, regenotyped across all 226 samples, and germline-filtered on the 81-sample unrelated subset.

DELLY INS calling was omitted from the pipeline because DELLY's insertion detection is unreliable for short-read data at ~5× coverage. Insertions are instead captured by Manta (MODULE_4F) for large events and Clair3 (MODULE_4A) for small indels.

### Annotation and summary (A04–A05)

Breakend events were annotated with gene overlap and classified as inter-chromosomal or intra-chromosomal. Intra-chromosomal BND pairs flanking the same genomic interval may indicate inversions, duplications, or other complex rearrangements. Summary statistics included per-sample BND burden, inter- vs intra-chromosomal ratio, and chromosome-pair frequency distributions.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| DELLY version | 1.7.3 | A02, A03 |
| SV type | BND | A02 |
| Threads per call | 4 | A02 |
| Parallel samples | 30 | A02 |
| Shared exclude BED | from MODULE_4B | A02 |
| Shared markdup BAMs | from MODULE_4B | A02, A03 |
| Full cohort | 226 samples | A02–A03 |
| Unrelated subset | 81 samples | A03 |
