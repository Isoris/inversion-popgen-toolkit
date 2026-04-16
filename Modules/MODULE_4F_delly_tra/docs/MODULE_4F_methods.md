# MODULE 4F — Translocation Structural Variant Discovery

## Methods

### Translocation discovery and genotyping (A01–A03)

Inter-chromosomal translocations were called using DELLY v1.7.3 (`delly call -t TRA`) on the shared duplicate-marked BAMs with the callable-region exclusion BED. Per-sample discovery was parallelized across 30 concurrent calls with 4 threads each. DELLY detects translocations from discordant read pairs mapping to different chromosomes and split-read evidence spanning inter-chromosomal breakpoint junctions. Per-sample BCFs were merged into a unified site list, all 226 samples were regenotyped at merged sites, and the cohort BCF was produced. An 81-sample unrelated subset was extracted and germline-filtered.

### Annotation and summary (A04–A05)

Translocation events were annotated with gene overlap at both breakpoint ends and classified by chromosome pair. Summary statistics included per-sample TRA burden, chromosome-pair frequency distributions, and recurrent translocation hotspots.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| DELLY version | 1.7.3 | A02, A03 |
| SV type | TRA | A02 |
| Threads per call | 4 | A02 |
| Parallel samples | 30 | A02 |
| Shared exclude BED | from MODULE_4B | A02 |
| Shared markdup BAMs | from MODULE_4B | A02, A03 |
| Full cohort | 226 samples | A02–A03 |
| Unrelated subset | 81 samples | A03 |
