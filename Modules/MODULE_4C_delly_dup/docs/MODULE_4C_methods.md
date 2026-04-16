# MODULE 4C — Duplication Structural Variant Discovery

## Methods

### Duplication discovery and genotyping (A01–A03)

Tandem duplications were called using DELLY v1.7.3 (`delly call -t DUP`) on the shared duplicate-marked BAMs (produced by MODULE_4B) with the callable-region exclusion BED. Per-sample discovery was parallelized across 30 concurrent calls with 4 threads each. Per-sample BCFs were merged into a unified site list (`delly merge`), all 226 samples were regenotyped at merged sites (`delly call -v`), and the cohort BCF was produced with `bcftools merge`. An 81-sample unrelated subset was extracted and germline-filtered with `delly filter -f germline`.

### Annotation and downstream analysis (A04–A05, B01–B02)

Duplications were annotated with gene, exon, and CDS overlap, repeat content, and depth support using the same annotation framework as MODULE_4B. A master annotation table was constructed merging all annotation layers with carrier counts, frequency class, and size class. Per-sample DUP burden summaries and tiered marker selections were produced. Presence/absence tables were generated for population-level DUP sharing analysis.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| DELLY version | 1.7.3 | A02, A03 |
| SV type | DUP | A02 |
| Threads per call | 4 | A02 |
| Parallel samples | 30 | A02 |
| Shared exclude BED | from MODULE_4B | A02 |
| Shared markdup BAMs | from MODULE_4B | A02, A03 |
| Strict marker QUAL | > 500 | B01 |
| Strict marker PE | > 4 | B01 |
| Full cohort | 226 samples | A02–A03 |
| Unrelated subset | 81 samples | A03 |
