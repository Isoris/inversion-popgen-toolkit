# MODULE 4G — Manta Structural Variant Discovery

## Methods

### Manta per-sample discovery (A01–A02)

Structural variants were independently called using Manta on each of the 226 duplicate-marked BAMs. Call regions were defined by inverting the callable-based exclusion BED to produce a bgzipped and tabix-indexed regions file. A custom Manta configuration was used with `minCandidateVariantSize = 50` to capture small SVs at the 50 bp boundary. Discovery was parallelized across 30 samples concurrently with 4 threads per Manta workflow. Manta detects deletions, duplications, insertions, inversions (encoded as BND pairs with INV3/INV5 tags), and breakends in a single run using local assembly of breakpoint regions.

### Inversion conversion, merging, and splitting (A03)

Manta encodes inversions as paired BND records with INV3 and INV5 tags. These were converted to standard INV records using `convertInversion_py3.py` (adapted from Manta's libexec) on each per-sample VCF prior to merging, as the conversion algorithm requires single-sample VCF IDs. Converted per-sample VCFs were merged into a 226-sample cohort VCF with `bcftools merge`. An 81-sample unrelated subset was extracted. The merged VCFs were split into per-type catalogs: DEL, DUP, small INS (fully assembled, SVLEN present, 50–~200 bp), large INS (incompletely assembled, LEFT_SVINSSEQ/RIGHT_SVINSSEQ present), INV, and BND. GT matrices and BED files were produced for each type.

### Evidence-based tiered filtering (A04)

Per-type catalogs were filtered using five evidence axes: total alt evidence (PR_alt + SR_alt summed across carriers), per-carrier minimum evidence, QUAL score, size, and carrier count. IMPRECISE variants (lacking split-read breakpoint assembly) were flagged but retained, as real SVs at ~9× coverage are frequently IMPRECISE due to insufficient reads for junction assembly. Three tiers were defined: Tier 1 (lenient, minimum evidence), Tier 2 (publication-grade, good evidence with carrier breadth and QUAL support), and Tier 3 (strict, strong evidence across many carriers). Both 226-sample and 81-sample cohorts were processed at all tiers.

### Summary reporting (A05)

Per-sample SV burden was computed across type, size class, tier, and sharing spectrum. Size distributions, chromosome-level counts, tier composition breakdowns, pairwise sharing matrices, and 1-Mb window density tables were generated. A combined multi-type summary enabled direct comparison of Manta and DELLY catalogs.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| Manta minCandidateVariantSize | 50 bp | A02 |
| Manta threads per sample | 4 | A02 |
| Parallel samples | 30 | A02 |
| Insertion classes | small (assembled) / large (partial) | A03 |
| Filter axes | PR+SR alt, per-carrier min, QUAL, size, carriers | A04 |
| Tier 1 | lenient | A04 |
| Tier 2 | publication-grade | A04 |
| Tier 3 | strict high-confidence | A04 |
| Full cohort | 226 samples | A02–A04 |
| Unrelated subset | 81 samples | A03–A04 |
