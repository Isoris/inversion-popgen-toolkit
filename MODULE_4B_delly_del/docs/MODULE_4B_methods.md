# MODULE 4B — Deletion Structural Variant Discovery

## Methods

### Input preparation and duplicate marking (A01)

DELLY requires unfiltered BAMs with duplicate marking. Raw minimap2 alignments (MODULE_1 manifest column 2, prior to population-genomics filtering) were duplicate-marked using `samtools markdup` to produce per-sample markdup BAMs. An exclusion BED was constructed by identifying 50-kb genomic bins with fewer than 500 callable base pairs (from mosdepth coverage data) and merging adjacent uncallable bins into blocks of at least 50 kb. Additionally, the first and last 50 kb of every chromosome were unconditionally masked to exclude telomeric and centromeric artifacts. Sample lists for the full cohort (226 samples) and the unrelated subset (81 samples after NAToRA first-degree pruning) were generated.

### Deletion discovery (A02)

Per-sample deletion calling was performed using DELLY v1.7.3 (`delly call -t DEL`) with the exclusion BED, running 30 samples in parallel with 4 threads per call on a full compute node (120 cores utilized). DELLY uses paired-end mapping and split-read evidence to detect deletions, producing per-sample BCF files.

### Site merging, regenotyping, and filtering (A03)

Per-sample discovery BCFs were merged into a unified site list using `delly merge`. All 226 samples were then regenotyped at the merged site set using `delly call -v` (genotyping mode). Regenotyped BCFs were merged into a cohort-level BCF with `bcftools merge`. An 81-sample unrelated subset was extracted, and germline filtering was applied using `delly filter -f germline` to produce the final population-level DEL catalog. Final outputs included VCFs, GT matrices (genotype per sample per site), and BED files for both the full 226-sample and filtered 81-sample catalogs.

### Functional and quality annotation (A04)

Deletions were annotated with functional overlap using sorted BED files for genes, exons, CDS regions, and repeats via `bedtools intersect`. Depth support was assessed using mosdepth coverage data (500-bp windows, MAPQ ≥ 30). Mate-distance quality control flagged deletions with anomalous SVLEN relative to insert size distribution, using thresholds at 20 kb (warning), 50 kb (suspicious), and 100 kb (extreme).

### Summary statistics and reporting (A05)

Per-sample deletion burden was computed across size classes (small < 1 kb, medium 1–10 kb, large > 10 kb), functional categories (repeat vs non-repeat, CDS vs exonic vs intronic), and sharing spectrum (private, rare, common). Cohort-level summary tables were generated including total counts, size distributions, and quality metric distributions.

### Publication figures (B01–B03)

A publication figure suite was generated including: genome-wide deletion heatmap across all chromosomes and samples, per-sample burden barplots, PCA from the deletion GT matrix, pairwise sample sharing heatmap, size distribution histograms, and per-chromosome deletion density. Extended plots included ancestry-grouped overlays using MODULE_2B labels, detailed quality breakdowns, and marker selection summaries.

### Downstream analysis (B02)

A master annotation table was constructed by merging VCF-derived variant information with functional class, repeat overlap, depth support, mate-distance QC, carrier counts, frequency class, and size class. Per-sample summary tables, tiered marker selections (Tier 1: PASS + PRECISE + QUAL > 500 + PE > 4; Tier 2: relaxed; Tier 3: all germline), and gene/chromosome overlap tables were produced.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| DELLY version | 1.7.3 | A02, A03 |
| SV type | DEL | A02 |
| Threads per call | 4 | A02 |
| Parallel samples | 30 | A02 |
| Chromosome-end mask | 50 kb | A01 |
| Excl min callable_bp | 500 per 50-kb bin | A01 |
| Excl min block size | 50 kb | A01 |
| Depth QC window | 500 bp | A04 |
| Depth MAPQ | 30 | A04 |
| Mate warn threshold | 20 kb | A04 |
| Mate suspicious threshold | 50 kb | A04 |
| Strict marker QUAL | > 500 | B02 |
| Strict marker PE | > 4 | B02 |
| Full cohort | 226 samples | A02–A03 |
| Unrelated subset | 81 samples | A03 |
