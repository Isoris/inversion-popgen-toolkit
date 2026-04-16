# MODULE 4A — SNP and Small Indel Discovery

## Methods

### Variant calling with Clair3 (D01–D03)

Single-sample SNP and small indel discovery was performed using Clair3 with the built-in Illumina model on the population-genomics-filtered BAMs. Per-chromosome BED files were generated from the reference index to scope each Clair3 run. Variant calling was parallelized as SLURM array jobs (one task per sample × chromosome combination). Key parameters included: quality threshold 20 for PASS/LowQual classification (`--qual 20`), minimum allele frequency 0.08 for both SNPs and indels (`--snp_min_af 0.08`, `--indel_min_af 0.08`) to capture weak candidates for downstream rescue, minimum mapping quality 20, and minimum coverage 2 to accommodate the ~5× sequencing depth. WhatsHap was enabled for both intermediate and final output phasing (`--use_whatshap_for_final_output_phasing`), producing phased genotypes (pipe-separated `0|1`) and phase set (PS) tags in the output VCFs. GVCF output was additionally generated for potential future joint genotyping.

Note: overlapping read-pair allele depth correction (PEAD) was not applied because overlapping mate inflation was already handled upstream by BamUtil clipOverlap in MODULE_1.

### VCF parsing and local event detection (P01–P02)

Per-sample Clair3 VCFs were parsed and annotated with variant type (SNP, insertion, deletion, MNP, complex), quality tier (PASS vs LowQual), allele depth, genotype quality, and genomic context from the reference FASTA index (P01). Local event blocks were detected by identifying spatially clustered variants within configurable distance thresholds, flagging potential complex events, multinucleotide variants, and alignment artifact clusters (P02).

### Phase-aware haplotype blocks (P02B)

Local haplotype phase blocks were constructed using two tiers of evidence. TIER_1 (WhatsHap) extracted phase set assignments directly from the Clair3 WhatsHap output, grouping co-phased variants by their PS tag. TIER_2 (read-pair rescue) extended phase information to unphased heterozygous variants by performing a region-centric single-pass BAM scan with union-find clustering, linking variants observed on the same read or read pair. Phase blocks were characterized by size, variant count, and tier of evidence. At ~5× Illumina coverage, phase blocks are typically short (2–5 variants) but locally reliable. Phase block outputs feed directly into the inversion detection pipeline (MODULE_5A) for haplotype contrast scoring.

### Strong variant rescue and weak candidate export (P03–P04)

Filtered (LowQual) variants with strong individual-level evidence were rescued using criteria including allele depth, genotype quality, and local context (P03). Remaining low-confidence indel candidates were exported as weak candidates with read-level evidence metrics for population-level validation (P04).

### Population-level regenotyping (P05–P07)

Weak indel candidates were clustered across all 226 samples by genomic position and variant signature to identify recurrent candidates (P05). A shared regenotype catalog was constructed from clustered candidates exceeding sample-count thresholds (P06). Each sample was independently scanned against the shared catalog using a linear-pass BAM reader with binary-search candidate lookup, extracting allele support at each catalog site (P07A). Per-sample regenotype results were merged and cohort-level summary statistics (allele frequency, sample support count, genotype concordance) were computed (P07B).

### Final classification and marker packaging (P08–P10)

Each variant was assigned to one of six classes based on individual quality, population support, and regenotype evidence: (1) PASS-confirmed, (2) PASS-unvalidated, (3) rescued-strong, (4) rescued-population, (5) weak-supported, (6) rejected (P08). For confirmed variants, a marker handoff package was generated containing flanking sequences, BED coordinates, and a master table suitable for multiplex assay design (P09). Per-sample publication figures summarizing variant counts, quality distributions, classification breakdown, and phase block statistics were produced (P10).

### Cohort-wide downstream analysis (X01–X21)

A population-level variant catalog was constructed by merging per-sample classified variants across all chromosomes (X01). Phase block catalogs were compiled (X02). Variants were annotated with gene, exon, and CDS overlap using annotation BEDs (X03). Per-sample burden summaries, pairwise sample distance matrices, marker selection optimization, and gene/chromosome summary tables were computed (X04–X07). Variant effect prediction and breeding concern scoring were performed using SnpEff-style consequence annotation and a five-level breeding concern scale (BC0–BC4) based on variant type, location, predicted functional impact, and population frequency (X20–X21).

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| Clair3 platform | ilmn | D02 |
| Clair3 --qual | 20 | D02 |
| Clair3 --snp_min_af | 0.08 | D02 |
| Clair3 --indel_min_af | 0.08 | D02 |
| Clair3 --min_mq | 20 | D02 |
| Clair3 --min_coverage | 2 | D02 |
| WhatsHap phasing | enabled (intermediate + final) | D02, P02B |
| GVCF output | enabled | D02 |
| Phase block TIER_1 | WhatsHap PS tags | P02B |
| Phase block TIER_2 | read-pair union-find | P02B |
| Classification classes | 6 | P08 |
| Breeding concern scale | BC0–BC4 | X21 |
