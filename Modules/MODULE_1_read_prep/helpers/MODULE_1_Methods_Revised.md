# MODULE 1 — Read Preparation and BAM Processing: Methods (Revised)

## 1.1 Read quality control and adapter trimming

Raw paired-end Illumina sequencing reads (PE[TODO: 150 or 151 — check fastp log "read1_mean_length"]) from 226 *Clarias gariepinus* broodstock individuals sampled from a commercial Thai hatchery that maintains broodstock for F₁ hybrid catfish seed production (*C. macrocephalus* ♀ × *C. gariepinus* ♂) were generated on [TODO: platform, e.g. Illumina NovaSeq 6000] at [TODO: facility], producing approximately [TODO: X] Gbp of raw data across [TODO: Y] run-lane units (range [TODO: min]–[TODO: max] per sample, median [TODO: N]). Reads were processed with fastp [TODO: version] (Chen et al., 2018) using the following parameters: automatic adapter detection for paired-end reads, poly-G and poly-X 3′ tail trimming (consistent with two-colour chemistry error profiles on [TODO: NovaSeq/NextSeq]), sliding-window quality trimming from both the 5′ and 3′ ends, a minimum Phred-scaled base quality threshold of 20 (fastp -q 20), hard trimming of 5 bp from both ends of both reads (fastp -f 5 -F 5 -t 5 -T 5), zero tolerance for ambiguous bases (fastp -n 0), and a minimum post-trimming read length of 30 bp (fastp -l 30). The 5 bp end-trimming was applied to remove residual ligation and end-repair artefacts commonly observed in the first and last bases of Illumina reads, which showed elevated error rates in per-base quality profiles during preliminary QC. Each sample was processed using 8 fastp threads with up to 15 concurrent jobs on a single 127-core compute node. [TODO: N] run-lane units with corrupted gzip outputs (unexpected EOF errors) were identified and re-processed from original raw reads (S01–S02). After trimming, the mean Q30 rate was [TODO: X]% (range [TODO: min]–[TODO: max]%), and [TODO: X]% of read pairs passed all quality filters. Complete per-sample QC metrics are reported in Supplementary Table S[TODO: N].

## 1.2 Parentage species verification

Species-level identity of each sequencing run-lane unit was verified using Mash distance-based screening against the *C. gariepinus* (Gar) and *C. macrocephalus* (Mac) haplotype-resolved reference subgenomes (Ondov et al., 2016). Although all 226 samples were labelled by the supplying hatchery as *C. gariepinus*, commercial hatcheries that produce F₁ *C. macrocephalus* × *C. gariepinus* hybrid seed maintain both parental broodstocks on site, and several morphologically cryptic *Clarias* species (notably *C. batrachus*) co-occur in the region. A read-level species-identity QC step was therefore essential to detect sample mix-ups, contamination with non-target *Clarias* broodstock, or non-catfish contamination before downstream analyses.

Mash sketches were built for both parental subgenomes using a k-mer size of 31, sketch size of 50,000, and minimum k-mer copy count of 2. Each run-lane unit's fastp-trimmed read pair was sketched in read mode with the same parameters, and Mash distances were computed against both subgenome sketches. For classification purposes, a run-lane was recorded as Gar-closer if dist(Gar) + 0.002 < dist(Mac), Mac-closer if dist(Mac) + 0.002 < dist(Gar), and Ambiguous otherwise, using a margin threshold of 0.002. Per-sample consensus calls were derived by majority vote across all run-lane units belonging to the same biological sample. [TODO: report actual pattern — expected result for a pure *C. gariepinus* cohort is that [TODO: X of 226] samples were classified as Gar-closer.]

An independent verification using meryl k-mer species-specific sets (k = 21, Rhie et al., 2020) with a log₂ ratio threshold of |log₂ ratio| > 1.0 showed [TODO: X]% concordance with the Mash results. Species assignment results were visualised as scatter plots of median Mash distance to Gar versus Mac per sample (Supplementary Figure S[TODO: N]).

## 1.3 Read alignment

Trimmed reads were aligned to the *C. gariepinus* (Gar) subgenome reference (fClaHyb_Gar_LG, 28 pseudochromosomes, 963,905,721 bp — the *gariepinus* haplotype extracted from the haplotype-resolved F₁ hybrid assembly described in Section 1 of the manuscript) using minimap2 [TODO: version] (Li, 2018) in short-read mode with the following flags: -ax sr (short-read preset), -Y (soft clipping for supplementary alignments), -c (output CIGAR), --secondary-seq (output sequence for secondary alignments), and --sam-hit-only (omit unmapped reads from SAM output). Read group information was set per run-lane unit with SM (sample), LB (library, set to sample name), PU (platform unit, set to run.lane), and PL (platform, set to ILLUMINA). Alignments were coordinate-sorted with samtools [TODO: version] (Danecek et al., 2021) sort and indexed with samtools index, using 16 threads per sample. The mean mapping rate across the cohort was [TODO: X]% (range [TODO: min]–[TODO: max]%).

## 1.4 BAM merging, duplicate marking, and overlap clipping

Per-run BAMs were merged into per-sample BAMs using samtools merge, preserving all read group information. Prior to merging, RG/SM tag consistency was validated: the SM field of every @RG header line in each input BAM was checked to match the expected sample name, and any mismatch triggered an abort. The merged BAM was then processed through the following pipeline: (1) name-sort with samtools sort -n; (2) mate information repair with samtools fixmate -m (adding mate score tags required for duplicate marking); (3) coordinate-sort with samtools sort; (4) duplicate marking (mark only, not removal) with samtools markdup -s (writing statistics to stderr); and (5) overlap clipping with BamUtil [TODO: version] (Jun et al., 2015) clipOverlap (--storeOrig OC --stats --unmapped). Overlap clipping removes the redundant portion of read pairs whose inserts are shorter than twice the read length, preventing double-counting of bases in overlapping regions that would inflate apparent depth and bias genotype likelihood estimation in downstream GL-based analyses (ANGSD, PCAngsd, NGSadmix). Memory per sorting thread was dynamically computed as 80% of the SLURM allocation divided by the thread count, clamped between 256 MB and 1,500 MB per thread.

Each final BAM was validated with samtools quickcheck, and QC metrics were generated with samtools flagstat and samtools stats. The mean PCR duplicate rate across the cohort was [TODO: X]% (range [TODO: min]–[TODO: max]%). [TODO: N] of 226 samples had duplicate rates above [TODO: threshold], but none were excluded as all remained within acceptable ranges for low-coverage population genomics.

## 1.5 Insert size characterisation

Insert size distributions were characterised using two complementary approaches. First, approximate p95 and p99 absolute TLEN percentiles were computed from properly paired, primary, non-duplicate, non-supplementary reads with MAPQ ≥ 60 (samtools view -f 0x2 -F 0xF0C -q 60). The p99 value of 514 bp and p95 value of 446 bp were derived from [TODO: a single representative sample / the median across all 226 samples / a pooled computation — clarify, and if single, state which sample and why representative]. These values were used as the upper TLEN bound for downstream BAM filtering.

Second, detailed insert size statistics were computed for each sample using reads filtered with MAPQ ≥ 30 (-f 2 -F 3844 -q 30), restricted to same-chromosome mate pairs (RNEXT = "="). Raw statistics (mean, standard deviation, median, MAD) were computed from all valid absolute TLEN values, followed by robust trimming using a median ± 10 × MAD window. The standard Illumina PE layout was confirmed as inward-facing (FR orientation) with a read length of [TODO: 150 or 151] bp. The cohort-wide trimmed mean insert size was [TODO: X] bp (SD [TODO: X] bp, range [TODO: min]–[TODO: max] bp).

## 1.6 Population genomics BAM filtering

Final BAMs for all downstream population genomics analyses (ANGSD, PCAngsd, NGSadmix, ngsRelate) were produced by applying the following filters to the merged, duplicate-marked, overlap-clipped BAMs: mapping quality ≥ 60 (samtools -q 60); proper pairs only (-f 0x2); exclusion of unmapped, mate-unmapped, secondary, QC-failed, duplicate, and supplementary reads (-F 0xF0C, corresponding to SAM flags 0x4 + 0x8 + 0x100 + 0x200 + 0x400 + 0x800); same-chromosome mate pairs only (enforced via awk filter on the SAM stream); and absolute TLEN within the empirical p99 band (0–514 bp, enforced via awk).

A stringent MAPQ threshold of 60 was chosen because the haplotype-resolved reference genome retains both parental subgenomes as separate chromosome sequences, creating extensive homeologous regions where reads of parental origin could map ambiguously between the Gar and Mac subgenomes. At MAPQ 60, only reads with essentially unique placement are retained, ensuring that population-level allele frequency estimates are not confounded by cross-subgenome mismapping. This threshold retained approximately [TODO: X]% of mapped reads genome-wide; by comparison, a MAPQ ≥ 30 threshold would have retained [TODO: Y]%, with the additional reads concentrated in repetitive and homeologous regions where mapping confidence is inherently lower. Filtered BAMs were indexed and validated with samtools flagstat.

## 1.7 Depth and coverage quality control

Per-sample and per-chromosome depth and breadth-of-coverage statistics were computed using mosdepth [TODO: version] (Pedersen & Quinlan, 2018) across five genomic region classes: whole-chromosome, repeat (soft-masked), non-repeat, callable (normal ACGT), and non-callable. Region class BEDs were derived from the reference FASTA mask annotation and their complements computed using bedtools (Quinlan & Hall, 2010). All five region BEDs were combined into a single annotated BED for a single mosdepth run per sample, using the --thresholds 1,5,10 option to compute breadth of coverage at 1×, 5×, and 10× depth. Per-chromosome MAPQ distributions were computed by counting primary mapped reads (excluding secondary and supplementary, -F 0x904) with MQ = 0 and MQ ≥ 30.

The mean sequencing depth across the 226 samples was [TODO: X]× (median [TODO: X]×, range [TODO: min]–[TODO: max]×, computed as bases mapped by CIGAR divided by the genome size of 963,905,721 bp). At the 1× threshold, breadth of coverage was [TODO: X]% of the genome; at 5×, [TODO: X]%; at 10×, [TODO: X]%. Heuristic ANGSD per-individual depth cutoffs were derived from the depth distribution: setMinDepthInd = 3 (based on p10 depth [TODO: X]× ≥ 6×) and setMaxDepthInd = 57 (approximately 5× the 90th-percentile per-sample depth of [TODO: X]×). Complete per-sample depth statistics are reported in Supplementary Table S[TODO: N].

## 1.8 BAM provenance and final sample manifest

A provenance table was constructed linking each of the 226 samples to its raw minimap2 per-run BAM path and its final population-genomics-filtered BAM path. All 226 samples produced valid filtered BAMs and entered downstream analyses; no samples were lost during the MODULE 1 pipeline. The provenance table and filtered BAM list (bamlist.pp.samechr.tlenP99.filtered.txt, 226 entries) served as the definitive sample manifests consumed by all downstream modules (MODULE 2 population structure, MODULE 3 heterozygosity, and MODULE 4 variant calling). Exact parameters and tool versions for every processing step were recorded in per-script .arg sidecar files, which are available in the data repository.

## Key parameters summary

**Table 1a.** MODULE 1 processing parameters.

| Step | Parameter | Value |
|------|-----------|-------|
| Trimming (S01) | fastp -q | 20 |
| Trimming (S01) | fastp -f/-F/-t/-T | 5 bp each end |
| Trimming (S01) | fastp -n (max N) | 0 |
| Trimming (S01) | fastp -l (min length) | 30 bp |
| Trimming (S01) | poly-G/poly-X trim | yes |
| Species (S04) | Mash k / sketch / min-copies | 31 / 50,000 / 2 |
| Species (S04) | Distance margin | 0.002 |
| Meryl (S05) | k / log₂ threshold | 21 / 1.0 |
| Mapping (S09) | minimap2 preset + flags | -ax sr -Y -c --secondary-seq --sam-hit-only |
| Merge/dup (S10) | fixmate / markdup | -m (mate scores) / mark only, not removed |
| Merge/dup (S10) | clipOverlap | --storeOrig OC --stats --unmapped |
| Filter (S13) | MAPQ | ≥ 60 |
| Filter (S13) | SAM flags require/exclude | 0x2 / 0xF0C |
| Filter (S13) | same-chr mate / TLEN | enforced / 0–514 bp |
| Depth QC (S14) | mosdepth thresholds | 1×, 5×, 10× |
| Depth QC (S14) | region classes | chr_all, repeat, nonrepeat, callable, noncallable |

**Table 1b.** MODULE 1 observed metrics (derived from data*).

| Metric | Value |
|--------|-------|
| Samples | 226 |
| Genome size | 963,905,721 bp (28 chromosomes) |
| Run-lane units | [TODO: N total] (range [TODO: min–max] per sample) |
| Raw data | [TODO: X] Gbp |
| Post-trim Q30 | [TODO: X]% (mean) |
| Mapping rate | [TODO: X]% (mean) |
| Duplicate rate | [TODO: X]% (mean) |
| Mean depth | [TODO: X]× (range [TODO: min–max]×) |
| TLEN p95 / p99* | 446 / 514 bp |
| Insert size (trimmed mean) | [TODO: X] bp |
| MAPQ ≥ 60 retention | [TODO: X]% of mapped reads |
| Breadth ≥ 1× / ≥ 5× / ≥ 10× | [TODO: X / Y / Z]% |
| Samples entering Module 2 | 226 (no attrition) |

*Empirically derived thresholds used as filtering parameters.
