# MODULE 1 — Read Preparation and BAM Processing

## Methods

### Read quality control and trimming (S01–S03)

Raw paired-end Illumina reads (PE151) were processed with fastp (v0.23.x) using the following parameters: adapter auto-detection for paired-end reads (`--detect_adapter_for_pe`), poly-G and poly-X tail trimming (`--trim_poly_g --trim_poly_x`), 5′ and 3′ quality trimming enabled (`-5 -3`), minimum Phred quality score 20 (`-q 20`), 5 bp trimmed from both ends of both reads (`-f 5 -F 5 -t 5 -T 5`), zero N bases allowed (`-n 0`), and minimum read length 30 bp (`-l 30`). Samples with corrupted gzip outputs (unexpected EOF) were identified and re-processed from original raw reads (S01–S02). Per-sample QC metrics including Q20/Q30 percentages, adapter trimming rates, duplication rates, and insert size peaks were extracted from fastp log files (S03) and summarized across the cohort.

### Species verification (S04–S06)

Parentage species identity of each sequencing run-lane unit was verified using Mash distance-based assignment against the *C. gariepinus* (Gar) and *C. macrocephalus* (Mac) haplotype-resolved reference genomes. Mash sketches were built with k=31, sketch size=50,000, minimum k-mer copies=2. Each run-lane unit was assigned as Gariepinus, Macrocephalus, or Ambiguous using a distance margin of 0.002. Per-CGA (per-sample) consensus calls were derived by majority vote across run-lane units (S04). An independent verification using meryl k-mer species-specific sets was also available (S05), using k=21 and a log₂ ratio threshold of 1.0.

### Read alignment (S07–S09)

Trimmed reads were aligned to the *C. gariepinus* haplotype reference (fClaHyb_Gar_LG) using minimap2 in short-read mode (`-ax sr`) with supplementary alignment soft-clipping (`-Y`), secondary sequence output, and SAM hit-only mode. Read groups were set per run-lane unit with SM=sample, LB=sample, PU=run.lane, PL=ILLUMINA. Alignments were coordinate-sorted and indexed with samtools.

### BAM merging and duplicate marking (S10)

Per-run BAMs were merged into per-sample BAMs with samtools merge, preserving read group information. RG/SM tag consistency was validated before merging. The merged BAM was name-sorted, processed with samtools fixmate (-m), coordinate-sorted, and duplicate-marked (mark-only, not removed) with samtools markdup. Overlapping read pairs were clipped using BamUtil clipOverlap (`--storeOrig OC --stats --unmapped`). Each final BAM was indexed and QC'd with samtools quickcheck, flagstat, and stats.

### Insert size characterization (S11–S12)

TLEN percentiles (p95, p99) were computed from properly paired, primary, non-duplicate, non-supplementary reads with MAPQ ≥ 60 and same-chromosome mate pairs (S11). Detailed insert size statistics including raw and MAD-trimmed (median ± 10×MAD) distributions were computed with MAPQ ≥ 30 (S12).

### Population genomics BAM filtering (S13)

Final BAMs for all downstream population genomics analyses were produced by filtering the merged markdup+clip BAMs with the following criteria: MAPQ ≥ 60, proper pairs only (`-f 0x2`), exclusion of unmapped, mate-unmapped, secondary, QC-fail, duplicate, and supplementary reads (`-F 0xF0C`), same-chromosome mate pairs only, and absolute TLEN within the empirical p99 band (0–514 bp). Filtered BAMs were indexed and flagstat was computed.

### Coverage and depth QC (S14–S15)

Per-sample and per-chromosome depth and breadth-of-coverage statistics were computed using mosdepth across five region classes: whole-chromosome, repeat, non-repeat, callable, and non-callable (as defined by the assembly mask BEDs). Thresholds at 1×, 5×, and 10× were recorded. MAPQ distributions (MQ0, MQ30+) were computed per chromosome. Summary statistics (min, p10, p25, p50, mean, p75, p90, max) were generated across samples for all metrics. Heuristic ANGSD depth cutoffs were suggested based on the depth distribution (minDepthInd=3 if p10 ≥ 6×, else 2; maxDepthInd=min(60, round(p90×5))).

### BAM provenance (S16)

A final provenance table was constructed linking each sample to its raw minimap2 BAM and popgen-filtered BAM, serving as the definitive sample manifest for MODULE_2 (biSNP discovery) and MODULE_3 (ancestry and relatedness).

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| fastp -q | 20 | S01 |
| fastp -f/-F/-t/-T | 5 bp each | S01 |
| fastp -n | 0 | S01 |
| fastp -l | 30 | S01 |
| Mash k | 31 | S04 |
| Mash sketch | 50,000 | S04 |
| Mash min copies | 2 | S04 |
| Mash margin | 0.002 | S04 |
| minimap2 preset | sr | S09 |
| Filter MAPQ | ≥ 60 | S13 |
| Filter flags exclude | 0xF0C | S13 |
| Filter TLEN max | 514 (p99) | S13 |
| Genome size | 963,905,721 bp | S15 |
