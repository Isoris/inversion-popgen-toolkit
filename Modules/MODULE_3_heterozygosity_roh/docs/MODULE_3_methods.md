# MODULE 3 — Heterozygosity & Runs of Homozygosity

## Methods

### Input preparation and quality control (A01)

Sample BAM paths were extracted from the MODULE_1 provenance manifest and validated for file existence, index presence, and read group consistency. A QC-pass BAM list and corresponding sample list were generated. BEAGLE genotype likelihood files from MODULE_2A were reformatted for ngsF-HMM input, producing a positions file and sample index file alongside the genotype likelihoods.

### Per-sample heterozygosity estimation (A02)

Individual-level heterozygosity was estimated from the population-genomics-filtered BAMs using ANGSD and realSFS. For each sample, site allele frequency (SAF) likelihoods were computed using ANGSD (`-doSaf 1`, `-GL 1`) with the *C. gariepinus* haplotype reference as ancestral sequence, restricted to callable regions. Filters included minimum base quality 20, minimum mapping quality 30, and excessive mismatch downweight (`-C 50`). The per-sample folded 1D site frequency spectrum was estimated with realSFS (`-fold 1`, maximum 2,000 iterations, convergence tolerance 1 × 10⁻¹⁶).

Per-site theta estimates (Watterson's θ_W and pairwise θ_π) were computed with `thetaStat do_stat` in non-overlapping 500 kb windows (`-type 2`, physical coordinates anchored to chromosome start). Three additional window scales were computed for local diversity landscape characterization: 5 kb windows with 1 kb step, 10 kb windows with 2 kb step, and 50 kb windows with 10 kb step.

### Runs of homozygosity detection (A03)

Runs of homozygosity were called using ngsF-HMM, which operates directly on genotype likelihoods to avoid hard-call genotyping bias at low sequencing coverage. The model was run with 10 independent random restarts (seeds 42–51), and the replicate with the highest log-likelihood was retained as the final ROH call set. ngsF-HMM classifies each site along the genome as belonging to an inbred (homozygous) or non-inbred state via a hidden Markov model, outputting per-site posterior probabilities of IBD (identity by descent) status.

### ROH parsing and FROH computation (A04)

The ngsF-HMM `.ibd` posterior output was converted to BED-format ROH intervals using a posterior probability threshold. ROH were classified into length bins reflecting different demographic origins: short ROH (< 1 Mb, reflecting ancestral population bottlenecks), medium ROH (1–5 Mb, reflecting historical inbreeding), and long ROH (> 5 Mb, reflecting recent inbreeding or consanguinity). The genomic inbreeding coefficient FROH was computed as the proportion of the callable genome covered by ROH, both genome-wide and per-chromosome. Heterozygosity was separately computed within and outside ROH intervals to characterize the contrast in diversity between inbred and non-inbred genomic segments.

### Visualization and statistical analysis (B01)

A comprehensive figure suite was generated including: genome-wide heterozygosity distributions across samples, ROH count and length distributions, FROH distributions, per-chromosome ROH heatmaps, local theta diversity ideograms, scatterplots of heterozygosity versus ROH burden, depth, and FROH, and ancestry-grouped overlay versions of all core plots using the MODULE_2B ancestry labels. Statistical summaries at both sample and chromosome levels were computed, including means, medians, ranges, and correlations between heterozygosity, ROH, and ancestry assignment. An automated report generator produced manuscript-ready Methods and Results sections with all parameter values and summary statistics filled in from the pipeline outputs.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| ANGSD -GL | 1 (SAMtools) | A02 |
| ANGSD -minQ | 20 | A02 |
| ANGSD -minMapQ | 30 | A02 |
| ANGSD -C | 50 | A02 |
| realSFS -fold | 1 (folded) | A02 |
| realSFS maxIter | 2000 | A02 |
| realSFS tolerance | 1e-16 | A02 |
| thetaStat main window | 500 kb non-overlapping | A02 |
| thetaStat multiscale | 5kb/1kb, 10kb/2kb, 50kb/10kb | A02 |
| ngsF-HMM replicates | 10 | A03 |
| ngsF-HMM seed range | 42–51 | A03 |
| ROH short bin | < 1 Mb | A04 |
| ROH medium bin | 1–5 Mb | A04 |
| ROH long bin | > 5 Mb | A04 |
| Callable genome | 963,905,721 bp | A02, A04 |
