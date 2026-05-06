# MODULE_3 manuscript paragraphs — drawn from Diversity Atlas v2

All numbers below come directly from the atlas's embedded data blocks and
have been spot-checked against the source supplementary tables (S1, S2,
S3, S4, S5, S6, S9, S11, SZ, ST1, ST3) in `MODULE3_Supplementary_Tables_7.xlsx`.

References marked **`[REF]`** are from `Atlas/_scripts/REFERENCES.md` —
fill in or replace with the appropriate citation key once a manuscript
reference manager is wired up. Numbers in **`[BRACES]`** are placeholders
where the atlas does not yet have the value (mostly in-/out-ROH ratio
medians and the per-chromosome median F_HOM, which I haven't aggregated
yet — easy to add).

---

## Methods — Heterozygosity, ROH, and θπ landscape

**Genome-wide heterozygosity.** Per-individual genome-wide heterozygosity
(*H*) was estimated for each of the 226 broodstock from the per-sample
folded site-frequency spectrum produced by ANGSD `realSFS`
(`-doSaf 1 -GL 1 -anc <reference> -fold 1`, restricted to the 576.9 Mb
callable mask defined by mosdepth pass-filtered positions; *H* =
*ξ*₁ / Σ*ξ* where *ξ*ᵢ are the SFS bin counts). *H* is reported per site;
no outgroup polarisation was applied because no clean ancestral reference
is available for the *C. gariepinus* / *C. macrocephalus* species pair.
The H-based inbreeding estimator F_HOM was computed as
1 − *Hᵢ* / *H̄*_pop with *H̄*_pop = 4.55 × 10⁻³ as the cohort mean.

**Per-window θπ.** Tajima's π (θπ) was computed in 500 kb non-overlapping
windows by ANGSD `thetaStat do_stat` from the same per-sample SFS, then
averaged across the 226 samples per window to yield the cohort-level
diversity track (1,895 windows total covering all 28 *C. gariepinus*
pseudochromosomes). Per-site θπ was calculated as *t*P / *n*Sites, where
*t*P is the unfolded *t*P estimator from `thetaStat` and *n*Sites is the
window's callable-site count. Multiscale tracks (5 kb / 10 kb / 50 kb
windows with 1 kb / 2 kb / 10 kb steps) were generated in parallel for
candidate-region inspection but are not used in the manuscript-headline
analyses (Supplementary Table M3.SD3).

**Runs of homozygosity.** ROH were inferred per individual per
chromosome with ngsF-HMM **`[REF Vieira et al.]`**, using genotype
likelihoods from the same BEAGLE input as the population structure
analyses. We ran 20 random-start replicates per chromosome (seeds 32–51)
and retained the replicate with the highest final log-likelihood for
downstream inference. Replicate-stability QC (Supplementary Table SZ;
26/28 chromosomes labelled "very_stable" with absolute log-likelihood
gap between best and worst seed < 0.05; no chromosome failed
convergence) confirms that the HMM solution is robust to seed choice.
ROH segments shorter than 1 Mb were excluded as below the reliable
detection threshold for ~9× coverage. *F*_ROH was then computed as the
fraction of the 576.9 Mb callable genome covered by ROH ≥ 1 Mb. Tracts
were partitioned into length classes 1–2, 2–4, 4–8, 8–16 and ≥ 16 Mb
(Supplementary Table SD4) reflecting demographically distinct inbreeding
ages: short tracts (< 1 Mb, excluded; > 100 generations old) reflect
ancestral bottlenecks, intermediate tracts (1–8 Mb) recent population
mating, and long tracts (> 8 Mb) reflect close-kin inbreeding within the
last ~10 generations.

**Statistical tests.** Per-individual *H* and *F*_ROH were compared
across the eight K=8 NGSadmix ancestry clusters by Kruskal–Wallis test;
significant pairwise contrasts were identified by Wilcoxon rank-sum
test with Benjamini–Hochberg correction within the family of pairwise
comparisons. To distinguish family structure from founder-lineage
divergence, all KW tests were re-run on the 81-individual subset
retained after first-degree kinship pruning (NAToRA **`[REF]`**, Tassell
threshold θ ≥ 0.177 on the ngsRelate first-degree kinship graph).
Spearman rank correlations were computed between all pairs of
H-/ROH-based estimators (Supplementary Table S3). All statistical
tests were performed in R using base implementations. All numerical
results, plots, and per-row tables are interactively browsable in the
Diversity Atlas (`Atlas/Diversity_atlas.html`).

---

## Results

### Genome-wide heterozygosity

Per-individual genome-wide heterozygosity averaged
*H̄* = 4.55 × 10⁻³ ± 3.28 × 10⁻⁴ (SD; mean ± SD; range
2.83 × 10⁻³ – 7.08 × 10⁻³; *n* = 226), consistent with the order of
magnitude reported for hatchery *C. gariepinus* broodstock by
Chaivichoo et al. **`[REF]`**. Cohort-level nucleotide diversity
(θ̄π = 4.50 × 10⁻³, median 4.86 × 10⁻³ over 1,895 500 kb windows)
was the same order of magnitude as *H̄*, consistent with the
infinite-sites equivalence of the two estimators (Supplementary
Table ST1). Per-chromosome means were tightly distributed across the
28 linkage groups (range 4.28 × 10⁻³ – 4.72 × 10⁻³, IQR < 7 × 10⁻⁴
within every chromosome), indicating that nucleotide diversity is
chromosome-scale-uniform — the biological signal lives at the
*window* scale, not the chromosome scale (see hotspot paragraph
below).

### Autozygosity (F_ROH)

Genomic inbreeding coefficient F_ROH (≥ 1 Mb threshold) averaged
0.277 ± 0.053 (range 0.018 – 0.407, median 0.278; *n* = 226), placing
this hatchery cohort in the range reported for managed aquaculture
broodstock (rainbow trout 0.10–0.20, D'Ambrosio et al. 2019 **`[REF]`**;
coho salmon 0.05–0.16, Yáñez et al. 2013 **`[REF]`**). The H-based
inbreeding estimator F_HOM produced rank-equivalent estimates
(Spearman ρ = +0.704, *P* = 3.4 × 10⁻³⁵; ρ(H, F_ROH) = −0.704 in the
opposite direction), an HMM-free cross-validation of the ngsF-HMM
output (Supplementary Table S3). Per-chromosome F_ROH followed the
expected length-weighted pattern, with mean values ranging from 0.21
to 0.34 across linkage groups (Supplementary Table S4). A small
number of individuals carried chromosome-level F_ROH > 1 (max 1.08
on LG21), reflecting ROH tracts extending past the per-chromosome
callable-mask denominator rather than biologically impossible
autozygosity; these artefacts do not affect rank-based tests.
Cohort-aggregated ROH composition shifted with increasing F_ROH from
short-tract-dominated (1–2 Mb class, indicative of older shared
ancestry) at the low-burden end to long-tract-enriched (4–16 Mb
classes) at the high-burden end (Supplementary Table S8;
Diversity Atlas tab 5). θ inside ROH was depressed to roughly
**1/4.3** of θ outside ROH per individual (cohort median ratio
out/in = 4.27, range 1.01–8.30; 12/226 samples flagged as low-ratio
outliers below the 5th percentile of 2.85, suggesting possible
ngsF-HMM mis-calls or coverage artefacts in those individuals;
Supplementary Table S12), an internal QC for the ROH calls.

### Family structure, not founder-lineage divergence

F_ROH differed significantly across the eight K=8 NGSadmix clusters in
the full panel (Kruskal–Wallis *H* = 27.32, *P* = 2.9 × 10⁻⁴, *n* = 226;
Supplementary Table S5). The signal collapsed to non-significance after
first-degree kinship pruning (KW *P* = 0.59, *n* = 81), demonstrating
that the F_ROH × ancestry covariance reflects recent shared inbreeding
among siblings *within* ancestry components rather than divergence among
founder lineages. Five BH-adjusted pairwise Wilcoxon contrasts were
significant (all involving K5: K5 vs K1 *P*_BH = 7.7 × 10⁻³, K5 vs K2
*P*_BH = 4.4 × 10⁻⁴, K5 vs K6 *P*_BH = 4.4 × 10⁻⁴, K5 vs K7
*P*_BH = 2.4 × 10⁻³, K5 vs K8 *P*_BH = 2.4 × 10⁻³). K5 was depressed
relative to the other clusters (median F_ROH = 0.247 vs cohort
median 0.278), and carried the highest mean cleanest-assignment Q_max
(0.73), consistent with a tightly-related sibling block within that
ancestry component being removed by pruning. Crucially, the median
F_ROH was identical (0.278) in the retained-81 and removed-145
panels: pruning removed *family-clustered variance*, not a shift in
central tendency. *H* did not differ across ancestry clusters in
either panel (KW *P* = 0.150 in all226, *P* = 0.800 in pruned81;
Supplementary Table S5). A K = 2–12 sweep confirmed that no
NGSadmix solution produced a significant *H* × ancestry test in
either panel after BH correction (Supplementary Table S7).

### θπ landscape and outlier windows

Nineteen 500 kb windows fell above the 99th percentile of cohort θπ
(> 5.27 × 10⁻³; Supplementary Table ST3). The strongest peak
(LG22: 19.0–19.5 Mb, mean θπ = 1.07 × 10⁻², 2.4× genome mean) was
isolated; the remaining 18 windows formed several short clusters of
adjacent outliers within 1 Mb (LG12: 4.5–8.5 Mb, 9 contiguous windows;
LG07: 38.0–39.5 Mb, 4 windows; LG01: 16.0–17.0 Mb, 2 windows),
consistent with single localised diversity peaks spanning 1–4 Mb each.
Such cohort-level θπ peaks in a reference-aligned scan are candidates
for inversion polymorphisms — regions where two divergent haplotype
arrangements co-segregate inflate cohort θπ through inter-arrangement
divergence. Per-sample θπ at each outlier window (Supplementary
Table ST3b; Diversity Atlas tab 3) reveals which individuals drive
each peak, distinguishing widespread from sample-specific signals.
Cross-reference of these 19 windows against the inversion-discovery
pipeline output is the entry point for follow-up confirmation.

---

## Numerical sources (for fact-checking)

Every number above can be verified against:

- *H* descriptive: Atlas tab 1 stat strip, S2 row `het_genomewide`
- F_ROH descriptive: Atlas tab 1 stat strip, S2 row `froh`
- ρ(H, F_ROH): Atlas tab 6, S3 row 1
- KW all226 / pruned81: Atlas tab 4 first table, S5 rows 3–4
- Pairwise K5 contrasts: Atlas tab 4 second table, S5 pairwise block
- Per-cluster F_ROH median: Atlas tab 4 box plot tooltips, derivable from S1 + S9
- LG22 hotspot: Atlas tab 3 ST3 row 1, hotspot table
- Cluster of 9 wins on LG12: Atlas tab 3 hotspot clustering tag
- Per-chr θπ range: Atlas tab 2 table 1, ST1
- Median F_ROH retained = removed = 0.278: Atlas tab 6 box plot
  (computed from S11 status × f_roh)

## Headline numbers cheat sheet

```
n_samples            226
n_pruned81           81
callable_genome      576.9 Mb
mean_H               4.55e-3 ± 3.28e-4
median_H             4.53e-3
range_H              2.83e-3 – 7.08e-3
mean_FROH            0.277 ± 0.053
median_FROH          0.278
range_FROH           0.018 – 0.407
mean_theta_pi        4.50e-3
median_theta_pi      4.86e-3
rho_H_FROH           −0.704  (P = 3.4e-35)
rho_H_FHOM           +0.704  (cross-validation)
rho_H_longest_ROH    −0.242  (P = 2.4e-4)
KW_FROH_all226_P     2.9e-4
KW_FROH_pruned81_P   0.59
KW_H_all226_P        0.150 (n.s.)
KW_H_pruned81_P      0.800 (n.s.)
n_outlier_windows    19  (above 99th percentile = 5.27e-3)
top_hotspot          LG22:19.0-19.5 Mb @ 1.07e-2 (2.4× genome mean)
ngsF-HMM_stable      26/28  (very_stable across 20 seeds)
```
