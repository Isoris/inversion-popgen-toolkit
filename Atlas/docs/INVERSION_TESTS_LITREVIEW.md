# Statistical tests for chromosomal inversions: AGE and INSIDE
## Literature-grounded table of tests, tiered by importance for the *C. gariepinus* manuscript

**Scope:** 226-sample hatchery cohort, 9× Illumina, ANGSD/genotype-likelihood pipeline, 81 NAToRA-pruned unrelated individuals available for population-genetic statistics, theta-ladder + VESM + SIFT4G already wired, target Nature Communications.

**Structure:**

1. Domain 1 — Tests that estimate the **AGE** of an inversion
2. Domain 2 — Tests that characterise what is **INSIDE** an inversion (cargo, diversity, selection, recombination, load, repeats, function)
3. Master summary table (all tests, tiered)
4. Minimum methods set for a Nat Comms-tier paper in 2024–2026
5. Gaps in the current Q5/MODULE_3/MODULE_CONSERVATION pipeline
6. Recommended workflow order
7. Controversies + how to navigate them
8. The "interesting cargo" question — what reviewers actually want

Tier definitions used throughout:
- **TIER 1** — community-standard, must report. Reviewers expect to see it.
- **TIER 2** — strongly recommended, adds rigor or independent triangulation.
- **TIER 3** — niche/specific contexts; report only if the data demand it.
- **TIER 4** — not feasible at 9× / 226 / no phasing / no outgroup. Skip.

---

## Domain 1 — AGE of the inversion

### 1.1 dXY between arrangements / 2μ — the workhorse split-time estimator

- **What it estimates:** time since the most recent common ancestor of the two arrangement-defined haplotype pools, assuming complete recombination suppression and a constant mutation rate. t̂ ≈ dXY / 2μ.
- **Canonical citations:** Hudson, Kreitman & Aguadé 1987 *Genetics*; Charlesworth & Charlesworth 2010 textbook; Hey 1991 *Mol. Biol. Evol.*; Charlesworth 2024 *Genetics* on neutral diversity in inversions.
- **Inversion applications:** Berg et al. 2017 *Heredity* (Atlantic cod, ~1.5 Mya for chr 1 and 2 inversions); Mérot et al. 2020 *Mol. Biol. Evol.* (Coelopa Cf-Inv(1) — combined dXY + coalescent ~0.2–1 Mya); Han et al. 2020 / Pettersson et al. 2024 (Atlantic herring inversions); Lamichhaney et al. 2016 (ruff Faeder ~3.8 Mya); Sun et al. 2018 (sparrow ZAL2m ~2 Mya); Jay et al. 2018 (Heliconius numata P1 ~2 Mya).
- **Assumes:** single origin; no gene flux, or gene flux only at the centre of the inversion that can be excluded; mutation rate known. **Gene flux mid-inversion homogenises the centre toward genome-wide diversity within ~10⁵ generations** (Korunes & Noor 2019 *Mol. Ecol.*, gene conversion 1×10⁻⁵–2.5×10⁻⁵ converted sites/bp/generation), so naive whole-inversion dXY underestimates age for old inversions.
- **Best practice:** restrict dXY to **near-breakpoint windows** (e.g. inner 100–500 kb at each end) where gene flux is rarest. Schaeffer 2008 *Evolution*; Kapun et al. 2023 bioRxiv on *D. melanogaster* In(3R)Payne; Berdan et al. 2023 *J. Evol. Biol.* explicitly recommend this. Plot a per-window dXY landscape across the inversion and report both whole-inversion and breakpoint-restricted estimates.
- **At 9× / ANGSD:** YES — ANGSD's `realSFS fst index` produces per-pair-of-population SFS from genotype likelihoods, dXY can be computed from per-site allele frequencies. Korunes-Noor's gene-conversion detection is harder without phasing but the *point estimate* of dXY does not need phased haplotypes.
- **At n=226:** 81 unrelated × 2 chromosomes × ~3 karyotype classes gives plenty of arrangement-pool sequences; with karyotype calls of reasonable confidence, easily sufficient. Better with ≥30 individuals per arrangement; you have that for the common karyotypes.
- **Phasing required:** No (you treat the karyotype calls as an arrangement label and pool reads).
- **Failure mode:** if karyotype calls are imperfect (mis-classified heterokaryotypes), dXY is contaminated by within-arrangement diversity → systematic underestimation of age.
- **Combine with:** Tajima's D inside per arrangement (1.4); MSMC2 cross-coalescence (1.2); Twisst-style topology weighting (1.6).
- **TIER 1.** Already in your `q5_dxy_*` registry. Make sure breakpoint-restricted dXY is wired explicitly.

### 1.2 MSMC2 / PSMC' / SMC++ on per-arrangement haplotype pairs

- **What it estimates:** Ne(t) trajectory and cross-coalescence rate for two haplotype groups; the time at which cross-coalescence rate drops to zero approximates the split.
- **Canonical citations:** Schiffels & Durbin 2014 *Nat. Genet.* (MSMC); Wang et al. 2020 (MSMC2 protocol — Schiffels & Wang 2020 *Methods Mol. Biol.*); Terhorst, Kamm & Song 2017 *Nat. Genet.* (SMC++); Li & Durbin 2011 *Nature* (PSMC).
- **Inversion applications:** less common because applying MSMC2 to an inversion arrangement (rather than to whole genomes) is awkward — you need diploid pseudo-genomes per arrangement. Done in: Tuttle/Sun white-throated sparrow work; Faria et al. 2019 Littorina; the Atlantic herring 2024 paper used per-arrangement coalescent inferences. Mostly applied genome-wide rather than per-arrangement.
- **Assumes:** SMC' approximation; recombination rate known (or estimated); mutation rate known. **Important:** standard PSMC strongly biases the most-recent-time bin (Patton et al. 2024 *Curr. Biol.* "Avoidable false PSMC peaks" — first time window must be split). MSMC2 on pairs is an exact SMC' model and avoids the multi-haplotype bias of MSMC.
- **At 9× / ANGSD:** PARTIAL. PSMC/MSMC2 strictly require diploid genotype calls. There is a **ngsPSMC** for genotype-likelihood input (Lyngsø et al.) and a beta-PSMC. For inversion-specific use, you would need to pseudo-haploidise per arrangement, which is borderline at 9×. Realistically: TIER 3 in your data.
- **At n=226:** OK in principle; you'd select a few high-coverage individuals per karyotype.
- **Phasing required:** PSMC no (single diploid); MSMC2 yes for ≥2 haplotypes from one individual; cross-coalescence between two unphased diploids works.
- **Failure mode:** time-bin artefacts; sensitivity to demographic structure; recombination map needed.
- **Combine with:** dXY (1.1) for cross-validation; Tajima's D (1.4).
- **TIER 3** for your dataset. Worth doing on one or two flagship inversions where you also have higher-coverage samples; not a default for all candidates.

### 1.3 Coalescent simulation + ABC (fastsimcoal2 / msprime / dadi)

- **What it estimates:** posterior on inversion age plus other demographic parameters by simulating SFS / linkage / divergence statistics under candidate models and matching to observed data.
- **Canonical citations:** Excoffier et al. 2013 *PLoS Genet.* (fastsimcoal2); Kelleher et al. 2016 *PLoS Comput. Biol.* (msprime); Gutenkunst et al. 2009 *PLoS Genet.* (∂a∂i); Beaumont et al. 2002 *Genetics* (ABC).
- **Inversion applications:** Lamichhaney et al. 2016 *Nat. Genet.* (ruff dating); Sun et al. 2018 *Nat. Genet.* (sparrow ZAL2m); Jay et al. 2018 *Curr. Biol.* (Heliconius P locus); Pettersson 2024 (herring); Joron et al. 2011 *Nature*. Standard for headline inversions in flagship papers.
- **Assumes:** the model is correctly specified; summary statistics are informative; mutation/recombination rates are known. Sensitive to demographic mis-specification.
- **At 9× / ANGSD:** YES — you compute the SFS via ANGSD's realSFS, then use it as the input to fastsimcoal2 (which fits to SFS-based summary statistics). dadi also takes folded SFS. msprime for forward simulation does not need the data; ABC compares simulated to observed summaries.
- **At n=226:** ideal — large sample improves SFS estimation.
- **Phasing required:** No (works on SFS).
- **Failure mode:** model mis-specification gives confident-but-wrong posteriors; difficult to evaluate without sensitivity analysis. The 95% CIs are typically much narrower than the true uncertainty.
- **Combine with:** dXY (1.1), Tajima's D (1.4), Twisst (1.6) — use these for model fit checks.
- **TIER 1 for flagship inversions, TIER 2 in general.** Worth running on the 3–6 inversions you're highlighting in the manuscript. Not necessary for every minor candidate.

### 1.4 Tajima's D inside the inversion (per arrangement)

- **What it estimates:** deviation of mean pairwise diversity from Watterson's θ. Positive D → excess of intermediate-frequency variants → balancing selection or population contraction. Negative D → excess rare variants → recent sweep or expansion.
- **Canonical citation:** Tajima 1989 *Genetics*.
- **Inversion-specific interpretation:** **for the OLDER arrangement under balancing selection, expect strongly positive Tajima's D inside the inversion; for the YOUNGER (derived) arrangement post-sweep, expect negative D**. Comparing D between arrangements is a coarse age polariser. Charlesworth 2024 *Genetics* gives the modern theoretical treatment of expected diversity statistics inside inversions.
- **Inversion applications:** Kapun et al. 2023 (In(3R)Payne); Berg et al. 2016; Mérot et al. 2020; the user's own theta-ladder pipeline already produces this.
- **At 9× / ANGSD:** YES — `thetaStat do_stat` from realSFS-derived per-site thetas. This is the user's MODULE_3 output.
- **At n=226:** great with 81 unrelated.
- **Phasing required:** No.
- **Failure mode:** demography (bottlenecks, expansion) confounds; **Tajima's D is not directly an age estimator**, only a polariser of relative age combined with selection regime.
- **Combine with:** dXY (1.1) for absolute time; Fu & Li's D, Fay & Wu's H for SFS-based triangulation.
- **TIER 1.** Already wired (`q5_tajima_D_std`, `q5_tajima_D_inv`). Make sure you report **per arrangement separately**, not pooled.

### 1.5 Fu & Li's D, Fay & Wu's H, Zeng's E — additional SFS neutrality statistics

- **What they estimate:** different shapes of the site-frequency spectrum.
  - Fu & Li's D / F (Fu & Li 1993 *Genetics*) is sensitive to singletons (recent variants).
  - Fay & Wu's H (Fay & Wu 2000 *Genetics*) is sensitive to high-frequency derived alleles → signature of recent positive sweep. **Requires polarisation by an outgroup**.
  - Zeng et al. 2006 E *Genetics* is sensitive to high-frequency derived alleles after sweep recovery.
- **Inversion application:** Faria et al. 2019 Littorina; Kapun et al. 2023; classical *Drosophila* In(3R)Payne work.
- **At 9× / ANGSD:** YES for D/F (folded SFS works). Fay & Wu's H needs ancestral-state polarisation, which requires an outgroup BAM at the same sites — feasible if you have a *C. macrocephalus* genome aligned. ANGSD supports this via `-anc` flag.
- **At n=226:** fine.
- **Phasing required:** No.
- **Failure mode:** demography confounds; Fay & Wu's H very sensitive to polarisation errors.
- **TIER 2.** Strong addition: report all of (πw vs π, Tajima's D, Fu & Li's D, Fay & Wu's H) in a single per-arrangement table. Reviewers like seeing the full SFS-based battery.

### 1.6 Twisst — topology weighting

- **What it estimates:** the proportion of windowed gene trees that match each possible taxon topology when the "taxa" are arrangement labels (HOMO_1, HOMO_2, optional outgroup). If the two arrangements are clearly monophyletic (each one a clean clade) across the inversion, that is the prerequisite condition for dXY-based dating to be valid.
- **Canonical citation:** Martin & Van Belleghem 2017 *Genetics*.
- **Inversion applications:** Heliconius wing-pattern supergene work; Kim et al. 2022 *Phil. Trans. R. Soc. B* (butterfly supergene); routine in modern inversion papers (Pettersson 2024 herring).
- **Assumes:** windowed gene-tree inference is reliable; sufficient SNPs per window.
- **At 9× / ANGSD:** PARTIAL — you'd build neighbour-joining or RAxML windowed trees from genotype calls or pseudo-haploid sequences. Doable but more work than the SFS-based statistics.
- **At n=226:** fine, possibly subsample.
- **Phasing required:** Strictly no, but works better with phased haplotypes (otherwise you collapse to "consensus" per individual which loses heterozygote information). For karyotype-based topology, you can pseudo-phase by arrangement.
- **Combine with:** dXY (1.1), age estimates from coalescent simulation (1.3).
- **TIER 2.** A high-impact figure for the manuscript: a **Twisst-weighted topology landscape across the inversion** showing arrangements monophyletic in the centre and not at the breakpoints (where double-crossovers create discordant topologies) is becoming standard. **Add this to your pipeline.**

### 1.7 Phylogenetic / Dollo dating with outgroups

- **What it estimates:** if the inversion is shared (polymorphic or fixed) in two or more species, the divergence time of those species sets an upper bound on inversion age (and Dollo parsimony assumes a single origin).
- **Canonical citations:** Farris 1977 (Dollo); modern implementations in PAUP*, MrBayes; for inversions specifically: Bartolomé & Charlesworth 2006 *Genetics*; Wallace et al. 2011 *Mol. Biol. Evol.* (D. pseudoobscura phylogeny of arrangements).
- **Inversion applications:** D. pseudoobscura third-chromosome arrangements (Wallace 2011 Bayesian phylogenetics on breakpoint sequences); estrildid finch four-inversion case 2024 *MBE* (Stryjewski-Sorenson group, msa092); Atlantic cod trans-Atlantic inversion sharing.
- **At 9× / ANGSD:** YES if you have outgroup sequences for breakpoint regions (which you do via *C. macrocephalus*).
- **At n=226:** fine.
- **Phasing required:** No.
- **Failure mode:** if the inversion has polyphyletic origin (homoplasy), Dollo gives wrong ancestral state.
- **TIER 1** because you have a sister species (*C. macrocephalus*) and a third planned outgroup. Already in your `q5_dollo_*` block.

### 1.8 Molecular-clock dating of breakpoint TE insertions

- **What it estimates:** if TE insertions are exclusively on one arrangement at the breakpoints, their pairwise divergence (LTR-LTR distance for retrotransposons; or family consensus distance for DNA TEs) dates the insertion event, which sets an upper bound on inversion age (or coincides with it if TE-mediated NAHR caused the inversion).
- **Canonical citations:** SanMiguel et al. 1998 *Nat. Genet.* (LTR-LTR dating); Lerat et al. 2003.
- **Inversion applications:** Atlantic cod 2025 *GBE* paper — MITE/hAT TEs at breakpoints exclusively on inverted haplotype, used to date the inversions. Cáceres et al. 1999 D. buzzatii Galileo TE.
- **At 9× / ANGSD:** PARTIAL — needs assembly-quality sequence at breakpoints to identify TE family and compute pairwise divergence. Possible with the F1 hybrid Cgar assembly you already have, less so for the resequencing samples.
- **TIER 2** if the F1 assembly resolves breakpoints and you can polarise insertion to one arrangement; otherwise TIER 3.

### 1.9 LD-decay-based dating (recent times only)

- **What it estimates:** time since admixture / divergence from the rate at which LD decays with physical distance. Calibrated for admixture pulses.
- **Canonical citations:** Kong et al. 2010 *Nature* (relate to recombination rate), Loh et al. 2013 *Genetics* (ALDER); Patterson et al. 2012 (ROLLOFF). Recent: Hapne-LD (Fournier et al. 2023 *Nat. Commun.*); IBDNe (Browning & Browning 2015).
- **Inversion application:** rare. The LD inside an inversion does not decay (it's suppressed). The relevant signal is in flanking regions just outside breakpoints — but interpretation is messy.
- **TIER 3.** Not a primary tool here; useful only if you separately want to date a recent admixture or population split.

### 1.10 Charlesworth balancing-selection age formula

- **What it estimates:** for an inversion at intermediate frequency under balancing selection, the **expected steady-state divergence between arrangements** is approximately π_B = 2 μ Ne (this is the result for a perfectly balanced biallelic system), and the time to reach that steady state is on the order of Ne generations. Under balancing selection, dXY between arrangements does NOT keep growing with time — it equilibrates at this value. Therefore, observing dXY > 2μNe is evidence the inversion is older than equilibrium ; observing dXY ~ 2μNe means it's been at balanced equilibrium for ≥Ne generations.
- **Canonical citation:** Charlesworth 1974 *Heredity*; Strobeck 1983 *Genetics*; recent treatment Charlesworth 2024 *Genetics*.
- **Inversion applications:** Kapun et al. 2023 (In(3R)Payne); Mérot et al. 2020 Coelopa.
- **TIER 2.** Useful when you already have dXY and an Ne estimate; reports as "this dXY is consistent with steady-state balancing selection at Ne ~ X" rather than a direct age. Pair with 1.1 and 1.4.

### 1.11 ARG-based methods (Relate, ARG-Weaver, tsinfer/tsdate)

- **What it estimates:** ancestral recombination graph for the whole region; gives per-coalescence-event ages. The Speidel/Myers Relate and Wong et al. tsdate methods give branch-length-calibrated ARGs.
- **Canonical citations:** Rasmussen et al. 2014 *PLoS Genet.* (ARG-Weaver); Speidel et al. 2019 *Nat. Genet.* (Relate); Wohns et al. 2022 *Science* (tsdate); Kelleher et al. 2019 *Nat. Genet.* (tsinfer).
- **Inversion applications:** very recent. Hejase et al. 2025 *PLoS Genet.* "Length of haplotype blocks signals SV in reconstructed genealogies" — explicitly uses Relate to detect inversion-related clade structure (chr 17 and chr 10 inversions detected from genealogy). State-of-the-art.
- **At 9× / ANGSD:** PARTIAL — Relate requires phased haplotypes; tsinfer can work on genotypes but performs better with phasing. At 9× phasing is uncertain unless reference-based imputation is available.
- **At n=226:** great.
- **TIER 3** for now in your data. **Worth piloting on one flagship inversion** as a methodological showcase if you have any phased haplotypes (e.g. from the F1 assembly samples).

---

## Domain 2 — INSIDE the inversion

### 2A. Cargo / gene content

#### 2A.1 Gene density inside vs outside

- **Test:** count of protein-coding genes per Mb inside the inversion vs (i) chromosome-wide median, (ii) length-and-GC-matched random intervals (regioneR-style, Gel et al. 2016 *Bioinformatics*).
- **Permutation framework:** **regioneR** in Bioconductor draws random matched intervals; recompute gene count per draw; observed-vs-null p-value. Alternative: bedtools shuffle with constraints.
- **TIER 1** — universal first descriptor. Trivial to implement and reviewers expect to see it.

#### 2A.2 GO-term enrichment (BP / MF / CC) of cargo

- **Tools:** topGO (Alexa et al. 2006 *Bioinformatics*) — uses elim/parent-child to handle GO hierarchy; goseq (Young et al. 2010 *Genome Biol.*) — corrects for gene-length bias; clusterProfiler (Yu et al. 2012 *OMICS*); g:Profiler (Raudvere et al. 2019 *NAR*).
- **Critical caveat:** **inversion regions have non-random gene density and gene composition** (e.g. higher TE density, often gene-poor near breakpoints). The standard GO-enrichment null is "all genes equally likely to be in the inversion", which is wrong. **Use the regioneR-style null** — sample random genome regions matched for size + chromosome + recombination-class, do GO enrichment on each, build empirical null distribution. Otherwise GO enrichment p-values are massively inflated.
- **Inversion applications:** ubiquitous (cod *GBE* 2026; sunflower Todesco 2020; Atlantic herring 2024; sardine; salmon).
- **TIER 1**, but **only if matched-null permutation is used**. Naive Fisher-exact GO is TIER 3 because reviewers will reject it.

#### 2A.3 KEGG / Reactome pathway enrichment

- Same caveats as 2A.2. Often more interpretable than GO for "is this an immune-related inversion?" type claims.
- **TIER 2.** Useful supplement but never a primary claim.

#### 2A.4 Manual candidate-gene curation

- The "is there a known sex-determination / growth / immunity gene in the inversion?" check, with the gene list interpreted in the literature. Reviewers respond more positively to this than to GO enrichment because of the inversion-cargo interpretive controversy (see §8).
- **Inversion applications:** the cod LG12 vision/development genes (Matschiner 2022; the "darter inversion shares vision genes with cod and trout" Evolution 2023 paper). Pearse 2019 trout Omy05 GREB1. Ruff Faeder/Satellite-defining genes HSD17B2 / SDR42E1 (Loveland et al. 2021 *G3*).
- **TIER 1.** Absolutely required for the manuscript narrative.

#### 2A.5 Enrichment for breeding-relevant categories (immunity, growth, reproduction, sex det)

- **Specific to your aquaculture context.** Curate a list of catfish / fish breeding-relevant genes (e.g. growth: GH, IGF1, MSTN; immunity: MHC, TLRs, IL-family; reproduction: cyp19, dmrt, sox9, foxl2, amh; sex det: amhy, sdY in salmon-type systems).
- Test by Fisher's exact for overlap with cargo, against background of all annotated genes in the genome.
- **TIER 1 for the breeding/manuscript story.** This is what makes the paper applied-relevant.

#### 2A.6 PPI-network density within cargo (STRING)

- **Test:** STRING analysis of cargo gene set; density of edges inside vs random matched gene set.
- **TIER 3.** Looks good but adds little. Skip unless you have a specific functional-module hypothesis.

#### 2A.7 Gene-family expansion/contraction (CAFE5)

- **Test:** for inversions enriched in lineage-specific expansions of a gene family (e.g. immune-receptor genes); requires multiple genome assemblies.
- **Inversion application:** Dahn et al. 2026 *Genome Biol.* on Arctic codfishes — CAFE5 on cargo of inversions vs chromosome background.
- **TIER 3** — useful only if you have several catfish genome assemblies (which you partly do).

### 2B. Within-inversion sequence diversity (per arrangement)

#### 2B.1 π (nucleotide diversity), Watterson's θ, π/θ ratio per arrangement

- **What:** standard SFS-based summaries; the **expectation under balancing selection inside an inversion is reduced π in the minor (younger) arrangement and chromosome-typical π in the major (older) arrangement in the long run, but elevated π in the centre when the two arrangements coalesce.**
- **Canonical citation:** Charlesworth 2024 *Genetics*; Berdan et al. 2023 *J. Evol. Biol.*; Wellenreuther & Bernatchez 2018 *TREE* table.
- **At 9× / ANGSD:** YES — your theta-ladder pipeline already computes this. Make sure it is reported **per arrangement separately** (HOMO_1 only, HOMO_2 only, plus optionally heterokaryotype "between" diversity).
- **TIER 1.**

#### 2B.2 Tajima's D / Fu & Li / Fay & Wu — see §1.4–1.5 above

- **TIER 1** for D, **TIER 2** for the others.

#### 2B.3 HKA test (polymorphism vs divergence)

- **What:** Hudson-Kreitman-Aguadé 1987 *Genetics* — a chi-square test that the ratio of polymorphism within species to divergence between species is uniform across loci. Loci with too much polymorphism relative to divergence → balancing selection; too little → directional selection.
- **Application to inversions:** test inversion-region against genome-wide background. Detects balancing selection when polymorphism is abnormally high relative to outgroup divergence.
- **At 9× / ANGSD:** PARTIAL — needs an outgroup. You have *C. macrocephalus*; doable.
- **TIER 2.** Strong evidence for balancing selection if positive. Often skipped because of outgroup-data hassle, but you have the outgroup.

#### 2B.4 LD decay inside vs outside

- **Test:** mean r² as a function of physical distance inside the inversion (computed per arrangement separately!) and chromosome-wide. Inversions show **slow LD decay across the inversion in heterokaryotypes** (no recombination) but normal LD decay within homokaryotypes. The contrast is the diagnostic.
- **At 9× / ANGSD:** YES — ngsLD (Fox et al. 2019) for genotype-likelihood-aware LD; or `plink --r2`. Phasing helps but is not strictly required.
- **TIER 1.** Universal supplementary figure.

#### 2B.5 Heterozygote-excess inside the inversion

- **What:** inside an inversion, heterokaryotypes have all the between-arrangement variation as heterozygous sites → predicted **elevated observed heterozygosity inside inversion in heterokaryotypes**, no change in homokaryotypes. π_AB ≈ π_A + π_B + 2 dXY (Charlesworth 2024).
- **TIER 1** — easy, diagnostic, often the cleanest visual evidence the inversion is real.

### 2C. Selection scans inside the inversion

**Critical preface:** Berdan et al. 2023 explicitly warn that selection scans inside inversions are **strongly confounded by recombination suppression itself** — low recombination inflates LD-based statistics, deflates SFS-based statistics, and produces false positives in many scans. Bioarxiv 2025 paper (Selection inference in a complex genomic landscape, Schiønning-style) shows that scan signals must be **interpreted in homokaryotypes only** to disentangle haplotype-specific sweeps from inversion-induced LD elevation.

#### 2C.1 iHS, nSL — within-population EHH-based scans

- **What:** Voight et al. 2006 *PLoS Biol.* (iHS); Ferrer-Admetlla et al. 2014 *Mol. Biol. Evol.* (nSL). nSL is more robust to recombination-rate variation and is preferred for inversion regions.
- **Canonical inversion application:** the European spruce bark beetle paper (bioRxiv 2025 "Selection inference in a complex genomic landscape: the impact of polymorphic inversions") explicitly recommends running iHS/nSL **on inversion-homozygote subsets only**, to detect adaptive sweeps within an arrangement.
- **At 9× / ANGSD:** PARTIAL — selscan and rehh require phased haplotypes. Klassmann & Gautier 2022 *PLoS One* extended EHH stats to unphased data (rehh's `unphased` mode), which is what you need.
- **At n=226 (~50–80 per homokaryotype):** OK but borderline; nSL/iHS usually want ≥40 phased haplotypes per pop.
- **Phasing required:** ideally yes. Doable with read-aware phasing (Beagle 5 / SHAPEIT4) at 9× if you have a reference panel; less reliable de novo at 9×.
- **TIER 2.** Run on homokaryotype subsets; report nSL preferred.

#### 2C.2 XP-EHH, XP-nSL — cross-population EHH

- **What:** Sabeti et al. 2007 *Nature* (XP-EHH); Szpiech et al. 2021 (XP-nSL). Tests for sweep in one arrangement vs the other.
- **Application to inversions:** comparing HOMO_1 vs HOMO_2 directly. **A clear sweep in one arrangement post-inversion** is the diagnostic for recent adaptive evolution within the inversion.
- **TIER 2.** Same phasing constraint as 2C.1.

#### 2C.3 SweepFinder2 / SweeD — composite likelihood sweeps

- **What:** Nielsen et al. 2005 (SweepFinder); Pavlidis et al. 2013 *Mol. Biol. Evol.* (SweeD); DeGiorgio et al. 2016 (SweepFinder2). Compares observed SFS to a sweep-shaped expected SFS.
- **At 9× / ANGSD:** PARTIAL — needs called or imputed genotypes. Doable.
- **Phasing required:** No.
- **TIER 2** if you can run on homokaryotype subsets.
- Note Crisci et al. 2013 *Front. Genet.* "impact of equilibrium assumptions" — bottlenecks easily produce false-positive SweepFinder/SweeD signals. Always run with appropriate demographic null.

#### 2C.4 OmegaPlus — LD-based sweeps

- **What:** Alachiotis et al. 2012 *Bioinformatics*. Detects the characteristic LD pattern around a sweep.
- **Inside inversions:** **strongly biased** because LD is already elevated. Use only on homokaryotype subsets and even then with caution.
- **TIER 3** — niche.

#### 2C.5 RAiSD — multi-signature scan

- **What:** Alachiotis & Pavlidis 2018 *Commun. Biol.*. Combines SFS, LD, and diversity drop. Fast, scalable.
- **Phasing:** No.
- **TIER 2.** Modern replacement for OmegaPlus + SweepFinder; nice single-tool approach.

#### 2C.6 PCAdapt — PC-based outlier detection

- **What:** Luu et al. 2017 *Mol. Ecol. Resour.*. PCA-based outliers, identifies SNPs strongly contributing to top PCs.
- **Inside inversions:** **karyotype IS one of the top PCs**, so naive PCAdapt will flag the entire inversion as an outlier — uninformative. Useful only **inside homokaryotype subsets**, where it can detect within-arrangement local adaptation.
- **At 9× / ANGSD:** YES — PCAngsd implements PCAdapt-style scans on genotype likelihoods (Meisner et al. 2021 *Bioinformatics*; PMC8480091).
- **TIER 2** — restricted to homokaryotype subsets.

#### 2C.7 BayPass / RDA / LFMM — genotype-environment association inside the inversion

- **What:** test whether allele frequency inside the inversion correlates with environmental variables. BayPass (Gautier 2015 *Genetics*); RDA (Forester et al. 2018 *Mol. Ecol.*); LFMM (Frichot et al. 2013 *Mol. Biol. Evol.*).
- **Application to inversions:** clean signal that the inversion is locally adapted is a clinal frequency change with environment. Without environmental variables (single hatchery cohort), **not applicable in your dataset**.
- **TIER 4 in current dataset.** TIER 1 in the future *C. macrocephalus* wild-cohort paper.

#### 2C.8 hapFLK / FLK — population-tree-aware test

- **What:** Bonhomme et al. 2010; Fariello et al. 2013 *Genetics*. F-statistic that accounts for population structure via a hierarchical population tree.
- **Inside inversions:** rarely applied; useful when you have multiple populations.
- **TIER 4 in current single-population dataset; TIER 2 in multi-population work.**

### 2D. Recombination, gene flux, within-inversion architecture

#### 2D.1 LD heatmap (R²) inside vs outside

- See 2B.4. The visual confirmation of recombination suppression. Do this **separately for HOMO_1, HOMO_2, and HET** subsets — only HET shows the suppression.
- **TIER 1.**

#### 2D.2 Population recombination rate ρ = 4Nec inferred (LDhat / pyrho / FastEPRR / ReLERNN)

- **What:** infer ρ along the chromosome; see the inversion-induced trough.
- **Tools:** LDhat (McVean et al. 2004); pyrho (Spence & Song 2019 *Sci. Adv.*) — handles unphased data + variable Ne; FastEPRR (Gao et al. 2016) — coalescent likelihood; ReLERNN (Adrion et al. 2020 *Mol. Biol. Evol.*) — deep learning, robust to inversions; LDhelmet (Chan et al. 2012).
- **Inversion applications:** ReLERNN paper specifically demonstrates inversion-detection from recombination-rate troughs; pyrho protocol on cichlid inversions.
- **At 9× / ANGSD:** pyrho is the right tool — handles unphased data and demography. ReLERNN trains a neural network on simulated data; doable.
- **TIER 2.** Adds a per-window ρ landscape figure that visually anchors the recombination-suppression claim. Adds methodological credibility.

#### 2D.3 Gene-conversion detection (Betran/Schaeffer-Anderson; GENECONV)

- **What:** identify SNPs that diagnostically polarise to the "wrong" arrangement → evidence of gene conversion track. Tract length 200–500 bp in *Drosophila* (Korunes & Noor 2019; Schaeffer & Anderson 2005).
- **At 9× / ANGSD:** PARTIAL — needs phased / arrangement-pooled haplotypes; reasonable with karyotype calls.
- **TIER 3** for routine reporting; **TIER 2 if you want to make a "selected genes within inversion" claim** (Korunes et al. 2024 *G3* — gene conversion reveals selected genes inside inversions). The Korunes 2024 method is the right modern reference.
- Worth running on flagship inversions only.

#### 2D.4 Four-gamete test

- **What:** Hudson & Kaplan 1985 *Genetics*; classic test for absence of recombination between two SNPs (presence of all 4 gametes is evidence for recombination). Inside an inversion in heterokaryotypes, the four-gamete test should fail almost everywhere except where gene conversion occurred.
- **TIER 3** — qualitative, mostly superseded by quantitative ρ inference.

#### 2D.5 Distance-from-breakpoint π profile (Schaeffer's "U-shape")

- **What:** plot π between arrangements as a function of distance from the breakpoint. Older inversions show a **U-shape**: high dXY near breakpoints (no gene flux) decaying toward genome-wide values in the centre (gene flux homogenises). Recent inversions show **flat dXY** profile.
- **Inversion applications:** Schaeffer 2008 *Evolution*; Kapun et al. 2023 *In(3R)Payne* — used as both age and gene-flux diagnostic.
- **TIER 1.** Easy from per-window dXY (you already have this in MODULE_3). One panel per inversion. Diagnostic at a glance.

### 2E. Deleterious load / mutation accumulation

#### 2E.1 Counts of derived deleterious variants per arrangement (SIFT, PolyPhen, VESM, GERP, CADD, PROVEAN)

- **What:** count putatively deleterious variants per arrangement, compare. Key: **per-haplotype** counts (not per-individual), to detect arrangement-specific accumulation. Equivalent to "load score per arrangement".
- **Best practice:** use multiple predictors (SIFT4G, VESM, CADD, GERP) and report concordance. Your MODULE_CONSERVATION (SIFT4G + VESM + bcftools csq + Cactus alignment) is the right setup.
- **Inversion applications:** Heliconius numata (Jay et al. 2021 *Nat. Genet.*) — heavy load on the derived inversion arrangements, frequency-dependent selection; Coelopa Cf-Inv(1) (Berdan et al. 2021 *Evol. Lett.*); Sunflower (Huang et al. 2022 *Mol. Biol. Evol.*) — load negatively correlated with inversion heterozygosity.
- **TIER 1** — central narrative test for any inversion paper that touches on supergene maintenance or breeding programmes.

#### 2E.2 R_xy statistic (Do et al. 2015 *Nat. Genet.*)

- **What:** a ratio of derived-allele frequencies of putatively deleterious vs putatively neutral variants, between two populations. R_xy > 1 means population X has higher deleterious load than Y. Adapted for arrangements: R_xy(HOMO_1 vs HOMO_2) tests for differential load.
- **Inversion application:** rare formal use, but conceptually applied in many supergene-load papers (e.g. Heliconius, Coelopa).
- **At 9× / ANGSD:** YES — needs derived-allele-frequency estimates per arrangement, which you can get by partitioning samples by karyotype and estimating per-pool AF.
- **TIER 2.** Strong: explicit ratio test, paired with confidence intervals via jackknife or bootstrap. **Add this — it's missing from your registry.**

#### 2E.3 dN/dS inside the inversion (PAML CodeML; HyPhy aBSREL/BUSTED/RELAX)

- **What:** ratio of nonsynonymous to synonymous substitutions in cargo genes between arrangements (or to outgroup). Long-term selective regime.
- **Canonical:** Yang 2007 *Mol. Biol. Evol.* (PAML); Smith et al. 2015 *Mol. Biol. Evol.* (HyPhy aBSREL); Wertheim et al. 2015 *Mol. Biol. Evol.* (RELAX).
- **At 9× / ANGSD:** needs codon-aware alignment, which your bcftools csq + Cactus pipeline can produce. Doable.
- **Caveat:** dN/dS over short timescales (within-species) is uninformative; meaningful only with outgroup divergence. *C. macrocephalus* gives that.
- **TIER 2.** Useful for supplementary table; not always the primary load test.

#### 2E.4 McDonald-Kreitman test

- **What:** McDonald & Kreitman 1991 *Nature*. Contingency table of polymorphism and divergence × synonymous and nonsynonymous. Detects recurrent positive selection (DN/DS > PN/PS).
- **Caveats:** very sensitive to demographic non-equilibrium; weakly deleterious mutations bias the test (use iMKT, Murga-Moreno et al. 2019 *NAR*). Power is low per gene; pool genes (e.g. all genes inside the inversion) for a single test.
- **At 9× / ANGSD:** YES with outgroup. iMKT framework gives α (proportion of adaptive substitutions).
- **TIER 2.** Good supplementary evidence for adaptive cargo evolution. Note Wen et al. 2020 bioRxiv on MKT vs PAML overlap and false-positive rates — be cautious.

#### 2E.5 F_ROH inside the inversion (homozygosity contrast)

- **What:** are HOMO_1 individuals more homozygous (ROH-rich) inside the inversion than expected? Indicates recessive-load purging or recent bottleneck specific to the arrangement.
- **At 9× / ANGSD:** YES via your existing F_ROH module restricted to inversion coordinates.
- **TIER 2.**

#### 2E.6 Heterozygote-deficit / homozygote-excess of putatively deleterious variants

- **What:** for putatively recessive deleterious variants inside the inversion, test for HWE deviation. If they're under purifying selection in homozygotes, the population shows heterozygote excess at those sites specifically.
- **TIER 3** — niche.

#### 2E.7 Number of segregating LoF (loss-of-function) variants per arrangement

- Similar to 2E.1 but restricted to LoF (stop-gain, frameshift, splice-disrupting) — the "high-confidence deleterious" subset. Use bcftools csq + LOFTEE-style filters.
- **TIER 2.** Often more interpretable than full SIFT/CADD because LoF is unambiguous.

### 2F. Repeats and structural variants inside

#### 2F.1 TE family counts inside vs outside (RepeatMasker / EarlGrey)

- **Test:** EarlGrey (Baril et al. 2024 *Mol. Biol. Evol.*) is the modern de novo + curated TE annotation pipeline; RepeatMasker for known libraries.
- **What:** test whether inversion is enriched for specific TE families relative to chromosome background, using regioneR-style permutation. The Cáceres et al. 1999 *Drosophila* / Atlantic cod 2025 *GBE* hypothesis is that breakpoints are enriched, not the interior — separately analyse interior and breakpoint flanks.
- **TIER 1.** Already wired in your TE-density pipeline.

#### 2F.2 Tandem-repeat density (TRF, GMATA)

- **Test:** Benson 1999 *NAR* (TRF). Microsatellite/minisatellite enrichment is a separate signal from interspersed-TE.
- **TIER 2.**

#### 2F.3 Structural-variant load inside (DELLY / Manta / GRIDSS)

- **What:** count of indels/CNVs/inversions-within-inversions per arrangement. Inversions with reduced recombination accumulate other SVs.
- **TIER 2.**

#### 2F.4 LTR-LTR pairwise divergence dating of TE insertions (see 1.8)

- **TIER 2** for breakpoint-flanking TEs.

### 2G. Functional-genomics signals (when phenotype/expression data exist)

#### 2G.1 Allele-specific expression in heterokaryotypes

- **Test:** compare expression of inversion-side alleles vs standard-side alleles in heterokaryotypes (RNA-seq from tissues). The Sun et al. 2018 sparrow paper showed dosage compensation; Loveland et al. 2021 G3 ruff paper showed allelic imbalance. Berdan et al. 2021 *Evol. Lett.* Coelopa.
- **Requires:** RNA-seq + heterokaryotype individuals.
- **TIER 4 in current dataset (no RNA-seq).** TIER 1 if you ever generate RNA-seq for the cohort.

#### 2G.2 eQTL enrichment

- **Test:** are inversion-tag SNPs enriched among genome-wide eQTLs? Requires expression data.
- **TIER 4 in current dataset.**

#### 2G.3 Hi-C / TAD-boundary disruption

- **TIER 4 in current dataset** (no Hi-C). TIER 2 if you generate Hi-C as part of the genome-assembly paper that's separate.

#### 2G.4 Methylation / chromatin

- **TIER 4 in current dataset.**

### 2H. Cohort/family-aware tests (specific to your aquaculture context)

#### 2H.1 Kinship-aware karyotype frequency (ngsRelate-pruned)

- Re-run all per-arrangement statistics on NAToRA-pruned 81 unrelated individuals, then on full 226. If discrepancies are large, the inversion signal is partly a family-pack artefact.
- **TIER 1.** This is your core defence against the "is this a real inversion or a family-clustering artefact?" reviewer attack.

#### 2H.2 Mixed-model test of karyotype-trait association

- **What:** if phenotype data are available, fit `phenotype ~ karyotype + (1|family)` using lme4. Random family effect absorbs pedigree confound.
- **TIER 1 if phenotypes exist.**

#### 2H.3 Cochran-Mantel-Haenszel test of karyotype-family independence

- **What:** stratified test for karyotype clumping in families. Tests whether some families have skewed karyotype frequencies (which would indicate parental segregation rather than population-level polymorphism).
- **TIER 2.** A Q6 family-linkage diagnostic.

---

## 3. Master summary table (all tests)

| Domain | Test | Key citation | ANGSD-OK | n=226-OK | Phasing | Tier | One-line interpretation if positive |
|---|---|---|---|---|---|---|---|
| Age | dXY/2μ between arrangements | Hudson 1987; Charlesworth 2024 | ✅ | ✅ | ❌ | 1 | Inversion split time in My/generations |
| Age | Breakpoint-restricted dXY | Schaeffer 2008; Berdan 2023 | ✅ | ✅ | ❌ | 1 | Less gene-flux-biased age |
| Age | Tajima's D per arrangement | Tajima 1989 | ✅ | ✅ | ❌ | 1 | Older balancing-selected arrangement → +D |
| Age | Fu&Li D, Fay&Wu H, Zeng E | Fu&Li 1993; Fay&Wu 2000 | ✅ | ✅ | ❌ | 2 | Triangulate SFS shape |
| Age | Charlesworth π_B equilibrium check | Charlesworth 1974, 2024 | ✅ | ✅ | ❌ | 2 | "dXY at 2μNe" → balanced equilibrium |
| Age | Coalescent + ABC (fastsimcoal2) | Excoffier 2013 | ✅ | ✅ | ❌ | 1 (flagship) / 2 | Posterior on age + demography |
| Age | MSMC2 / SMC++ cross-coalescence | Schiffels 2014 | ⚠️ | ✅ | partial | 3 | Time of cross-coalescence drop |
| Age | Twisst topology weighting | Martin & Van Belleghem 2017 | ✅ | ✅ | partial | 2 | Arrangement monophyly across inversion |
| Age | Phylogenetic / Dollo with outgroup | Bartolomé & Charlesworth 2006 | ✅ | ✅ | ❌ | 1 | Polarisation + upper bound on age |
| Age | TE molecular-clock at breakpoints | SanMiguel 1998; Atlantic cod 2025 | partial | ✅ | partial | 2 (if assembly resolves BP) |
| Age | LD-decay-based admixture/split | Loh 2013 ALDER | ✅ | ✅ | ❌ | 3 | Recent timescales only |
| Age | ARG-based (Relate, tsdate) | Speidel 2019; Wohns 2022 | ⚠️ | ✅ | yes | 3 | Per-coalescence-event ages |
| Cargo | Gene density (regioneR null) | Gel 2016 | ✅ | NA | ❌ | 1 | Inversion is gene-rich/poor |
| Cargo | GO enrichment (matched-null) | Alexa 2006 topGO | ✅ | NA | ❌ | 1 | Functional enrichment categories |
| Cargo | KEGG / Reactome | Yu 2012 | ✅ | NA | ❌ | 2 | Pathway enrichment |
| Cargo | Manual candidate-gene curation | — | ✅ | NA | ❌ | 1 | Specific functional gene present |
| Cargo | Breeding-relevant gene set test | — | ✅ | NA | ❌ | 1 | Aquaculture relevance |
| Cargo | PPI-network density (STRING) | — | ✅ | NA | ❌ | 3 | Functional module |
| Cargo | CAFE5 family expansion | Mendes 2020 | ⚠️ | NA | ❌ | 3 | Lineage-specific expansion in cargo |
| Diversity | π / θw per arrangement | Nei 1987 | ✅ | ✅ | ❌ | 1 | Diversity contrast |
| Diversity | π_AB heterokaryotype excess | Charlesworth 2024 | ✅ | ✅ | ❌ | 1 | Confirms inversion is real |
| Diversity | LD heatmap inside vs outside | — | ✅ | ✅ | partial | 1 | Recombination suppression |
| Diversity | HKA polymorphism vs divergence | Hudson 1987 | ✅ (with outgroup) | ✅ | ❌ | 2 | Balancing or directional selection |
| Selection | iHS / nSL on homokaryotypes | Voight 2006; Ferrer-Admetlla 2014 | partial | borderline | yes | 2 | Within-arrangement sweep |
| Selection | XP-EHH / XP-nSL HOMO1 vs HOMO2 | Sabeti 2007 | partial | borderline | yes | 2 | Differential sweep between arrangements |
| Selection | SweepFinder2 / SweeD | Pavlidis 2013 | ✅ | ✅ | ❌ | 2 | Composite-likelihood sweep |
| Selection | RAiSD multi-signature | Alachiotis & Pavlidis 2018 | ✅ | ✅ | ❌ | 2 | Combined sweep signature |
| Selection | OmegaPlus | Alachiotis 2012 | ✅ | ✅ | ❌ | 3 | LD-based sweep (biased inside inv) |
| Selection | PCAdapt on homokaryotypes | Luu 2017 (PCAngsd) | ✅ | ✅ | ❌ | 2 | Within-arrangement local adaptation |
| Selection | BayPass / RDA / LFMM (GEA) | Gautier 2015 | ✅ | ✅ | ❌ | 4 (current) / 1 (multi-pop future) | Environmental association |
| Selection | hapFLK / FLK | Fariello 2013 | ✅ | ✅ | partial | 4 (current) | Multi-population test |
| Recombination | ρ along chromosome (pyrho) | Spence & Song 2019 | ✅ | ✅ | ❌ | 2 | Recombination-rate trough |
| Recombination | ReLERNN deep learning | Adrion 2020 | ✅ | ✅ | ❌ | 3 | Modern alternative to pyrho |
| Recombination | Gene-conversion detection | Schaeffer & Anderson 2005; Korunes 2019, 2024 | partial | ✅ | partial | 2 (flagship) | Selected genes within inversion |
| Recombination | Distance-from-breakpoint π profile | Schaeffer 2008; Kapun 2023 | ✅ | ✅ | ❌ | 1 | U-shape = old + gene flux |
| Load | Per-arrangement deleterious counts (SIFT/VESM/CADD) | Multiple | ✅ | ✅ | ❌ | 1 | Differential load |
| Load | R_xy statistic | Do 2015 | ✅ | ✅ | ❌ | 2 | Formal load ratio with CI |
| Load | dN/dS (PAML CodeML, HyPhy) | Yang 2007 | ✅ (with outgroup) | ✅ | ❌ | 2 | Long-term selection regime |
| Load | McDonald-Kreitman / iMKT | McDonald & Kreitman 1991; Murga-Moreno 2019 | ✅ (with outgroup) | ✅ | ❌ | 2 | Recurrent positive selection |
| Load | F_ROH inside inversion | — | ✅ | ✅ | ❌ | 2 | Recessive-load purging |
| Load | LoF variant counts per arrangement | bcftools csq + LOFTEE | ✅ | ✅ | ❌ | 2 | High-confidence load asymmetry |
| Repeats | TE family counts (RepeatMasker / EarlGrey) | Smit; Baril 2024 | ✅ | NA | ❌ | 1 | TE accumulation |
| Repeats | Tandem repeats (TRF / GMATA) | Benson 1999 | ✅ | NA | ❌ | 2 | Microsatellite enrichment |
| Repeats | SV load (DELLY / Manta) | — | ✅ | ✅ | ❌ | 2 | Other SVs accumulate inside |
| Functional | Allele-specific expression | Sun 2018 sparrow | needs RNA-seq | NA | partial | 4 (current) | Inversion changes expression |
| Functional | eQTL enrichment | — | needs eQTL | NA | ❌ | 4 (current) | Functional cargo |
| Functional | Hi-C TAD disruption | — | needs Hi-C | NA | ❌ | 4 (current) | Chromatin reorganisation |
| Cohort | ngsRelate-pruned recomputation | Korneliussen 2014 | ✅ | ✅ | ❌ | 1 | Family-confound check |
| Cohort | GLMM phenotype ~ karyotype + (1|family) | lme4 | ✅ | ✅ | ❌ | 1 (if pheno) | Karyotype-trait association |
| Cohort | Cochran-Mantel-Haenszel | — | ✅ | ✅ | ❌ | 2 | Karyotype-family clumping |

---

## 4. Minimum methods set for a Nat Comms-tier inversion paper, 2024–2026

Reading the methods sections of recent inversion papers in Nat Commun, Nat Ecol Evol, Mol Biol Evol, Mol Ecol, GBE, JEB (Berdan 2023; Pettersson 2024 Atlantic herring; Harringmeyer & Hoekstra 2022; Reeve et al. 2023 Littorina; Jay 2018, 2021; Atlantic cod 2025/2026; Stryjewski-Sorenson estrildid 2024 MBE; Songbird inversion 2025 bioRxiv), the consensus minimum is:

1. **Local PCA discovery** (lostruct or equivalent) ← you have this.
2. **Three-cluster karyotype calling** with explicit cluster-quality metric ← you have this.
3. **LD heatmap inside vs outside** (per-karyotype) ← partially have.
4. **Per-arrangement π, θw, Tajima's D landscape** ← you have this.
5. **dXY between arrangements**, both whole-inversion AND breakpoint-restricted ← you have whole, **add breakpoint-restricted explicitly**.
6. **Distance-from-breakpoint dXY profile** (the U-shape figure) ← **add this**.
7. **Twisst-weighted topology landscape** showing arrangement monophyly ← **add this**.
8. **Polarisation by outgroup (Dollo)** with at least one outgroup ← you have *C. macrocephalus*; wire `q5_dollo_*` properly.
9. **Coalescent / ABC age estimate for flagship inversions** ← **add this for 3–6 highlighted inversions**.
10. **Cargo gene list + GO enrichment with matched-null permutation** ← partially have; **make the null matched, not all-genes**.
11. **Manual curation of candidate genes** ← do this.
12. **TE family enrichment at breakpoints** (NAHR mechanism evidence) ← you have this.
13. **Per-arrangement deleterious load** (SIFT + VESM + CADD concordance) ← you have this.
14. **R_xy or equivalent formal load ratio with CI** ← **add this**.
15. **Family-confounding check via ngsRelate-pruned re-run** ← you have this; make sure it's explicit in figures.
16. **Recombination-rate landscape (pyrho or LDhat)** showing the inversion trough ← **add pyrho**.

Items in **bold** are gaps in your current pipeline that should be added before submission.

---

## 5. Gaps in current Q5 / MODULE_3 / MODULE_CONSERVATION pipeline

### Critical gaps (TIER 1 missing or incomplete)

- **Breakpoint-restricted dXY** is not explicitly wired as a flat key. Your `q5_dxy_*` is presumably whole-inversion. Add `q5_dxy_breakpoint_left`, `q5_dxy_breakpoint_right`, `q5_dxy_centre`, plus `q5_dxy_profile_shape` (linear / flat / U-shape).
- **Twisst topology weighting** is missing. Add as `q5_twisst_arrangement_monophyly_weight` (the % of windows where arrangements are monophyletic) and the per-window TSV as a Layer-D attached file.
- **Coalescent / ABC age estimate** with confidence interval. Currently `q5_dollo_mya` is the upper-bound age; add `q5_age_abc_point`, `q5_age_abc_lo95`, `q5_age_abc_hi95`, `q5_age_abc_method` (fastsimcoal2 / msprime+ABC).
- **U-shape dXY classification** (`q5_dxy_profile_class` ∈ {flat, linear-decay, U-shape}). One-line summary that reviewers can read.
- **Matched-null GO enrichment**. Your current GO infrastructure is "planned"; make sure the null is **size-and-recombination-class-matched random intervals**, not all genes.
- **Per-arrangement R_xy load ratio** with jackknife CI. You have load counts; combine into a formal Do-2015-style ratio.
- **pyrho ρ landscape** showing recombination-rate trough at the inversion. Adds methodological credibility.

### Secondary gaps (TIER 2 missing)

- **Fay & Wu's H per arrangement** (needs *C. macrocephalus* polarisation, which you have).
- **HKA test** of inversion vs background (needs outgroup).
- **iHS/nSL on homokaryotype subsets** with rehh `unphased` mode.
- **Gene-conversion detection** on flagship inversions (Korunes 2024 method).
- **dN/dS PAML/HyPhy on cargo** (you have the alignments).
- **iMKT** for adaptive substitution proportion.
- **TE-LTR molecular-clock dating** at breakpoints if assembly resolves them.

### TIER 3 niche tests in your current setup that could be deprioritised

- `q5_gds_gap`, `q5_gds_fst_spearman_rho` — not standard literature. Either justify them or drop in favour of the more standard set.
- `q5_segregating_sites_std` and `q5_segregating_sites_inv` — useful for HKA but otherwise rarely the first reported quantity. Keep but de-emphasise.

### TIER 4 (not feasible / not applicable)

- MSMC2 / SMC++ — defer to a flagship inversion only.
- ARG-based (Relate, tsdate) — defer to a methodological appendix or skip.
- BayPass / RDA / LFMM — defer to the multi-population *C. macrocephalus* paper.
- Allele-specific expression / eQTL / Hi-C — only if you generate the data.

---

## 6. Recommended workflow order

The principle: **independent evidence first, derived/synthesis second**. Don't compute MK before you have the SFS; don't run pyrho before you've done LD heatmap; don't run ABC before you've computed dXY.

### Step 1 — Per-arrangement SFS (the foundation)
- Karyotype calls finalised (you have this).
- ANGSD `realSFS` per arrangement on the unrelated 81 subset.
- Per-arrangement π, θw, Tajima's D, Fu&Li D, Fay&Wu H, Zeng E.
- dXY between arrangements, both whole-inversion and breakpoint-restricted.

### Step 2 — Inversion structure descriptors
- LD heatmap per karyotype subset.
- Distance-from-breakpoint π profile (U-shape diagnostic).
- pyrho ρ landscape (recombination trough).
- Twisst arrangement monophyly weight.

### Step 3 — Polarisation and age
- Dollo with *C. macrocephalus* outgroup.
- HKA test (uses outgroup).
- TE breakpoint dating (if assembly resolves).
- Coalescent / ABC age for 3–6 flagship inversions.

### Step 4 — Selection signatures (homokaryotype subsets only!)
- iHS/nSL with rehh unphased mode.
- XP-nSL between HOMO_1 and HOMO_2.
- SweepFinder2 + RAiSD.
- PCAdapt within homokaryotypes.

### Step 5 — Cargo characterisation
- Gene list extraction.
- regioneR-style matched-null gene density and GO enrichment.
- Manual curation against breeding-relevant gene catalogue.
- TE family enrichment inside vs at breakpoints.

### Step 6 — Load
- Per-arrangement deleterious-variant counts (SIFT + VESM + CADD concordance).
- R_xy ratio with jackknife CI.
- LoF variant counts.
- F_ROH inside inversion.
- dN/dS / PAML on cargo (long-term).
- iMKT (recent positive selection).
- Gene-conversion detection on flagship inversions.

### Step 7 — Cohort validation
- ngsRelate-pruned re-run of all per-arrangement statistics.
- GLMM phenotype × karyotype with family random effect (if phenotypes exist).

### What can be parallelised vs sequential

- Steps 1, 2, 5, 6 are independent → run in parallel.
- Step 3 depends on Step 1 (needs SFS / dXY).
- Step 4 depends on Step 1 and Step 2 (needs karyotype + recombination).
- Step 7 is a sanity check across all of the above.

### Independent-evidence triangulation pairs (what NOT to report alone)

| Claim | Don't rely on | Triangulate with |
|---|---|---|
| Recombination is suppressed | LD heatmap alone | + pyrho ρ landscape + four-gamete test |
| Inversion is X My old | dXY alone | + Tajima's D pattern + ABC + Twisst monophyly + outgroup polarisation |
| Cargo is functionally enriched | GO enrichment with naive null | + matched-null permutation + manual curation + breeding-relevant set test |
| Sweep inside inversion | iHS alone | + nSL + RAiSD + SweepFinder2; on homokaryotypes only |
| Differential load | SIFT counts | + R_xy + LoF counts + dN/dS |
| Inversion is balanced-selected | HWE deviation alone | + π_B at Charlesworth equilibrium + Tajima's D positive in older arrangement + frequency stable across populations |

---

## 7. Controversies and how to navigate them

### 7.1 dXY age vs MSMC2 age can disagree by an order of magnitude

This is real and unresolved. dXY uses point divergence and assumes a clock; MSMC2 uses cross-coalescence and is sensitive to demographic structure. Berdan et al. 2023 explicitly notes that **the relevant time** is "time of split between arrangements" which is not the same as "time of inversion mutation" — the inversion may have arisen on a haplotype that already had within-population diversity, so dXY/2μ overestimates the age of the inversion mutation itself. **Recommendation:** report dXY/2μ as a "split time of arrangement haplotype pools" upper bound, ABC posterior as best estimate, and explicitly note the distinction in methods. This is how Lamichhaney 2016, Sun 2018, Jay 2018 all framed it.

### 7.2 FST inside vs outside an inversion is misleading (Cruickshank & Hahn 2014)

Settled. **Always report dXY in addition to (or instead of) FST** for between-arrangement comparison. FST inflates in low-diversity regions (which inversions often are post-sweep), giving a false impression of differentiation. Cruickshank & Hahn 2014 *Mol. Ecol.* is the canonical citation; Matthey-Doret & Whitlock 2019 *Mol. Ecol.* is the follow-up confirming background-selection inflates FST.

### 7.3 Selection scans inside inversions produce false positives because of recombination suppression alone

Settled. **Run iHS, nSL, OmegaPlus, RAiSD ONLY on homokaryotype subsets**, never on the pooled cohort. The bioRxiv 2025 European spruce bark beetle paper is explicit on this. Run within-HOMO_1 and within-HOMO_2 separately, then ask whether either subset shows local sweep signal.

### 7.4 "Inversions cause speciation" vs "inversions protect already-divergent regions"

Berdan et al. 2023 explicitly resolve this in favour of the latter for most documented cases. **Don't claim inversions are speciation drivers** unless you have explicit hybrid-fitness or RI-locus-mapping evidence; instead claim they "facilitate maintenance of divergent allele combinations under gene flow" (Kirkpatrick & Barton 2006 framing).

### 7.5 GO enrichment is weak evidence for "interesting cargo"

True. See §8 below.

### 7.6 dN/dS within species is often non-informative

McDonald & Kreitman 1991 is more powerful than dN/dS for within-species selection inference. Use dN/dS for between-species comparison only (i.e. *C. gariepinus* vs *C. macrocephalus*).

### 7.7 Tajima's D is not an age estimator

It's a polariser combined with a selection inference. Saying "the inversion is 1 My old because Tajima's D is +2.1" is wrong. Saying "Tajima's D is positive in arrangement A and approximately zero in arrangement B, consistent with arrangement A being the older / balancing-selected one" is right.

### 7.8 The "inversion frequency stability over time" claim

A claim that frequencies have been stable for thousands of generations requires temporal sampling (museum / archived specimens; Mérot 2025 Coelopa paper does this). Single-time-point cohorts cannot make this claim. Don't overclaim balance based on one snapshot.

### 7.9 Underdominance vs overdominance vs neutral models for intermediate-frequency inversions

Berdan, Blanckaert, Butlin & Bank 2021 *PLoS Genet.* show that **intermediate-frequency inversions can arise from accumulation of recessive deleterious load alone, without any positive selection**. So an intermediate-frequency inversion is NOT evidence for adaptive maintenance unless you separately demonstrate positive fitness effect or environmental association. Conservative framing is "inversion at intermediate frequency, consistent with multiple maintenance mechanisms including [list]".

---

## 8. The "interesting cargo" question — what reviewers actually want

GO enrichment alone does not satisfy reviewers in 2024–2026. Looking at recent successful inversion papers (Heliconius numata, Atlantic cod, deer mouse, ruff, sparrow, sunflower, Atlantic herring, Mimulus, sardine, bark beetle, songbird inversions), the components of a convincing "this cargo is interesting" claim are:

1. **A specific named gene with a published functional role**, ideally in a related species' breeding/biology literature. Reviewers respond to "the inversion contains *amh*, the master sex-determination locus in salmonids" much more than to "GO category 'reproductive process' is enriched, p=0.003".

2. **Convergent functional category across multiple inversions in the manuscript** — if 4 of your 6 inversions contain immune-receptor genes (TLRs, NLRs, MHC), that's a story. Single-inversion enrichments are weak.

3. **Cross-species cargo conservation** — if the same gene families recur in inversions in multiple fish species (cod LG12, trout Omy05, your catfish LG28 — the darter paper Evolution 2023 makes exactly this argument with vision genes), that's strong evidence the inversion captures functionally important variation. **This is your cross-species page 13 narrative.**

4. **Trait association** — if you have phenotypes (body weight, growth rate, survival), test karyotype × trait with family-aware GLMM. Even one significant association elevates the cargo claim from "potentially functional" to "functionally validated".

5. **Allele-specific expression / eQTL** — gold standard but requires RNA-seq.

6. **Permutation-based null** — every cargo claim must be referenced against a null of size-matched, recombination-class-matched random intervals. The Cabrera-Brufau et al. 2022 *Genome Res.* and Werme et al. 2022 *Nat. Genet.* approaches are the modern standard.

7. **Negative cargo statement** — sometimes the most convincing claim is "the inversion is gene-poor and contains no obvious candidates", which then redirects the story to neutral/load-driven maintenance (Berdan 2021 framework). **Don't force a positive cargo claim if the data don't support one.**

The structure of a strong cargo paragraph in a 2024+ paper looks like:

> "The inversion spans X.Y Mb and contains N genes, no different from chromosome-wide expectation (matched-null permutation p = 0.42). However, the cargo is enriched (matched-null p < 0.001) for genes annotated to GO category 'response to oxidative stress', driven by three named genes [G1, G2, G3] previously implicated in [trait]-related QTL in related species [refs]. Two of these genes are in the homologous inversion identified in [related species, ref], suggesting cross-species conservation of inversion content. In our cohort, karyotype-stratified GLMM finds significant association between karyotype and [trait] (β = X, 95% CI [Y, Z], p = 0.001), with random family effect explaining N% of variance. We interpret this as evidence that the inversion captures functionally relevant variation in [trait pathway], consistent with [hypothesis: local adaptation / breeding selection / etc]."

That's the template. Match it.

---

## End notes

Tier rankings here are my reading of the current literature and your data constraints. They are not absolute — for example, BayPass is TIER 4 in your single-population hatchery cohort but will become TIER 1 in the *C. macrocephalus* multi-population paper. Keep the framework axis-keyed (per the previous review's six-axis framework) and the test-tier table will rebalance automatically when cohort context changes.

The tests you are MOST likely missing and that would most strengthen the manuscript right now are, in priority order:
1. Breakpoint-restricted dXY + U-shape profile classification.
2. Twisst topology weighting per inversion.
3. ABC coalescent age estimate for 3–6 flagship inversions.
4. Matched-null GO + cargo enrichment.
5. R_xy formal load ratio with CI.
6. pyrho recombination landscape.

These six additions, together with what you already have, constitute the Nat Commons-tier minimum methods set.
