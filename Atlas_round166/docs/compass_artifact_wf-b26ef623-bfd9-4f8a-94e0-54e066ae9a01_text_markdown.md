# A Manuscript-Grade Classification Framework for Chromosomal Inversions in Population Genomics: Literature Synthesis and Unified Axis Hierarchy

## TL;DR

- **The literature converges on six orthogonal axes that describe every inversion** — *Discovery / Evidence*, *Structure*, *Mechanism / Origin*, *Evolutionary History*, *Population Dynamics*, and *Function & Application* — and these map cleanly onto your existing Q1–Q7 registry, your 14-axis tier view, and your 5-bucket stats profile, collapsing redundancy and exposing two genuine gaps (a true *Mechanism* axis and an explicit *Detectability/Evidence* layer).
- **Your three current organising schemes are not actually three theories — they are three projections of the same six-axis space**, mostly differing in granularity. Q1–Q3 are mostly Structure + Mechanism; Q4 is Mechanism (cargo + breakpoint substrate); Q5 is Evolutionary History; Q6 is Population Dynamics; Q7 (existence A–D + tier) is the *Discovery/Evidence* axis. The "stats profile" buckets are a presentation re-slicing and should be retired as an organizing layer; keep them only as figure-panel groupings.
- **For the *C. gariepinus* manuscript and the portable *Macrocephalus*/aquaculture export**, adopt a hierarchy of (i) Evidence Layer A–D as a *gate*, (ii) Structure & Mechanism as the *descriptive* axes (size, position, breakpoint substrate, simple/nested), (iii) History as the *origin* axis (age, polarisation, ancestral vs derived, recurrence), (iv) Population Dynamics as the *frequency* axis (HWE, family confound, cline, substructure), (v) Function as the *consequence* axis (cargo, GO enrichment, trait association), and (vi) Application as the *breeding-utility* axis (Tier 1–4 markers). Each axis gets one omnibus test plus one to two descriptors — no method is duplicated across axes.

---

## Key Findings

### 1. The literature has crystallised, since ~2018, around a small, recurrent set of descriptive dimensions

The modern conceptual scaffolding — Hoffmann & Rieseberg 2008 *Annu. Rev. Ecol. Evol. Syst.*; Wellenreuther & Bernatchez 2018 *Trends Ecol. Evol.*; Faria, Johannesson, Butlin & Westram 2019 *Trends Ecol. Evol.* ("Evolving Inversions"); Mérot, Oomen, Tigano & Wellenreuther 2020 *Trends Ecol. Evol.* (the SV roadmap); and especially **Berdan et al. 2023 *J. Evol. Biol.*** ("How chromosomal inversions reorient the evolutionary process") — uses a remarkably consistent toolbox of axes. Berdan et al. 2023 is the most useful template because it explicitly partitions inversion biology into **(a) origin**, **(b) selection regime**, **(c) life history of the inversion through time**, **(d) molecular consequences (recombination, drift, mutation accumulation)**, and **(e) phenotypic and species-level outcomes**. Wellenreuther & Bernatchez 2018 Table 1 is the canonical comparative table with axes: species, inversion name, **size**, **number of genes**, **estimated age (Mya)**, **selective process**. Their Figure 1 is the "paracentric vs pericentric" visual that almost every review reuses.

This convergent vocabulary means a reviewer trained on these papers expects to see, for every inversion: *evidence layer*, *size*, *centromeric position*, *breakpoint mechanism*, *cargo content*, *age estimate*, *frequency / cline information*, *signature of selection*, and *trait association*. Your six-axis collapse should mirror those expectations.

### 2. Foundational frameworks and what each foregrounds (the "axis emphasis" map)

- **Sturtevant (1917, 1921) and Dobzhansky (Genetics & the Origin of Species, 1937–1951; Krimbas & Powell 1992 *Drosophila Inversion Polymorphism*, CRC):** founded the cytological vocabulary — *paracentric vs pericentric*, *gene arrangement* (e.g. ST, AR, CH on Drosophila pseudoobscura 3rd chromosome), *cosmopolitan vs endemic*, *common vs rare*. Foregrounded **Structure** (paracentric/pericentric) and **Population Dynamics** (cline-tracking, frequency change over decades).
- **White 1973 *Animal Cytology and Evolution* and Stebbins 1958:** established that pericentric inversions in heterozygotes produce unbalanced gametes (underdominance) while paracentrics produce acentric/dicentric products that are inviable but do not depress fertility because the offending gamete is selected against at the polar body — this is the canonical "pericentric matters more for hybrid sterility" dichotomy. Foregrounds **Mechanism** (recombination consequences).
- **Kirkpatrick & Barton 2006 *Genetics*:** the local-adaptation-with-migration model — an inversion can spread by capturing locally adapted alleles and shielding them from gene flow, *without epistasis*. This is the dominant theoretical reference for *why* inversions are at intermediate frequency. Foregrounds **Evolutionary History** (origin of selective benefit) and **Function** (locally adaptive cargo). Their Figure 1 (selection-migration regimes) is the canonical schematic. Note the 2018 erratum for their Equation 3 — Charlesworth & Barton 2018 *Genetics* gives the corrected n-locus rate.
- **Hoffmann & Rieseberg 2008 *Annu. Rev. Ecol. Evol. Syst.*:** synthesised the shift from "inversions = underdominance speciation engines" to "inversions = recombination-suppression vehicles for adaptive divergence". Foregrounds **Function** and **Evolutionary History**. Provides the standard list of phenotypic traits associated with inversions.
- **Wellenreuther & Bernatchez 2018 *Trends Ecol. Evol.*:** the comparative-table review — every inversion ever characterised gets a row with size, age, selective regime. Foregrounds **Structure** (size), **History** (age), **Function** (selective regime). Their Table 1 should be your template for the manuscript's per-inversion summary table.
- **Faria et al. 2019 *Trends Ecol. Evol.* "Evolving Inversions":** adds time as an explicit axis — inversions are not static; their allelic content evolves via mutation, gene conversion, and rare double crossovers, and the *age* and *life-stage* of the inversion change which evolutionary forces dominate. Foregrounds **Mechanism** (gene flux) and **Evolutionary History** (lifecycle). Introduces the "Type I vs Type II polymorphism" distinction (transient/divergent vs balanced).
- **Mérot, Oomen, Tigano & Wellenreuther 2020 *Trends Ecol. Evol.*:** roadmap for *detection* of structural variants, including local-PCA, long-read-vs-short-read, pangenome-graph approaches. Foregrounds **Discovery / Evidence**.
- **Mérot 2020 *Mol. Ecol.* (commentary on Huang et al.):** lays out the genome-scan-then-test-association methodology that is now standard. Foregrounds **Discovery** + **Function**.
- **Berdan et al. 2023 *J. Evol. Biol.*** (19-author consortium): the explicit decomposition of inversion evolution into mechanism categories that can each be **disentangled** by a specific signature. Their Figure 2 (coalescent predictions for genes across a chromosomal region under different evolutionary models) is the schematic to reproduce mentally for every inversion: neutral, locally adapted breakpoints, locally adapted alleles within, balancing selection, and so on. Foregrounds **all six axes** but most directly **Mechanism** and **History**.
- **Connallon & Olito 2022 *Mol. Ecol.*:** length distribution as a diagnostic — neutral and underdominant inversions trend small, locally adaptive ones trend large. Foregrounds **Structure** as a *signature* of **Evolutionary History**. Useful because it gives you a principled way to interpret your size_class.
- **Schwander, Libbrecht & Keller 2014 *Curr. Biol.*; Thompson & Jiggins 2014 *Heredity*; Charlesworth 2016 (*Heredity*; *PLoS Biol.* commentary 2016 on supergenes); Gutiérrez-Valencia, Hughes, Berdan & Slotte 2021 *Genome Biol. Evol.*:** the supergene-as-special-case-of-inversion literature. Foregrounds **Function** (multi-trait coordinated phenotype), **Mechanism** (recombination suppression by an inversion as the supergene-forming agent), and the *fates* axis (degeneration, balanced lethals, AOD) which Berdan, Blanckaert, Butlin & Bank 2021 *PLoS Genet.* models explicitly.
- **Huang & Rieseberg 2020 (*Mol. Ecol.*)** and **Todesco et al. 2020 *Nature*** (sunflower haploblocks — 37 large 1–100 Mb non-recombining haplotype blocks): foreground **Discovery** (local-PCA / haplotype-block scans on RAD or whole-genome data) and **Function** (genotype-environment association).

### 3. The discovery / evidence axis is universally underspecified in the literature but central to your pipeline

There is no single agreed evidence ladder, but a defensible four-tier scheme emerges from the methods literature:

- **A — Direct sequence evidence of inversion**: assembly-based confirmation (PacBio HiFi / ONT contigs spanning the breakpoint with split alignment in opposite orientation; pangenome graph bubble; Hi-C contact-map flip; karyotype/FISH). Sniffles, SVIM, SVIM-asm, pbsv, npInv (NAHR-specialised), Sawfish, and minigraph/PGGB/Cactus are the standard tools (see SVIM, Heller & Vingron 2019 *Bioinformatics*; PacBio benchmarks; the *Nat. Commun.* 2024 SV-detection benchmark paper).
- **B — Strong indirect signal across population**: local-PCA producing three discrete clusters (homo/hetero/homo) at expected 1:2:1 or skewed proportions (Li & Ralph 2019 *Genetics* "lostruct"; the lostruct R package). Three-cluster PCA is the **single most diagnostic indirect signature** and is what Huang et al. 2020 *Mol. Ecol.* (sunflower), Todesco et al. 2020 *Nature*, Mérot et al. 2021 *Mol. Ecol.* (Coregonus), Harringmeyer & Hoekstra 2022 *Nat. Ecol. Evol.* (deer mouse), Reeve, Butlin, Koch, Stankowski & Faria 2023 *Mol. Ecol.* (Littorina), and Le Moan et al. 2024 *Evol. Lett.* (Littorina fabalis) all use.
- **C — Linkage/LD signal alone**: long-distance high-LD blocks (Kemppainen et al. 2015 in Littorina; the Faria et al. 2019 *Mol. Ecol.* "evolutionary chunks" approach), reduced recombination on a linkage map.
- **D — Phenotypic / cline signal alone**: differentiation outlier or cline that is later interpreted as inversion-driven; this is the weakest evidence level and corresponds to historical "islands of divergence" first documented in Atlantic cod (Hemmer-Hansen et al. 2013; Berg et al. 2016 *Heredity*).

This A–D ladder is exactly your `q7_existence_layer`, and it is the right gate for the manuscript: report all candidates at layer A or B, treat C as supporting, treat D as exploratory.

### 4. Structural axes — converging numerical thresholds

- **Size**: the literature has no consensus, but a workable trichotomy from synthesising plant/animal surveys is **small <1 Mb, intermediate 1–5 Mb, large >5 Mb**, with the additional sometimes-used cut at **>10 Mb = "supergene-scale"**. Empirical anchor points: Drosophila In(3R)Payne ~8 Mb; Anopheles 2La ~22 Mb (~10% of genome); ruff ~4.5 Mb; white-throated sparrow ZAL2m ~100 Mb; Heliconius numata P1+P2+P3 nested ~1.7 Mb total; Coelopa frigida Cf-Inv(1) ~25 Mb; sunflower haploblocks 1–100 Mb (Todesco et al. 2020 *Nature*); deer mouse 1.5–43.8 Mb (Harringmeyer & Hoekstra 2022 *Nat. Ecol. Evol.*); Atlantic cod 4 inversions of several Mb each (Berg et al. 2016 *GBE*; Kirubakaran et al. 2016; Barth et al. 2017 *Heredity*); Atlantic salmon SVs (Bertolotti et al. 2020 *Nat. Commun.*); rainbow trout Omy05 ~55 Mb pericentric. Connallon & Olito 2022 *Mol. Ecol.* gives a theoretical reason for the upper-tail bias under local adaptation: large inversions are favoured when they capture migration-protected locally adapted alleles. Crop-genomics groups (Huang et al. 2024 *Plant Biotechnol. J.*) propose a relative measure (% of chromosome) but absolute Mb thresholds dominate animal literature. Your `size_class` (small <1Mb / medium 1–5Mb / large >5Mb) is essentially correct; consider adding a "supergene-scale ≥10 Mb" tier for the *Macrocephalus* portability.
- **Centromeric position**: paracentric (one arm) vs pericentric (spanning centromere). Functionally pericentric inversions in heterozygotes more often produce duplication-deficiency gametes and so contribute more strongly to underdominance; paracentrics largely escape that cost in animals (acentric/dicentric products are eliminated at meiosis I or II) but in plants single crossovers within a paracentric inversion can produce inviable gametes (Stebbins 1958 effects). White 1973 and Hoffmann & Rieseberg 2008 are the citations.
- **Position relative to telomere**: Harringmeyer & Hoekstra 2022 found that deer mouse inversion breakpoints **frequently occur in centromeric and telomeric regions**, reflecting both the recombinational landscape and the substrate for ectopic recombination at repeat-rich ends. Your sub-telomeric LG28 case fits this pattern. The relevant citations are (i) Harringmeyer & Hoekstra 2022, (ii) Pevzner & Tesler 2003 *PNAS* and Peng, Pevzner & Tesler 2006 *PLoS Comput. Biol.* for the "fragile breakage" model with reusable breakpoint regions enriched at chromosome ends, and (iii) Murphy et al. 2005 for mammalian breakpoint reuse. Sub-telomeric vs interstitial vs peri-centromeric is a defensible three-way classification.
- **Single vs nested vs complex**: Drosophila pseudoobscura has nested inversions on the third chromosome (Schaeffer 2008 *Evolution*; Wallace, Detweiler & Schaeffer 2011 *Mol. Biol. Evol.*); Anopheles 2La/2R complexes likewise; Heliconius numata has three sequentially accreted inversions P1+P2+P3 (Jay et al. 2018 *Curr. Biol.*; Jay, Chouteau et al. 2021 *Nat. Genet.*). Coelopa frigida Cf-Inv(1) is itself three overlapping inversions (Aziz 1975; Mérot et al. 2020 *Mol. Biol. Evol.*).
- **Cargo content**: gene-rich / gene-poor / contains-immune-or-MHC / contains-reproductive-genes. Corbett-Detig & Hartl 2012 *PLoS Genet.* showed that cosmopolitan Drosophila inversions are enriched for breakpoints in *intergenic* regions — i.e. selection or mutational bias against breakpoints that disrupt gene sequences. The "gene-rich vs gene-poor" axis is meaningful because gene-rich inversions provide more selective targets and more potential trait associations.

### 5. Mechanism / origin axes — the breakpoint-substrate signature

There are two canonical mechanisms producing chromosomal inversions, and they leave different molecular signatures at breakpoints:

- **Ectopic recombination / NAHR between inverted repeats** (often LTR retrotransposons, LINEs, or DNA-transposon TIRs). This produces breakpoints flanked by *chimeric* repeat copies (Cáceres et al. 1999 *PNAS* and Casals et al. 2003 on Drosophila Galileo transposon; Delprat et al. 2009 *PLoS One* in *D. buzzatii*; Lakich et al. 1993 *Nat. Genet.* — the factor VIII inversion in haemophilia A; Lupski reviews on NAHR). In Atlantic cod, a recent paper (cited by you in 2024 work — Matschiner / Tørresen / Star group, *GBE* 2025 in PMC12368964) finds **tandem accumulation of MITEs / hAT transposons** at the breakpoints of all four large cod inversions, exclusively on the inverted haplotypes — strong evidence for TE-mediated ectopic recombination. The deer mouse case (Harringmeyer & Hoekstra 2022) shows breakpoints flanked by 0.5–50 kb inverted repeats. NAHR between LINE elements has been characterised genome-wide in humans (Startek et al. 2015 *NAR*).
- **Non-homologous end-joining / FoSTeS / MMBIR**: produces breakpoints with no homology, often with short microhomologies, sometimes with flanking inverted duplications (Hastings, Lupski, Rosenberg & Ira 2009 *Nat. Rev. Genet.*).

Your Q4 already captures TE-family enrichment, GO enrichment, and tandem-repeat context at breakpoints — this is the right Mechanism axis. The Kuang et al. 2024 catfish chromosome-evolution paper you already cite fits cleanly here as the species-specific anchor. Add explicit reporting of: (a) repeat density on each side of each breakpoint, (b) whether the flanking elements are the **same family** in **opposite orientation** (the diagnostic for NAHR), (c) whether short microhomology or inverted duplications are present (the diagnostic for NHEJ/MMBIR).

### 6. Evolutionary history axes — age, polarisation, recurrence

- **Age estimation**: use **dXY / 2μ between arrangements** measured at sites near breakpoints (because gene flux via gene conversion erodes the signal in inversion centres — Schaeffer & Anderson 2005 *Genetics*; Schaeffer 2008; Korunes & Noor 2019 in *D. pseudoobscura* show conversion tract length ~205 bp, rate ~3.4×10⁻⁶, two orders of magnitude above mutation rate). **Tajima's D** *inside* the inversion is expected to be elevated under balancing selection (older arrangements) and depressed under recent sweeps. **MSMC2 / IBDseq / coalescent simulation**-based methods give absolute age but are sensitive to the demographic model. ABC has been used for ruff (Lamichhaney et al. 2016 *Nat. Genet.* — Faeder ~3.8 Mya, Satellite ~0.5 Mya), white-throated sparrow ZAL2m (~2–3 Mya per Tuttle et al. 2016 *Curr. Biol.*; Sun et al. 2018 *Nat. Genet.*), Heliconius numata P locus (~2 Mya for P1, younger for nested P2/P3 — Jay et al. 2018), Atlantic cod inversions (>100 ka, predating trans-Atlantic separation — Berg et al. 2017 *Heredity*), Coelopa Cf-Inv(1) (~0.2–1 Mya, Mérot et al. 2020 *Mol. Biol. Evol.*), Mimulus DIV1 (Lowry & Willis 2010 *PLoS Biol.*; Twyford & Friedman 2015 *Evolution*), seaweed fly clinal inversions (Mérot et al. 2018 *Proc. R. Soc. B*).
- **Polarisation (ancestral vs derived)**: requires outgroup orientation and Dollo-parsimony reconstruction with **two or more outgroups**. Critical for interpretation: if the derived arrangement is at high frequency, balancing selection or a recent sweep must be invoked; if ancestral, the inversion may simply be the standard. Your Q5 Dollo block is the right slot for this.
- **Recurrence / hotspots**: Pevzner & Tesler 2003 *PNAS* / Peng, Pevzner & Tesler 2006 fragile-breakage; Murphy et al. 2005 mammalian; Harringmeyer & Hoekstra 2022 deer mouse breakpoints in centro/telomeric regions; Ranz et al. 2007 *Drosophila*; the Atlantic cod TE-tandem accumulation paper.
- **Type I vs Type II polymorphism (Faria et al. 2019)**: Type I = transient (selection driving towards fixation in one population); Type II = balanced (polymorphism maintained by balancing selection of any flavor). This is a useful evolutionary-state classification.

### 7. Population-dynamics axes

- **Frequency class**: rare (<5%), intermediate (5–50%), common (>50%) — the operational thresholds used by Wellenreuther & Bernatchez 2018 Table 1 and most subsequent papers.
- **Cline / GEA**: latitudinal or environmental clines are the strongest evidence of selection. Examples: Drosophila In(3R)Payne / In(3L)P (Kapun & Flatt 2019 *Mol. Ecol.*; Rane et al. 2015; Kapun et al. 2023 *Mol. Biol. Evol.*); Anopheles 2La aridity cline (Coluzzi et al. 1979; Ayala & Coluzzi 2005 *Mol. Ecol.*; Gray et al. 2009 *Malar. J.*; Fouet et al. 2012 *PLoS One*; Cheng et al. 2012 *Genetics*); seaweed fly Cf-Inv(1) and Cf-Inv(4.1) cline (Mérot et al. 2018 *Proc. R. Soc. B*; Mérot et al. 2025 bioRxiv).
- **Hardy–Weinberg dynamics**: deviation from HWE under heterokaryotype advantage (e.g. Coelopa Cf-Inv(1), ruff Satellite homozygote-lethal). Your `q6_hwe_*` block is well-placed.
- **Family / pedigree confound**: this is the *aquaculture-specific* axis the literature largely ignores. In a hatchery cohort with ~226 samples and known family structure, an inversion-like local-PCA signal can be a family-pack artefact. The defensible fix is (i) ngsRelate or KING relatedness matrix, (ii) re-running the local-PCA on a relatedness-pruned subset, and (iii) testing for genotype–family association (Cochran–Mantel–Haenszel or mixed-model LRT). Cite Andrews et al. 2018 / Robledo et al. 2018 aquaculture-genomics reviews; for the principle, Cridland et al. 2013 / Huang et al. 2014 on hidden relatedness in DGRP (Mackay et al. 2012 *Nature*).
- **Substructure within karyotype**: the Schaeffer *D. pseudoobscura* tradition — sub-arrangements within an arrangement, gene conversion creating mosaic haplotypes within the inverted background. Berdan, Blanckaert, Butlin & Bank 2021 *PLoS Genet.* models exactly this fates dynamic and finds that inversions can branch into multiple highly differentiated within-arrangement haplotypes that halt fitness decay. Methodologically: silhouette + gap statistic + bootstrap hierarchical clustering on within-karyotype genotype matrices, plus DAPC-cross-validation (Jombart, Devillard & Balloux 2010 *BMC Genet.*; Thia 2023 *Mol. Ecol. Resour.* guidelines).

### 8. Function / phenotype axes

- **Trait association catalog** (canonical, organise by selection regime in the manuscript table):
  - *Mating-system/morph supergenes*: ruff (Lamichhaney 2016; Küpper 2016 *Nat. Genet.*); white-throated sparrow ZAL2m (Tuttle 2016; Sun et al. 2018; Maney et al. 2020 *PNAS*); fire ant Gp-9 (Wang et al. 2013 *Nature*; Y.-C. Huang et al. 2018 *Proc. R. Soc. B*).
  - *Mimicry/wing pattern*: Heliconius numata P locus (Joron et al. 2011 *Nature*; Jay et al. 2018 *Curr. Biol.*; Jay, Chouteau et al. 2021 *Nat. Genet.*); Papilio polytes (Nishikawa et al. 2015 *Nat. Genet.*).
  - *Migratory/life-history ecotypes*: Atlantic cod (Berg et al. 2016, 2017); Atlantic herring (Han et al. 2020); deer mouse forest/prairie (Hager et al. 2022 *Science*).
  - *Local adaptation along environmental gradient*: Drosophila In(3R)Payne; Anopheles 2La; sunflower haploblocks (Todesco 2020); Mimulus DIV1 (Lowry & Willis 2010); stickleback marine/freshwater inversions on chr I, XI, XXI (Jones et al. 2012 *Nature*).
  - *Sexual antagonism / sex chromosomes*: white-throated sparrow ZAL2m; Littorina sex-linked inversions (Hearn et al. 2022 *Mol. Ecol.*).
- **GO enrichment** at cargo content level is the standard reporting approach. Your Q4 block already encodes this.
- **Reproductive isolation contributions**: the speciation-island debate. Noor et al. 2001 *PNAS*; Rieseberg 2001 *Trends Ecol. Evol.*; Faria & Navarro 2010 *Trends Ecol. Evol.*; Mérot 2020. Empirically settled view (Berdan et al. 2023): inversions contribute to RI mostly by *protecting* divergence rather than *causing* underdominance, except in some plant systems (Stebbins).

### 9. Applied / breeding axes (highly relevant for your aquaculture context)

The aquaculture-genomics literature is starting to engage seriously with inversions:
- The 492-genome Atlantic salmon SV catalogue (Bertolotti et al. 2020 *Nat. Commun.*) is the methodological template for population-scale SV discovery in a farmed fish.
- Rainbow trout Omy05 (Pearse et al. 2014, 2019) and Omy20 (Leitwein et al. 2024 *G3*) are the canonical salmonid examples; Arctic charr LG12 (Christensen et al. 2021 *PMC8473973*) is another.
- Pacific salmon GREB1/ROCK1 spawning-time region (Prince et al. 2017 *Sci. Adv.*; Waples et al. 2022).
- Conservation-management framing: Wellenreuther et al. 2024 / Matschiner et al. 2025 *PMC12684361* "Putting Structural Variants Into Practice" (marine management).

For your Tier-1–4 marker scheme:
- **Tier 1** = breakpoint-spanning PCR amplicon with diagnostic indel/junction (gold standard; e.g. 2La PCR-RFLP — Coulibaly et al. 2016).
- **Tier 2** = within-inversion private indel that distinguishes the two arrangements with near-perfect LD.
- **Tier 3** = SNP-tag panel (5–20 SNPs in highest-LD with karyotype, validated with classification accuracy ≥98%).
- **Tier 4** = exploratory / population-scale local-PCA-based call with no validated marker.

This maps cleanly to literature on "structural-haplotype tagging marker" frameworks in dairy cattle (Daetwyler et al.), pig QTL, salmonid marker-assisted selection (Yáñez et al. — Cooke Aquaculture / SeqSNP work).

### 10. The statistical battery — what the field actually does

The Wellenreuther 2018 / Mérot 2020 / Berdan 2023 papers all use a similar minimum methods set, which I recommend as your **community-accepted minimum**:

1. **Local PCA (lostruct, Li & Ralph 2019)** in 100 SNP or 100 kb windows + MDS clustering of windows = primary inversion-discovery scan.
2. **Per-window PCA with explicit 3-cluster test** (Gaussian mixture model with k=3 or simple K-means + silhouette) = karyotype calling.
3. **Linkage disequilibrium heatmap** across the inversion region with chromosome-wide background = recombination-suppression confirmation. Network-based LD analyses (Kemppainen et al. 2015) for chromosome-scale screens.
4. **dXY / FST / Tajima's D** in sliding windows comparing arrangements = age and selection signature.
5. **Genotype–environment association (GEA)** or genotype–phenotype association = selection regime / function. The Huang et al. 2020 sunflower workflow is the template.
6. **Family-aware GLMM** (lme4 or BLUPF90 style, family random effect) for any per-feature trait test in a hatchery cohort.
7. **PERMANOVA on the multivariate feature space** (Anderson 2001 *Austral Ecology*; PERMANOVA-S; vegan::adonis2) = omnibus test that the karyotypes differ jointly across features. **PERMDISP** (vegan::betadisper) accompanies it to distinguish centroid-shift from spread-shift.
8. **Mantel / partial-Mantel** for distance-based confound checks.
9. **DAPC** (Jombart et al. 2010 *BMC Genet.*; Thia 2023 *Mol. Ecol. Resour.* guidelines; DAPCy 2025 for scale) = supervised classification with cross-validation. Use this *instead of* LDA when n is large relative to p, but **prune sex-linked and inversion-linked loci before** computing the global DAPC, or fold inversion karyotype as a single "super-locus" — this is the explicit Thia 2023 recommendation.
10. **Random Forest** = classification when feature interactions matter and you want feature-importance ranking; sPCA (spatial PCA) when geography matters.
11. **Substructure-discovery ladder**: silhouette → gap statistic → bootstrap hierarchical clustering with multi-method concordance (your plan) is exactly correct and matches recent best practice.

---

## Details — The Recommended Six-Axis Hierarchy

### Top-level structure

```
Axis 1. EVIDENCE LAYER (gate, A–D) — Discovery
Axis 2. STRUCTURE — descriptive
Axis 3. MECHANISM / ORIGIN — molecular
Axis 4. EVOLUTIONARY HISTORY — temporal
Axis 5. POPULATION DYNAMICS — frequency
Axis 6. FUNCTION & APPLICATION — phenotypic / breeding
```

The carving is principled because:
- Each axis answers a *different* reviewer question (see §"Reviewer key questions" below).
- The axes are *biologically independent*: a young inversion (Axis 4) can be at any size (Axis 2) with any cargo (Axis 6); a sub-telomeric inversion (Axis 2) can have any age and any frequency.
- Each axis maps to a *distinct evidentiary base* (sequence vs SNP cluster vs LD vs cline vs phenotype).
- This six-axis split matches the implicit decomposition in Berdan et al. 2023 and Wellenreuther & Bernatchez 2018 Table 1, so reviewers will recognise the framework.

### Per-axis canonical question, literature precedent, statistical test, and your data input

#### Axis 1 — Evidence Layer (A–D)

- **Canonical question**: How confident are we that this is a real inversion, and what kind of evidence supports it?
- **Literature precedent**: Mérot, Oomen, Tigano & Wellenreuther 2020 *Trends Ecol. Evol.* (SV roadmap); Mérot 2020 *Mol. Ecol.* commentary; Li & Ralph 2019 *Genetics*; Reeve et al. 2023 *Mol. Ecol.* for the Littorina detection-method ladder.
- **Sub-axes**:
  - 1.1 Layer A: assembly/long-read/Hi-C/karyotype confirmation
  - 1.2 Layer B: three-cluster local-PCA + heterozygote-excess LD pattern
  - 1.3 Layer C: long-distance LD block alone
  - 1.4 Layer D: phenotypic/cline outlier alone
- **Test/descriptor**: cluster-quality metrics (silhouette, BIC of k=3 GMM); LD-block boundary clarity; presence/absence of long-read split-alignment.
- **Your data input**: Q7 existence A–D + tier; bundle export already has this.

#### Axis 2 — Structure

- **Canonical question**: What does this inversion look like physically?
- **Sub-axes**:
  - 2.1 Size class (small <1 Mb / intermediate 1–5 Mb / large 5–10 Mb / supergene-scale ≥10 Mb)
  - 2.2 Centromeric position (paracentric / pericentric)
  - 2.3 Chromosomal-arm position (sub-telomeric / interstitial / peri-centromeric)
  - 2.4 Topology (single / nested / overlapping / part of a complex)
  - 2.5 Boundary clarity (sharp breakpoint coordinates / repeat-blurred / ambiguous)
- **Literature precedent**: Wellenreuther & Bernatchez 2018 Table 1; Connallon & Olito 2022 *Mol. Ecol.* (size as a signature); Harringmeyer & Hoekstra 2022 *Nat. Ecol. Evol.* (breakpoint positions); Schaeffer 2008 / Wallace, Detweiler & Schaeffer 2011 (nesting); Corbett-Detig & Hartl 2012.
- **Test/descriptor**: physical size in Mb and as % of chromosome; coordinate of midpoint relative to centromere/telomere; nesting graph; boundary CI from breakpoint reads.
- **Your data input**: Q1 (shape/composite of local-PCA), Q3 (boundaries L/R repeat density, GC, repeat class).

#### Axis 3 — Mechanism / Origin

- **Canonical question**: What molecular mechanism produced this inversion, and what is at the breakpoints?
- **Sub-axes**:
  - 3.1 Breakpoint substrate (NAHR-via-TE / NAHR-via-segmental-duplication / NHEJ / FoSTeS-MMBIR / unknown)
  - 3.2 Repeat-class enrichment at breakpoints (LTR / LINE / DNA-TE / tandem-repeat)
  - 3.3 Evidence of recurrent or single origin (chimeric repeats vs unique sequence)
  - 3.4 Recombination-suppression profile (full / partial / breakpoint-leaky); gene-flux signature
- **Literature precedent**: Cáceres et al. 1999 *PNAS*; Casals et al. 2003; Delprat et al. 2009 *PLoS One* (Drosophila Galileo); Lakich et al. 1993 *Nat. Genet.* (factor VIII NAHR); Hastings et al. 2009 *Nat. Rev. Genet.* (MMBIR); Atlantic cod TE-tandem-mediated inversions (*GBE* 2025, Tørresen group); Harringmeyer & Hoekstra 2022 (deer mouse inverted-repeat flanks); Schaeffer & Anderson 2005 *Genetics* (gene-flux); Korunes & Noor 2019; Berdan, Blanckaert, Butlin & Bank 2021 *PLoS Genet.*; Kuang et al. 2024 (your catfish citation).
- **Test/descriptor**: repeat-density Z-score on each side; same-family/opposite-orientation flanking-element check; microhomology length at junction; per-breakpoint gene-disruption test (Corbett-Detig & Hartl 2012's permutation against random breakpoint model).
- **Your data input**: Q4 (TE family / GO enrichment / tandem-repeat at breakpoints); Q3 (repeat density, GC, repeat class).

#### Axis 4 — Evolutionary History

- **Canonical question**: How old is this inversion, which arrangement is ancestral, and is it shared across populations/species?
- **Sub-axes**:
  - 4.1 Age (point estimate ± CI in My or generations; method: dXY, MSMC2, ABC, IBDseq)
  - 4.2 Polarisation (ancestral vs derived arrangement; outgroup-based Dollo)
  - 4.3 Polymorphism class (Type I transient / Type II balanced — Faria et al. 2019)
  - 4.4 Recurrence / shared-with-sister-species (e.g. Atlantic cod inversions trans-Atlantic; Littorina ancestral polymorphism)
- **Literature precedent**: Hudson-Kreitman-Aguadé 1987; Charlesworth & Charlesworth 2010 textbook; Faria et al. 2019 *Trends Ecol. Evol.*; Berg et al. 2017 *Heredity* (cod trans-Atlantic); Reeve et al. 2023 *Mol. Ecol.* (Littorina); Fuller, Haynes, Richards & Schaeffer 2017–2018; Stankowski et al. 2024 *Science* (Littorina live-bearing); Lamichhaney et al. 2016 (ruff dating); Tuttle et al. 2016 / Sun et al. 2018 (sparrow); Jay et al. 2018 (Heliconius numata).
- **Test/descriptor**: dXY between arrangements at near-breakpoint windows / 2μ; Tajima's D inside vs outside; Dollo parsimony with ≥2 outgroups; population-tree topology weighting (Twisst).
- **Your data input**: Q5 (age / origin / Tajima D / dXY / Dollo).

#### Axis 5 — Population Dynamics

- **Canonical question**: How is the inversion distributed in the population, and is its frequency consistent with neutrality or selection?
- **Sub-axes**:
  - 5.1 Frequency class (rare <5% / intermediate 5–50% / common >50%)
  - 5.2 HWE deviation and direction (heterozygote excess → balancing/overdominance; deficit → assortment/drift/structure)
  - 5.3 Family/pedigree confound class (cleared via relatedness pruning / contaminated / unresolved)
  - 5.4 Geographic/environmental cline (present / absent / parallel-across-populations)
  - 5.5 Within-arrangement substructure (single / multiple haplotype lineages within one arrangement)
- **Literature precedent**: Charlesworth 2016 (balancing selection theory); Kapun & Flatt 2019 (clinal); Mérot et al. 2018 *Proc. R. Soc. B* (Coelopa cline); Fouet et al. 2012 (Anopheles 2La); Wellenreuther et al. 2017 *J. Evol. Biol.* (clinal damselfly); Mérot et al. 2020 *Mol. Biol. Evol.* (Coelopa within-arrangement substructure); Reeve et al. 2023 *Mol. Ecol.* (Littorina); Cridland et al. 2013 / Huang et al. 2014 (DGRP relatedness as warning template). For aquaculture-specific family confound see Robledo et al. 2018 *Trends Genet.* on aquaculture genomic selection.
- **Test/descriptor**: HWE chi-square per population; F_IS within karyotype; ngsRelate-pruned re-test; latitudinal/environmental GLMM with karyotype as response; silhouette + gap-statistic + bootstrap-hierarchical-clustering ladder for within-karyotype substructure.
- **Your data input**: Q6 (frequency / HWE / family-linkage / polymorphism class); ngsRelate; NGSadmix Q-matrix; family pedigree.

#### Axis 6 — Function & Application

- **Canonical question**: What does this inversion do phenotypically, and is it useful as a breeding marker?
- **Sub-axes**:
  - 6.1 Cargo class (gene-rich / gene-poor / contains-immune-or-MHC / contains-reproductive / contains-growth-or-feed-conversion)
  - 6.2 GO/KEGG enrichment of cargo
  - 6.3 Trait association (presence/absence; effect size; sex-specific)
  - 6.4 Selection-regime inference (locally adaptive / overdominant / sweep / neutral / unresolved)
  - 6.5 Marker tier (Tier 1 breakpoint-PCR / Tier 2 private-indel / Tier 3 SNP-tag / Tier 4 exploratory)
- **Literature precedent**: Hoffmann & Rieseberg 2008; Wellenreuther & Bernatchez 2018; Schwander, Libbrecht & Keller 2014; Gutiérrez-Valencia et al. 2021; Lowry & Willis 2010 (Mimulus); Lamichhaney 2016; Tuttle 2016; Sun 2018; Joron 2011; Jay 2018, 2021; Todesco 2020; Hager et al. 2022 (deer mouse multi-trait); Berg et al. 2016/2017 (cod); Han et al. 2020 (herring); Pearse 2019 (rainbow trout Omy05); Bertolotti 2020 (Atlantic salmon); Robledo et al. 2018; Yáñez Cooke Aquaculture imputation; Wellenreuther et al. 2024 marine-management framing.
- **Test/descriptor**: GO-enrichment hypergeometric test (Bonferroni or FDR-corrected); GLMM with karyotype as fixed effect, family as random; per-trait omnibus PERMANOVA on multi-trait phenotype; classification accuracy of Tier-3 SNP panel via leave-one-out cross-validation.
- **Your data input**: Q1 cargo descriptors; GHSL gene lists; phenotype tables; family pedigree; the 14-axis tier view's "burden" and "confidence tier" entries.

---

## How your three existing schemes map onto the proposed unified framework

### Q1–Q7 registry mapping

| Q-block | Proposed axis | Status in unified scheme |
| --- | --- | --- |
| Q1 (shape/composite of local-PCA) | Axis 1 sub-axis 1.2 + Axis 2 sub-axis 2.5 | Split — the *signal-quality* part is Evidence; the *boundary-clarity* part is Structure |
| Q2 (internal dynamics / decomposition / recombinant map) | Axis 5.5 + Axis 3.4 | Mostly substructure (5.5); small overlap with recombination-suppression profile (3.4) |
| Q3 (boundaries L/R repeat density, GC, repeat class) | Axis 2.5 + Axis 3.1–3.2 | Split — boundary clarity is Structure; repeat substrate is Mechanism |
| Q4 (TE family / GO enrichment / tandem-repeat at breakpoints) | Axis 3.2 + Axis 6.2 | Split — breakpoint-TE is Mechanism; cargo-GO is Function |
| Q5 (age / origin / Tajima D / dXY / Dollo) | Axis 4.1–4.2 + Axis 6.4 | Mostly Evolutionary History; selection-regime inference (Tajima D) is Function |
| Q6 (frequency / HWE / family-linkage / polymorphism class) | Axis 5.1–5.4 + Axis 4.3 | Mostly Population Dynamics; polymorphism class (Type I/II) is Evolutionary History |
| Q7 (existence A–D + independence + tier) | Axis 1 (entirely) | Becomes the gate |

**What is redundant** in your current registry:
- Q1 "shape/composite" overlaps with Q7's "existence layer" — both are evidence/quality metrics.
- Q3 "repeat density at boundaries" appears again as Q4 "TE family at breakpoints" — these are the same evidence; merge into Axis 3.2.
- The 14-axis tier view's "boundary quality" and "group validation" duplicate Q7 layer logic and Q1 signal-quality.
- The "stats profile" 5-bucket grouping (Breakpoint architecture / Genomic composition / Functional cargo / Population variation / Breeding utility) re-organises the Q-blocks but adds no new content; retire as an organising layer, keep only as figure-panel groupings.

**What is missing** from your current scheme:
- A clean **Mechanism axis** distinct from Structure — Q3 and Q4 currently smear breakpoint-substrate (mechanism) with size/position (structure) and cargo (function).
- An explicit **polarisation / ancestral-vs-derived** sub-axis at Axis 4.2 — your Q5 has Dollo but it's not surfaced as a first-class sub-axis.
- An explicit **recurrent / shared-with-sister-species** flag for the *C. macrocephalus* portability — Axis 4.4.
- An explicit **selection-regime inference** field at Axis 6.4 distinct from cargo content.

**What should be merged**:
- Merge Q1 signal-quality + Q7 existence layer + 14-axis-tier "boundary quality" + "group validation" → Axis 1 alone, with a single per-inversion **Evidence Card** that lists all of: long-read assembly, three-cluster PCA, LD block, GEA outlier.
- Merge Q3 repeats + Q4 breakpoint-TE → Axis 3 (Mechanism) with a single **Breakpoint Architecture Card**.
- Merge Q4 cargo-GO + 5-bucket "Functional cargo" + 14-axis "burden" → Axis 6.1–6.2 (Function/Cargo Card).
- Merge Q6 frequency + HWE + family + 14-axis "polymorphism class" → Axis 5 (Population Dynamics Card).

### 14-axis tier view mapping

The 14 axes collapse cleanly:
- Existence layers A–D → Axis 1.
- Boundary quality → Axis 1.2 + Axis 2.5.
- Group validation → Axis 1.2 (cluster-quality metrics).
- Internal structure → Axis 5.5.
- Recombinant class → Axis 3.4.
- Family linkage → Axis 5.3.
- Polymorphism class → Axis 4.3.
- Mechanism → Axis 3.
- Age → Axis 4.1.
- Burden → Axis 6.1 (gene-rich/gene-poor) and Axis 3.4 (deleterious-load context).
- Confidence tier → Axis 1 final score (an aggregate descriptor).

### 5-bucket stats profile mapping

- Breakpoint architecture → Axis 3.
- Genomic composition → Axis 2 + Axis 6.1.
- Functional cargo → Axis 6.1–6.3.
- Population variation → Axis 5.
- Breeding utility → Axis 6.5.

The 5-bucket scheme is essentially Axes 3, (2+6.1), 6, 5, 6.5 — i.e. it conflates Structure and Cargo and double-counts Function. **Retire it as the primary scheme; keep it as the figure-panel layout** in the bundle export so each per-inversion result panel has the same five visual groupings.

---

## Reviewer "key questions" — every inversion gets answered on these eight points

1. **Is it real?** → Axis 1 (Evidence layer A–D + cluster quality)
2. **How big and where?** → Axis 2 (size class, paracentric/pericentric, sub-telomeric/interstitial, simple/nested)
3. **What's at the breakpoints?** → Axis 3 (NAHR/NHEJ; TE family; repeat density)
4. **How old, and is it ancestral or derived?** → Axis 4 (age estimate ± CI; outgroup-based polarisation)
5. **What is its frequency and is it neutral?** → Axis 5 (rare/intermediate/common; HWE; cline)
6. **What does it carry, and does it associate with a trait?** → Axis 6.1–6.3 (cargo content; GO; trait association)
7. **What evolutionary force is maintaining it?** → Axis 6.4 (locally adaptive / overdominant / sweep / neutral)
8. **Is it useful as a breeding marker?** → Axis 6.5 (Tier 1–4)

These eight questions are precisely the eight reporting fields in Wellenreuther & Bernatchez 2018 Table 1 (with #1 implicit and #8 added for breeding context). If your manuscript Table 1 has these eight columns, reviewers will read it as canonical.

---

## Species-specific vs species-portable axes — the bundle-export design point

**Fully portable (always reportable for any species)**: Axis 1, 2, 3, 4, 5.1–5.2, 6.1–6.4. These are all sequence-derivable and need only a reference assembly + population resequencing.

**Conditionally portable (require species-specific extras)**:
- Axis 5.3 (family confound) — only relevant in pedigreed/hatchery cohorts. For wild-cohort species (future *C. macrocephalus* wild population study), this sub-axis is N/A and should be auto-suppressed.
- Axis 5.4 (cline) — requires multi-population sampling; flag as N/A when the dataset is single-population.
- Axis 6.3 (trait association) — requires phenotyping; the bundle export should auto-flag as "no phenotype data" when the species has none.
- Axis 6.5 (marker tier) — only relevant when the species has a breeding programme; for wild population genomics studies, the panel should be re-labeled "Diagnostic-marker portability" with the same Tier 1–4 ladder.

**Truly species-specific**: nothing in the framework. The framework is fully generalisable; the *content* of cargo (e.g. catfish-specific immunity gene families) is species-specific but slots into the same Axis 6.1 schema.

The bundle-export design implication: each per-inversion panel has six sub-panels (one per axis) plus an optional seventh (cohort-context: family/pedigree visualisation) that auto-suppresses when the cohort has no pedigree. The panel template is fixed; the data filling each sub-panel is species-specific.

---

## Recommendations

### Stage 1 — Immediate (manuscript-critical, do in next 1–2 weeks)

1. **Adopt the six-axis framework as the manuscript's organising structure**, with Wellenreuther & Bernatchez 2018 Table 1 as the per-inversion summary template (8 columns: ID, size, position, breakpoint mechanism, age, frequency-cline class, cargo+trait, marker tier).
2. **Collapse Q1–Q7, the 14-axis tier view, and the 5-bucket stats profile into Axis 1–6** using the mapping table above. Document this as a single internal "axis dictionary" markdown file so future-you and *Macrocephalus* future-you don't re-fragment.
3. **Replace the layered Q-block reporting in figures with Evidence Card + Breakpoint Card + Population Dynamics Card + Function Card + Marker Card** — exactly five panels per inversion, mirroring the canonical reviewer questions.
4. **Promote polarisation (ancestral vs derived, Dollo with ≥2 outgroups) to a first-class sub-axis 4.2** in the export.
5. **Make family-confound resolution (Axis 5.3) explicit** — every panel reports both the raw and ngsRelate-pruned local-PCA so the reviewer can see that the inversion signal is not a family-pack artefact. This is a *specific* aquaculture vulnerability that no Drosophila / Anopheles paper has had to defend; pre-empt it.

### Stage 2 — Next month (portability for *Macrocephalus* and beyond)

6. **Refactor the bundle export to be axis-keyed rather than Q-keyed**. The export YAML/JSON should have six top-level keys (`evidence`, `structure`, `mechanism`, `history`, `population`, `function_application`) with sub-axis fields under each, and a `cohort_context` block that records which sub-axes are N/A for the current dataset.
7. **Auto-suppress N/A panels** based on `cohort_context` (e.g. no pedigree → suppress Axis 5.3; single population → suppress Axis 5.4) and replace with a "not applicable for this study design" stamp.
8. **Build a one-page "Inversion Spec Sheet" (8 fields, A4-half-page) for each catalogued inversion**, identical in layout for catfish and any future species. Reviewers will recognise this as analogous to the JTTLR / GenBank-style per-feature sheet.

### Stage 3 — Strategic (for the *Nature Communications* readability)

9. **Lead the manuscript with one canonical worked example** (e.g. the LG28 sub-telomeric case) walking the reader through Axis 1–6 in order. This establishes the framework before the reader meets the catalogue.
10. **Place the comparative table (Wellenreuther & Bernatchez 2018-style, 8 columns) as Table 1**; relegate the Q-registry to supplementary material and explicitly state in Methods that the registry is the implementation back-end of the six-axis framework.
11. **For figures with the LD block visualisation, use the Wellenreuther & Bernatchez 2018 Figure 1 conventions** (paracentric/pericentric schematic top, LD heatmap below, three-cluster PCA inset) — reviewers know this layout.

### Benchmarks / thresholds that would change these recommendations

- If ngsRelate-pruning eliminates >30% of an inversion signal: re-classify that inversion to Evidence layer C (LD-block-only) and treat as exploratory in the manuscript.
- If <90% of inversions can be assigned to Layer A or B: do not report a comprehensive catalogue; report only the validated subset and frame the rest as preliminary.
- If size_class shows the inversion is <500 kb and inside a single gene: reconsider whether it is an "inversion" in the supergene-tradition sense or a small SV closer to the structural-variant call set; the literature treats sub-Mb inversions less canonically (cf. Yano et al. 2025 *Trends Ecol. Evol.* "Beyond supergenes").
- If polarisation is unresolved (no ≥2 outgroups available): explicitly downgrade Axis 4.2 to "unpolarised" rather than guessing.

---

## Caveats — where the literature is unsettled, contested, or weak

1. **Size thresholds are not consensus** — small/intermediate/large boundaries vary between papers (Hoffmann & Rieseberg use one set, Connallon & Olito 2022 use relative size, Huang et al. 2024 *Plant Biotechnol. J.* propose percentage-of-chromosome thresholds). State your thresholds explicitly.
2. **The "inversion vs supergene" boundary is fuzzy** — some authors reserve "supergene" for inversions controlling discrete morphs (Schwander et al. 2014); others use it for any large recombination-suppressed haplotype block (Todesco et al. 2020 "haploblocks"; Berdan et al. 2023). Pick one convention and say so.
3. **Underdominance vs locally-adaptive vs overdominance is overdetermined by data alone** — Berdan et al. 2023 explicitly warn that "most patterns are overdetermined". Be careful claiming a specific selection regime without convergent multiple lines of evidence (cline + heterozygote excess + GO enrichment + age, ideally all four).
4. **Age estimation from dXY is biased by gene flux** — Schaeffer & Anderson 2005 / Korunes & Noor 2019 / Berdan et al. 2021 show gene conversion homogenises mid-inversion regions; restrict dXY to near-breakpoint windows and note the assumption.
5. **The fragile-breakage model is not universally accepted** — Sankoff & Trinh 2004 / 2007 disputed the original Pevzner & Tesler 2003 finding; the consensus has shifted toward "fragile regions exist but turnover" (Alekseyev & Pevzner 2010). Frame breakpoint hotspots as "consistent with fragile-breakage regions in this lineage" rather than asserting universality.
6. **Local-PCA can miss inversions whose population frequency is too low or too imbalanced** — Li & Ralph 2019 explicitly state this; rare or near-fixed inversions may be invisible. Report your detection-power estimate.
7. **Local-PCA can also produce false positives from family structure or population structure**; ngsRelate-pruned re-runs are essential. Your awareness of this in the cohort design is a strength to call out.
8. **The "sub-telomeric enrichment" of inversion breakpoints** (Harringmeyer & Hoekstra 2022) is so far documented in mammals; the catfish case is consistent but the n is small and the literature is not yet settled on whether it generalises across vertebrates.
9. **The TE-mediated NAHR mechanism is robustly documented in Drosophila (Galileo), Atlantic cod (hAT/MITE tandems), deer mouse (long inverted repeats)**, but for many species the breakpoint resolution is too coarse to discriminate NAHR from NHEJ. State your confidence honestly.
10. **The Tier 1–4 marker scheme is your invention** — the literature does not have a standardised inversion-marker tier system. This is fine and indeed a contribution; cite the conceptual lineage (cattle/pig QTL marker design; salmonid marker-assisted selection — e.g. Yáñez et al. work; Robledo et al. 2018 *Trends Genet.*; Bertolotti et al. 2020 *Nat. Commun.*) and present the four-tier scheme as a generalised framework adapted from those, not as field-standard.
11. **Hardy–Weinberg deviation interpretation is multi-causal** — heterozygote excess can arise from heterokaryotype advantage *or* from secondary contact / admixture *or* from family structure *or* from null alleles in a SNP panel. Always rule out the alternatives.
12. **The 14-author Berdan et al. 2023 paper itself notes** that disentangling adaptive from non-adaptive inversion evolution is currently the central open problem in the field. Your manuscript should acknowledge that even after all six axes are populated, attribution of selection regime remains hypothesis, not certainty, for any individual inversion.

---

This six-axis framework is the literature-grounded backbone you need. It is orthogonal across axes, exhaustive across the eight reviewer questions, mappable from your existing Q1–Q7 + 14-axis + 5-bucket schemes without losing information, portable from *gariepinus* to *macrocephalus* to other aquaculture species, and explicit about where the literature is unsettled. Adopt it as the manuscript's organising structure and the bundle-export's data schema in parallel.