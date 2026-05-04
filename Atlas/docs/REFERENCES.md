# REFERENCES — MS_Inversions_North_african_catfish

Flat list of citations the manuscript methods will pull from, organized by
pipeline stream. **Single source of truth** — when writing methods or
building the eventual atlas references page, read from this file.

Format: `Author Year. Title. Journal vol:pages. doi:xxx.`

Entries marked `[CITED]` are confirmed in the working draft. Entries
marked `[CANDIDATE]` are tools/methods used in the pipeline that should
be cited when methods are written. Entries marked `[BACKGROUND]` are
foundational references — cite only if the methods section defines the
underlying concept.

Cohort: 226 pure *C. gariepinus* hatchery samples on LANTA. Citations
specific to the F1 hybrid genome assembly paper or the *C. macrocephalus*
wild cohort do not belong in this file.

---

## Variant calling, genotype likelihoods, low-coverage framework

[CANDIDATE] **Korneliussen TS, Albrechtsen A, Nielsen R. 2014.**
ANGSD: Analysis of Next Generation Sequencing Data.
*BMC Bioinformatics* 15:356. doi:10.1186/s12859-014-0356-4.

[CANDIDATE] **Meisner J, Albrechtsen A. 2018.**
Inferring population structure and admixture proportions in low-depth NGS
data. *Genetics* 210:719–731. doi:10.1534/genetics.118.301336.
*(PCAngsd)*

[CANDIDATE] **Skotte L, Korneliussen TS, Albrechtsen A. 2013.**
Estimating individual admixture proportions from next generation
sequencing data. *Genetics* 195:693–702.
doi:10.1534/genetics.113.154138.
*(NGSadmix)*

[CANDIDATE] **Korneliussen TS, Moltke I. 2015.**
NgsRelate: a software tool for estimating pairwise relatedness from
next-generation sequencing data. *Bioinformatics* 31:4009–4011.
doi:10.1093/bioinformatics/btv509.

[CANDIDATE] **Garcia-Erill G, Albrechtsen A. 2020.**
Evaluation of model fit of inferred admixture proportions.
*Mol Ecol Resour* 20:936–949. doi:10.1111/1755-0998.13171.
*(evalAdmix)*

[CANDIDATE] **Andrade P, Lopes S, Gomes I, Boschiero C, et al. 2019.**
NAToRA, a relatedness-pruning method to enhance the inference of
population history. *Bioinformatics* 35:5079–5081.
*(NAToRA — verify exact citation; double-check year/journal)*

---

## Heterozygosity, ROH, nucleotide diversity (θπ track)

[CITED] **Korunes KL, Samuk K. 2021.**
pixy: Unbiased estimation of nucleotide diversity and divergence in the
presence of missing data. *Mol Ecol Resour* 21:1359–1368.
doi:10.1111/1755-0998.13326.
*(Per-site θπ normalization in the presence of missing data. The
pestPG `tP / nSites` correction in STEP_TR_A follows this convention.)*

[BACKGROUND] **Nei M, Li WH. 1979.**
Mathematical model for studying genetic variation in terms of restriction
endonucleases. *PNAS* 76:5269–5273.
*(Foundational π. Cite only if the methods section defines π formally.)*

[CANDIDATE] **Vieira FG, Albrechtsen A, Nielsen R. 2016.**
Estimating IBD tracts from low coverage NGS data.
*Bioinformatics* 32:2096–2102. doi:10.1093/bioinformatics/btw117.
*(ngsF-HMM, used in MODULE_3 ROH track.)*

---

## Phasing, haplotypes (founder-pack track)

[CANDIDATE] **Patterson M, Marschall T, Pisanti N, et al. 2015.**
WhatsHap: weighted haplotype assembly for future-generation sequencing
reads. *J Comput Biol* 22:498–509. doi:10.1089/cmb.2014.0157.
*(WhatsHap phase blocks: 44–819 bp on our short-read data.)*

---

## Local PCA, inversion discovery (dosage + GHSL streams)

[CANDIDATE] **Li H, Ralph P. 2019.**
Local PCA shows how the effect of population structure differs along
the genome. *Genetics* 211:289–304. doi:10.1534/genetics.118.301747.
*(Foundational local-PCA-for-inversions framework. Cite for the
methodological lineage of MODULE_2A and MODULE_2E approaches.)*

[CANDIDATE — TBD] **GHSL methodology citation.**
*(If a published method underlies the Group-level Haplotype Similarity
by Locus formulation, add here. Otherwise this is novel methodology
described in the manuscript itself.)*

---

## Structural variants

[CANDIDATE] **Rausch T, Zichner T, Schlattl A, Stütz AM, Benes V,
Korbel JO. 2012.**
DELLY: structural variant discovery by integrated paired-end and
split-read analysis. *Bioinformatics* 28:i333–i339.
doi:10.1093/bioinformatics/bts378.

[CANDIDATE] **Chen X, Schulz-Trieglaff O, Shaw R, et al. 2016.**
Manta: rapid detection of structural variants and indels for germline
and cancer sequencing applications. *Bioinformatics* 32:1220–1222.
doi:10.1093/bioinformatics/btv710.

---

## Conservation, deleterious variants (MODULE_CONSERVATION)

[CANDIDATE] **Cingolani P, Platts A, Wang LL, et al. 2012.**
A program for annotating and predicting the effects of single nucleotide
polymorphisms, SnpEff. *Fly (Austin)* 6:80–92.
doi:10.4161/fly.19695.

[CANDIDATE] **Vaser R, Adusumalli S, Leng SN, Sikic M, Ng PC. 2016.**
SIFT missense predictions for genomes. *Nat Protoc* 11:1–9.
doi:10.1038/nprot.2015.123.
*(SIFT4G.)*

[CANDIDATE] **Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. 2010.**
Detection of nonneutral substitution rates on mammalian phylogenies.
*Genome Res* 20:110–121. doi:10.1101/gr.097857.109.
*(GERP++.)*

[CANDIDATE] **Armstrong J, Hickey G, Diekhans M, et al. 2020.**
Progressive Cactus is a multiple-genome aligner for the thousand-genome
era. *Nature* 587:246–251. doi:10.1038/s41586-020-2871-y.

[CANDIDATE — TBD] **VESM_650M citation.**
*(Add when settled which VESM/ESM variant the pipeline uses.)*

---

## Aquaculture / hatchery context (manuscript framing)

[CANDIDATE — TBD] **C. gariepinus aquaculture / pedigree references.**
*(Standard SE Asian catfish aquaculture references go here. Pull from
prior chapter drafts when methods is written.)*

---

## Pipeline / tooling (cite if Methods describes software stack)

[CANDIDATE] **Li H. 2018.** Minimap2: pairwise alignment for
nucleotide sequences. *Bioinformatics* 34:3094–3100.
doi:10.1093/bioinformatics/bty191.

[CANDIDATE] **Danecek P, Bonfield JK, Liddle J, et al. 2021.**
Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
doi:10.1093/gigascience/giab008.

---

## NOT to cite in manuscript

These came up during pipeline development but are not citable sources
for a methods section. Useful in the lab notebook / code comments only.

- ANGSD GitHub issue #329 (the θπ scaling thread).
  *Trace of how the pestPG `tP` sum-vs-per-site bug was identified.
  Lives in the STEP_TR_A code comment, not the bibliography.*

- Any biorxiv preprint that has since been published.
  *Always cite the peer-reviewed version. The pixy preprint
  (`bioRxiv 10.1101/2020.06.27.175091`) became Korunes & Samuk 2021,
  cited above.*

- StackOverflow / Biostars / pop-genomics-discord threads.

---

## Working notes

- Five `[CANDIDATE — TBD]` entries flagged above need verification or
  filling in before submission. Do this pass after LG28 results land,
  not now.
- DOIs are present where verified. Entries without DOI need a quick
  PubMed check during the methods-writing pass.
- This file is plain Markdown — no rendering, no linking to a bib
  manager. When a manuscript draft starts, export to BibTeX with
  whichever tool the journal expects (most use either RIS or
  ENW import from Zotero / Mendeley).
